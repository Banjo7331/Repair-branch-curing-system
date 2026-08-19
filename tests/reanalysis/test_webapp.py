"""The review server, and the one thing it is allowed to change.

Acknowledging is the only write in this package that makes the system *quieter*
— the surfacing policy suppresses every future change carrying an acknowledged
fingerprint. So the tests here are mostly about the conditions on that write:
it has to be attributed, it has to hit an event that is actually open, and it
must not be replayable by a refresh.

The rest are the ordinary obligations of anything that renders somebody else's
strings into HTML and calls itself a review surface.
"""

from __future__ import annotations

import http.client
import logging
import threading
import urllib.error
import urllib.parse
import urllib.request
from datetime import UTC, datetime
from http.server import ThreadingHTTPServer
from pathlib import Path

import pytest

from repairbench.reanalysis import operations, webapp
from repairbench.reanalysis.dashboard import escape
from repairbench.reanalysis.ledger import EventStatus
from repairbench.reanalysis.operations import OperationError
from repairbench.reanalysis.routing import ReviewQueue, Urgency
from repairbench.reanalysis.store import (
    JsonAssessmentStore,
    JsonCaseRepository,
    StoredEvent,
    StoreError,
)
from repairbench.reanalysis.webapp import ReviewApp
from repairbench.reanalysis.world import DriftAxis, Pin

NOW = datetime(2026, 8, 17, 12, 0, tzinfo=UTC)


def phenotype() -> Pin:
    return Pin(axis=DriftAxis.PHENOTYPE, version="hpo-day-1", digest="0" * 64)


def event(event_id: str = "evt-1", *, status: EventStatus = EventStatus.OPEN) -> StoredEvent:
    return StoredEvent(
        event_id=event_id,
        variant_key="PLUSG-c158",
        fingerprint=f"fp-{event_id}",
        status=status,
        summary="loss_of_function → undetermined: the mechanism no longer resolves",
        urgency=Urgency.HIGH,
        queue=ReviewQueue.CLINICAL_SIGNOUT,
    )


@pytest.fixture
def repository(tmp_path: Path) -> JsonCaseRepository:
    return JsonCaseRepository(tmp_path)


@pytest.fixture
def app(repository: JsonCaseRepository) -> ReviewApp:
    return ReviewApp(repository)


def registered(
    repository: JsonCaseRepository,
    case_id: str = "NICU-014",
    *,
    events: list[StoredEvent] | None = None,
) -> None:
    ledger = repository.register(case_id, ["PLUSG-c158"], phenotype())
    ledger.events = list(events or [])
    ledger.last_examined_at = NOW
    repository.save(ledger)


def reload_event(repository: JsonCaseRepository, case_id: str = "NICU-014") -> StoredEvent:
    return repository.get(case_id).events[0]  # type: ignore[return-value]


# --------------------------------------------------------------------------
# Acknowledging: the one write
# --------------------------------------------------------------------------


def test_an_acknowledgement_records_who_made_it(app: ReviewApp, repository: JsonCaseRepository):
    """Not a flag. A suppression that nobody can be asked about later is one an
    incident review cannot reconstruct."""
    registered(repository, events=[event()])

    app.acknowledge("NICU-014", "evt-1", "B. Kowalski", "confirmed against ClinGen")

    stored = reload_event(repository)
    assert stored.status is EventStatus.ACKNOWLEDGED
    assert stored.acknowledged_by == "B. Kowalski"
    assert stored.acknowledged_note == "confirmed against ClinGen"
    assert stored.acknowledged_at is not None


def test_an_anonymous_acknowledgement_is_refused(app: ReviewApp, repository: JsonCaseRepository):
    registered(repository, events=[event()])

    with pytest.raises(ValueError, match="needs the name"):
        app.acknowledge("NICU-014", "evt-1", "   ", "")

    assert reload_event(repository).status is EventStatus.OPEN


def test_a_note_is_optional_but_a_reviewer_is_not(app: ReviewApp, repository: JsonCaseRepository):
    """A reviewer may have nothing to add. Requiring prose would produce prose
    like "ok", which is worse than an empty field because it looks like content."""
    registered(repository, events=[event()])

    app.acknowledge("NICU-014", "evt-1", "B. Kowalski", "")

    assert reload_event(repository).status is EventStatus.ACKNOWLEDGED


def test_acknowledging_an_unknown_event_is_refused_rather_than_ignored(
    app: ReviewApp, repository: JsonCaseRepository
):
    """A form post that silently did nothing would leave the reviewer believing
    they had signed something off."""
    registered(repository, events=[event()])

    with pytest.raises(StoreError, match="not an open event"):
        app.acknowledge("NICU-014", "evt-does-not-exist", "B. Kowalski", "")


def test_a_superseded_event_cannot_be_acknowledged(
    app: ReviewApp, repository: JsonCaseRepository
):
    """Signing for a transition a later run already overtook is the exact mistake
    the ledger supersedes events to prevent."""
    registered(repository, events=[event(status=EventStatus.SUPERSEDED)])

    with pytest.raises(StoreError, match="not an open event"):
        app.acknowledge("NICU-014", "evt-1", "B. Kowalski", "")


def test_acknowledging_twice_is_refused_the_second_time(
    app: ReviewApp, repository: JsonCaseRepository
):
    registered(repository, events=[event()])
    app.acknowledge("NICU-014", "evt-1", "B. Kowalski", "")

    with pytest.raises(StoreError):
        app.acknowledge("NICU-014", "evt-1", "Somebody Else", "")

    assert reload_event(repository).acknowledged_by == "B. Kowalski"


def test_an_acknowledgement_survives_being_written_and_read_back(
    app: ReviewApp, tmp_path: Path, repository: JsonCaseRepository
):
    """The whole point is that the next scheduled run, in a different process,
    does not raise this change again."""
    registered(repository, events=[event()])
    app.acknowledge("NICU-014", "evt-1", "B. Kowalski", "seen")

    fresh = JsonCaseRepository(tmp_path).get("NICU-014")

    assert fresh.acknowledged_fingerprints == frozenset({"fp-evt-1"})
    assert fresh.open_events() == ()


# --------------------------------------------------------------------------
# The pages
# --------------------------------------------------------------------------


def test_the_queue_links_into_each_case(app: ReviewApp, repository: JsonCaseRepository):
    registered(repository, events=[event()])

    assert '<a href="/case/NICU-014">' in app.queue()


def test_a_case_page_shows_the_form_only_for_open_events(
    app: ReviewApp, repository: JsonCaseRepository
):
    registered(repository, events=[event()])
    assert "<form" in app.case("NICU-014")

    app.acknowledge("NICU-014", "evt-1", "B. Kowalski", "")
    assert "<form" not in app.case("NICU-014")


def test_a_closed_event_keeps_its_reviewer_and_note_on_the_page(
    app: ReviewApp, repository: JsonCaseRepository
):
    """An append-only ledger whose notes nobody can find later is a ledger
    nobody had a reason to append to."""
    registered(repository, events=[event()])
    app.acknowledge("NICU-014", "evt-1", "B. Kowalski", "plan paused pending review")

    page = app.case("NICU-014")

    assert "B. Kowalski" in page
    assert "plan paused pending review" in page


def test_the_page_warns_that_acknowledging_suppresses_future_alerts(
    app: ReviewApp, repository: JsonCaseRepository
):
    """A decision with a longer reach than "mark as read" suggests, said at the
    moment somebody makes it."""
    registered(repository, events=[event()])

    assert "suppresses every future alert" in app.case("NICU-014")


def test_the_page_says_the_name_is_not_authenticated(
    app: ReviewApp, repository: JsonCaseRepository
):
    registered(repository, events=[event()])

    assert "does not authenticate" in app.case("NICU-014")


def test_an_unknown_case_is_refused(app: ReviewApp, repository: JsonCaseRepository):
    registered(repository)

    with pytest.raises(StoreError, match="no case"):
        app.case("NOT-A-CASE")


def test_a_reviewer_name_cannot_inject_markup(app: ReviewApp, repository: JsonCaseRepository):
    """Typed into a box by whoever has the page open."""
    registered(repository, events=[event()])

    app.acknowledge("NICU-014", "evt-1", "<img src=x onerror=alert(1)>", "")

    page = app.case("NICU-014").split("</style>", 1)[1]
    assert "<img src=x" not in page
    assert "&lt;img src=x onerror=alert(1)&gt;" in page


def test_no_page_carries_script_or_reaches_the_network(
    app: ReviewApp, repository: JsonCaseRepository
):
    registered(repository, events=[event()])

    for page in (app.queue(), app.case("NICU-014")):
        assert "<script" not in page.lower()
        assert "http://" not in page and "https://" not in page


# --------------------------------------------------------------------------
# Over a real socket
# --------------------------------------------------------------------------


@pytest.fixture
def server(repository: JsonCaseRepository):
    """The handler, on an ephemeral port, so the routing is exercised too.

    Worth the fixture: the routes are string surgery on a path, and a unit test
    of ``ReviewApp`` would pass with every one of them mis-wired.
    """
    holder: dict[str, ThreadingHTTPServer] = {}

    def run() -> None:
        original = ThreadingHTTPServer.serve_forever

        def capture(self: ThreadingHTTPServer, *args: object, **kwargs: object) -> None:
            holder["server"] = self
            started.set()
            original(self)

        ThreadingHTTPServer.serve_forever = capture  # type: ignore[method-assign]
        try:
            webapp.serve("127.0.0.1:0", repository, logging.getLogger("test"))
        finally:
            ThreadingHTTPServer.serve_forever = original  # type: ignore[method-assign]

    started = threading.Event()
    thread = threading.Thread(target=run, daemon=True)
    thread.start()
    assert started.wait(5), "server did not start"
    yield f"http://127.0.0.1:{holder['server'].server_address[1]}"
    holder["server"].shutdown()


def post(base: str, path: str, fields: dict[str, str]) -> int:
    body = urllib.parse.urlencode(fields).encode()
    request = urllib.request.Request(base + path, data=body, method="POST")
    try:
        with urllib.request.urlopen(request) as response:
            return int(response.status)
    except urllib.error.HTTPError as failure:
        return int(failure.code)


def test_the_routes_are_wired(server: str, repository: JsonCaseRepository):
    registered(repository, events=[event()])

    with urllib.request.urlopen(server + "/") as response:
        assert response.status == 200
    with urllib.request.urlopen(server + "/case/NICU-014") as response:
        assert b"PLUSG-c158" in response.read()


def test_a_successful_post_redirects_so_a_refresh_does_not_replay_it(
    server: str, repository: JsonCaseRepository
):
    """Redirect after post, so pressing F5 does not try to sign the change
    off twice — the second attempt would be refused, but with an error page
    that reads as though something went wrong."""
    registered(repository, events=[event()])
    host, _, port = server.removeprefix("http://").partition(":")

    # http.client rather than urllib: urllib follows a 303 and hands back the
    # 200 from the page after it, which hides the thing under test.
    connection = http.client.HTTPConnection(host, int(port))
    connection.request(
        "POST",
        "/case/NICU-014/acknowledge",
        body=urllib.parse.urlencode({"event_id": "evt-1", "by": "B. Kowalski"}),
        headers={"Content-Type": "application/x-www-form-urlencoded"},
    )
    response = connection.getresponse()

    assert response.status == 303
    # The banner rides along in the query string, because a redirect that lands
    # on an identical-looking page cannot say whether anything happened.
    assert response.getheader("Location").startswith("/case/NICU-014?said=")
    assert reload_event(repository).status is EventStatus.ACKNOWLEDGED


def test_a_blank_reviewer_over_http_is_a_refusal_not_a_crash(
    server: str, repository: JsonCaseRepository
):
    registered(repository, events=[event()])

    assert post(server, "/case/NICU-014/acknowledge", {"event_id": "evt-1", "by": ""}) == 400
    assert reload_event(repository).status is EventStatus.OPEN


def test_acknowledging_something_that_is_not_open_is_a_conflict(
    server: str, repository: JsonCaseRepository
):
    registered(repository, events=[event()])
    fields = {"event_id": "evt-nope", "by": "B. Kowalski"}

    assert post(server, "/case/NICU-014/acknowledge", fields) == 409


def test_an_unknown_path_is_a_404_rather_than_a_stack_trace(server: str):
    request = urllib.request.Request(server + "/admin")
    with pytest.raises(urllib.error.HTTPError) as refusal:
        urllib.request.urlopen(request)

    assert refusal.value.code == 404


# --------------------------------------------------------------------------
# Registering and running, from the browser
# --------------------------------------------------------------------------

DEPLOYMENT = Path(__file__).parents[1] / "data" / "deployment"


@pytest.fixture
def operating(repository: JsonCaseRepository) -> ReviewApp:
    """A server that can start runs, because it was given a catalogue."""
    return ReviewApp(repository, catalogue=DEPLOYMENT / "catalogue.yaml")


def test_a_case_can_be_registered_from_the_form(operating: ReviewApp):
    case_id = operating.register(
        "NICU-021", "PLUSG:stop_gained:158:heterozygous", "hpo-day-1", ""
    )

    assert case_id == "NICU-021"
    assert operating.run("NICU-021") is not None


def test_several_variants_come_from_several_lines(operating: ReviewApp, tmp_path: Path):
    operating.register(
        "NICU-021",
        "PLUSG:stop_gained:158\n# a comment, ignored\nDEMOG:missense_variant:306",
        "",
        "",
    )

    ledger = JsonCaseRepository(tmp_path).get("NICU-021")
    assert ledger.variant_keys == ("PLUSG-c158", "DEMOG-c306")


def test_an_unknown_consequence_is_named_rather_than_rejected_generically(
    operating: ReviewApp,
):
    """The person typed something into a box. "Bad request" tells them nothing;
    the list of what the rules understand tells them what to type instead."""
    with pytest.raises(OperationError, match="unknown consequence"):
        operating.register("X-1", "PLUSG:nonsens:158", "", "")


def test_a_case_with_no_variants_is_refused(operating: ReviewApp):
    with pytest.raises(OperationError, match="watching nothing"):
        operating.register("X-1", "   \n\n", "", "")


def test_registering_the_same_case_twice_is_refused(operating: ReviewApp):
    """Re-registering would replace the watched variants while keeping the
    history, so the case could carry events about variants it no longer
    watches — a system that looks like it lost track."""
    operating.register("NICU-021", "PLUSG:stop_gained:158", "", "")

    with pytest.raises(OperationError, match="already registered"):
        operating.register("NICU-021", "DEMOG:missense_variant:306", "", "")


def test_a_case_identifier_cannot_escape_the_state_directory(operating: ReviewApp):
    """The identifier becomes a file name, and this form faces a browser."""
    with pytest.raises(OperationError, match="cannot contain a path"):
        operating.register("../../etc/passwd", "PLUSG:stop_gained:158", "", "")


def test_a_run_records_that_the_case_was_examined(operating: ReviewApp, tmp_path: Path):
    """Even when it finds nothing — that is the whole point of the field."""
    operating.register("NICU-021", "PLUSG:stop_gained:158", "", "")

    operating.run("NICU-021")

    assert JsonCaseRepository(tmp_path).get("NICU-021").last_examined_at is not None


def test_a_server_without_a_catalogue_says_so_instead_of_hiding_the_button(
    app: ReviewApp, repository: JsonCaseRepository
):
    """Absence of a catalogue is a state to report, not an error to raise: such
    a server can still show and sign off a queue that cron fills."""
    registered(repository, events=[event()])

    assert app.can_run is False
    assert "cannot run a comparison" in app.queue()
    with pytest.raises(OperationError, match="without --catalogue"):
        app.run("NICU-014")


def test_running_every_case_reports_the_ones_it_could_not_run(
    operating: ReviewApp, repository: JsonCaseRepository
):
    """One case with a missing variant file must not stop the other forty, and
    a caller told "1 examined" without being told which failed cannot act."""
    operating.register("GOOD-1", "PLUSG:stop_gained:158", "", "")
    registered(repository, "BROKEN-1")  # a ledger with no variants sidecar

    examined, failures = operating.run_all()

    assert examined == 1
    assert "BROKEN-1" in failures
    assert "variant file" in failures["BROKEN-1"]


def test_the_queue_offers_both_actions_when_it_can_run(operating: ReviewApp):
    operating.register("NICU-021", "PLUSG:stop_gained:158", "", "")

    page = operating.queue()

    assert 'action="/run-all"' in page
    assert 'action="/register"' in page


def test_an_action_at_an_unknown_address_is_a_lookup_failure(operating: ReviewApp):
    """So the handler answers 404 rather than inventing a redirect."""
    with pytest.raises(LookupError):
        operating.perform("/wipe-everything", {})


def test_the_banner_reports_what_the_last_action_did(operating: ReviewApp):
    """Every action redirects, and a page that comes back identical cannot say
    whether anything happened. "Nothing moved" is the commonest result."""
    operating.register("NICU-021", "PLUSG:stop_gained:158", "", "")
    _, said = operating.perform("/case/NICU-021/run", {})

    assert "re-examined" in said
    assert escape(said) in operating.queue(said)


# --------------------------------------------------------------------------
# The baseline, and why registration takes it
# --------------------------------------------------------------------------


def test_registering_lays_down_the_baseline(operating: ReviewApp, tmp_path: Path):
    """Without this the first run after registration can never report drift,
    however much moved in between: it would be seeing every variant for the
    first time and would quietly adopt today's answer as the reference."""
    operating.register("NICU-021", "PLUSG:stop_gained:158", "", "")

    stored = JsonAssessmentStore(tmp_path).latest_for("NICU-021", "PLUSG-c158")

    assert stored is not None
    assert JsonCaseRepository(tmp_path).get("NICU-021").last_world is not None


def test_a_run_right_after_registering_against_newer_releases_reports_the_drift(
    tmp_path: Path,
):
    """The whole loop, in the two steps a browser actually performs. Registered
    when ClinGen still supported dosage sensitivity for the gene; run once the
    refutation has landed."""
    repository = JsonCaseRepository(tmp_path)
    old = ReviewApp(repository, catalogue=DEPLOYMENT / "catalogue-old.yaml")
    old.register("NICU-021", "PLUSG:stop_gained:158:heterozygous", "hpo-day-1", "")

    current = ReviewApp(repository, catalogue=DEPLOYMENT / "catalogue.yaml")
    report = current.run("NICU-021")

    assert report.events  # type: ignore[attr-defined]
    raised = report.events[0]  # type: ignore[attr-defined]
    assert "gene_curation" in raised.summary()
    assert raised.urgency.at_least(Urgency.HIGH)


def test_registering_without_a_catalogue_still_registers(
    app: ReviewApp, tmp_path: Path
):
    """A server that cannot run is still allowed to record who is being watched;
    the baseline is simply taken by whatever runs next."""
    app.register("NICU-021", "PLUSG:stop_gained:158", "", "")

    assert "NICU-021" in JsonCaseRepository(tmp_path).case_ids()


def test_a_baseline_that_could_not_be_taken_does_not_lose_the_registration(
    tmp_path: Path,
):
    """The variants are the part somebody typed. Discarding them because a
    release file was unreadable would make the person enter them again to fix
    something that is not their problem."""
    broken = tmp_path / "not-a-catalogue.yaml"
    broken.write_text("gene_curation: [{version: '1', path: nowhere.tsv}]\n")

    registration = operations.register(
        tmp_path, "NICU-021", [operations.parse_variant("PLUSG:stop_gained:158")],
        catalogue_path=broken,
    )

    assert registration.baseline.startswith("no baseline taken")
    assert "NICU-021" in JsonCaseRepository(tmp_path).case_ids()


def test_seeding_the_demo_is_idempotent(tmp_path: Path):
    """Restarting `repairbench demo` must not stack duplicate cases, and must
    not wipe the acknowledgements somebody made while trying it."""
    older, log = DEPLOYMENT / "catalogue-old.yaml", logging.getLogger("t")
    case_id = operations.seed_demo(tmp_path, older, log)
    again = operations.seed_demo(tmp_path, older, log)

    assert case_id == again
    assert JsonCaseRepository(tmp_path).case_ids() == [case_id]


def test_the_seeded_case_is_named_so_nobody_mistakes_it_for_a_patient():
    assert operations.DEMO_CASE.startswith("DEMO")
