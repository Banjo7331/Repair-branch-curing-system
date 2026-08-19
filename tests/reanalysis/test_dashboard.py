"""The page a human opens, and the two mistakes it must not make.

A dashboard is where a system's honesty goes to die. It is the one surface
nobody diffs, it summarises by deleting, and every deletion is a judgement made
by whoever wrote the template rather than by the policy that was reviewed.

So these tests are almost all about what the page is forbidden to do:

* it may not change an urgency, in either direction;
* it may not present a case nobody has examined as a case with nothing wrong;
* it may not drop a quiet case, because a page of three rows has to mean
  "three cases exist", not "three cases need attention";
* and it may not interpolate a case identifier that came out of somebody
  else's file without escaping it.
"""

from __future__ import annotations

from datetime import UTC, datetime, timedelta
from pathlib import Path

import pytest

from repairbench.reanalysis.dashboard import (
    STALE_AFTER,
    collect,
    render_html,
    render_text,
    write,
)
from repairbench.reanalysis.ledger import EventStatus
from repairbench.reanalysis.routing import ReviewQueue, Urgency
from repairbench.reanalysis.store import JsonCaseRepository, StoredEvent
from repairbench.reanalysis.world import DriftAxis, Pin, World

NOW = datetime(2026, 8, 17, 12, 0, tzinfo=UTC)


def phenotype() -> Pin:
    return Pin(axis=DriftAxis.PHENOTYPE, version="hpo-day-1", digest="0" * 64)


def world() -> World:
    return World.of(
        Pin(axis=axis, version="v1", digest=f"{index}" * 64)
        for index, axis in enumerate(DriftAxis)
    )


def event(
    event_id: str = "e1",
    *,
    urgency: Urgency = Urgency.HIGH,
    queue: ReviewQueue = ReviewQueue.CLINICAL_SIGNOUT,
    status: EventStatus = EventStatus.OPEN,
    variant_key: str = "PLUSG-c158",
) -> StoredEvent:
    return StoredEvent(
        event_id=event_id,
        variant_key=variant_key,
        fingerprint=f"fp-{event_id}",
        status=status,
        summary="loss_of_function → undetermined: the mechanism no longer resolves",
        urgency=urgency,
        queue=queue,
    )


@pytest.fixture
def repository(tmp_path: Path) -> JsonCaseRepository:
    return JsonCaseRepository(tmp_path)


def registered(
    repository: JsonCaseRepository,
    case_id: str,
    *,
    events: list[StoredEvent] | None = None,
    examined: datetime | None = NOW,
) -> None:
    ledger = repository.register(case_id, ["PLUSG-c158"], phenotype())
    ledger.events = list(events or [])
    ledger.last_world = world()
    ledger.last_examined_at = examined
    repository.save(ledger)


# --------------------------------------------------------------------------
# What it reads
# --------------------------------------------------------------------------


def test_an_empty_state_directory_is_a_page_that_says_nothing_is_watched(
    repository: JsonCaseRepository,
):
    view = collect(repository, now=NOW)

    assert view.rows == ()
    assert "Nothing is being watched." in render_html(view).replace("\n", " ")


def test_a_quiet_case_still_gets_a_row(repository: JsonCaseRepository):
    """Three rows must mean three cases exist. Hiding the quiet ones would make
    the same page mean "three cases need attention" on a different day."""
    registered(repository, "QUIET-1")

    view = collect(repository, now=NOW)

    assert [row.case_id for row in view.rows] == ["QUIET-1"]
    assert view.waiting == ()
    assert "nothing open" in render_html(view)


def test_sidecar_files_beside_a_ledger_are_not_read_as_cases(
    repository: JsonCaseRepository, tmp_path: Path
):
    """The command line writes `<case>.variants.json` into the same directory,
    and the first thing this page ever did was crash on it."""
    registered(repository, "NICU-014")
    (tmp_path / "cases" / "NICU-014.variants.json").write_text('[{"key": "PLUSG-c158"}]')

    assert [row.case_id for row in collect(repository, now=NOW).rows] == ["NICU-014"]


def test_an_acknowledged_event_leaves_the_queue(repository: JsonCaseRepository):
    registered(repository, "NICU-014", events=[event(status=EventStatus.ACKNOWLEDGED)])

    assert collect(repository, now=NOW).waiting == ()


def test_a_superseded_event_leaves_the_queue(repository: JsonCaseRepository):
    """Two unread alerts about one variant is one alert too many: the older
    describes a world that no longer exists."""
    registered(repository, "NICU-014", events=[event(status=EventStatus.SUPERSEDED)])

    assert collect(repository, now=NOW).waiting == ()


# --------------------------------------------------------------------------
# What it must not decide for itself
# --------------------------------------------------------------------------


def test_the_urgency_shown_is_the_urgency_that_was_recorded(repository: JsonCaseRepository):
    """The page renders a decision somebody's policy made and that was reviewed.
    A dashboard that scored events itself would be a second opinion nobody saw."""
    registered(repository, "NICU-014", events=[event(urgency=Urgency.LOW)])

    view = collect(repository, now=NOW)
    # Split off the stylesheet, which names every urgency class by definition.
    body = render_html(view).split("</style>", 1)[1]

    assert view.rows[0].worst is Urgency.LOW
    assert "u-low" in body
    assert "u-critical" not in body
    assert "u-high" not in body


def test_events_are_filed_under_the_queue_the_run_chose(repository: JsonCaseRepository):
    registered(
        repository,
        "NICU-014",
        events=[
            event("e1", queue=ReviewQueue.VALIDATION),
            event("e2", queue=ReviewQueue.CLINICAL_SIGNOUT, variant_key="DEMOG-c306"),
        ],
    )

    view = collect(repository, now=NOW)

    assert len(view.in_queue(ReviewQueue.VALIDATION)) == 1
    assert len(view.in_queue(ReviewQueue.CLINICAL_SIGNOUT)) == 1


def test_the_most_pressing_change_is_first(repository: JsonCaseRepository):
    registered(
        repository,
        "NICU-014",
        events=[
            event("e1", urgency=Urgency.ROUTINE),
            event("e2", urgency=Urgency.CRITICAL, variant_key="DEMOG-c306"),
        ],
    )

    urgencies = [entry.urgency for _, entry in collect(repository, now=NOW).waiting]

    assert urgencies == [Urgency.CRITICAL, Urgency.ROUTINE]


# --------------------------------------------------------------------------
# Silence, which is the failure mode
# --------------------------------------------------------------------------


def test_a_case_examined_recently_with_nothing_open_is_not_stale(
    repository: JsonCaseRepository,
):
    registered(repository, "QUIET-1", examined=NOW - timedelta(hours=6))

    assert collect(repository, now=NOW).stale == ()


def test_a_case_nobody_has_examined_is_stale_even_with_an_empty_queue(
    repository: JsonCaseRepository,
):
    """The distinction the whole page exists for. A case examined nightly for a
    year with nothing to report and a case whose scheduler died in March both
    show an empty queue, and only one of them is fine.

    The first version of this module measured the age of the last *event*, which
    made every healthy quiet case look dead — the reason `last_examined_at` is
    recorded on every run, including the ones that find nothing."""
    registered(repository, "STALE-1", examined=NOW - STALE_AFTER - timedelta(hours=1))

    view = collect(repository, now=NOW)

    assert view.waiting == ()
    assert [row.case_id for row in view.stale] == ["STALE-1"]


def test_a_case_registered_and_never_examined_is_stale(repository: JsonCaseRepository):
    registered(repository, "NEW-1", examined=None)

    assert [row.case_id for row in collect(repository, now=NOW).stale] == ["NEW-1"]
    assert "never" in render_html(collect(repository, now=NOW))


def test_a_healthy_page_says_an_empty_queue_means_nothing_moved(
    repository: JsonCaseRepository,
):
    """Without this sentence an empty page is ambiguous, and the ambiguity
    always resolves in the reassuring direction."""
    registered(repository, "QUIET-1")

    assert "means nothing moved, rather than nothing ran" in render_html(
        collect(repository, now=NOW)
    ).replace("\n", " ")


# --------------------------------------------------------------------------
# The page as a file
# --------------------------------------------------------------------------


def test_the_page_is_self_contained(repository: JsonCaseRepository):
    """No script, no network. A dashboard that needs a server running is a
    dashboard that is down exactly when the pipeline is."""
    registered(repository, "NICU-014", events=[event()])

    page = render_html(collect(repository, now=NOW))

    assert "<script" not in page.lower()
    assert "http://" not in page and "https://" not in page


def test_a_case_identifier_from_a_file_cannot_inject_markup(
    repository: JsonCaseRepository,
):
    """Case identifiers come from whoever ran `watch`, which is a script reading
    somebody's spreadsheet."""
    # No slash in the payload: it also has to be a legal file name, which is
    # exactly the constraint that makes people assume such a value is safe.
    registered(repository, "<img src=x onerror=alert(1)>")

    page = render_html(collect(repository, now=NOW)).split("</style>", 1)[1]

    assert "<img src=x" not in page
    assert "&lt;img src=x onerror=alert(1)&gt;" in page


def test_the_page_names_the_state_directory_it_read(repository: JsonCaseRepository):
    registered(repository, "NICU-014")

    assert str(repository.root) in render_html(collect(repository, now=NOW))


def test_writing_the_page_returns_the_path_it_wrote(
    repository: JsonCaseRepository, tmp_path: Path
):
    registered(repository, "NICU-014")

    written = write(collect(repository, now=NOW), tmp_path / "queue.html")

    assert written.exists()
    assert written.read_text().startswith("<!doctype html>")


def test_the_text_form_leads_with_the_thing_that_needs_a_person(
    repository: JsonCaseRepository,
):
    """This is what a cron job mails. A digest whose first line is not the
    headline is a digest nobody reads to the end."""
    registered(repository, "NICU-014", events=[event()])

    first = render_text(collect(repository, now=NOW)).splitlines()[0]

    assert "1 case(s) watched" in first
    assert "1 change(s) waiting" in first


def test_the_text_form_says_when_nothing_is_waiting(repository: JsonCaseRepository):
    registered(repository, "QUIET-1")

    assert "nothing waiting" in render_text(collect(repository, now=NOW))
