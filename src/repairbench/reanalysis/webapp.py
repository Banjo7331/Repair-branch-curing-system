"""A local review server: the queue, one case at a time, and one button.

The dashboard writes a page. This serves one, and the difference is not
convenience — it is that a page cannot *do* anything, and there is exactly one
thing a reviewer needs to do that nothing in this package could until now.

**It is the operating surface for the reanalysis half.** Register a case, run a
comparison, read what moved, sign it off. What it deliberately cannot do is edit
a conclusion: ``run`` performs exactly the comparison a scheduled process would,
from the same pinned files, and nothing here lets a caller assert a mechanism,
an urgency or a queue. A button that started a run is a person doing what the
clock does; a button that changed a verdict would be a different product.

**Acknowledging is the gap this closes.** The surfacing policy suppresses any
change whose fingerprint has already been acknowledged, so acknowledgement is
what stops a queue from being a list that only grows. ``CaseLedger`` has had the
method since the beginning and nothing exposed it, which meant a reviewer who
read a change, decided it needed nothing, and closed the tab would be shown the
same change at the next release, forever.

Three decisions, and the first two are the reason this is not a web framework.

**Standard library only.** ``http.server``, the same thing ``observability.py``
already runs. A package whose whole argument is that its clinical judgements
live in readable files should not acquire a framework, a template engine, a
migration tool and an ORM in order to render two tables and accept one form.

**Server-rendered, no client state.** Every page is a full document; every
action is a form post followed by a redirect. Nothing here depends on
JavaScript, so the review surface degrades to *readable* rather than to *blank*.

**It binds to the loopback address and says what that is worth.** This is a
single-operator tool for a laboratory workstation. The name typed into the
acknowledgement box is **attribution, not authentication** — it records who
*said* they reviewed a change, which is what the ledger needs and is not the
same as proving it. A deployment where that distinction matters needs an
identity provider in front of this, and this module is not going to pretend
otherwise by adding a password box.

One thing it deliberately will not do: change anything except an event's
acknowledged status. It cannot re-run an analysis, edit a rule, or alter an
urgency. Every other number on these pages was decided by a scheduled run and is
reproduced unchanged, so the worst a mistake here can do is mark one event read.
"""

from __future__ import annotations

import logging
from datetime import UTC, datetime
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, quote, unquote

from repairbench.reanalysis import operations
from repairbench.reanalysis.dashboard import (
    STYLE,
    collect,
    describe_world,
    escape,
    render_html,
    summary_of,
)
from repairbench.reanalysis.ledger import CaseLedger, EventStatus
from repairbench.reanalysis.operations import OperationError
from repairbench.reanalysis.store import JsonCaseRepository, StoreError

#: Loopback by default, and it is a default with an argument behind it. See the
#: module docstring: there is no authentication here, so the only binding this
#: module can honestly offer out of the box is one that does not leave the
#: machine.
DEFAULT_ADDRESS = "127.0.0.1:8080"


class ReviewApp:
    """The routes, separated from the socket so they can be tested without one."""

    def __init__(
        self,
        repository: JsonCaseRepository,
        *,
        catalogue: Path | None = None,
        logger: logging.Logger | None = None,
        version: str = "",
    ) -> None:
        self._cases = repository
        #: Where the releases live. Optional, and its absence is a *state* the
        #: pages report rather than an error they raise: a server started
        #: without one can still show and sign off a queue that cron fills, and
        #: saying "no catalogue was given, so nothing here can start a run" is
        #: more use than a missing button.
        self._catalogue = catalogue
        self._logger = logger or logging.getLogger("repairbench")
        self._version = version

    @property
    def can_run(self) -> bool:
        return self._catalogue is not None

    # -- reading ---------------------------------------------------------

    def queue(self, message: str = "") -> str:
        """The same page ``dashboard`` writes, with links and the two actions."""
        view = collect(self._cases)
        page = render_html(view, version=self._version)
        # The dashboard is the authority on what the queue says; this only makes
        # the case identifiers reachable. Rewriting its markup here rather than
        # teaching it about URLs keeps the written page free of links that go
        # nowhere when it is mailed to somebody.
        for row in view.rows:
            plain = f"<code>{escape(row.case_id)}</code>"
            linked = f'<a href="/case/{quote(row.case_id)}">{plain}</a>'
            page = page.replace(plain, linked)
        controls = _run_all_control(self.can_run) + _register_form()
        # The banner goes above the counts, not down with the actions: it says
        # what the press the reader just made did, and a report of an action
        # placed below the thing it changed is read after the fact it explains.
        page = page.replace('<div class="cards">', _banner(message) + '<div class="cards">', 1)
        # The written page carries only the dashboard's stylesheet, because it
        # has no forms in it. Served, it does — so the form rules are appended
        # here rather than pushed back into the module that writes files.
        page = page.replace("</style>", _EXTRA_STYLE + "</style>", 1)
        return page.replace("<footer>", controls + "<footer>", 1)

    def case(self, case_id: str, message: str = "") -> str:
        ledger = self._cases.get(case_id)
        return _case_page(ledger, self._version, can_run=self.can_run, message=message)

    # -- the one thing it can change -------------------------------------

    def register(self, case_id: str, variants: str, phenotype: str, tissue: str) -> str:
        """Register a case from the form, and return its identifier.

        The variant text is one specification per line, in the same grammar the
        command line takes. A textarea rather than a repeating widget because
        the format is the one already documented, and a second syntax for the
        same thing is a second thing to get wrong.
        """
        parsed = [
            operations.parse_variant(line)
            for line in variants.splitlines()
            if line.strip() and not line.strip().startswith("#")
        ]
        registration = operations.register(
            self._cases.root,
            case_id,
            parsed,
            phenotype=phenotype.strip() or "unrecorded",
            tissue=tissue.strip(),
            # Registering through the page takes the baseline immediately, so
            # the next run is a comparison rather than a first sighting.
            catalogue_path=self._catalogue,
            logger=self._logger,
        )
        return registration.case_id

    def run(self, case_id: str) -> object:
        """Run one comparison — the same one a scheduled process would."""
        if self._catalogue is None:
            raise OperationError(
                "this server was started without --catalogue, so it cannot run a "
                "comparison. It can still show the queue and sign changes off; the "
                "runs themselves are whatever fills it."
            )
        return operations.run(self._cases.root, self._catalogue, case_id, self._logger)

    def run_all(self) -> tuple[int, dict[str, str]]:
        if self._catalogue is None:
            raise OperationError(
                "this server was started without --catalogue, so it cannot run a "
                "comparison."
            )
        reports, failures = operations.run_all(self._cases.root, self._catalogue, self._logger)
        return len(reports), failures

    def perform(self, path: str, form: dict[str, list[str]]) -> tuple[str, str]:
        """Carry out one posted action, and say where to go and what to report.

        A table rather than a ladder of branches in the request handler, for the
        same reason the rule files are tables: the set of things this server can
        change should be readable in one place. There are four, and three of
        them are here — the fourth, the scheduled run, is not a route at all.

        Raises ``LookupError`` for an address that is not an action, so the
        handler answers 404 rather than inventing a redirect.
        """
        def first(name: str) -> str:
            return form.get(name, [""])[0]

        if path == "/register":
            case_id = self.register(
                first("case_id"), first("variants"), first("phenotype"), first("tissue")
            )
            return f"/case/{quote(case_id)}", f"{case_id} registered"

        if path == "/run-all":
            examined, failures = self.run_all()
            return "/", _run_all_summary(examined, failures)

        if path.startswith("/case/") and path.endswith("/run"):
            case_id = unquote(path[len("/case/") : -len("/run")])
            report = self.run(case_id)
            return f"/case/{quote(case_id)}", report.headline()  # type: ignore[attr-defined]

        if path.startswith("/case/") and path.endswith("/acknowledge"):
            case_id = unquote(path[len("/case/") : -len("/acknowledge")])
            self.acknowledge(
                case_id, first("event_id"), first("by"), first("note")
            )
            return f"/case/{quote(case_id)}", "acknowledged"

        raise LookupError(path)

    def acknowledge(self, case_id: str, event_id: str, by: str, note: str) -> None:
        """Record that a named person read one event.

        Raises rather than reporting success for an event that is not there or
        is not open. A form post that silently did nothing would leave the
        reviewer believing they had signed something off.
        """
        ledger = self._cases.get(case_id)
        if not ledger.acknowledge(event_id, by=by, note=note, at=datetime.now(UTC)):
            raise StoreError(
                f"{event_id} is not an open event on {case_id}. It may have been "
                "acknowledged already, or superseded by a later change to the same "
                "variant — in which case the later one is what needs reading."
            )
        self._cases.save(ledger)


# --------------------------------------------------------------------------
# The case page
# --------------------------------------------------------------------------


def _case_page(
    ledger: CaseLedger, version: str, *, can_run: bool = False, message: str = ""
) -> str:
    open_events = ledger.open_events()
    history = [event for event in ledger.events if event.status is not EventStatus.OPEN]
    return "\n".join(
        [
            "<!doctype html>",
            '<html lang="en"><head><meta charset="utf-8">',
            '<meta name="viewport" content="width=device-width,initial-scale=1">',
            f"<title>{escape(ledger.case_id)} — repairbench</title>",
            f"<style>{STYLE}{_EXTRA_STYLE}</style></head><body><main>",
            '<p class="sub"><a href="/">← queue</a></p>',
            f"<h1>{escape(ledger.case_id)}</h1>",
            f'<p class="sub">{len(ledger.variant_keys)} variant(s) watched · '
            f"phenotype {escape(ledger.phenotype.version)}"
            f"{' · ' + escape(version) if version else ''}</p>",
            _banner(message),
            _run_one_control(ledger.case_id, can_run),
            _waiting_section(ledger, open_events),
            _history_section(history),
            _world_section(ledger),
            _footer(),
            "</main></body></html>",
        ]
    )


def _waiting_section(ledger: CaseLedger, open_events: tuple[Any, ...]) -> str:
    if not open_events:
        return (
            '<h2>Waiting</h2><div class="note">Nothing open for this case. '
            "That means every change raised so far has been read — not that nothing "
            "has changed, and not that the case was examined recently; the queue page "
            "says which.</div>"
        )
    return "<h2>Waiting</h2>" + "".join(
        _event_card(ledger.case_id, event) for event in open_events
    )


def _event_card(case_id: str, event: Any) -> str:
    """One change, and the form that signs it off.

    The warning under the button is not decoration. Acknowledging suppresses
    every future alert carrying this fingerprint, which is a decision with a
    longer reach than "mark as read" suggests, and the person clicking should
    be told so at the moment they click.
    """
    return (
        '<div class="card wide">'
        f'<div class="head"><span class="u u-{escape(str(event.urgency))}">'
        f"{escape(str(event.urgency))}</span> <code>{escape(event.variant_key)}</code>"
        f'<span class="quiet"> → {escape(str(event.queue))}</span></div>'
        f'<div class="reason">{escape(summary_of(event))}</div>'
        f'<form method="post" action="/case/{quote(case_id)}/acknowledge">'
        f'<input type="hidden" name="event_id" value="{escape(event.event_id)}">'
        '<label>Reviewer<input name="by" required autocomplete="name" '
        'placeholder="who is signing this off"></label>'
        '<label>Note<input name="note" placeholder="why it needs nothing, or what was done">'
        "</label>"
        '<button type="submit">Acknowledge</button>'
        '<p class="warn">This suppresses every future alert with the same fingerprint. '
        "The name is recorded as attribution — this server does not authenticate it.</p>"
        "</form></div>"
    )


def _history_section(history: list[Any]) -> str:
    """Everything that is no longer open, and why it is not.

    Kept on the page because the ledger is append-only and the point of an
    append-only ledger is that somebody can read it. An acknowledgement whose
    note nobody can find later is a note that was never worth writing.
    """
    if not history:
        return '<h2>History</h2><div class="note">Nothing has been closed yet.</div>'
    rows = "".join(
        f"<tr><td><code>{escape(event.variant_key)}</code></td>"
        f"<td>{escape(str(event.status))}</td>"
        f"<td>{escape(_closed_by(event))}</td>"
        f'<td class="reason">{escape(summary_of(event))}</td></tr>'
        for event in history
    )
    return (
        "<h2>History</h2><table><tr><th>variant</th><th>status</th>"
        f"<th>by</th><th>change</th></tr>{rows}</table>"
    )


def _closed_by(event: Any) -> str:
    who = getattr(event, "acknowledged_by", "")
    when = getattr(event, "acknowledged_at", None)
    if not who:
        superseded = getattr(event, "superseded_by", None)
        return f"superseded by {superseded}" if superseded else "—"
    note = getattr(event, "acknowledged_note", "")
    stamp = when.strftime("%Y-%m-%d") if when else "date not recorded"
    return f"{who}, {stamp}" + (f" — {note}" if note else "")


def _world_section(ledger: CaseLedger) -> str:
    return (
        "<h2>Last compared against</h2>"
        f'<div class="note pins">{escape(describe_world(ledger.last_world))}</div>'
    )


def _footer() -> str:
    return (
        "<footer>Nothing on this page was computed here. Every urgency, queue and "
        "transition was decided by a scheduled run and is reproduced unchanged; the only "
        "thing this server can alter is whether an event is marked as read, and by whom. "
        "Not a medical device.</footer>"
    )


_EXTRA_STYLE = """
a { color:#1d4ed8 }
.card.wide { min-width:0; margin-bottom:.75rem }
.card .head { margin-bottom:.35rem }
form { margin-top:.9rem; display:flex; flex-wrap:wrap; gap:.6rem; align-items:flex-end }
label { display:flex; flex-direction:column; gap:.2rem; font-size:.78rem;
        text-transform:uppercase; letter-spacing:.04em; color:var(--muted); flex:1 1 14rem }
input { font:inherit; text-transform:none; letter-spacing:normal; color:var(--ink);
        padding:.45rem .6rem; border:1px solid var(--line); border-radius:.35rem }
button { font:inherit; font-weight:600; padding:.5rem 1.1rem; border-radius:.35rem;
         border:1px solid var(--ink); background:var(--ink); color:#fff; cursor:pointer }
.warn { flex:1 1 100%; margin:.1rem 0 0; font-size:.8rem; color:var(--muted) }
form.stack { flex-direction:column; align-items:stretch }
form.stack label { flex:1 1 auto }
form.inline { margin:0 0 1rem }
textarea { font:13px/1.5 ui-monospace,SFMono-Regular,Menlo,monospace; text-transform:none;
           letter-spacing:normal; color:var(--ink); padding:.5rem .6rem;
           border:1px solid var(--line); border-radius:.35rem; resize:vertical }
.note.said { border-left-color:var(--ink); color:var(--ink) }
"""


# --------------------------------------------------------------------------
# The socket
# --------------------------------------------------------------------------


def serve(
    address: str,
    repository: JsonCaseRepository,
    logger: logging.Logger,
    *,
    catalogue: Path | None = None,
    version: str = "",
) -> None:
    """Serve the review pages until interrupted."""
    app = ReviewApp(repository, catalogue=catalogue, logger=logger, version=version)
    host, _, port = address.rpartition(":")

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:
            path, _, query = self.path.partition("?")
            # The one thing carried in a query string: what the previous action
            # did. Redirect-after-post loses it otherwise, and a page that comes
            # back identical is indistinguishable from one where nothing ran.
            said = parse_qs(query).get("said", [""])[0]
            if path == "/":
                self._html(200, app.queue(said))
            elif path.startswith("/case/"):
                self._case(unquote(path[len("/case/") :]), said)
            else:
                self._html(404, _message("Not found", "No page at this address."))

        def do_POST(self) -> None:
            try:
                location, said = app.perform(self.path.split("?", 1)[0], self._form())
            except LookupError:
                self._html(404, _message("Not found", "No action at this address."))
            except OperationError as refusal:
                # A refusal the operations layer made — an unknown consequence, a
                # duplicate case, a server started without a catalogue. Repeated
                # to the person because the reason is the useful part.
                self._html(400, _message("Not done", str(refusal)))
            except ValueError as refusal:
                # A blank reviewer. The refusal is the ledger's, and it is
                # repeated rather than translated into a generic "bad request".
                self._html(400, _message("Not recorded", str(refusal)))
            except StoreError as missing:
                self._html(409, _message("Not recorded", str(missing)))
            else:
                # Redirect after post, so a refresh does not replay the action.
                self._redirect(location, said)

        def _case(self, case_id: str, said: str = "") -> None:
            try:
                self._html(200, app.case(case_id, said))
            except StoreError as missing:
                self._html(404, _message("No such case", str(missing)))

        def _redirect(self, location: str, said: str = "") -> None:
            target = f"{location}?said={quote(said)}" if said else location
            self.send_response(303)
            self.send_header("Location", target)
            self.end_headers()

        def _form(self) -> dict[str, list[str]]:
            length = int(self.headers.get("Content-Length", "0") or 0)
            body = self.rfile.read(length).decode("utf-8", errors="replace")
            return parse_qs(body, keep_blank_values=True)

        def _html(self, status: int, body: str) -> None:
            payload = body.encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(payload)))
            self.end_headers()
            self.wfile.write(payload)

        def log_message(self, fmt: str, *args: Any) -> None:
            logger.info("review %s", fmt % args)

    server = ThreadingHTTPServer((host or "127.0.0.1", int(port)), Handler)
    logger.info("review server on http://%s", address)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


def _run_all_summary(examined: int, failures: dict[str, str]) -> str:
    """What a sweep did, including what it could not do.

    The failures are named rather than counted. "39 of 40 examined" tells a
    reader that something is wrong and not which case to go and look at, which
    is the half of the sentence that would have been useful.
    """
    said = f"{examined} case(s) examined"
    if failures:
        said += "; could not run " + ", ".join(
            f"{case} ({reason})" for case, reason in sorted(failures.items())
        )
    return said


def _message(title: str, detail: str) -> str:
    return (
        "<!doctype html>"
        '<html lang="en"><head><meta charset="utf-8">'
        f"<title>{escape(title)}</title><style>{STYLE}{_EXTRA_STYLE}</style></head>"
        f'<body><main><p class="sub"><a href="/">← queue</a></p>'
        f'<h1>{escape(title)}</h1><div class="note">{escape(detail)}</div>'
        "</main></body></html>"
    )


# --------------------------------------------------------------------------
# The two controls on the queue page
# --------------------------------------------------------------------------


def _banner(message: str) -> str:
    """What the last action did, carried across the redirect in the URL.

    Shown rather than swallowed because every action here is followed by a
    redirect, and a page that comes back looking identical is indistinguishable
    from one where nothing happened. "3 cases examined, nothing moved" is the
    commonest result and the one most worth saying out loud.
    """
    return f'<div class="note said">{escape(message)}</div>' if message else ""


def _run_all_control(can_run: bool) -> str:
    if not can_run:
        return (
            '<h2>Runs</h2><div class="note">This server was started without a release '
            "catalogue, so it cannot run a comparison — it can show the queue and sign "
            "changes off. Restart it with <code>--catalogue</code> to run from here.</div>"
        )
    return (
        "<h2>Runs</h2>"
        '<form method="post" action="/run-all">'
        '<button type="submit">Re-examine every case</button>'
        '<p class="warn">Runs exactly the comparison a scheduled process would, against '
        "the releases the catalogue currently names. It cannot reach a conclusion cron "
        "would not have reached — it is the same run, started by a person.</p></form>"
    )


def _register_form() -> str:
    """One case, one textarea, the grammar the command line already documents."""
    return (
        "<h2>Watch a new case</h2>"
        '<form method="post" action="/register" class="stack">'
        '<label>Case identifier<input name="case_id" required '
        'placeholder="NICU-014"></label>'
        '<label>Variants, one per line<textarea name="variants" rows="3" required '
        'placeholder="PLUSG:stop_gained:158:heterozygous"></textarea></label>'
        '<label>Phenotype snapshot<input name="phenotype" placeholder="unrecorded"></label>'
        '<label>Affected tissue, as GTEx names it<input name="tissue" '
        'placeholder="left blank: nothing is checked against where the gene is on"></label>'
        '<button type="submit">Register</button>'
        '<p class="warn">Each variant is <code>gene:consequence:cds_position[:zygosity]</code>. '
        "Zygosity left out means unknown rather than heterozygous — guessing it would offer "
        "every modality that needs an intact allele.</p>"
        "</form>"
    )


def _run_one_control(case_id: str, can_run: bool) -> str:
    if not can_run:
        return ""
    return (
        f'<form method="post" action="/case/{quote(case_id)}/run" class="inline">'
        '<button type="submit">Re-examine now</button></form>'
    )
