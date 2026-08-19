"""The queue, as a page somebody actually opens.

Everything else in this package renders one answer for one question asked now.
Reanalysis is not that: it runs at three in the morning, compares today with
last month, and exits. Nobody is watching it, which is the whole point — and it
means the output has to survive being read a week late by someone who was not
there when it ran.

Two audiences, and the page serves both because they fail together.

**The reviewer** needs the work list: which cases moved, what the change was,
how loudly it is asking, and which queue it went to. That is the ledger, filed
by urgency, and this module does not compute one number that is not already in
it. A dashboard that derived its own severity would be a second opinion nobody
reviewed.

**The operator** needs the thing no error rate detects: the job quietly not
running. The failure mode of a scheduled reanalysis is not a stack trace, it is
silence that looks exactly like *nothing changed*. So the loudest element on the
page is not an event — it is a case whose last run is old. An empty queue under
a fresh run is the system working; an empty queue under a run from March is the
system dead, and the two must not look alike.

Three refusals.

**No JavaScript, no network, no build step.** One file, inline style, opens from
disk. A dashboard that needs a server running is a dashboard that is down
exactly when the pipeline is.

**Nothing is recomputed.** Urgency and queue are read from the ledger as the run
recorded them. The page cannot promote an event, and cannot quietly demote one
either.

**Silence is rendered, not omitted.** A case with no open events gets a row
saying so. Leaving it out would make a page of five rows mean either "five cases
need attention" or "five cases exist", and a reviewer cannot tell which.
"""

from __future__ import annotations

import html
from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from pathlib import Path

from repairbench.reanalysis.ledger import CaseLedger, LedgerEntry
from repairbench.reanalysis.routing import ReviewQueue, Urgency
from repairbench.reanalysis.store import JsonCaseRepository
from repairbench.reanalysis.world import World

#: Past this, a case is reported as stale rather than quiet. Two days, matching
#: the alert the metrics module documents — the two should not disagree about
#: what "overdue" means, because they are read by the same person.
STALE_AFTER = timedelta(days=2)


@dataclass(frozen=True, slots=True)
class CaseRow:
    """One watched case, as the page shows it."""

    case_id: str
    variants: int
    open_events: tuple[LedgerEntry, ...]
    world: World | None
    #: When a run last examined this case. Not when it last *found* something —
    #: the first version of this page measured the second, and every healthy
    #: quiet case came out looking dead.
    last_examined_at: datetime | None

    @property
    def worst(self) -> Urgency | None:
        """The most pressing thing waiting, or ``None`` when nothing is."""
        urgencies = [event.urgency for event in self.open_events]  # type: ignore[attr-defined]
        return min(urgencies, key=lambda u: u.rank) if urgencies else None

    def age(self, now: datetime) -> timedelta | None:
        return None if self.last_examined_at is None else now - self.last_examined_at

    def is_stale(self, now: datetime) -> bool:
        """Has this case gone quiet for longer than a run interval?

        ``None`` — registered by ``watch`` and never re-examined since — counts
        as stale. It is the state a case is in for the hour between registration
        and the first scheduled run, and also the state it is in forever if that
        run never happens.
        """
        age = self.age(now)
        return age is None or age > STALE_AFTER


@dataclass(frozen=True, slots=True)
class QueueView:
    """Every watched case, and what is waiting across all of them."""

    rows: tuple[CaseRow, ...]
    generated_at: datetime
    state_directory: str

    @property
    def waiting(self) -> tuple[tuple[CaseRow, LedgerEntry], ...]:
        """Open events across every case, most pressing first.

        Paired with the case rather than flattened, because "which patient" is
        the first thing a reviewer needs and the event does not carry it.
        """
        pairs = [(row, event) for row in self.rows for event in row.open_events]
        pairs.sort(key=lambda pair: (pair[1].urgency.rank, pair[0].case_id))  # type: ignore[attr-defined]
        return tuple(pairs)

    def in_queue(self, queue: ReviewQueue) -> tuple[tuple[CaseRow, LedgerEntry], ...]:
        return tuple(pair for pair in self.waiting if pair[1].queue is queue)  # type: ignore[attr-defined]

    @property
    def stale(self) -> tuple[CaseRow, ...]:
        return tuple(row for row in self.rows if row.is_stale(self.generated_at))


def collect(repository: JsonCaseRepository, *, now: datetime | None = None) -> QueueView:
    """Read every case in the state directory into one view."""
    moment = now or datetime.now(UTC)
    rows = []
    for case_id in sorted(repository.case_ids()):
        ledger = repository.get(case_id)
        rows.append(_row(ledger))
    return QueueView(
        rows=tuple(rows),
        generated_at=moment,
        state_directory=str(repository.root),
    )


def _row(ledger: CaseLedger) -> CaseRow:
    return CaseRow(
        case_id=ledger.case_id,
        variants=len(ledger.variant_keys),
        open_events=ledger.open_events(),
        world=ledger.last_world,
        last_examined_at=ledger.last_examined_at,
    )


# --------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------

STYLE = """
:root { --ink:#16181d; --muted:#6b7280; --line:#e5e7eb; --paper:#fff; --wash:#f7f8fa;
        --critical:#b42318; --high:#c2410c; --routine:#a16207; --low:#3f6212; --silent:#6b7280; }
* { box-sizing:border-box }
body { margin:0; padding:2.5rem 1.5rem; background:var(--wash); color:var(--ink);
       font:15px/1.55 ui-sans-serif,-apple-system,"Segoe UI",Roboto,sans-serif }
main { max-width:70rem; margin:0 auto }
h1 { font-size:1.4rem; margin:0 0 .2rem }
h2 { font-size:1rem; margin:2.2rem 0 .6rem; text-transform:uppercase;
     letter-spacing:.06em; color:var(--muted) }
.sub { color:var(--muted); margin:0 0 1.6rem; font-size:.9rem }
.cards { display:flex; flex-wrap:wrap; gap:.75rem; margin-bottom:.5rem }
.card { background:var(--paper); border:1px solid var(--line); border-radius:.5rem;
        padding:.8rem 1.1rem; min-width:9rem }
.card .n { font-size:1.7rem; font-weight:600; line-height:1.1 }
.card .k { color:var(--muted); font-size:.82rem }
.card.alarm { border-color:var(--critical); background:#fef3f2 }
.card.alarm .n { color:var(--critical) }
table { width:100%; border-collapse:collapse; background:var(--paper);
        border:1px solid var(--line); border-radius:.5rem; overflow:hidden }
th { text-align:left; font-size:.75rem; text-transform:uppercase; letter-spacing:.05em;
     color:var(--muted); padding:.6rem .9rem; border-bottom:1px solid var(--line) }
td { padding:.7rem .9rem; border-bottom:1px solid var(--line); vertical-align:top }
tr:last-child td { border-bottom:none }
.u { font-weight:600; font-size:.78rem; text-transform:uppercase; letter-spacing:.04em }
.u-critical{color:var(--critical)} .u-high{color:var(--high)} .u-routine{color:var(--routine)}
.u-low{color:var(--low)} .u-silent{color:var(--silent)}
code { font:13px/1.4 ui-monospace,SFMono-Regular,Menlo,monospace; background:var(--wash);
       padding:.1rem .3rem; border-radius:.2rem; white-space:nowrap }
.reason { color:var(--muted); font-size:.88rem; margin-top:.25rem }
.quiet { color:var(--muted) }
.stale { color:var(--critical); font-weight:600 }
.note { background:var(--paper); border:1px solid var(--line); border-left:3px solid var(--muted);
        border-radius:.4rem; padding:.9rem 1.1rem; color:var(--muted); font-size:.88rem }
.pins { font:12px/1.5 ui-monospace,SFMono-Regular,Menlo,monospace; color:var(--muted) }
footer { margin-top:2.5rem; color:var(--muted); font-size:.82rem;
         border-top:1px solid var(--line); padding-top:1rem }
"""


def render_html(view: QueueView, *, version: str = "") -> str:
    """The whole page, as one self-contained document."""
    waiting = view.waiting
    clinical = view.in_queue(ReviewQueue.CLINICAL_SIGNOUT)
    validation = view.in_queue(ReviewQueue.VALIDATION)
    stale = view.stale

    parts = [
        "<!doctype html>",
        '<html lang="en"><head><meta charset="utf-8">',
        '<meta name="viewport" content="width=device-width,initial-scale=1">',
        "<title>repairbench — reanalysis queue</title>",
        f"<style>{STYLE}</style></head><body><main>",
        "<h1>Reanalysis queue</h1>",
        f'<p class="sub">{escape(view.state_directory)} · generated '
        f"{escape(view.generated_at.strftime('%Y-%m-%d %H:%M UTC'))}"
        f"{' · ' + escape(version) if version else ''}</p>",
        _summary_cards(view, waiting, clinical, validation, stale),
        _staleness_section(view, stale),
        _queue_section("Waiting for clinical sign-out", clinical),
        _queue_section("Waiting for validation", validation),
        _cases_section(view),
        _footer(),
        "</main></body></html>",
    ]
    return "\n".join(parts)


def _summary_cards(
    view: QueueView,
    waiting: tuple[tuple[CaseRow, LedgerEntry], ...],
    clinical: tuple[tuple[CaseRow, LedgerEntry], ...],
    validation: tuple[tuple[CaseRow, LedgerEntry], ...],
    stale: tuple[CaseRow, ...],
) -> str:
    """Four numbers, and the stale one is styled to be the one you see first."""
    cards = [
        ("cases watched", len(view.rows), False),
        ("open changes", len(waiting), False),
        ("clinical sign-out", len(clinical), False),
        ("validation", len(validation), False),
        ("not run recently", len(stale), bool(stale)),
    ]
    rendered = "".join(
        f'<div class="card{" alarm" if alarm else ""}">'
        f'<div class="n">{count}</div><div class="k">{escape(label)}</div></div>'
        for label, count, alarm in cards
    )
    return f'<div class="cards">{rendered}</div>'


def _staleness_section(view: QueueView, stale: tuple[CaseRow, ...]) -> str:
    """The section that exists because silence is the failure mode.

    A queue with nothing in it is the commonest correct output this system
    produces, and it is indistinguishable from a scheduler that stopped. This
    says which of the two it is, before anything else on the page.
    """
    if not stale:
        return (
            '<h2>Scheduler</h2><div class="note">Every watched case has been '
            f"re-examined within the last {STALE_AFTER.days} days. An empty queue below "
            "therefore means nothing moved, rather than nothing ran.</div>"
        )
    rows = "".join(
        f"<tr><td><code>{escape(row.case_id)}</code></td>"
        f'<td class="stale">{_age(row.age(view.generated_at))}</td>'
        f'<td class="quiet">{row.variants} variant(s) watched</td></tr>'
        for row in stale
    )
    return (
        "<h2>Not run recently</h2>"
        '<div class="note">A scheduled reanalysis fails by not happening. Nothing has '
        f"examined these cases for longer than {STALE_AFTER.days} days, so an empty queue "
        "for them is not evidence that nothing moved.</div>"
        f"<table><tr><th>case</th><th>last examined</th><th></th></tr>{rows}</table>"
    )


def _queue_section(title: str, pairs: tuple[tuple[CaseRow, LedgerEntry], ...]) -> str:
    if not pairs:
        return f'<h2>{escape(title)}</h2><div class="note">Nothing waiting.</div>'
    rows = "".join(
        f"<tr><td><code>{escape(row.case_id)}</code></td>"
        f"<td>{_urgency_badge(event)}</td>"
        f"<td><code>{escape(event.variant_key)}</code>"
        f'<div class="reason">{escape(summary_of(event))}</div></td></tr>'
        for row, event in pairs
    )
    return (
        f"<h2>{escape(title)}</h2>"
        f"<table><tr><th>case</th><th>urgency</th><th>change</th></tr>{rows}</table>"
    )


def _urgency_badge(event: LedgerEntry) -> str:
    urgency = escape(str(event.urgency))  # type: ignore[attr-defined]
    return f'<span class="u u-{urgency}">{urgency}</span>'


def _cases_section(view: QueueView) -> str:
    """Every case, including the quiet ones — see the module docstring."""
    if not view.rows:
        return (
            '<h2>Cases</h2><div class="note">No cases are registered in this state '
            "directory. Nothing is being watched.</div>"
        )
    rows = "".join(
        f"<tr><td><code>{escape(row.case_id)}</code></td>"
        f"<td>{row.variants}</td>"
        f"<td>{_worst_cell(row)}</td>"
        f'<td class="pins">{escape(describe_world(row.world))}</td></tr>'
        for row in view.rows
    )
    return (
        "<h2>Cases</h2>"
        "<table><tr><th>case</th><th>variants</th><th>waiting</th>"
        f"<th>world last compared against</th></tr>{rows}</table>"
    )


def _worst_cell(row: CaseRow) -> str:
    worst = row.worst
    if worst is None:
        return '<span class="quiet">nothing open</span>'
    return (
        f'<span class="u u-{escape(str(worst))}">{escape(str(worst))}</span> '
        f'<span class="quiet">({len(row.open_events)})</span>'
    )


def _footer() -> str:
    return (
        "<footer>Every urgency and queue on this page was decided by the run that raised "
        "the event and is reproduced here unchanged — this page ranks nothing and promotes "
        "nothing. Not a medical device: a change reaching clinical sign-out is a prompt to "
        "look, not a finding.</footer>"
    )


def summary_of(event: LedgerEntry) -> str:
    """The sentence a reviewer reads, from either kind of ledger entry."""
    stored = getattr(event, "summary", None)
    if callable(stored):
        return str(stored())
    return str(stored) if stored else str(event.event_id)


def describe_world(world: World | None) -> str:
    if world is None:
        return "not recorded"
    return " ".join(f"{pin.axis.value}@{pin.version}" for pin in world.pins)


def _age(age: timedelta | None) -> str:
    if age is None:
        return "never"
    days = age.days
    if days >= 1:
        return f"{days} day(s) ago"
    return f"{age.seconds // 3600} hour(s) ago"


def escape(text: object) -> str:
    """Escape. Case identifiers and gene symbols come from files somebody else
    wrote, and a page that interpolates them raw is a page that can be made to
    say anything.

    Public, along with ``STYLE``, ``describe_world`` and ``summary_of``, because
    the review server renders pages that have to look and behave like these. Two
    copies of an escaping helper is how one of them ends up not being called.
    """
    return html.escape(str(text), quote=True)


def write(view: QueueView, path: Path, *, version: str = "") -> Path:
    path = Path(path)
    path.write_text(render_html(view, version=version), encoding="utf-8")
    return path


def render_text(view: QueueView) -> str:
    """The same view for a terminal, because a page is not always what is wanted.

    Kept deliberately short: this is the thing a cron job can mail, and a
    fifty-line digest is a digest nobody reads to the end.
    """
    lines = [f"{len(view.rows)} case(s) watched, {len(view.waiting)} change(s) waiting"]
    if view.stale:
        lines.append(
            f"  ! {len(view.stale)} case(s) not examined for over {STALE_AFTER.days} "
            "days — an empty queue for them means nothing"
        )
    for row, event in view.waiting:
        lines.append(f"  [{event.urgency}] {row.case_id}  {event.variant_key}")  # type: ignore[attr-defined]
        lines.append(f"      {summary_of(event)}")
    if not view.waiting:
        lines.append("  nothing waiting")
    return "\n".join(lines)
