"""The case as a thing that lives in time.

A single assessment is a photograph. A watched case is the film: the same
variants, re-examined at every release, accumulating a record of what changed,
who was told and what they decided. The ledger is append-only — an event is
never edited, only acknowledged or superseded — because "what did we know, and
when" has to survive being asked years later.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from enum import StrEnum
from typing import Protocol

from repairbench.reanalysis.drift import Attribution
from repairbench.reanalysis.routing import ReviewQueue, SurfacingDecision, Urgency
from repairbench.reanalysis.world import Pin, World


class LedgerEntry(Protocol):
    """What the ledger asks of an event, and nothing more.

    Two kinds of thing satisfy it. A ``DriftEvent`` from the run that raised it,
    carrying the full attribution. And a ``StoredEvent`` read back from disk,
    which carries the fingerprint and the status — the two things a later run
    consults — and refuses to pretend it still holds the causal claim.

    Writing the requirement down as a protocol is what lets the second exist
    without the first being weakened into something rehydratable.
    """

    @property
    def event_id(self) -> str: ...

    @property
    def variant_key(self) -> str: ...

    @property
    def fingerprint(self) -> str: ...

    @property
    def status(self) -> EventStatus: ...

    def superseded(self, by_event_id: str) -> LedgerEntry: ...

    def acknowledged(self, by: str, note: str, at: datetime) -> LedgerEntry: ...


class EventStatus(StrEnum):
    OPEN = "open"
    ACKNOWLEDGED = "acknowledged"
    SUPERSEDED = "superseded"
    """A later run found a further change to the same variant before anyone read
    this one."""


@dataclass(frozen=True, slots=True)
class DriftEvent:
    """One surfaced change, with its cause and its routing decision attached."""

    event_id: str
    attribution: Attribution
    decision: SurfacingDecision
    raised_at: datetime | None = None
    status: EventStatus = EventStatus.OPEN
    superseded_by: str | None = None
    #: Who said they had read this, what they wrote, and when. Recorded rather
    #: than merely flipping a flag: "somebody dealt with it" is not something a
    #: laboratory can answer an audit with, and the note is where a reviewer
    #: says *why* a change needed nothing — which is the part a later reviewer
    #: needs and the part a boolean throws away.
    acknowledged_by: str = ""
    acknowledged_note: str = ""
    acknowledged_at: datetime | None = None

    @property
    def variant_key(self) -> str:
        return self.attribution.delta.variant_key

    @property
    def fingerprint(self) -> str:
        return self.attribution.delta.fingerprint

    @property
    def urgency(self) -> Urgency:
        return self.decision.urgency

    @property
    def queue(self) -> ReviewQueue:
        return self.decision.queue

    def acknowledged(self, by: str, note: str, at: datetime) -> DriftEvent:
        """Mark as read, by somebody, at a time.

        The fingerprint will not be raised again — which is exactly why the
        person has to be named. This is the one action in the package that makes
        the system *quieter*, and an anonymous switch that suppresses future
        alerts is the thing an incident review cannot reconstruct.
        """
        return _replace(
            self,
            status=EventStatus.ACKNOWLEDGED,
            acknowledged_by=by,
            acknowledged_note=note,
            acknowledged_at=at,
        )

    def superseded(self, by_event_id: str) -> DriftEvent:
        return _replace(self, status=EventStatus.SUPERSEDED, superseded_by=by_event_id)

    def summary(self) -> str:
        delta = self.attribution.delta
        return (
            f"{delta.variant_key}  {delta.before.mechanism} → {delta.after.mechanism}"
            f"  [{delta.kind}, {self.urgency}]  {self.decision.reason}"
        )


def _replace(event: DriftEvent, **changes: object) -> DriftEvent:
    fields = {
        "event_id": event.event_id,
        "attribution": event.attribution,
        "decision": event.decision,
        "raised_at": event.raised_at,
        "status": event.status,
        "superseded_by": event.superseded_by,
        "acknowledged_by": event.acknowledged_by,
        "acknowledged_note": event.acknowledged_note,
        "acknowledged_at": event.acknowledged_at,
    }
    fields.update(changes)
    return DriftEvent(**fields)  # type: ignore[arg-type]


@dataclass(slots=True)
class CaseLedger:
    """Every variant a case is watching, and everything that has happened."""

    case_id: str
    variant_keys: tuple[str, ...]
    #: The patient's own coordinate, which moves on a ward round rather than on
    #: a release schedule.
    phenotype: Pin
    events: list[LedgerEntry] = field(default_factory=list)
    last_world: World | None = None
    #: When a run last *looked* at this case, which is not when it last found
    #: something. The two were conflated until a dashboard tried to report on
    #: them, and the conflation is the dangerous direction: a case examined
    #: nightly for a year with nothing to report has no events, and reading
    #: that as "never run" cries wolf, while reading "never run" as "quiet"
    #: hides a scheduler that died. Only this field can tell them apart.
    last_examined_at: datetime | None = None

    @property
    def acknowledged_fingerprints(self) -> frozenset[str]:
        return frozenset(
            event.fingerprint for event in self.events if event.status is EventStatus.ACKNOWLEDGED
        )

    def open_events(self) -> tuple[LedgerEntry, ...]:
        return tuple(event for event in self.events if event.status is EventStatus.OPEN)

    def record(self, event: LedgerEntry) -> None:
        """Append, superseding anything still open for the same variant.

        Two unread alerts about one variant is one alert too many: the older
        describes a world that no longer exists, and leaving it in the queue
        invites a reviewer to sign out a transition already overtaken. It stays
        in the ledger — marked superseded, still auditable — but leaves the queue.
        """
        self.events = [
            existing.superseded(event.event_id)
            if existing.status is EventStatus.OPEN and existing.variant_key == event.variant_key
            else existing
            for existing in self.events
        ]
        self.events.append(event)

    def acknowledge(self, event_id: str, *, by: str, note: str, at: datetime) -> bool:
        """Record that a named person read this event.

        ``by`` is required and must not be blank. An acknowledgement suppresses
        every future alert carrying the same fingerprint, so it is the one
        gesture here that removes information from somebody's screen — and an
        unattributed one is a suppression nobody can be asked about.

        Only an open event can be acknowledged. Acknowledging a superseded one
        would sign for a transition that has already been overtaken, which is
        the specific mistake ``record`` supersedes events to prevent.
        """
        if not by.strip():
            raise ValueError(
                "an acknowledgement needs the name of whoever made it: it suppresses "
                "every future alert with this fingerprint, and an anonymous suppression "
                "is one nobody can be asked about later"
            )
        found = False
        updated = []
        for event in self.events:
            if event.event_id == event_id and event.status is EventStatus.OPEN:
                updated.append(event.acknowledged(by.strip(), note.strip(), at))
                found = True
            else:
                updated.append(event)
        self.events = updated
        return found


@dataclass(frozen=True, slots=True)
class ReanalysisReport:
    """What one run produced — the unit a scheduler logs and a UI renders."""

    case_id: str
    baseline: World
    candidate: World
    events: tuple[DriftEvent, ...]
    variants_examined: int
    rule_evaluations: int

    @property
    def moved_axes(self) -> tuple[str, ...]:
        return tuple(axis.value for axis in self.candidate.axes_differing_from(self.baseline))

    def by_queue(self, queue: ReviewQueue) -> tuple[DriftEvent, ...]:
        return tuple(event for event in self.events if event.queue is queue)

    def at_least(self, floor: Urgency) -> tuple[DriftEvent, ...]:
        return tuple(event for event in self.events if event.urgency.at_least(floor))

    @property
    def needs_a_human(self) -> bool:
        return any(event.decision.reaches_a_human for event in self.events)

    def headline(self) -> str:
        if not self.events:
            return f"{self.case_id}: {self.variants_examined} variants re-examined, nothing moved"
        return (
            f"{self.case_id}: {len(self.events)} change(s) across {self.variants_examined} "
            f"variants, {len(self.at_least(Urgency.HIGH))} needing prompt review "
            f"(axes moved: {', '.join(self.moved_axes) or 'none'})"
        )
