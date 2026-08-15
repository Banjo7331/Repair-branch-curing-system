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

    def acknowledged(self) -> LedgerEntry: ...


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

    def acknowledged(self) -> DriftEvent:
        """Mark as read. Its fingerprint will not be raised again."""
        return _replace(self, status=EventStatus.ACKNOWLEDGED)

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

    def acknowledge(self, event_id: str) -> bool:
        found = False
        updated = []
        for event in self.events:
            if event.event_id == event_id:
                updated.append(event.acknowledged())
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
