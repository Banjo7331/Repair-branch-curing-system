"""Where a change goes once it has been detected and explained.

The value objects only. The rules that produce a decision live in
``surfacing.py``, which keeps the ledger able to record a decision without
depending on the policy that made it.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum


class Urgency(StrEnum):
    """How quickly a change needs a human."""

    CRITICAL = "critical"
    HIGH = "high"
    ROUTINE = "routine"
    LOW = "low"
    SILENT = "silent"

    @property
    def rank(self) -> int:
        """Position on the scale; lower is more pressing."""
        return list(Urgency).index(self)

    def at_least(self, floor: Urgency) -> bool:
        return self.rank <= floor.rank


class ReviewQueue(StrEnum):
    """Where a change is filed."""

    CLINICAL_SIGNOUT = "clinical_signout"
    """A qualified reviewer looks at it against the patient's record."""

    VALIDATION = "validation"
    """A change in our own rules. Bioinformatics reviews it before anyone
    clinical does."""

    WATCHLIST = "watchlist"
    """No action; visible so the next move is not a surprise."""

    NONE = "none"
    """Nobody is told."""


@dataclass(frozen=True, slots=True)
class SurfacingDecision:
    """Where a change goes, and the sentence explaining why it went there."""

    urgency: Urgency
    queue: ReviewQueue
    reason: str
    suppressed: bool = False
    rule_change_caveat: bool = False

    @property
    def reaches_a_human(self) -> bool:
        return self.queue is not ReviewQueue.NONE and not self.suppressed
