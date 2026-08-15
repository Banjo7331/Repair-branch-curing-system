"""Why the assessment moved — established by experiment, not by heuristics.

The tempting way to attribute a change is to guess: ClinVar moved and the
mechanism moved, so ClinVar must be the reason. That guess is wrong often enough
to matter, because releases arrive in batches — a gnomAD refresh, a ClinGen
recuration, a new transcript version and a rule-file edit can all land in the
same week, and only one of them moved the needle.

So this does not guess. It re-runs the rules on **counterfactual worlds**: the
old world with exactly one axis advanced, and the new world with exactly one
axis held back. Two questions get an answer per axis — is it *sufficient*, is it
*necessary* — and the cost is at most two re-evaluations per moved axis.

What it buys is a claim that survives review: not "we think the rule edit did
this", but "we re-ran it with the old rule file and the change still happened".
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

from repairbench.modality import Modality
from repairbench.model import Confidence, Mechanism
from repairbench.reanalysis.drift import (
    Assessment,
    AssessmentDelta,
    Attribution,
    AttributionPattern,
    AxisRole,
)
from repairbench.reanalysis.world import DriftAxis, World

#: Re-run everything for one variant as it would have been assessed in a given
#: world. Implementations live behind ``ports.AssessmentEngine``.
Evaluator = Callable[[str, World], Assessment]

#: What "the same result" means when comparing a counterfactual to the endpoints.
#: The modality sets are included because a change that only moves them is still
#: a change somebody has to act on.
Outcome = tuple[Mechanism, Confidence, frozenset[Modality], frozenset[Modality]]


class NothingToAttributeError(ValueError):
    """Attribution was asked for a delta in which nothing changed.

    It costs a re-evaluation per moved axis. Asking for it when there is no
    movement is a caller bug, and answering "nothing caused nothing" would hide
    it.
    """


def _outcome(assessment: Assessment) -> Outcome:
    return (
        assessment.mechanism,
        assessment.confidence,
        assessment.indicated,
        assessment.contraindicated,
    )


@dataclass(slots=True)
class _Memo:
    """Wraps an evaluator so a world is only ever assessed once, and counts the
    calls that actually reached the rules."""

    evaluate: Evaluator
    variant_key: str
    calls: int = 0
    cache: dict[str, Assessment] = field(default_factory=dict)

    def seed(self, assessment: Assessment) -> None:
        """Register an assessment already in hand.

        Both endpoints are always seeded, so when a single axis moved every
        counterfactual world *is* one of the endpoints and the attribution costs
        nothing at all.
        """
        self.cache[assessment.world.digest] = assessment

    def __call__(self, world: World) -> Assessment:
        cached = self.cache.get(world.digest)
        if cached is not None:
            return cached
        self.calls += 1
        result = self.evaluate(self.variant_key, world)
        self.cache[world.digest] = result
        return result


class Attributor:
    """Assigns causal roles to the axes that moved between two worlds."""

    def attribute(self, delta: AssessmentDelta, evaluate: Evaluator) -> Attribution:
        if not delta.is_material:
            raise NothingToAttributeError(
                f"{delta.variant_key} did not change between "
                f"{delta.before.world.short_digest} and {delta.after.world.short_digest}"
            )

        moved = delta.moved_axes
        if not moved:
            # The assessment changed while every pin stood still. Either the
            # rules are non-deterministic or a pin is lying about its contents.
            # Both are defects, and neither is a clinical finding.
            return Attribution(delta=delta, roles={}, pattern=AttributionPattern.UNATTRIBUTED)

        probe = _Memo(evaluate, delta.variant_key)
        probe.seed(delta.before)
        probe.seed(delta.after)
        baseline, candidate = _outcome(delta.before), _outcome(delta.after)

        sufficient = self._sufficient(delta, moved, baseline, candidate, probe)
        necessary = self._necessary(delta, moved, candidate, probe, sufficient)

        roles: dict[DriftAxis, AxisRole] = {}
        for axis in moved:
            is_sufficient, is_necessary = axis in sufficient, axis in necessary
            if is_sufficient and is_necessary:
                roles[axis] = AxisRole.DECISIVE
            elif is_sufficient:
                roles[axis] = AxisRole.SUFFICIENT
            elif is_necessary:
                roles[axis] = AxisRole.NECESSARY
            else:
                roles[axis] = AxisRole.CONTRIBUTING

        if len(sufficient) == 1:
            pattern = AttributionPattern.SOLE
        elif len(sufficient) > 1:
            pattern = AttributionPattern.REDUNDANT
        else:
            pattern = AttributionPattern.INTERACTION

        return Attribution(
            delta=delta, roles=roles, pattern=pattern, evaluations_performed=probe.calls
        )

    def _sufficient(
        self,
        delta: AssessmentDelta,
        moved: tuple[DriftAxis, ...],
        baseline: Outcome,
        candidate: Outcome,
        probe: _Memo,
    ) -> set[DriftAxis]:
        """Advance one axis at a time from the old world; see who reproduces it."""
        found: set[DriftAxis] = set()
        for axis in moved:
            world = delta.before.world.with_pin(delta.after.world.pin_for(axis))
            result = _outcome(probe(world))
            if result == candidate and result != baseline:
                found.add(axis)
        return found

    def _necessary(
        self,
        delta: AssessmentDelta,
        moved: tuple[DriftAxis, ...],
        candidate: Outcome,
        probe: _Memo,
        sufficient: set[DriftAxis],
    ) -> set[DriftAxis]:
        """Hold one axis back from the new world; see whose absence prevents it.

        With a single moved axis the answer is arithmetic rather than empirical —
        reverting the only thing that moved must restore the only other world we
        have — so the rules are not re-run.
        """
        if len(moved) == 1:
            return set(moved) if sufficient else set()
        found: set[DriftAxis] = set()
        for axis in moved:
            world = delta.after.world.with_pin(delta.before.world.pin_for(axis))
            if _outcome(probe(world)) != candidate:
                found.add(axis)
        return found
