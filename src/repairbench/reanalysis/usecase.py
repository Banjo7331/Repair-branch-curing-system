"""Re-examining one watched case against the current world.

The shape of it is the argument of the whole project in a dozen lines: assemble
today's world, re-assess each variant under it, diff against the last
assessment, attribute anything that moved, and let the surfacing policy decide
who — if anyone — hears about it.

Note what is *not* here. No thresholds, no "if ClinVar changed then alert", no
notification rules. Those are domain decisions and live in the rule files and in
``surfacing.py``, where they can be tested one at a time without a database.
"""

from __future__ import annotations

from datetime import datetime
from typing import Protocol

from repairbench.reanalysis.attribution import Attributor
from repairbench.reanalysis.drift import Assessment, compare
from repairbench.reanalysis.ledger import CaseLedger, DriftEvent, EventStatus, ReanalysisReport
from repairbench.reanalysis.surfacing import SurfacingPolicy
from repairbench.reanalysis.world import Pin, World


class AssessmentEngine(Protocol):
    """Runs the mechanism and modality rules for one variant, in one world.

    The contract that makes everything above it honest: the engine must *honour
    the world it is handed*. An implementation that quietly reads today's gene
    curation when asked for last year's makes every attribution in this module a
    lie — which is why the merged project is better than the two halves were.
    Before, this port had to be satisfied by somebody else's classifier and the
    contract was a hope; now the engine is the rule files in this same package,
    and honouring a world means loading the rule file that pin names.
    """

    def assess(self, variant_key: str, world: World) -> Assessment: ...


class SnapshotCatalog(Protocol):
    """Knows which releases exist and which are current."""

    def latest_global_pins(self) -> list[Pin]: ...


class CaseRepository(Protocol):
    def get(self, case_id: str) -> CaseLedger: ...

    def save(self, ledger: CaseLedger) -> None: ...


class AssessmentStore(Protocol):
    """Append-only history of everything ever concluded about a case."""

    def latest_for(self, case_id: str, variant_key: str) -> Assessment | None: ...

    def record(self, case_id: str, assessment: Assessment) -> None: ...


class Clock(Protocol):
    def now(self) -> datetime: ...


class IdentifierFactory(Protocol):
    def next_id(self) -> str: ...


class Notifier(Protocol):
    def publish(self, report: ReanalysisReport) -> None: ...


class _Counter:
    """Counts what a run cost in rule evaluations."""

    def __init__(self, engine: AssessmentEngine) -> None:
        self._engine = engine
        self.calls = 0

    def __call__(self, variant_key: str, world: World) -> Assessment:
        self.calls += 1
        return self._engine.assess(variant_key, world)


class ReanalyseCase:
    """Run one case against the current world and file whatever moved."""

    def __init__(
        self,
        engine: AssessmentEngine,
        catalog: SnapshotCatalog,
        cases: CaseRepository,
        assessments: AssessmentStore,
        clock: Clock,
        identifiers: IdentifierFactory,
        notifier: Notifier,
        policy: SurfacingPolicy | None = None,
        attributor: Attributor | None = None,
    ) -> None:
        self._engine = engine
        self._catalog = catalog
        self._cases = cases
        self._assessments = assessments
        self._clock = clock
        self._identifiers = identifiers
        self._notifier = notifier
        self._policy = policy or SurfacingPolicy()
        self._attributor = attributor or Attributor()

    def execute(self, case_id: str) -> ReanalysisReport:
        ledger = self._cases.get(case_id)
        candidate = self._assemble_world(ledger)
        counter = _Counter(self._engine)
        acknowledged = ledger.acknowledged_fingerprints

        events: list[DriftEvent] = []
        for variant_key in ledger.variant_keys:
            event = self._examine(ledger, variant_key, candidate, counter, acknowledged)
            if event is not None:
                ledger.record(event)
                events.append(event)

        baseline = ledger.last_world or candidate
        ledger.last_world = candidate
        # Recorded even when nothing moved — especially then. A run that found
        # nothing is the commonest correct outcome and the one indistinguishable
        # from no run at all, and this is the only thing that distinguishes them.
        ledger.last_examined_at = self._clock.now()
        self._cases.save(ledger)

        report = ReanalysisReport(
            case_id=case_id,
            baseline=baseline,
            candidate=candidate,
            events=tuple(events),
            variants_examined=len(ledger.variant_keys),
            rule_evaluations=counter.calls,
        )
        if report.needs_a_human:
            self._notifier.publish(report)
        return report

    def _assemble_world(self, ledger: CaseLedger) -> World:
        """Global releases plus this patient's phenotype.

        ``World`` refuses to be built with an axis missing, so a catalogue that
        forgets to publish one fails here rather than three steps later inside an
        attribution that would have blamed the wrong thing.
        """
        return World.of([*self._catalog.latest_global_pins(), ledger.phenotype])

    def _examine(
        self,
        ledger: CaseLedger,
        variant_key: str,
        candidate: World,
        counter: _Counter,
        acknowledged: frozenset[str],
    ) -> DriftEvent | None:
        previous = self._assessments.latest_for(ledger.case_id, variant_key)

        if previous is None:
            # First sight of this variant: establish the baseline, claim nothing.
            self._assessments.record(ledger.case_id, counter(variant_key, candidate))
            return None

        if previous.world.digest == candidate.digest:
            # Nothing this variant depends on has moved. The rules are not run at
            # all, which is what makes a nightly sweep over a cohort affordable.
            return None

        current = counter(variant_key, candidate)
        self._assessments.record(ledger.case_id, current)

        delta = compare(previous, current)
        if not delta.is_material:
            return None

        attribution = self._attributor.attribute(delta, counter)
        decision = self._policy.decide(attribution, acknowledged)
        if not decision.reaches_a_human:
            return None

        return DriftEvent(
            event_id=self._identifiers.next_id(),
            attribution=attribution,
            decision=decision,
            raised_at=self._clock.now(),
            status=EventStatus.OPEN,
        )
