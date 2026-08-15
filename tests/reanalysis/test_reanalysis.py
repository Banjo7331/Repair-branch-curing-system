"""The reanalysis half, end to end and rule by rule.

The engine here is real: it runs the shipped mechanism and modality rule files
against a gene context that varies with the world. That is the payoff of merging
the two halves — before, this port had to be satisfied by somebody else's
classifier and "the engine honours the world it is handed" was a hope written in
a docstring. Here it is arranged and asserted.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pytest

from repairbench.engine import resolve
from repairbench.features import MechanismQuery, Variant
from repairbench.modality import Modality
from repairbench.modality_rules import load_modality_ruleset
from repairbench.model import (
    Consequence,
    DosageScore,
    Gene,
    Mechanism,
    MissenseDistribution,
    Zygosity,
)
from repairbench.reanalysis.attribution import Attributor, NothingToAttributeError
from repairbench.reanalysis.drift import (
    Assessment,
    AssessmentDelta,
    Attribution,
    AttributionPattern,
    AxisRole,
    DeltaKind,
    compare,
)
from repairbench.reanalysis.ledger import CaseLedger, DriftEvent, EventStatus
from repairbench.reanalysis.routing import ReviewQueue, Urgency
from repairbench.reanalysis.surfacing import SurfacingPolicy
from repairbench.reanalysis.usecase import ReanalyseCase
from repairbench.reanalysis.world import (
    DriftAxis,
    IncompleteWorldError,
    Pin,
    PinConflictError,
    World,
    all_axes,
)
from repairbench.ruleset import load_ruleset
from repairbench.selector import select
from repairbench.transcript import Transcript

RULES = Path(__file__).parents[2] / "rules"
VARIANT = "GRCh38-2-166046000-C-T"


# --------------------------------------------------------------------------
# Fixtures
# --------------------------------------------------------------------------


def pin(axis: DriftAxis, version: str) -> Pin:
    return Pin(axis=axis, version=version, digest=f"{axis.value}-{version}")


def world(**versions: str) -> World:
    return World.of(pin(axis, versions.get(axis.value, "1")) for axis in all_axes())


class RuleDrivenEngine:
    """Runs the shipped rule files against a case whose facts depend on the world.

    Two axes are wired to something real, and both are the kind of change that
    looks like nothing and inverts an answer.

    ``annotation`` at version 2 is a new transcript release that moves the exon
    boundaries. The same stop codon that sat comfortably upstream of the last
    junction now sits within fifty nucleotides of it, so the transcript escapes
    decay, a truncated protein is made, and a collagen-like gene goes from
    haploinsufficiency to dominant-negative. Nobody learned anything about the
    patient; the therapy inverts anyway. That is why annotation is an axis.

    ``gene_curation`` at version 2 refutes dosage sensitivity, which takes a
    settled loss-of-function call back to undetermined.
    """

    def __init__(self) -> None:
        self.mechanism_rules = load_ruleset(RULES / "mechanism-v1.yaml")
        self.modality_rules = load_modality_ruleset(RULES / "modality-v1.yaml")
        self.calls = 0

    def assess(self, variant_key: str, world_: World) -> Assessment:
        self.calls += 1
        # v1: the stop at c.400 sits 200 nt upstream of the last junction at 600.
        # v2: the junction moves to 400, so the same stop escapes decay.
        exons = (
            (300, 100, 500)
            if world_.pin_for(DriftAxis.ANNOTATION).version == "2"
            else (300, 300, 300)
        )
        dosage = (
            DosageScore.UNLIKELY
            if world_.pin_for(DriftAxis.GENE_CURATION).version == "2"
            else DosageScore.SUFFICIENT_EVIDENCE
        )
        gene = Gene(
            symbol="TESTG",
            haploinsufficiency=dosage,
            loeuf=0.2,
            forms_multimer=True,
            truncating_variants_are_milder=True,
            distribution=MissenseDistribution(
                pathogenic_missense_total=60,
                pathogenic_missense_in_hotspot=15,
                pathogenic_truncating_total=40,
            ),
        )
        query = MechanismQuery(
            variant=Variant(
                gene="TESTG",
                consequence=Consequence.NONSENSE,
                cds_position=400,
                zygosity=Zygosity.HETEROZYGOUS,
            ),
            transcript=Transcript("NM_000001.1", "TESTG", exons, mane_select=True),
            gene=gene,
        )
        call = resolve(query, self.mechanism_rules)
        return Assessment.of(variant_key, world_, call, select(call, query, self.modality_rules))


@pytest.fixture
def engine() -> RuleDrivenEngine:
    return RuleDrivenEngine()


# --------------------------------------------------------------------------
# The world
# --------------------------------------------------------------------------


def test_a_world_missing_an_axis_is_rejected():
    """An assessment that cannot name everything it read is not comparable to
    any other, so it must not be constructible."""
    partial = [pin(axis, "1") for axis in all_axes() if axis is not DriftAxis.ANNOTATION]

    with pytest.raises(IncompleteWorldError, match="annotation"):
        World.of(partial)


def test_the_same_version_with_two_digests_is_a_conflict():
    contradiction = Pin(axis=DriftAxis.CLINVAR, version="1", digest="something-else")

    with pytest.raises(PinConflictError, match="clinvar"):
        World.of([*[pin(axis, "1") for axis in all_axes()], contradiction])


def test_the_digest_ignores_pin_order():
    forwards = World.of(pin(axis, "1") for axis in all_axes())
    backwards = World.of(reversed([pin(axis, "1") for axis in all_axes()]))

    assert forwards.digest == backwards.digest


def test_only_the_moved_axes_are_reported():
    after = world(clinvar="2", rules="2")

    assert after.axes_differing_from(world()) == (DriftAxis.CLINVAR, DriftAxis.RULES)


def test_rules_is_the_only_non_clinical_axis():
    """The whole reporting layer leans on this distinction."""
    assert [axis for axis in all_axes() if not axis.is_clinical] == [DriftAxis.RULES]


def test_annotation_is_an_axis_because_a_transcript_version_can_invert_an_answer():
    """A new transcript version moves the exon boundaries, which moves the NMD
    boundary, which can change the mechanism without anybody learning anything."""
    assert DriftAxis.ANNOTATION in all_axes()
    assert DriftAxis.ANNOTATION.is_clinical


# --------------------------------------------------------------------------
# What a change means
# --------------------------------------------------------------------------


def test_a_new_transcript_version_inverts_the_mechanism(engine: RuleDrivenEngine):
    """The signal this project exists to catch, and the case that justifies the
    annotation axis. Nothing was learned about the patient — the exon boundaries
    moved — and supplementation goes from being the answer to being the hazard."""
    before = engine.assess(VARIANT, world())
    after = engine.assess(VARIANT, world(annotation="2"))

    delta = compare(before, after)

    assert before.mechanism is Mechanism.LOSS_OF_FUNCTION
    assert after.mechanism is Mechanism.DOMINANT_NEGATIVE
    assert delta.kind is DeltaKind.MECHANISM_INVERTED


def test_a_refuted_dosage_curation_takes_a_settled_answer_away(engine: RuleDrivenEngine):
    """The other direction, and it is not a lesser event: we knew, and now we do
    not, so everything downstream is unsupported until we do again."""
    delta = compare(
        engine.assess(VARIANT, world()), engine.assess(VARIANT, world(gene_curation="2"))
    )

    assert delta.kind is DeltaKind.MECHANISM_LOST


def test_a_withdrawn_route_is_distinguished_from_one_that_merely_dropped_off(
    engine: RuleDrivenEngine,
):
    """A route that is now *contraindicated* is a warning; one that quietly left
    the list is only a change."""
    before = engine.assess(VARIANT, world())
    after = engine.assess(VARIANT, world(annotation="2"))

    delta = compare(before, after)

    assert set(delta.withdrawn_modalities) <= before.indicated
    for modality in delta.withdrawn_modalities:
        assert modality in after.contraindicated


def test_diffing_two_different_variants_is_refused(engine: RuleDrivenEngine):
    first = engine.assess(VARIANT, world())
    other = engine.assess("GRCh38-16-9800000-T-C", world())

    with pytest.raises(ValueError, match="cannot diff"):
        compare(first, other)


def test_a_fingerprint_identifies_the_transition_not_the_variant(engine: RuleDrivenEngine):
    once = compare(engine.assess(VARIANT, world()), engine.assess(VARIANT, world(annotation="2")))
    again = compare(
        engine.assess(VARIANT, world(clinvar="9")),
        engine.assess(VARIANT, world(clinvar="9", annotation="2")),
    )

    assert once.fingerprint == again.fingerprint


# --------------------------------------------------------------------------
# Attribution
# --------------------------------------------------------------------------


def test_one_axis_moved_and_it_is_the_cause(engine: RuleDrivenEngine):
    """Costs no rule evaluations at all — both endpoints are already in hand."""
    delta = compare(engine.assess(VARIANT, world()), engine.assess(VARIANT, world(annotation="2")))

    result = Attributor().attribute(delta, engine.assess)

    assert result.pattern is AttributionPattern.SOLE
    assert result.roles[DriftAxis.ANNOTATION] is AxisRole.DECISIVE
    assert result.evaluations_performed == 0
    assert result.explain() == "caused by annotation"


def test_an_axis_that_moved_without_mattering_is_not_blamed(engine: RuleDrivenEngine):
    """Releases arrive in batches; being in the batch is not being the cause."""
    delta = compare(
        engine.assess(VARIANT, world()),
        engine.assess(VARIANT, world(annotation="2", clinvar="2", population_frequency="2")),
    )

    result = Attributor().attribute(delta, engine.assess)

    assert result.roles[DriftAxis.ANNOTATION] is AxisRole.DECISIVE
    assert result.roles[DriftAxis.CLINVAR] is AxisRole.CONTRIBUTING
    assert result.roles[DriftAxis.POPULATION_FREQUENCY] is AxisRole.CONTRIBUTING
    assert result.primary_axes == (DriftAxis.ANNOTATION,)


def test_a_rule_edit_landing_in_the_same_week_does_not_borrow_the_credit(
    engine: RuleDrivenEngine,
):
    """Half of the "our corrections are not discoveries" rule, and the half that
    is easy to get wrong: the rules axis moved, and it still did nothing."""
    delta = compare(
        engine.assess(VARIANT, world()),
        engine.assess(VARIANT, world(annotation="2", rules="2")),
    )

    result = Attributor().attribute(delta, engine.assess)

    assert result.roles[DriftAxis.RULES] is AxisRole.CONTRIBUTING
    assert not result.is_rule_change_implicated
    assert not result.is_purely_our_rules


def test_a_change_with_no_moved_axis_is_a_defect_not_a_finding(engine: RuleDrivenEngine):
    """Same world, different answer: the rules are not deterministic."""
    before = engine.assess(VARIANT, world())
    moved = engine.assess(VARIANT, world(annotation="2"))
    # Same world, different conclusions: the only way to reach this is for the
    # rules to be non-deterministic, which is a defect rather than a finding.
    after = Assessment(
        variant_key=VARIANT,
        world=before.world,
        mechanism=moved.mechanism,
        confidence=moved.confidence,
        indicated=moved.indicated,
        contraindicated=moved.contraindicated,
    )

    result = Attributor().attribute(compare(before, after), engine.assess)

    assert result.pattern is AttributionPattern.UNATTRIBUTED
    assert result.roles == {}


def test_attribution_is_refused_when_nothing_changed(engine: RuleDrivenEngine):
    delta = compare(engine.assess(VARIANT, world()), engine.assess(VARIANT, world(panel="2")))

    with pytest.raises(NothingToAttributeError):
        Attributor().attribute(delta, engine.assess)


def test_cost_is_bounded_by_the_number_of_moved_axes(engine: RuleDrivenEngine):
    delta = compare(
        engine.assess(VARIANT, world()),
        engine.assess(VARIANT, world(annotation="2", clinvar="2", population_frequency="2")),
    )

    result = Attributor().attribute(delta, engine.assess)

    assert result.evaluations_performed <= 2 * len(delta.moved_axes)


# --------------------------------------------------------------------------
# Surfacing
# --------------------------------------------------------------------------


def attribution_for(kind: DeltaKind, roles=None, engine: RuleDrivenEngine | None = None):
    """Build an attribution with a chosen kind, for testing the policy alone."""
    engine = engine or RuleDrivenEngine()
    before = engine.assess(VARIANT, world())
    after = engine.assess(VARIANT, world(annotation="2"))
    delta = AssessmentDelta(variant_key=VARIANT, before=before, after=after, kind=kind)
    return Attribution(
        delta=delta,
        roles=roles if roles is not None else {DriftAxis.ANNOTATION: AxisRole.DECISIVE},
        pattern=AttributionPattern.SOLE,
    )


def test_an_inverted_mechanism_is_the_most_urgent_thing_the_system_can_say():
    decision = SurfacingPolicy().decide(attribution_for(DeltaKind.MECHANISM_INVERTED))

    assert decision.urgency is Urgency.CRITICAL
    assert decision.queue is ReviewQueue.CLINICAL_SIGNOUT
    assert "inverted" in decision.reason


def test_a_withdrawn_route_is_as_urgent_as_an_inversion():
    """It may already be in a plan, which an opened route cannot be."""
    withdrawn = SurfacingPolicy().decide(attribution_for(DeltaKind.MODALITY_WITHDRAWN))
    opened = SurfacingPolicy().decide(attribution_for(DeltaKind.MODALITY_OPENED))

    assert withdrawn.urgency is Urgency.CRITICAL
    assert opened.urgency is Urgency.ROUTINE


def test_our_own_rule_change_never_reaches_the_clinical_queue_first():
    roles = {DriftAxis.RULES: AxisRole.DECISIVE}

    decision = SurfacingPolicy().decide(attribution_for(DeltaKind.MECHANISM_INVERTED, roles))

    assert decision.queue is ReviewQueue.VALIDATION
    assert decision.urgency is Urgency.CRITICAL
    assert decision.rule_change_caveat


def test_a_mixed_cause_goes_to_the_clinic_but_carries_the_caveat():
    roles = {DriftAxis.RULES: AxisRole.NECESSARY, DriftAxis.CLINVAR: AxisRole.NECESSARY}

    decision = SurfacingPolicy().decide(attribution_for(DeltaKind.MECHANISM_INVERTED, roles))

    assert decision.queue is ReviewQueue.CLINICAL_SIGNOUT
    assert decision.rule_change_caveat


def test_an_already_signed_out_transition_is_not_raised_again():
    attribution = attribution_for(DeltaKind.MECHANISM_INVERTED)

    decision = SurfacingPolicy().decide(attribution, frozenset({attribution.delta.fingerprint}))

    assert decision.suppressed
    assert not decision.reaches_a_human


def test_every_delta_kind_is_routed_deliberately():
    """A guard against a future kind quietly defaulting into the clinical queue."""
    policy = SurfacingPolicy()
    clinical = {}
    for kind in DeltaKind:
        decision = policy.decide(attribution_for(kind))
        assert decision.reason, f"{kind} was routed without a reason"
        clinical[kind] = decision.queue is ReviewQueue.CLINICAL_SIGNOUT

    assert {kind for kind, is_clinical in clinical.items() if is_clinical} == {
        DeltaKind.MECHANISM_INVERTED,
        DeltaKind.MODALITY_WITHDRAWN,
        DeltaKind.MECHANISM_RESOLVED,
        DeltaKind.MECHANISM_LOST,
        DeltaKind.MODALITY_OPENED,
    }


# --------------------------------------------------------------------------
# The run
# --------------------------------------------------------------------------


class Catalog:
    def __init__(self) -> None:
        self.versions = {
            axis.value: "1" for axis in all_axes() if not axis.is_case_scoped
        }

    def latest_global_pins(self) -> list[Pin]:
        return [
            pin(axis, self.versions[axis.value])
            for axis in all_axes()
            if not axis.is_case_scoped
        ]

    def advance(self, axis: DriftAxis, version: str) -> None:
        self.versions[axis.value] = version


class Cases:
    def __init__(self, ledger: CaseLedger) -> None:
        self.ledger = ledger

    def get(self, case_id: str) -> CaseLedger:
        return self.ledger

    def save(self, ledger: CaseLedger) -> None:
        self.ledger = ledger


class Store:
    def __init__(self) -> None:
        self.latest: dict[str, Assessment] = {}

    def latest_for(self, case_id: str, variant_key: str) -> Assessment | None:
        return self.latest.get(f"{case_id}|{variant_key}")

    def record(self, case_id: str, assessment: Assessment) -> None:
        self.latest[f"{case_id}|{assessment.variant_key}"] = assessment


class Ids:
    def __init__(self) -> None:
        self.n = 0

    def next_id(self) -> str:
        self.n += 1
        return f"evt-{self.n:04d}"


class FixedClock:
    def now(self) -> datetime:
        return datetime(2026, 8, 15, 3, 0, 0)


class Publisher:
    def __init__(self) -> None:
        self.reports = []

    def publish(self, report) -> None:
        self.reports.append(report)


@pytest.fixture
def harness(engine: RuleDrivenEngine):
    catalog = Catalog()
    cases = Cases(
        CaseLedger(
            case_id="NICU-014",
            variant_keys=(VARIANT,),
            phenotype=pin(DriftAxis.PHENOTYPE, "hpo-day-1"),
        )
    )
    store, publisher = Store(), Publisher()
    usecase = ReanalyseCase(engine, catalog, cases, store, FixedClock(), Ids(), publisher)
    return usecase, engine, catalog, cases, publisher


def test_the_first_run_establishes_a_baseline_and_claims_nothing(harness):
    """A variant seen for the first time has nothing to compare against, and
    reporting it as changed would be the easiest possible false positive."""
    usecase, engine, _, _, _ = harness

    report = usecase.execute("NICU-014")

    assert report.events == ()
    assert engine.calls == 1


def test_a_run_against_an_unchanged_world_costs_nothing(harness):
    """What makes a nightly sweep over a cohort affordable: if no pin moved, the
    rules are not run at all."""
    usecase, engine, _, _, _ = harness
    usecase.execute("NICU-014")
    before = engine.calls

    report = usecase.execute("NICU-014")

    assert engine.calls == before
    assert report.events == ()


def test_a_recuration_raises_an_attributed_critical_event(harness):
    usecase, _, catalog, _, publisher = harness
    usecase.execute("NICU-014")

    catalog.advance(DriftAxis.ANNOTATION, "2")
    report = usecase.execute("NICU-014")

    assert len(report.events) == 1
    event = report.events[0]
    assert event.attribution.delta.kind is DeltaKind.MECHANISM_INVERTED
    assert event.urgency is Urgency.CRITICAL
    assert event.attribution.roles[DriftAxis.ANNOTATION] is AxisRole.DECISIVE
    assert len(publisher.reports) == 1


def test_the_report_names_the_world_it_ran_in(harness):
    """Without it the report is uncitable: a reviewer reading it in a year cannot
    tell which releases produced it."""
    usecase, _, _, cases, _ = harness

    report = usecase.execute("NICU-014")

    assert "gene_curation@1" in report.candidate.describe()
    assert cases.ledger.last_world is not None


def test_a_second_change_supersedes_an_unread_first(harness):
    """Two unread alerts about one variant is one alert too many: the older
    describes a world that no longer exists."""
    usecase, _, catalog, cases, _ = harness
    usecase.execute("NICU-014")
    catalog.advance(DriftAxis.ANNOTATION, "2")
    usecase.execute("NICU-014")
    catalog.advance(DriftAxis.ANNOTATION, "1")
    usecase.execute("NICU-014")

    statuses = [event.status for event in cases.ledger.events]

    assert statuses.count(EventStatus.SUPERSEDED) == 1
    assert len(cases.ledger.open_events()) == 1


def test_an_acknowledged_transition_does_not_resurface(harness):
    usecase, _, catalog, cases, _ = harness
    usecase.execute("NICU-014")
    catalog.advance(DriftAxis.ANNOTATION, "2")
    first = usecase.execute("NICU-014")
    cases.ledger.acknowledge(first.events[0].event_id)

    catalog.advance(DriftAxis.ANNOTATION, "1")
    usecase.execute("NICU-014")
    catalog.advance(DriftAxis.ANNOTATION, "2")
    again = usecase.execute("NICU-014")

    assert again.events == ()


def test_an_event_carries_its_cause_into_the_ledger(harness):
    usecase, _, catalog, _, _ = harness
    usecase.execute("NICU-014")
    catalog.advance(DriftAxis.ANNOTATION, "2")

    event: DriftEvent = usecase.execute("NICU-014").events[0]

    assert "caused by annotation" in event.decision.reason
    assert Modality.GENE_ADDITION in event.attribution.delta.withdrawn_modalities
