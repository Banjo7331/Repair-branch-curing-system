"""The whole thing, against files on disk.

Everything else in the suite runs the rules against objects built in memory.
This runs them against a miniature deployment — two rule revisions, two
curation releases, an annotation, a state directory — because the claim that
matters most in this project is only testable that way.

The claim: *"we re-ran it with January's curation and the change did not
happen"* has to be literally what happened. A counterfactual that quietly fell
back to today's files would produce a confident causal story about an experiment
that was never run, and nothing in an in-memory test would catch it.
"""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

import pytest

from repairbench.model import Consequence, MissenseDistribution, Zygosity
from repairbench.observability import Metrics
from repairbench.reanalysis.attribution import Attributor
from repairbench.reanalysis.catalogue import SourceCatalogue
from repairbench.reanalysis.drift import DeltaKind, compare
from repairbench.reanalysis.engine import RepairbenchEngine, WatchedVariant
from repairbench.reanalysis.routing import ReviewQueue, Urgency
from repairbench.reanalysis.store import JsonAssessmentStore, JsonCaseRepository, StoreError
from repairbench.reanalysis.surfacing import SurfacingPolicy
from repairbench.reanalysis.usecase import ReanalyseCase
from repairbench.reanalysis.world import DriftAxis, Pin, World

DEPLOYMENT = Path(__file__).parents[1] / "data" / "deployment"
VARIANT = "PLUSG-c158"

# c.158 is chosen to sit between the two NMD boundaries the literature quotes.
# Under a 50 nt rule the stop is far enough from the last junction to trigger
# decay; under 55 nt it escapes. Same variant, same evidence, different rule.
WATCHED = WatchedVariant(
    key=VARIANT,
    gene="PLUSG",
    consequence=Consequence.NONSENSE,
    zygosity=Zygosity.HETEROZYGOUS,
    cds_position=158,
    distribution=MissenseDistribution(60, 15, 40),
)


@pytest.fixture
def catalogue() -> SourceCatalogue:
    return SourceCatalogue.load(DEPLOYMENT / "catalogue.yaml")


@pytest.fixture
def engine(catalogue: SourceCatalogue) -> RepairbenchEngine:
    return RepairbenchEngine(catalogue, {VARIANT: WATCHED})


def world_at(catalogue: SourceCatalogue, **versions: str) -> World:
    pins = [
        catalogue.pin_for(axis, versions.get(axis.value, default))
        for axis, default in (
            (DriftAxis.CLINVAR, "2026-06"),
            (DriftAxis.POPULATION_FREQUENCY, "v4.1"),
            (DriftAxis.GENE_CURATION, "2026-01"),
            (DriftAxis.PANEL, "epilepsy-3.2"),
            (DriftAxis.EXPRESSION, "gtex-v10"),
            (DriftAxis.ANNOTATION, "r1"),
            (DriftAxis.RULES, "v1"),
        )
    ]
    return World.of([*pins, Pin(axis=DriftAxis.PHENOTYPE, version="day-1", digest="hpo-day-1")])


# --------------------------------------------------------------------------
# The catalogue
# --------------------------------------------------------------------------


def test_a_version_the_deployment_no_longer_holds_is_refused(catalogue: SourceCatalogue):
    """The refusal that makes attribution trustworthy. Falling back to the
    current release would make every counterfactual a claim about an experiment
    nobody ran."""
    with pytest.raises(Exception, match="cannot reproduce"):
        catalogue.pin_for(DriftAxis.GENE_CURATION, "2019-01")


def test_a_directory_release_is_digested_as_a_whole(catalogue: SourceCatalogue):
    """The rule files and our own curation version together, because the same
    people edit them in the same review. One pin covers the combination, so a
    report cannot cite a mixture nobody reviewed."""
    v1 = catalogue.pin_for(DriftAxis.RULES, "v1")
    v2 = catalogue.pin_for(DriftAxis.RULES, "v2")

    assert v1.digest != v2.digest
    assert len(v1.digest) == 64


def test_every_non_case_scoped_axis_must_have_a_release(tmp_path: Path):
    partial = tmp_path / "catalogue.yaml"
    partial.write_text("clinvar:\n  - {version: x, path: c.txt}\n")
    (tmp_path / "c.txt").write_text("x")

    with pytest.raises(Exception, match="no releases for"):
        SourceCatalogue.load(partial)


# --------------------------------------------------------------------------
# The distinction the whole project exists for
# --------------------------------------------------------------------------


def test_a_rule_edit_and_a_recuration_reach_the_same_answer_by_different_routes(
    engine: RepairbenchEngine, catalogue: SourceCatalogue
):
    """The sharpest test in the suite.

    Two changes produce an identical outcome — a settled mechanism becomes
    unsettled — and they must not be reported alike. One is ClinGen refuting
    dosage sensitivity, which is the field learning something. The other is us
    moving the NMD boundary from 50 nt to 55, which is a defensible reading of
    the same literature and tells nobody anything about the patient.
    """
    before = engine.assess(VARIANT, world_at(catalogue))
    after_rule_edit = engine.assess(VARIANT, world_at(catalogue, rules="v2"))
    after_recuration = engine.assess(VARIANT, world_at(catalogue, gene_curation="2026-04"))

    assert before.mechanism.value == "loss_of_function"
    assert after_rule_edit.mechanism.value == "undetermined"
    assert after_recuration.mechanism.value == "undetermined"

    policy, attributor = SurfacingPolicy(), Attributor()
    rule_edit = policy.decide(
        attributor.attribute(compare(before, after_rule_edit), engine.assess)
    )
    recuration = policy.decide(
        attributor.attribute(compare(before, after_recuration), engine.assess)
    )

    assert rule_edit.queue is ReviewQueue.VALIDATION
    assert rule_edit.rule_change_caveat
    assert "not because the evidence did" in rule_edit.reason

    assert recuration.queue is ReviewQueue.CLINICAL_SIGNOUT
    assert recuration.urgency is Urgency.HIGH
    assert "caused by gene_curation" in recuration.reason


def test_the_counterfactual_actually_loads_the_older_files(
    engine: RepairbenchEngine, catalogue: SourceCatalogue
):
    """Not a mock. The probe reads January's TSV off the disk."""
    delta = compare(
        engine.assess(VARIANT, world_at(catalogue)),
        engine.assess(VARIANT, world_at(catalogue, gene_curation="2026-04", rules="v2")),
    )

    attribution = Attributor().attribute(delta, engine.assess)

    # Both axes independently produce the same outcome, so neither is necessary
    # and the honest report is that either would have done it.
    assert attribution.roles[DriftAxis.GENE_CURATION].is_causal
    assert attribution.roles[DriftAxis.RULES].is_causal
    assert "independently caused by" in attribution.explain()


# --------------------------------------------------------------------------
# Surviving between runs
# --------------------------------------------------------------------------


class Ids:
    def __init__(self) -> None:
        self.n = 0

    def next_id(self) -> str:
        self.n += 1
        return f"evt-{self.n:04d}"


class Clock:
    def now(self) -> datetime:
        return datetime(2026, 8, 15, 3, 0, tzinfo=UTC)


class Quiet:
    def publish(self, report: object) -> None:
        """Deliberately empty."""


def build_run(tmp_path: Path, engine: RepairbenchEngine, catalogue: SourceCatalogue, version: str):
    class PinnedCatalog:
        def latest_global_pins(self) -> list[Pin]:
            return [p for p in world_at(catalogue, gene_curation=version).pins
                    if not p.axis.is_case_scoped]

    cases = JsonCaseRepository(tmp_path)
    store = JsonAssessmentStore(tmp_path)
    return (
        ReanalyseCase(engine, PinnedCatalog(), cases, store, Clock(), Ids(), Quiet()),
        cases,
        store,
    )


def test_a_case_must_be_registered_before_it_can_be_run(tmp_path: Path, engine, catalogue):
    """A run will not invent the list of variants it is meant to be watching."""
    usecase, _, _ = build_run(tmp_path, engine, catalogue, "2026-01")

    with pytest.raises(StoreError, match="Register it first"):
        usecase.execute("NICU-014")


def test_state_survives_between_processes(tmp_path: Path, engine, catalogue):
    """The property that makes this runnable from cron at all: the second
    invocation compares against the first, which is on disk."""
    usecase, cases, _ = build_run(tmp_path, engine, catalogue, "2026-01")
    cases.register("NICU-014", [VARIANT], Pin(DriftAxis.PHENOTYPE, "day-1", "hpo-day-1"))

    first = usecase.execute("NICU-014")
    assert first.events == ()

    # A fresh set of objects, as a second cron invocation would have.
    later, cases_again, _ = build_run(tmp_path, engine, catalogue, "2026-04")
    report = later.execute("NICU-014")

    assert len(report.events) == 1
    assert report.events[0].attribution.delta.kind is DeltaKind.MECHANISM_LOST
    assert cases_again.get("NICU-014").open_events()


def test_an_assessment_is_never_overwritten(tmp_path: Path, engine, catalogue):
    """"What did we know, and when" has to survive being asked years later."""
    usecase, cases, _ = build_run(tmp_path, engine, catalogue, "2026-01")
    cases.register("NICU-014", [VARIANT], Pin(DriftAxis.PHENOTYPE, "day-1", "hpo-day-1"))
    usecase.execute("NICU-014")
    later, _, store_again = build_run(tmp_path, engine, catalogue, "2026-04")
    later.execute("NICU-014")

    history = store_again.history_for("NICU-014", VARIANT)

    assert len(history) == 2
    assert {entry["mechanism"] for entry in history} == {"loss_of_function", "undetermined"}


def test_a_stored_event_refuses_to_hand_back_a_causal_claim(tmp_path: Path, engine, catalogue):
    """A year-old attribution is not reconstructed, because re-presenting it as
    if it had just been established would misrepresent its footing."""
    usecase, cases, _ = build_run(tmp_path, engine, catalogue, "2026-01")
    cases.register("NICU-014", [VARIANT], Pin(DriftAxis.PHENOTYPE, "day-1", "hpo-day-1"))
    usecase.execute("NICU-014")
    later, cases_again, _ = build_run(tmp_path, engine, catalogue, "2026-04")
    later.execute("NICU-014")

    restored = cases_again.get("NICU-014").events[0]

    assert restored.fingerprint
    with pytest.raises(StoreError, match="not reconstructed"):
        _ = restored.attribution.delta


def test_an_acknowledged_transition_stays_acknowledged_across_processes(
    tmp_path: Path, engine, catalogue
):
    usecase, cases, _ = build_run(tmp_path, engine, catalogue, "2026-01")
    cases.register("NICU-014", [VARIANT], Pin(DriftAxis.PHENOTYPE, "day-1", "hpo-day-1"))
    usecase.execute("NICU-014")
    later, cases_again, _ = build_run(tmp_path, engine, catalogue, "2026-04")
    raised = later.execute("NICU-014").events[0]

    ledger = cases_again.get("NICU-014")
    ledger.acknowledge(raised.event_id)
    cases_again.save(ledger)

    # Back to January's curation and forward again: the identical transition.
    back, _, _ = build_run(tmp_path, engine, catalogue, "2026-01")
    back.execute("NICU-014")
    forward, _, _ = build_run(tmp_path, engine, catalogue, "2026-04")

    assert forward.execute("NICU-014").events == ()


# --------------------------------------------------------------------------
# Being watchable
# --------------------------------------------------------------------------


def test_metrics_expose_the_one_signal_that_matters(tmp_path: Path, engine, catalogue):
    """A scheduled job that stops running is invisible unless something measures
    its absence. Alert on the age of this gauge, not on an error rate."""
    usecase, cases, _ = build_run(tmp_path, engine, catalogue, "2026-01")
    cases.register("NICU-014", [VARIANT], Pin(DriftAxis.PHENOTYPE, "day-1", "hpo-day-1"))
    metrics = Metrics()

    metrics.run_completed(usecase.execute("NICU-014"), elapsed_seconds=0.4)
    exposed = metrics.expose()

    assert "repairbench_last_run_timestamp_seconds" in exposed
    assert 'repairbench_runs_total{outcome="quiet"} 1' in exposed
    assert "# TYPE repairbench_rule_evaluations_total counter" in exposed


def test_events_are_counted_by_transition_queue_and_urgency(tmp_path: Path, engine, catalogue):
    usecase, cases, _ = build_run(tmp_path, engine, catalogue, "2026-01")
    cases.register("NICU-014", [VARIANT], Pin(DriftAxis.PHENOTYPE, "day-1", "hpo-day-1"))
    usecase.execute("NICU-014")
    later, _, _ = build_run(tmp_path, engine, catalogue, "2026-04")
    metrics = Metrics()

    metrics.run_completed(later.execute("NICU-014"), elapsed_seconds=1.0)

    assert 'kind="mechanism_lost"' in metrics.expose()
    assert 'queue="clinical_signout"' in metrics.expose()


def test_the_state_directory_is_readable_without_this_code(tmp_path: Path, engine, catalogue):
    """Plain JSON on purpose. An operator debugging a scheduled run at 3am should
    not need the library to find out what it last concluded."""
    usecase, cases, _ = build_run(tmp_path, engine, catalogue, "2026-01")
    cases.register("NICU-014", [VARIANT], Pin(DriftAxis.PHENOTYPE, "day-1", "hpo-day-1"))
    usecase.execute("NICU-014")

    written = json.loads((tmp_path / "cases" / "NICU-014.json").read_text())

    assert written["variant_keys"] == [VARIANT]
    assert written["last_world"]
