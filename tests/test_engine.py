"""How the engine adds up what the rules say.

These tests use invented rulesets rather than the shipped one, because the
question here is arithmetic — curation outranking inference, cautions capping
confidence, conflicts being kept — and mixing that with real clinical rules
would make a failure ambiguous about which of the two broke.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from repairbench.engine import resolve
from repairbench.features import MechanismQuery, Variant
from repairbench.model import (
    Confidence,
    Consequence,
    DosageScore,
    Gene,
    Mechanism,
    MissenseDistribution,
)
from repairbench.ruleset import load_ruleset
from repairbench.transcript import Transcript

TRANSCRIPT = Transcript("NM_000001.1", "TESTG", (300, 300, 300, 300))


def query(**gene_kwargs) -> MechanismQuery:
    return MechanismQuery(
        variant=Variant(gene="TESTG", consequence=Consequence.NONSENSE, cds_position=500),
        transcript=TRANSCRIPT,
        gene=Gene(symbol="TESTG", **gene_kwargs),
    )


def ruleset_from(tmp_path: Path, body: str, name: str = "rules.yaml"):
    path = tmp_path / name
    path.write_text("version: test\nthresholds: {}\nrules:\n" + body)
    return load_ruleset(path)


def test_a_decisive_rule_ends_the_matter(tmp_path: Path):
    rules = ruleset_from(
        tmp_path,
        """
  - id: CURATED
    supports: loss_of_function
    strength: decisive
    when: {feature: gene.has_curated_mechanism, is: true}
    because: curation
  - id: INFERRED_OTHERWISE
    supports: gain_of_function
    strength: strong
    when: {feature: consequence.is_predicted_null, is: true}
    because: inference that disagrees
""",
    )

    call = resolve(
        query(curated_mechanism="loss_of_function", curated_mechanism_source="ClinGen"), rules
    )

    assert call.mechanism is Mechanism.LOSS_OF_FUNCTION
    assert call.confidence is Confidence.ESTABLISHED


def test_the_disagreeing_rule_is_kept_as_a_conflict_not_discarded(tmp_path: Path):
    """A curation describes the gene; a variant can be an exception to it, and a
    reviewer needs to see both facts."""
    rules = ruleset_from(
        tmp_path,
        """
  - id: CURATED
    supports: loss_of_function
    strength: decisive
    when: {feature: gene.has_curated_mechanism, is: true}
    because: curation
  - id: DISAGREES
    supports: dominant_negative
    strength: strong
    when: {feature: consequence.is_predicted_null, is: true}
    because: this variant may be an exception
""",
    )

    call = resolve(
        query(curated_mechanism="loss_of_function", curated_mechanism_source="ClinGen"), rules
    )

    assert [e.rule_id for e in call.conflicts] == ["DISAGREES"]
    assert call.needs_review


def test_one_supporting_rule_is_not_an_answer(tmp_path: Path):
    """Below the minimum, the engine says it does not know rather than laundering
    a single weak rule into a mechanism."""
    rules = ruleset_from(
        tmp_path,
        """
  - id: WEAK
    supports: loss_of_function
    strength: supporting
    when: {feature: consequence.is_predicted_null, is: true}
    because: barely anything
""",
    )

    call = resolve(query(), rules)

    assert call.mechanism is Mechanism.UNDETERMINED
    assert call.confidence is Confidence.NONE


def test_a_strong_caution_caps_confidence_at_possible(tmp_path: Path):
    """The rule file's way of saying 'not so fast' without saying 'no'."""
    body = """
  - id: STRONG_FOR_LOF
    supports: loss_of_function
    strength: strong
    when: {feature: consequence.is_predicted_null, is: true}
    because: strong
  - id: MORE_FOR_LOF
    supports: loss_of_function
    strength: moderate
    when: {feature: consequence.is_predicted_null, is: true}
    because: more
"""
    without_caution = resolve(query(), ruleset_from(tmp_path, body, "a.yaml"))
    with_caution = resolve(
        query(),
        ruleset_from(
            tmp_path,
            body
            + """
  - id: CAUTION
    supports: undetermined
    strength: strong
    when: {feature: consequence.is_predicted_null, is: true}
    because: the evidence does not settle it
""",
            "b.yaml",
        ),
    )

    assert without_caution.confidence is Confidence.PROBABLE
    assert with_caution.confidence is Confidence.POSSIBLE
    assert with_caution.mechanism is Mechanism.LOSS_OF_FUNCTION


def test_a_close_contest_does_not_reach_probable(tmp_path: Path):
    """Two mechanisms within one point of each other is not a confident answer."""
    rules = ruleset_from(
        tmp_path,
        """
  - id: FOR_LOF
    supports: loss_of_function
    strength: strong
    when: {feature: consequence.is_predicted_null, is: true}
    because: one side
  - id: FOR_DN
    supports: dominant_negative
    strength: moderate
    when: {feature: consequence.is_predicted_null, is: true}
    because: the other side
""",
    )

    call = resolve(query(), rules)

    assert call.mechanism is Mechanism.LOSS_OF_FUNCTION
    assert call.confidence is Confidence.POSSIBLE
    assert call.conflicts


def test_no_rule_firing_is_an_answer_of_its_own(tmp_path: Path):
    rules = ruleset_from(
        tmp_path,
        """
  - id: NEVER
    supports: loss_of_function
    strength: strong
    when: {feature: consequence.is_missense, is: true}
    because: not this variant
""",
    )

    call = resolve(query(), rules)

    assert call.mechanism is Mechanism.UNDETERMINED
    assert call.evidence == ()
    assert call.needs_review


def test_an_undetermined_call_warns_against_using_its_feasibility_flags(tmp_path: Path):
    """Every downstream modality depends on the mechanism; without one, the flags
    are not usable and the call says so."""
    rules = ruleset_from(
        tmp_path,
        """
  - id: NEVER
    supports: loss_of_function
    strength: strong
    when: {feature: consequence.is_missense, is: true}
    because: not this variant
""",
    )

    call = resolve(query(), rules)

    assert any("mechanism is unresolved" in note for note in call.feasibility.notes)


def test_gene_addition_is_coherent_only_for_loss_of_function(tmp_path: Path):
    """The flag that protects against the worst mistake this module could enable:
    supplying a normal copy where the mutant one is the problem."""
    for mechanism, coherent in (
        (Mechanism.LOSS_OF_FUNCTION, True),
        (Mechanism.GAIN_OF_FUNCTION, False),
        (Mechanism.DOMINANT_NEGATIVE, False),
        (Mechanism.SPLICING_DISRUPTION, False),
        (Mechanism.UNDETERMINED, False),
    ):
        assert mechanism.tolerates_gene_addition is coherent, mechanism


def test_dosage_vocabulary_separates_refutation_from_absence():
    """"No evidence yet" and "actively refuted" are different findings, and a
    predicted null variant means something different under each."""
    assert DosageScore.SUFFICIENT_EVIDENCE.supports_haploinsufficiency
    assert not DosageScore.NO_EVIDENCE.supports_haploinsufficiency
    assert not DosageScore.NO_EVIDENCE.refutes_haploinsufficiency
    assert DosageScore.UNLIKELY.refutes_haploinsufficiency
    assert DosageScore.AUTOSOMAL_RECESSIVE.refutes_haploinsufficiency


def test_a_hotspot_count_above_the_total_is_rejected():
    with pytest.raises(Exception, match="more hotspot missense"):
        MissenseDistribution(pathogenic_missense_total=3, pathogenic_missense_in_hotspot=5)
