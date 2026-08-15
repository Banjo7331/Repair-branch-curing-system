"""The validation that matters: does the rule file reproduce the literature?

Every case in ``reference/mechanisms.yaml`` is a gene whose mechanism the field
already agrees on. If a case fails, the first question is not "which assertion
do I relax" but "which rule is wrong" — the reference set is the specification
and the rule file is the implementation, not the other way round.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from repairbench.engine import resolve
from repairbench.features import MechanismQuery, SplicePrediction, Variant
from repairbench.model import (
    Confidence,
    Consequence,
    DosageScore,
    Gene,
    Imprinting,
    Mechanism,
    MissenseDistribution,
)
from repairbench.ruleset import load_ruleset
from repairbench.transcript import Transcript

REFERENCE = Path(__file__).parent / "reference" / "mechanisms.yaml"
RULES = Path(__file__).parents[1] / "rules" / "mechanism-v1.yaml"


def load_cases() -> list[dict]:
    return yaml.safe_load(REFERENCE.read_text())["cases"]


def build_query(case: dict) -> MechanismQuery:
    variant_spec = case["variant"]
    transcript_spec = case["transcript"]
    gene_spec = dict(case["gene"])

    exon_lengths = transcript_spec.get("exon_lengths")
    if exon_lengths is None:
        exon_lengths = uniform_exons(
            transcript_spec["exon_count"], transcript_spec["coding_length"]
        )
    transcript = Transcript(
        accession=transcript_spec["accession"],
        gene=variant_spec["gene"],
        coding_exon_lengths=tuple(exon_lengths),
        mane_select=True,
    )

    distribution = MissenseDistribution(**gene_spec.pop("distribution", {}))
    gene = Gene(
        symbol=variant_spec["gene"],
        haploinsufficiency=DosageScore(gene_spec.pop("haploinsufficiency", "no_evidence")),
        triplosensitivity=DosageScore(gene_spec.pop("triplosensitivity", "no_evidence")),
        loeuf=gene_spec.pop("loeuf", None),
        forms_multimer=gene_spec.pop("forms_multimer", False),
        truncating_variants_are_milder=gene_spec.pop("truncating_variants_are_milder", False),
        imprinting=Imprinting(gene_spec.pop("imprinting", "not_imprinted")),
        silenced_allele_intact=gene_spec.pop("silenced_allele_intact", False),
        distribution=distribution,
        curated_mechanism=gene_spec.pop("curated_mechanism", None),
        curated_mechanism_source=gene_spec.pop("curated_mechanism_source", None),
    )
    assert not gene_spec, f"unused gene fields in fixture: {sorted(gene_spec)}"

    return MechanismQuery(
        variant=Variant(
            gene=variant_spec["gene"],
            consequence=Consequence(variant_spec["consequence"]),
            cds_position=variant_spec["cds_position"],
            protein_change=variant_spec.get("protein_change", ""),
            hgvs_c=variant_spec.get("hgvs_c", ""),
        ),
        transcript=transcript,
        gene=gene,
        splice=SplicePrediction(max_delta=case.get("splice_max_delta")),
    )


def uniform_exons(count: int, coding_length: int) -> list[int]:
    """Spread a coding length over a number of exons, remainder on the last."""
    base = coding_length // count
    lengths = [base] * count
    lengths[-1] += coding_length - base * count
    return lengths


@pytest.mark.parametrize("case", load_cases(), ids=lambda c: c["name"])
def test_reference_case(case: dict) -> None:
    ruleset = load_ruleset(RULES)

    call = resolve(build_query(case), ruleset)

    expected = case["expect"]
    assert call.mechanism is Mechanism(expected["mechanism"]), (
        f"{case['name']}\n"
        f"  expected {expected['mechanism']}, got {call.mechanism}\n"
        f"  evidence: {[e.rule_id for e in call.evidence]}\n"
        f"  conflicts: {[e.rule_id for e in call.conflicts]}"
    )

    if "confidence" in expected:
        assert call.confidence is Confidence(expected["confidence"]), (
            f"{case['name']}: confidence {call.confidence}, expected {expected['confidence']}\n"
            f"  evidence: {[(e.rule_id, e.strength.value) for e in call.evidence]}"
        )

    if expected.get("has_conflicts"):
        assert call.conflicts, f"{case['name']}: expected the call to record a conflict"

    for flag in (
        "gene_addition_coherent",
        "fits_viral_payload",
        "exon_skipping_preserves_frame",
        "silenced_allele_available",
        "allele_specific_silencing_indicated",
    ):
        if flag in expected:
            assert getattr(call.feasibility, flag) == expected[flag], (
                f"{case['name']}: feasibility.{flag} = "
                f"{getattr(call.feasibility, flag)}, expected {expected[flag]}"
            )


def test_every_case_carries_its_reasoning() -> None:
    """A reference case without a note is a case nobody can review."""
    for case in load_cases():
        assert case.get("note", "").strip(), f"{case['name']} has no note explaining why"


def test_a_determined_call_always_cites_a_rule() -> None:
    """No mechanism without evidence — the discipline the whole module rests on."""
    ruleset = load_ruleset(RULES)
    for case in load_cases():
        call = resolve(build_query(case), ruleset)
        if call.mechanism.is_determined:
            assert call.evidence, f"{case['name']}: determined mechanism with no supporting rule"
            assert call.ruleset_version, f"{case['name']}: call does not name the ruleset"
