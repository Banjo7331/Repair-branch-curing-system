"""Does the module pick the route that was actually taken?

Every case in ``reference/modalities.yaml`` is a disease where a therapy exists
or where the field has settled which route is available. The test runs the full
pipeline — mechanism rules, then modality rules on their output — so a failure
can be in either file. That is deliberate: M6 is only ever as good as the
mechanism it is handed.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from repairbench.cli import build_query
from repairbench.engine import resolve
from repairbench.modality import Modality, Verdict
from repairbench.modality_rules import load_modality_ruleset
from repairbench.model import Mechanism
from repairbench.ruleset import load_ruleset
from repairbench.selector import select

REFERENCE = Path(__file__).parent / "reference" / "modalities.yaml"
MECHANISM_RULES = Path(__file__).parents[1] / "rules" / "mechanism-v1.yaml"
MODALITY_RULES = Path(__file__).parents[1] / "rules" / "modality-v1.yaml"


def load_cases() -> list[dict]:
    return yaml.safe_load(REFERENCE.read_text())["cases"]


def run(case: dict):
    query = build_query(case)
    call = resolve(query, load_ruleset(MECHANISM_RULES))
    return call, select(call, query, load_modality_ruleset(MODALITY_RULES))


@pytest.mark.parametrize("case", load_cases(), ids=lambda c: c["name"])
def test_reference_case(case: dict) -> None:
    call, selection = run(case)
    expected = case["expect"]

    assert call.mechanism is Mechanism(expected["mechanism"]), (
        f"{case['name']}: mechanism came out as {call.mechanism}, "
        f"expected {expected['mechanism']} — the modality assessment below is moot"
    )

    if expected.get("blocked"):
        assert selection.is_blocked, f"{case['name']}: expected every modality to be blocked"
        assert all(
            a.verdict is Verdict.BLOCKED_BY_UNRESOLVED_MECHANISM for a in selection.assessments
        )
        return

    indicated = {a.modality for a in selection.indicated}
    for name in expected.get("indicated", []):
        modality = Modality(name)
        assert modality in indicated, (
            f"{case['name']}: expected {name} to be indicated, got "
            f"{selection.verdict_for(modality)}\n"
            f"  indicated: {[a.modality.value for a in selection.indicated]}"
        )

    if "caveat_mentions" in expected:
        assert any(expected["caveat_mentions"] in caveat for caveat in selection.caveats), (
            f"{case['name']}: expected a caveat mentioning "
            f"{expected['caveat_mentions']!r}, got {selection.caveats}"
        )

    for name in expected.get("contraindicated", []):
        modality = Modality(name)
        assert selection.verdict_for(modality) is Verdict.CONTRAINDICATED, (
            f"{case['name']}: expected {name} to be contraindicated, got "
            f"{selection.verdict_for(modality)}"
        )


def test_gene_addition_is_never_indicated_against_a_harmful_product() -> None:
    """The module's worst available mistake, asserted across the whole set.

    Supplying a normal copy where the mutant one is the problem leaves every
    poisoning subunit exactly where it was. No case, and no future rule, may
    produce that recommendation.
    """
    for case in load_cases():
        call, selection = run(case)
        if call.mechanism in {Mechanism.DOMINANT_NEGATIVE, Mechanism.GAIN_OF_FUNCTION}:
            assert selection.verdict_for(Modality.GENE_ADDITION) is Verdict.CONTRAINDICATED, (
                f"{case['name']}: gene addition was not contraindicated for a "
                f"{call.mechanism} mechanism"
            )


def test_a_known_zygosity_never_leaves_a_zygosity_caveat() -> None:
    """The caveat is for missing data, not for data the module dislikes.

    Asserted specifically rather than as "no caveats at all": these cases supply
    no affected tissue, so they correctly carry a tissue caveat, and a blanket
    assertion here would have to be weakened every time an honest caveat is
    added. Naming the one under test keeps it sharp.
    """
    for case in load_cases():
        if "zygosity" not in case["variant"]:
            continue
        _, selection = run(case)
        assert not [c for c in selection.caveats if "zygosity" in c], (
            f"{case['name']}: zygosity caveat raised despite a known zygosity"
        )


def test_no_modality_needing_a_wild_type_allele_is_offered_without_one() -> None:
    """The defect this whole change closes, asserted across the set."""
    ruleset = load_modality_ruleset(MODALITY_RULES)
    for case in load_cases():
        query = build_query(case)
        if query.variant.zygosity.leaves_a_wild_type_allele is not False:
            continue
        _, selection = run(case)
        for assessment in selection.indicated:
            assert assessment.modality not in ruleset.require_wild_type_allele, (
                f"{case['name']}: {assessment.modality} was offered to a "
                f"{query.variant.zygosity} patient with no intact allele"
            )


def test_every_indication_carries_its_reasoning() -> None:
    for case in load_cases():
        _, selection = run(case)
        for assessment in selection.indicated:
            assert assessment.indications, (
                f"{case['name']}: {assessment.modality} indicated with no rule"
            )
            for evidence in assessment.indications:
                assert evidence.because.strip()
                assert evidence.citation.strip()


def test_every_selection_names_the_ruleset() -> None:
    """A modality list that cannot name the rules behind it cannot be compared
    with the one produced next month."""
    for case in load_cases():
        _, selection = run(case)
        assert selection.ruleset_version.startswith("modality-v1@")


def test_every_case_carries_its_reasoning() -> None:
    for case in load_cases():
        assert case.get("note", "").strip(), f"{case['name']} has no note explaining why"
