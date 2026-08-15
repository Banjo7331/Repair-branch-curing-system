"""The rule file is a document a geneticist edits, so the loader has to fail
loudly and specifically when they get something wrong.

The rule these tests protect: a mistake in the rule file must never produce a
rule that silently never fires. That failure mode is indistinguishable from a
rule that correctly never fires, and only one of the two is a bug.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from repairbench.features import FeatureRecord
from repairbench.model import Mechanism, Strength
from repairbench.ruleset import Rule, RulesetError, load_ruleset
from repairbench.transcript import NMDOutcome, NMDPrediction

SHIPPED_RULES = Path(__file__).parents[1] / "rules" / "mechanism-v1.yaml"


def features(**values) -> FeatureRecord:
    return FeatureRecord(
        values=values,
        nmd=NMDPrediction(NMDOutcome.NOT_APPLICABLE, "fixture"),
    )


def rule(when: dict) -> Rule:
    return Rule(
        id="TEST",
        supports=Mechanism.LOSS_OF_FUNCTION,
        strength=Strength.MODERATE,
        when=when,
        because="fixture",
    )


def test_the_shipped_ruleset_loads_and_pins_itself():
    ruleset = load_ruleset(SHIPPED_RULES)

    assert ruleset.version == "mechanism-v1"
    assert len(ruleset.digest) == 64
    assert ruleset.pin.startswith("mechanism-v1@")


def test_the_digest_changes_when_the_file_does(tmp_path: Path):
    """A call names the ruleset that produced it; two rulesets must not share a name."""
    original = load_ruleset(SHIPPED_RULES)
    edited = tmp_path / "edited.yaml"
    edited.write_text(SHIPPED_RULES.read_text() + "\n# a curator's comment\n")

    assert load_ruleset(edited).digest != original.digest


def test_every_shipped_rule_cites_something():
    """A rule without a citation is an assertion nobody can check."""
    for shipped in load_ruleset(SHIPPED_RULES).rules:
        assert shipped.citation.strip(), f"{shipped.id} has no citation"
        assert shipped.because.strip(), f"{shipped.id} has no reasoning"


def test_a_typo_in_a_feature_name_fails_the_run():
    with pytest.raises(RulesetError, match="unknown feature"):
        rule({"feature": "gene.multimer", "is": True}).fires(
            features(**{"gene.forms_multimer": True})
        )


def test_an_unknown_operator_fails_the_run():
    with pytest.raises(RulesetError, match="unknown operator"):
        rule({"feature": "x", "roughly": 3}).fires(features(x=3))


def test_a_leaf_with_two_operators_is_ambiguous_and_refused():
    with pytest.raises(RulesetError, match="exactly one operator"):
        rule({"feature": "x", "gt": 1, "lt": 5}).fires(features(x=3))


def test_combinators_compose():
    record = features(a=True, b=False, n=5)

    assert rule({"all": [{"feature": "a", "is": True}, {"feature": "n", "gte": 5}]}).fires(record)
    assert rule({"any": [{"feature": "b", "is": True}, {"feature": "n", "eq": 5}]}).fires(record)
    assert rule({"not": {"feature": "b", "is": True}}).fires(record)
    both_true = {"all": [{"feature": "a", "is": True}, {"feature": "b", "is": True}]}
    assert not rule(both_true).fires(record)


def test_a_missing_value_does_not_satisfy_a_comparison():
    """gnomAD constraint is often absent; a rule about it must not fire on None."""
    record = features(**{"gene.loeuf": None})

    assert not rule({"feature": "gene.loeuf", "lt": 0.35}).fires(record)
    assert rule({"feature": "gene.loeuf", "present": False}).fires(record)


def test_an_unknown_mechanism_in_a_rule_is_rejected(tmp_path: Path):
    broken = tmp_path / "broken.yaml"
    broken.write_text(
        "version: t\nrules:\n"
        "  - id: R\n    supports: telepathy\n    strength: strong\n"
        "    when: {feature: x, is: true}\n    because: nonsense\n"
    )

    with pytest.raises(RulesetError, match="unknown mechanism"):
        load_ruleset(broken)


def test_duplicate_rule_ids_are_rejected(tmp_path: Path):
    """Two rules with one id makes the evidence trail ambiguous."""
    duplicated = tmp_path / "dup.yaml"
    duplicated.write_text(
        "version: t\nrules:\n"
        "  - id: R\n    supports: loss_of_function\n    strength: strong\n"
        "    when: {feature: x, is: true}\n    because: one\n"
        "  - id: R\n    supports: gain_of_function\n    strength: strong\n"
        "    when: {feature: x, is: false}\n    because: two\n"
    )

    with pytest.raises(RulesetError, match="duplicate rule id"):
        load_ruleset(duplicated)
