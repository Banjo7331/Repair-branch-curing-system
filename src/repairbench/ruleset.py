"""Loading, validating and pinning the rule file.

The rule file is the product. This module is the smallest thing that can read it
safely, and its design goals are in that order: safe, then small.

*Safe* means the predicate language is interpreted, not executed. There is no
``eval``, no expression string, and no way for a rule to reach anything the
feature record does not name — so a rule file is a document a geneticist can be
handed without it also being code they are running.

*Pinned* means a ruleset carries the digest of the bytes it was loaded from. A
mechanism call names the ruleset that produced it, for the same reason a
reanalysis names its knowledge snapshot: two calls made under different rules
are not comparable, and a system that cannot tell them apart will eventually
report a rule change as a biological finding.
"""

from __future__ import annotations

import hashlib
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Protocol

import yaml

from repairbench.model import Mechanism, RepairbenchError, Strength


class Features(Protocol):
    """Anything a rule can be evaluated against.

    A Protocol rather than the mechanism feature record, because three rule
    files now share this interpreter — mechanism, modality and off-target risk —
    and only the first two describe a variant. What they have in common is the
    contract that makes the language safe: a named lookup that *refuses* unknown
    names, and the ability to list what it does know so the refusal is useful.
    """

    def get(self, name: str) -> Any: ...

    def known_names(self) -> list[str]: ...


class RulesetError(RepairbenchError):
    """The rule file is malformed, or a rule refers to something that does not exist."""


#: Leaf operators. Each takes the feature value and the operand.
_OPERATORS: dict[str, Callable[[Any, Any], bool]] = {
    "eq": lambda value, operand: value == operand,
    "ne": lambda value, operand: value != operand,
    "is": lambda value, operand: value is operand,
    "gt": lambda value, operand: value is not None and value > operand,
    "gte": lambda value, operand: value is not None and value >= operand,
    "lt": lambda value, operand: value is not None and value < operand,
    "lte": lambda value, operand: value is not None and value <= operand,
    "in": lambda value, operand: value in operand,
    "present": lambda value, operand: (value is not None) is operand,
}


@dataclass(frozen=True, slots=True)
class Rule:
    """One clinical claim, in a form the engine can apply and a human can check."""

    id: str
    supports: Mechanism
    strength: Strength
    when: dict[str, Any]
    because: str
    citation: str = ""

    @property
    def is_decisive(self) -> bool:
        """A decisive rule ends the matter — it is how curation outranks inference."""
        return self.strength is Strength.DECISIVE

    def fires(self, features: Features) -> bool:
        return evaluate_predicate(self.when, features, self.id)


@dataclass(frozen=True, slots=True)
class Thresholds:
    """The numbers the engine and the feature builder need, all declared in the
    rule file rather than scattered through the code."""

    nmd_boundary_nt: int = 50
    reinitiation_window_nt: int = 150
    aav_payload_kb: float = 4.4
    #: Median TPM below which a gene counts as not transcribed. A convention
    #: rather than a measurement, which is why it is a threshold a laboratory
    #: can move rather than a constant in the code.
    expressed_above_tpm: float = 1.0
    #: Total points below which no mechanism is claimed at all.
    minimum_points: int = 2
    #: Total points, and margin over the runner-up, required for "probable".
    probable_points: int = 4
    probable_margin: int = 2


@dataclass(frozen=True, slots=True)
class Ruleset:
    """A versioned, digested set of rules."""

    version: str
    description: str
    thresholds: Thresholds
    rules: tuple[Rule, ...]
    digest: str

    @property
    def short_digest(self) -> str:
        return self.digest[:12]

    @property
    def pin(self) -> str:
        """The citation form that goes into a mechanism call."""
        return f"{self.version}@{self.short_digest}"


def load_ruleset(path: str | Path) -> Ruleset:
    """Read and validate a rule file."""
    raw_bytes = Path(path).read_bytes()
    digest = hashlib.sha256(raw_bytes).hexdigest()
    document = yaml.safe_load(raw_bytes)

    if not isinstance(document, dict):
        raise RulesetError(f"{path}: rule file must be a mapping at the top level")
    for required in ("version", "rules"):
        if required not in document:
            raise RulesetError(f"{path}: rule file has no {required!r}")

    thresholds = Thresholds(**(document.get("thresholds") or {}))

    rules: list[Rule] = []
    seen_ids: set[str] = set()
    for index, entry in enumerate(document["rules"], start=1):
        rule = _parse_rule(entry, index, path)
        if rule.id in seen_ids:
            raise RulesetError(f"{path}: duplicate rule id {rule.id!r}")
        seen_ids.add(rule.id)
        rules.append(rule)

    if not rules:
        raise RulesetError(f"{path}: rule file declares no rules")

    return Ruleset(
        version=str(document["version"]),
        description=str(document.get("description", "")),
        thresholds=thresholds,
        rules=tuple(rules),
        digest=digest,
    )


def _parse_rule(entry: Any, index: int, path: str | Path) -> Rule:
    if not isinstance(entry, dict):
        raise RulesetError(f"{path}: rule {index} is not a mapping")
    for required in ("id", "supports", "strength", "when", "because"):
        if required not in entry:
            raise RulesetError(f"{path}: rule {index} has no {required!r}")

    try:
        supports = Mechanism(entry["supports"])
    except ValueError as exc:
        raise RulesetError(
            f"{path}: rule {entry['id']!r} supports unknown mechanism {entry['supports']!r}"
        ) from exc
    try:
        strength = Strength(entry["strength"])
    except ValueError as exc:
        raise RulesetError(
            f"{path}: rule {entry['id']!r} has unknown strength {entry['strength']!r}"
        ) from exc

    return Rule(
        id=str(entry["id"]),
        supports=supports,
        strength=strength,
        when=entry["when"],
        because=" ".join(str(entry["because"]).split()),
        citation=str(entry.get("citation", "")),
    )


def evaluate_predicate(node: Any, features: Features, rule_id: str) -> bool:
    """Interpret one predicate node.

    Public because the modality rules in M6 and the off-target risk rules in M7
    are written in the same language and evaluated by the same interpreter. One
    predicate dialect across every rule file is worth more than tidy module
    boundaries: a curator learns it once.

    Combinators are ``all``, ``any`` and ``not``; a leaf names a feature and one
    operator. Anything else is an error rather than a false, because a rule that
    never fires because of a typo is indistinguishable from a rule that
    correctly never fires — and only one of those is a bug.
    """
    if not isinstance(node, dict):
        raise RulesetError(
            f"rule {rule_id!r}: predicate must be a mapping, got {type(node).__name__}"
        )

    if "all" in node:
        return all(evaluate_predicate(child, features, rule_id) for child in node["all"])
    if "any" in node:
        return any(evaluate_predicate(child, features, rule_id) for child in node["any"])
    if "not" in node:
        return not evaluate_predicate(node["not"], features, rule_id)

    if "feature" not in node:
        raise RulesetError(f"rule {rule_id!r}: predicate names no feature and is not a combinator")

    name = node["feature"]
    try:
        value = features.get(name)
    except KeyError as exc:
        raise RulesetError(
            f"rule {rule_id!r} tests unknown feature {name!r}; "
            f"known features are: {', '.join(features.known_names())}"
        ) from exc

    operators = [key for key in node if key not in {"feature"}]
    if len(operators) != 1:
        raise RulesetError(
            f"rule {rule_id!r}: feature {name!r} needs exactly one operator, "
            f"got {operators or 'none'}"
        )
    operator = operators[0]
    if operator not in _OPERATORS:
        raise RulesetError(
            f"rule {rule_id!r}: unknown operator {operator!r}; "
            f"available: {', '.join(sorted(_OPERATORS))}"
        )
    return bool(_OPERATORS[operator](value, node[operator]))
