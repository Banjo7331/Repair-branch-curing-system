"""Loading the modality rule file.

Same predicate dialect as the mechanism rules — a curator learns it once — but a
different conclusion shape: a modality rule names an intervention class and
takes a *stance* on it rather than arguing for one answer among several.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from repairbench.features import FeatureRecord
from repairbench.modality import Modality, Stance
from repairbench.model import Strength
from repairbench.ruleset import RulesetError, evaluate_predicate


@dataclass(frozen=True, slots=True)
class ModalityRule:
    """One claim about whether a class of intervention fits a mechanism."""

    id: str
    modality: Modality
    stance: Stance
    strength: Strength
    when: dict[str, Any]
    because: str
    citation: str = ""

    def fires(self, features: FeatureRecord) -> bool:
        return evaluate_predicate(self.when, features, self.id)


@dataclass(frozen=True, slots=True)
class ModalityThresholds:
    """Indication points required before a modality is called indicated.

    Set above one on purpose: a single supporting rule is a hint, and this module
    exists downstream of a mechanism call that may itself be only *possible*.
    """

    minimum_indication_points: int = 2


@dataclass(frozen=True, slots=True)
class ModalityRuleset:
    version: str
    description: str
    thresholds: ModalityThresholds
    rules: tuple[ModalityRule, ...]
    digest: str
    #: Modalities that only work if the patient has an unaffected copy. Declared
    #: in the rule file rather than hard-coded, because which interventions
    #: depend on a wild-type allele is a clinical claim like any other — and it
    #: is what the engine consults when zygosity is unknown.
    require_wild_type_allele: frozenset[Modality] = frozenset()

    @property
    def short_digest(self) -> str:
        return self.digest[:12]

    @property
    def pin(self) -> str:
        return f"{self.version}@{self.short_digest}"


def load_modality_ruleset(path: str | Path) -> ModalityRuleset:
    """Read and validate a modality rule file."""
    raw_bytes = Path(path).read_bytes()
    document = yaml.safe_load(raw_bytes)

    if not isinstance(document, dict):
        raise RulesetError(f"{path}: modality rule file must be a mapping")
    for required in ("version", "rules"):
        if required not in document:
            raise RulesetError(f"{path}: modality rule file has no {required!r}")

    rules: list[ModalityRule] = []
    seen: set[str] = set()
    for index, entry in enumerate(document["rules"], start=1):
        rule = _parse(entry, index, path)
        if rule.id in seen:
            raise RulesetError(f"{path}: duplicate rule id {rule.id!r}")
        seen.add(rule.id)
        rules.append(rule)

    if not rules:
        raise RulesetError(f"{path}: modality rule file declares no rules")

    requires: set[Modality] = set()
    for name in document.get("modalities_requiring_wild_type_allele", []):
        try:
            requires.add(Modality(name))
        except ValueError as exc:
            raise RulesetError(
                f"{path}: modalities_requiring_wild_type_allele names unknown modality {name!r}"
            ) from exc

    return ModalityRuleset(
        version=str(document["version"]),
        description=str(document.get("description", "")),
        thresholds=ModalityThresholds(**(document.get("thresholds") or {})),
        rules=tuple(rules),
        digest=hashlib.sha256(raw_bytes).hexdigest(),
        require_wild_type_allele=frozenset(requires),
    )


def _parse(entry: Any, index: int, path: str | Path) -> ModalityRule:
    if not isinstance(entry, dict):
        raise RulesetError(f"{path}: rule {index} is not a mapping")
    for required in ("id", "modality", "stance", "strength", "when", "because"):
        if required not in entry:
            raise RulesetError(f"{path}: rule {index} has no {required!r}")

    try:
        modality = Modality(entry["modality"])
    except ValueError as exc:
        raise RulesetError(
            f"{path}: rule {entry['id']!r} names unknown modality {entry['modality']!r}"
        ) from exc
    try:
        stance = Stance(entry["stance"])
    except ValueError as exc:
        raise RulesetError(
            f"{path}: rule {entry['id']!r} takes unknown stance {entry['stance']!r}"
        ) from exc
    try:
        strength = Strength(entry["strength"])
    except ValueError as exc:
        raise RulesetError(
            f"{path}: rule {entry['id']!r} has unknown strength {entry['strength']!r}"
        ) from exc
    if strength is Strength.DECISIVE:
        raise RulesetError(
            f"{path}: rule {entry['id']!r} is marked decisive. Modality rules do not "
            "have a decisive tier — a contraindication already outranks every "
            "indication, and nothing here should be able to end the discussion "
            "in the other direction."
        )

    return ModalityRule(
        id=str(entry["id"]),
        modality=modality,
        stance=stance,
        strength=strength,
        when=entry["when"],
        because=" ".join(str(entry["because"]).split()),
        citation=str(entry.get("citation", "")),
    )
