"""Ranking off-target hits by where they land, not by how many mismatches they have.

This is the part the search tools leave to you, and the part that decides
whether a hit list is usable. Cas-OFFinder returns every site within *n*
mismatches, sorted by mismatch count if it sorts at all — which puts a
four-mismatch hit in the coding sequence of a tumour suppressor below a
two-mismatch hit in an intergenic desert. The ranking is arithmetically correct
and clinically backwards.

What makes a hit matter is where it lands:

* Is it in coding sequence, or in an intron, or nowhere near a gene?
* Is the gene one a cell cannot lose, or one whose disruption drives cancer?
* Is the gene even transcribed in the tissue this therapy is aimed at? A hit in
  a gene switched off there is a different finding from the same hit in the one
  tissue where the gene does its job — and the expression data for asking is
  already in this package.

The tiers and the reasons for them live in ``rules/offtarget-v1.yaml``, in the
same predicate language as the mechanism and modality rules, for the same
reason: which off-target hits are unacceptable is a clinical judgement, and a
clinical judgement belongs in a file a reviewer can read.

The sequence-level score is a separate matter and is deliberately absent. CFD
weights every mismatch by its position and identity, from a published table this
package does not carry — so a hit here is ranked by context, the mismatch count
is reported as the search gave it, and the report says which of the two you are
looking at.
"""

from __future__ import annotations

import hashlib
from collections.abc import Callable
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any, TypeAlias

import yaml

from repairbench.annotation.store import Placement
from repairbench.context.expression import Tissue, median_tpm
from repairbench.context.genelists import GeneList, GeneLists
from repairbench.design.offtarget import OffTargetHit
from repairbench.ruleset import RulesetError, evaluate_predicate

#: Anything that can say what is at a coordinate — ``TranscriptStore.locate``
#: in practice, and a lambda in a test.
Locator: TypeAlias = Callable[[str, int], Placement]


class RiskTier(StrEnum):
    """How much a hit matters, worst first."""

    PROHIBITIVE = "prohibitive"
    """Would stop this guide being used. Coding sequence of a gene the patient
    cannot afford to lose."""

    SERIOUS = "serious"
    """Needs answering before anything is ordered — usually by sequencing the
    site, which is a bench task rather than a software one."""

    MODERATE = "moderate"
    """Worth recording and worth watching."""

    LOW = "low"
    """Nothing here says it is safe. It says nothing in the rule file fired."""

    UNASSESSED = "unassessed"
    """No annotation covered this coordinate, so nothing was asked about it. The
    tier exists so that "we did not look" cannot be read as "we looked and it
    was fine"."""

    @classmethod
    def worst_first(cls) -> tuple[RiskTier, ...]:
        return tuple(cls)

    @property
    def rank(self) -> int:
        return list(RiskTier).index(self)


@dataclass(frozen=True, slots=True)
class RiskRule:
    id: str
    tier: RiskTier
    when: dict[str, Any]
    because: str
    citation: str = ""

    def fires(self, features: HitFeatures) -> bool:
        return evaluate_predicate(self.when, features, self.id)


@dataclass(frozen=True, slots=True)
class RiskRuleset:
    version: str
    description: str
    rules: tuple[RiskRule, ...]
    digest: str

    @property
    def pin(self) -> str:
        return f"{self.version}@{self.digest[:12]}"


@dataclass(frozen=True, slots=True)
class HitFeatures:
    """One hit, flattened into the names the rules use."""

    values: dict[str, Any]

    def get(self, name: str) -> Any:
        if name not in self.values:
            raise KeyError(name)
        return self.values[name]

    def known_names(self) -> list[str]:
        return sorted(self.values)


@dataclass(frozen=True, slots=True)
class RuleFired:
    rule_id: str
    tier: RiskTier
    because: str
    citation: str = ""


@dataclass(frozen=True, slots=True)
class RankedHit:
    hit: OffTargetHit
    tier: RiskTier
    where: str
    evidence: tuple[RuleFired, ...] = ()


@dataclass(frozen=True, slots=True)
class OffTargetAssessment:
    guide: str
    hits: tuple[RankedHit, ...]
    ruleset_pin: str
    source_pin: str
    scoring: str = ""
    unranked: str = ""

    def in_tier(self, tier: RiskTier) -> tuple[RankedHit, ...]:
        return tuple(hit for hit in self.hits if hit.tier is tier)

    @property
    def worst(self) -> RiskTier:
        tiers = [hit.tier for hit in self.hits]
        return min(tiers, key=lambda tier: tier.rank, default=RiskTier.LOW)


def build_hit_features(
    hit: OffTargetHit,
    placement: Placement | None,
    *,
    lists: GeneLists | None = None,
    expression: dict[str, float] | None = None,
    tissue: Tissue | None = None,
) -> HitFeatures:
    """Everything the risk rules may test about one hit.

    A feature nobody could compute is ``None``, never a default. "Not on the
    essential list" and "no list was supplied" are different claims, and a rule
    that fired on the second because it was written for the first would rank a
    hit as safe on the strength of a missing file.
    """
    gene = placement.gene if placement else None
    memberships = lists.lists_for(gene) if lists and gene else ()
    tpm = median_tpm(expression, tissue)

    values: dict[str, Any] = {
        # --- what the search found
        "hit.mismatches": hit.mismatches,
        "hit.bulge_size": hit.bulge_size,
        "hit.has_bulge": hit.has_bulge,
        # --- where it landed
        "hit.annotated": placement is not None,
        "hit.gene": gene,
        "hit.in_a_gene": gene is not None,
        "hit.in_coding_sequence": placement.in_coding_sequence if placement else None,
        "hit.in_transcript_span": placement.in_transcript_span if placement else None,
        # --- what kind of gene it is
        "hit.gene_lists_known": lists is not None,
        "hit.gene_is_essential": (GeneList.ESSENTIAL in memberships) if lists and gene else None,
        "hit.gene_is_oncogene": (GeneList.ONCOGENE in memberships) if lists and gene else None,
        "hit.gene_is_tumour_suppressor": (
            (GeneList.TUMOUR_SUPPRESSOR in memberships) if lists and gene else None
        ),
        # --- whether the gene is even on where the therapy is aimed
        "hit.tpm_in_target_tissue": tpm,
        "hit.expression_measured": tpm is not None,
        "tissue.name": tissue.name if tissue else None,
    }
    return HitFeatures(values=values)


def assess(
    hits: tuple[OffTargetHit, ...],
    ruleset: RiskRuleset,
    *,
    source_pin: str,
    locate: Locator | None = None,
    lists: GeneLists | None = None,
    expression: dict[str, dict[str, float]] | None = None,
    tissue: Tissue | None = None,
) -> OffTargetAssessment:
    """Rank every hit, worst first, with the rule that put it there.

    ``locate`` is any callable taking ``(chromosome, position)`` and returning a
    ``Placement`` — ``TranscriptStore.locate`` in practice. Without one, every
    hit comes back ``unassessed``, which is the honest answer: a hit list with no
    annotation behind it has not been reviewed, it has only been read.
    """
    ranked: list[RankedHit] = []

    for hit in hits:
        placement = locate(hit.chromosome, hit.position) if locate else None
        gene_expression = (
            expression.get(placement.gene) if expression and placement and placement.gene else None
        )
        features = build_hit_features(
            hit, placement, lists=lists, expression=gene_expression, tissue=tissue
        )

        fired = tuple(
            RuleFired(rule.id, rule.tier, rule.because, rule.citation)
            for rule in ruleset.rules
            if rule.fires(features)
        )
        tier = (
            min((rule.tier for rule in fired), key=lambda entry: entry.rank)
            if fired
            else (RiskTier.LOW if placement is not None else RiskTier.UNASSESSED)
        )
        ranked.append(
            RankedHit(
                hit=hit,
                tier=tier,
                where=placement.describe if placement else "no annotation was supplied",
                evidence=fired,
            )
        )

    ranked.sort(key=lambda entry: (entry.tier.rank, -entry.hit.mismatches))

    return OffTargetAssessment(
        guide=hits[0].guide if hits else "",
        hits=tuple(ranked),
        ruleset_pin=ruleset.pin,
        source_pin=source_pin,
        scoring=(
            "hits are ranked by where they land, not by a sequence score. No CFD table is "
            "attached, so the mismatch count is reported as the search gave it and is not "
            "weighted by position or by which base swapped for which"
        ),
        unranked=(
            ""
            if locate
            else (
                "no annotation was supplied, so no hit was placed in a gene and none of the "
                "context rules could be asked. This is a hit list that has been read, not one "
                "that has been reviewed"
            )
        ),
    )


def load_risk_rules(path: str | Path) -> RiskRuleset:
    """Read and validate the off-target risk rules."""
    raw = Path(path).read_bytes()
    document = yaml.safe_load(raw)

    if not isinstance(document, dict):
        raise RulesetError(f"{path}: risk rule file must be a mapping")
    for required in ("version", "rules"):
        if required not in document:
            raise RulesetError(f"{path}: risk rule file has no {required!r}")

    rules: list[RiskRule] = []
    seen: set[str] = set()
    for index, entry in enumerate(document["rules"], start=1):
        if not isinstance(entry, dict):
            raise RulesetError(f"{path}: rule {index} is not a mapping")
        for required in ("id", "tier", "when", "because"):
            if required not in entry:
                raise RulesetError(f"{path}: rule {index} has no {required!r}")
        try:
            tier = RiskTier(entry["tier"])
        except ValueError as error:
            raise RulesetError(
                f"{path}: rule {entry['id']!r} has unknown tier {entry['tier']!r}; "
                f"known tiers are {', '.join(tier.value for tier in RiskTier)}"
            ) from error
        if entry["id"] in seen:
            raise RulesetError(f"{path}: duplicate rule id {entry['id']!r}")
        seen.add(str(entry["id"]))
        rules.append(
            RiskRule(
                id=str(entry["id"]),
                tier=tier,
                when=entry["when"],
                because=" ".join(str(entry["because"]).split()),
                citation=str(entry.get("citation", "")),
            )
        )

    if not rules:
        raise RulesetError(f"{path}: risk rule file declares no rules")

    return RiskRuleset(
        version=str(document["version"]),
        description=str(document.get("description", "")),
        rules=tuple(rules),
        digest=hashlib.sha256(raw).hexdigest(),
    )
