"""Applying the rules and adding up what they say.

The engine is deliberately dull. Every clinical judgement lives in the rule
file; what is left here is arithmetic over the rules that fired, plus one
structural decision worth defending.

**Curation outranks inference.** A rule marked ``decisive`` — in practice, the
one that reads an expert determination from ClinGen — ends the matter. But the
other rules are still evaluated, and any that disagree are reported as conflicts
rather than discarded. A curated loss-of-function gene whose variant escapes
nonsense-mediated decay is exactly the case where a reviewer needs to see both
facts, because the curation describes the gene and the arithmetic describes this
variant, and the second can be an exception to the first.
"""

from __future__ import annotations

from repairbench.features import FeatureRecord, MechanismQuery, build_features
from repairbench.model import (
    Confidence,
    Evidence,
    Feasibility,
    Mechanism,
    MechanismCall,
    Strength,
)
from repairbench.ruleset import Ruleset


def resolve(query: MechanismQuery, ruleset: Ruleset) -> MechanismCall:
    """Determine the mechanism of one variant."""
    features = build_features(
        query,
        nmd_boundary_nt=ruleset.thresholds.nmd_boundary_nt,
        reinitiation_window_nt=ruleset.thresholds.reinitiation_window_nt,
    )

    fired = [
        Evidence(
            rule_id=rule.id,
            supports=rule.supports,
            strength=rule.strength,
            because=rule.because,
            citation=rule.citation,
        )
        for rule in ruleset.rules
        if rule.fires(features)
    ]

    mechanism, confidence = _aggregate(fired, ruleset)
    supporting = tuple(e for e in fired if e.supports is mechanism)
    # A caution argues that the evidence does not settle the question. When the
    # answer is "undetermined" it *is* the answer, and when the answer is
    # anything else it is the reason a reviewer should look twice — so it is
    # never simply dropped.
    conflicts = tuple(e for e in fired if e.supports is not mechanism)

    return MechanismCall(
        gene=query.gene.symbol,
        transcript=query.transcript.accession,
        mechanism=mechanism,
        confidence=confidence,
        evidence=supporting,
        feasibility=_feasibility(query, features, mechanism, ruleset),
        conflicts=conflicts,
        ruleset_version=ruleset.pin,
    )


def _aggregate(fired: list[Evidence], ruleset: Ruleset) -> tuple[Mechanism, Confidence]:
    """Turn the rules that fired into one mechanism and a confidence.

    Three kinds of rule meet here. A *decisive* rule ends the matter. Rules that
    support a mechanism add up. And *cautions* — rules whose ``supports`` is
    ``undetermined`` — argue that the evidence does not settle the question; a
    strong one caps confidence at "possible", which is how the rule file says
    "not so fast" without having to say "no".
    """
    for evidence in fired:
        if evidence.strength is Strength.DECISIVE:
            return evidence.supports, Confidence.ESTABLISHED

    totals: dict[Mechanism, int] = {}
    strong_caution = False
    for evidence in fired:
        if not evidence.supports.is_determined:
            strong_caution = strong_caution or evidence.strength is Strength.STRONG
            continue
        totals[evidence.supports] = totals.get(evidence.supports, 0) + evidence.strength.points

    if not totals:
        return Mechanism.UNDETERMINED, Confidence.NONE

    ranked = sorted(totals.items(), key=lambda item: item[1], reverse=True)
    winner, winning_points = ranked[0]
    runner_up_points = ranked[1][1] if len(ranked) > 1 else 0
    margin = winning_points - runner_up_points

    thresholds = ruleset.thresholds
    if winning_points < thresholds.minimum_points:
        # Something fired, but not enough to point anywhere. Reporting the
        # mechanism anyway would launder one supporting rule into an answer.
        return Mechanism.UNDETERMINED, Confidence.NONE

    reaches_probable = (
        winning_points >= thresholds.probable_points and margin >= thresholds.probable_margin
    )
    if reaches_probable and not strong_caution:
        return winner, Confidence.PROBABLE
    return winner, Confidence.POSSIBLE


def _feasibility(
    query: MechanismQuery,
    features: FeatureRecord,
    mechanism: Mechanism,
    ruleset: Ruleset,
) -> Feasibility:
    """What this mechanism and this transcript make possible.

    Derivation, not judgement — which is why it is in Python rather than in the
    rule file. Each flag follows from a fact already established; none of them
    is a clinical claim about whether an intervention would work.
    """
    notes: list[str] = []
    transcript, gene = query.transcript, query.gene

    fits_payload = transcript.cdna_kilobases <= ruleset.thresholds.aav_payload_kb
    if not fits_payload:
        notes.append(
            f"coding sequence is {transcript.cdna_kilobases:.1f} kb, above the "
            f"{ruleset.thresholds.aav_payload_kb} kb practical AAV ceiling — a shortened "
            "construct would be protein engineering, not gene replacement"
        )

    exon_frame: bool | None = None
    if features.get("consequence.is_splice_affecting") or features.get(
        "consequence.is_predicted_null"
    ):
        exon_frame = features.get("transcript.variant_exon_preserves_frame")
        if exon_frame is False:
            notes.append(
                f"exon {features.get('transcript.variant_exon')} is not a multiple of three, "
                "so skipping it would shift the reading frame rather than rescue it"
            )

    silenced_available = gene.imprinting.is_imprinted and gene.silenced_allele_intact
    if gene.imprinting.is_imprinted and not gene.silenced_allele_intact:
        notes.append(
            "the gene is imprinted but the silenced allele is not recorded as intact — "
            "reactivation has nothing to reactivate unless that is established"
        )

    silencing_indicated = mechanism in {Mechanism.GAIN_OF_FUNCTION, Mechanism.DOMINANT_NEGATIVE}
    if silencing_indicated:
        notes.append(
            "the mutant allele is doing harm rather than being absent, so supplying a normal "
            "copy does not address the mechanism and may worsen it"
        )
    if mechanism is Mechanism.UNDETERMINED:
        notes.append(
            "no feasibility flag below should be acted on: the mechanism is unresolved, and "
            "every downstream modality depends on it"
        )

    return Feasibility(
        gene_addition_coherent=mechanism.tolerates_gene_addition,
        fits_viral_payload=fits_payload,
        exon_skipping_preserves_frame=exon_frame,
        silenced_allele_available=silenced_available,
        allele_specific_silencing_indicated=silencing_indicated,
        notes=tuple(notes),
    )
