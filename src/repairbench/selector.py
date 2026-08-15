"""Applying the modality rules to a mechanism call.

Two structural decisions carry this module, and both are here rather than in the
rule file because they are not clinical claims — they are statements about what
this software is allowed to conclude.

**An unresolved mechanism blocks everything.** If M5 returned ``undetermined``,
no modality is assessed at all. The alternative is worse than useless: rules
would fire on the transcript facts alone and produce a plausible-looking list of
interventions resting on nothing.

**A contraindication outranks any number of indications.** Reasons to try
something accumulate; one reason not to ends the matter. This is the ordering
that stops the module's worst available mistake — proposing gene addition for a
dominant-negative variant, where a normal copy leaves every poisoning subunit
exactly where it was.
"""

from __future__ import annotations

from typing import Any

from repairbench.features import FeatureRecord, MechanismQuery, build_features
from repairbench.modality import (
    ModalityAssessment,
    ModalityEvidence,
    ModalitySelection,
    Stance,
    Verdict,
    all_modalities,
)
from repairbench.modality_rules import ModalityRuleset
from repairbench.model import Consequence, MechanismCall


def select(
    call: MechanismCall,
    query: MechanismQuery,
    ruleset: ModalityRuleset,
) -> ModalitySelection:
    """Assess every modality against one mechanism call."""
    if not call.mechanism.is_determined:
        return ModalitySelection(
            gene=call.gene,
            mechanism=call.mechanism,
            assessments=tuple(
                ModalityAssessment(modality, Verdict.BLOCKED_BY_UNRESOLVED_MECHANISM)
                for modality in all_modalities()
            ),
            ruleset_version=ruleset.pin,
            blocked_reason=(
                "the mechanism is unresolved, and every modality below depends on it — "
                "assessing them anyway would produce a list of interventions resting on "
                "transcript facts alone"
            ),
        )

    features = _features_for(call, query)

    fired = [
        ModalityEvidence(
            rule_id=rule.id,
            modality=rule.modality,
            stance=rule.stance,
            strength=rule.strength,
            because=rule.because,
            citation=rule.citation,
        )
        for rule in ruleset.rules
        if rule.fires(features)
    ]

    assessments = []
    for modality in all_modalities():
        indications = tuple(
            e for e in fired if e.modality is modality and e.stance is Stance.INDICATES
        )
        contraindications = tuple(
            e for e in fired if e.modality is modality and e.stance is Stance.CONTRAINDICATES
        )
        assessments.append(
            ModalityAssessment(
                modality=modality,
                verdict=_verdict(indications, contraindications, ruleset),
                indications=indications,
                contraindications=contraindications,
            )
        )

    return ModalitySelection(
        gene=call.gene,
        mechanism=call.mechanism,
        assessments=tuple(assessments),
        ruleset_version=ruleset.pin,
        caveats=_caveats(features, tuple(assessments), ruleset),
    )


def _caveats(
    features: FeatureRecord,
    assessments: tuple[ModalityAssessment, ...],
    ruleset: ModalityRuleset,
) -> tuple[str, ...]:
    """Things that do not rule a modality out but change how to read its verdict.

    Unknown zygosity is the case this exists for. Treating it as "no wild-type
    allele" would rule out real options on missing data; treating it as "yes"
    would offer options that may not exist. Neither is acceptable, so the
    verdict stands and the gap is stated.
    """
    caveats: list[str] = []

    tissue = features.get("tissue.name")
    if tissue is None:
        caveats.append(
            "no affected tissue was given, so nothing below was checked against where the "
            "gene is actually switched on — a variant in a gene silent in the affected "
            "tissue is unlikely to be the cause of anything"
        )
    elif any(assessment.verdict.is_available for assessment in assessments):
        caveats.append(
            f"delivery to {features.get('tissue.system')} is not assessed anywhere in this "
            f"package. Whether a vector, an oligonucleotide or an editor reaches {tissue} at "
            "a useful dose is the question that decides most of these in practice, and "
            "nothing here answers it"
        )

    if features.get("variant.has_wild_type_allele") is not None:
        return tuple(caveats)

    affected = [
        assessment.modality.value
        for assessment in assessments
        if assessment.verdict.is_available
        and assessment.modality in ruleset.require_wild_type_allele
    ]
    if affected:
        caveats.append(
            "zygosity was not supplied, and these modalities only work if the patient has an "
            f"unaffected copy of the gene: {', '.join(affected)}. Their verdicts above assume "
            "one exists and must not be read as established until it is confirmed."
        )
    return tuple(caveats)


def _verdict(
    indications: tuple[ModalityEvidence, ...],
    contraindications: tuple[ModalityEvidence, ...],
    ruleset: ModalityRuleset,
) -> Verdict:
    if contraindications:
        return Verdict.CONTRAINDICATED
    points = sum(e.strength.points for e in indications)
    if points >= ruleset.thresholds.minimum_indication_points:
        return Verdict.INDICATED
    return Verdict.NOT_INDICATED


def _features_for(call: MechanismCall, query: MechanismQuery) -> FeatureRecord:
    """The mechanism features, plus what M5 concluded.

    Recomputing the mechanism-level features rather than passing a separate
    record keeps one vocabulary across both rule files: a curator who knows
    ``gene.forms_multimer`` from the mechanism rules can use it in the modality
    rules without learning a second set of names.
    """
    base = build_features(query)
    feasibility = call.feasibility

    extra: dict[str, Any] = {
        "mechanism": call.mechanism.value,
        "mechanism.confidence": call.confidence.value,
        "mechanism.needs_review": call.needs_review,
        "mechanism.tolerates_gene_addition": call.mechanism.tolerates_gene_addition,
        "feasibility.gene_addition_coherent": feasibility.gene_addition_coherent,
        "feasibility.fits_viral_payload": feasibility.fits_viral_payload,
        "feasibility.exon_skipping_preserves_frame": feasibility.exon_skipping_preserves_frame,
        "feasibility.silenced_allele_available": feasibility.silenced_allele_available,
        "feasibility.allele_specific_silencing_indicated": (
            feasibility.allele_specific_silencing_indicated
        ),
        # Whether a mutant transcript can be told apart from the normal one at
        # all. A single-nucleotide change gives an antisense oligonucleotide
        # something to bind that the wild-type allele does not have; a whole-gene
        # deletion gives it nothing.
        "variant.provides_allele_specific_target": query.variant.consequence
        in _ALLELE_DISTINGUISHING,
    }
    return FeatureRecord(values={**base.values, **extra}, nmd=base.nmd)


#: Consequences that leave a sequence difference an oligonucleotide could bind.
#: Whole-gene and multi-exon deletions do not appear here, and they are the case
#: where allele-specific silencing has nothing to aim at.
_ALLELE_DISTINGUISHING = frozenset(
    {
        Consequence.MISSENSE,
        Consequence.NONSENSE,
        Consequence.FRAMESHIFT,
        Consequence.INFRAME_DELETION,
        Consequence.INFRAME_INSERTION,
        Consequence.SPLICE_ACCEPTOR,
        Consequence.SPLICE_DONOR,
        Consequence.SPLICE_REGION,
        Consequence.START_LOST,
        Consequence.STOP_LOST,
    }
)
