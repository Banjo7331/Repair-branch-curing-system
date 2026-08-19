"""Turning a case into the flat record the rules are written against.

Every rule in ``rules/mechanism-*.yaml`` tests a named feature. This module is
the single place those names are produced, which means the rule file can be read
as a document — a geneticist checking ``gene.nulls_are_milder`` does not need to
know that it came from a dataclass field.

Two decisions worth stating. Features are computed once and are pure data: no
rule can reach past this record into an object and pull out something the record
does not name. And a feature that cannot be computed is ``None`` rather than a
default, so a rule that depends on missing data does not quietly fire on a zero.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from repairbench.context.expression import Tissue, median_tpm
from repairbench.model import Consequence, Gene, Zygosity
from repairbench.transcript import (
    NMD_BOUNDARY_NT,
    REINITIATION_WINDOW_NT,
    NMDOutcome,
    NMDPrediction,
    Transcript,
)


@dataclass(frozen=True, slots=True)
class SplicePrediction:
    """A splice-effect prediction from an external tool, kept as a raw score.

    The threshold is not applied here. SpliceAI's authors publish two of them —
    0.2 for recall, 0.5 for precision — and which one a laboratory uses is a
    policy choice, so it lives in the rule file where a reviewer can see it.
    """

    max_delta: float | None = None
    tool: str = "SpliceAI"


@dataclass(frozen=True, slots=True)
class Variant:
    """The variant as annotated against one transcript."""

    gene: str
    consequence: Consequence
    #: CDS offset of the variant. For a predicted-null variant this is where the
    #: premature stop lands, which is what the NMD arithmetic needs.
    cds_position: int
    protein_change: str = ""
    hgvs_c: str = ""
    #: How many copies the patient has left. Half the interventions in M6 depend
    #: on there being an unaffected allele, so this is not optional detail.
    zygosity: Zygosity = Zygosity.UNKNOWN


@dataclass(frozen=True, slots=True)
class MechanismQuery:
    """Everything the resolver is allowed to look at."""

    variant: Variant
    transcript: Transcript
    gene: Gene
    splice: SplicePrediction = field(default_factory=SplicePrediction)
    #: The tissue the disease affects. Case-scoped, like the phenotype: it is a
    #: fact about the patient rather than a release.
    tissue: Tissue | None = None
    #: Median TPM per tissue for this gene, from the expression source. Absent
    #: when no expression release was loaded, which the rules handle by not
    #: firing rather than by assuming silence.
    expression: dict[str, float] | None = None


@dataclass(frozen=True, slots=True)
class FeatureRecord:
    """A flat, immutable view of one case, keyed by the names rules use."""

    values: dict[str, Any]
    nmd: NMDPrediction

    def get(self, name: str) -> Any:
        """Look up a feature, refusing unknown names.

        A typo in a rule file must fail the run rather than evaluate to false —
        a rule that silently never fires is indistinguishable from a rule that
        correctly never fires, and only one of those is a bug.
        """
        if name not in self.values:
            raise KeyError(name)
        return self.values[name]

    def known_names(self) -> list[str]:
        return sorted(self.values)


def build_features(
    query: MechanismQuery,
    *,
    nmd_boundary_nt: int = NMD_BOUNDARY_NT,
    reinitiation_window_nt: int = REINITIATION_WINDOW_NT,
) -> FeatureRecord:
    """Compute every feature the rules may test.

    The two NMD distances are parameters rather than constants because they are
    rule values, not arithmetic: the junction boundary is quoted as 50–55 nt in
    the literature, and a laboratory that prefers 55 should be able to say so in
    the rule file rather than in a patch.
    """
    variant, transcript, gene = query.variant, query.transcript, query.gene

    if variant.consequence.is_predicted_null:
        nmd = transcript.predict_nmd_with(
            variant.cds_position,
            boundary_nt=nmd_boundary_nt,
            reinitiation_window_nt=reinitiation_window_nt,
        )
    else:
        nmd = NMDPrediction(
            NMDOutcome.NOT_APPLICABLE,
            f"{variant.consequence.value} does not introduce a premature termination codon",
        )

    variant_exon = transcript.exon_at(variant.cds_position)
    distribution = gene.distribution

    values: dict[str, Any] = {
        # --- what the variant does to the sequence
        "consequence": variant.consequence.value,
        "consequence.is_predicted_null": variant.consequence.is_predicted_null,
        "consequence.is_splice_affecting": variant.consequence.is_splice_affecting,
        "consequence.preserves_reading_frame": variant.consequence.preserves_reading_frame,
        "consequence.is_missense": variant.consequence is Consequence.MISSENSE,
        # --- what the transcript does with it
        "nmd.outcome": nmd.outcome.value,
        "nmd.certain": nmd.certain,
        "transcript.exon_count": transcript.exon_count,
        "transcript.cdna_kb": round(transcript.cdna_kilobases, 3),
        "transcript.mane_select": transcript.mane_select,
        "transcript.variant_exon": variant_exon,
        "transcript.variant_in_last_exon": variant_exon == transcript.exon_count,
        "transcript.variant_exon_preserves_frame": transcript.exon_preserves_frame(variant_exon),
        # --- what is known about the gene
        "gene.symbol": gene.symbol,
        "gene.haploinsufficiency": gene.haploinsufficiency.value,
        "gene.haploinsufficiency.supported": gene.haploinsufficiency.supports_haploinsufficiency,
        "gene.haploinsufficiency.refuted": gene.haploinsufficiency.refutes_haploinsufficiency,
        #: Whether anybody has looked, as distinct from what they found. A rule
        #: that wants to hedge on an uncurated gene can say so; none does yet,
        #: and the feature exists so that the choice is available in the rule
        #: file rather than requiring a code change.
        "gene.haploinsufficiency.curated": gene.haploinsufficiency.is_curated,
        "gene.triplosensitivity": gene.triplosensitivity.value,
        "gene.loeuf": gene.loeuf,
        "gene.forms_multimer": gene.forms_multimer,
        "gene.nulls_are_milder": gene.truncating_variants_are_milder,
        "gene.imprinted": gene.imprinting.is_imprinted,
        "gene.imprinting": gene.imprinting.value,
        "gene.silenced_allele_intact": gene.silenced_allele_intact,
        "gene.curated_mechanism": gene.curated_mechanism,
        "gene.has_curated_mechanism": gene.curated_mechanism is not None,
        # --- where pathogenic variation sits
        "gene.missense_total": distribution.pathogenic_missense_total,
        "gene.missense_in_hotspot": distribution.pathogenic_missense_in_hotspot,
        "gene.missense_clustering_ratio": round(distribution.clustering_ratio, 3),
        "gene.truncating_total": distribution.pathogenic_truncating_total,
        "gene.has_pathogenic_truncating": distribution.pathogenic_truncating_total > 0,
        # A share rather than a yes/no. The boolean above is kept because a rule
        # may legitimately want it, but nothing in mechanism-v1 reads it any
        # more: four truncating variants among eighty is what PIK3CA actually
        # looks like in ClinVar, and the boolean calls that "truncation causes
        # disease here".
        "gene.truncating_fraction": (
            None
            if distribution.truncating_fraction is None
            else round(distribution.truncating_fraction, 3)
        ),
        "gene.variants_counted": distribution.counted,
        # --- what the patient has left
        "variant.zygosity": variant.zygosity.value,
        "variant.has_wild_type_allele": variant.zygosity.leaves_a_wild_type_allele,
        "variant.zygosity_known": variant.zygosity is not Zygosity.UNKNOWN,
        # --- where the disease is, and whether the gene is on there
        "tissue.name": query.tissue.name if query.tissue else None,
        "tissue.system": query.tissue.system.value if query.tissue else None,
        "tissue.known": query.tissue is not None,
        "expression.tpm_in_affected_tissue": median_tpm(query.expression, query.tissue),
        "expression.measured": median_tpm(query.expression, query.tissue) is not None,
        # --- external predictions
        "splice.max_delta": query.splice.max_delta,
        "splice.tool": query.splice.tool,
    }
    return FeatureRecord(values=values, nmd=nmd)
