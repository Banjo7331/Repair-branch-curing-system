"""ClinVar submissions, reduced to the two things this package asks of them.

**Where pathogenic variation sits in a gene**, which is the input to the
clustering rule — the one that separates gain of function from loss of function
by noticing that pathogenic missense variants pile into one stretch of protein
while truncating variants are absent. Those counts were typed in by hand until
this module existed, which meant the most inferential rule in the package was
resting on numbers this project made up.

**Where one variant actually is**, so a reference case can name a position
somebody reported rather than a position that seemed plausible.

Three refusals shape everything below.

**Review status is not decoration.** A submission from one laboratory with no
assertion criteria and a reviewed expert-panel classification are not the same
evidence, and averaging them is how a count becomes confident nonsense. The
minimum accepted level is a threshold in the rule file; the star rating is kept
per variant so a report can say what it counted.

**Conflicting is not pathogenic.** ClinVar's own aggregate says when submitters
disagree, and a parser that matched the substring "Pathogenic" would count
"Conflicting classifications of pathogenicity" as support.

**A hotspot here means a dense window, not a domain.** Clustering is computed
from protein positions with a window width declared in the rule file, and named
for what it is. A curated functional domain is a different claim and this file
has no access to one.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass
from enum import StrEnum

from repairbench.context.source import ContextError, Fact, Provenance, Source, read_tsv
from repairbench.model import MissenseDistribution

GENE_COLUMN = "GeneSymbol"
SIGNIFICANCE_COLUMN = "ClinicalSignificance"
REVIEW_COLUMN = "ReviewStatus"
ASSEMBLY_COLUMN = "Assembly"
CHROMOSOME_COLUMN = "Chromosome"
START_COLUMN = "PositionVCF"
REFERENCE_COLUMN = "ReferenceAlleleVCF"
ALTERNATE_COLUMN = "AlternateAlleleVCF"
NAME_COLUMN = "Name"
TYPE_COLUMN = "Type"

REQUIRED = {
    GENE_COLUMN,
    SIGNIFICANCE_COLUMN,
    REVIEW_COLUMN,
    ASSEMBLY_COLUMN,
    CHROMOSOME_COLUMN,
    NAME_COLUMN,
}

#: ``NM_000088.4(COL1A1):c.2461G>A (p.Gly821Ser)`` — the protein change is what
#: says whether a variant truncates, and the c. part is what a case cites.
_HGVS = re.compile(
    r"(?P<transcript>N[MR]_\d+\.\d+)"
    r"(?:\((?P<gene>[^)]+)\))?"
    r":(?P<coding>c\.[^ ]+)"
    r"(?:\s+\((?P<protein>p\.[^)]+)\))?"
)
#: A three-letter amino acid, its position, and what replaces it.
_PROTEIN = re.compile(r"p\.(?P<from>[A-Z][a-z]{2})(?P<position>\d+)(?P<to>[A-Z][a-z]{2}|Ter|=)")

#: Aggregate classifications this package counts as pathogenic. Written out
#: rather than matched as a substring, because "Conflicting classifications of
#: pathogenicity" contains the word and means the opposite of agreement.
PATHOGENIC = frozenset(
    {
        "Pathogenic",
        "Likely pathogenic",
        "Pathogenic/Likely pathogenic",
    }
)


class ReviewStatus(StrEnum):
    """ClinVar's star ratings, in its own words.

    Kept as a scale rather than a boolean because the threshold is a policy: a
    laboratory building a variant list may want two stars, and one counting
    variation to find a hotspot may accept one.
    """

    PRACTICE_GUIDELINE = "practice guideline"
    REVIEWED_BY_EXPERT_PANEL = "reviewed by expert panel"
    MULTIPLE_SUBMITTERS_NO_CONFLICTS = "criteria provided, multiple submitters, no conflicts"
    SINGLE_SUBMITTER = "criteria provided, single submitter"
    CONFLICTING = "criteria provided, conflicting classifications"
    NO_CRITERIA = "no assertion criteria provided"
    UNKNOWN = "unknown"

    @property
    def stars(self) -> int:
        return {
            ReviewStatus.PRACTICE_GUIDELINE: 4,
            ReviewStatus.REVIEWED_BY_EXPERT_PANEL: 3,
            ReviewStatus.MULTIPLE_SUBMITTERS_NO_CONFLICTS: 2,
            ReviewStatus.SINGLE_SUBMITTER: 1,
        }.get(self, 0)

    @classmethod
    def parse(cls, raw: str) -> ReviewStatus:
        """Map a review-status string, tolerating wording drift.

        Unknown wording becomes ``UNKNOWN`` with zero stars rather than raising:
        ClinVar rephrases these, and a release that renamed one should reduce
        what is counted rather than stop the run. What it must not do is score
        an unrecognised status as though it were reviewed.
        """
        text = raw.strip().lower()
        for status in cls:
            if status is not cls.UNKNOWN and status.value == text:
                return status
        return cls.UNKNOWN


class VariantKind(StrEnum):
    """What the protein change says the variant does."""

    MISSENSE = "missense"
    TRUNCATING = "truncating"
    SYNONYMOUS = "synonymous"
    OTHER = "other"


@dataclass(frozen=True, slots=True)
class ClinVarVariant:
    """One pathogenic submission, reduced to what this package reads."""

    gene: str
    kind: VariantKind
    review: ReviewStatus
    transcript: str = ""
    coding: str = ""
    protein: str = ""
    #: Protein position, when the name carries one. The clustering arithmetic
    #: needs it and roughly a fifth of records do not have it.
    protein_position: int | None = None
    chromosome: str = ""
    position: int | None = None
    reference: str = ""
    alternate: str = ""

    @property
    def is_placed(self) -> bool:
        return self.position is not None and bool(self.reference and self.alternate)


def classify(protein: str, variant_type: str) -> VariantKind:
    """What kind of change this is, from the protein consequence.

    Read from the HGVS protein name rather than from ClinVar's ``Type`` column,
    because ``Type`` says how the *sequence* changed — deletion, duplication,
    single nucleotide variant — and the rules ask what happened to the
    *product*. A single nucleotide variant can be missense or nonsense, and
    those argue for opposite mechanisms.
    """
    if not protein:
        return VariantKind.TRUNCATING if "deletion" in variant_type.lower() else VariantKind.OTHER
    if "Ter" in protein or "*" in protein or "fs" in protein:
        return VariantKind.TRUNCATING
    if protein.endswith("=") or "Sil" in protein:
        return VariantKind.SYNONYMOUS
    return VariantKind.MISSENSE if _PROTEIN.search(protein) else VariantKind.OTHER


def read_variants(
    source: Source,
    *,
    genes: set[str],
    assembly: str = "GRCh38",
    minimum_stars: int = 1,
) -> list[ClinVarVariant]:
    """Every pathogenic variant in these genes, at or above a review threshold.

    ``genes`` is required rather than optional. The file has millions of rows
    and a call without a filter would read all of them to answer a question
    about a handful — and, worse, would look like it was working.
    """
    if not genes:
        raise ContextError(
            "read_variants needs the genes to look for: the file is millions of rows, and a "
            "call without a filter would read all of them and look like it was working"
        )

    found: list[ClinVarVariant] = []
    for row in read_tsv(source, required=REQUIRED):
        if row[ASSEMBLY_COLUMN] != assembly:
            continue
        gene = row[GENE_COLUMN]
        if gene not in genes:
            continue
        if row[SIGNIFICANCE_COLUMN] not in PATHOGENIC:
            continue
        review = ReviewStatus.parse(row[REVIEW_COLUMN])
        if review.stars < minimum_stars:
            continue

        match = _HGVS.search(row[NAME_COLUMN])
        protein = match.group("protein") or "" if match else ""
        protein_match = _PROTEIN.search(protein)

        found.append(
            ClinVarVariant(
                gene=gene,
                kind=classify(protein, row.get(TYPE_COLUMN, "")),
                review=review,
                transcript=match.group("transcript") if match else "",
                coding=match.group("coding") if match else "",
                protein=protein,
                protein_position=int(protein_match.group("position")) if protein_match else None,
                chromosome=row[CHROMOSOME_COLUMN],
                position=_maybe_int(row.get(START_COLUMN, "")),
                reference=row.get(REFERENCE_COLUMN, ""),
                alternate=row.get(ALTERNATE_COLUMN, ""),
            )
        )
    return found


def _maybe_int(raw: str) -> int | None:
    try:
        return int(raw)
    except (TypeError, ValueError):
        return None


def distribution_for(
    variants: list[ClinVarVariant], *, hotspot_window_aa: int = 20
) -> MissenseDistribution:
    """Counts of pathogenic variation, and how tightly the missense clusters.

    The hotspot count is the number of missense variants inside the densest
    window of ``hotspot_window_aa`` residues. That is a measurement of
    clustering and not a claim about a functional domain — the rule that reads
    it says so, and the window is a rule-file threshold because where the line
    falls is a judgement rather than a fact.
    """
    missense = [v for v in variants if v.kind is VariantKind.MISSENSE]
    truncating = [v for v in variants if v.kind is VariantKind.TRUNCATING]
    positions = sorted(v.protein_position for v in missense if v.protein_position is not None)

    return MissenseDistribution(
        pathogenic_missense_total=len(missense),
        pathogenic_missense_in_hotspot=_densest_window(positions, hotspot_window_aa),
        pathogenic_truncating_total=len(truncating),
    )


def _densest_window(positions: list[int], width: int) -> int:
    """How many of these positions fall inside the busiest window of that width.

    A two-pointer sweep over sorted positions: for each start, advance the end
    while it stays inside the window. Linear, and exact — approximating it by
    binning would make the answer depend on where the bin edges happened to
    fall, which for a gene whose hotspot straddles one would understate the
    clustering that matters most.
    """
    if not positions:
        return 0
    best = start = 0
    for end, position in enumerate(positions):
        while position - positions[start] >= width:
            start += 1
        best = max(best, end - start + 1)
    return best


def review_summary(variants: list[ClinVarVariant]) -> str:
    """What was counted, by review status — for the provenance line.

    A distribution assembled mostly from single-submitter records and one
    assembled from expert-panel records support the same rule to very different
    degrees, and the count alone hides which one this is.
    """
    counts = Counter(variant.review for variant in variants)
    parts = [
        f"{count}×{status.stars}★"
        for status, count in sorted(counts.items(), key=lambda item: -item[0].stars)
    ]
    return ", ".join(parts) or "nothing counted"


def group_by_gene(variants: list[ClinVarVariant]) -> dict[str, list[ClinVarVariant]]:
    """The same variants, filed under the gene each one names."""
    grouped: dict[str, list[ClinVarVariant]] = {}
    for variant in variants:
        grouped.setdefault(variant.gene, []).append(variant)
    return grouped


def commonest_transcript(variants: list[ClinVarVariant]) -> str:
    """Which transcript these submissions are mostly written against.

    Not a claim about which transcript is *right* — that question belongs to the
    annotation store, which answers it from MANE and curation status. This is
    only what the submitters used, and it is worth knowing because a c. position
    read off one transcript and asserted against another is off by the length of
    whatever exons differ.
    """
    counts = Counter(variant.transcript for variant in variants if variant.transcript)
    return counts.most_common(1)[0][0] if counts else ""


def exemplars(
    variants: list[ClinVarVariant],
    kind: VariantKind,
    *,
    transcript: str = "",
    limit: int = 5,
) -> list[ClinVarVariant]:
    """Reported variants of one kind, best-reviewed first.

    For the one job of letting a case cite a position somebody actually
    submitted instead of a position that seemed plausible. Restricted to a
    single transcript when one is given, because c. numbering is only meaningful
    against the transcript it was written on.
    """
    pool = [
        variant
        for variant in variants
        if variant.kind is kind
        and variant.coding
        and (not transcript or variant.transcript == transcript)
    ]
    pool.sort(key=lambda variant: (-variant.review.stars, variant.protein_position or 0))
    return pool[:limit]


def ingest(
    source: Source,
    into: dict[str, Provenance],
    *,
    genes: set[str],
    minimum_stars: int = 1,
    hotspot_window_aa: int = 20,
) -> int:
    """Record each gene's distribution of pathogenic variation."""
    variants = read_variants(source, genes=genes, minimum_stars=minimum_stars)
    by_gene = group_by_gene(variants)

    for gene, found in by_gene.items():
        into.setdefault(gene, Provenance(gene=gene)).record(
            Fact(
                field="distribution",
                value=distribution_for(found, hotspot_window_aa=hotspot_window_aa),
                source=source,
                citation=(
                    f"{len(found)} pathogenic submissions at ≥{minimum_stars}★ "
                    f"({review_summary(found)}); hotspot = densest "
                    f"{hotspot_window_aa}-residue window"
                ),
            )
        )
    return len(by_gene)
