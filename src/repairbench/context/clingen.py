"""ClinGen dosage sensitivity, as published.

The file is the gene curation list ClinGen ships for each assembly. What this
module takes from it is two columns and one refusal.

The refusal is the interesting part. ClinGen's scale is not a scale: 0 through 3
rank evidence for haploinsufficiency, but 30 means the gene acts recessively and
40 means dosage sensitivity was *actively refuted*. Flattening those into "low"
would erase the distinction the mechanism rules most need — "no evidence yet" is
a gap and "refuted" is a finding, and a predicted null variant means something
different under each.
"""

from __future__ import annotations

from repairbench.context.source import ContextError, Fact, Provenance, Source, read_tsv
from repairbench.model import DosageScore

GENE_COLUMN = "Gene Symbol"
HAPLO_COLUMN = "Haploinsufficiency Score"
TRIPLO_COLUMN = "Triplosensitivity Score"

#: ClinGen's numeric codes, in their own vocabulary rather than flattened.
_SCORES: dict[str, DosageScore] = {
    "0": DosageScore.NO_EVIDENCE,
    "1": DosageScore.LITTLE_EVIDENCE,
    "2": DosageScore.EMERGING_EVIDENCE,
    "3": DosageScore.SUFFICIENT_EVIDENCE,
    "30": DosageScore.AUTOSOMAL_RECESSIVE,
    "40": DosageScore.UNLIKELY,
}

#: Values that mean "nobody has looked". Distinct from a score of 0, which means
#: somebody looked and found nothing.
_NOT_EVALUATED = {"", "Not yet evaluated", "Not evaluated", "-"}


def parse_score(raw: str, *, gene: str, column: str) -> DosageScore | None:
    """Map one cell, refusing values that are neither a score nor an absence."""
    if raw in _NOT_EVALUATED:
        return None
    score = _SCORES.get(raw)
    if score is None:
        raise ContextError(
            f"{gene}: {column} is {raw!r}, which is not a ClinGen dosage code. "
            f"Known codes: {', '.join(sorted(_SCORES))}"
        )
    return score


def ingest(source: Source, into: dict[str, Provenance], *, genes: set[str] | None = None) -> int:
    """Read dosage curation into a provenance map, one fact at a time.

    ``genes`` restricts the parse, which is what makes this usable against the
    full list: the file covers the genome and a case concerns a handful of genes.

    A gene that is present but unevaluated contributes nothing rather than a
    default. That is deliberate — recording ``no_evidence`` where ClinGen has
    simply not looked would put a claim in the provenance that ClinGen never
    made.
    """
    recorded = 0
    for row in read_tsv(source, required={GENE_COLUMN, HAPLO_COLUMN, TRIPLO_COLUMN}):
        gene = row[GENE_COLUMN]
        if not gene or (genes is not None and gene not in genes):
            continue

        # Parse before touching the map. Creating the entry first would make a
        # gene ClinGen has not evaluated *present with nothing in it*, which
        # reads downstream as "we have context for this gene" — the exact
        # failure this module exists to avoid.
        scores = {
            field_name: parse_score(row[column], gene=gene, column=column)
            for column, field_name in (
                (HAPLO_COLUMN, "haploinsufficiency"),
                (TRIPLO_COLUMN, "triplosensitivity"),
            )
        }
        supplied = {name: score for name, score in scores.items() if score is not None}
        if not supplied:
            continue

        provenance = into.setdefault(gene, Provenance(gene=gene))
        for field_name, score in supplied.items():
            provenance.record(Fact(field=field_name, value=score, source=source))
            recorded += 1
    return recorded
