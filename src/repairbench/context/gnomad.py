"""gnomAD constraint, as published.

One number is taken from this file: LOEUF, the upper bound of the confidence
interval on observed-over-expected loss-of-function variation. Low means the
population carries fewer null alleles in this gene than chance would predict,
which is consistent with dosage sensitivity.

Two things this module is careful about.

**It reads the MANE Select row.** The constraint file has one row per
transcript, and the numbers differ between them. Taking whichever row came first
would make the value depend on file ordering.

**It keeps LOEUF and does not derive a verdict from it.** Constraint is a
statement about a *gene* across a population, not about the variant in front of
you, and the rule that reads it is deliberately marked supporting for exactly
that reason. Turning it into a boolean here would move a judgement out of the
rule file, which is the one place this project keeps them.
"""

from __future__ import annotations

from repairbench.context.source import ContextError, Fact, Provenance, Source, read_tsv

GENE_COLUMN = "gene"
LOEUF_COLUMN = "lof.oe_ci.upper"
MANE_COLUMN = "mane_select"

_TRUE = {"true", "TRUE", "True", "1", "yes"}


def ingest(source: Source, into: dict[str, Provenance], *, genes: set[str] | None = None) -> int:
    """Read constraint into a provenance map."""
    recorded = 0
    for row in read_tsv(source, required={GENE_COLUMN, LOEUF_COLUMN, MANE_COLUMN}):
        gene = row[GENE_COLUMN]
        if not gene or (genes is not None and gene not in genes):
            continue
        if row[MANE_COLUMN] not in _TRUE:
            # One row per transcript; the numbers differ between them, and
            # taking whichever came first would make the value depend on file
            # ordering rather than on biology.
            continue

        raw = row[LOEUF_COLUMN]
        if raw in {"", "NA", "NaN"}:
            # Constraint is undefined for short genes: too few expected variants
            # for the ratio to mean anything. Absent is the correct answer, and
            # the rules already refuse to fire on a missing value.
            continue
        try:
            loeuf = float(raw)
        except ValueError as error:
            raise ContextError(
                f"{gene}: {LOEUF_COLUMN} is {raw!r}, which is not a number"
            ) from error

        into.setdefault(gene, Provenance(gene=gene)).record(
            Fact(field="loeuf", value=loeuf, source=source)
        )
        recorded += 1
    return recorded
