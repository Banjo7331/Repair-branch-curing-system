"""gnomAD constraint, as published — in whichever schema the release uses.

One number is taken from this file: LOEUF, the upper bound of the confidence
interval on observed-over-expected loss-of-function variation. Low means the
population carries fewer null alleles in this gene than chance would predict,
which is consistent with dosage sensitivity.

**The column names are not stable between releases, and this is where that
stops being a footnote.** The v2.1.1 constraint file — the one the downloads
page has offered longest — calls the number ``oe_lof_upper`` and marks the row
to use with ``canonical``, Ensembl's canonical transcript. The v4 files call the
same number ``lof.oe_ci.upper`` and mark the row with ``mane_select``. This
module was written against the second set of names, and pointed at the file the
site actually links to it would fail with "no such column" — the correct
failure, and a useless one to a reader who has the right file in hand.

So both schemas are declared, detected from the header, and the refusal names
every spelling it looked for. What is *not* done is matching a column by
position or by resemblance: a constraint value read out of the wrong column is a
number that looks entirely reasonable.

**Which row is used differs with the schema, and the difference is real.** MANE
Select is an agreement between RefSeq and Ensembl about the clinical reference
transcript; Ensembl canonical is Ensembl's own pick. They usually agree, they
are not the same claim, and the provenance records which one a value came
through.

**LOEUF is kept and no verdict is derived from it.** Constraint is a statement
about a *gene* across a population, not about the variant in front of you, and
the rule that reads it is deliberately marked supporting for exactly that
reason. Turning it into a boolean here would move a judgement out of the rule
file, which is the one place this project keeps them.
"""

from __future__ import annotations

from dataclasses import dataclass

from repairbench.context.source import ContextError, Fact, Provenance, Source, read_tsv

GENE_COLUMN = "gene"

_TRUE = {"true", "TRUE", "True", "1", "yes"}
_ABSENT = {"", "NA", "NaN", "."}


@dataclass(frozen=True, slots=True)
class Schema:
    """One release's spelling of the columns this module needs."""

    release: str
    loeuf: str
    preferred_row: str
    #: What the preferred-row column means, for the provenance line.
    row_meaning: str

    @property
    def columns(self) -> set[str]:
        return {GENE_COLUMN, self.loeuf, self.preferred_row}


#: Newest first, so a deployment holding both reads the newer one.
SCHEMAS: tuple[Schema, ...] = (
    Schema(
        release="v4",
        loeuf="lof.oe_ci.upper",
        preferred_row="mane_select",
        row_meaning="MANE Select",
    ),
    Schema(
        release="v2.1.1",
        loeuf="oe_lof_upper",
        preferred_row="canonical",
        row_meaning="Ensembl canonical",
    ),
)


def detect_schema(source: Source) -> Schema:
    """Which release's column names this file uses.

    Read from the header rather than from the file name, because a file renamed
    on the way to disk is ordinary and its columns are not.
    """
    header = next(read_tsv(source, required=set()), None)
    if header is None:
        raise ContextError(f"{source.path.name}: no rows to read a schema from")

    for schema in SCHEMAS:
        if schema.columns <= set(header):
            return schema

    wanted = "; ".join(
        f"{schema.release} wants {', '.join(sorted(schema.columns))}" for schema in SCHEMAS
    )
    raise ContextError(
        f"{source.path.name} matches no constraint schema this package knows ({wanted}). "
        f"It has: {', '.join(sorted(header))}"
    )


def ingest(source: Source, into: dict[str, Provenance], *, genes: set[str] | None = None) -> int:
    """Read constraint into a provenance map."""
    schema = detect_schema(source)
    recorded = 0

    for row in read_tsv(source, required=schema.columns):
        gene = row[GENE_COLUMN]
        if not gene or (genes is not None and gene not in genes):
            continue
        if row[schema.preferred_row] not in _TRUE:
            # One row per transcript; the numbers differ between them, and
            # taking whichever came first would make the value depend on file
            # ordering rather than on biology.
            continue

        raw = row[schema.loeuf]
        if raw in _ABSENT:
            # Constraint is undefined for short genes: too few expected variants
            # for the ratio to mean anything. Absent is the correct answer, and
            # the rules already refuse to fire on a missing value.
            continue
        try:
            loeuf = float(raw)
        except ValueError as error:
            raise ContextError(
                f"{gene}: {schema.loeuf} is {raw!r}, which is not a number"
            ) from error

        into.setdefault(gene, Provenance(gene=gene)).record(
            Fact(
                field="loeuf",
                value=loeuf,
                source=source,
                citation=f"gnomAD {schema.release}, {schema.row_meaning} transcript",
            )
        )
        recorded += 1
    return recorded
