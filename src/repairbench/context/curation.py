"""The facts nobody publishes as a table.

Whether a gene product assembles into a complex, and whether its null alleles
are milder than its missense ones, are judgements read out of the literature.
There is no download for them. That is not a gap in this project's ingest — it
is the state of the field, and the honest response is to keep them somewhere
that looks nothing like an ingested table.

So this file is different in three ways, each of them enforced rather than
suggested:

* **Every entry demands a citation.** A claim that *COL1A1* null alleles are
  milder is a clinical assertion; without a pointer to where that was
  established it is folklore with a version number.
* **It is pinned like everything else**, so a report can say that this fact came
  from our own curation at a particular revision — and a reviewer who disagrees
  knows exactly which file to argue with.
* **It supplies only the two fields that have no public table.** Anything ClinGen
  or gnomAD publishes must come from ClinGen or gnomAD, and the loader refuses
  local overrides of them outright.
"""

from __future__ import annotations

from typing import Any

import yaml

from repairbench.context.source import ContextError, Fact, Provenance, Source

#: The only fields local curation may supply. Everything else has a public
#: source, and a local override of a published fact is a way to be quietly wrong.
CURATABLE = frozenset({"forms_multimer", "truncating_variants_are_milder", "curated_mechanism"})


def ingest(source: Source, into: dict[str, Provenance], *, genes: set[str] | None = None) -> int:
    """Read local curation into a provenance map."""
    document = yaml.safe_load(source.path.read_text()) or {}
    if not isinstance(document, dict) or "genes" not in document:
        raise ContextError(f"{source.path.name}: expected a mapping with a 'genes' key")

    recorded = 0
    for gene, entry in (document["genes"] or {}).items():
        if genes is not None and gene not in genes:
            continue
        if not isinstance(entry, dict):
            raise ContextError(f"{source.path.name}: entry for {gene} is not a mapping")

        provenance = into.setdefault(gene, Provenance(gene=gene))
        for field_name, claim in entry.items():
            recorded += _record(provenance, gene, field_name, claim, source)
    return recorded


def _record(
    provenance: Provenance,
    gene: str,
    field_name: str,
    claim: Any,
    source: Source,
) -> int:
    if field_name not in CURATABLE:
        raise ContextError(
            f"{gene}: local curation may not supply {field_name!r}. "
            f"It is published by a source of its own; only {', '.join(sorted(CURATABLE))} "
            "have no public table and belong here."
        )
    if not isinstance(claim, dict) or "value" not in claim:
        raise ContextError(
            f"{gene}.{field_name}: expected a mapping with 'value' and 'citation'"
        )

    citation = str(claim.get("citation", "")).strip()
    if not citation:
        raise ContextError(
            f"{gene}.{field_name}: no citation. A clinical assertion without a pointer to "
            "where it was established is folklore with a version number."
        )

    provenance.record(
        Fact(field=field_name, value=claim["value"], source=source, citation=citation)
    )
    return 1
