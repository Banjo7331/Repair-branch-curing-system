"""Assembling a gene's context out of several pinned sources.

The registry is the join. It reads each source once, keeps the provenance of
every fact, and hands back a ``Gene`` the rules can run against — plus, for each
gene, the record of where every field in it came from.

What it will not do is fill a gap. A gene missing from ClinGen gets no dosage
score rather than a default of "no evidence", because those are different
claims: one is what ClinGen said, the other is what ClinGen has not looked at.
The mechanism rules already refuse to fire on a missing value, so a gap
propagates into an honest "undetermined" instead of a confident wrong answer.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from repairbench.context import clingen, clinvar, curation, expression, gnomad
from repairbench.context.source import ContextError, Provenance, Source
from repairbench.model import DosageScore, Gene, MissenseDistribution


@dataclass(frozen=True, slots=True)
class SourcedGene:
    """A gene context, and the record of where each fact in it came from."""

    gene: Gene
    provenance: Provenance

    @property
    def expression(self) -> dict[str, float] | None:
        """Median TPM per tissue, when an expression source supplied it.

        Kept off ``Gene`` on purpose. Everything on ``Gene`` is a property of the
        gene alone; expression is a property of the gene *in a tissue*, and the
        tissue is a fact about the case. Folding it in would invite a rule to
        read "expression" as though it were a single number.
        """
        fact = self.provenance.facts.get("expression")
        return dict(fact.value) if fact else None

    def cite(self, field_name: str) -> str:
        source = self.provenance.source_for(field_name)
        return source.pin if source else "not supplied by any source"


class GeneContextRegistry:
    """Gene-level facts, read from pinned files."""

    def __init__(self, sources: dict[str, Provenance], pins: dict[str, Source]) -> None:
        self._genes = sources
        self._pins = pins

    @classmethod
    def load(
        cls,
        *,
        dosage: Path | None = None,
        constraint: Path | None = None,
        local: Path | None = None,
        expression_matrix: Path | None = None,
        variant_summary: Path | None = None,
        dosage_version: str = "unversioned",
        constraint_version: str = "unversioned",
        local_version: str = "unversioned",
        expression_version: str = "unversioned",
        clinvar_version: str = "unversioned",
        minimum_stars: int = 1,
        hotspot_window_aa: int = 20,
        genes: set[str] | None = None,
    ) -> GeneContextRegistry:
        """Read whichever sources are supplied.

        All of them are optional and none is defaulted to an empty stand-in: a
        registry loaded without gnomAD simply has no constraint for any gene,
        which the rules handle, rather than a fabricated one they would act on.
        """
        collected: dict[str, Provenance] = {}
        pins: dict[str, Source] = {}

        for path, version, name, reader in (
            (dosage, dosage_version, "clingen_dosage", clingen.ingest),
            (constraint, constraint_version, "gnomad_constraint", gnomad.ingest),
            (local, local_version, "local_curation", curation.ingest),
            (expression_matrix, expression_version, "expression", expression.ingest),
        ):
            if path is None:
                continue
            source = Source.of(name, path, version)
            reader(source, collected, genes=genes)
            pins[name] = source

        if variant_summary is not None:
            # Read outside the loop because it is the one source that cannot be
            # read whole. The other four are gene-per-row tables of a few
            # hundred thousand lines; ClinVar's summary is millions of
            # submissions, and a registry loaded without a gene filter would
            # spend minutes reading them to answer a question about nine.
            if not genes:
                raise ContextError(
                    "loading ClinVar needs the genes to look for: the file is millions of "
                    "submissions, and reading all of them to build context for a handful "
                    "would take minutes and look like it was working"
                )
            source = Source.of("clinvar", variant_summary, clinvar_version)
            clinvar.ingest(
                source,
                collected,
                genes=genes,
                minimum_stars=minimum_stars,
                hotspot_window_aa=hotspot_window_aa,
            )
            pins["clinvar"] = source

        if not pins:
            raise ContextError("a registry with no sources cannot supply anything")
        return cls(collected, pins)

    @property
    def pins(self) -> tuple[Source, ...]:
        return tuple(self._pins.values())

    @property
    def genes(self) -> list[str]:
        return sorted(self._genes)

    def describe(self) -> str:
        return "\n".join(source.pin for source in self.pins)

    def provenance_for(self, symbol: str) -> Provenance:
        provenance = self._genes.get(symbol)
        if provenance is not None and not provenance.facts:
            provenance = None
        if provenance is None:
            raise ContextError(
                f"{symbol} is not present in any loaded source "
                f"({', '.join(source.pin for source in self.pins)})"
            )
        return provenance

    def gene(self, symbol: str, *, distribution: MissenseDistribution | None = None) -> SourcedGene:
        """Build the gene context the rules read.

        ``distribution`` — where pathogenic variation sits in the gene — comes
        from ClinVar when a summary was loaded, and the parameter now overrides
        it rather than being the only way to supply it. The override stays
        because a curated count from a disease-specific database is better
        evidence than a submission tally, and the reference set has to be able
        to hold a gene still while a rule is being examined.

        What it does *not* do is default to zero when no source supplied one.
        An empty distribution reads to the clustering rule as "no pattern",
        which is the correct thing for it to conclude from no data — but the
        provenance says nothing was counted, so a report cannot present it as a
        finding.
        """
        provenance = self.provenance_for(symbol)
        facts = provenance.facts

        def value(field_name: str, fallback: object = None) -> object:
            fact = facts.get(field_name)
            return fact.value if fact else fallback

        gene = Gene(
            symbol=symbol,
            # A gene ClinGen has not evaluated gets *not evaluated*, not a
            # score. The registry already refuses to invent a value; this is
            # the same refusal spelled where the Gene is built.
            haploinsufficiency=value("haploinsufficiency", DosageScore.NOT_EVALUATED),  # type: ignore[arg-type]
            triplosensitivity=value("triplosensitivity", DosageScore.NOT_EVALUATED),  # type: ignore[arg-type]
            loeuf=value("loeuf"),  # type: ignore[arg-type]
            forms_multimer=bool(value("forms_multimer", False)),
            truncating_variants_are_milder=bool(value("truncating_variants_are_milder", False)),
            distribution=distribution or value("distribution", MissenseDistribution()),  # type: ignore[arg-type]
            curated_mechanism=value("curated_mechanism"),  # type: ignore[arg-type]
            curated_mechanism_source=(
                facts["curated_mechanism"].citation if "curated_mechanism" in facts else None
            ),
        )
        return SourcedGene(gene=gene, provenance=provenance)
