"""Reading a VCF: the patient's own data, as a caller emitted it.

Everything in ``annotation/`` and ``context/`` is reference material — files that
describe the genome or the literature. This reads the one file that describes the
patient, and it is the only input the project has that nobody else versions.

Four things it is careful about, and each of them is a way to be quietly wrong
rather than loudly broken.

**Genotype is where zygosity comes from**, and zygosity decides roughly half the
modalities. ``1/2`` is the case worth knowing: two different alternate alleles at
one site means the sample carries *no reference allele*, so every intervention
that works by raising output from an intact copy is off the table — the same
conclusion as a homozygote, reached by a different route.

**Multi-allelic records are decomposed, and decomposition without normalisation
is wrong.** Splitting ``A>C,CTT`` gives two records whose representations are not
canonical, and a variant written non-canonically fails to match the same variant
written properly. So a reader handed a sequence provider normalises on the way
out, and one without a reference says so rather than pretending.

**Consequence is not predicted here.** This project reads consequences, it does
not derive them; that is VEP's job. If the VCF carries a CSQ or ANN field the
reader takes it, and if it does not, the caller has to supply one. Guessing
"missense because the alleles are the same length" would be a prediction dressed
as a parse.

**The assembly is checked when the file declares one.** A VCF called against a
different build than the annotation produces coordinates that are all plausible
and all about the wrong part of the genome.
"""

from __future__ import annotations

import gzip
import re
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from repairbench.annotation.fasta import SequenceProvider
from repairbench.annotation.normalise import left_align
from repairbench.model import Consequence, RepairbenchError, Zygosity

PASS_FILTERS = frozenset({"PASS", ".", ""})

#: Sequence Ontology terms a VEP/SnpEff annotation may carry, mapped to the
#: coarse vocabulary the rules read. Terms outside this map are not errors —
#: a VCF carries many — they simply do not answer the question this project asks.
_CONSEQUENCE_TERMS: dict[str, Consequence] = {
    "stop_gained": Consequence.NONSENSE,
    "frameshift_variant": Consequence.FRAMESHIFT,
    "splice_acceptor_variant": Consequence.SPLICE_ACCEPTOR,
    "splice_donor_variant": Consequence.SPLICE_DONOR,
    "splice_region_variant": Consequence.SPLICE_REGION,
    "missense_variant": Consequence.MISSENSE,
    "inframe_deletion": Consequence.INFRAME_DELETION,
    "inframe_insertion": Consequence.INFRAME_INSERTION,
    "start_lost": Consequence.START_LOST,
    "stop_lost": Consequence.STOP_LOST,
    "synonymous_variant": Consequence.SYNONYMOUS,
    "intron_variant": Consequence.INTRONIC,
}

#: Consequences ordered by how much they matter to this project, so that a
#: multi-consequence annotation resolves to the one the rules act on rather than
#: to whichever VEP happened to list first.
_SEVERITY: tuple[Consequence, ...] = (
    Consequence.NONSENSE,
    Consequence.FRAMESHIFT,
    Consequence.SPLICE_ACCEPTOR,
    Consequence.SPLICE_DONOR,
    Consequence.START_LOST,
    Consequence.STOP_LOST,
    Consequence.MISSENSE,
    Consequence.INFRAME_DELETION,
    Consequence.INFRAME_INSERTION,
    Consequence.SPLICE_REGION,
    Consequence.SYNONYMOUS,
    Consequence.INTRONIC,
)

_ASSEMBLY_HINTS = {"GRCh38": ("GRCh38", "hg38"), "GRCh37": ("GRCh37", "hg19", "b37")}


class VcfError(RepairbenchError):
    """The VCF is malformed, or does not carry what was asked of it."""


@dataclass(frozen=True, slots=True)
class ObservedVariant:
    """One alternate allele, as one sample carries it.

    A record straight out of a VCF, before any interpretation. It knows nothing
    about mechanisms — that is the point of the boundary.
    """

    chromosome: str
    position: int
    reference: str
    alternate: str
    zygosity: Zygosity
    sample: str
    filters: tuple[str, ...] = ()
    gene: str | None = None
    consequence: Consequence | None = None
    #: True when the reader left-aligned this record against a reference.
    normalised: bool = False

    @property
    def key(self) -> str:
        return f"{self.chromosome}-{self.position}-{self.reference}-{self.alternate}"

    @property
    def passed_filters(self) -> bool:
        return not self.filters or set(self.filters) <= PASS_FILTERS

    @property
    def is_interpretable(self) -> bool:
        """Can the rules do anything with this?

        Both a gene and a consequence are needed, and neither is derived here.
        A record missing them is not an error — most records in a genome are —
        it simply is not something this project can reason about.
        """
        return self.gene is not None and self.consequence is not None


def zygosity_from_genotype(genotype: str) -> Zygosity:
    """Read zygosity from a GT field.

    The cases, and why each answer is what it is:

    * ``0/1``, ``0|1`` — heterozygous. One reference allele remains.
    * ``1/1`` — homozygous. No reference allele.
    * ``1`` — a haploid call: hemizygous, as on a male X outside the
      pseudoautosomal regions.
    * ``1/2`` — two *different* alternate alleles. The sample is heterozygous
      for each of them and carries no reference allele at all, which for every
      rule that asks "is there an intact copy" is the same answer as a
      homozygote. Reported as compound heterozygous so the distinction survives
      in the record even though the consequence is identical.
    * ``./.``, ``.`` — no call. Unknown, and unknown is not a guess.

    A note on the case this cannot see: many callers emit ``1/1`` for a male X
    chromosome rather than a haploid ``1``. That reads here as homozygous rather
    than hemizygous — and the two agree on the only question the rules ask of
    zygosity, so the misreading is inert. Distinguishing them properly needs the
    sample's sex, which a VCF does not carry.
    """
    alleles = [allele for allele in re.split(r"[/|]", genotype.strip()) if allele != ""]
    if not alleles or all(allele == "." for allele in alleles):
        return Zygosity.UNKNOWN
    if len(alleles) == 1:
        return Zygosity.HEMIZYGOUS

    called = {allele for allele in alleles if allele != "."}
    if "0" in called:
        return Zygosity.HETEROZYGOUS
    if len(called) == 1:
        return Zygosity.HOMOZYGOUS
    return Zygosity.COMPOUND_HETEROZYGOUS


def _open(path: Path) -> Iterator[str]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as handle:
            yield from handle
    else:
        with path.open() as handle:
            yield from handle


@dataclass(frozen=True, slots=True)
class VcfReader:
    """Reads one VCF, optionally normalising and filtering as it goes."""

    path: Path
    sample: str | None = None
    #: When set, records whose FILTER is not PASS are dropped. A caller that
    #: wants them must ask for them, because a failed filter is the caller
    #: saying it does not believe its own call.
    require_pass: bool = True
    expected_assembly: str | None = None

    def read(self, sequences: SequenceProvider | None = None) -> list[ObservedVariant]:
        """Parse the file into one record per carried alternate allele."""
        header_samples: list[str] = []
        observed: list[ObservedVariant] = []
        seen_header = False

        for number, raw_line in enumerate(_open(self.path), start=1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                self._check_meta(line)
                continue
            if line.startswith("#CHROM"):
                header_samples = line.split("\t")[9:]
                seen_header = True
                continue
            if not seen_header:
                raise VcfError(f"{self.path.name}:{number}: data before the #CHROM header")
            observed.extend(self._read_record(line, header_samples, number, sequences))

        if not seen_header:
            raise VcfError(f"{self.path.name}: no #CHROM header line")
        return observed

    def _check_meta(self, line: str) -> None:
        if not line.startswith("##reference") or self.expected_assembly is None:
            return
        hints = _ASSEMBLY_HINTS.get(self.expected_assembly, (self.expected_assembly,))
        if not any(hint.lower() in line.lower() for hint in hints):
            raise VcfError(
                f"{self.path.name}: {line.strip()} does not look like "
                f"{self.expected_assembly}. A VCF called against a different build gives "
                "coordinates that are all plausible and all about the wrong part of the genome."
            )

    def _read_record(
        self,
        line: str,
        header_samples: list[str],
        number: int,
        sequences: SequenceProvider | None,
    ) -> list[ObservedVariant]:
        columns = line.split("\t")
        minimum_columns = 8
        if len(columns) < minimum_columns:
            raise VcfError(f"{self.path.name}:{number}: expected at least 8 columns")

        chromosome, position, _, reference, alternates, _, filter_field, info = columns[:8]
        filters = tuple(f for f in filter_field.split(";") if f)
        if self.require_pass and not set(filters) <= PASS_FILTERS:
            return []

        genotype = self._genotype(columns, header_samples, number)
        sample_name = self._sample_name(header_samples)
        annotations = _parse_annotations(info)

        records: list[ObservedVariant] = []
        for index, alternate in enumerate(alternates.split(","), start=1):
            if alternate in {".", "*"}:
                # A spanning deletion placeholder. It is a statement about
                # another record, not a variant of its own.
                continue
            zygosity = _zygosity_for_allele(genotype, index)
            if zygosity is None:
                continue

            chrom, pos, ref, alt = chromosome, int(position), reference, alternate
            normalised = False
            if sequences is not None:
                aligned = left_align(chrom, pos, ref, alt, sequences)
                chrom, pos, ref, alt = (
                    aligned.chromosome,
                    aligned.position,
                    aligned.reference,
                    aligned.alternate,
                )
                normalised = True

            gene, consequence = annotations.get(alternate, (None, None))
            records.append(
                ObservedVariant(
                    chromosome=chrom,
                    position=pos,
                    reference=ref,
                    alternate=alt,
                    zygosity=zygosity,
                    sample=sample_name,
                    filters=filters,
                    gene=gene,
                    consequence=consequence,
                    normalised=normalised,
                )
            )
        return records

    def _genotype(self, columns: list[str], header_samples: list[str], number: int) -> str:
        sample_columns = columns[9:]
        if not sample_columns:
            return ""
        index = 0
        if self.sample is not None:
            if self.sample not in header_samples:
                raise VcfError(
                    f"{self.path.name}: no sample {self.sample!r}; it has "
                    f"{', '.join(header_samples) or 'none'}"
                )
            index = header_samples.index(self.sample)
        elif len(header_samples) > 1:
            raise VcfError(
                f"{self.path.name} has {len(header_samples)} samples "
                f"({', '.join(header_samples)}). Name the one to read — picking the first "
                "would silently interpret a parent's genotype as the child's."
            )

        format_keys = columns[8].split(":") if len(columns) > 8 else []
        if "GT" not in format_keys:
            raise VcfError(f"{self.path.name}:{number}: FORMAT carries no GT")
        values = sample_columns[index].split(":")
        return values[format_keys.index("GT")]

    def _sample_name(self, header_samples: list[str]) -> str:
        if self.sample is not None:
            return self.sample
        return header_samples[0] if header_samples else "unnamed"


def _zygosity_for_allele(genotype: str, allele_index: int) -> Zygosity | None:
    """Zygosity for one alternate allele, or ``None`` when the sample does not
    carry it at all."""
    if not genotype:
        return Zygosity.UNKNOWN
    alleles = [allele for allele in re.split(r"[/|]", genotype.strip()) if allele != ""]
    if str(allele_index) not in alleles:
        return None
    return zygosity_from_genotype(genotype)


def _parse_annotations(info: str) -> dict[str, tuple[str | None, Consequence | None]]:
    """Read gene and consequence out of a VEP ``CSQ`` or SnpEff ``ANN`` field.

    Both formats put the alternate allele first and a pipe-separated list after
    it; the fields differ, so this looks for the terms it knows rather than
    counting positions. Positional parsing would break the moment somebody ran
    VEP with a different ``--fields``.
    """
    found: dict[str, tuple[str | None, Consequence | None]] = {}
    for entry in info.split(";"):
        key, _, value = entry.partition("=")
        if key not in {"CSQ", "ANN"} or not value:
            continue
        for annotation in value.split(","):
            parts = annotation.split("|")
            if len(parts) < 2:
                continue
            allele = parts[0]
            consequence = _most_severe(parts)
            gene = _gene_symbol(parts)
            previous = found.get(allele)
            if previous is None or _outranks(consequence, previous[1]):
                found[allele] = (gene or (previous[0] if previous else None), consequence)
    return found


def _most_severe(parts: list[str]) -> Consequence | None:
    terms: set[Consequence] = set()
    for part in parts:
        for term in part.split("&"):
            mapped = _CONSEQUENCE_TERMS.get(term.strip())
            if mapped is not None:
                terms.add(mapped)
    for candidate in _SEVERITY:
        if candidate in terms:
            return candidate
    return None


def _gene_symbol(parts: list[str]) -> str | None:
    """The gene symbol, found by shape rather than by position.

    A symbol is an uppercase-ish token that is neither a consequence term nor an
    identifier. This is a heuristic and is written down as one — the alternative
    is hard-coding a field index that changes with the annotator's options.
    """
    for part in parts[1:]:
        candidate = part.strip()
        if not candidate or candidate in _CONSEQUENCE_TERMS:
            continue
        if candidate.startswith(("ENSG", "ENST", "NM_", "NR_", "XM_")):
            continue
        if re.fullmatch(r"[A-Z][A-Z0-9\-]{1,14}", candidate):
            return candidate
    return None


def _outranks(candidate: Consequence | None, existing: Consequence | None) -> bool:
    if candidate is None:
        return False
    if existing is None:
        return True
    return _SEVERITY.index(candidate) < _SEVERITY.index(existing)
