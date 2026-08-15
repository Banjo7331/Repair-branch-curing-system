"""Where a gene is switched on, as GTEx measured it.

The question this answers is narrow and worth stating precisely, because tissue
is a dimension where it is easy to claim more than the data supports.

**What it answers:** is this gene transcribed in the tissue the patient's disease
affects? If it is not, a variant in it is unlikely to be the cause, and every
intervention that works by supplying or silencing its product *there* has nothing
to act on.

**What it does not answer:** whether a delivery route reaches that tissue,
whether the tissue's dominant transcript is the one the rules were run against,
or whether the gene was expressed at the developmental moment that mattered.
GTEx measures adult post-mortem bulk tissue. A gene switched off at fifty may
have been essential at four weeks, and nothing in this file can tell.

The file format is GTEx's median-TPM matrix: two identifier columns and one
column per tissue, parsed by header name because the tissue set changes between
releases. GTEx ships it as GCT, which carries a two-line preamble before the
header; strip that first, or export as TSV. The parser does not guess at it,
because a preamble it mis-skipped would shift every column by one and produce
expression values that are all plausible and all for the wrong tissue.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

from repairbench.context.source import ContextError, Fact, Provenance, Source, read_tsv

GENE_COLUMN = "Description"
#: The Ensembl identifier column. Present, unused: this project joins on symbol
#: because that is what ClinGen and the rule files speak.
ID_COLUMN = "Name"


class TissueSystem(StrEnum):
    """Coarse groupings, for the questions delivery asks.

    Deliberately much coarser than GTEx's fifty-odd tissues. The fine vocabulary
    is right for "is this gene on here"; the coarse one is right for "could
    anything get there", and conflating them would let a rule about the brain
    key off *Brain - Nucleus accumbens (basal ganglia)*.
    """

    CNS = "central_nervous_system"
    MUSCLE = "skeletal_muscle"
    HEART = "heart"
    LIVER = "liver"
    KIDNEY = "kidney"
    EYE = "eye"
    BLOOD = "blood"
    SKIN = "skin"
    LUNG = "lung"
    OTHER = "other"


def all_tissue_systems() -> tuple[TissueSystem, ...]:
    return tuple(TissueSystem)


#: Prefixes GTEx uses, mapped to the coarse system. Order matters: the first
#: match wins, so more specific prefixes come first.
_SYSTEM_PREFIXES: tuple[tuple[str, TissueSystem], ...] = (
    ("Brain", TissueSystem.CNS),
    ("Nerve", TissueSystem.CNS),
    ("Pituitary", TissueSystem.CNS),
    ("Muscle", TissueSystem.MUSCLE),
    ("Heart", TissueSystem.HEART),
    ("Liver", TissueSystem.LIVER),
    ("Kidney", TissueSystem.KIDNEY),
    ("Whole Blood", TissueSystem.BLOOD),
    ("Spleen", TissueSystem.BLOOD),
    ("Skin", TissueSystem.SKIN),
    ("Lung", TissueSystem.LUNG),
)


def system_for(tissue: str) -> TissueSystem:
    """Map a GTEx tissue name onto its coarse system.

    Anything unrecognised is ``OTHER`` rather than an error. A new GTEx release
    adding a tissue should not stop a run; it should produce an answer that is
    honest about being unclassified.
    """
    for prefix, system in _SYSTEM_PREFIXES:
        if tissue.startswith(prefix):
            return system
    return TissueSystem.OTHER


@dataclass(frozen=True, slots=True)
class Tissue:
    """The tissue a case is about.

    Case-scoped, like the phenotype: a patient's affected tissue is a fact about
    the patient, not a release. The GTEx name is kept verbatim so a lookup is
    exact, and the coarse system is derived so delivery rules have something
    they can key on without matching on a basal ganglia subregion.
    """

    name: str

    @property
    def system(self) -> TissueSystem:
        return system_for(self.name)


def ingest(source: Source, into: dict[str, Provenance], *, genes: set[str] | None = None) -> int:
    """Read median expression into a provenance map, one gene at a time.

    The whole per-tissue row is kept rather than a single number, because which
    tissue matters is a property of the case and is not known here. Reducing it
    early would mean re-reading the file for every case.
    """
    recorded = 0
    for row in read_tsv(source, required={GENE_COLUMN}):
        gene = row[GENE_COLUMN]
        if not gene or (genes is not None and gene not in genes):
            continue

        profile: dict[str, float] = {}
        for column, value in row.items():
            if column in {GENE_COLUMN, ID_COLUMN} or not value:
                continue
            try:
                profile[column] = float(value)
            except ValueError as error:
                raise ContextError(
                    f"{gene}: expression in {column!r} is {value!r}, which is not a number"
                ) from error

        if not profile:
            continue
        into.setdefault(gene, Provenance(gene=gene)).record(
            Fact(field="expression", value=profile, source=source)
        )
        recorded += 1
    return recorded


def median_tpm(profile: dict[str, float] | None, tissue: Tissue | None) -> float | None:
    """Expression of one gene in one tissue, or ``None`` when either is unknown.

    ``None`` rather than zero, and the distinction carries: zero means GTEx
    measured this gene here and found nothing, which is evidence. ``None`` means
    nobody looked, which is not.
    """
    if profile is None or tissue is None:
        return None
    return profile.get(tissue.name)
