"""Tiling antisense oligonucleotides along a target, and saying what is unknown.

An ASO is the simplest molecule in this package to design and the hardest to
design *well*, and the gap between those two is the reason this module is shaped
the way it is.

Simple: slide a window along the target, take the reverse complement, check the
composition. Every window of the right length is a candidate, and a tiling run
produces hundreds of them.

Hard: which of them works depends on whether the site is **accessible** — on
whether that stretch of transcript is paired up inside the molecule's own fold,
occupied by a protein, or exposed. Two candidates identical on every composition
rule can differ tenfold for that reason. Predicting it needs a folding model,
and no folding model is attached here.

So this module does what it can do honestly. It tiles, it applies the
composition and motif rules from a pinned file, it places each window against
the exon boundaries when an annotation is supplied — a splice acceptor and an
exon interior are different targets — and it refuses the one confusion that
inverts the therapy: a cleaving chemistry aimed at a splice site destroys the
transcript that the splice redirection was supposed to rescue.

Everything it does not know is a first-class part of the output, because a list
of two hundred oligonucleotides is exactly the kind of result that reads as an
answer.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from typing import Protocol

from repairbench.design.editors import DesignError
from repairbench.design.flags import (
    Flag,
    FlagRuleset,
    FlatFeatures,
    Severity,
    sort_weight,
    worst_of,
)
from repairbench.design.sequence import is_resolved, reverse_complement


class Action(StrEnum):
    """What a chemistry does to the transcript it binds."""

    CLEAVES = "cleaves"
    """Recruits RNase H. Destroys the transcript — the knockdown instrument."""

    BLOCKS = "blocks"
    """Occupies a site and recruits nothing. Redirects splicing or translation."""


@dataclass(frozen=True, slots=True)
class Chemistry:
    """One oligonucleotide chemistry, as the rule file declares it."""

    id: str
    action: Action
    length: int
    description: str = ""
    citation: str = ""

    @property
    def cleaves(self) -> bool:
        return self.action is Action.CLEAVES


class AccessibilityModel(Protocol):
    """Whatever can say how exposed a target site is inside the folded transcript.

    Declared and left unimplemented for the same reason the efficiency model in
    the base editor is: this is the factor that decides which ASO works, it
    needs a folding model this package does not carry, and inventing a number
    for it would turn a tiling list into a ranked list that looks like knowledge.
    """

    @property
    def name(self) -> str: ...

    @property
    def availability(self) -> str: ...

    def accessibility(self, start: int, end: int) -> float | None: ...


@dataclass(frozen=True, slots=True)
class NoStructureModel:
    """The default. Says what is missing rather than guessing at it."""

    @property
    def name(self) -> str:
        return "none"

    @property
    def availability(self) -> str:
        return (
            "no structure model is attached, so nothing below is ranked by target "
            "accessibility — which is the factor that most decides whether an antisense "
            "oligonucleotide works. Two windows identical on every rule here can differ "
            "tenfold in a cell because one of them is paired up inside the transcript's own "
            "fold. Folding the target is a job for RNAfold or its equivalents, and none is "
            "running here"
        )

    def accessibility(self, start: int, end: int) -> float | None:
        return None


class Region(StrEnum):
    """Where a window sits relative to the annotated exon."""

    ACCEPTOR_SITE = "acceptor site"
    DONOR_SITE = "donor site"
    EXON_INTERIOR = "exon interior"
    INTRON = "intron"
    UNANNOTATED = "unannotated"


@dataclass(frozen=True, slots=True)
class Exon:
    """One exon's genomic bounds, and how far into the intron a splice site reaches."""

    start: int
    end: int
    #: How many bases either side of a boundary count as the splice signal. The
    #: canonical acceptor and donor motifs are shorter than this; the number is
    #: the practical footprint an oligonucleotide needs to cover to occupy one.
    splice_site_nt: int = 12

    def region_for(self, start: int, end: int) -> Region:
        covers_acceptor = start <= self.start + 2 and end >= self.start - self.splice_site_nt
        covers_donor = end >= self.end - 2 and start <= self.end + self.splice_site_nt
        if covers_acceptor:
            return Region.ACCEPTOR_SITE
        if covers_donor:
            return Region.DONOR_SITE
        if start >= self.start and end <= self.end:
            return Region.EXON_INTERIOR
        return Region.INTRON


@dataclass(frozen=True, slots=True)
class AsoCandidate:
    """One window, the oligonucleotide against it, and what is wrong with it."""

    chemistry: Chemistry
    sequence: str
    chromosome: str
    #: 1-based inclusive genomic span of the *target*, not of the oligonucleotide,
    #: which is its reverse complement and has no coordinates of its own.
    span: tuple[int, int]
    target: str
    gc_fraction: float
    melting_temperature_c: float
    cpg_count: int
    region: Region = Region.UNANNOTATED
    flags: tuple[Flag, ...] = ()

    @property
    def severity(self) -> Severity | None:
        return worst_of(self.flags)

    @property
    def is_blocked(self) -> bool:
        return self.severity is Severity.BLOCKING

    def describe(self) -> str:
        lines = [
            f"{self.sequence}  {self.chromosome}:{self.span[0]}-{self.span[1]}  "
            f"{self.chemistry.id}, {self.region}",
            f"    {len(self.sequence)} nt  GC {self.gc_fraction:.0%}  "
            f"Tm ~{self.melting_temperature_c:.0f} °C  CpG {self.cpg_count}",
        ]
        lines.extend(f"    ! {flag.describe()}" for flag in self.flags)
        return "\n".join(lines)


@dataclass(frozen=True, slots=True)
class AsoOutcome:
    gene: str
    chromosome: str
    span: tuple[int, int]
    chemistry: Chemistry
    candidates: tuple[AsoCandidate, ...] = ()
    refusals: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()
    ruleset_pin: str = ""
    ranking: str = ""
    tiled: int = 0

    @property
    def has_candidates(self) -> bool:
        return bool(self.candidates)

    @property
    def usable(self) -> tuple[AsoCandidate, ...]:
        return tuple(candidate for candidate in self.candidates if not candidate.is_blocked)


def chemistries(rules: FlagRuleset) -> dict[str, Chemistry]:
    """The chemistries the rule file declares, by id."""
    declared = rules.extra.get("chemistries") or []
    if not declared:
        raise DesignError(f"{rules.pin}: the ASO rule file declares no chemistry")

    catalogue: dict[str, Chemistry] = {}
    for entry in declared:
        try:
            action = Action(entry["action"])
        except (KeyError, ValueError) as error:
            raise DesignError(
                f"{rules.pin}: chemistry {entry.get('id')!r} must declare an action of "
                f"{' or '.join(item.value for item in Action)} — an oligonucleotide that "
                "cleaves and one that blocks are not interchangeable"
            ) from error
        catalogue[str(entry["id"])] = Chemistry(
            id=str(entry["id"]),
            action=action,
            length=int(entry["length"]),
            description=" ".join(str(entry.get("description", "")).split()),
            citation=str(entry.get("citation", "")),
        )
    return catalogue


def melting_temperature(sequence: str) -> float:
    """Wallace again, and approximate for the same reasons as in ``prime``."""
    sequence = sequence.upper()
    at = sequence.count("A") + sequence.count("T")
    gc = sequence.count("G") + sequence.count("C")
    return 2.0 * at + 4.0 * gc


def longest_self_complement(sequence: str) -> int:
    """The longest stretch that could pair with another part of the same molecule.

    A crude proxy — it finds the longest substring whose reverse complement also
    occurs, which is a hairpin's necessary condition and not its sufficient one.
    It is here to catch the obvious cases, and the rule file says so.
    """
    sequence = sequence.upper()
    for length in range(len(sequence) // 2, 2, -1):
        for start in range(len(sequence) - length + 1):
            window = sequence[start : start + length]
            complement = reverse_complement(window)
            if complement in sequence[:start] or complement in sequence[start + length :]:
                return length
    return 0


def tile(
    gene: str,
    chromosome: str,
    start: int,
    target: str,
    rules: FlagRuleset,
    *,
    chemistry: str,
    exon: Exon | None = None,
    must_cover: int | None = None,
    model: AccessibilityModel | None = None,
) -> AsoOutcome:
    """Every window of one chemistry's length along a target region.

    ``target`` is the genomic sequence of the region, sense strand, beginning at
    ``start``. Sense strand rather than transcript strand on purpose: the caller
    knows the gene's orientation, and an ASO for a minus-strand gene is the
    *forward* sequence of that region — which is the mistake worth making
    impossible rather than documenting.

    ``must_cover`` drops every window that does not span one coordinate, and it
    exists for allele-specific silencing. The only thing distinguishing the
    mutant transcript from the normal one is the variant itself, so an
    oligonucleotide that does not cover it cannot discriminate — it knocks down
    both alleles, which for a dominant-negative variant removes the good product
    along with the bad. Whether one that *does* cover it discriminates enough is
    a different question, and not one this package answers.
    """
    model = model or NoStructureModel()
    catalogue = chemistries(rules)
    if chemistry not in catalogue:
        raise DesignError(
            f"{chemistry!r} is not a chemistry in {rules.pin}; it declares "
            f"{', '.join(sorted(catalogue))}"
        )
    selected = catalogue[chemistry]
    target = target.upper()

    if len(target) < selected.length:
        raise DesignError(
            f"the target region is {len(target)} nt and {selected.id} is {selected.length} nt "
            "long — there is no window to tile"
        )

    step = int(rules.threshold("step_nt", 1))
    homopolymer = int(rules.threshold("homopolymer_run_nt", 4))
    guanine_run = int(rules.threshold("guanine_run_nt", 4))

    candidates: list[AsoCandidate] = []
    tiled = 0
    for offset in range(0, len(target) - selected.length + 1, step):
        window = target[offset : offset + selected.length]
        tiled += 1
        if not is_resolved(window):
            continue

        span = (start + offset, start + offset + selected.length - 1)
        if must_cover is not None and not span[0] <= must_cover <= span[1]:
            continue
        sequence = reverse_complement(window)
        region = exon.region_for(*span) if exon else Region.UNANNOTATED

        features = FlatFeatures(
            values={
                "aso.length_nt": len(sequence),
                "aso.gc_fraction": _gc(sequence),
                "aso.tm_c": melting_temperature(sequence),
                "aso.cpg_count": sequence.count("CG"),
                "aso.has_guanine_run": "G" * guanine_run in sequence,
                "aso.has_homopolymer_run": any(
                    base * homopolymer in sequence for base in "ACGT"
                ),
                "aso.self_complement_nt": longest_self_complement(sequence),
                "aso.chemistry": selected.id,
                "aso.chemistry_cleaves": selected.cleaves,
                "aso.chemistry_blocks": not selected.cleaves,
                "aso.region": region.value,
                "aso.covers_splice_site": region in {Region.ACCEPTOR_SITE, Region.DONOR_SITE},
                "aso.exon_annotated": exon is not None,
                "aso.accessibility": model.accessibility(*span),
            }
        )

        candidates.append(
            AsoCandidate(
                chemistry=selected,
                sequence=sequence,
                chromosome=chromosome,
                span=span,
                target=window,
                gc_fraction=_gc(sequence),
                melting_temperature_c=melting_temperature(sequence),
                cpg_count=sequence.count("CG"),
                region=region,
                flags=rules.raise_flags(features),
            )
        )

    ordered = sorted(
        candidates,
        key=lambda candidate: (
            sort_weight(candidate.severity),
            len(candidate.flags),
            candidate.span[0],
        ),
    )

    notes: list[str] = []
    if must_cover is not None:
        notes.append(
            f"only windows covering {chromosome}:{must_cover} were kept. That base is the one "
            "thing telling the two transcripts apart, so an oligonucleotide missing it silences "
            "the healthy allele as well as the affected one. Whether covering it is enough to "
            "discriminate — a single mismatch in a 20-mer often is not — is not answered here"
        )
    if exon is None:
        notes.append(
            "no exon was supplied, so no window below is placed against a splice boundary. For "
            "a steric blocker that is the whole question — it does nothing unless it covers "
            "something worth covering — and for a gapmer it is the difference between cutting "
            "an intron and cutting the message"
        )

    return AsoOutcome(
        gene=gene,
        chromosome=chromosome,
        span=(start, start + len(target) - 1),
        chemistry=selected,
        candidates=tuple(ordered),
        refusals=(),
        notes=tuple(notes),
        ruleset_pin=rules.pin,
        ranking=model.availability,
        tiled=tiled,
    )


def _gc(sequence: str) -> float:
    if not sequence:
        return 0.0
    return round((sequence.count("G") + sequence.count("C")) / len(sequence), 3)
