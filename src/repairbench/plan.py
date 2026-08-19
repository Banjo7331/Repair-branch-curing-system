"""One case, end to end: mechanism, then modality, then the molecule.

Until this module existed the package had three good halves and no seam. You
could run ``assess`` and learn that gene addition is contraindicated, then run
``design`` on the same variant and get a page of protospacers, and nothing
anywhere connected the two. The commands were independent, so the safety
property M6 exists for — *a contraindication outranks any number of
indications* — stopped at the edge of M6.

This is the seam, and its one rule is the one that could not be enforced before:

**A modality M6 ruled out is never designed.** Not designed and marked, not
designed with a warning at the bottom — not designed. A page of pegRNAs for a
contraindicated intervention is worse than no output at all, because a sequence
is a thing somebody can order and a caveat is a thing somebody can skim.

Two weaker rules follow it. An unresolved mechanism designs nothing, because
every modality below it is resting on transcript facts alone. And a modality
that was merely *not indicated* is not designed either, but is listed with its
verdict rather than omitted — the reader should be able to see that it was
considered.

What comes out is one document: why the variant causes disease, which classes
of intervention that admits, and — for those with a designer in this package —
the actual molecules, each with what is wrong with it. Every rule file that
contributed is pinned at the bottom, because a plan that cannot name the six
files it was produced under cannot be reproduced or compared with a later one.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from repairbench.annotation.fasta import SequenceProvider
from repairbench.annotation.store import TranscriptStore
from repairbench.design.aso import AsoOutcome, Exon
from repairbench.design.aso import tile as tile_asos
from repairbench.design.candidate import DesignOutcome
from repairbench.design.designer import CorrectionRequest
from repairbench.design.designer import design as design_base_edits
from repairbench.design.editors import DesignError, EditorCatalogue
from repairbench.design.flags import FlagRuleset
from repairbench.design.prime import EditRequest, PrimeOutcome, design_pegrnas
from repairbench.features import MechanismQuery
from repairbench.modality import Modality, ModalitySelection, Verdict
from repairbench.modality_rules import ModalityRuleset
from repairbench.model import MechanismCall, RepairbenchError
from repairbench.ruleset import Ruleset, RulesetError

Outcome = DesignOutcome | PrimeOutcome | AsoOutcome


class PlanError(RepairbenchError):
    """A plan cannot be assembled from what was supplied."""


@dataclass(frozen=True, slots=True)
class Route:
    """Where one modality goes, and why it goes there."""

    modality: Modality
    designer: str
    because: str
    chemistry: str = ""
    target: str = ""

    @property
    def has_designer(self) -> bool:
        return self.designer != "none"


@dataclass(frozen=True, slots=True)
class RoutingTable:
    """The seam between M6 and M7, pinned like every other rule file."""

    version: str
    description: str
    routes: dict[Modality, Route]
    digest: str
    thresholds: dict[str, Any] = field(default_factory=dict)

    @property
    def pin(self) -> str:
        return f"{self.version}@{self.digest[:12]}"

    def threshold(self, name: str, default: Any = None) -> Any:
        return self.thresholds.get(name, default)

    def route_for(self, modality: Modality) -> Route:
        route = self.routes.get(modality)
        if route is None:
            raise PlanError(
                f"{self.pin} does not say what to do with {modality}. Every modality must be "
                "listed, including the ones with no designer — an unlisted one is a modality "
                "nobody decided about, and it would silently produce nothing"
            )
        return route


def load_routing(path: str | Path) -> RoutingTable:
    """Read and validate the routing table.

    Refuses a table that does not mention every modality. A missing entry would
    look exactly like "no candidates found" in the output, which is the failure
    mode this whole file exists to prevent.
    """
    raw = Path(path).read_bytes()
    document = yaml.safe_load(raw)
    if not isinstance(document, dict) or "designers" not in document:
        raise RulesetError(f"{path}: routing file must be a mapping with a 'designers' list")

    routes: dict[Modality, Route] = {}
    for index, entry in enumerate(document["designers"], start=1):
        try:
            modality = Modality(entry["modality"])
        except (KeyError, ValueError) as error:
            raise RulesetError(
                f"{path}: entry {index} names an unknown modality {entry.get('modality')!r}"
            ) from error
        routes[modality] = Route(
            modality=modality,
            designer=str(entry.get("designer", "none")),
            because=" ".join(str(entry.get("because", "")).split()),
            chemistry=str(entry.get("chemistry", "")),
            target=str(entry.get("target", "")),
        )

    missing = [modality.value for modality in Modality if modality not in routes]
    if missing:
        raise RulesetError(f"{path}: no routing entry for {', '.join(missing)}")

    return RoutingTable(
        version=str(document["version"]),
        description=str(document.get("description", "")),
        routes=routes,
        digest=hashlib.sha256(raw).hexdigest(),
        thresholds=dict(document.get("thresholds") or {}),
    )


@dataclass(frozen=True, slots=True)
class DesignedFor:
    """What happened to one modality on the way to a molecule."""

    modality: Modality
    verdict: Verdict
    designer: str
    outcome: Outcome | None = None
    #: Why nothing was designed, when nothing was. Never empty when ``outcome``
    #: is ``None``, because "no output" is the sentence a reader most needs
    #: explained.
    refusal: str = ""

    @property
    def designed(self) -> bool:
        return self.outcome is not None


@dataclass(frozen=True, slots=True)
class Plan:
    """One case, from why to what."""

    gene: str
    call: MechanismCall
    selection: ModalitySelection
    designs: tuple[DesignedFor, ...] = ()
    pins: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()

    @property
    def designed(self) -> tuple[DesignedFor, ...]:
        return tuple(design for design in self.designs if design.designed)

    @property
    def ruled_out(self) -> tuple[DesignedFor, ...]:
        return tuple(
            design for design in self.designs if design.verdict is Verdict.CONTRAINDICATED
        )

    @property
    def has_designs(self) -> bool:
        return bool(self.designed)


@dataclass(frozen=True, slots=True)
class Designers:
    """The rule files each designer runs under, gathered so a plan can pin them."""

    editors: EditorCatalogue
    prime: FlagRuleset
    aso: FlagRuleset
    routing: RoutingTable


@dataclass(frozen=True, slots=True)
class Locus:
    """Where the variant is, in the terms the designers need.

    Separate from ``MechanismQuery`` on purpose. M5 and M6 reason in CDS
    offsets, which is right for them and useless for a designer: a protospacer
    is placed in genomic coordinates, against reference sequence. Assembling
    this is where a plan finds out that the case file never carried a genomic
    coordinate at all, and that is worth its own refusal.
    """

    chromosome: str
    position: int
    reference: str
    alternate: str
    #: Genomic bounds of the coding exon the variant sits in, when an annotation
    #: was supplied. The exon-skipping designer needs it and nothing else does.
    exon: tuple[int, int] | None = None
    #: Which strand the gene is transcribed from. ``None`` without an
    #: annotation, and that is a refusal rather than a default: an antisense
    #: oligonucleotide copies one genomic strand for a plus-strand gene and the
    #: other for a minus-strand one, so guessing produces a molecule with the
    #: sequence of its own target.
    strand: str | None = None

    @classmethod
    def from_case(cls, case: dict[str, Any], store: TranscriptStore | None) -> Locus | None:
        genomic = case.get("genomic")
        if not genomic:
            return None
        for required in ("chromosome", "position", "reference", "alternate"):
            if required not in genomic:
                raise PlanError(
                    f"the genomic block has no {required!r}. A designer needs the reference and "
                    "alternate alleles as written on the genome, not only a CDS offset"
                )

        chromosome = str(genomic["chromosome"])
        position = int(genomic["position"])
        exon = None
        if store is not None:
            gene = case["variant"]["gene"]
            record, _ = store.preferred_for(gene)
            exon = next(
                (
                    (start, end)
                    for start, end in record.cds_blocks
                    if start <= position <= end
                ),
                None,
            )
            strand = record.strand
        else:
            strand = None
        return cls(
            chromosome=chromosome,
            position=position,
            reference=str(genomic["reference"]).upper(),
            alternate=str(genomic["alternate"]).upper(),
            exon=exon,
            strand=strand,
        )


def plan(
    query: MechanismQuery,
    call: MechanismCall,
    selection: ModalitySelection,
    designers: Designers,
    *,
    locus: Locus | None = None,
    sequences: SequenceProvider | None = None,
    mechanism_rules: Ruleset | None = None,
    modality_rules: ModalityRuleset | None = None,
) -> Plan:
    """Take a finished assessment through to the molecules it admits."""
    designs: list[DesignedFor] = []
    notes: list[str] = []

    if selection.is_blocked:
        notes.append(
            "the mechanism is unresolved, so nothing below was designed. Every modality "
            "depends on it, and a molecule designed against an undetermined mechanism is a "
            "molecule designed against nothing"
        )

    for assessment in selection.assessments:
        route = designers.routing.route_for(assessment.modality)
        designs.append(
            _design_one(
                assessment.modality,
                assessment.verdict,
                route,
                designers,
                locus,
                sequences,
                query.gene.symbol,
            )
        )

    if locus is None:
        notes.append(
            "the case file carries no genomic coordinate, so no molecule could be designed "
            "for any modality. M5 and M6 reason in CDS offsets; a protospacer is placed on the "
            "genome, and the two are not interconvertible without an annotation"
        )
    elif sequences is None:
        notes.append(
            "no reference sequence was supplied, so nothing was designed. Every designer here "
            "reads bases"
        )

    pins = [
        f"mechanism  {mechanism_rules.pin if mechanism_rules else call.ruleset_version}",
        f"modality   {modality_rules.pin if modality_rules else selection.ruleset_version}",
        f"editors    {designers.editors.pin}",
        f"prime      {designers.prime.pin}",
        f"aso        {designers.aso.pin}",
        f"routing    {designers.routing.pin}",
    ]

    return Plan(
        gene=query.gene.symbol,
        call=call,
        selection=selection,
        designs=tuple(designs),
        pins=tuple(pins),
        notes=tuple(notes),
    )


def _design_one(  # noqa: PLR0911 — one return per reason nothing was designed, which
    # is the point: collapsing them would produce one refusal text for six
    # different situations, and the text is the part a reader acts on.
    modality: Modality,
    verdict: Verdict,
    route: Route,
    designers: Designers,
    locus: Locus | None,
    sequences: SequenceProvider | None,
    gene: str,
) -> DesignedFor:
    """One modality, and the several reasons it may produce no molecule.

    The order of these checks is the safety property. The contraindication is
    tested first and returns first, so no combination of missing inputs, absent
    designers or later failures can route around it.
    """
    if verdict is Verdict.CONTRAINDICATED:
        return DesignedFor(
            modality=modality,
            verdict=verdict,
            designer=route.designer,
            refusal=(
                "ruled out by the modality rules, so nothing was designed. A sequence is "
                "something a reader can order and a caveat is something a reader can skim, "
                "which is why this is a refusal rather than a warning"
            ),
        )

    if verdict is Verdict.BLOCKED_BY_UNRESOLVED_MECHANISM:
        return DesignedFor(
            modality=modality,
            verdict=verdict,
            designer=route.designer,
            refusal="the mechanism is unresolved, so no modality was assessed and none designed",
        )

    if verdict is not Verdict.INDICATED:
        return DesignedFor(
            modality=modality,
            verdict=verdict,
            designer=route.designer,
            refusal=(
                "not indicated by the modality rules. Nothing here says it is impossible — "
                "only that no rule argued for it, which is a different claim and a weaker one"
            ),
        )

    if not route.has_designer:
        return DesignedFor(
            modality=modality, verdict=verdict, designer=route.designer, refusal=route.because
        )

    if locus is None or sequences is None:
        return DesignedFor(
            modality=modality,
            verdict=verdict,
            designer=route.designer,
            refusal=(
                "indicated, and not designed: this run had no "
                f"{'genomic coordinate' if locus is None else 'reference sequence'} to design "
                "against"
            ),
        )

    try:
        outcome = _run_designer(route, designers, locus, sequences, gene)
    except (DesignError, PlanError) as error:
        return DesignedFor(
            modality=modality,
            verdict=verdict,
            designer=route.designer,
            refusal=f"the {route.designer} designer declined: {error}",
        )
    return DesignedFor(
        modality=modality, verdict=verdict, designer=route.designer, outcome=outcome
    )


def _run_designer(
    route: Route,
    designers: Designers,
    locus: Locus,
    sequences: SequenceProvider,
    gene: str,
) -> Outcome:
    """Call the designer this modality routes to, in its own terms.

    The allele convention is worth stating once: the patient carries the
    *alternate* allele and the *reference* is what it should read. Swapping them
    produces a design that installs the disease.
    """
    if route.designer == "base_editor":
        return design_base_edits(
            CorrectionRequest(
                gene=gene,
                chromosome=locus.chromosome,
                position=locus.position,
                patient_base=locus.alternate,
                wild_type_base=locus.reference,
            ),
            sequences,
            designers.editors,
        )

    if route.designer == "pegrna":
        return design_pegrnas(
            EditRequest(
                gene=gene,
                chromosome=locus.chromosome,
                position=locus.position,
                patient_allele=locus.alternate,
                wild_type_allele=locus.reference,
            ),
            sequences,
            designers.prime,
        )

    if route.designer == "aso":
        if locus.strand is None:
            raise PlanError(
                "no annotation was supplied, so the gene's strand is unknown and no "
                "antisense oligonucleotide can be designed. Which genomic strand the "
                "molecule copies is decided by the gene's orientation, and a guess "
                "produces one with the sequence of the transcript it was meant to bind."
            )
        start, end, exon = _aso_target(route, designers.routing, locus)
        region = sequences.fetch(locus.chromosome, start, end)
        discriminating = route.target != "affected_exon"
        if discriminating:
            # Tiled against the *patient's* sequence, not the reference. An
            # oligonucleotide complementary to the reference base is
            # complementary to the healthy transcript — which is the allele it
            # was supposed to spare.
            offset = locus.position - start
            region = region[:offset] + locus.alternate + region[offset + len(locus.reference) :]
        return tile_asos(
            gene,
            locus.chromosome,
            start,
            region,
            designers.aso,
            chemistry=route.chemistry,
            strand=locus.strand,
            exon=exon,
            must_cover=locus.position if discriminating else None,
        )

    raise PlanError(
        f"{designers.routing.pin} names a designer this package does not have: "
        f"{route.designer!r}"
    )


def _aso_target(
    route: Route, routing: RoutingTable, locus: Locus
) -> tuple[int, int, Exon | None]:
    """Which stretch of sequence to tile, which is a different answer per modality.

    Skipping and splice correction aim at the exon and its boundaries, so the
    target is the exon plus enough flanking intron to put the acceptor and donor
    sites inside the tiled region rather than at its edge.

    Allele-specific silencing aims at the variant itself, and this is the honest
    part: the only thing distinguishing the mutant transcript from the normal one
    is that base, so an oligonucleotide that does not cover it cannot
    discriminate. Whether one that does cover it discriminates *enough* is a
    question this package does not answer.
    """
    if route.target == "affected_exon":
        if locus.exon is None:
            raise PlanError(
                "no annotation placed this variant in a coding exon, so there is no exon to "
                "tile across. Skipping is defined by which exon comes out, and guessing at its "
                "boundaries would put the oligonucleotide over the wrong splice site"
            )
        flank = int(routing.threshold("splice_flank_nt", 30))
        start, end = locus.exon
        return max(1, start - flank), end + flank, Exon(start=start, end=end)

    window = int(routing.threshold("variant_window_nt", 60))
    return max(1, locus.position - window), locus.position + window, None
