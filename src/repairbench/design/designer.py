"""Placing a protospacer so the base to be corrected lands in an editor's window.

The arithmetic is short and every step of it is a way to be wrong quietly.

**The design is against the patient's sequence, not the reference.** The base to
be corrected is by definition the one the reference does not have. Scanning the
reference window would look for a PAM around a base that is not there — and
where the variant itself creates or destroys a PAM, the reference scan produces
guides that do not exist in this patient, or misses the ones that do.

**Which strand is decided by the conversion, not by preference.** A deaminase
makes A→G and C→T and nothing else. Correcting a patient's T to a C is not a
C-editor's job — it is an A-editor's job on the *other* strand, where the same
base pair reads A and needs to read G. Half the correctable variants in a
genome are only correctable that way, and a designer that scans one strand
misses them and reports "no candidates" instead of "look at the minus strand".

**Position numbering runs from the PAM-distal end.** Protospacer position 1 is
the far end and 20 sits against the PAM, which is the convention every window in
the literature is quoted in. Numbering the other way would put every window at
the wrong end of every guide and still produce plausible-looking output.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, replace

from repairbench.annotation.fasta import SequenceProvider
from repairbench.design.candidate import Bystander, DesignOutcome, EditCandidate
from repairbench.design.editors import Conversion, DesignError, Editor, EditorCatalogue
from repairbench.design.efficiency import EfficiencyModel, NoModelAttached, ordered
from repairbench.design.sequence import (
    complement_base,
    is_resolved,
    matches_pam,
    reverse_complement,
)
from repairbench.model import Zygosity

PLUS, MINUS = "+", "-"


@dataclass(frozen=True, slots=True)
class CorrectionRequest:
    """One base to put back, and everything needed to look for a way to do it."""

    gene: str
    chromosome: str
    #: 1-based genomic coordinate of the base to correct.
    position: int
    #: What the patient has there. The allele the disease is attributed to.
    patient_base: str
    #: What it should read instead. Usually the reference base, and named
    #: separately because "reference" and "wild type" are not always the same
    #: claim — a reference allele can itself be the minor one.
    wild_type_base: str
    zygosity: Zygosity = Zygosity.UNKNOWN

    def __post_init__(self) -> None:
        for name, base in (("patient", self.patient_base), ("wild type", self.wild_type_base)):
            if len(base) != 1 or base.upper() not in "ACGT":
                raise DesignError(
                    f"{name} allele {base!r} is not a single unambiguous base. Base editing "
                    "corrects single-base substitutions; an insertion, a deletion or a "
                    "multi-base change is prime editing's problem, and this module will not "
                    "pretend to address it"
                )
        if self.patient_base.upper() == self.wild_type_base.upper():
            raise DesignError(
                f"{self.gene} {self.chromosome}:{self.position}: the patient base and the "
                "wild-type base are the same — there is nothing to correct"
            )


def design(
    request: CorrectionRequest,
    sequences: SequenceProvider,
    catalogue: EditorCatalogue,
    *,
    model: EfficiencyModel | None = None,
    coding: Callable[[int], bool] | None = None,
) -> DesignOutcome:
    """Every protospacer in the catalogue that puts this base in a window.

    ``coding`` — if supplied — answers whether a genomic position sits in coding
    sequence, and is used only to annotate bystanders. It is a callable rather
    than a transcript so that this module does not acquire an opinion about
    which transcript matters; the caller has already decided that.
    """
    model = model or NoModelAttached()
    patient, wild_type = request.patient_base.upper(), request.wild_type_base.upper()

    routes = _routes(patient, wild_type)
    if not routes:
        return DesignOutcome(
            gene=request.gene,
            chromosome=request.chromosome,
            position=request.position,
            patient_base=patient,
            wild_type_base=wild_type,
            refusals=(
                f"{patient}→{wild_type} is a transversion, and no deaminase makes one. Base "
                "editing converts A→G and C→T on the strand it acts on, which covers the four "
                "transitions and nothing else. This correction needs prime editing or a "
                "template-driven repair, neither of which this module designs",
            ),
            catalogue_pin=catalogue.pin,
        )

    padding = catalogue.thresholds.search_padding_nt
    window_start = max(1, request.position - padding)
    window_end = request.position + padding
    reference = sequences.fetch(request.chromosome, window_start, window_end).upper()

    observed = reference[request.position - window_start]
    notes: list[str] = []
    if observed != wild_type:
        notes.append(
            f"the reference reads {observed} at {request.chromosome}:{request.position} while the "
            f"correction targets {wild_type}. That is not necessarily wrong — a reference allele "
            "is not always the wild-type one — but a mismatch here is also what a wrong assembly "
            "or a wrong coordinate looks like"
        )

    # The one substitution that makes this the patient's genome rather than
    # anybody's. Everything below is scanned against this string.
    index = request.position - window_start
    patient_sequence = reference[:index] + patient + reference[index + 1 :]

    candidates: list[EditCandidate] = []
    refusals: list[str] = []
    considered: list[str] = []

    for strand, conversion in routes.items():
        editors = catalogue.making(conversion)
        if not editors:
            refusals.append(
                f"correcting this base on the {strand} strand needs a {conversion} editor, and "
                f"the catalogue ({catalogue.pin}) has none"
            )
            continue
        for editor in editors:
            considered.append(f"{editor.id} on the {strand} strand")
            found = _scan(request, editor, patient_sequence, window_start, window_end, strand)
            candidates.extend(found)

    if not candidates and not refusals:
        refusals.append(
            f"no PAM in the catalogue sits at a distance that puts g.{request.position} inside "
            f"an editing window, within {padding} nt either side. A relaxed-PAM nuclease is the "
            "usual next thing to try, and the catalogue's SpRY entry says what that costs"
        )

    annotated = tuple(_annotate(candidate, request, coding) for candidate in candidates)
    notes.extend(_notes_for(request))

    return DesignOutcome(
        gene=request.gene,
        chromosome=request.chromosome,
        position=request.position,
        patient_base=patient,
        wild_type_base=wild_type,
        candidates=ordered(annotated, model),
        refusals=tuple(refusals),
        notes=tuple(notes),
        catalogue_pin=catalogue.pin,
        ranking=model.availability,
        considered=tuple(considered),
    )


def _routes(patient: str, wild_type: str) -> dict[str, Conversion]:
    """Which strand admits which conversion.

    On the minus strand the same base pair reads as its complement, so a
    patient's T needing to become C is an A needing to become G — the commonest
    correction there is, and invisible to anything that only looks at the plus
    strand.
    """
    routes: dict[str, Conversion] = {}
    forward = Conversion.between(patient, wild_type)
    if forward is not None:
        routes[PLUS] = forward
    reverse = Conversion.between(complement_base(patient), complement_base(wild_type))
    if reverse is not None:
        routes[MINUS] = reverse
    return routes


def _scan(
    request: CorrectionRequest,
    editor: Editor,
    patient_sequence: str,
    window_start: int,
    window_end: int,
    strand: str,
) -> list[EditCandidate]:
    """Every PAM on one strand that places the target in this editor's window.

    Minus-strand scanning is done by reverse-complementing the whole window and
    running the identical plus-strand logic over it, with one function mapping
    an index back to a genomic coordinate. Two sets of index arithmetic for the
    two strands is where off-by-ones live.
    """
    if strand == PLUS:
        strand_sequence = patient_sequence
        def to_genomic(index: int) -> int:
            return window_start + index
    else:
        strand_sequence = reverse_complement(patient_sequence)
        def to_genomic(index: int) -> int:
            return window_end - index

    target_index = (
        request.position - window_start if strand == PLUS else window_end - request.position
    )
    length, pam_length = editor.protospacer_length, len(editor.pam)
    found: list[EditCandidate] = []

    for pam_index in range(length, len(strand_sequence) - pam_length + 1):
        pam = strand_sequence[pam_index : pam_index + pam_length]
        if not matches_pam(pam, editor.pam):
            continue

        start = pam_index - length
        position_in_protospacer = target_index - start + 1
        if position_in_protospacer not in editor.window:
            continue

        protospacer = strand_sequence[start:pam_index]
        if not is_resolved(protospacer):
            continue

        bystanders = tuple(
            Bystander(
                position_in_protospacer=offset + 1,
                genomic_position=to_genomic(start + offset),
                becomes=editor.conversion.product_base,
            )
            for offset in range(editor.window_start - 1, editor.window_end)
            if strand_sequence[start + offset] == editor.conversion.source_base
            and start + offset != target_index
        )

        span = sorted((to_genomic(start), to_genomic(pam_index - 1)))
        pam_span = sorted((to_genomic(pam_index), to_genomic(pam_index + pam_length - 1)))
        found.append(
            EditCandidate(
                editor=editor,
                chromosome=request.chromosome,
                strand=strand,
                protospacer=protospacer,
                pam=pam,
                span=(span[0], span[1]),
                pam_span=(pam_span[0], pam_span[1]),
                target_position_in_protospacer=position_in_protospacer,
                target_genomic_position=request.position,
                bystanders=bystanders,
            )
        )
    return found


def _annotate(
    candidate: EditCandidate,
    request: CorrectionRequest,
    coding: Callable[[int], bool] | None,
) -> EditCandidate:
    """Attach the warnings that belong to a candidate rather than to a run."""
    warnings: list[str] = []

    if len(candidate.bystanders) >= 2:
        warnings.append(
            f"{len(candidate.bystanders)} other editable bases sit in the window; each is an "
            "unintended change, and this package does not predict what any of them does to the "
            "protein"
        )

    edge = candidate.target_position_in_protospacer in {
        candidate.editor.window_start,
        candidate.editor.window_end,
    }
    if edge:
        warnings.append(
            f"the target sits at position {candidate.target_position_in_protospacer}, the edge of "
            f"this editor's stated window ({candidate.editor.window_start}-"
            f"{candidate.editor.window_end}) — the position where the window is least certain"
        )

    if coding is None:
        return replace(candidate, warnings=tuple(warnings))

    bystanders = tuple(
        replace(bystander, in_coding_sequence=coding(bystander.genomic_position))
        for bystander in candidate.bystanders
    )
    if any(bystander.in_coding_sequence for bystander in bystanders):
        warnings.append(
            "a bystander above is in coding sequence, where a silent change and a missense "
            "change look identical from here"
        )
    return replace(candidate, bystanders=bystanders, warnings=tuple(warnings))


def _notes_for(request: CorrectionRequest) -> list[str]:
    """What is true of the run rather than of any one candidate.

    The heterozygote note is the reason this function exists. Attached to each
    candidate it was repeated six times in one report, which is how a real
    caution gets skimmed past; it is a fact about the patient, so it belongs
    once, at the top.
    """
    if not request.zygosity.leaves_a_wild_type_allele:
        return []
    return [
        "the patient has an unaffected copy, and the one base that distinguishes the two "
        "alleles is the target itself — which sits in the PAM-distal half of the protospacer, "
        "where a single mismatch discriminates poorly. Every guide below should be assumed to "
        "bind both alleles. The intended edit is harmless on the healthy one, since that base "
        "already reads correctly, but every bystander applies to it as well"
    ]
