"""Designing pegRNAs: the edit written into the guide that installs it.

Prime editing is the modality base editing refuses to be. A deaminase makes four
transitions; a prime editor writes whatever the template says — every
substitution, and small insertions and deletions besides. What it costs is a
design space large enough that picking badly out of it is the normal outcome.

The geometry, once, because everything below is arithmetic over it:

* A Cas9 nickase cuts **one** strand — the one the protospacer matches — three
  bases 5' of the PAM, between protospacer positions 17 and 18.
* The 3' end that cut creates is a primer. The **primer binding site** at the
  far end of the pegRNA anneals to it, which is why the PBS is the reverse
  complement of the bases *immediately upstream* of the nick.
* Reverse transcriptase then copies the **reverse transcription template**
  backwards off the pegRNA, writing new sequence 3' of the nick. The edit lives
  in that template, so the edit must sit **downstream of the nick** — and a
  protospacer whose nick falls past the edit is useless no matter how good its
  PAM is. That single constraint eliminates about half of the PAMs near any
  given variant, and it is the first thing a naive implementation misses.
* The template must continue past the edit far enough for the new flap to
  anneal against the genome. That tail is the homology arm.

Two consequences worth stating because they are where a first implementation
goes wrong. The pegRNA extension is written 5'→3' as *template then primer* —
RTT first, PBS last — which is the reverse of the order they are used in and the
reverse of the order they are usually explained in. And both are reverse
complements of the strand the protospacer matches, so a minus-strand protospacer
means two complementations rather than one; getting either wrong produces a
sequence of exactly the right length that installs nothing.
"""

from __future__ import annotations

from dataclasses import dataclass

from repairbench.annotation.fasta import SequenceProvider
from repairbench.design.editors import DesignError
from repairbench.design.flags import (
    Flag,
    FlagRuleset,
    FlatFeatures,
    Severity,
    sort_weight,
    worst_of,
)
from repairbench.design.sequence import is_resolved, matches_pam, reverse_complement

PLUS, MINUS = "+", "-"


@dataclass(frozen=True, slots=True)
class EditRequest:
    """What should be at a coordinate instead of what is.

    Written as two alleles rather than as a "change type", because prime editing
    does not care which of the three it is: a substitution is a one-for-one
    template, an insertion a longer one, a deletion a shorter one, and the
    arithmetic below is identical. Both alleles are written VCF-style with a
    shared anchor base, so neither is ever empty.
    """

    gene: str
    chromosome: str
    #: 1-based coordinate of the first base of the patient's allele.
    position: int
    patient_allele: str
    wild_type_allele: str

    def __post_init__(self) -> None:
        for name, allele in (
            ("patient", self.patient_allele),
            ("wild type", self.wild_type_allele),
        ):
            if not allele or not is_resolved(allele.upper()):
                raise DesignError(
                    f"{name} allele {allele!r} is not unambiguous sequence. Write an insertion "
                    "or a deletion the way a VCF does, with the anchor base included in both "
                    "alleles, rather than as an empty string"
                )
        if self.patient_allele.upper() == self.wild_type_allele.upper():
            raise DesignError(
                f"{self.gene} {self.chromosome}:{self.position}: the two alleles are the same — "
                "there is nothing to correct"
            )

    @property
    def kind(self) -> str:
        if len(self.patient_allele) == len(self.wild_type_allele):
            return "substitution"
        return "insertion" if len(self.wild_type_allele) > len(self.patient_allele) else "deletion"


@dataclass(frozen=True, slots=True)
class NickingGuide:
    """A second nick, for PE3 — and whether it is the PE3b kind."""

    protospacer: str
    pam: str
    strand: str
    nick_position: int
    distance_nt: int
    #: True when this guide's protospacer matches the *edited* sequence and not
    #: the original. It cannot fire until the edit is installed, which is what
    #: makes PE3b produce fewer indels than PE3.
    edit_dependent: bool = False

    def describe(self) -> str:
        kind = "PE3b" if self.edit_dependent else "PE3 "
        return (
            f"{kind}  {self.protospacer} {self.pam.lower()} ({self.strand})  "
            f"nick at g.{self.nick_position}, {self.distance_nt} nt away"
        )


@dataclass(frozen=True, slots=True)
class PegRna:
    """One complete pegRNA: spacer, template, primer, and what is wrong with it."""

    spacer: str
    pam: str
    strand: str
    chromosome: str
    protospacer_span: tuple[int, int]
    nick_position: int
    nick_to_edit_nt: int
    pbs: str
    rtt: str
    homology_arm_nt: int
    pam_disrupted_by_edit: bool = False
    nicking_guides: tuple[NickingGuide, ...] = ()
    flags: tuple[Flag, ...] = ()

    @property
    def extension(self) -> str:
        """The 3' extension as it is ordered: template first, primer last."""
        return f"{self.rtt}{self.pbs}"

    @property
    def severity(self) -> Severity | None:
        return worst_of(self.flags)

    @property
    def is_blocked(self) -> bool:
        return self.severity is Severity.BLOCKING

    @property
    def pe3b(self) -> tuple[NickingGuide, ...]:
        return tuple(guide for guide in self.nicking_guides if guide.edit_dependent)

    def describe(self) -> str:
        lines = [
            f"{self.spacer} {self.pam.lower()}  {self.chromosome}:"
            f"{self.protospacer_span[0]}-{self.protospacer_span[1]} ({self.strand})",
            f"    nick at g.{self.nick_position}, {self.nick_to_edit_nt} nt from the edit"
            f"{'' if self.pam_disrupted_by_edit else '; the PAM survives the edit'}",
            f"    PBS {len(self.pbs):>2} nt  {self.pbs}",
            f"    RTT {len(self.rtt):>2} nt  {self.rtt}  "
            f"({self.homology_arm_nt} nt of homology past the edit)",
            f"    3' extension  {self.extension}",
        ]
        lines.extend(f"      · {guide.describe()}" for guide in self.nicking_guides[:2])
        lines.extend(f"    ! {flag.describe()}" for flag in self.flags)
        return "\n".join(lines)


@dataclass(frozen=True, slots=True)
class PrimeOutcome:
    """Every pegRNA, and what could not be designed."""

    gene: str
    chromosome: str
    position: int
    edit: str
    candidates: tuple[PegRna, ...] = ()
    refusals: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()
    ruleset_pin: str = ""
    ranking: str = ""

    @property
    def has_candidates(self) -> bool:
        return bool(self.candidates)

    @property
    def usable(self) -> tuple[PegRna, ...]:
        """Those nothing blocking was said about. Not "the good ones"."""
        return tuple(candidate for candidate in self.candidates if not candidate.is_blocked)

    @property
    def protospacers(self) -> int:
        """How many distinct placements the candidates came from.

        Reported because the candidate count on its own is misleading: ten
        pegRNAs differing only in primer length are one choice with ten
        settings, not ten choices.
        """
        return len({(candidate.spacer, candidate.strand) for candidate in self.candidates})


def melting_temperature(sequence: str) -> float:
    """The Wallace approximation, and named as an approximation on purpose.

    Two degrees per A or T, four per G or C. It is the rule of thumb for short
    oligonucleotides, it is not nearest-neighbour, and the primer binding site
    is an RNA:DNA hybrid for which neither model is exactly right. It is here
    because the threshold in the rule file was fitted against a similar
    approximation, and computing a more precise number to compare with a rough
    threshold would be false precision rather than more accuracy.
    """
    sequence = sequence.upper()
    at = sequence.count("A") + sequence.count("T")
    gc = sequence.count("G") + sequence.count("C")
    return 2.0 * at + 4.0 * gc


@dataclass(frozen=True, slots=True)
class _Nuclease:
    """The nickase, as the rule file declares it."""

    id: str
    pam: str
    protospacer_length: int
    nick_after_position: int
    citation: str = ""

    @property
    def pam_offset_after_nick(self) -> int:
        """How many bases of new sequence sit between the nick and the PAM.

        Three, for SpCas9 — which is the same statement as "the nick falls three
        bases 5' of the PAM", read from the other side. Derived rather than
        declared so that a rule file moving the nick cannot leave this behind.
        """
        return self.protospacer_length - self.nick_after_position


def _nuclease(rules: FlagRuleset) -> _Nuclease:
    declared = rules.extra.get("nucleases") or []
    if not declared:
        raise DesignError(f"{rules.pin}: the prime rule file declares no nuclease")
    entry = declared[0]
    return _Nuclease(
        id=str(entry["id"]),
        pam=str(entry["pam"]).upper(),
        protospacer_length=int(entry["protospacer_length"]),
        nick_after_position=int(entry["nick_after_position"]),
        citation=str(entry.get("citation", "")),
    )


@dataclass(frozen=True, slots=True)
class _Scan:
    """One strand's view of the locus, with the coordinate map back.

    Both strands are scanned in *protospacer space* — the sequence as a
    protospacer on that strand reads it, 5' to 3'. For the minus strand that is
    the reverse complement of the window. Carrying one mapping function beats
    two sets of index arithmetic, which is where off-by-ones live.

    The patient view and the edited view are aligned wherever it matters: they
    share everything up to the edit, the nick is upstream of the edit by
    construction, so an insertion or a deletion shifts only the region past the
    edit — which is inside the template and nowhere else.
    """

    strand: str
    patient_view: str
    edited_view: str
    edit_index: int
    window_start: int
    window_end: int

    def to_genomic(self, index: int) -> int:
        if self.strand == PLUS:
            return self.window_start + index
        return self.window_end - index


def _views(
    strand: str,
    patient_sequence: str,
    edited_sequence: str,
    edit_index: int,
    patient_length: int,
    window_start: int,
    window_end: int,
) -> _Scan:
    if strand == PLUS:
        return _Scan(
            strand=strand,
            patient_view=patient_sequence,
            edited_view=edited_sequence,
            edit_index=edit_index,
            window_start=window_start,
            window_end=window_end,
        )
    return _Scan(
        strand=strand,
        patient_view=reverse_complement(patient_sequence),
        edited_view=reverse_complement(edited_sequence),
        # The allele's last base on the plus strand is its first on the minus.
        edit_index=len(patient_sequence) - (edit_index + patient_length),
        window_start=window_start,
        window_end=window_end,
    )


def design_pegrnas(
    request: EditRequest,
    sequences: SequenceProvider,
    rules: FlagRuleset,
    *,
    search_padding_nt: int = 150,
) -> PrimeOutcome:
    """Every pegRNA in the rule file's design space that installs this edit."""
    nuclease = _nuclease(rules)
    patient = request.patient_allele.upper()
    wild_type = request.wild_type_allele.upper()

    window_start = max(1, request.position - search_padding_nt)
    window_end = request.position + search_padding_nt
    reference = sequences.fetch(request.chromosome, window_start, window_end).upper()

    index = request.position - window_start
    notes = _reference_note(reference, index, patient, wild_type, request)

    # The patient's genome, and the genome as it should read. Protospacers are
    # found in the first; the template writes the second.
    #
    # Both splice out ``len(wild_type)`` bases, because that is what the
    # *reference* carries at this index — and both used to splice out
    # ``len(patient)`` instead. For a substitution the two are equal and nothing
    # showed; for an insertion or a deletion, which is the entire reason this
    # module exists rather than the base editor, both sequences came out wrong.
    # An insertion of four bases silently consumed four reference bases, so the
    # "patient" sequence was the reference unchanged and the "edited" sequence
    # was the reference with four bases deleted — a template that writes a
    # deletion nobody asked for. Reproducing a published pegRNA is what showed
    # it: the primer binding site matched the paper exactly, and the reverse
    # transcription template encoded a different product.
    patient_sequence = reference[:index] + patient + reference[index + len(wild_type) :]
    edited_sequence = reference[:index] + wild_type + reference[index + len(wild_type) :]

    candidates: list[PegRna] = []
    for strand in (PLUS, MINUS):
        scan = _views(
            strand, patient_sequence, edited_sequence, index, len(patient), window_start, window_end
        )
        candidates.extend(_scan_strand(request, rules, nuclease, scan, len(wild_type)))

    refusals: list[str] = []
    if not candidates:
        refusals.append(
            f"no protospacer within {search_padding_nt} nt puts a nick upstream of "
            f"g.{request.position} and within {rules.threshold('nick_to_edit_max_nt', 30)} nt of "
            "it. The nick must fall *before* the edit on the strand it cuts, which rules out "
            "every PAM downstream of the edit however good the PAM is — and it is why prime "
            "editing at some positions needs a PAM-relaxed nickase rather than a longer template"
        )

    # Ordered by stated criteria only, in the order the rule file states them:
    # nothing blocking first, then the primer binding site closest to the
    # melting temperature the file asks for, then the shortest reach. None of
    # these is a prediction of how well a pegRNA edits.
    target_tm = float(rules.threshold("pbs_target_tm_c", 30))
    ordered = sorted(
        candidates,
        key=lambda candidate: (
            sort_weight(candidate.severity),
            abs(melting_temperature(candidate.pbs) - target_tm),
            candidate.nick_to_edit_nt,
            len(candidate.rtt),
        ),
    )

    return PrimeOutcome(
        gene=request.gene,
        chromosome=request.chromosome,
        position=request.position,
        edit=f"{patient}→{wild_type} ({request.kind})",
        candidates=tuple(ordered),
        refusals=tuple(refusals),
        notes=tuple(notes),
        ruleset_pin=rules.pin,
        ranking=(
            "no efficiency model is attached. These are ordered by the severity of what is "
            "wrong with them and then by distance from the nick — a statement about design "
            "rules, not a prediction of how much editing any of them would produce. That "
            "prediction is what PRIDICT and its successors make from thousands of measured "
            "outcomes, and none of them is running here"
        ),
    )


def _reference_note(
    reference: str,
    index: int,
    patient: str,
    wild_type: str,
    request: EditRequest,
) -> list[str]:
    """Which of the two alleles the reference carries, when it is worth saying.

    Three cases, and only one of them is silent. Normally the reference carries
    the wild-type allele and there is nothing to report. Sometimes it carries
    the patient's — a pathogenic allele can be the reference allele, and the
    design is still correct, but a reader should know that "reference" and
    "wild type" have come apart here. And sometimes it carries neither, which is
    what a wrong coordinate or a wrong assembly looks like and is the only one
    of the three that usually means a mistake.
    """
    if reference[index : index + len(wild_type)] == wild_type:
        return []
    if reference[index : index + len(patient)] == patient:
        return [
            f"the reference carries the patient's allele at {request.chromosome}:"
            f"{request.position}, not the wild-type one. That happens — a pathogenic allele "
            "can be the reference allele — and the design below is unaffected, but nothing "
            "downstream should treat 'reference' and 'wild type' as the same word here"
        ]
    observed = reference[index : index + max(len(patient), len(wild_type))]
    return [
        f"the reference reads {observed} at {request.chromosome}:{request.position} and matches "
        f"neither allele ({patient} or {wild_type}). Every pegRNA below would still be "
        "synthesisable and none of them would edit what you think — this is what a wrong "
        "coordinate or a wrong assembly looks like"
    ]


def _scan_strand(
    request: EditRequest,
    rules: FlagRuleset,
    nuclease: _Nuclease,
    scan: _Scan,
    wild_type_length: int,
) -> list[PegRna]:
    """Every pegRNA on one strand, one template length per placement."""
    length, pam_length = nuclease.protospacer_length, len(nuclease.pam)
    max_distance = int(rules.threshold("nick_to_edit_max_nt", 30))
    homology_min = int(rules.threshold("rtt_homology_min_nt", 5))
    homology_max = int(rules.threshold("rtt_homology_max_nt", 20))
    rtt_max = int(rules.threshold("rtt_max_nt", 40))
    pbs_lengths = range(
        int(rules.threshold("pbs_min_nt", 8)), int(rules.threshold("pbs_max_nt", 17)) + 1
    )

    found: list[PegRna] = []
    for pam_index in range(length, len(scan.patient_view) - pam_length + 1):
        pam = scan.patient_view[pam_index : pam_index + pam_length]
        if not matches_pam(pam, nuclease.pam):
            continue

        start = pam_index - length
        protospacer = scan.patient_view[start:pam_index]
        if not is_resolved(protospacer):
            continue

        # The last base retained on the nicked strand. Everything after it is
        # re-synthesised from the template.
        nick_index = start + nuclease.nick_after_position - 1
        distance = scan.edit_index - nick_index
        if distance < 1 or distance > max_distance:
            continue

        # Every template that carries the edit, from the minimum homology arm up
        # to the maximum the rule file allows. Both this and the primer binding
        # site are scanned, because both are scanned at the bench — emitting one
        # template length is what kept the published HEXA pegRNA out of this
        # module's design space entirely.
        for homology_arm in range(homology_min, homology_max + 1):
            rtt_length = distance + wild_type_length + homology_arm - 1
            if rtt_length > rtt_max or nick_index + 1 + rtt_length > len(scan.edited_view):
                break
            template = scan.edited_view[nick_index + 1 : nick_index + 1 + rtt_length]
            if not is_resolved(template):
                break
            rtt = reverse_complement(template)

            for pbs_length in pbs_lengths:
                if nick_index - pbs_length + 1 < 0:
                    continue
                pbs = reverse_complement(
                    scan.patient_view[nick_index - pbs_length + 1 : nick_index + 1]
                )
                span = sorted((scan.to_genomic(start), scan.to_genomic(pam_index - 1)))
                found.append(
                    PegRna(
                        spacer=protospacer,
                        pam=pam,
                        strand=scan.strand,
                        chromosome=request.chromosome,
                        protospacer_span=(span[0], span[1]),
                        nick_position=scan.to_genomic(nick_index),
                        nick_to_edit_nt=distance,
                        pbs=pbs,
                        rtt=rtt,
                        homology_arm_nt=homology_arm,
                        pam_disrupted_by_edit=_pam_disrupted(rtt, pam, nuclease),
                    )
                )

    return [_reviewed(candidate, request, rules, nuclease, scan) for candidate in found]


def _pam_disrupted(rtt: str, pam: str, nuclease: _Nuclease) -> bool:
    """Does the edit destroy the PAM this pegRNA used?

    Read out of the template rather than off the genome, which is what makes it
    correct for insertions and deletions: the reverse complement of the template
    *is* the new sequence downstream of the nick, so the PAM's new reading is at
    a fixed offset in it however much the edit shifted the coordinates.

    The answer matters because a surviving PAM means the editor can re-engage
    the allele it has just corrected, and the second pass installs an indel
    rather than the edit.
    """
    new_sequence = reverse_complement(rtt)
    offset = nuclease.pam_offset_after_nick
    return new_sequence[offset : offset + len(pam)].upper() != pam.upper()


def _reviewed(
    candidate: PegRna,
    request: EditRequest,
    rules: FlagRuleset,
    nuclease: _Nuclease,
    scan: _Scan,
) -> PegRna:
    """One candidate, with its second nicks found and its rule file applied."""
    guides = _nicking_guides(candidate, rules, nuclease, scan)
    features = _features(candidate, request, rules, guides)
    return PegRna(
        spacer=candidate.spacer,
        pam=candidate.pam,
        strand=candidate.strand,
        chromosome=candidate.chromosome,
        protospacer_span=candidate.protospacer_span,
        nick_position=candidate.nick_position,
        nick_to_edit_nt=candidate.nick_to_edit_nt,
        pbs=candidate.pbs,
        rtt=candidate.rtt,
        homology_arm_nt=candidate.homology_arm_nt,
        pam_disrupted_by_edit=candidate.pam_disrupted_by_edit,
        nicking_guides=guides,
        flags=rules.raise_flags(features),
    )


def _features(
    candidate: PegRna,
    request: EditRequest,
    rules: FlagRuleset,
    guides: tuple[NickingGuide, ...],
) -> FlatFeatures:
    """One pegRNA, flattened into the names the rule file uses."""
    run = "T" * int(rules.threshold("polyt_run_nt", 4))
    extension = candidate.extension

    return FlatFeatures(
        values={
            "pegrna.spacer_length_nt": len(candidate.spacer),
            "pegrna.pbs_length_nt": len(candidate.pbs),
            "pegrna.pbs_tm_c": melting_temperature(candidate.pbs),
            "pegrna.pbs_gc_fraction": _gc(candidate.pbs),
            "pegrna.rtt_length_nt": len(candidate.rtt),
            "pegrna.homology_arm_nt": candidate.homology_arm_nt,
            "pegrna.nick_to_edit_nt": candidate.nick_to_edit_nt,
            "pegrna.extension_begins_with_c": extension.startswith("C"),
            "pegrna.has_polyt_run": run in candidate.spacer or run in extension,
            "pegrna.pam_disrupted_by_edit": candidate.pam_disrupted_by_edit,
            "pegrna.edit_kind": request.kind,
            "pegrna.edit_length_nt": max(
                len(request.patient_allele), len(request.wild_type_allele)
            ),
            "pegrna.pe3_options": len(guides),
            "pegrna.pe3b_available": any(guide.edit_dependent for guide in guides),
            "pegrna.strand": candidate.strand,
        }
    )


def _gc(sequence: str) -> float:
    if not sequence:
        return 0.0
    return round((sequence.count("G") + sequence.count("C")) / len(sequence), 3)


def _nicking_guides(
    candidate: PegRna,
    rules: FlagRuleset,
    nuclease: _Nuclease,
    scan: _Scan,
) -> tuple[NickingGuide, ...]:
    """Second nicks on the opposite strand, in the window where PE3 helps.

    The PE3b test is the part worth reading closely. A nicking guide is PE3b
    when its protospacer matches the *edited* sequence and not the patient's:
    it cannot cut until the edit is already installed, so the two nicks are
    never open at the same time. That simultaneity is the mechanism behind PE3's
    indels, and avoiding it costs nothing when the geometry allows it.
    """
    opposite = MINUS if candidate.strand == PLUS else PLUS
    minimum = int(rules.threshold("pe3_nick_min_nt", 40))
    maximum = int(rules.threshold("pe3_nick_max_nt", 90))
    pe3b_maximum = int(rules.threshold("pe3b_search_nt", 60))

    # The opposite strand's view of the same locus, built from the same two
    # sequences so the PE3b comparison is like for like.
    patient_view = reverse_complement(scan.patient_view)
    edited_view = reverse_complement(scan.edited_view)
    length_delta = len(edited_view) - len(patient_view)

    def to_genomic(index: int) -> int:
        mirrored = len(patient_view) - 1 - index
        return scan.to_genomic(mirrored)

    length, pam_length = nuclease.protospacer_length, len(nuclease.pam)
    guides: list[NickingGuide] = []

    for pam_index in range(length, len(patient_view) - pam_length + 1):
        start = pam_index - length
        nick_index = start + nuclease.nick_after_position - 1
        distance = abs(to_genomic(nick_index) - candidate.nick_position)
        if distance > max(maximum, pe3b_maximum) or distance == 0:
            continue

        patient_pam = patient_view[pam_index : pam_index + pam_length]
        patient_protospacer = patient_view[start:pam_index]
        # On this strand the edit sits upstream in index terms, so an indel
        # shifts the edited view by its length difference. Comparing without
        # that shift would call every guide past the edit "edit-dependent".
        edited_start, edited_pam_index = start + length_delta, pam_index + length_delta
        edited_pam = edited_view[edited_pam_index : edited_pam_index + pam_length]
        edited_protospacer = edited_view[edited_start:edited_pam_index]

        matches_patient = matches_pam(patient_pam, nuclease.pam) and is_resolved(
            patient_protospacer
        )
        matches_edited = matches_pam(edited_pam, nuclease.pam) and is_resolved(edited_protospacer)

        if matches_edited and edited_protospacer != patient_protospacer:
            if distance > pe3b_maximum:
                continue
            guides.append(
                NickingGuide(
                    protospacer=edited_protospacer,
                    pam=edited_pam,
                    strand=opposite,
                    nick_position=to_genomic(nick_index),
                    distance_nt=distance,
                    edit_dependent=True,
                )
            )
        elif matches_patient:
            if not minimum <= distance <= maximum:
                continue
            guides.append(
                NickingGuide(
                    protospacer=patient_protospacer,
                    pam=patient_pam,
                    strand=opposite,
                    nick_position=to_genomic(nick_index),
                    distance_nt=distance,
                )
            )

    # PE3b first, then closest.
    guides.sort(key=lambda guide: (not guide.edit_dependent, guide.distance_nt))
    return tuple(guides)
