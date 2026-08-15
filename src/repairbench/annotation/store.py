"""Choosing which transcript to interpret a variant against.

A gene has many transcripts and they disagree about where the exons are, so the
choice is not a detail — it decides the NMD prediction and therefore the
therapeutic direction. The rule here is the one the field settled on: prefer
MANE Select, the transcript RefSeq and Ensembl agreed to call the reference for
clinical reporting. Where there is none, take the longest coding sequence and
say so, because a silent fallback would hide the fact that a choice was made.
"""

from __future__ import annotations

from dataclasses import dataclass

from repairbench.annotation.gff import Annotation, AnnotationError, TranscriptRecord
from repairbench.transcript import Transcript


@dataclass(frozen=True, slots=True)
class Placement:
    """What the annotation says is at a coordinate — possibly nothing.

    All four fields being empty is a real answer and a common one: most of a
    genome is not in any transcript. It is kept distinct from "we did not look",
    which is the absence of a ``Placement`` altogether.
    """

    gene: str | None = None
    accession: str | None = None
    in_coding_sequence: bool = False
    in_transcript_span: bool = False

    @property
    def describe(self) -> str:
        if self.gene is None:
            return "not in any transcript in this annotation"
        where = "coding sequence" if self.in_coding_sequence else "intron or untranslated region"
        return f"{self.gene} ({self.accession}), {where}"


@dataclass(frozen=True, slots=True)
class ResolvedVariant:
    """A variant placed on a transcript: which one, where, and how it was chosen."""

    record: TranscriptRecord
    cds_position: int
    #: Why this transcript and not another — carried so a report can say.
    selection_reason: str

    @property
    def transcript(self) -> Transcript:
        return self.record.to_transcript()


class TranscriptStore:
    """Lookup over a parsed annotation."""

    def __init__(self, annotation: Annotation) -> None:
        self._annotation = annotation
        self._by_gene: dict[str, list[TranscriptRecord]] = {}
        for record in annotation.transcripts.values():
            self._by_gene.setdefault(record.gene, []).append(record)

    @property
    def pin(self) -> str:
        """The annotation release this store speaks for."""
        return self._annotation.pin

    @property
    def genes(self) -> list[str]:
        return sorted(self._by_gene)

    def by_accession(self, accession: str) -> TranscriptRecord:
        record = self._annotation.transcripts.get(accession)
        if record is None:
            raise AnnotationError(
                f"{accession} is not in {self._annotation.pin}. Note that the version matters: "
                "exon boundaries move between transcript versions."
            )
        return record

    def preferred_for(self, gene: str) -> tuple[TranscriptRecord, str]:
        """The transcript to interpret this gene against, and why."""
        candidates = self._by_gene.get(gene)
        if not candidates:
            raise AnnotationError(f"{gene} has no coding transcript in {self._annotation.pin}")

        mane = [record for record in candidates if record.mane_select]
        if len(mane) == 1:
            return mane[0], "MANE Select"
        if len(mane) > 1:
            # Two MANE Select transcripts for one gene should be impossible. If
            # the file says otherwise, the file is wrong and guessing would bury
            # that under a plausible answer.
            raise AnnotationError(
                f"{gene} has {len(mane)} transcripts tagged MANE Select in "
                f"{self._annotation.pin}: {', '.join(sorted(r.accession for r in mane))}"
            )

        longest = max(candidates, key=lambda record: (record.coding_length, record.accession))
        return longest, (
            f"longest coding sequence ({longest.coding_length} nt of "
            f"{len(candidates)} transcripts) — no MANE Select transcript is annotated"
        )

    def locate(self, chromosome: str, position: int) -> Placement:
        """What sits at a genomic coordinate, for a position nobody chose.

        A different question from ``resolve``, and the difference is why both
        exist. ``resolve`` places a variant somebody is asking about, and refuses
        when it is not in coding sequence. ``locate`` is asked about an
        off-target hit — a coordinate a search picked rather than a clinician —
        where "intronic" and "nowhere near a gene" are answers rather than
        errors, and are the answers that decide how much the hit matters.
        """
        best: Placement | None = None
        for record in self._annotation.transcripts.values():
            if record.seqid != chromosome or not record.cds_blocks:
                continue
            starts = [start for start, _ in record.cds_blocks]
            ends = [end for _, end in record.cds_blocks]
            if not min(starts) <= position <= max(ends):
                continue
            placement = Placement(
                gene=record.gene,
                accession=record.accession,
                in_coding_sequence=record.cds_offset(position) is not None,
                in_transcript_span=True,
            )
            # A coding placement outranks an intronic one at the same
            # coordinate. Transcripts overlap, and reporting the intronic
            # reading of a position that is exonic in another transcript would
            # understate the finding.
            if placement.in_coding_sequence:
                return placement
            best = best or placement
        return best or Placement()

    def resolve(
        self,
        gene: str,
        chromosome: str,
        position: int,
        *,
        accession: str | None = None,
    ) -> ResolvedVariant:
        """Place a genomic coordinate on a transcript as a CDS offset.

        Refuses rather than guesses in the two cases that matter: a coordinate on
        the wrong chromosome, and a coordinate that is not in coding sequence at
        all. The second is not an error in general — introns exist — but it is an
        error *here*, because every rule downstream reads a CDS offset, and
        inventing one would put the variant in an exon it is not in.
        """
        if accession is not None:
            record, reason = self.by_accession(accession), "requested explicitly"
        else:
            record, reason = self.preferred_for(gene)

        if record.seqid != chromosome:
            raise AnnotationError(
                f"{record.accession} is on {record.seqid}, but the variant is on {chromosome}"
            )

        offset = record.cds_offset(position)
        if offset is None:
            raise AnnotationError(
                f"{chromosome}:{position} is not inside the coding sequence of "
                f"{record.accession} — it is intronic or untranslated, and this module "
                "reasons in CDS coordinates"
            )
        return ResolvedVariant(record=record, cds_position=offset, selection_reason=reason)
