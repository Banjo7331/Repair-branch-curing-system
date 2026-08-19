"""Parsing GFF3 into transcripts the rest of the package can reason about.

The parser is deliberately narrow. It reads ``mRNA`` and ``CDS`` records and
ignores everything else, because everything else is either derivable or
irrelevant to the one question this package asks a transcript: where does a
coding position sit relative to the exon-exon junctions.

Three things it gets right that a naive parser gets wrong, and each of them
silently corrupts an NMD prediction rather than raising:

* **Strand.** On the minus strand the first coding exon is the one with the
  *highest* genomic coordinate. Sorting CDS blocks by start position and calling
  it transcript order inverts the exon numbering for half the genome.
* **Coding, not exonic.** UTRs are not in the CDS records, which is what makes
  the offsets here CDS offsets. Mixing exon and CDS coordinates shifts every
  position by the length of the 5′ UTR.
* **Versioned accessions.** Exon boundaries move between transcript versions, so
  the version is part of the identity and is kept.
"""

from __future__ import annotations

import gzip
import hashlib
import re
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

from repairbench.annotation.naming import NO_ALIASES, Aliases
from repairbench.model import RepairbenchError
from repairbench.transcript import Transcript


class AnnotationError(RepairbenchError):
    """The annotation file is malformed, or describes something unusable."""


_ATTRIBUTE = re.compile(r"([^=;]+)=([^;]*)")


def _parse_attributes(field_text: str) -> dict[str, str]:
    return {key.strip(): value.strip() for key, value in _ATTRIBUTE.findall(field_text)}


@dataclass(slots=True)
class TranscriptRecord:
    """One transcript as the annotation describes it.

    ``cds_blocks`` are 1-based inclusive genomic intervals **in transcript
    order** — already reversed for minus-strand genes, so index 0 is always the
    first coding exon translated.
    """

    accession: str
    gene: str
    seqid: str
    strand: str
    cds_blocks: list[tuple[int, int]] = field(default_factory=list)
    mane_select: bool = False
    #: Empty when the transcript is usable; otherwise why it is not. Set by
    #: ``finalise`` and read by the parser, which drops it and counts it.
    unusable: str = ""

    @property
    def is_curated(self) -> bool:
        """Did a human curate this transcript, or did a pipeline model it?

        RefSeq says so in the accession: ``NM_`` and ``NR_`` are curated,
        ``XM_`` and ``XR_`` are computed predictions. The distinction is invisible
        in a synthetic fixture — everything there is invented equally — and it
        decides the answer for any gene in a duplicated region, where the curated
        transcript often sits on a scaffold and the chromosome carries only
        models.
        """
        return self.accession.startswith(("NM_", "NR_"))

    @property
    def coding_exon_lengths(self) -> tuple[int, ...]:
        return tuple(end - start + 1 for start, end in self.cds_blocks)

    @property
    def coding_length(self) -> int:
        return sum(self.coding_exon_lengths)

    def to_transcript(self) -> Transcript:
        """The domain object the NMD arithmetic runs on."""
        return Transcript(
            accession=self.accession,
            gene=self.gene,
            coding_exon_lengths=self.coding_exon_lengths,
            mane_select=self.mane_select,
        )

    def cds_offset(self, genomic_position: int) -> int | None:
        """Map a genomic coordinate onto a 1-based CDS offset.

        Returns ``None`` for a position outside the coding sequence — an
        intronic or untranslated coordinate is not an error, it simply has no
        CDS offset, and the caller decides what that means.
        """
        offset = 0
        for start, end in self.cds_blocks:
            length = end - start + 1
            if start <= genomic_position <= end:
                within = (
                    genomic_position - start + 1
                    if self.strand == "+"
                    else end - genomic_position + 1
                )
                return offset + within
            offset += length
        return None

    def genomic_position(self, cds_offset: int) -> int:
        """The inverse. Raises rather than guessing on an out-of-range offset."""
        if not 1 <= cds_offset <= self.coding_length:
            raise AnnotationError(
                f"c.{cds_offset} is outside {self.accession} (coding length {self.coding_length})"
            )
        remaining = cds_offset
        for start, end in self.cds_blocks:
            length = end - start + 1
            if remaining <= length:
                return start + remaining - 1 if self.strand == "+" else end - remaining + 1
            remaining -= length
        raise AnnotationError(f"c.{cds_offset} not located in {self.accession}")

    def finalise(self) -> None:
        """Put the blocks in transcript order, and say if they cannot be used.

        Returns nothing and raises nothing: the caller decides what to do with
        an unusable transcript, because in a real annotation the answer is "drop
        this one" rather than "reject the file".
        """
        self.cds_blocks.sort()
        self.unusable = self._overlap_reason()
        if self.strand == "-":
            self.cds_blocks.reverse()

    def _overlap_reason(self) -> str:
        """Why this transcript's blocks cannot carry CDS arithmetic, if they cannot.

        Overlapping CDS blocks are not a corrupt file. *PEG10* is the case that
        taught this package so: a retrotransposon-derived gene translated
        through a programmed ribosomal frameshift, where the ribosome slips back
        one base and reads on in another frame, so two CDS blocks legitimately
        share a coordinate. RefSeq annotates it correctly and this package
        cannot use it — a position inside the overlap has two CDS offsets, and
        every calculation downstream assumes it has one.

        So the transcript is dropped and the reason is kept. Refusing the whole
        file, which is what this used to do, made one gene's real biology deny
        access to the other hundred and eighty thousand transcripts in it.
        """
        for (_, previous_end), (next_start, _) in zip(
            self.cds_blocks, self.cds_blocks[1:], strict=False
        ):
            if next_start <= previous_end:
                return (
                    f"CDS blocks overlap at {next_start}, which is how a programmed "
                    "ribosomal frameshift is annotated — a coordinate inside the overlap has "
                    "two CDS offsets, and this package assumes one"
                )
        return ""


@dataclass(frozen=True, slots=True)
class Annotation:
    """A parsed annotation file, pinned by the digest of its bytes.

    The pin is not decoration. A transcript structure is an input to the NMD
    calculation exactly as a ClinVar release is an input to a classification,
    and a prediction that cannot name the annotation it was made against cannot
    be compared with one made a year later.
    """

    source: str
    digest: str
    transcripts: dict[str, TranscriptRecord]
    #: Transcripts the file describes and this package cannot use, with the
    #: reason for each. Never silently empty-by-omission: a real annotation has
    #: a handful, and a reader deciding whether a missing gene is missing or
    #: rejected needs to be able to tell.
    rejected: dict[str, str] = field(default_factory=dict)
    #: Every name each chromosome goes by, read from the file's own region
    #: records. Empty for an annotation that declares none, which is every
    #: hand-written fixture and is fine — the ``chr`` prefix still resolves.
    aliases: Aliases = NO_ALIASES

    @property
    def short_digest(self) -> str:
        return self.digest[:12]

    @property
    def pin(self) -> str:
        return f"{Path(self.source).name}@{self.short_digest}"

    def __len__(self) -> int:
        return len(self.transcripts)


def _open(path: Path) -> Iterator[str]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as handle:
            yield from handle
    else:
        with path.open() as handle:
            yield from handle


def parse_gff3(path: str | Path, *, genes: set[str] | None = None) -> Annotation:
    """Read transcripts from a GFF3 file.

    ``genes`` restricts the parse to a set of symbols, which is what makes this
    usable against a full RefSeq annotation: the file is gigabytes, the question
    concerns a handful of genes, and there is no reason to hold the rest.
    """
    path = Path(path)
    digest = _digest(path)

    records: dict[str, TranscriptRecord] = {}
    by_feature_id: dict[str, TranscriptRecord] = {}
    aliases: list[tuple[str, str]] = []

    for line in _open(path):
        if line.startswith("#") or not line.strip():
            continue
        columns = line.rstrip("\n").split("\t")
        if len(columns) != 9:
            raise AnnotationError(f"{path}: expected 9 GFF3 columns, got {len(columns)}")

        seqid, _, feature, start, end, _, strand, _, attribute_text = columns
        if feature not in {"mRNA", "CDS", "region"}:
            continue

        if feature == "region":
            alias = _chromosome_alias(seqid, attribute_text)
            if alias is not None:
                aliases.append(alias)
            continue

        attributes = _parse_attributes(attribute_text)

        if feature == "mRNA":
            record = _transcript_from(attributes, seqid, strand, genes)
            if record is None:
                continue
            records[record.accession] = record
            if "ID" in attributes:
                by_feature_id[attributes["ID"]] = record
            continue

        # CDS
        parent = by_feature_id.get(attributes.get("Parent", ""))
        if parent is None:
            continue
        parent.cds_blocks.append((int(start), int(end)))

    usable, rejected = _coding_only(records)

    # An annotation whose transcripts were all *rejected* is a different thing
    # from one that parsed nothing at all, and it is returned rather than
    # raised: the reasons are the useful part, and a caller asking for a gene
    # that lost every transcript then gets a refusal naming the gene rather
    # than one blaming the file.
    if not usable and not rejected:
        raise AnnotationError(
            f"{path}: no coding transcripts parsed"
            + (f" for genes {sorted(genes)}" if genes else "")
        )

    return Annotation(
        source=str(path),
        digest=digest,
        transcripts=usable,
        rejected=rejected,
        aliases=Aliases.of(aliases),
    )


def _coding_only(
    records: dict[str, TranscriptRecord],
) -> tuple[dict[str, TranscriptRecord], dict[str, str]]:
    """Drop what has no CDS, and put what remains in transcript order.

    An mRNA with no CDS records is a non-coding transcript or a truncated file.
    Either way it cannot answer the question this package asks, and keeping it
    would produce a transcript of coding length zero — which reads downstream as
    a real structure rather than as an absent one.
    """
    usable: dict[str, TranscriptRecord] = {}
    rejected: dict[str, str] = {}
    for accession, record in records.items():
        if not record.cds_blocks:
            continue
        record.finalise()
        if record.unusable:
            rejected[accession] = record.unusable
            continue
        usable[accession] = record
    return usable, rejected


def _transcript_from(
    attributes: dict[str, str],
    seqid: str,
    strand: str,
    genes: set[str] | None,
) -> TranscriptRecord | None:
    """One mRNA record, or ``None`` when it is not one this parse wants.

    Real annotation carries the gene symbol under ``gene`` and the accession
    under ``transcript_id``; hand-written files often use ``gene_name`` and
    ``Name``. Both are read, in that order, because a parser that only knew one
    convention would return an empty annotation and no explanation.
    """
    gene = attributes.get("gene") or attributes.get("gene_name", "")
    if genes is not None and gene not in genes:
        return None
    accession = (
        attributes.get("transcript_id") or attributes.get("Name") or attributes.get("ID", "")
    )
    if not accession:
        return None
    return TranscriptRecord(
        accession=accession,
        gene=gene,
        seqid=seqid,
        strand=strand,
        mane_select="MANE Select" in attributes.get("tag", ""),
    )


def _chromosome_alias(seqid: str, attribute_text: str) -> tuple[str, str] | None:
    """The ordinary name of a chromosome, out of the record that declares it.

    NCBI opens each chromosome with a ``region`` record carrying
    ``chromosome=17``, so the mapping from accession to ordinary name is in the
    file, written by the people who assigned both. Reading it beats a hard-coded
    table, which goes stale with every assembly patch.

    Scaffolds and patches carry the attribute too, so the ``genome=chromosome``
    check is what keeps an unplaced contig from claiming a chromosome's name.
    """
    attributes = _parse_attributes(attribute_text)
    name = attributes.get("chromosome")
    if not name or attributes.get("genome") != "chromosome":
        return None
    return (seqid, name)


def _digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            hasher.update(chunk)
    return hasher.hexdigest()
