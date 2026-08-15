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
        """Put the blocks in transcript order and check they do not overlap."""
        self.cds_blocks.sort()
        for (_, previous_end), (next_start, _) in zip(
            self.cds_blocks, self.cds_blocks[1:], strict=False
        ):
            if next_start <= previous_end:
                raise AnnotationError(
                    f"{self.accession} has overlapping CDS blocks near {previous_end}"
                )
        if self.strand == "-":
            self.cds_blocks.reverse()


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

    for line in _open(path):
        if line.startswith("#") or not line.strip():
            continue
        columns = line.rstrip("\n").split("\t")
        if len(columns) != 9:
            raise AnnotationError(f"{path}: expected 9 GFF3 columns, got {len(columns)}")

        seqid, _, feature, start, end, _, strand, _, attribute_text = columns
        if feature not in {"mRNA", "CDS"}:
            continue

        attributes = _parse_attributes(attribute_text)

        if feature == "mRNA":
            gene = attributes.get("gene") or attributes.get("gene_name", "")
            if genes is not None and gene not in genes:
                continue
            accession = (
                attributes.get("transcript_id")
                or attributes.get("Name")
                or attributes.get("ID", "")
            )
            if not accession:
                continue
            record = TranscriptRecord(
                accession=accession,
                gene=gene,
                seqid=seqid,
                strand=strand,
                mane_select="MANE Select" in attributes.get("tag", ""),
            )
            records[accession] = record
            if "ID" in attributes:
                by_feature_id[attributes["ID"]] = record
            continue

        # CDS
        parent = by_feature_id.get(attributes.get("Parent", ""))
        if parent is None:
            continue
        parent.cds_blocks.append((int(start), int(end)))

    usable: dict[str, TranscriptRecord] = {}
    for accession, record in records.items():
        if not record.cds_blocks:
            # An mRNA with no CDS records is a non-coding transcript or a
            # truncated file. Either way it cannot answer the question this
            # package asks, and silently keeping it would produce a transcript
            # of coding length zero.
            continue
        record.finalise()
        usable[accession] = record

    if not usable:
        raise AnnotationError(
            f"{path}: no coding transcripts parsed"
            + (f" for genes {sorted(genes)}" if genes else "")
        )

    return Annotation(source=str(path), digest=digest, transcripts=usable)


def _digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            hasher.update(chunk)
    return hasher.hexdigest()
