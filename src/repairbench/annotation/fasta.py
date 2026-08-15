"""Random access to a reference FASTA, without loading a genome into memory.

Left-aligning an indel needs the reference bases immediately upstream of it, one
at a time, and a human genome is three gigabytes. The standard answer is the
``.fai`` index that ``samtools faidx`` writes: five numbers per sequence that
turn a coordinate into a file offset. Reading it is a dozen lines, so this
package reads it rather than taking a dependency for the privilege.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

from repairbench.annotation.gff import AnnotationError


class SequenceProvider(Protocol):
    """Whatever can hand back reference bases for a genomic interval."""

    def fetch(self, chromosome: str, start: int, end: int) -> str:
        """Return the reference sequence for a 1-based inclusive interval."""


@dataclass(frozen=True, slots=True)
class _IndexEntry:
    length: int
    offset: int
    bases_per_line: int
    bytes_per_line: int


class IndexedFasta:
    """A FASTA read through its ``.fai`` index.

    The index must exist. Building one means reading the whole file, which is
    exactly what the index exists to avoid, and doing it silently would turn a
    missing-file problem into a mysterious pause.
    """

    def __init__(self, path: str | Path) -> None:
        self._path = Path(path)
        index_path = Path(f"{self._path}.fai")
        if not index_path.exists():
            raise AnnotationError(
                f"{self._path} has no .fai index — run `samtools faidx` on it first"
            )
        self._index = _read_index(index_path)
        self._handle = self._path.open("rb")

    def close(self) -> None:
        self._handle.close()

    def __enter__(self) -> IndexedFasta:
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    @property
    def chromosomes(self) -> list[str]:
        return sorted(self._index)

    def fetch(self, chromosome: str, start: int, end: int) -> str:
        """Reference sequence for a 1-based inclusive interval, upper-cased.

        Soft-masked genomes store repeats in lower case, and a repeat region is
        precisely where indel left-alignment matters most — so case is
        normalised here rather than left as a trap for a string comparison
        three modules away.
        """
        entry = self._index.get(chromosome)
        if entry is None:
            raise AnnotationError(
                f"{chromosome!r} is not in {self._path.name}; it has {', '.join(self.chromosomes)}"
            )
        if start < 1 or end > entry.length or start > end:
            raise AnnotationError(
                f"{chromosome}:{start}-{end} is outside the sequence (length {entry.length})"
            )

        wanted = end - start + 1
        line, column = divmod(start - 1, entry.bases_per_line)
        self._handle.seek(entry.offset + line * entry.bytes_per_line + column)

        # Read generously: the interval spans line breaks, which are bytes in
        # the file and not bases in the sequence.
        newline_overhead = entry.bytes_per_line - entry.bases_per_line
        span = wanted + (wanted // entry.bases_per_line + 2) * newline_overhead
        raw = self._handle.read(span).decode()
        return "".join(character for character in raw if not character.isspace())[:wanted].upper()


def _read_index(path: Path) -> dict[str, _IndexEntry]:
    index: dict[str, _IndexEntry] = {}
    for number, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            raise AnnotationError(f"{path}:{number}: expected 5 columns in a .fai record")
        name, length, offset, bases, width = fields[:5]
        index[name] = _IndexEntry(int(length), int(offset), int(bases), int(width))
    if not index:
        raise AnnotationError(f"{path}: index is empty")
    return index


@dataclass(frozen=True, slots=True)
class InMemorySequences:
    """A sequence provider backed by a dictionary. For tests and small fixtures."""

    sequences: dict[str, str]

    def fetch(self, chromosome: str, start: int, end: int) -> str:
        sequence = self.sequences.get(chromosome)
        if sequence is None:
            raise AnnotationError(f"{chromosome!r} is not present")
        if start < 1 or end > len(sequence) or start > end:
            raise AnnotationError(
                f"{chromosome}:{start}-{end} is outside the sequence (length {len(sequence)})"
            )
        return sequence[start - 1 : end].upper()
