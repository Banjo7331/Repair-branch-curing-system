"""A source file, its identity, and where each fact it supplied came from.

Two ideas carry this module.

A **source** is a file plus the digest of its bytes. That digest is what turns a
pin from a declaration into something checkable: before this package, a world
could claim ``gene_curation@2026-01`` with nothing behind the label.

**Provenance** is per fact, not per gene. A gene's context is assembled from
several files, and a report that says "dominant negative, because null alleles
are milder" is only reviewable if a reader can find out that the milder-nulls
claim came from our own curation with a PMID attached, while the dosage score
came from ClinGen. Recording the source at the level of the gene would lose
exactly the distinction that matters.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import itertools
from collections.abc import Iterator
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

from repairbench.model import RepairbenchError


class ContextError(RepairbenchError):
    """A source file is missing, malformed, or does not say what it claims to."""


@dataclass(frozen=True, slots=True)
class Source:
    """One ingested file, pinned by the digest of its bytes."""

    name: str
    path: Path
    digest: str
    #: The release label the file carries or was given, e.g. "2026-01".
    version: str

    @property
    def short_digest(self) -> str:
        return self.digest[:12]

    @property
    def pin(self) -> str:
        """The citation form. This is what makes a world's pin *earned*."""
        return f"{self.name}@{self.version}/{self.short_digest}"

    @classmethod
    def of(cls, name: str, path: str | Path, version: str) -> Source:
        path = Path(path)
        if not path.exists():
            raise ContextError(f"{name}: no such file: {path}")
        return cls(name=name, path=path, digest=_digest(path), version=version)


def _digest(path: Path) -> str:
    """The content digest of a file, or of a directory taken as a whole.

    Some releases are a directory rather than a file — the rule files and our own
    curation version together, because the same people edit them in the same
    review. Digesting the directory means its pin covers the combination, so a
    report cannot cite a mixture that nobody reviewed as one. Names are included
    in the hash, so adding a file changes the pin even if no existing byte moved.
    """
    hasher = hashlib.sha256()
    if path.is_dir():
        for child in sorted(p for p in path.rglob("*") if p.is_file()):
            hasher.update(str(child.relative_to(path)).encode())
            hasher.update(b"\0")
            hasher.update(_file_digest(child).encode())
        return hasher.hexdigest()
    return _file_digest(path)


def _file_digest(path: Path) -> str:
    hasher = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


@dataclass(frozen=True, slots=True)
class Fact:
    """One value, and the source that supplied it."""

    field: str
    value: Any
    source: Source
    #: A citation, where the source carries one. Required for local curation and
    #: absent for downloaded tables, which cite themselves.
    citation: str = ""

    def describe(self) -> str:
        line = f"{self.field} = {self.rendered}  [{self.source.pin}]"
        return f"{line}  {self.citation}" if self.citation else line

    @property
    def rendered(self) -> str:
        """The value as a report shows it.

        A per-tissue expression profile is a fact with fifty numbers in it.
        Printing all of them would bury the four single-valued facts beside it,
        so it is summarised — and the summary says how many tissues, because
        "measured in 54 tissues" and "measured in one" are different claims
        about how much the value is worth.
        """
        if isinstance(self.value, dict):
            return f"{len(self.value)} tissues measured"
        # Enums before strings, because every enum here is a StrEnum and would
        # otherwise be printed by its repr — ``<DosageScore.SUFFICIENT_EVIDENCE:
        # 'sufficient_evidence'>``, which is Python talking about itself in a
        # line meant to tell a reviewer what ClinGen said.
        if isinstance(self.value, Enum):
            return str(self.value.value)
        # Strings keep their quotes, so a value that was supplied and is empty
        # reads as blank rather than as a line nobody wrote.
        if isinstance(self.value, str):
            return repr(self.value)
        return str(self.value)


@dataclass(slots=True)
class Provenance:
    """Every fact that went into one gene's context."""

    gene: str
    facts: dict[str, Fact] = field(default_factory=dict)

    def record(self, fact: Fact) -> None:
        existing = self.facts.get(fact.field)
        if existing is not None and existing.value != fact.value:
            raise ContextError(
                f"{self.gene}: {fact.field} was given as {existing.value!r} by "
                f"{existing.source.pin} and {fact.value!r} by {fact.source.pin}. "
                "Two sources disagree, and picking one silently would hide it."
            )
        self.facts[fact.field] = fact

    def source_for(self, field_name: str) -> Source | None:
        fact = self.facts.get(field_name)
        return fact.source if fact else None

    @property
    def sources(self) -> tuple[Source, ...]:
        seen: dict[str, Source] = {}
        for fact in self.facts.values():
            seen.setdefault(fact.source.pin, fact.source)
        return tuple(seen.values())

    def describe(self) -> str:
        lines = [f"{self.gene}"]
        ordered = sorted(self.facts.values(), key=lambda fact: fact.field)
        lines.extend(f"  {fact.describe()}" for fact in ordered)
        return "\n".join(lines)


def _open_text(path: Path) -> Iterator[str]:
    """Open a source file, transparently if it is gzipped.

    Every one of these releases ships compressed and several ship *only*
    compressed. Requiring the caller to gunzip first meant either a decompressed
    copy on disk next to the original, or a helpful error at the point of use —
    and the second is what happened: the parser read the gzip magic bytes as a
    header and complained about missing columns.
    """
    if path.suffix in _COMPRESSED:
        with gzip.open(path, "rt", newline="") as handle:
            yield from handle
    else:
        with path.open(newline="") as handle:
            yield from handle


#: Suffixes that mean gzip. ``.bgz`` is bgzip — block-gzip, which every gzip
#: reader accepts and which is what gnomAD and most genomics releases ship.
_COMPRESSED = {".gz", ".bgz", ".bgzf"}


#: The first line of a GCT file. The format declares itself, so skipping its
#: preamble is reading the format rather than guessing at it.
_GCT_MAGIC = "#1."


def _without_preamble(handle: Iterator[str]) -> Iterator[str]:
    """Drop the comment lines above the header, and keep the header itself.

    Real files put several. ClinGen's gene curation list opens with five lines
    of provenance and then a header that is *also* commented — ``#Gene
    Symbol``, tab-separated — and GTEx's GCT opens with two lines of counts
    before a header that is not commented at all.

    The rule that reads both: while lines begin with ``#``, keep only the last
    one that contains a tab, because a header has columns and a sentence does
    not. If none does, the first uncommented line is the header. Skipping every
    ``#`` line unconditionally — which this used to do — silently made
    ClinGen's first line of provenance into the column names, and the failure
    surfaced as "this file has no Gene Symbol column" about a file whose second
    column is Gene Symbol.
    """
    first = next(handle, None)
    if first is None:
        return
    if first.startswith(_GCT_MAGIC):
        # A GCT: version line, then a line of dimensions, then the header. The
        # dimensions line is tab-separated and uncommented, so anything looking
        # for "the first line with tabs" would take it for the columns and read
        # 56200 as a gene name.
        next(handle, None)
        yield from handle
        return
    handle = itertools.chain([first], handle)

    header: str | None = None
    for line in handle:
        if line.startswith("#"):
            if "\t" in line:
                header = line
            continue
        if header is not None:
            yield header
            header = None
        yield line
        break
    else:
        if header is not None:
            yield header
        return
    yield from handle


def read_tsv(source: Source, *, required: set[str]) -> Iterator[dict[str, str]]:
    """Stream a tab-separated file, keyed by header name.

    Header-keyed rather than positional on purpose: column *order* changes
    between releases of both files this package reads, and a positional parser
    would keep working while silently reading the wrong column.
    """
    reader = csv.DictReader(_without_preamble(_open_text(source.path)), delimiter="\t")
    fieldnames = [name.lstrip("#").strip() for name in (reader.fieldnames or [])]
    missing = sorted(required - set(fieldnames))
    if missing:
        raise ContextError(
            f"{source.name} ({source.path.name}) has no column(s): {', '.join(missing)}. "
            f"It has: {', '.join(fieldnames)}"
        )
    for row in reader:
        yield {
            key.lstrip("#").strip(): (value or "").strip()
            for key, value in row.items()
            if key
        }
