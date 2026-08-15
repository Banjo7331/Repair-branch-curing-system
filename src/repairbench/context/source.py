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
import hashlib
from collections.abc import Iterator
from dataclasses import dataclass, field
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
        return repr(self.value)


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


def read_tsv(source: Source, *, required: set[str]) -> Iterator[dict[str, str]]:
    """Stream a tab-separated file, keyed by header name.

    Header-keyed rather than positional on purpose: column *order* changes
    between releases of both files this package reads, and a positional parser
    would keep working while silently reading the wrong column.
    """
    with source.path.open(newline="") as handle:
        lines = (line for line in handle if not line.startswith("##"))
        reader = csv.DictReader(lines, delimiter="\t")
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
