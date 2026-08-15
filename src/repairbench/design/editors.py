"""The editor catalogue, read from a file and pinned like every other rule file.

A base editor is three facts — which conversion its deaminase makes, which PAM
its nuclease needs, and where in the protospacer it acts — and all three are
published measurements that get revised. ABE8e's window is not a property of
arithmetic; it is what a group measured, and a later group will measure it
differently.

So the catalogue is data, digested, and cited in every candidate this module
produces. A design made under ``editors-v1@abc123`` and a design made after
somebody widened a window are not the same design, and a package that could not
tell them apart would eventually report a file edit as a better candidate.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any

import yaml

from repairbench.model import RepairbenchError


class DesignError(RepairbenchError):
    """A design cannot be produced, or the catalogue describing one is malformed."""


class Conversion(StrEnum):
    """What a deaminase does to a base on the protospacer strand.

    Two, and only two, are on offer. Every other single-base change — every
    transversion — is outside what base editing does at all, and the designer
    says so rather than returning an empty list that reads like "none found".
    """

    A_TO_G = "A>G"
    C_TO_T = "C>T"

    @property
    def source_base(self) -> str:
        return self.value[0]

    @property
    def product_base(self) -> str:
        return self.value[-1]

    @classmethod
    def between(cls, source: str, product: str) -> Conversion | None:
        """The conversion that turns one base into another, if a deaminase does it."""
        for conversion in cls:
            if (conversion.source_base, conversion.product_base) == (source, product):
                return conversion
        return None


@dataclass(frozen=True, slots=True)
class Editor:
    """One editor: a deaminase, a nuclease, and the window between them."""

    id: str
    conversion: Conversion
    nuclease: str
    #: IUPAC PAM, 3' of the protospacer. Every nuclease in this catalogue places
    #: it there; a 5' PAM nuclease would need the geometry generalising, and
    #: pretending otherwise would silently place protospacers backwards.
    pam: str
    protospacer_length: int
    window_start: int
    window_end: int
    because: str = ""
    citation: str = ""

    def __post_init__(self) -> None:
        if not 1 <= self.window_start <= self.window_end <= self.protospacer_length:
            raise DesignError(
                f"{self.id}: editing window {self.window_start}-{self.window_end} does not fit "
                f"in a protospacer of {self.protospacer_length} nt"
            )

    @property
    def window(self) -> range:
        """Protospacer positions the deaminase acts on, 1-based inclusive."""
        return range(self.window_start, self.window_end + 1)

    @property
    def label(self) -> str:
        return f"{self.id} ({self.conversion}, PAM {self.pam})"

    def edits(self, base: str) -> bool:
        return base.upper() == self.conversion.source_base


@dataclass(frozen=True, slots=True)
class CatalogueThresholds:
    search_padding_nt: int = 40
    crowded_window_bystanders: int = 2


@dataclass(frozen=True, slots=True)
class EditorCatalogue:
    """Every editor this package will design against, and the digest of the file."""

    version: str
    description: str
    thresholds: CatalogueThresholds
    editors: tuple[Editor, ...]
    digest: str

    @property
    def short_digest(self) -> str:
        return self.digest[:12]

    @property
    def pin(self) -> str:
        return f"{self.version}@{self.short_digest}"

    def making(self, conversion: Conversion) -> tuple[Editor, ...]:
        return tuple(editor for editor in self.editors if editor.conversion is conversion)

    def __iter__(self) -> Any:
        return iter(self.editors)

    def __len__(self) -> int:
        return len(self.editors)


def load_editors(path: str | Path) -> EditorCatalogue:
    """Read and validate the editor catalogue."""
    raw = Path(path).read_bytes()
    document = yaml.safe_load(raw)

    if not isinstance(document, dict):
        raise DesignError(f"{path}: editor catalogue must be a mapping at the top level")
    for required in ("version", "editors"):
        if required not in document:
            raise DesignError(f"{path}: editor catalogue has no {required!r}")

    thresholds = CatalogueThresholds(**(document.get("thresholds") or {}))

    editors: list[Editor] = []
    seen: set[str] = set()
    for index, entry in enumerate(document["editors"], start=1):
        editor = _parse_editor(entry, index, path)
        if editor.id in seen:
            raise DesignError(f"{path}: duplicate editor id {editor.id!r}")
        seen.add(editor.id)
        editors.append(editor)

    if not editors:
        raise DesignError(f"{path}: catalogue declares no editors")

    return EditorCatalogue(
        version=str(document["version"]),
        description=str(document.get("description", "")),
        thresholds=thresholds,
        editors=tuple(editors),
        digest=hashlib.sha256(raw).hexdigest(),
    )


def _parse_editor(entry: Any, index: int, path: str | Path) -> Editor:
    if not isinstance(entry, dict):
        raise DesignError(f"{path}: editor {index} is not a mapping")
    for required in ("id", "conversion", "pam", "protospacer_length", "window"):
        if required not in entry:
            raise DesignError(f"{path}: editor {index} has no {required!r}")

    try:
        conversion = Conversion(entry["conversion"])
    except ValueError as error:
        raise DesignError(
            f"{path}: editor {entry['id']!r} claims conversion {entry['conversion']!r}; "
            f"base editing makes {', '.join(c.value for c in Conversion)} and nothing else"
        ) from error

    window = entry["window"]
    if not isinstance(window, dict) or "start" not in window or "end" not in window:
        raise DesignError(f"{path}: editor {entry['id']!r} needs a window with start and end")

    return Editor(
        id=str(entry["id"]),
        conversion=conversion,
        nuclease=str(entry.get("nuclease", "unspecified")),
        pam=str(entry["pam"]).upper(),
        protospacer_length=int(entry["protospacer_length"]),
        window_start=int(window["start"]),
        window_end=int(window["end"]),
        because=" ".join(str(entry.get("because", "")).split()),
        citation=str(entry.get("citation", "")),
    )
