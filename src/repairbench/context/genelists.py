"""Which genes it is worse to hit by accident.

An off-target hit is a change to a gene somebody did not intend to change, and
how much that matters depends almost entirely on which gene. Two published kinds
of list answer most of it: genes a cell cannot lose (DepMap's common essentials)
and genes whose disruption drives cancer (COSMIC's Cancer Gene Census, split into
oncogenes and tumour suppressors — a distinction that matters, because a
loss-of-function hit in a tumour suppressor is the dangerous direction).

Neither is redistributable, so this package reads them rather than shipping
them: a two-column file of symbol and membership, pinned by content digest like
every other source here. A gene absent from the file is absent from the list,
not absent from the genome, and the risk rules are written to know the
difference.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path

from repairbench.context.source import ContextError, Source, read_tsv

SYMBOL_COLUMN = "symbol"
LIST_COLUMN = "list"
CITATION_COLUMN = "citation"


class GeneList(StrEnum):
    """The memberships the risk rules can ask about."""

    ESSENTIAL = "essential"
    """Loss is lethal to the cell in most backgrounds — DepMap common essentials."""

    ONCOGENE = "oncogene"
    """Activated by mutation. A bystander edit here is the wrong kind of luck."""

    TUMOUR_SUPPRESSOR = "tumour_suppressor"
    """Disrupted by loss, which is exactly what an unintended indel or stop does."""


@dataclass(frozen=True, slots=True)
class GeneLists:
    """Membership by symbol, with the source that supplied it."""

    source: Source
    members: dict[GeneList, frozenset[str]] = field(default_factory=dict)
    citations: dict[str, str] = field(default_factory=dict)

    @property
    def pin(self) -> str:
        return self.source.pin

    def lists_for(self, symbol: str | None) -> tuple[GeneList, ...]:
        if symbol is None:
            return ()
        return tuple(
            name for name, symbols in self.members.items() if symbol.upper() in symbols
        )

    def contains(self, symbol: str | None, name: GeneList) -> bool:
        return name in self.lists_for(symbol)

    def cite(self, symbol: str) -> str:
        return self.citations.get(symbol.upper(), self.source.pin)


def load_gene_lists(path: str | Path, version: str = "unversioned") -> GeneLists:
    """Read a symbol/list file and pin it."""
    source = Source.of("gene_lists", Path(path), version)
    members: dict[GeneList, set[str]] = {}
    citations: dict[str, str] = {}

    for row in read_tsv(source, required={SYMBOL_COLUMN, LIST_COLUMN}):
        symbol = row[SYMBOL_COLUMN].strip().upper()
        if not symbol:
            continue
        try:
            name = GeneList(row[LIST_COLUMN].strip())
        except ValueError as error:
            raise ContextError(
                f"{symbol}: {row[LIST_COLUMN]!r} is not a list this package knows; "
                f"it knows {', '.join(entry.value for entry in GeneList)}"
            ) from error
        members.setdefault(name, set()).add(symbol)
        if row.get(CITATION_COLUMN):
            citations[symbol] = row[CITATION_COLUMN]

    if not members:
        raise ContextError(f"{Path(path).name} lists no genes")
    return GeneLists(
        source=source,
        members={name: frozenset(symbols) for name, symbols in members.items()},
        citations=citations,
    )
