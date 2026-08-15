"""Reading somebody else's off-target search.

Finding every site in a genome within six mismatches of a 20-mer, with bulges,
is a job for an indexed search over three gigabytes — Cas-OFFinder does it on a
GPU, and reimplementing it here in Python would be slower by orders of magnitude
and wrong in ways that are hard to see. So this module reads its output.

What this package adds is the part those tools leave out, and it is in
``risk.py``: a hit ranked by mismatch count alone is ranked by the wrong thing.
Three mismatches into the intron of a gene that is silent in the target tissue
is not the same finding as four mismatches into the coding sequence of a tumour
suppressor, and a list sorted by mismatch count puts them the wrong way round.

The file format is Cas-OFFinder's own, tab separated:

    crRNA  chromosome  position  matched DNA  strand  mismatches

with two extra columns (bulge type, bulge size) in the bulge-enabled output.
Both are accepted; a row that is neither is refused rather than guessed at,
because a column read as the wrong field produces coordinates that look fine.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from repairbench.design.editors import DesignError

#: Column counts Cas-OFFinder emits. Anything else is refused.
_CLASSIC, _WITH_BULGES = 6, 8


@dataclass(frozen=True, slots=True)
class OffTargetHit:
    """One site the search says resembles the guide."""

    guide: str
    chromosome: str
    #: 0-based in Cas-OFFinder's output; kept as the tool wrote it and converted
    #: once, here, rather than at every use.
    position: int
    sequence: str
    strand: str
    mismatches: int
    bulge_size: int = 0

    @property
    def key(self) -> str:
        return f"{self.chromosome}:{self.position}:{self.strand}"

    @property
    def has_bulge(self) -> bool:
        return self.bulge_size > 0


def _same_guide(crrna: str, guide: str) -> bool:
    """Is this row's crRNA the guide we were asked about?

    Compared by prefix in either direction, because a Cas-OFFinder query carries
    the PAM pattern appended in IUPAC code — ``...AAAANRG`` for a guide written
    ``...AAAA`` — and an exact match would never fire while a substring test
    would quietly match a guide that merely contains another.
    """
    left, right = crrna.upper(), guide.upper()
    return left.startswith(right) or right.startswith(left)


def read_casoffinder(path: str | Path, *, guide: str | None = None) -> tuple[OffTargetHit, ...]:
    """Parse a Cas-OFFinder output file.

    ``guide`` restricts to one crRNA. A file produced for several guides at once
    is the normal case, and silently mixing hits for two different guides into
    one risk assessment would attribute one guide's worst site to the other.
    """
    source = Path(path)
    hits: list[OffTargetHit] = []

    for number, line in enumerate(source.read_text().splitlines(), start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) not in {_CLASSIC, _WITH_BULGES}:
            raise DesignError(
                f"{source.name}:{number}: expected {_CLASSIC} columns (or {_WITH_BULGES} with "
                f"bulges) from Cas-OFFinder, got {len(fields)}. A file in another layout would "
                "parse into plausible coordinates for the wrong sites"
            )

        crrna, chromosome, position, sequence, strand, mismatches = fields[:6]
        if guide is not None and not _same_guide(crrna, guide):
            continue
        try:
            hits.append(
                OffTargetHit(
                    guide=crrna,
                    chromosome=chromosome,
                    position=int(position),
                    sequence=sequence.upper(),
                    strand=strand,
                    mismatches=int(mismatches),
                    bulge_size=int(fields[7]) if len(fields) == _WITH_BULGES else 0,
                )
            )
        except ValueError as error:
            raise DesignError(f"{source.name}:{number}: {error}") from error

    if not hits:
        raise DesignError(
            f"{source.name} contains no hits for "
            f"{guide if guide else 'any guide'}. An empty search is not the same as a safe "
            "guide — check that the file is the one the search wrote"
        )
    return tuple(hits)
