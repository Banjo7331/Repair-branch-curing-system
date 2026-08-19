"""Reconciling the two names every chromosome has.

Point this package at real files and the first thing that happens is a name
collision. NCBI's GFF3 calls chromosome 17 ``NC_000017.11``. UCSC's FASTA calls
it ``chr17``. A VCF from a clinical pipeline calls it ``17`` or ``chr17``
depending on who wrote it. All four are the same sequence.

The tempting fix is a hard-coded table, and it is wrong twice over: it goes
stale with every assembly patch, and it silently invents an answer for anything
it does not know. The better source is the annotation itself. NCBI's GFF3 opens
each chromosome with a ``region`` record carrying ``chromosome=17`` in its
attributes, so the mapping between the accession and the ordinary name is *in
the file*, written by the people who assigned both.

What this module will not do is guess. ``chr17`` and ``17`` and
``NC_000017.11`` resolve to one another because the annotation says so or
because the prefix rule is unambiguous; anything else raises. A coordinate
looked up under a name that resolved to the wrong contig would succeed, return
plausible sequence, and be wrong in a way nothing downstream could detect.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field

_CHR_PREFIX = "chr"


def bare(name: str) -> str:
    """A chromosome name with the ``chr`` prefix removed, if it had one.

    The one normalisation safe to do without consulting anything: ``chr17`` and
    ``17`` differ by a prefix that carries no information, and every reference
    source disagrees about whether to write it.
    """
    return name[len(_CHR_PREFIX) :] if name.lower().startswith(_CHR_PREFIX) else name


@dataclass(frozen=True, slots=True)
class Aliases:
    """Every name one sequence goes by, read from the annotation that named it.

    ``by_alias`` maps each known spelling to the accession the annotation uses,
    so a lookup can be answered without knowing which convention the caller
    prefers. It is deliberately not a two-way dictionary: the accession is the
    canonical form here because it is the one the GFF3's own records carry.
    """

    by_alias: Mapping[str, str] = field(default_factory=dict)

    @classmethod
    def of(cls, pairs: Iterable[tuple[str, str]]) -> Aliases:
        """Build from (accession, ordinary name) pairs — ``NC_000017.11``, ``17``."""
        table: dict[str, str] = {}
        for accession, name in pairs:
            for alias in (accession, name, f"{_CHR_PREFIX}{bare(name)}"):
                table[alias.lower()] = accession
        return cls(by_alias=table)

    def is_chromosome(self, seqid: str) -> bool:
        """Is this accession one of the assembled chromosomes?

        The alternative is a scaffold, a patch or an unplaced contig, and real
        annotation is full of them: a gene in a segmentally duplicated region
        has transcripts on several at once. Which of those a caller should be
        given is a judgement, and it cannot be made without first knowing which
        is which.
        """
        return seqid in set(self.by_alias.values())

    def canonical(self, name: str) -> str:
        """The accession this name refers to, or the name itself when unknown.

        Returning the name unchanged rather than raising is deliberate: an
        annotation of a single synthetic contig has no region records and needs
        no aliases, and that case must keep working exactly as it did.
        """
        return self.by_alias.get(name.lower(), name)

    def same_sequence(self, one: str, other: str) -> bool:
        """Do two names refer to the same sequence?

        Three ways they can: identical, aliases of one accession, or differing
        only by the ``chr`` prefix. The third is needed because a file may use a
        convention the annotation never mentioned.
        """
        if one == other:
            return True
        if self.canonical(one) == self.canonical(other):
            return True
        return bare(one).lower() == bare(other).lower()

    def __len__(self) -> int:
        return len(self.by_alias)


#: An empty table, for annotations that name only one contig and need no
#: reconciling. Its ``same_sequence`` still handles the ``chr`` prefix, which is
#: the one difference that needs no evidence to resolve.
NO_ALIASES = Aliases()
