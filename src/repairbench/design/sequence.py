"""Sequence arithmetic: complements, IUPAC PAMs, and which strand is which.

Small and boring, and every function in it is somewhere a design can be quietly
wrong rather than loudly broken. A protospacer placed on the wrong strand still
looks like a protospacer; a PAM matched without expanding its IUPAC codes still
looks matched; an N in the reference still looks like a base.
"""

from __future__ import annotations

from repairbench.design.editors import DesignError

_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

#: IUPAC ambiguity codes, expanded. PAMs are written in this alphabet — NGG,
#: NNGRRT, NRN — and matching them by string equality would find nothing.
IUPAC: dict[str, str] = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "AG",
    "Y": "CT",
    "S": "CG",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}

UNAMBIGUOUS = frozenset("ACGT")


def reverse_complement(sequence: str) -> str:
    """The other strand, read 5' to 3'.

    Half of every PAM scan happens here: a protospacer on the minus strand runs
    against increasing genomic coordinates, and its PAM sits at the *lower*
    coordinates. Working in reverse-complement space rather than with two sets
    of index arithmetic is what keeps that from being a source of off-by-ones.
    """
    return sequence.translate(_COMPLEMENT)[::-1]


def complement_base(base: str) -> str:
    return base.translate(_COMPLEMENT)


def matches_pam(sequence: str, pam: str) -> bool:
    """Does this sequence satisfy an IUPAC PAM pattern?

    An unknown base (``N`` in the *reference*, not in the pattern) never
    matches. A design against unresolved reference sequence would be a design
    against a guess, and the caller is told the site is unusable rather than
    handed a candidate resting on it.
    """
    if len(sequence) != len(pam):
        return False
    for base, code in zip(sequence, pam, strict=True):
        allowed = IUPAC.get(code)
        if allowed is None:
            raise DesignError(f"PAM {pam!r} contains {code!r}, which is not an IUPAC code")
        if base not in allowed or base not in UNAMBIGUOUS:
            return False
    return True


def is_resolved(sequence: str) -> bool:
    """True when every base is A, C, G or T.

    Reference genomes carry N runs, and a protospacer overlapping one cannot be
    ordered as an oligonucleotide, let alone scored.
    """
    return bool(sequence) and set(sequence) <= UNAMBIGUOUS
