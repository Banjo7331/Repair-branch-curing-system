"""Left-alignment: giving a variant the same name everybody else gives it.

An insertion of ``A`` into a run of eight ``A``s can be written nine ways, all
of them describing the same change. Every database picks the leftmost, and a
variant written any other way silently fails to match ClinVar, fails to match
the patient's earlier report, and — in this package — fails to land on the exon
the NMD arithmetic thinks it lands on.

M5 refused to do this and said so, because trimming shared affixes is safe with
no external input while left-shifting is not: shifting requires the reference
sequence. Here the reference is in hand, so the refusal is lifted, and the
algorithm is the one bcftools and vt implement.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from repairbench.annotation.fasta import SequenceProvider
from repairbench.annotation.gff import AnnotationError

_DNA = re.compile(r"^[ACGTN]+$")

#: A guard against a pathological or corrupt reference sending the shift loop
#: leftwards forever. Real indels shift by a handful of bases; a homopolymer run
#: long enough to exhaust this is a data problem, not a variant.
MAX_SHIFT = 1000


@dataclass(frozen=True, slots=True)
class NormalisedVariant:
    """A variant in its canonical, leftmost, minimal representation."""

    chromosome: str
    position: int
    reference: str
    alternate: str
    #: How far the representation moved left. Zero means the variant was already
    #: canonical. A negative value means trimming a shared prefix moved it
    #: *right*, which is still normalisation. Anything non-zero means an earlier
    #: record of the same change, written the other way, would not have matched.
    shifted_by: int = 0

    @property
    def key(self) -> str:
        return f"{self.chromosome}-{self.position}-{self.reference}-{self.alternate}"

    @property
    def is_indel(self) -> bool:
        return len(self.reference) != len(self.alternate)


def left_align(
    chromosome: str,
    position: int,
    reference: str,
    alternate: str,
    sequences: SequenceProvider,
) -> NormalisedVariant:
    """Trim and left-shift a variant against the reference.

    The algorithm is the one described by Tan, Abecasis and Kang (2015) and
    implemented by bcftools and vt, and the detail that makes it correct is the
    representation: internally an allele may be *empty*. An insertion is an
    empty reference and a one-base alternate, not a padded pair, and the loop
    then reads as two moves that alternate cleanly —

    * if both alleles are non-empty and end with the same base, drop it;
    * if either allele is empty, borrow the reference base to the left.

    A substitution satisfies neither and falls straight through, which is why
    normalising a whole VCF costs nothing on the 99% of records that are already
    canonical. An insertion inside a homopolymer run walks left one base at a
    time until the borrowed base finally differs, which is precisely the leftmost
    representation every database stores.

    Padding is restored at the end, so what comes out is a valid VCF record
    rather than the internal form.
    """
    _require_dna("reference", reference)
    _require_dna("alternate", alternate)
    if position < 1:
        raise AnnotationError(f"position must be 1-based, got {position}")
    if reference == alternate:
        raise AnnotationError(f"reference and alternate are identical at {chromosome}:{position}")

    original_position = position

    # Drop a shared prefix first. This is the cheap half of normalisation and is
    # the only half M5 was willing to do without a reference in hand.
    while reference and alternate and reference[0] == alternate[0]:
        reference, alternate = reference[1:], alternate[1:]
        position += 1

    shifts = 0
    while True:
        if reference and alternate and reference[-1] == alternate[-1]:
            reference, alternate = reference[:-1], alternate[:-1]
            continue
        if reference and alternate:
            break
        if position <= 1:
            # Nothing further left to borrow. This is the canonical form.
            break
        preceding = sequences.fetch(chromosome, position - 1, position - 1)
        if not preceding or preceding == "N":
            # An unplaced or masked base is not something to shift into.
            break
        reference, alternate = preceding + reference, preceding + alternate
        position -= 1
        shifts += 1
        if shifts > MAX_SHIFT:
            raise AnnotationError(
                f"{chromosome}:{original_position} shifted more than {MAX_SHIFT} bases; "
                "the reference or the variant is wrong"
            )

    if not reference or not alternate:
        # VCF alleles cannot be empty: pad with the base to the left, or with
        # the base to the right at the very start of a chromosome.
        if position > 1:
            pad = sequences.fetch(chromosome, position - 1, position - 1)
            reference, alternate = pad + reference, pad + alternate
            position -= 1
        else:
            pad = sequences.fetch(chromosome, position, position)
            reference, alternate = reference + pad, alternate + pad

    return NormalisedVariant(
        chromosome=chromosome,
        position=position,
        reference=reference,
        alternate=alternate,
        shifted_by=original_position - position,
    )


def _require_dna(label: str, allele: str) -> None:
    if not _DNA.match(allele):
        raise AnnotationError(f"{label} allele is not a DNA string: {allele!r}")


def verify_reference(
    chromosome: str,
    position: int,
    reference: str,
    sequences: SequenceProvider,
) -> None:
    """Check that the reference allele is what the reference genome says.

    The cheapest possible detection of the most damaging possible error: a VCF
    called against a different assembly than the annotation. Every coordinate
    would be plausible, every lookup would succeed, and every answer would be
    about the wrong part of the genome.
    """
    actual = sequences.fetch(chromosome, position, position + len(reference) - 1)
    if actual != reference.upper():
        raise AnnotationError(
            f"{chromosome}:{position} claims reference {reference!r} but the genome has "
            f"{actual!r} — the variant and the reference are from different assemblies"
        )
