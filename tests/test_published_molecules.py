"""The validation that is not self-referential: reproduce a molecule somebody made.

Every other test in this suite checks that the package agrees with itself. The
reference sets check that the rules reproduce what the field already concluded —
which is a real check, and it is still a check against *statements*. This file
checks the package against an *object*: a drug that exists, whose sequence is
printed on an FDA label, and which somebody manufactured and gave to patients.

Two subjects, one per designer.

**Eteplirsen** (EXONDYS 51), a thirty-nucleotide phosphorodiamidate
morpholino that makes the spliceosome skip *DMD* exon 51. Three properties make
it the right one to spend this on:

* its sequence is public and unambiguous, quoted verbatim in the label rather
  than reconstructed from a figure;
* *DMD* is on the **minus strand**, which is where an antisense design goes
  wrong invisibly;
* it targets an exon interior rather than a splice site, so it also tests that
  the rules caution rather than refuse when a molecule does something the rule
  file did not anticipate.

**It failed the first time it was run, and that is the point.** ``tile`` was
reverse-complementing every window unconditionally, so for a minus-strand gene
it returned the sequence of the messenger instead of something complementary to
it — a molecule identical to its own target, which hybridises with nothing.
Every antisense oligonucleotide this package had ever printed for *DMD*,
*COL1A1* or *MECP2* was that. Nothing caught it: the synthetic fixtures have no
orientation, and ``test_real_locus`` had written the expected oligonucleotide
down as a constant and then never asserted it.

**FAH c.1062+5G>A** does the same for the base editor, and its result is
stranger and better. Given nothing but a patient allele and a reference
chromosome, this package independently reaches the published correction guide's
editor class, PAM, strand and position-in-protospacer — and differs from the
printed sequence at exactly one base. That base turns out to be a bystander edit
in the authors' cell line, and it is the bystander *this package predicts* for
the guide they used to build that cell line. The disagreement was the package
being right about something nobody had told it.

Skips without ``refdata/``. Needs the real chromosomes and the real annotation —
a fixture cannot test this, because the whole claim is that the answer matches
something outside this repository.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from repairbench.annotation.fasta import IndexedFasta
from repairbench.annotation.gff import parse_gff3
from repairbench.annotation.store import TranscriptStore
from repairbench.design.aso import Exon, tile
from repairbench.design.designer import CorrectionRequest, design
from repairbench.design.editors import load_editors
from repairbench.design.flags import load_flag_rules
from repairbench.design.prime import EditRequest, design_pegrnas
from repairbench.design.sequence import reverse_complement
from repairbench.design.sequence import reverse_complement as revcomp

ROOT = Path(__file__).parents[1]
GENOME = ROOT / "refdata" / "chrX.fa"
ANNOTATION = ROOT / "refdata" / "GRCh38_latest_genomic.gff.gz"
RULES = ROOT / "rules" / "aso-v1.yaml"

#: EXONDYS 51 (eteplirsen), FDA label 206488, DESCRIPTION section, quoted
#: verbatim: the 5'→3' base sequence of the thirty subunits.
ETEPLIRSEN = "CTCCAACATCAAGGAAGATGGCATTTCTAG"

#: Where those thirty bases turn out to sit in GRCh38. Not taken from any
#: paper — found by searching the assembly for the label's sequence, then
#: checked against the exon the annotation gives. Written down so that a
#: reference or annotation change that moves it fails here loudly.
EXON_51 = (31773960, 31774192)
ETEPLIRSEN_SPAN = (31774098, 31774127)

pytestmark = pytest.mark.skipif(
    not (GENOME.exists() and ANNOTATION.exists()),
    reason="needs refdata/ — run scripts/fetch-reference-data.sh",
)


@pytest.fixture(scope="module")
def annotation() -> TranscriptStore:
    """Both genes in one parse.

    The real annotation is 1.5 GB and takes the better part of twenty seconds to
    read; parsing it once per test turned this file into the slowest thing in
    the suite, which is how a test stops being run."""
    return TranscriptStore(parse_gff3(ANNOTATION, genes={"DMD", "HEXA"}))


@pytest.fixture(scope="module")
def exon_51(annotation: TranscriptStore) -> tuple[int, int]:
    """DMD exon 51, as the real annotation places it."""
    record, _ = annotation.preferred_for("DMD")
    assert record.accession == "NM_004006.3"
    assert record.strand == "-"
    return record.cds_blocks[50]


@pytest.fixture(scope="module")
def candidates(exon_51: tuple[int, int]):
    start, end = exon_51
    rules = load_flag_rules(RULES)
    with IndexedFasta(GENOME) as genome:
        target = genome.fetch("chrX", start, end).upper()
    outcome = tile(
        "DMD",
        "chrX",
        start,
        target,
        rules,
        chemistry="steric-PMO-30",
        strand="-",
        exon=Exon(start=start, end=end),
    )
    return outcome.candidates


def test_the_annotation_puts_exon_51_where_the_drug_needs_it(exon_51):
    """233 nucleotides, not a multiple of three — which is why skipping it
    restores the reading frame for the deletions eteplirsen is indicated for."""
    start, end = exon_51

    assert (start, end) == EXON_51
    assert end - start + 1 == 233
    assert (end - start + 1) % 3 != 0


def test_the_drugs_sequence_is_in_the_reference_where_expected(exon_51):
    """Before asking whether the designer finds it: is it there at all?"""
    with IndexedFasta(GENOME) as genome:
        found = genome.fetch("chrX", *ETEPLIRSEN_SPAN).upper()

    assert found == ETEPLIRSEN
    start, end = exon_51
    assert start <= ETEPLIRSEN_SPAN[0] and ETEPLIRSEN_SPAN[1] <= end


def test_the_designer_produces_eteplirsen_base_for_base(candidates):
    """The whole file, in one assertion."""
    matches = [candidate for candidate in candidates if candidate.sequence == ETEPLIRSEN]

    assert len(matches) == 1, (
        "the approved molecule is not among the windows this package would design "
        f"for its own target; {len(candidates)} windows were produced"
    )
    assert matches[0].span == ETEPLIRSEN_SPAN


def test_the_molecule_is_not_the_transcript_it_binds(candidates):
    """The defect this file was written to catch, stated directly.

    For a minus-strand gene the messenger is the reverse complement of the
    forward sequence. An oligonucleotide carrying *that* is a copy of its own
    target: perfectly designed, and inert."""
    sequences = {candidate.sequence for candidate in candidates}

    assert ETEPLIRSEN in sequences
    assert reverse_complement(ETEPLIRSEN) not in sequences


def test_the_rules_caution_about_the_approved_drug_rather_than_refusing_it(candidates):
    """An honest result that is worth keeping rather than tuning away.

    Eteplirsen binds an exon-internal enhancer, not a splice site, and the rule
    file says a steric blocker away from any splice signal is worth a second
    look. That caution is correct as written and the drug works anyway — which
    is exactly why it is a caution and not a veto, and why nothing in this
    package ranks candidates for a customer."""
    drug = next(c for c in candidates if c.sequence == ETEPLIRSEN)
    flags = {flag.rule_id for flag in drug.flags}

    assert "A_STERIC_BLOCKER_AWAY_FROM_ANY_SPLICE_SIGNAL" in flags
    assert not any(flag.severity.value == "blocking" for flag in drug.flags)


def test_the_melting_temperature_says_it_is_out_of_range(candidates):
    """Reproducing a thirty-mer is what exposed this: the Wallace rule returns
    something near 90 °C for one, which is not a melting temperature."""
    drug = next(c for c in candidates if c.sequence == ETEPLIRSEN)

    assert "TM_IS_OUTSIDE_WHAT_THIS_APPROXIMATION_COVERS" in {
        flag.rule_id for flag in drug.flags
    }


def test_a_strand_must_be_given(exon_51):
    """No default, because the defect was a silent assumption of one."""
    start, end = exon_51
    with IndexedFasta(GENOME) as genome:
        target = genome.fetch("chrX", start, end).upper()

    with pytest.raises(Exception, match="not a strand"):
        tile(
            "DMD",
            "chrX",
            start,
            target,
            load_flag_rules(RULES),
            chemistry="steric-PMO-30",
            strand="unknown",
        )


# --------------------------------------------------------------------------
# Base editing: a published correction for a real disease variant
# --------------------------------------------------------------------------

CHR15 = ROOT / "refdata" / "chr15.fa"
EDITORS = ROOT / "rules" / "editors-v1.yaml"

#: Hereditary tyrosinemia type 1, FAH c.1062+5G>A — the fifth base of the
#: intron 12 donor. Both guides below are quoted from a published worked
#: example: one installs the variant in a cell line with BE3, the other
#: corrects it with ABE7.10, and both put the site at protospacer position 5.
FAH_VARIANT = 80180230
INSTALLING_GUIDE = "GATACTCACCGGCCCGCTGA"      # BE3, PAM TGG, minus strand
PUBLISHED_CORRECTION = "GTAAATATCTGGCTGCACTG"  # ABE7.10, PAM AGG, plus strand

fah = pytest.mark.skipif(
    not (CHR15.exists() and ANNOTATION.exists()),
    reason="needs refdata/chr15.fa — run scripts/fetch-reference-data.sh",
)


@pytest.fixture(scope="module")
def genome15():
    with IndexedFasta(CHR15) as fasta:
        yield fasta


@pytest.fixture(scope="module")
def editors():
    return load_editors(EDITORS)


def correct(genome15, editors, patient: str, wild_type: str):
    return design(
        CorrectionRequest(
            gene="FAH",
            chromosome="chr15",
            position=FAH_VARIANT,
            patient_base=patient,
            wild_type_base=wild_type,
        ),
        genome15,
        editors,
    )


@fah
def test_the_annotation_places_the_variant_at_the_donor(genome15):
    """c.1062+5 should be the fifth base of an intron, and the intron should
    open with the canonical GT."""
    store = TranscriptStore(parse_gff3(ANNOTATION, genes={"FAH"}))
    record, _ = store.preferred_for("FAH")

    assert record.accession == "NM_000137.4"
    assert record.strand == "+"
    last_coding = record.genomic_position(1062)
    assert last_coding + 5 == FAH_VARIANT
    assert genome15.fetch("chr15", last_coding + 1, last_coding + 2).upper() == "GT"
    assert genome15.fetch("chr15", FAH_VARIANT, FAH_VARIANT).upper() == "G"


@fah
def test_the_designer_reaches_the_published_correction(genome15, editors):
    """Editor class, PAM, strand and position in the protospacer, independently.

    Nothing about the published guide is given to the designer — it is handed a
    patient allele and a reference chromosome, and asked what would put that
    base in an adenine editor's window."""
    outcome = correct(genome15, editors, patient="A", wild_type="G")
    at_five = [c for c in outcome.candidates if c.target_position_in_protospacer == 5]

    assert at_five, "nothing places the variant at the position the published guide uses"
    ours = at_five[0]
    assert ours.pam.upper() == "AGG"
    assert ours.strand == "+"
    assert any(c.editor.id.startswith("ABE7.10") for c in at_five)
    assert len(ours.protospacer) == len(PUBLISHED_CORRECTION)


@fah
def test_the_one_base_of_disagreement_is_a_bystander_this_package_predicts(
    genome15, editors
):
    """The result worth the whole exercise.

    Our correction guide differs from the published one at exactly one base:
    position 3, where the paper has A and the reference has G. That is not an
    error in either. The published guide was written against a *cell line* in
    which BE3 had installed the variant — and BE3's window covers a second
    cytosine, which it edits too. The authors' correction guide therefore
    matches an allele carrying the variant plus a bystander; a patient carries
    only the variant.

    The check: run the *installing* guide through this package and see whether
    the bystander it predicts is the base in question."""
    ours = next(
        c
        for c in correct(genome15, editors, patient="A", wild_type="G").candidates
        if c.target_position_in_protospacer == 5
    )
    differences = [
        (index + 1, published, mine)
        for index, (published, mine) in enumerate(
            zip(PUBLISHED_CORRECTION, ours.protospacer, strict=True)
        )
        if published != mine
    ]
    assert differences == [(3, "A", "G")]

    installing = next(
        c
        for c in correct(genome15, editors, patient="G", wild_type="A").candidates
        if c.protospacer == INSTALLING_GUIDE
    )
    assert installing.pam.upper() == "TGG"
    assert installing.strand == "-"
    assert installing.target_position_in_protospacer == 5

    predicted = {(b.position_in_protospacer, b.genomic_position) for b in installing.bystanders}
    assert predicted == {(7, FAH_VARIANT - 2)}, (
        "the bystander this package predicts for the installing guide is exactly the base "
        "by which the published correction guide differs from ours"
    )


# --------------------------------------------------------------------------
# Prime editing: the pegRNA from the paper that introduced it
# --------------------------------------------------------------------------

#: Anzalone et al. 2019, Nature 576:149, Supplementary Information, "Figure 5
#: sequences": pegRNA ``HEXA_5b_correct``, which reverts the Ashkenazi
#: Tay-Sachs allele HEXA c.1274_1277dupTATC (long written 1278insTATC).
#:
#: The spacers in that table carry a 5' G added for U6 transcription; the
#: protospacer is the remaining 20.
HEXA_SPACER = "ATCCTTCCAGTCAGGGCCAT"
HEXA_EXTENSION = "ACCTGAACCGTATATCCTATGGCCCTGACTG"
HEXA_PBS = HEXA_EXTENSION[-10:]        # the table's "PBS length (nt)" column: 10
HEXA_RTT = HEXA_EXTENSION[:21]         # the table's "RT template length (nt)": 21
HEXA_NICKING = "TACCTGAACCGTATATCCTA"  # nicking gRNA GTACCTGAACCGTATATCCTA, less its G

#: The duplication, as a VCF-style insertion on the forward strand. HEXA is on
#: the minus strand, so the four coding bases TATC are GATA here.
HEXA_ANCHOR = 72346579
HEXA_PATIENT = "GGATA"
HEXA_WILD_TYPE = "G"

prime = pytest.mark.skipif(
    not (CHR15.exists() and ANNOTATION.exists()),
    reason="needs refdata/chr15.fa — run scripts/fetch-reference-data.sh",
)


@pytest.fixture(scope="module")
def pegrnas(genome15):
    return design_pegrnas(
        EditRequest(
            gene="HEXA",
            chromosome="chr15",
            position=HEXA_ANCHOR,
            patient_allele=HEXA_PATIENT,
            wild_type_allele=HEXA_WILD_TYPE,
        ),
        genome15,
        load_flag_rules(ROOT / "rules" / "prime-v1.yaml"),
    )


@prime
def test_the_duplication_is_where_the_transcript_says_it_is(genome15, annotation):
    """c.1274_1277 should read TATC in the transcript, which is what makes the
    Ashkenazi allele a tandem duplication rather than a plain insertion."""
    record, _ = annotation.preferred_for("HEXA")
    assert record.accession == "NM_000520.6"
    assert record.strand == "-"

    first, last = record.genomic_position(1277), record.genomic_position(1274)
    assert revcomp(genome15.fetch("chr15", first, last).upper()) == "TATC"
    assert record.genomic_position(1278) == HEXA_ANCHOR


@prime
def test_the_designer_reaches_the_published_pegrna_exactly(pegrnas):
    """Spacer, PAM, strand, nick, primer binding site and template — all of it.

    Nothing from the paper is supplied: the module is given a patient allele
    and chromosome 15, and every component below is derived from the assembly."""
    exact = [
        c
        for c in pegrnas.candidates
        if c.spacer == HEXA_SPACER and c.pbs == HEXA_PBS and c.rtt == HEXA_RTT
    ]

    assert len(exact) == 1
    found = exact[0]
    assert found.pam == "AGG"
    assert found.strand == "+"
    assert found.extension == HEXA_EXTENSION
    assert found.nick_position == 72346574
    assert found.nick_to_edit_nt == 5


@prime
def test_the_published_nicking_guide_is_offered_and_called_pe3b(pegrnas):
    """The paper's second nick is the install pegRNA's own spacer, which matches
    the *corrected* allele — so it cannot fire until the edit is made. That is
    the definition of PE3b, and this package has to reach it unaided."""
    found = next(
        c
        for c in pegrnas.candidates
        if c.spacer == HEXA_SPACER and c.pbs == HEXA_PBS and c.rtt == HEXA_RTT
    )
    published = [g for g in found.nicking_guides if g.protospacer == HEXA_NICKING]

    assert len(published) == 1
    assert published[0].edit_dependent, "the published nicking guide is the PE3b kind"


@prime
def test_an_insertion_does_not_silently_consume_reference_bases(genome15, pegrnas):
    """The defect this reproduction found, stated as a property.

    Both the patient and the edited sequence used to splice out ``len(patient)``
    reference bases instead of ``len(wild_type)``. For a substitution those are
    equal and nothing showed; for an insertion — which is what prime editing is
    for — the "patient" sequence came out as the untouched reference and the
    "edited" sequence carried a deletion nobody asked for. The primer binding
    site still matched the paper, because it is read upstream of the nick; the
    template did not."""
    found = next(c for c in pegrnas.candidates if c.rtt == HEXA_RTT)
    written = revcomp(found.rtt)
    reference = genome15.fetch("chr15", found.nick_position + 1, found.nick_position + 21).upper()

    assert written == reference, (
        "the template must write the reference sequence back, since the correction "
        "removes the duplicated bases"
    )
