"""Reading the one file that describes the patient.

Most of these tests are about genotype, because genotype is where zygosity comes
from and zygosity decides roughly half the modalities. The rest are about the
ways a VCF reader is quietly wrong rather than loudly broken: the wrong sample,
an unnormalised indel, a consequence invented rather than read.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from repairbench.annotation.fasta import IndexedFasta, InMemorySequences
from repairbench.annotation.normalise import left_align
from repairbench.model import Consequence, Zygosity
from repairbench.vcf import VcfError, VcfReader, zygosity_from_genotype

DATA = Path(__file__).parent / "data"
CASE = DATA / "deployment" / "case.vcf"
FASTA = DATA / "mini.fa"


@pytest.fixture
def genome():
    with IndexedFasta(FASTA) as fasta:
        yield fasta


def read(sample: str = "CHILD", **kwargs):
    return VcfReader(path=CASE, sample=sample, **kwargs).read()


# --------------------------------------------------------------------------
# Genotype
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("genotype", "expected"),
    [
        ("0/1", Zygosity.HETEROZYGOUS),
        ("0|1", Zygosity.HETEROZYGOUS),
        ("1|0", Zygosity.HETEROZYGOUS),
        ("1/1", Zygosity.HOMOZYGOUS),
        ("1", Zygosity.HEMIZYGOUS),
        ("./.", Zygosity.UNKNOWN),
        (".", Zygosity.UNKNOWN),
    ],
)
def test_the_ordinary_genotypes(genotype: str, expected: Zygosity):
    assert zygosity_from_genotype(genotype) is expected


def test_two_different_alternate_alleles_leave_no_reference_copy():
    """The case worth knowing about. ``1/2`` is heterozygous for each allele and
    carries no reference allele at all — which, for every rule that asks whether
    an intact copy exists, is the same answer as a homozygote."""
    zygosity = zygosity_from_genotype("1/2")

    assert zygosity is Zygosity.COMPOUND_HETEROZYGOUS
    assert zygosity.leaves_a_wild_type_allele is False


def test_a_partial_no_call_is_read_from_what_was_called():
    """``./1`` says one allele is alternate and the other is unknown. Reading it
    as heterozygous would invent a reference allele nobody observed."""
    assert zygosity_from_genotype("./1") is Zygosity.HOMOZYGOUS


def test_unknown_is_not_a_guess():
    """Zygosity that was not called stays uncalled: the modality layer raises a
    caveat for it, which is the honest handling of missing data."""
    assert zygosity_from_genotype("./.").leaves_a_wild_type_allele is None


# --------------------------------------------------------------------------
# Records
# --------------------------------------------------------------------------


def test_only_alleles_the_sample_carries_are_returned():
    """A multi-allelic site lists what the cohort has, not what this child has."""
    variants = read()

    at_900 = [v for v in variants if v.position == 900]
    assert len(at_900) == 1
    assert at_900[0].alternate == "T"  # the sample is 0/2; allele 1 is somebody else's


def test_a_failed_filter_is_dropped_by_default():
    """A failed filter is the caller saying it does not believe its own call."""
    assert all(v.passed_filters for v in read())
    assert any(not v.passed_filters for v in read(require_pass=False))


def test_reading_the_wrong_sample_is_refused_rather_than_guessed():
    """Picking the first column would silently interpret a parent's genotype as
    the child's, and every downstream answer would be about the wrong person."""
    with pytest.raises(VcfError, match="Name the one to read"):
        VcfReader(path=CASE).read()


def test_an_unknown_sample_name_lists_what_is_there():
    with pytest.raises(VcfError, match="MOTHER"):
        VcfReader(path=CASE, sample="FATHER").read()


def test_the_mother_and_the_child_are_different_people():
    child = {v.key: v.zygosity for v in read("CHILD")}
    mother = {v.key: v.zygosity for v in read("MOTHER")}

    assert child != mother


# --------------------------------------------------------------------------
# Normalisation
# --------------------------------------------------------------------------


def test_indels_are_left_aligned_when_a_reference_is_supplied(genome):
    """Decomposition without normalisation is wrong: a variant written
    non-canonically fails to match the same variant written properly, including
    the patient's own earlier report."""
    normalised = VcfReader(path=CASE, sample="CHILD").read(genome)

    keys = {v.key for v in normalised}
    assert "17-999-G-GA" in keys  # written at 1005 in the file
    assert "17-3299-G-GT" in keys  # written at 3305
    assert all(v.normalised for v in normalised)


def test_without_a_reference_the_reader_says_so_rather_than_pretending():
    """The flag is on the record, so a caller can tell whether the key it is
    about to join on is canonical."""
    assert all(not v.normalised for v in read())


def test_normalisation_does_not_disturb_a_substitution(genome):
    read_back = VcfReader(path=CASE, sample="CHILD").read(genome)
    substitutions = [v for v in read_back if v.position == 560]

    assert substitutions[0].key == "17-560-A-T"


def test_a_variant_at_the_start_of_a_contig_still_normalises():
    """A homopolymer at position 1 has nothing to its left to borrow."""
    tiny = InMemorySequences({"1": "AAAAAAAA"})

    assert left_align("1", 4, "A", "AA", tiny).position == 1


# --------------------------------------------------------------------------
# Annotation
# --------------------------------------------------------------------------


def test_gene_and_consequence_are_read_not_derived():
    """This project reads consequences; it does not predict them. Guessing
    'missense because the alleles are the same length' would be a prediction
    dressed as a parse."""
    by_key = {v.key: v for v in read()}

    assert by_key["17-560-A-T"].consequence is Consequence.NONSENSE
    assert by_key["17-560-A-T"].gene == "PLUSG"


def test_the_annotation_is_matched_to_the_allele_it_describes():
    """A multi-allelic CSQ carries one entry per allele, and attaching the wrong
    one would give the right gene with the wrong consequence."""
    by_key = {v.key: v for v in read()}

    assert by_key["17-900-C-T"].consequence is Consequence.MISSENSE


def test_the_most_severe_consequence_wins(tmp_path: Path):
    """VEP lists several per transcript and the order is its own. Taking the
    first would make the answer depend on how VEP was invoked."""
    path = tmp_path / "multi.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        "17\t100\t.\tA\tT\t50\tPASS\t"
        "CSQ=T|splice_region_variant&stop_gained|GENEX|ENSG1\tGT\t0/1\n"
    )

    variant = VcfReader(path=path, sample="S").read()[0]

    assert variant.consequence is Consequence.NONSENSE


def test_a_record_without_an_annotation_is_not_an_error(tmp_path: Path):
    """Most records in a genome carry nothing this project can use. That is not
    a parse failure; it is a record the rules have no question about."""
    path = tmp_path / "bare.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        "17\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\n"
    )

    variant = VcfReader(path=path, sample="S").read()[0]

    assert variant.consequence is None
    assert not variant.is_interpretable


# --------------------------------------------------------------------------
# The assembly
# --------------------------------------------------------------------------


def test_a_mismatched_assembly_is_caught_from_the_header(tmp_path: Path):
    """Every coordinate would be plausible and every answer would be about the
    wrong part of the genome."""
    path = tmp_path / "wrong.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n##reference=file:///ref/hg19.fa\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    )

    with pytest.raises(VcfError, match="does not look like GRCh38"):
        VcfReader(path=path, sample="S", expected_assembly="GRCh38").read()


def test_the_declared_assembly_is_accepted():
    VcfReader(path=CASE, sample="CHILD", expected_assembly="GRCh38").read()


def test_data_before_the_header_is_refused(tmp_path: Path):
    path = tmp_path / "headerless.vcf"
    path.write_text("17\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\n")

    with pytest.raises(VcfError, match="before the #CHROM header"):
        VcfReader(path=path, sample="S").read()


def test_a_spanning_deletion_placeholder_is_not_a_variant(tmp_path: Path):
    """``*`` is a statement about another record, not a variant of its own."""
    path = tmp_path / "star.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        "17\t100\t.\tA\t*,T\t50\tPASS\t.\tGT\t1/2\n"
    )

    variants = VcfReader(path=path, sample="S").read()

    assert [v.alternate for v in variants] == ["T"]
