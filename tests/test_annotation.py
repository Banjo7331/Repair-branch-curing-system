"""Reading real annotation, and the arithmetic that depends on getting it right.

The fixtures are synthetic and labelled as such — what is under test is the
parser, the strand handling and the coordinate mapping, none of which care
whether the coordinates correspond to a real locus. What they do care about is
being wrong in a way that raises rather than in a way that produces a plausible
number, and that is most of what these tests check.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from repairbench.annotation.fasta import IndexedFasta, InMemorySequences
from repairbench.annotation.gff import AnnotationError, parse_gff3
from repairbench.annotation.normalise import left_align, verify_reference
from repairbench.annotation.store import TranscriptStore

DATA = Path(__file__).parent / "data"
GFF = DATA / "mini.gff3"
FASTA = DATA / "mini.fa"


@pytest.fixture
def annotation():
    return parse_gff3(GFF)


@pytest.fixture
def store(annotation):
    return TranscriptStore(annotation)


@pytest.fixture
def genome():
    with IndexedFasta(FASTA) as fasta:
        yield fasta


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------


def test_coding_exons_are_read_in_transcript_order(annotation):
    record = annotation.transcripts["NM_000001.1"]

    assert record.coding_exon_lengths == (120, 90, 300)
    assert record.coding_length == 510
    assert record.mane_select


def test_a_minus_strand_transcript_is_reversed(annotation):
    """The first coding exon on the minus strand has the highest coordinate.

    Sorting CDS blocks by start and calling that transcript order inverts the
    exon numbering for half the genome — and inverts the NMD prediction with it,
    because the last junction ends up at the wrong end of the gene.
    """
    record = annotation.transcripts["NM_000002.2"]

    assert record.strand == "-"
    assert record.cds_blocks[0] == (841, 900)
    assert record.coding_exon_lengths == (60, 90, 150)


def test_a_non_coding_transcript_is_dropped(annotation):
    """An mRNA with no CDS records would otherwise become a transcript of
    coding length zero, and every downstream calculation would divide into it."""
    assert "NR_000005.1" not in annotation.transcripts


def test_the_annotation_pins_itself(annotation):
    assert annotation.pin.startswith("mini.gff3@")
    assert len(annotation.digest) == 64


def test_parsing_can_be_restricted_to_the_genes_that_matter():
    """The real file is gigabytes and the question concerns a handful of genes."""
    restricted = parse_gff3(GFF, genes={"PLUSG"})

    assert set(restricted.transcripts) == {"NM_000001.1"}


def test_asking_for_a_gene_that_is_not_there_fails_loudly():
    with pytest.raises(AnnotationError, match="no coding transcripts"):
        parse_gff3(GFF, genes={"NOSUCHGENE"})


# --------------------------------------------------------------------------
# Coordinates
# --------------------------------------------------------------------------


def test_genomic_to_cds_on_the_plus_strand(annotation):
    record = annotation.transcripts["NM_000001.1"]

    assert record.cds_offset(500) == 1
    assert record.cds_offset(619) == 120
    assert record.cds_offset(900) == 121  # first base of the second coding exon
    assert record.cds_offset(1500) == 211


def test_genomic_to_cds_on_the_minus_strand(annotation):
    """Counting runs the other way, and the first coding base is the last coordinate."""
    record = annotation.transcripts["NM_000002.2"]

    assert record.cds_offset(900) == 1
    assert record.cds_offset(841) == 60
    assert record.cds_offset(489) == 61


def test_an_intronic_coordinate_has_no_cds_offset(annotation):
    """Not an error in general; an error here, because every rule downstream
    reads a CDS offset and inventing one would place the variant in an exon it
    is not in."""
    assert annotation.transcripts["NM_000001.1"].cds_offset(700) is None


def test_the_coordinate_mapping_round_trips(annotation):
    for accession in ("NM_000001.1", "NM_000002.2"):
        record = annotation.transcripts[accession]
        for offset in range(1, record.coding_length + 1):
            assert record.cds_offset(record.genomic_position(offset)) == offset


# --------------------------------------------------------------------------
# Choosing a transcript
# --------------------------------------------------------------------------


def test_mane_select_is_preferred(store):
    record, reason = store.preferred_for("PLUSG")

    assert record.accession == "NM_000001.1"
    assert reason == "MANE Select"


def test_without_mane_the_longest_coding_sequence_wins_and_says_so(store):
    """A silent fallback would hide that a choice was made at all."""
    record, reason = store.preferred_for("ALTG")

    assert record.accession == "NM_000004.1"
    assert "no MANE Select" in reason


def test_resolving_places_a_genomic_variant_on_a_transcript(store):
    resolved = store.resolve("PLUSG", "17", 1500)

    assert resolved.record.accession == "NM_000001.1"
    assert resolved.cds_position == 211
    assert resolved.transcript.coding_exon_lengths == (120, 90, 300)


def test_a_variant_on_the_wrong_chromosome_is_refused(store):
    with pytest.raises(AnnotationError, match="is on 17"):
        store.resolve("PLUSG", "11", 1500)


def test_an_intronic_variant_is_refused_rather_than_placed(store):
    with pytest.raises(AnnotationError, match="not inside the coding sequence"):
        store.resolve("PLUSG", "17", 700)


def test_an_unversioned_accession_is_not_silently_matched(store):
    with pytest.raises(AnnotationError, match="version matters"):
        store.by_accession("NM_000001")


# --------------------------------------------------------------------------
# The reference genome
# --------------------------------------------------------------------------


def test_fasta_fetch_spans_line_breaks(genome):
    """Line breaks are bytes in the file and not bases in the sequence."""
    assert genome.fetch("17", 59, 62) == "GTAC"
    assert len(genome.fetch("17", 1, 300)) == 300


def test_fasta_knows_where_the_homopolymer_is(genome):
    assert genome.fetch("17", 1000, 1009) == "A" * 10
    assert genome.fetch("17", 999, 999) == "G"


def test_reading_past_the_end_is_refused(genome):
    with pytest.raises(AnnotationError, match="outside the sequence"):
        genome.fetch("17", 4999, 5100)


def test_normalisation_can_move_a_variant_into_a_different_answer(genome, store):
    """The case the whole package needed this module for.

    An insertion written at 17:3305 sits in one place; written canonically it
    sits six bases earlier. Both spellings describe the same change, both land
    in coding sequence, and only one of them matches what a database stores —
    so before this module existed, half the spellings would have failed to join
    against the patient's own earlier report.
    """
    normalised = left_align("17", 3305, "T", "TT", genome)

    assert normalised.position == 3299
    assert store.resolve("DEMOG", "17", 3305).cds_position == 306
    assert store.resolve("DEMOG", "17", normalised.position).cds_position == 300


def test_a_missing_index_is_reported_rather_than_built(tmp_path):
    """Building one means reading the whole file, which is what the index exists
    to avoid; doing it silently would turn a missing file into a mystery pause."""
    unindexed = tmp_path / "bare.fa"
    unindexed.write_text(">1\nACGT\n")

    with pytest.raises(AnnotationError, match=r"no \.fai index"):
        IndexedFasta(unindexed)


# --------------------------------------------------------------------------
# Normalisation
# --------------------------------------------------------------------------


def test_an_insertion_walks_to_the_left_of_its_homopolymer(genome):
    """The whole point. Inserting an A anywhere in a run of ten A's is one
    change with ten possible spellings, and databases store the leftmost."""
    result = left_align("17", 1005, "A", "AA", genome)

    assert (result.position, result.reference, result.alternate) == (999, "G", "GA")
    assert result.shifted_by == 6


def test_every_spelling_of_the_same_insertion_normalises_to_one_key(genome):
    """The property that makes a join on variant identity trustworthy."""
    keys = {left_align("17", pos, "A", "AA", genome).key for pos in range(1000, 1010)}

    assert keys == {"17-999-G-GA"}


def test_a_deletion_in_the_same_run_normalises_too(genome):
    result = left_align("17", 1004, "AA", "A", genome)

    assert (result.position, result.reference, result.alternate) == (999, "GA", "G")


def test_a_substitution_is_left_alone(genome):
    result = left_align("17", 1500, "A", "T", genome)

    assert (result.position, result.reference, result.alternate) == (1500, "A", "T")
    assert result.shifted_by == 0


def test_a_shared_prefix_is_trimmed(genome):
    """Trimming can move a variant right. That is still normalisation."""
    result = left_align("17", 1500, "AC", "AG", genome)

    assert (result.position, result.reference, result.alternate) == (1501, "C", "G")


def test_normalisation_is_idempotent(genome):
    once = left_align("17", 1005, "A", "AA", genome)
    twice = left_align(once.chromosome, once.position, once.reference, once.alternate, genome)

    assert once.key == twice.key


def test_the_shift_cannot_run_off_the_start_of_a_chromosome():
    """A homopolymer at position 1 has nothing to its left to borrow."""
    sequences = InMemorySequences({"1": "AAAAAAAAAA"})

    result = left_align("1", 5, "A", "AA", sequences)

    assert result.position == 1


def test_a_mismatched_reference_allele_is_caught(genome):
    """The cheapest detection of the most damaging error: a VCF called against a
    different assembly. Every coordinate would be plausible and every answer
    would be about the wrong part of the genome."""
    verify_reference("17", 1000, "A", genome)

    with pytest.raises(AnnotationError, match="different assemblies"):
        verify_reference("17", 1000, "T", genome)


def test_a_non_dna_allele_is_refused(genome):
    with pytest.raises(AnnotationError, match="not a DNA string"):
        left_align("17", 1000, "A", "<DEL>", genome)
