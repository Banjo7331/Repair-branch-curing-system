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

from repairbench.annotation.fasta import IndexedFasta, InMemorySequences, write_index
from repairbench.annotation.gff import AnnotationError, parse_gff3
from repairbench.annotation.naming import NO_ALIASES
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


# --------------------------------------------------------------------------
# Writing the index, so a real genome does not also need samtools
# --------------------------------------------------------------------------


def test_the_index_written_here_is_the_one_samtools_would_have_written(tmp_path: Path):
    """Checked against the committed fixture indexes, which were written by
    hand from the format's definition. If these ever diverge, every coordinate
    read through the index moves."""
    for fixture in (DATA / "mini.fa", DATA / "design" / "target.fa"):
        copy = tmp_path / fixture.name
        copy.write_bytes(fixture.read_bytes())

        assert write_index(copy).read_text() == Path(f"{fixture}.fai").read_text()


def test_an_indexed_file_reads_back_the_sequence_it_indexed(tmp_path: Path):
    path = tmp_path / "two.fa"
    path.write_text(">chr1 with a description\nACGTACGTAC\nGTACGTACGT\nAC\n>chr2\nTTTTTTTT\n")
    write_index(path)

    with IndexedFasta(path) as fasta:
        assert fasta.fetch("chr1", 1, 22) == "ACGTACGTACGTACGTACGTAC"
        assert fasta.fetch("chr1", 11, 12) == "GT"
        assert fasta.fetch("chr2", 8, 8) == "T"


def test_a_ragged_sequence_is_refused_rather_than_indexed(tmp_path: Path):
    """The format cannot describe one, and an index written anyway returns
    plausible sequence from the wrong coordinates."""
    path = tmp_path / "ragged.fa"
    path.write_text(">chr1\nACGTAC\nAC\nACGTAC\n")

    with pytest.raises(AnnotationError, match="ragged"):
        write_index(path)


def test_something_that_is_not_a_fasta_is_refused(tmp_path: Path):
    path = tmp_path / "notes.txt"
    path.write_text("this is not a FASTA\n")

    with pytest.raises(AnnotationError, match="is this a FASTA"):
        write_index(path)


# --------------------------------------------------------------------------
# The two names every chromosome has
# --------------------------------------------------------------------------


def test_the_alias_table_is_read_from_the_annotations_own_records(tmp_path: Path):
    """NCBI opens each chromosome with a region record carrying its ordinary
    name. Reading the mapping out of the file beats a hard-coded table, which
    goes stale with every assembly patch."""
    path = tmp_path / "real.gff3"
    path.write_text(
        "##gff-version 3\n"
        "NC_000017.11\tRefSeq\tregion\t1\t83257441\t.\t+\t.\t"
        "ID=NC_000017.11:1..83257441;chromosome=17;genome=chromosome\n"
        "NC_000017.11\tBestRefSeq\tmRNA\t500\t2000\t.\t+\t.\t"
        "ID=rna-NM_000088.4;gene=COL1A1;transcript_id=NM_000088.4;tag=MANE Select\n"
        "NC_000017.11\tBestRefSeq\tCDS\t500\t619\t.\t+\t0\tID=cds-1;Parent=rna-NM_000088.4\n"
    )

    annotation = parse_gff3(path)

    assert annotation.aliases.canonical("chr17") == "NC_000017.11"
    assert annotation.aliases.canonical("17") == "NC_000017.11"
    assert annotation.aliases.same_sequence("NC_000017.11", "chr17")


def test_a_variant_written_the_ucsc_way_resolves_against_a_refseq_annotation(tmp_path: Path):
    """The collision a real run meets immediately: the annotation says
    NC_000017.11 and the FASTA says chr17. Without this, every real coordinate
    is refused as being on the wrong chromosome."""
    path = tmp_path / "real.gff3"
    path.write_text(
        "##gff-version 3\n"
        "NC_000017.11\tRefSeq\tregion\t1\t83257441\t.\t+\t.\t"
        "ID=NC_000017.11:1..83257441;chromosome=17;genome=chromosome\n"
        "NC_000017.11\tBestRefSeq\tmRNA\t500\t2000\t.\t+\t.\t"
        "ID=rna-NM_000088.4;gene=COL1A1;transcript_id=NM_000088.4;tag=MANE Select\n"
        "NC_000017.11\tBestRefSeq\tCDS\t500\t619\t.\t+\t0\tID=cds-1;Parent=rna-NM_000088.4\n"
    )
    store = TranscriptStore(parse_gff3(path))

    for spelling in ("chr17", "17", "NC_000017.11"):
        assert store.resolve("COL1A1", spelling, 550).cds_position == 51


def test_an_unplaced_scaffold_does_not_claim_a_chromosomes_name(tmp_path: Path):
    """Scaffolds carry a chromosome attribute too. Without the genome=chromosome
    check, a patch contig would alias itself to the chromosome it patches, and
    coordinates would resolve against the wrong sequence."""
    path = tmp_path / "scaffold.gff3"
    path.write_text(
        "##gff-version 3\n"
        "NT_187614.1\tRefSeq\tregion\t1\t100000\t.\t+\t.\t"
        "ID=NT_187614.1:1..100000;chromosome=17;genome=genomic\n"
        "NT_187614.1\tBestRefSeq\tmRNA\t500\t2000\t.\t+\t.\t"
        "ID=rna-NM_1.1;gene=GENEX;transcript_id=NM_1.1\n"
        "NT_187614.1\tBestRefSeq\tCDS\t500\t619\t.\t+\t0\tID=cds-1;Parent=rna-NM_1.1\n"
    )

    assert parse_gff3(path).aliases.canonical("17") == "17"


def test_the_chr_prefix_alone_resolves_without_any_annotation():
    """A fixture that declares no region records still has to work, and the
    prefix is the one difference that needs no evidence."""
    assert NO_ALIASES.same_sequence("17", "chr17")
    assert not NO_ALIASES.same_sequence("17", "18")


def test_a_fasta_written_the_other_way_still_answers(tmp_path: Path):
    """UCSC writes chr17, Ensembl writes 17, and a run mixing the two sources
    should not fail on a prefix."""
    path = tmp_path / "one.fa"
    path.write_text(">chr17\nACGTACGTAC\n")
    write_index(path)

    with IndexedFasta(path) as fasta:
        assert fasta.fetch("17", 1, 4) == fasta.fetch("chr17", 1, 4) == "ACGT"


def test_an_unknown_chromosome_still_names_what_the_file_has(tmp_path: Path):
    path = tmp_path / "one.fa"
    path.write_text(">chr17\nACGTACGTAC\n")
    write_index(path)

    with IndexedFasta(path) as fasta, pytest.raises(AnnotationError, match="chr17"):
        fasta.fetch("chr9", 1, 4)


# --------------------------------------------------------------------------
# One transcript's biology must not take the file down
# --------------------------------------------------------------------------


def test_overlapping_cds_blocks_drop_the_transcript_not_the_file(tmp_path: Path):
    """*PEG10* is the case that taught this package the difference. It is
    translated through a programmed ribosomal frameshift — the ribosome slips
    back a base and reads on in another frame — so two of its CDS blocks
    legitimately share a coordinate. RefSeq annotates that correctly and this
    package cannot use it, because a position inside the overlap has two CDS
    offsets. Refusing the whole file, which is what this used to do, let one
    gene deny access to the other hundred and thirty-six thousand transcripts."""
    path = tmp_path / "frameshift.gff3"
    path.write_text(
        "##gff-version 3\n"
        "1\ttest\tmRNA\t100\t900\t.\t+\t.\tID=rna-NM_000111.1;gene=SLIPPY;transcript_id=NM_000111.1\n"
        "1\ttest\tCDS\t100\t300\t.\t+\t0\tID=c1;Parent=rna-NM_000111.1\n"
        "1\ttest\tCDS\t300\t600\t.\t+\t0\tID=c2;Parent=rna-NM_000111.1\n"
        "1\ttest\tmRNA\t100\t900\t.\t+\t.\tID=rna-NM_000222.1;gene=ORDINARY;transcript_id=NM_000222.1\n"
        "1\ttest\tCDS\t100\t300\t.\t+\t0\tID=c3;Parent=rna-NM_000222.1\n"
        "1\ttest\tCDS\t400\t600\t.\t+\t0\tID=c4;Parent=rna-NM_000222.1\n"
    )

    annotation = parse_gff3(path)

    assert "NM_000222.1" in annotation.transcripts
    assert "NM_000111.1" not in annotation.transcripts
    assert "ribosomal frameshift" in annotation.rejected["NM_000111.1"]


def test_a_rejected_transcript_is_not_silently_missing(tmp_path: Path):
    """A reader deciding whether a gene is absent or unusable has to be able to
    tell, so the reason travels with the annotation rather than to a log."""
    path = tmp_path / "frameshift.gff3"
    path.write_text(
        "##gff-version 3\n"
        "1\ttest\tmRNA\t100\t900\t.\t-\t.\tID=rna-NM_000111.1;gene=SLIPPY;transcript_id=NM_000111.1\n"
        "1\ttest\tCDS\t100\t300\t.\t-\t0\tID=c1;Parent=rna-NM_000111.1\n"
        "1\ttest\tCDS\t250\t600\t.\t-\t0\tID=c2;Parent=rna-NM_000111.1\n"
    )

    annotation = parse_gff3(path)

    assert not annotation.transcripts
    assert list(annotation.rejected) == ["NM_000111.1"]

    with pytest.raises(AnnotationError, match="no coding transcript"):
        TranscriptStore(annotation).preferred_for("SLIPPY")
