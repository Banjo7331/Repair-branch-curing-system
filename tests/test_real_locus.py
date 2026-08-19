"""One real variant, at a real coordinate, through the whole package.

Every other test in this suite runs against sequence this project invented. That
is the right way to test arithmetic — a fixture can be built so the answer is
knowable — and it is exactly why it cannot test the thing this file does: that
the arithmetic still holds when the input is a chromosome somebody else made.

The case is *COL1A1* p.(Gly821Ser), c.2461G>A. Three properties make it the
right one to spend a real genome on:

* The gene is on the **minus strand**, so every coordinate, every complement and
  every exon ordering is exercised in the direction where mistakes hide. A
  plus-strand gene would pass with half the code wrong.
* The change is G>A on the coding strand and therefore **C>T on the plus
  strand**, which is where a designer that reasons in genomic coordinates
  without thinking about strand produces a confident, wrong molecule.
* It is a glycine substitution in the triple helix — the textbook dominant
  negative, and the mechanism the whole package was built to tell apart from
  haploinsufficiency in the same gene.

The expected values below were verified by hand against GRCh38: the codon at the
mapped coordinate reads GGC on the coding strand, which is a glycine, and the
protospacer and oligonucleotide were read back out of the reference
independently of the code that produced them.

The test skips without ``refdata/``. Running it needs 100 MB of downloads that
are not ours to redistribute — ``scripts/fetch-reference-data.sh`` fetches them —
so CI runs everything else and this runs where the data is.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from repairbench.annotation.fasta import IndexedFasta
from repairbench.annotation.gff import parse_gff3
from repairbench.annotation.store import TranscriptStore
from repairbench.design.aso import tile
from repairbench.design.flags import load_flag_rules
from repairbench.design.sequence import reverse_complement
from repairbench.modality import Modality
from repairbench.model import Mechanism

ROOT = Path(__file__).parents[1]
REFDATA = ROOT / "refdata"
CASE = Path(__file__).parent / "data" / "real" / "col1a1-gly821ser.yaml"

GENOME = REFDATA / "chr17.fa"
ANNOTATION = REFDATA / "GRCh38_latest_genomic.gff.gz"

#: The variant, and what was verified by hand at that coordinate.
POSITION = 50190099
CDS_POSITION = 2461
CODON = "GGC"
PROTOSPACER = "GACAGCCAACCTGGTGCTAA"
#: The discriminating gapmer over the variant, carrying the patient's T where
#: the reference has C. COL1A1 is on the minus strand, so the molecule *is* the
#: forward sequence of that window.
#:
#: This constant read ``AGCCAACCTGGTGCTAAAGG`` — the reverse complement, the
#: sequence of the transcript rather than of anything that binds it — and was
#: never asserted. It sat in the block headed "verified by hand" while the one
#: check that would have caught the strand defect was the check nobody wrote.
OLIGONUCLEOTIDE = "CCTTTAGCACCAGGTTGGCT"
OLIGONUCLEOTIDE_SPAN = (50190080, 50190099)

pytestmark = pytest.mark.skipif(
    not (GENOME.exists() and ANNOTATION.exists()),
    reason="needs refdata/ — run scripts/fetch-reference-data.sh",
)


@pytest.fixture(scope="module")
def record():
    store = TranscriptStore(parse_gff3(ANNOTATION, genes={"COL1A1"}))
    chosen, _ = store.preferred_for("COL1A1")
    return chosen


@pytest.fixture(scope="module")
def genome():
    with IndexedFasta(GENOME) as fasta:
        yield fasta


@pytest.fixture(scope="module")
def case():
    return yaml.safe_load(CASE.read_text())


def test_the_transcript_is_the_one_the_clinic_would_use(record):
    assert record.accession == "NM_000088.4"
    assert record.strand == "-"
    assert len(record.cds_blocks) == 51


def test_the_coordinate_and_the_hgvs_agree(record, case):
    """The case names both, and they were written independently: the genomic
    position came from the annotation, the c. position from the literature. If
    the CDS arithmetic were wrong on a minus-strand gene, they would disagree."""
    assert case["genomic"]["position"] == POSITION
    assert record.cds_offset(POSITION) == CDS_POSITION
    assert record.genomic_position(CDS_POSITION) == POSITION


def test_the_codon_at_that_coordinate_is_a_glycine(record, genome):
    """The claim the case has been making since before there was a genome to
    check it against. c.2461 is the first base of its codon, the gene reads
    right to left, and the codon comes out GGC."""
    codon = reverse_complement(genome.fetch("chr17", POSITION - 2, POSITION))

    assert codon == CODON
    assert codon.startswith("GG")  # every glycine codon does


def test_the_reference_carries_the_wild_type_base(genome, case):
    assert genome.fetch("chr17", POSITION, POSITION) == case["genomic"]["reference"]


def test_the_protospacer_reads_back_out_of_the_genome(genome):
    """Read independently of the designer: take the span, put the patient's
    allele in, complement it, and the result must be the guide."""
    start, end = 50190083, 50190102
    span = genome.fetch("chr17", start, end)
    patient = span[: POSITION - start] + "T" + span[POSITION - start + 1 :]

    assert reverse_complement(patient) == PROTOSPACER


def test_the_whole_plan_runs_and_reaches_the_right_mechanism(case):
    from repairbench.cli import build_query_with_provenance  # noqa: PLC0415
    from repairbench.engine import resolve  # noqa: PLC0415
    from repairbench.modality_rules import load_modality_ruleset  # noqa: PLC0415
    from repairbench.ruleset import load_ruleset  # noqa: PLC0415
    from repairbench.selector import select  # noqa: PLC0415

    query, notes = build_query_with_provenance(case, ANNOTATION, GENOME)
    call = resolve(query, load_ruleset(ROOT / "rules" / "mechanism-v1.yaml"))
    selection = select(call, query, load_modality_ruleset(ROOT / "rules" / "modality-v1.yaml"))

    assert call.mechanism is Mechanism.DOMINANT_NEGATIVE
    assert Modality.GENE_ADDITION in {a.modality for a in selection.contraindicated}
    assert Modality.ALLELE_SPECIFIC_SILENCING in {a.modality for a in selection.indicated}
    assert any("NM_000088.4" in note for note in notes)


def test_the_chromosome_name_in_the_case_is_not_the_one_in_the_annotation(case, record):
    """The case says chr17; RefSeq says NC_000017.11. They are the same sequence
    and nothing reconciles them by accident — the aliases come from the file's
    own region records, and this is the first place that mattered."""
    assert case["genomic"]["chromosome"] == "chr17"
    assert record.seqid == "NC_000017.11"


def test_an_annotation_this_size_rejects_a_few_transcripts_and_says_which():
    """PEG10 is translated through a programmed ribosomal frameshift, so two of
    its CDS blocks share a coordinate. That is real biology, correctly
    annotated, and unusable here — a position in the overlap has two CDS
    offsets. It is dropped with a reason rather than taking the file down."""
    annotation = parse_gff3(ANNOTATION, genes={"PEG10"})

    assert annotation.rejected
    assert any("ribosomal frameshift" in reason for reason in annotation.rejected.values())
    assert annotation.transcripts, "the gene's other transcripts must survive"


def test_the_discriminating_oligonucleotide_is_the_forward_strand(genome, record):
    """The value above, asserted at last.

    COL1A1 is transcribed from the minus strand, so an oligonucleotide
    complementary to its messenger carries the *forward* genomic sequence. The
    window must also end on the variant, because that base is the only thing
    telling the two alleles apart.
    """
    start, end = 50190039, 50190159
    region = genome.fetch("chr17", start, end).upper()
    offset = POSITION - start
    patient = region[:offset] + "T" + region[offset + 1 :]

    outcome = tile(
        "COL1A1",
        "chr17",
        start,
        patient,
        load_flag_rules(ROOT / "rules" / "aso-v1.yaml"),
        chemistry="gapmer-2MOE",
        strand=record.strand,
        must_cover=POSITION,
    )

    match = [c for c in outcome.candidates if c.span == OLIGONUCLEOTIDE_SPAN]
    assert len(match) == 1
    assert match[0].sequence == OLIGONUCLEOTIDE
    assert match[0].sequence.endswith("T"), "the patient's allele, not the reference's"
