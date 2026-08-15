"""Ranking somebody else's hit list by where the hits land.

The claim under test is one sentence long: a hit list sorted by mismatch count
is sorted by the wrong thing. The fixture is built to make that fail loudly if
it is not true — its worst hit has *more* mismatches than several of its
harmless ones, so any ranking that reduces to arithmetic over mismatch counts
puts them in the wrong order and a test here goes red.

The second thing under test is the shape of not knowing. No annotation, no gene
lists, no expression release: each absence has to produce a different and
visibly weaker answer than the presence would, rather than a confident "low".
"""

from __future__ import annotations

from pathlib import Path

import pytest

from repairbench.annotation.gff import parse_gff3
from repairbench.annotation.store import Placement, TranscriptStore
from repairbench.cli import main
from repairbench.context.expression import Tissue, ingest
from repairbench.context.genelists import GeneList, load_gene_lists
from repairbench.context.source import ContextError, Provenance, Source
from repairbench.design.editors import DesignError
from repairbench.design.offtarget import OffTargetHit, read_casoffinder
from repairbench.design.report import render_offtarget
from repairbench.design.risk import RiskTier, assess, build_hit_features, load_risk_rules

DATA = Path(__file__).parent / "data"
HITS = DATA / "design" / "hits.txt"
LISTS = DATA / "design" / "gene_lists.tsv"
GFF = DATA / "mini.gff3"
MATRIX = DATA / "context" / "gtex_median_tpm.tsv"
RULES = Path(__file__).parents[1] / "rules" / "offtarget-v1.yaml"

MUSCLE = Tissue("Muscle - Skeletal")


@pytest.fixture
def hits():
    return read_casoffinder(HITS)


@pytest.fixture
def rules():
    return load_risk_rules(RULES)


@pytest.fixture
def lists():
    return load_gene_lists(LISTS, "fixture")


@pytest.fixture
def locate():
    return TranscriptStore(parse_gff3(GFF)).locate


@pytest.fixture
def expression():
    collected: dict[str, Provenance] = {}
    ingest(Source.of("expression", MATRIX, "fixture"), collected)
    return {
        gene: dict(provenance.facts["expression"].value)
        for gene, provenance in collected.items()
    }


def review(hits, rules, **kwargs):
    return assess(hits, rules, source_pin="fixture", **kwargs)


def tier_of(assessment, chromosome: str, position: int) -> RiskTier:
    return next(
        hit.tier
        for hit in assessment.hits
        if hit.hit.chromosome == chromosome and hit.hit.position == position
    )


# --------------------------------------------------------------------------
# Reading the file
# --------------------------------------------------------------------------


def test_the_hits_are_read_as_the_search_wrote_them(hits):
    assert len(hits) == 6
    assert hits[0].mismatches == 0
    assert {hit.strand for hit in hits} == {"+", "-"}


def test_a_file_in_another_layout_is_refused(tmp_path: Path):
    """A column read as the wrong field produces coordinates that look fine and
    point somewhere else."""
    path = tmp_path / "wrong.txt"
    path.write_text("GUIDE\t17\t100\n")

    with pytest.raises(DesignError, match="expected 6 columns"):
        read_casoffinder(path)


def test_an_empty_search_is_not_read_as_a_safe_guide(tmp_path: Path):
    path = tmp_path / "empty.txt"
    path.write_text("# nothing here\n")

    with pytest.raises(DesignError, match="not the same as a safe guide"):
        read_casoffinder(path)


def test_one_guide_can_be_selected_from_a_multi_guide_file(tmp_path: Path):
    """Mixing two guides' hits into one review would attribute one guide's worst
    site to the other."""
    path = tmp_path / "two.txt"
    path.write_text(
        "AAAACCCCGGGGTTTTAAAANRG\t17\t100\tAAAACCCCGGGGTTTTAAAAAGG\t+\t0\n"
        "TTTTGGGGCCCCAAAATTTTNRG\t17\t500\tTTTTGGGGCCCCAAAATTTTAGG\t+\t1\n"
    )

    selected = read_casoffinder(path, guide="AAAACCCCGGGGTTTTAAAA")

    assert len(selected) == 1
    assert selected[0].position == 100


# --------------------------------------------------------------------------
# The gene lists
# --------------------------------------------------------------------------


def test_membership_is_read_per_list(lists):
    assert lists.contains("PLUSG", GeneList.ESSENTIAL)
    assert lists.contains("DEMOG", GeneList.TUMOUR_SUPPRESSOR)
    assert not lists.contains("PLUSG", GeneList.ONCOGENE)


def test_a_gene_absent_from_the_file_is_absent_from_the_list_only(lists):
    """Not from the genome, and not from mattering. The rules are written to
    know the difference, and this is the assertion that keeps them honest."""
    assert lists.lists_for("SOMETHING_ELSE") == ()


def test_a_list_this_package_does_not_know_stops_the_run(tmp_path: Path):
    path = tmp_path / "lists.tsv"
    path.write_text("symbol\tlist\nGENEX\tinteresting\n")

    with pytest.raises(ContextError, match="not a list this package knows"):
        load_gene_lists(path)


# --------------------------------------------------------------------------
# The ranking
# --------------------------------------------------------------------------


def test_a_distant_hit_in_an_essential_gene_outranks_a_close_one_in_nothing(
    hits, rules, lists, locate
):
    """The whole argument, in one assertion. Four mismatches into the coding
    sequence of an essential gene, against five mismatches into intergenic
    space — and the first is the one that stops the guide."""
    assessment = review(hits, rules, locate=locate, lists=lists)

    assert tier_of(assessment, "17", 1550) is RiskTier.PROHIBITIVE
    assert tier_of(assessment, "11", 5000) is RiskTier.LOW
    assert assessment.hits[0].hit.position in {545, 1550, 4100}


def test_a_perfect_match_elsewhere_is_prohibitive_whatever_it_hits(hits, rules, locate):
    """The guide is not unique in this genome, and no property of the
    surrounding sequence redeems that."""
    assessment = review(hits, rules, locate=locate)

    assert tier_of(assessment, "17", 545) is RiskTier.PROHIBITIVE
    assert any(
        fired.rule_id == "A_HIT_WITH_NO_MISMATCHES"
        for hit in assessment.hits
        for fired in hit.evidence
    )


def test_a_tumour_suppressor_outranks_an_oncogene(hits, rules, lists, locate):
    """Deliberate, and the reason is directional: a random indel disrupts, and
    disruption is how a tumour suppressor causes harm and not how an oncogene
    does."""
    assessment = review(hits, rules, locate=locate, lists=lists)

    assert tier_of(assessment, "17", 4100) is RiskTier.PROHIBITIVE
    assert tier_of(assessment, "11", 450) is RiskTier.SERIOUS


def test_the_worst_tier_a_hit_earns_is_the_one_it_gets(hits, rules, lists, locate):
    """Several rules fire on one hit and they disagree about severity. Reasons
    to worry do not average out."""
    assessment = review(hits, rules, locate=locate, lists=lists)
    perfect = next(hit for hit in assessment.hits if hit.hit.position == 545)

    assert len(perfect.evidence) > 1
    assert perfect.tier is RiskTier.PROHIBITIVE


def test_every_ranked_hit_carries_the_rule_that_ranked_it(hits, rules, lists, locate):
    assessment = review(hits, rules, locate=locate, lists=lists)

    for hit in assessment.hits:
        if hit.tier is not RiskTier.LOW:
            assert hit.evidence
            assert all(fired.because for fired in hit.evidence)


# --------------------------------------------------------------------------
# Where the tissue dimension earns its place
# --------------------------------------------------------------------------


def test_a_non_coding_hit_in_a_gene_silent_in_the_target_tissue_is_downgraded(
    hits, rules, lists, locate, expression
):
    """PLUSG is measured at 0.1 TPM in skeletal muscle. An intronic hit there is
    a different finding from the same hit in the tissue where the gene works —
    downgraded, and deliberately not dismissed, because an edit is permanent and
    an expression programme is not."""
    assessment = review(hits, rules, locate=locate, lists=lists, expression=expression,
                        tissue=MUSCLE)
    intronic = next(hit for hit in assessment.hits if hit.hit.position == 700)

    assert intronic.tier is RiskTier.MODERATE
    assert any(
        fired.rule_id == "HIT_IN_A_GENE_SILENT_IN_THE_TARGET_TISSUE"
        for fired in intronic.evidence
    )


def test_without_a_tissue_the_expression_rule_cannot_fire(hits, rules, lists, locate, expression):
    assessment = review(hits, rules, locate=locate, lists=lists, expression=expression)

    assert not any(
        fired.rule_id == "HIT_IN_A_GENE_SILENT_IN_THE_TARGET_TISSUE"
        for hit in assessment.hits
        for fired in hit.evidence
    )


# --------------------------------------------------------------------------
# The shape of not knowing
# --------------------------------------------------------------------------


def test_without_an_annotation_nothing_is_ranked_and_the_report_says_so(hits, rules):
    """A hit list with no annotation behind it has been read, not reviewed. The
    tier is a separate word from "low" so that nobody can read one as the
    other."""
    assessment = review(hits, rules)

    unplaced = [hit for hit in assessment.hits if not hit.evidence]
    assert unplaced and all(hit.tier is RiskTier.UNASSESSED for hit in unplaced)
    assert "has been read, not one that has been reviewed" in assessment.unranked


def test_a_missing_gene_list_is_not_read_as_clearance(hits, rules, locate):
    """Without the lists, the essential-gene rule cannot fire — and the hit still
    ranks, on the generic coding rule, rather than falling through to low."""
    assessment = review(hits, rules, locate=locate)

    assert tier_of(assessment, "17", 1550) is not RiskTier.LOW
    assert not any(
        fired.rule_id == "CODING_HIT_IN_AN_ESSENTIAL_GENE"
        for hit in assessment.hits
        for fired in hit.evidence
    )


def test_features_that_could_not_be_computed_are_none_not_false():
    hit = OffTargetHit("GUIDE", "17", 100, "ACGT", "+", 2)
    features = build_hit_features(hit, Placement(gene="GENEX", in_transcript_span=True))

    assert features.get("hit.gene_is_essential") is None
    assert features.get("hit.expression_measured") is False
    assert features.get("hit.in_coding_sequence") is False


def test_no_sequence_score_is_invented(hits, rules, locate):
    """CFD weights a mismatch by its position and identity, from a table this
    package does not carry. The report says what the ranking is, so nobody reads
    a context tier as a binding score."""
    assessment = review(hits, rules, locate=locate)

    assert "No CFD table is attached" in assessment.scoring
    assert "not weighted by position" in assessment.scoring


def test_the_report_leads_with_the_worst_tier(hits, rules, lists, locate):
    rendered = render_offtarget(review(hits, rules, locate=locate, lists=lists))

    assert rendered.index("prohibitive") < rendered.index("serious")
    assert "CODING_HIT_IN_AN_ESSENTIAL_GENE" in rendered


def test_the_command_fails_the_build_on_a_prohibitive_hit(capsys):
    """The exit code is the point: a guide with a prohibitive off-target should
    stop a pipeline rather than be recorded in a log nobody reads."""
    assert (
        main(
            [
                "offtarget",
                str(HITS),
                "--annotation", str(GFF),
                "--gene-lists", str(LISTS),
            ]
        )
        == 1
    )
    assert "prohibitive" in capsys.readouterr().out
