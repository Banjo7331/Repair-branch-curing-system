"""Counting pathogenic variation, which is where the most inferential rule feeds.

The clustering rule — pathogenic missense piling into one stretch of protein
while truncating variants are absent — is how this package tells gain of
function from loss of function without being told. Until this module existed,
the counts it reads were typed in by hand, so the least direct inference in the
package rested on the least direct evidence there is: numbers somebody
estimated.

The tests here are mostly about what must *not* be counted. A count is a single
number that hides everything that went into it, and three different mistakes all
produce a plausible one: matching "Pathogenic" as a substring so that
"Conflicting classifications of pathogenicity" counts as support, taking
one-submitter records and expert-panel records as the same evidence, and reading
ClinVar's Type column — which describes the sequence change — as though it
described the protein.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from repairbench.cli import main
from repairbench.context.clinvar import (
    ReviewStatus,
    VariantKind,
    classify,
    commonest_transcript,
    distribution_for,
    exemplars,
    group_by_gene,
    ingest,
    read_variants,
    review_summary,
)
from repairbench.context.source import ContextError, Provenance, Source

DATA = Path(__file__).parent / "data" / "context"
SUMMARY = DATA / "clinvar_variant_summary.txt"


def source() -> Source:
    return Source.of("clinvar", SUMMARY, "2026-08")


def variants(**kwargs):
    return read_variants(source(), genes={"COL1A1"}, **kwargs)


# --------------------------------------------------------------------------
# What must not be counted
# --------------------------------------------------------------------------


def test_conflicting_is_not_pathogenic():
    """ClinVar's own aggregate says when submitters disagree, and the phrase for
    it contains the word "pathogenicity". A substring match counts a
    disagreement as support."""
    assert not any(variant.coding == "c.500G>A" for variant in variants())


def test_a_submission_with_no_criteria_is_below_the_floor():
    """Zero stars. Whether that is acceptable is a policy, and the policy is a
    threshold rather than a hard-coded filter — but it is not the default."""
    assert not any(variant.coding == "c.600G>A" for variant in variants())
    assert any(variant.coding == "c.600G>A" for variant in variants(minimum_stars=0))


def test_the_other_assembly_is_not_a_second_variant():
    """The file carries both builds, so the same variant appears twice. Counting
    both would double every gene's evidence."""
    assert len([v for v in variants() if v.coding == "c.2461G>A"]) == 1
    assert all(v.chromosome == "17" for v in variants())


def test_the_gene_filter_is_required():
    """The real file is millions of rows. A call without a filter would read all
    of them to answer a question about a handful, and would look like it was
    working."""
    with pytest.raises(ContextError, match="millions of rows"):
        read_variants(source(), genes=set())


# --------------------------------------------------------------------------
# What kind of change it is
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("protein", "kind"),
    [
        ("p.Gly821Ser", VariantKind.MISSENSE),
        ("p.Arg334Ter", VariantKind.TRUNCATING),
        ("p.Gly1000fs", VariantKind.TRUNCATING),
        ("p.Ala33=", VariantKind.SYNONYMOUS),
        ("", VariantKind.OTHER),
    ],
)
def test_the_kind_is_read_from_the_protein_change(protein: str, kind: VariantKind):
    """Not from ClinVar's Type column, which says how the *sequence* changed. A
    single nucleotide variant can be missense or nonsense, and those two argue
    for opposite mechanisms."""
    assert classify(protein, "single nucleotide variant") is kind


def test_a_frameshift_counts_as_truncating_whatever_its_sequence_type():
    frameshift = next(v for v in variants() if "fs" in v.protein)

    assert frameshift.kind is VariantKind.TRUNCATING


def test_a_synonymous_submission_is_neither_missense_nor_truncating():
    distribution = distribution_for(variants())
    counted = distribution.pathogenic_missense_total + distribution.pathogenic_truncating_total

    assert counted == len([v for v in variants() if v.kind is not VariantKind.SYNONYMOUS])


# --------------------------------------------------------------------------
# The hotspot, which is a window and not a domain
# --------------------------------------------------------------------------


def test_clustering_is_the_densest_window_not_a_bin():
    """Three missense variants at residues 821, 824 and 827 sit inside twenty
    residues of each other. Binning by fixed edges would split them whenever the
    edge happened to fall between."""
    distribution = distribution_for(variants(), hotspot_window_aa=20)

    assert distribution.pathogenic_missense_total == 3
    assert distribution.pathogenic_missense_in_hotspot == 3


def test_a_narrower_window_finds_less_clustering():
    """The width is a judgement, so it is a parameter, and moving it must
    visibly move the answer rather than quietly not."""
    assert distribution_for(variants(), hotspot_window_aa=3).pathogenic_missense_in_hotspot < 3


def test_a_gene_with_no_missense_has_no_hotspot():
    truncating_only = [v for v in variants() if v.kind is VariantKind.TRUNCATING]

    assert distribution_for(truncating_only).pathogenic_missense_in_hotspot == 0


# --------------------------------------------------------------------------
# What the provenance says
# --------------------------------------------------------------------------


def test_the_citation_says_what_was_counted_and_at_what_quality():
    """A distribution built from single-submitter records and one built from
    expert-panel records support the same rule to very different degrees, and
    the count alone hides which this is."""
    into: dict[str, Provenance] = {}
    ingest(source(), into, genes={"COL1A1", "PIK3CA"})

    citation = into["COL1A1"].facts["distribution"].citation
    assert "★" in citation
    assert "densest 20-residue window" in citation


def test_the_star_breakdown_is_ordered_best_first():
    assert review_summary(variants()).startswith("1×3★")


def test_review_wording_this_release_does_not_use_scores_zero():
    """ClinVar rephrases these. A renamed status should reduce what is counted,
    never be scored as though it were reviewed."""
    assert ReviewStatus.parse("reviewed by a committee of one").stars == 0
    assert ReviewStatus.parse("reviewed by expert panel").stars == 3


def test_each_gene_gets_its_own_distribution():
    into: dict[str, Provenance] = {}
    ingest(source(), into, genes={"COL1A1", "PIK3CA"})

    assert into["PIK3CA"].facts["distribution"].value.pathogenic_missense_total == 2
    assert into["PIK3CA"].facts["distribution"].value.pathogenic_truncating_total == 0


# --------------------------------------------------------------------------
# Citing a position somebody reported
# --------------------------------------------------------------------------


def test_examples_are_restricted_to_one_transcript():
    """A c. position is only meaningful against the transcript it was written
    on, so a list mixing accessions would invite exactly the error this exists
    to prevent."""
    cited = exemplars(variants(), VariantKind.MISSENSE, transcript="NM_000088.4")

    assert cited
    assert {variant.transcript for variant in cited} == {"NM_000088.4"}


def test_examples_are_best_reviewed_first():
    cited = exemplars(variants(), VariantKind.TRUNCATING)

    assert [variant.review.stars for variant in cited] == sorted(
        (variant.review.stars for variant in cited), reverse=True
    )


def test_an_example_without_a_coding_name_is_not_offered():
    """A case cites c. positions. A record that carries none cannot supply one,
    and offering it with an empty string would put `c.` nowhere in a YAML file
    that reads as though it were somewhere."""
    assert all(variant.coding for variant in exemplars(variants(), VariantKind.MISSENSE))


def test_the_commonest_transcript_is_the_submitters_choice_not_ours():
    assert commonest_transcript(variants()) == "NM_000088.4"


def test_no_variants_names_no_transcript():
    assert commonest_transcript([]) == ""


def test_grouping_files_each_variant_under_its_own_gene():
    grouped = group_by_gene(read_variants(source(), genes={"COL1A1", "PIK3CA"}))

    assert set(grouped) == {"COL1A1", "PIK3CA"}
    assert all(v.gene == gene for gene, found in grouped.items() for v in found)


# --------------------------------------------------------------------------
# The command that puts the counts where the reference set can read them
# --------------------------------------------------------------------------


def test_the_command_prints_a_paste_ready_distribution(capsys):
    """The point of the command: updating the reference set after a release is
    transcription rather than judgement."""
    assert main(["clinvar", str(SUMMARY), "COL1A1", "--release", "fixture"]) == 0

    printed = capsys.readouterr().out
    assert "pathogenic_missense_total: 3" in printed
    assert "pathogenic_truncating_total: 2" in printed
    assert yaml.safe_load(
        printed.split("distribution:", 1)[1].split("\n", 1)[0]
    ) == {
        "pathogenic_missense_total": 3,
        "pathogenic_missense_in_hotspot": 3,
        "pathogenic_truncating_total": 2,
    }


def test_the_command_pins_the_file_it_counted(capsys):
    main(["clinvar", str(SUMMARY), "COL1A1", "--release", "fixture"])

    assert "clinvar@fixture/" in capsys.readouterr().out


def test_a_gene_that_matched_nothing_says_so_rather_than_reporting_zero(capsys):
    """A symbol ClinVar files under another name returns nothing, and printing
    that as zero pathogenic variation would invent a finding — one that argues
    for gain of function, since absent truncating variants is half that rule."""
    main(["clinvar", str(SUMMARY), "NOTAGENE"])

    printed = capsys.readouterr().out
    assert "nothing matched" in printed
    assert "pathogenic_missense_total" not in printed
