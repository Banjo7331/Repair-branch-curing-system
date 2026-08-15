"""Tiling antisense oligonucleotides, and the one confusion that inverts a therapy.

Most of this module is composition arithmetic, and most of these tests are about
the two places where an antisense design goes wrong in a way nothing else
catches.

The first is chemistry. A gapmer recruits RNase H and destroys what it binds; a
steric blocker occupies a site and destroys nothing. Aim a gapmer at a splice
site to redirect splicing and it degrades the transcript the redirection was
meant to rescue — the design is correct in every other respect and the therapy
is inverted.

The second is what a tiling run *is*. Two hundred oligonucleotides that pass
every composition rule are two hundred starting points, not two hundred
candidates, because the factor that decides which ones work — whether the site
is accessible inside the folded transcript — is not computed anywhere here.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from repairbench.annotation.fasta import IndexedFasta
from repairbench.cli import main
from repairbench.design.aso import (
    Action,
    Exon,
    NoStructureModel,
    Region,
    chemistries,
    longest_self_complement,
    tile,
)
from repairbench.design.editors import DesignError
from repairbench.design.flags import Severity, load_flag_rules
from repairbench.design.sequence import reverse_complement

DATA = Path(__file__).parent / "data" / "design"
FASTA = DATA / "target.fa"
RULES = Path(__file__).parents[1] / "rules" / "aso-v1.yaml"

#: The fixture's annotated exon, from tests/data/design/target.gff3.
EXON = Exon(start=250, end=350)


@pytest.fixture
def rules():
    return load_flag_rules(RULES)


@pytest.fixture
def target():
    with IndexedFasta(FASTA) as fasta:
        return fasta.fetch("17", 240, 360)


def run(target, rules, chemistry="steric-PMO", **kwargs):
    return tile("TARG", "17", 240, target, rules, chemistry=chemistry, **kwargs)


# --------------------------------------------------------------------------
# The chemistries
# --------------------------------------------------------------------------


def test_the_chemistries_are_data_with_an_action(rules):
    """Which one destroys and which one occupies is the fact everything else
    here depends on, so it is declared per chemistry rather than inferred from
    the name."""
    catalogue = chemistries(rules)

    assert catalogue["gapmer-2MOE"].action is Action.CLEAVES
    assert catalogue["steric-PMO"].action is Action.BLOCKS
    assert catalogue["steric-PMO"].length == 25


def test_a_chemistry_without_an_action_is_refused(tmp_path: Path):
    path = tmp_path / "aso.yaml"
    path.write_text(
        "version: bad\nchemistries:\n  - {id: mystery, length: 20}\n"
        "rules:\n  - {id: R, severity: note, when: {feature: aso.length_nt, gt: 0}, because: x}\n"
    )

    with pytest.raises(DesignError, match="not interchangeable"):
        chemistries(load_flag_rules(path))


def test_an_unknown_chemistry_lists_what_there_is(target, rules):
    with pytest.raises(DesignError, match="steric-PMO"):
        run(target, rules, chemistry="something-else")


def test_a_target_shorter_than_the_oligonucleotide_is_refused(rules):
    with pytest.raises(DesignError, match="no window to tile"):
        tile("TARG", "17", 1, "ACGT", rules, chemistry="steric-PMO")


# --------------------------------------------------------------------------
# The tiling
# --------------------------------------------------------------------------


def test_the_oligonucleotide_is_the_reverse_complement_of_its_target(target, rules):
    """Antisense. Reporting the target sequence as though it were the molecule
    to order is the most embarrassing possible bug here and the easiest to
    make."""
    outcome = run(target, rules)

    for candidate in outcome.candidates:
        assert candidate.sequence == reverse_complement(candidate.target)


def test_the_span_belongs_to_the_target_not_the_molecule(target, rules):
    """An oligonucleotide has no coordinates. The span is where it binds."""
    outcome = run(target, rules)
    candidate = outcome.candidates[0]

    assert candidate.span[1] - candidate.span[0] + 1 == candidate.chemistry.length
    assert 240 <= candidate.span[0] <= 360


def test_every_window_of_the_region_is_tiled(target, rules):
    outcome = run(target, rules)
    expected = len(target) - chemistries(rules)["steric-PMO"].length + 1

    assert outcome.tiled == expected


def test_the_step_comes_from_the_rule_file(target, rules, tmp_path: Path):
    """Tiling density is a trade between coverage and a list somebody can read,
    which makes it a policy rather than a constant."""
    coarse = tmp_path / "aso.yaml"
    coarse.write_text(RULES.read_text().replace("step_nt: 1", "step_nt: 5"))

    dense = run(target, rules).tiled
    sparse = run(target, load_flag_rules(coarse)).tiled

    assert sparse < dense


# --------------------------------------------------------------------------
# Where the window sits
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("start", "end", "region"),
    [
        (240, 264, Region.ACCEPTOR_SITE),
        (330, 354, Region.DONOR_SITE),
        (280, 304, Region.EXON_INTERIOR),
        (100, 124, Region.INTRON),
    ],
)
def test_a_window_is_placed_against_the_exon_boundaries(start, end, region):
    """A splice acceptor and an exon interior are different targets for the same
    molecule, and for a steric blocker the difference is whether it does
    anything at all."""
    assert EXON.region_for(start, end) is region


def test_without_an_exon_nothing_is_placed_and_the_report_says_so(target, rules):
    outcome = run(target, rules)

    assert all(candidate.region is Region.UNANNOTATED for candidate in outcome.candidates)
    assert any("no exon was supplied" in note for note in outcome.notes)


# --------------------------------------------------------------------------
# The confusion that inverts the therapy
# --------------------------------------------------------------------------


def test_a_gapmer_at_a_splice_site_is_blocked(target, rules):
    """RNase H cuts the transcript the gapmer binds. Aimed at a splice site to
    redirect splicing, it degrades the transcript that redirection was meant to
    rescue — correct in every other respect, and the opposite of the therapy."""
    outcome = run(target, rules, chemistry="gapmer-2MOE", exon=EXON)

    at_splice_sites = [
        candidate
        for candidate in outcome.candidates
        if candidate.region in {Region.ACCEPTOR_SITE, Region.DONOR_SITE}
    ]
    assert at_splice_sites
    for candidate in at_splice_sites:
        assert candidate.is_blocked
        assert any(
            flag.rule_id == "A_SPLICE_SITE_TARGET_NEEDS_A_STERIC_BLOCKER"
            for flag in candidate.flags
        )


def test_a_steric_blocker_at_the_same_site_is_not(target, rules):
    outcome = run(target, rules, chemistry="steric-PMO", exon=EXON)

    at_splice_sites = [
        candidate
        for candidate in outcome.candidates
        if candidate.region in {Region.ACCEPTOR_SITE, Region.DONOR_SITE}
    ]
    assert at_splice_sites
    assert not any(candidate.is_blocked for candidate in at_splice_sites)


def test_a_steric_blocker_in_the_middle_of_an_exon_is_cautioned(target, rules):
    """It binds and changes nothing, unless it happens to cover a silencer or an
    enhancer — and this package does not read those, which is why the flag is a
    caution rather than a block."""
    outcome = run(target, rules, chemistry="steric-PMO", exon=EXON)
    interior = [c for c in outcome.candidates if c.region is Region.EXON_INTERIOR]

    assert interior
    assert all(
        any(
            flag.rule_id == "A_STERIC_BLOCKER_AWAY_FROM_ANY_SPLICE_SIGNAL"
            for flag in candidate.flags
        )
        for candidate in interior
    )


# --------------------------------------------------------------------------
# Composition
# --------------------------------------------------------------------------


def test_a_run_of_four_guanines_blocks_a_window(rules):
    """G-quadruplex: the molecule folds on itself instead of on the target, and
    the window next door usually does not have the run."""
    sequence = "AT" * 20 + "CCCC" + "AT" * 20
    outcome = tile("G", "17", 1, sequence, rules, chemistry="steric-2MOE")
    containing = [c for c in outcome.candidates if "GGGG" in c.sequence]

    assert containing
    assert all(candidate.is_blocked for candidate in containing)


def test_cpg_dinucleotides_are_counted_and_flagged(rules):
    """Unmethylated CpG in a phosphorothioate backbone is a TLR9 agonist, and
    the fix is usually to slide the window rather than change the chemistry."""
    sequence = "ATCGATCGATCGATCGATCGATCGATCGATCG"
    outcome = tile("G", "17", 1, sequence, rules, chemistry="steric-2MOE")

    assert all(candidate.cpg_count > 1 for candidate in outcome.candidates)
    assert all(
        any(flag.rule_id == "TOO_MANY_CPG_MOTIFS" for flag in candidate.flags)
        for candidate in outcome.candidates
    )


def test_gc_content_outside_the_bounds_is_a_caution(rules):
    sequence = "AT" * 40
    outcome = tile("G", "17", 1, sequence, rules, chemistry="steric-2MOE")

    assert all(candidate.gc_fraction == 0.0 for candidate in outcome.candidates)
    assert all(
        any(flag.rule_id == "GC_CONTENT_OUT_OF_RANGE" for flag in candidate.flags)
        for candidate in outcome.candidates
    )


def test_self_complementarity_is_measured_crudely_and_labelled_as_crude():
    """It finds the longest stretch whose reverse complement also occurs, which
    is a hairpin's necessary condition and not its sufficient one."""
    assert longest_self_complement("AAAATTTT") >= 4
    assert longest_self_complement("AAAAAAAA") == 0


def test_a_window_with_an_unresolved_base_is_dropped_not_designed(rules):
    """An N in the reference is a base nobody knows. An oligonucleotide ordered
    against it is ordered against a guess."""
    sequence = "ACGT" * 5 + "N" + "ACGT" * 5
    outcome = tile("G", "17", 1, sequence, rules, chemistry="steric-2MOE")

    assert outcome.tiled > len(outcome.candidates)
    assert all("N" not in candidate.sequence for candidate in outcome.candidates)


# --------------------------------------------------------------------------
# The model that is not attached
# --------------------------------------------------------------------------


def test_accessibility_is_not_predicted_and_the_report_says_so(target, rules):
    outcome = run(target, rules)

    assert "no structure model is attached" in outcome.ranking
    assert NoStructureModel().accessibility(1, 20) is None


def test_a_model_that_is_attached_is_asked(target, rules):
    """The Protocol exists so that attaching RNAfold is a constructor argument
    rather than a rewrite."""

    class Flat:
        name = "flat"
        availability = "a test model that calls everything equally accessible"

        def accessibility(self, start: int, end: int) -> float:
            return 0.5

    outcome = tile("TARG", "17", 240, target, rules, chemistry="steric-PMO", model=Flat())

    assert "equally accessible" in outcome.ranking


def test_blocked_candidates_are_kept_and_separated(target, rules):
    """Not filtered out. A reader deciding between chemistries needs to see that
    a window was refused and why, and a shorter list would just look like fewer
    options."""
    outcome = run(target, rules, chemistry="gapmer-2MOE", exon=EXON)

    assert outcome.usable
    assert len(outcome.candidates) > len(outcome.usable)
    assert all(
        candidate.severity is Severity.BLOCKING
        for candidate in outcome.candidates
        if candidate.is_blocked
    )


# --------------------------------------------------------------------------
# What the operator sees
# --------------------------------------------------------------------------


def test_the_report_leads_with_how_many_were_tiled(capsys):
    assert (
        main(
            [
                "aso",
                "--gene", "TARG",
                "--at", "17:240-360",
                "--fasta", str(FASTA),
                "--chemistry", "steric-PMO",
                "--exon", "250-350",
                "--limit", "3",
            ]
        )
        == 0
    )
    out = capsys.readouterr().out

    assert "tiled" in out
    assert "windows" in out
    assert "no structure model is attached" in out
    assert "aso-v1@" in out


def test_a_span_in_the_wrong_shape_is_refused():
    assert (
        main(
            [
                "aso",
                "--gene", "TARG",
                "--at", "17:240",
                "--fasta", str(FASTA),
                "--chemistry", "steric-PMO",
            ]
        )
        == 1
    )
