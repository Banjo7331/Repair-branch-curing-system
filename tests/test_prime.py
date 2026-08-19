"""pegRNA geometry, which is four reverse complements and one inequality.

The inequality is the important one: the nick must fall *before* the edit on the
strand it cuts. A protospacer that fails it looks perfect — good PAM, right
distance, clean sequence — and its template writes over sequence the polymerase
never reaches. Half the tests here exist to keep that constraint honest.

The rest check the reverse complements. The primer binding site, the template
and the minus strand are each one, and a wrong one produces a pegRNA of exactly
the right length that installs nothing. So the template is verified against the
genome from the outside: reverse-complement it back and it must equal the
edited sequence downstream of the nick, base for base.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from repairbench.annotation.fasta import IndexedFasta, InMemorySequences
from repairbench.cli import main
from repairbench.design.editors import DesignError
from repairbench.design.flags import Severity, load_flag_rules, sort_weight
from repairbench.design.prime import (
    EditRequest,
    design_pegrnas,
    melting_temperature,
)
from repairbench.design.sequence import reverse_complement

DATA = Path(__file__).parent / "data" / "design"
FASTA = DATA / "target.fa"
RULES = Path(__file__).parents[1] / "rules" / "prime-v1.yaml"

#: The site the fixture was built around: a transversion, which base editing
#: refuses outright, with an NGG placed so the nick lands 14 nt upstream.
SITE = 500


@pytest.fixture
def rules():
    return load_flag_rules(RULES)


@pytest.fixture
def genome():
    with IndexedFasta(FASTA) as fasta:
        yield fasta


def outcome(genome, rules, patient: str = "A", wild_type: str = "C", position: int = SITE):
    return design_pegrnas(
        EditRequest("TARG", "17", position, patient, wild_type), genome, rules
    )


# --------------------------------------------------------------------------
# What prime editing is for
# --------------------------------------------------------------------------


def test_a_transversion_is_designed_rather_than_refused(genome, rules):
    """The whole reason this module exists beside the base editor. A→C is a
    transversion; no deaminase makes one, and a prime editor does not care."""
    result = outcome(genome, rules, "A", "C")

    assert result.has_candidates
    assert "substitution" in result.edit


@pytest.mark.parametrize(
    ("patient", "wild_type", "kind"),
    [("A", "C", "substitution"), ("A", "AGG", "insertion"), ("ACT", "A", "deletion")],
)
def test_all_three_edit_kinds_are_designed_the_same_way(
    genome, rules, patient, wild_type, kind
):
    """A substitution is a one-for-one template, an insertion a longer one, a
    deletion a shorter one. The arithmetic does not branch, and these assert
    that it does not need to."""
    result = outcome(genome, rules, patient, wild_type)

    assert kind in result.edit
    assert result.has_candidates


def test_an_allele_written_as_an_empty_string_is_refused():
    """VCF writes an insertion with an anchor base in both alleles. An empty
    allele makes the coordinate ambiguous, and this module will not guess which
    convention the caller meant."""
    with pytest.raises(DesignError, match="the way a VCF does"):
        EditRequest("G", "17", 500, "", "GG")


def test_correcting_an_allele_to_itself_is_refused():
    with pytest.raises(DesignError, match="nothing to correct"):
        EditRequest("G", "17", 500, "A", "A")


# --------------------------------------------------------------------------
# The inequality
# --------------------------------------------------------------------------


def test_every_nick_falls_before_the_edit(genome, rules):
    """The constraint that eliminates about half the PAMs near any variant.
    Reverse transcription runs from the nick forwards, so an edit upstream of it
    is never written."""
    result = outcome(genome, rules)

    for candidate in result.candidates:
        if candidate.strand == "+":
            assert candidate.nick_position < SITE
        else:
            assert candidate.nick_position > SITE
        assert 1 <= candidate.nick_to_edit_nt <= 30


def test_a_site_with_no_upstream_pam_says_which_constraint_failed(rules):
    """Not "no candidates" — the refusal names the geometry, because the fix is
    a different nuclease rather than a longer template."""
    flat = InMemorySequences({"17": "CT" * 300})
    result = design_pegrnas(EditRequest("G", "17", 300, "A", "C"), flat, rules)

    assert not result.has_candidates
    assert "must fall *before* the edit" in result.refusals[0]


# --------------------------------------------------------------------------
# The reverse complements
# --------------------------------------------------------------------------


def test_the_template_reads_back_as_the_edited_genome(genome, rules):
    """The strongest check available without a laboratory: reverse-complement
    the template and it must equal the genome as it should read, base for base,
    starting one base after the nick."""
    result = outcome(genome, rules, "A", "C")
    candidate = next(c for c in result.candidates if c.strand == "+")

    written = reverse_complement(candidate.rtt)
    start = candidate.nick_position + 1
    genomic = genome.fetch("17", start, start + len(written) - 1)
    expected = (
        genomic[: SITE - start] + "C" + genomic[SITE - start + 1 :]
    )  # the reference with the edit installed

    assert written == expected
    assert written[SITE - start] == "C"


def test_the_template_carries_the_inserted_bases(genome, rules):
    """An insertion is longer than the sequence it replaces, so the template is
    longer than the span it covers — which is the only way the arithmetic can
    tell an insertion from a substitution."""
    substitution = outcome(genome, rules, "A", "C").candidates[0]
    insertion = outcome(genome, rules, "A", "AGG").candidates[0]

    assert len(insertion.rtt) == len(substitution.rtt) + 2
    assert "AGG" in reverse_complement(insertion.rtt)


def test_the_primer_binding_site_anneals_to_what_is_left_of_the_nick(genome, rules):
    """The PBS pairs with the nicked strand's 3' end — the bases immediately
    upstream of the nick, not downstream. Getting the side wrong produces a
    pegRNA that anneals to the template it is about to write."""
    result = outcome(genome, rules)
    candidate = next(c for c in result.candidates if c.strand == "+")

    upstream = genome.fetch(
        "17", candidate.nick_position - len(candidate.pbs) + 1, candidate.nick_position
    )
    assert candidate.pbs == reverse_complement(upstream)


def test_the_extension_is_template_then_primer(genome, rules):
    """Written in the order it is synthesised in, which is the reverse of the
    order it is used in. A pegRNA ordered the other way round is a real and
    silent failure mode."""
    candidate = outcome(genome, rules).candidates[0]

    assert candidate.extension == candidate.rtt + candidate.pbs
    assert candidate.extension.endswith(candidate.pbs)


def test_the_homology_arm_is_scanned_across_the_range_the_rule_file_allows(genome, rules):
    """Both design parameters are scanned, not just the primer binding site.

    This asserted a single length until a published pegRNA turned out to use a
    longer template than the minimum — which meant the molecule somebody made
    was not expressible in this module's design space at all."""
    minimum = int(rules.threshold("rtt_homology_min_nt", 5))
    maximum = int(rules.threshold("rtt_homology_max_nt", 20))
    arms = {candidate.homology_arm_nt for candidate in outcome(genome, rules).candidates}

    assert arms, "no candidates to inspect"
    assert min(arms) == minimum
    assert len(arms) > 1, "the template length is a parameter, not a constant"
    assert all(minimum <= arm <= maximum for arm in arms)
    for candidate in outcome(genome, rules).candidates:
        assert len(candidate.rtt) >= candidate.homology_arm_nt


# --------------------------------------------------------------------------
# PE3 and PE3b
# --------------------------------------------------------------------------


def test_a_pe3b_nick_is_found_and_named(genome, rules):
    """The fixture places a protospacer on the opposite strand that overlaps the
    edit, so its spacer matches the corrected allele and not the original."""
    candidate = outcome(genome, rules).candidates[0]

    assert candidate.pe3b
    guide = candidate.pe3b[0]
    assert guide.edit_dependent
    assert guide.strand != candidate.strand


def test_pe3b_is_not_held_to_the_pe3_distance_window(genome, rules):
    """Mechanical rather than empirical: a PE3b spacer has to overlap the edit,
    so its nick lands wherever the edit is — usually much closer than the 40 nt
    PE3 floor. Holding it to that window would rule out the configuration with
    the fewest indels."""
    candidate = outcome(genome, rules).candidates[0]

    assert candidate.pe3b[0].distance_nt < int(rules.threshold("pe3_nick_min_nt", 40))


def test_pe3b_sorts_ahead_of_a_plain_second_nick(genome, rules):
    candidate = outcome(genome, rules).candidates[0]

    assert candidate.nicking_guides[0].edit_dependent


def test_no_second_nick_anywhere_is_a_caution_not_a_silence(rules):
    """PE2 alone works. It is also several-fold less efficient, so the absence
    is a fact to plan around and is stated as one."""
    sequence = "CTG" * 60 + "AGG" + "CTG" * 60
    flat = InMemorySequences({"17": sequence})
    # The edit sits 14 nt past the nick of the single PAM in this sequence.
    result = design_pegrnas(EditRequest("G", "17", 195, "A", "C"), flat, rules)

    if result.has_candidates:
        flagged = {
            flag.rule_id for candidate in result.candidates for flag in candidate.flags
        }
        assert "NO_SECOND_NICK_IN_RANGE" in flagged


# --------------------------------------------------------------------------
# The rule file
# --------------------------------------------------------------------------


def test_a_surviving_pam_is_flagged(genome, rules):
    """With the PAM intact the editor re-engages what it has just corrected, and
    the second pass installs an indel instead of the edit."""
    candidate = outcome(genome, rules).candidates[0]

    assert not candidate.pam_disrupted_by_edit
    assert any(flag.rule_id == "THE_PAM_SURVIVES_THE_EDIT" for flag in candidate.flags)


def test_a_poly_t_run_blocks_a_candidate_outright(rules, tmp_path: Path):
    """Four Ts is a Pol III terminator, so the molecule is truncated at the run
    rather than made as designed. That is not a preference."""
    blocked = load_flag_rules(RULES)
    sequence = "TTTT" + "CTG" * 100
    flat = InMemorySequences({"17": sequence})
    result = design_pegrnas(EditRequest("G", "17", 60, "A", "C"), flat, blocked)

    for candidate in result.candidates:
        if "TTTT" in candidate.spacer or "TTTT" in candidate.extension:
            assert candidate.is_blocked


def test_the_severity_of_a_candidate_is_its_worst_flag(genome, rules):
    """Reasons a design will not work do not average out against reasons it
    might."""
    for candidate in outcome(genome, rules).candidates:
        if candidate.flags:
            assert candidate.severity is min(
                (flag.severity for flag in candidate.flags), key=lambda level: level.rank
            )


def test_unflagged_candidates_sort_ahead_of_flagged_ones():
    """The one ordering mistake that inverts a whole report, and it is easy to
    make: the severity ranks count down from the worst, and a candidate with no
    flags has no severity at all."""
    assert sort_weight(None) < sort_weight(Severity.NOTE) < sort_weight(Severity.BLOCKING)


def test_the_ruleset_is_pinned(genome, rules):
    assert outcome(genome, rules).ruleset_pin.startswith("prime-v1@")


# --------------------------------------------------------------------------
# Arithmetic that is only approximate, and says so
# --------------------------------------------------------------------------


def test_the_melting_temperature_is_the_wallace_rule():
    """Two degrees per A or T, four per G or C. Named as an approximation in the
    docstring because the threshold it is compared against was fitted with one."""
    assert melting_temperature("AAAA") == 8.0
    assert melting_temperature("GCGC") == 16.0


def test_nothing_is_ranked_by_predicted_efficiency(genome, rules):
    assert "no efficiency model is attached" in outcome(genome, rules).ranking
    assert "PRIDICT" in outcome(genome, rules).ranking


# --------------------------------------------------------------------------
# The reference, which is normally the wild type and sometimes is not
# --------------------------------------------------------------------------


def test_the_ordinary_case_says_nothing(genome, rules):
    """The fixture's reference carries the wild-type base, which is the normal
    arrangement and not worth a line of output."""
    assert outcome(genome, rules, "A", "C").notes == ()


def test_a_reference_carrying_the_patients_allele_is_reported_not_refused(genome, rules):
    """A pathogenic allele can be the reference allele. The design is unaffected;
    what changes is that "reference" and "wild type" have come apart, and
    nothing downstream should treat them as one word."""
    result = outcome(genome, rules, "C", "A")

    assert any("carries the patient's allele" in note for note in result.notes)
    assert result.has_candidates


def test_a_reference_matching_neither_allele_is_called_out(genome, rules):
    result = outcome(genome, rules, "G", "T")

    assert any("matches neither allele" in note for note in result.notes)


# --------------------------------------------------------------------------
# What the operator sees
# --------------------------------------------------------------------------


def test_the_report_says_how_many_were_not_shown(capsys):
    assert (
        main(
            [
                "pegrna",
                "--gene", "TARG",
                "--at", f"17:{SITE}",
                "--patient", "A",
                "--wild-type", "C",
                "--fasta", str(FASTA),
                "--limit", "2",
            ]
        )
        == 0
    )
    out = capsys.readouterr().out

    assert "more, not shown" in out
    assert "prime-v1@" in out
    assert "PE3b" in out


def test_the_command_exits_non_zero_when_nothing_is_usable(capsys, tmp_path: Path):
    """A pipeline step that designed nothing should not look successful."""
    flat = tmp_path / "flat.fa"
    sequence = "CT" * 300
    flat.write_text(f">17 flat\n{sequence}\n")
    (tmp_path / "flat.fa.fai").write_text(f"17\t{len(sequence)}\t{len('>17 flat')+1}\t600\t601\n")

    assert (
        main(
            [
                "pegrna",
                "--gene", "TARG",
                "--at", "17:300",
                "--patient", "A",
                "--wild-type", "C",
                "--fasta", str(flat),
            ]
        )
        == 2
    )
    assert "not designed" in capsys.readouterr().out
