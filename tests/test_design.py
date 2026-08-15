"""Placing a base editor over a substitution.

Almost every test here is about a way this module could be wrong quietly. A
protospacer placed on the wrong strand is still twenty plausible bases; a window
numbered from the wrong end still puts the target inside it sometimes; a design
run against the reference instead of the patient produces guides for a genome
nobody has. None of those raises.

The fixture is synthetic and says so — ``tests/data/design/target.fa`` is a
repeating trimer with PAMs placed by hand at known distances, because what is
under test is the arithmetic, and the arithmetic does not care whether the
coordinates are a real locus.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from repairbench.annotation.fasta import IndexedFasta, InMemorySequences
from repairbench.annotation.gff import parse_gff3
from repairbench.annotation.store import TranscriptStore
from repairbench.cli import main
from repairbench.design.candidate import EditCandidate
from repairbench.design.designer import CorrectionRequest, design
from repairbench.design.editors import Conversion, DesignError, load_editors
from repairbench.design.efficiency import EfficiencyScore, NoModelAttached
from repairbench.design.sequence import matches_pam, reverse_complement
from repairbench.model import Zygosity

DATA = Path(__file__).parent / "data" / "design"
FASTA = DATA / "target.fa"
GFF = DATA / "target.gff3"
CATALOGUE = Path(__file__).parents[1] / "rules" / "editors-v1.yaml"

#: The two sites the fixture was built around. At 301 the patient has an A that
#: should read G, correctable on the plus strand. At 401 the patient has a T
#: that should read C — the same correction seen from the other strand, which is
#: the case a single-strand scanner reports as "no candidates".
PLUS_SITE, MINUS_SITE = 301, 401


@pytest.fixture
def catalogue():
    return load_editors(CATALOGUE)


@pytest.fixture
def genome():
    with IndexedFasta(FASTA) as fasta:
        yield fasta


def outcome_at(genome, catalogue, position: int, patient: str, wild_type: str, **kwargs):
    return design(
        CorrectionRequest(
            gene="TARG",
            chromosome="17",
            position=position,
            patient_base=patient,
            wild_type_base=wild_type,
            **kwargs,
        ),
        genome,
        catalogue,
    )


# --------------------------------------------------------------------------
# The catalogue
# --------------------------------------------------------------------------


def test_the_catalogue_is_pinned_by_its_contents(catalogue):
    """A design made before somebody widened an editing window and one made
    after are different designs, and the pin is what says so."""
    assert catalogue.pin.startswith("editors-v1@")
    assert len(catalogue.pin.split("@")[1]) == 12


def test_a_deaminase_that_makes_a_transversion_is_refused(tmp_path: Path):
    """No enzyme in this class does it, so a catalogue claiming one is a typo or
    a misunderstanding, and either way the run should stop."""
    path = tmp_path / "editors.yaml"
    path.write_text(
        "version: bad\neditors:\n"
        "  - {id: X, conversion: 'A>C', pam: NGG, protospacer_length: 20, "
        "window: {start: 4, end: 8}}\n"
    )

    with pytest.raises(DesignError, match="base editing makes"):
        load_editors(path)


def test_a_window_outside_the_protospacer_is_refused(tmp_path: Path):
    path = tmp_path / "editors.yaml"
    path.write_text(
        "version: bad\neditors:\n"
        "  - {id: X, conversion: 'A>G', pam: NGG, protospacer_length: 20, "
        "window: {start: 4, end: 25}}\n"
    )

    with pytest.raises(DesignError, match="does not fit"):
        load_editors(path)


# --------------------------------------------------------------------------
# Which correction is even a base edit
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("patient", "wild_type", "conversion"),
    [("A", "G", Conversion.A_TO_G), ("C", "T", Conversion.C_TO_T)],
)
def test_the_two_conversions_a_deaminase_makes(patient, wild_type, conversion):
    assert Conversion.between(patient, wild_type) is conversion


def test_a_transversion_is_refused_with_the_reason(genome, catalogue):
    """The commonest thing a user will ask for that this module cannot do. An
    empty candidate list would read as "none found near here", which is a
    different and much more encouraging claim."""
    outcome = outcome_at(genome, catalogue, PLUS_SITE, "A", "C")

    assert not outcome.has_candidates
    assert "transversion" in outcome.refusals[0]
    assert "prime editing" in outcome.refusals[0]


def test_an_indel_is_refused_before_anything_is_read():
    with pytest.raises(DesignError, match="prime editing's problem"):
        CorrectionRequest("G", "17", 301, patient_base="AT", wild_type_base="A")


def test_correcting_a_base_to_itself_is_refused():
    with pytest.raises(DesignError, match="nothing to correct"):
        CorrectionRequest("G", "17", 301, patient_base="A", wild_type_base="A")


# --------------------------------------------------------------------------
# The geometry
# --------------------------------------------------------------------------


def test_the_target_always_lands_inside_the_stated_window(genome, catalogue):
    """The one invariant this module exists to keep. Everything else is
    presentation."""
    outcome = outcome_at(genome, catalogue, PLUS_SITE, "A", "G")

    assert outcome.has_candidates
    for candidate in outcome.candidates:
        assert candidate.target_position_in_protospacer in candidate.editor.window


def test_the_protospacer_really_sits_where_the_candidate_says_it_does(genome, catalogue):
    """Reads the reference back and checks the arithmetic from the outside:
    twenty bases at the stated span, the PAM immediately 3' of them, and the
    target where the candidate claims it is."""
    outcome = outcome_at(genome, catalogue, PLUS_SITE, "A", "G")
    candidate = next(c for c in outcome.candidates if c.strand == "+")

    start, end = candidate.span
    assert end - start + 1 == candidate.editor.protospacer_length
    assert candidate.pam_span[0] == end + 1
    assert matches_pam(candidate.pam, candidate.editor.pam)
    assert (
        start + candidate.target_position_in_protospacer - 1 == candidate.target_genomic_position
    )


def test_a_minus_strand_correction_is_found_at_all(genome, catalogue):
    """T→C is not a cytosine editor's job — it is an adenine editor's job on the
    other strand, where the same base pair reads A and needs to read G. A
    designer that scans one strand reports nothing here, and reporting nothing
    is a wrong answer, not a missing one."""
    outcome = outcome_at(genome, catalogue, MINUS_SITE, "T", "C")

    assert outcome.has_candidates
    assert all(candidate.strand == "-" for candidate in outcome.candidates)
    assert all(candidate.conversion is Conversion.A_TO_G for candidate in outcome.candidates)


def test_the_minus_strand_protospacer_is_the_reverse_complement_of_the_reference(
    genome, catalogue
):
    outcome = outcome_at(genome, catalogue, MINUS_SITE, "T", "C")
    candidate = outcome.candidates[0]

    start, end = candidate.span
    plus_strand = genome.fetch("17", start, end)
    # The patient's base differs from the reference at exactly one position, so
    # compare everything except the target.
    expected = reverse_complement(plus_strand)
    index = candidate.target_position_in_protospacer - 1
    assert candidate.protospacer[:index] == expected[:index]
    assert candidate.protospacer[index + 1 :] == expected[index + 1 :]
    assert candidate.protospacer[index] == "A"


def test_the_design_is_against_the_patients_sequence_not_the_reference(genome, catalogue):
    """The reference does not have the base being corrected — that is what makes
    it the base being corrected. Every protospacer must carry the patient's
    allele, or the guide is for somebody else's genome."""
    outcome = outcome_at(genome, catalogue, PLUS_SITE, "A", "G")

    for candidate in outcome.candidates:
        assert candidate.protospacer[candidate.target_position_in_protospacer - 1] == "A"
    assert genome.fetch("17", PLUS_SITE, PLUS_SITE) == "G"


def test_a_reference_that_disagrees_with_the_wild_type_base_is_reported(genome, catalogue):
    """Not an error — a reference allele is not always the wild-type one — but it
    is also what a wrong coordinate or a wrong assembly looks like."""
    outcome = outcome_at(genome, catalogue, PLUS_SITE, "C", "T")

    assert any("the reference reads G" in note for note in outcome.notes)


def test_no_pam_at_a_usable_distance_says_so(catalogue):
    """A window with no G in it at all: nothing to make an NGG from, and the
    refusal names the relaxed-PAM option rather than leaving a bare empty list."""
    flat = InMemorySequences({"17": "CT" * 200})
    outcome = design(
        CorrectionRequest("TARG", "17", 100, patient_base="A", wild_type_base="G"), flat, catalogue
    )

    assert not outcome.has_candidates
    assert "no PAM in the catalogue" in outcome.refusals[0]


# --------------------------------------------------------------------------
# Bystanders
# --------------------------------------------------------------------------


def test_every_other_editable_base_in_the_window_is_listed(genome, catalogue):
    """The fixture puts a second A at 17:303, two positions from the target. It
    is inside the window for every editor that reaches the target, and it must
    appear as a coordinate rather than as a count."""
    outcome = outcome_at(genome, catalogue, PLUS_SITE, "A", "G")
    candidate = next(c for c in outcome.candidates if c.editor.id == "ABE8e-SpCas9")

    assert [bystander.genomic_position for bystander in candidate.bystanders] == [303]
    assert candidate.bystanders[0].becomes == "G"
    assert not candidate.is_clean


def test_a_base_outside_the_window_is_not_a_bystander(genome, catalogue):
    """ABE7.10's window stops at position 7 and ABE8e's at 8. A base at position
    8 is a bystander for one and not for the other, and collapsing the two would
    make the window in the catalogue decorative."""
    outcome = outcome_at(genome, catalogue, PLUS_SITE, "A", "G")

    for candidate in outcome.candidates:
        for bystander in candidate.bystanders:
            assert bystander.position_in_protospacer in candidate.editor.window


def test_the_target_is_never_listed_as_its_own_bystander(genome, catalogue):
    outcome = outcome_at(genome, catalogue, PLUS_SITE, "A", "G")

    for candidate in outcome.candidates:
        positions = [b.genomic_position for b in candidate.bystanders]
        assert candidate.target_genomic_position not in positions


def test_a_bystander_in_coding_sequence_is_marked_when_an_annotation_is_supplied(
    genome, catalogue
):
    """And only then. Without a transcript the field is ``None`` — nobody looked
    — rather than False, which would read as "checked, not coding"."""
    store = TranscriptStore(parse_gff3(GFF))
    request = CorrectionRequest("TARG", "17", PLUS_SITE, patient_base="A", wild_type_base="G")

    annotated = design(
        request,
        genome,
        catalogue,
        coding=lambda position: store.locate("17", position).in_coding_sequence,
    )
    bare = design(request, genome, catalogue)

    with_bystanders = next(c for c in annotated.candidates if c.bystanders)
    assert with_bystanders.bystanders[0].in_coding_sequence is True
    assert any("is in coding sequence" in warning for warning in with_bystanders.warnings)
    assert next(c for c in bare.candidates if c.bystanders).bystanders[0].in_coding_sequence is None


def test_an_unaffected_copy_earns_a_warning_about_the_healthy_allele(genome, catalogue):
    """A heterozygote's two alleles differ by one base, and it sits in the
    PAM-distal half of the guide where it discriminates poorly. The intended
    edit is harmless on the healthy allele — it is already correct — but every
    bystander applies to it too."""
    outcome = outcome_at(
        genome, catalogue, PLUS_SITE, "A", "G", zygosity=Zygosity.HETEROZYGOUS
    )

    assert any("bind both alleles" in note for note in outcome.notes)


def test_the_heterozygote_note_is_made_once_and_not_per_candidate(genome, catalogue):
    """It is a fact about the patient, not about a placement. Repeated under six
    candidates it is a caution nobody finishes reading."""
    outcome = outcome_at(genome, catalogue, PLUS_SITE, "A", "G", zygosity=Zygosity.HETEROZYGOUS)

    assert len([note for note in outcome.notes if "bind both alleles" in note]) == 1
    assert not any(
        "bind both alleles" in warning
        for candidate in outcome.candidates
        for warning in candidate.warnings
    )


def test_an_unknown_zygosity_raises_no_such_note(genome, catalogue):
    outcome = outcome_at(genome, catalogue, PLUS_SITE, "A", "G")

    assert not any("bind both alleles" in note for note in outcome.notes)


# --------------------------------------------------------------------------
# The model that is not attached
# --------------------------------------------------------------------------


def test_nothing_is_ranked_by_efficiency_and_the_report_says_so(genome, catalogue):
    outcome = outcome_at(genome, catalogue, PLUS_SITE, "A", "G")

    assert "no efficiency model is attached" in outcome.ranking
    assert NoModelAttached().score(outcome.candidates[0]) is None


def test_clean_candidates_come_first_because_that_is_a_stated_criterion(genome, catalogue):
    """Not a quality ranking — a count of what the editor can also reach. The
    ordering is declared rather than implied, which is the difference between a
    criterion and a heuristic."""
    outcome = outcome_at(genome, catalogue, MINUS_SITE, "T", "C")
    counts = [len(candidate.bystanders) for candidate in outcome.candidates]

    assert counts == sorted(counts)


def test_a_model_that_is_attached_decides_the_order(genome, catalogue):
    """The point of the Protocol: attaching a real model is a constructor
    argument, not a rewrite. This one scores by span so the test can assert the
    order changed."""

    class Backwards:
        name = "backwards"
        availability = "a test model that prefers whatever comes last"

        def score(self, candidate: EditCandidate) -> EfficiencyScore:
            return EfficiencyScore(value=candidate.span[0] / 10_000, model=self.name)

    request = CorrectionRequest("TARG", "17", PLUS_SITE, patient_base="A", wild_type_base="G")
    scored = design(request, genome, catalogue, model=Backwards())

    spans = [candidate.span[0] for candidate in scored.candidates]
    assert spans == sorted(spans, reverse=True)
    assert "prefers whatever comes last" in scored.ranking


def test_a_score_outside_zero_to_one_is_refused():
    with pytest.raises(ValueError, match="not a fraction"):
        EfficiencyScore(value=42.0, model="overconfident")


# --------------------------------------------------------------------------
# What the operator sees
# --------------------------------------------------------------------------


def test_the_command_exits_non_zero_when_it_designed_nothing(capsys):
    """A pipeline step that produces no candidate should not look successful.
    Two is used rather than one so a refusal is distinguishable from a crash."""
    assert (
        main(
            [
                "design",
                "--gene", "TARG",
                "--at", f"17:{PLUS_SITE}",
                "--patient", "A",
                "--wild-type", "C",
                "--fasta", str(FASTA),
            ]
        )
        == 2
    )
    assert "transversion" in capsys.readouterr().out


def test_the_report_prints_the_catalogue_it_was_made_under(capsys):
    main(
        [
            "design",
            "--gene", "TARG",
            "--at", f"17:{PLUS_SITE}",
            "--patient", "A",
            "--wild-type", "G",
            "--fasta", str(FASTA),
        ]
    )
    out = capsys.readouterr().out

    assert "editors-v1@" in out
    assert "no efficiency model is attached" in out
    assert "not a therapy" in out


def test_a_coordinate_in_the_wrong_shape_is_refused(capsys):
    assert main(
        [
            "design",
            "--gene", "TARG",
            "--at", "17-301",
            "--patient", "A",
            "--wild-type", "G",
            "--fasta", str(FASTA),
        ]
    ) == 1
