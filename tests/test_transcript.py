"""The NMD boundary, tested at the boundary.

Every case here is one nucleotide's difference from its neighbour, because that
is where this calculation is wrong if it is wrong — and being wrong here inverts
the therapy rather than degrading it.
"""

from __future__ import annotations

import pytest

from repairbench.model import InvalidTranscriptError
from repairbench.transcript import NMDOutcome, Transcript

# Four exons of 300 coding nucleotides. The last junction sits at c.900, so the
# escape window opens at c.851.
FOUR_EXONS = Transcript("NM_000001.1", "TESTG", (300, 300, 300, 300))


def test_a_stop_well_upstream_of_the_last_junction_triggers_decay():
    prediction = FOUR_EXONS.predict_nmd(500)

    assert prediction.outcome is NMDOutcome.PREDICTED
    assert prediction.certain
    assert "c.900" in prediction.reason


def test_the_last_nucleotide_inside_the_decay_window():
    """c.850 is 50 nt from the junction and still degraded."""
    assert FOUR_EXONS.predict_nmd(850).outcome is NMDOutcome.PREDICTED


def test_one_nucleotide_later_escapes():
    """c.851 is within 50 nt of the junction — the therapy inverts here."""
    prediction = FOUR_EXONS.predict_nmd(851)

    assert prediction.outcome is NMDOutcome.ESCAPE
    assert prediction.certain


def test_a_stop_in_the_final_exon_escapes():
    prediction = FOUR_EXONS.predict_nmd(1000)

    assert prediction.outcome is NMDOutcome.ESCAPE
    assert "final exon" in prediction.reason


def test_a_single_exon_transcript_has_no_junction_to_trigger_decay():
    single = Transcript("NM_000002.1", "TESTG", (900,))

    prediction = single.predict_nmd(100)

    assert prediction.outcome is NMDOutcome.ESCAPE
    assert "single-exon" in prediction.reason


def test_a_start_proximal_stop_escapes_but_not_certainly():
    """Reinitiation is a weaker rule, and the prediction says so."""
    prediction = FOUR_EXONS.predict_nmd(120)

    assert prediction.outcome is NMDOutcome.ESCAPE
    assert not prediction.certain


def test_the_boundary_is_a_parameter_not_a_constant():
    """A laboratory preferring 55 nt should not need a patch."""
    at_fifty = FOUR_EXONS.predict_nmd_with(851, boundary_nt=50)
    at_forty = FOUR_EXONS.predict_nmd_with(851, boundary_nt=40)

    assert at_fifty.outcome is NMDOutcome.ESCAPE
    assert at_forty.outcome is NMDOutcome.PREDICTED


def test_exon_lookup_is_one_based_and_inclusive():
    assert FOUR_EXONS.exon_at(1) == 1
    assert FOUR_EXONS.exon_at(300) == 1
    assert FOUR_EXONS.exon_at(301) == 2
    assert FOUR_EXONS.exon_at(1200) == 4


def test_frame_preservation_is_exon_length_modulo_three():
    """The arithmetic behind exon skipping, and behind Becker versus Duchenne."""
    mixed = Transcript("NM_000003.1", "TESTG", (300, 301, 299))

    assert mixed.exon_preserves_frame(1)
    assert not mixed.exon_preserves_frame(2)
    assert not mixed.exon_preserves_frame(3)


def test_an_unversioned_accession_is_refused():
    """Exon boundaries move between versions; an unnamed structure is unreproducible."""
    with pytest.raises(InvalidTranscriptError, match="unversioned"):
        Transcript("NM_000004", "TESTG", (300, 300))


def test_a_position_outside_the_coding_sequence_is_refused():
    with pytest.raises(InvalidTranscriptError, match="outside"):
        FOUR_EXONS.predict_nmd(5000)


def test_coding_length_drives_the_payload_question():
    assert FOUR_EXONS.cdna_kilobases == pytest.approx(1.2)
