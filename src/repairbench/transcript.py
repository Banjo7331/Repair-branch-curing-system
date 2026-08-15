"""Transcript arithmetic — the most consequential calculation in the package.

Whether a premature termination codon destroys the transcript or merely
truncates the protein decides which of two opposite therapies makes sense.

A transcript degraded by nonsense-mediated decay produces nothing. One allele's
worth of protein is missing, the mechanism is haploinsufficiency, and the
therapeutic question is how to supply more.

A transcript that escapes decay produces a shortened protein that is still
there, still folding, still binding its partners. If it poisons the complex it
belongs to, adding a working copy makes things worse, and the question inverts:
how to remove the mutant allele.

Same nonsense variant, same gene, opposite intervention. The difference is where
the stop codon sits relative to the last exon-exon junction — so that
calculation is written here, alone, with its boundaries named and tested.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import StrEnum

from repairbench.model import InvalidTranscriptError

#: Distance upstream of the final exon-exon junction within which a premature
#: stop escapes decay. The literature quotes 50–55 nt; the ClinGen SVI PVS1
#: decision tree uses 50, and so does this.
NMD_BOUNDARY_NT = 50

#: Distance from the start codon within which a premature stop may be bypassed
#: by reinitiation at a downstream methionine. Weaker than the junction rule and
#: treated as such — it lowers certainty rather than deciding the outcome.
REINITIATION_WINDOW_NT = 150

_VERSIONED_ACCESSION = re.compile(r"^[A-Z]{2}_\d+\.\d+$")


class NMDOutcome(StrEnum):
    """What the decay machinery is predicted to do."""

    PREDICTED = "predicted"
    """The transcript is degraded; this allele produces nothing."""

    ESCAPE = "escape"
    """The transcript survives; a truncated protein is produced."""

    NOT_APPLICABLE = "not_applicable"
    """The variant introduces no premature stop."""


@dataclass(frozen=True, slots=True)
class NMDPrediction:
    """An outcome with the reason it was reached and how sure the rule is.

    ``certain`` separates the junction rule, which is well established, from the
    reinitiation rule, which is not. A prediction that changes the therapeutic
    direction has to carry that distinction rather than flatten it.
    """

    outcome: NMDOutcome
    reason: str
    certain: bool = True


@dataclass(frozen=True, slots=True)
class Transcript:
    """The coding structure a variant is interpreted against.

    Exon lengths are *coding* nucleotides. Untranslated regions are excluded on
    purpose: every position here is a CDS offset, and mixing genomic with CDS
    coordinates is how NMD predictions silently go wrong.
    """

    accession: str
    gene: str
    coding_exon_lengths: tuple[int, ...]
    mane_select: bool = False

    def __post_init__(self) -> None:
        if not _VERSIONED_ACCESSION.match(self.accession):
            raise InvalidTranscriptError(
                f"transcript accession {self.accession!r} is unversioned — exon boundaries "
                "move between versions, and a prediction made against an unnamed structure "
                "cannot be reproduced"
            )
        if not self.coding_exon_lengths:
            raise InvalidTranscriptError(f"{self.accession} has no coding exons")
        for index, length in enumerate(self.coding_exon_lengths, start=1):
            if length <= 0:
                raise InvalidTranscriptError(
                    f"{self.accession} exon {index} has non-positive coding length {length}"
                )

    @property
    def coding_length(self) -> int:
        return sum(self.coding_exon_lengths)

    @property
    def cdna_kilobases(self) -> float:
        """Coding sequence in kilobases — what decides whether a gene fits in a vector."""
        return self.coding_length / 1000.0

    @property
    def exon_count(self) -> int:
        return len(self.coding_exon_lengths)

    @property
    def last_junction(self) -> int:
        """CDS offset of the final exon-exon junction.

        A single-exon transcript has none and reports 0.
        """
        if self.exon_count < 2:
            return 0
        return sum(self.coding_exon_lengths[:-1])

    def exon_at(self, cds_position: int) -> int:
        """Which exon contains a CDS position, 1-based."""
        self._require_in_range(cds_position)
        offset = 0
        for index, length in enumerate(self.coding_exon_lengths, start=1):
            offset += length
            if cds_position <= offset:
                return index
        raise InvalidTranscriptError(f"c.{cds_position} not located in {self.accession}")

    def exon_preserves_frame(self, exon_index: int) -> bool:
        """Does removing this exon leave the reading frame intact?

        The question behind exon skipping, and the reason Becker muscular
        dystrophy is milder than Duchenne: the same gene, deletions of similar
        size, and the ones that happen to be a multiple of three leave a protein
        that still works.
        """
        if not 1 <= exon_index <= self.exon_count:
            raise InvalidTranscriptError(f"{self.accession} has no exon {exon_index}")
        return self.coding_exon_lengths[exon_index - 1] % 3 == 0

    def predict_nmd(self, ptc_cds_position: int) -> NMDPrediction:
        """Apply the ClinGen SVI rules with the documented default boundaries."""
        return self.predict_nmd_with(ptc_cds_position)

    def predict_nmd_with(
        self,
        ptc_cds_position: int,
        *,
        boundary_nt: int = NMD_BOUNDARY_NT,
        reinitiation_window_nt: int = REINITIATION_WINDOW_NT,
    ) -> NMDPrediction:
        """Apply the rules with explicit boundaries.

        The distances are parameters because they are rule values rather than
        arithmetic: the junction boundary is quoted as 50–55 nt across the
        literature, and a laboratory that prefers 55 should be able to say so in
        the rule file instead of in a patch.
        """
        self._require_in_range(ptc_cds_position)

        if self.exon_count == 1:
            return NMDPrediction(
                NMDOutcome.ESCAPE,
                "single-exon transcript — no exon junction complex is deposited, "
                "so decay is not triggered",
            )

        junction = self.last_junction
        if ptc_cds_position > junction:
            return NMDPrediction(
                NMDOutcome.ESCAPE,
                f"stop at c.{ptc_cds_position} lies in the final exon, downstream of "
                f"the last junction at c.{junction}",
            )
        if ptc_cds_position > junction - boundary_nt:
            return NMDPrediction(
                NMDOutcome.ESCAPE,
                f"stop at c.{ptc_cds_position} is within {boundary_nt} nt of "
                f"the last junction at c.{junction}",
            )
        if ptc_cds_position <= reinitiation_window_nt:
            return NMDPrediction(
                NMDOutcome.ESCAPE,
                f"stop at c.{ptc_cds_position} sits within {reinitiation_window_nt} nt "
                "of the start codon; "
                "reinitiation at a downstream methionine may rescue translation",
                certain=False,
            )
        return NMDPrediction(
            NMDOutcome.PREDICTED,
            f"stop at c.{ptc_cds_position} lies more than {boundary_nt} nt upstream of "
            f"the last junction at c.{junction}",
        )

    def _require_in_range(self, cds_position: int) -> None:
        if not 1 <= cds_position <= self.coding_length:
            raise InvalidTranscriptError(
                f"CDS position {cds_position} is outside {self.accession} "
                f"(coding length {self.coding_length})"
            )
