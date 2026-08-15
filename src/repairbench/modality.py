"""What could be done about a mechanism — and, more importantly, what must not.

M5 answers *why this variant causes disease*. M6 answers *what class of
intervention that mechanism admits*. It stops there. Which molecule, at what
dose, in which tissue, is M7 and beyond.

The asymmetry in this module is deliberate and is its main safety property: a
contraindication outranks any number of indications. Reasons to try something
accumulate; one reason not to ends the matter. That is the correct ordering for
a system whose worst available mistake is proposing gene addition for a
dominant-negative variant — where supplying a normal copy leaves every poisoning
subunit exactly where it was.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

from repairbench.model import Mechanism, Strength


class Modality(StrEnum):
    """Classes of intervention, at the level of what they do to the biology."""

    GENE_ADDITION = "gene_addition"
    """Deliver a working copy. The old gene stays where it is."""

    TRUNCATED_CONSTRUCT = "truncated_construct"
    """Deliver a shortened but functional version, for genes too large to carry.
    Micro-dystrophin is the worked example, and it is protein engineering rather
    than gene replacement."""

    WILD_TYPE_UPREGULATION = "wild_type_upregulation"
    """Raise output from the intact allele — for instance by blocking a
    non-productive splice event with an antisense oligonucleotide, so more
    transcript reaches translation. Applicable where the gene is too large to
    deliver and one allele still works."""

    ALLELE_SPECIFIC_SILENCING = "allele_specific_silencing"
    """Remove the mutant transcript while leaving the normal one. The answer when
    the product is doing harm rather than being absent."""

    EXON_SKIPPING = "exon_skipping"
    """Splice out the affected exon, when what is left is still in frame and
    still works."""

    SPLICE_CORRECTION = "splice_correction"
    """Restore normal splicing disrupted by the variant."""

    SILENCED_ALLELE_REACTIVATION = "silenced_allele_reactivation"
    """Wake an intact copy that is switched off — by genomic imprinting, or by X
    inactivation in a heterozygous female. Different biology, same opportunity:
    the sequence is already there and only needs turning on."""

    BASE_EDITING = "base_editing"
    """Chemically convert one base in place. Narrow in what it can change, and
    it changes every eligible base in its window, not only the intended one."""

    PRIME_EDITING = "prime_editing"
    """Write a small arbitrary edit in place. The closest thing to literally
    correcting the sequence, and the first such therapy reached a patient in
    2025."""

    READ_THROUGH = "read_through"
    """Persuade the ribosome past a premature stop. Mechanistically appealing,
    clinically contested."""


def all_modalities() -> tuple[Modality, ...]:
    return tuple(Modality)


class Stance(StrEnum):
    """What a rule says about a modality."""

    INDICATES = "indicates"
    CONTRAINDICATES = "contraindicates"


class Verdict(StrEnum):
    """Where a modality ended up."""

    INDICATED = "indicated"
    """Coherent with the mechanism, nothing rules it out. Not a recommendation."""

    NOT_INDICATED = "not_indicated"
    """Nothing argues for it here. Silence, not prohibition."""

    CONTRAINDICATED = "contraindicated"
    """Something argues actively against it. This is the output that matters."""

    BLOCKED_BY_UNRESOLVED_MECHANISM = "blocked_by_unresolved_mechanism"
    """No mechanism was established, so no modality can be assessed at all."""

    @property
    def is_available(self) -> bool:
        return self is Verdict.INDICATED


@dataclass(frozen=True, slots=True)
class ModalityEvidence:
    """One rule that fired about one modality."""

    rule_id: str
    modality: Modality
    stance: Stance
    strength: Strength
    because: str
    citation: str = ""

    @property
    def label(self) -> str:
        mark = "+" if self.stance is Stance.INDICATES else "-"
        return f"{mark}{self.rule_id}"


@dataclass(frozen=True, slots=True)
class ModalityAssessment:
    """One modality, its verdict, and everything said about it."""

    modality: Modality
    verdict: Verdict
    indications: tuple[ModalityEvidence, ...] = ()
    contraindications: tuple[ModalityEvidence, ...] = ()

    @property
    def points(self) -> int:
        return sum(e.strength.points for e in self.indications)

    def summary(self) -> str:
        line = f"{self.modality:<32} {self.verdict}"
        if self.contraindications:
            line += f"  ({len(self.contraindications)} against)"
        return line


@dataclass(frozen=True, slots=True)
class ModalitySelection:
    """The module's output for one case."""

    gene: str
    mechanism: Mechanism
    assessments: tuple[ModalityAssessment, ...]
    ruleset_version: str = ""
    #: Set when the mechanism was not established and everything was blocked.
    blocked_reason: str = ""
    #: Things that do not rule a modality out but change how its verdict should
    #: be read — chiefly missing patient data that half the modalities depend on.
    caveats: tuple[str, ...] = ()

    @property
    def indicated(self) -> tuple[ModalityAssessment, ...]:
        """Ranked by accumulated indication strength, strongest first."""
        available = [a for a in self.assessments if a.verdict.is_available]
        return tuple(sorted(available, key=lambda a: a.points, reverse=True))

    @property
    def contraindicated(self) -> tuple[ModalityAssessment, ...]:
        return tuple(a for a in self.assessments if a.verdict is Verdict.CONTRAINDICATED)

    @property
    def is_blocked(self) -> bool:
        return bool(self.blocked_reason)

    def verdict_for(self, modality: Modality) -> Verdict:
        for assessment in self.assessments:
            if assessment.modality is modality:
                return assessment.verdict
        return Verdict.NOT_INDICATED
