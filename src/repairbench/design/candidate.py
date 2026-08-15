"""What a designed edit is, and what has to be said alongside it.

The shape of these types is the argument of the module. A candidate is not a
guide sequence — it is a guide sequence, the editor it belongs to, every
unintended base the editor can also reach, and the reasons this package cannot
tell you how well it would work. Returning the sequence alone would be the
easiest thing here to get wrong, because a bare 20-mer looks finished.
"""

from __future__ import annotations

import textwrap
from dataclasses import dataclass, field

from repairbench.design.editors import Conversion, Editor


@dataclass(frozen=True, slots=True)
class Bystander:
    """A base the editor can also convert, inside the same window.

    Every one is an unintended change to the patient's genome, and the point of
    this type existing is that they are *listed*, individually, with their
    coordinates — never summarised as a count and never left out because the
    intended edit is the interesting one.
    """

    position_in_protospacer: int
    genomic_position: int
    becomes: str
    #: Whether this base sits in coding sequence, when a transcript was supplied.
    #: ``None`` means nobody looked, which is not the same as "no".
    in_coding_sequence: bool | None = None

    def describe(self) -> str:
        where = ""
        if self.in_coding_sequence is True:
            where = ", in coding sequence"
        elif self.in_coding_sequence is False:
            where = ", outside coding sequence"
        return (
            f"position {self.position_in_protospacer} (g.{self.genomic_position}) "
            f"→ {self.becomes}{where}"
        )


@dataclass(frozen=True, slots=True)
class EditCandidate:
    """One protospacer that places the target in one editor's window."""

    editor: Editor
    chromosome: str
    strand: str
    protospacer: str
    pam: str
    #: 1-based inclusive genomic span of the protospacer, low coordinate first,
    #: regardless of which strand it reads from.
    span: tuple[int, int]
    pam_span: tuple[int, int]
    target_position_in_protospacer: int
    target_genomic_position: int
    bystanders: tuple[Bystander, ...] = ()
    warnings: tuple[str, ...] = ()

    @property
    def conversion(self) -> Conversion:
        return self.editor.conversion

    @property
    def is_clean(self) -> bool:
        """No other editable base in the window.

        The only property here worth sorting on without a model attached, and
        even it is a statement about what the editor *can* reach rather than
        about what it would do.
        """
        return not self.bystanders

    @property
    def guide(self) -> str:
        """The protospacer as it would be ordered, 5' to 3'."""
        return self.protospacer

    def describe(self) -> str:
        lines = [
            f"{self.editor.id}  {self.protospacer} {self.pam.lower()}  "
            f"{self.chromosome}:{self.span[0]}-{self.span[1]} ({self.strand})",
            f"    target at position {self.target_position_in_protospacer} "
            f"(g.{self.target_genomic_position}), {self.conversion}",
        ]
        if self.bystanders:
            lines.append(f"    bystanders ({len(self.bystanders)}):")
            lines.extend(f"      · {bystander.describe()}" for bystander in self.bystanders)
        else:
            lines.append("    no other editable base in the window")
        lines.extend(
            textwrap.fill(
                warning, width=88, initial_indent="    ! ", subsequent_indent="      "
            )
            for warning in self.warnings
        )
        return "\n".join(lines)


@dataclass(frozen=True, slots=True)
class DesignOutcome:
    """Every candidate, every refusal, and the pins both were produced under.

    Refusals are first-class and carry their own text. "No candidates" has at
    least four different meanings here — the change is a transversion, no PAM
    sits at a usable distance, the reference is unresolved, no editor in the
    catalogue makes the required conversion — and collapsing them into an empty
    list would lose the only part a reader can act on.
    """

    gene: str
    chromosome: str
    position: int
    patient_base: str
    wild_type_base: str
    candidates: tuple[EditCandidate, ...] = ()
    refusals: tuple[str, ...] = ()
    notes: tuple[str, ...] = ()
    catalogue_pin: str = ""
    #: What the efficiency model said, or why there is none. Never a number this
    #: package invented.
    ranking: str = ""
    considered: tuple[str, ...] = field(default_factory=tuple)

    @property
    def has_candidates(self) -> bool:
        return bool(self.candidates)

    @property
    def clean(self) -> tuple[EditCandidate, ...]:
        return tuple(candidate for candidate in self.candidates if candidate.is_clean)
