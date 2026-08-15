"""The world a mechanism call was made in.

Eight coordinates. Seven describe what the field knows; one — ``rules`` —
describes what *we* decided, and keeping those apart is the point of the whole
module.

Adding an axis is deliberately expensive: ``World`` refuses to be built without
every one of them, so a new axis breaks every existing world until each is
re-pinned. That cost is the feature. An axis nobody had to acknowledge would be
an input the system reads without naming.
"""

from __future__ import annotations

import hashlib
from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from datetime import date
from enum import StrEnum

from repairbench.model import RepairbenchError


class WorldError(RepairbenchError):
    """A world is incomplete, or contradicts itself."""


class IncompleteWorldError(WorldError):
    """A world was assembled without pinning every axis the rules read.

    A call that cannot name the version of everything it consulted cannot be
    reproduced, and an irreproducible call cannot be diffed against a later one —
    the difference would be unattributable by construction.
    """


class PinConflictError(WorldError):
    """The same release was presented with two different content digests.

    Releases are immutable by convention. When ``clinvar@2026-07`` has one digest
    today and another tomorrow, one of the two ingests is wrong, and silently
    preferring either would corrupt every historical comparison.
    """


class DriftAxis(StrEnum):
    """An independent coordinate along which a mechanism call can move."""

    CLINVAR = "clinvar"
    """Submitted assertions and their review status."""

    POPULATION_FREQUENCY = "population_frequency"
    """gnomAD-style frequencies and constraint. Moves the LOEUF rules."""

    GENE_CURATION = "gene_curation"
    """ClinGen dosage sensitivity and curated mechanisms. Moves almost everything."""

    PANEL = "panel"
    """Which genes a panel contains, at what confidence."""

    PHENOTYPE = "phenotype"
    """The patient's own HPO terms. The one axis that moves on a ward round
    rather than on a release schedule."""

    EXPRESSION = "expression"
    """Where genes are switched on. A new GTEx release can move the rule that
    asks whether this gene is transcribed in the affected tissue at all — which
    is a way for an answer to change without anything about the patient
    changing, so it earns an axis rather than riding on another one."""

    ANNOTATION = "annotation"
    """The transcript structures. A new transcript version moves the exon
    boundaries, which moves the NMD boundary, which can invert the mechanism —
    this axis exists because that is a real way for an answer to change without
    anybody learning anything."""

    RULES = "rules"
    """Our own rule files. Not evidence about the world — evidence about us."""

    @property
    def is_clinical(self) -> bool:
        """Does movement here mean the field learned something?

        ``RULES`` is the odd one out and the reason this property exists. When a
        mechanism changes because we edited a rule file, no new fact about the
        patient exists, and reporting it as a finding would inflate the yield of
        the software with its own corrections.
        """
        return self is not DriftAxis.RULES

    @property
    def is_case_scoped(self) -> bool:
        """True when the axis moves per patient rather than per release."""
        return self is DriftAxis.PHENOTYPE


def all_axes() -> tuple[DriftAxis, ...]:
    """Every axis, in the order a report presents them."""
    return tuple(DriftAxis)


@dataclass(frozen=True, slots=True, order=True)
class Pin:
    """One axis fixed at one version, with the digest that proves which.

    ``version`` is the human-facing label a report cites; ``digest`` is the
    content hash that makes the citation checkable. Both are required.
    """

    axis: DriftAxis
    version: str
    digest: str
    released_at: date | None = None

    def __post_init__(self) -> None:
        if not self.version.strip():
            raise WorldError(f"pin for {self.axis} has no version label")
        if not self.digest.strip():
            raise WorldError(f"pin for {self.axis} has no content digest")

    @property
    def label(self) -> str:
        return f"{self.axis.value}@{self.version}"

    def same_source(self, other: Pin) -> bool:
        """Do two pins name the same release? Release dates are metadata about
        when we learned of a release, not part of its identity."""
        return (self.axis, self.version, self.digest) == (other.axis, other.version, other.digest)


@dataclass(frozen=True, slots=True)
class World:
    """A complete, immutable set of pins — one per axis."""

    pins: tuple[Pin, ...]

    def __post_init__(self) -> None:
        seen: set[DriftAxis] = set()
        for pin in self.pins:
            if pin.axis in seen:
                raise PinConflictError(f"{pin.axis} is pinned twice in one world")
            seen.add(pin.axis)
        missing = [axis.value for axis in all_axes() if axis not in seen]
        if missing:
            raise IncompleteWorldError(f"world does not pin: {', '.join(missing)}")
        if tuple(sorted(self.pins)) != self.pins:
            raise WorldError("pins must be sorted — use World.of() rather than the constructor")

    @classmethod
    def of(cls, pins: Iterable[Pin]) -> World:
        """Assemble from pins in any order, rejecting contradictions."""
        collected: dict[DriftAxis, Pin] = {}
        for pin in pins:
            previous = collected.get(pin.axis)
            if previous is not None and not previous.same_source(pin):
                raise PinConflictError(
                    f"{pin.axis} pinned as {previous.version}/{previous.digest[:8]} "
                    f"and {pin.version}/{pin.digest[:8]} in the same world"
                )
            collected[pin.axis] = pin
        return cls(tuple(sorted(collected.values())))

    def __iter__(self) -> Iterator[Pin]:
        return iter(self.pins)

    @property
    def by_axis(self) -> Mapping[DriftAxis, Pin]:
        return {pin.axis: pin for pin in self.pins}

    def pin_for(self, axis: DriftAxis) -> Pin:
        for pin in self.pins:
            if pin.axis is axis:
                return pin
        raise WorldError(f"world does not pin {axis}")

    def with_pin(self, pin: Pin) -> World:
        """The same world with one axis moved — the counterfactual primitive."""
        return World.of([*(p for p in self.pins if p.axis is not pin.axis), pin])

    def axes_differing_from(self, other: World) -> tuple[DriftAxis, ...]:
        """Which coordinates moved. Attribution never looks further."""
        mine, theirs = self.by_axis, other.by_axis
        return tuple(axis for axis in all_axes() if not mine[axis].same_source(theirs[axis]))

    @property
    def digest(self) -> str:
        """A stable hash of the whole world. Two equal digests must reproduce."""
        canonical = "\n".join(f"{p.axis.value}\t{p.version}\t{p.digest}" for p in self.pins)
        return hashlib.sha256(canonical.encode()).hexdigest()

    @property
    def short_digest(self) -> str:
        return self.digest[:12]

    def describe(self) -> str:
        """The citation line for a report footer."""
        return " ".join(pin.label for pin in self.pins)
