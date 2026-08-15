"""What changed between two assessments of the same variant, and what caused it.

Two vocabularies, answering different questions.

``DeltaKind`` answers *what happened to the therapeutic direction*. This is
where the merge of the two halves of this project pays off: in isolation, a
reanalysis system watches an ACMG verdict move between five tiers, which is a
weaker signal than it looks. Here what is watched is the *mechanism* — and a
mechanism moving from loss-of-function to dominant-negative does not shift a
tier, it inverts the therapy. Adding a working copy goes from being the answer
to being the thing that makes it worse.

``AxisRole`` and ``AttributionPattern`` answer *why*: which coordinate of the
world had to move. Without that second vocabulary the system cannot separate
"ClinGen recurated this gene" from "we edited our own rule file", and those must
never be reported alike.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import StrEnum

from repairbench.modality import Modality, ModalitySelection
from repairbench.model import Confidence, Mechanism, MechanismCall
from repairbench.reanalysis.world import DriftAxis, World, all_axes


@dataclass(frozen=True, slots=True)
class Assessment:
    """Everything this project concluded about one variant, in one world.

    The unit that drifts. Pairing the mechanism with the modality sets is
    deliberate: a mechanism can hold steady while a modality is withdrawn — for
    instance when zygosity arrives and removes every route that needed an intact
    allele — and a system watching only the mechanism would report that week as
    quiet.

    The conclusions are held as plain values, and the ``call`` and ``selection``
    that produced them are optional detail. That split came from a real need
    rather than from taste: a scheduled run has to compare today's assessment
    with one made a month ago, which means writing it to disk and reading it
    back, and requiring the full objects for a comparison that reads four fields
    would have made persistence hard for no gain.
    """

    variant_key: str
    world: World
    mechanism: Mechanism
    confidence: Confidence
    indicated: frozenset[Modality]
    contraindicated: frozenset[Modality]
    #: What produced the conclusions above, when it is still in hand. Absent for
    #: an assessment read back from disk.
    call: MechanismCall | None = None
    selection: ModalitySelection | None = None

    @classmethod
    def of(
        cls,
        variant_key: str,
        world: World,
        call: MechanismCall,
        selection: ModalitySelection,
    ) -> Assessment:
        """Build from a live mechanism call and modality selection."""
        return cls(
            variant_key=variant_key,
            world=world,
            mechanism=call.mechanism,
            confidence=call.confidence,
            indicated=frozenset(a.modality for a in selection.indicated),
            contraindicated=frozenset(a.modality for a in selection.contraindicated),
            call=call,
            selection=selection,
        )


class DeltaKind(StrEnum):
    """The therapeutic meaning of a change."""

    MECHANISM_INVERTED = "mechanism_inverted"
    """One determined mechanism became a different one. The most consequential
    thing this system can report: loss-of-function to dominant-negative flips
    supplementation from the treatment to the hazard."""

    MODALITY_WITHDRAWN = "modality_withdrawn"
    """A route previously offered is now ruled out. The retraction analogue —
    something may already have been planned around it."""

    MECHANISM_RESOLVED = "mechanism_resolved"
    """Undetermined became determined. A direction now exists where none did."""

    MECHANISM_LOST = "mechanism_lost"
    """Determined became undetermined. We knew, and now we do not."""

    MODALITY_OPENED = "modality_opened"
    """A route that was not available now is."""

    CONFIDENCE_CHANGED = "confidence_changed"
    """Same mechanism, different footing. Invisible in a plan, visible before a
    reviewer signs one."""

    NONE = "none"
    """Nothing moved."""

    @property
    def inverts_direction(self) -> bool:
        """Does this change what somebody would already have been doing?"""
        return self in {DeltaKind.MECHANISM_INVERTED, DeltaKind.MODALITY_WITHDRAWN}


def all_delta_kinds() -> tuple[DeltaKind, ...]:
    return tuple(DeltaKind)


@dataclass(frozen=True, slots=True)
class AssessmentDelta:
    """The difference between two assessments of the same variant."""

    variant_key: str
    before: Assessment
    after: Assessment
    kind: DeltaKind

    @property
    def withdrawn_modalities(self) -> tuple[Modality, ...]:
        """Routes that were indicated and are now actively ruled out.

        Not merely 'no longer indicated'. A route that quietly dropped off the
        list is a change; a route that is now *contraindicated* is a warning,
        and only the second belongs in this list.
        """
        return tuple(
            sorted(self.before.indicated & self.after.contraindicated, key=lambda m: m.value)
        )

    @property
    def opened_modalities(self) -> tuple[Modality, ...]:
        return tuple(sorted(self.after.indicated - self.before.indicated, key=lambda m: m.value))

    @property
    def moved_axes(self) -> tuple[DriftAxis, ...]:
        return self.after.world.axes_differing_from(self.before.world)

    @property
    def is_material(self) -> bool:
        return self.kind is not DeltaKind.NONE

    @property
    def fingerprint(self) -> str:
        """Identity of *this change*, not of this variant.

        Two runs that detect the same transition produce the same fingerprint,
        which is how an acknowledged finding stays acknowledged instead of
        resurfacing at every release.
        """
        return (
            f"{self.variant_key}:{self.before.mechanism}->{self.after.mechanism}"
            f":{self.before.confidence}->{self.after.confidence}"
            f":{'|'.join(m.value for m in self.withdrawn_modalities)}"
            f":{'|'.join(m.value for m in self.opened_modalities)}"
        )


def compare(before: Assessment, after: Assessment) -> AssessmentDelta:
    """Diff two assessments, refusing to diff two different variants.

    The refusal matters more than it looks. The commonest way a reanalysis
    pipeline produces nonsense is joining records on a key that shifted
    underneath it — a lifted-over coordinate, a re-normalised indel, a new
    transcript version. Comparing on the identity itself turns that class of bug
    into an exception instead of a false "reclassified" alert.
    """
    if before.variant_key != after.variant_key:
        raise ValueError(f"cannot diff {before.variant_key} against {after.variant_key}")

    delta = AssessmentDelta(
        variant_key=before.variant_key, before=before, after=after, kind=DeltaKind.NONE
    )
    return AssessmentDelta(
        variant_key=before.variant_key,
        before=before,
        after=after,
        kind=_classify(before, after, delta),
    )


def _classify(before: Assessment, after: Assessment, delta: AssessmentDelta) -> DeltaKind:
    """Name the transition, most consequential first.

    The order *is* the logic, so it is written as an ordered table rather than
    buried in a ladder of conditionals. A run where the mechanism inverted and a
    modality was withdrawn is reported as an inversion, because that is the fact
    a reviewer has to act on and the withdrawal is its consequence.
    """
    determined_before = before.mechanism.is_determined
    determined_after = after.mechanism.is_determined

    ordered: tuple[tuple[bool, DeltaKind], ...] = (
        (
            determined_before and determined_after and before.mechanism is not after.mechanism,
            DeltaKind.MECHANISM_INVERTED,
        ),
        (bool(delta.withdrawn_modalities), DeltaKind.MODALITY_WITHDRAWN),
        (not determined_before and determined_after, DeltaKind.MECHANISM_RESOLVED),
        (determined_before and not determined_after, DeltaKind.MECHANISM_LOST),
        (bool(delta.opened_modalities), DeltaKind.MODALITY_OPENED),
        (before.confidence is not after.confidence, DeltaKind.CONFIDENCE_CHANGED),
    )
    for applies, kind in ordered:
        if applies:
            return kind
    return DeltaKind.NONE


class AxisRole(StrEnum):
    """How one axis relates to the change, established by re-running the rules."""

    DECISIVE = "decisive"
    """Moving it alone reproduces the change, and holding it back prevents it."""

    SUFFICIENT = "sufficient"
    """Moving it alone reproduces the change, but another axis would have too."""

    NECESSARY = "necessary"
    """Holding it back prevents the change, but alone it is not enough."""

    CONTRIBUTING = "contributing"
    """Moved, but changed nothing — it arrived in the same release window."""

    @property
    def is_causal(self) -> bool:
        return self is not AxisRole.CONTRIBUTING


def all_axis_roles() -> tuple[AxisRole, ...]:
    return tuple(AxisRole)


class AttributionPattern(StrEnum):
    """The shape of the causal story."""

    SOLE = "sole"
    REDUNDANT = "redundant"
    INTERACTION = "interaction"
    UNATTRIBUTED = "unattributed"


def all_attribution_patterns() -> tuple[AttributionPattern, ...]:
    return tuple(AttributionPattern)


@dataclass(frozen=True, slots=True)
class Attribution:
    """A delta together with the evidence for what caused it."""

    delta: AssessmentDelta
    roles: Mapping[DriftAxis, AxisRole]
    pattern: AttributionPattern
    evaluations_performed: int = 0

    @property
    def primary_axes(self) -> tuple[DriftAxis, ...]:
        ranked: list[DriftAxis] = []
        for wanted in (AxisRole.DECISIVE, AxisRole.SUFFICIENT, AxisRole.NECESSARY):
            ranked.extend(axis for axis in all_axes_in(self.roles) if self.roles[axis] is wanted)
        return tuple(dict.fromkeys(ranked)) or tuple(all_axes_in(self.roles))

    @property
    def rules_role(self) -> AxisRole | None:
        return self.roles.get(DriftAxis.RULES)

    @property
    def is_rule_change_implicated(self) -> bool:
        role = self.rules_role
        return role is not None and role.is_causal

    @property
    def is_purely_our_rules(self) -> bool:
        """No clinical axis played a causal role — we changed, the world did not.

        Axes that merely *moved* do not count. Releases arrive together, and a
        gnomAD refresh landing in the same week as a rule edit did not become a
        cause by being present; the counterfactuals already established that it
        changed nothing.
        """
        causal = [axis for axis, role in self.roles.items() if role.is_causal]
        return bool(causal) and all(not axis.is_clinical for axis in causal)

    def explain(self) -> str:
        """One sentence, in the register a reviewer reads."""
        if self.pattern is AttributionPattern.UNATTRIBUTED:
            return "no axis of the world explains this change — the pipeline is not reproducible"
        named = ", ".join(axis.value for axis in self.primary_axes)
        if self.pattern is AttributionPattern.SOLE:
            return f"caused by {named}"
        if self.pattern is AttributionPattern.REDUNDANT:
            return f"independently caused by any of: {named}"
        moved = ", ".join(axis.value for axis in self.delta.moved_axes)
        return f"no single axis suffices; required the combination of {moved}"


def all_axes_in(roles: Mapping[DriftAxis, AxisRole]) -> tuple[DriftAxis, ...]:
    """The moved axes in canonical order, so output is stable across runs."""
    return tuple(axis for axis in all_axes() if axis in roles)
