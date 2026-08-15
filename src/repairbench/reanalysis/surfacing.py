"""Which changes reach a human, with what urgency, and to which queue.

Reanalysis has one predictable failure mode, and it is not missing something —
it is finding so many immaterial changes that the laboratory stops reading the
output. A monthly ClinVar release moves hundreds of records; almost none of them
change what anyone should do.

Three rules are worth stating out loud:

* **An inverted mechanism outranks everything.** Loss-of-function becoming
  dominant-negative does not shift a tier — it turns supplementation from the
  treatment into the hazard, and anything planned on the old reading is now
  wrong.
* **A change our own rule file caused is not a clinical finding.** It goes to
  validation. It may well *become* one once someone confirms the new rule is
  right, but the confirmation comes first — otherwise our corrections get
  counted as discoveries.
* **A withdrawn route is louder than a new one.** An opened modality is an
  opportunity nobody has acted on yet; a withdrawn one may already be in a plan.
"""

from __future__ import annotations

from dataclasses import dataclass

from repairbench.reanalysis.drift import Attribution, AttributionPattern, DeltaKind
from repairbench.reanalysis.routing import ReviewQueue, SurfacingDecision, Urgency


@dataclass(frozen=True, slots=True)
class SurfacingPolicy:
    """The rules, in one place, testable one by one."""

    def decide(
        self,
        attribution: Attribution,
        acknowledged: frozenset[str] = frozenset(),
    ) -> SurfacingDecision:
        delta = attribution.delta

        if not delta.is_material:
            return SurfacingDecision(Urgency.SILENT, ReviewQueue.NONE, "nothing changed")

        if delta.fingerprint in acknowledged:
            return SurfacingDecision(
                Urgency.SILENT,
                ReviewQueue.NONE,
                "this exact transition has already been reviewed and signed out",
                suppressed=True,
            )

        if attribution.pattern is AttributionPattern.UNATTRIBUTED:
            return SurfacingDecision(
                Urgency.HIGH,
                ReviewQueue.VALIDATION,
                "the assessment moved while every pin stood still — the pipeline is not "
                "reproducible and no clinical conclusion may be drawn from this run",
            )

        if attribution.is_purely_our_rules:
            return self._rule_change_only(attribution)
        return self._clinical(attribution)

    def _rule_change_only(self, attribution: Attribution) -> SurfacingDecision:
        kind = attribution.delta.kind
        if kind.inverts_direction:
            return SurfacingDecision(
                Urgency.CRITICAL,
                ReviewQueue.VALIDATION,
                "our own rule change reverses a direction that may already have been "
                "acted on; the rule must be confirmed before anything is amended",
                rule_change_caveat=True,
            )
        if kind in {DeltaKind.MECHANISM_RESOLVED, DeltaKind.MECHANISM_LOST}:
            return SurfacingDecision(
                Urgency.ROUTINE,
                ReviewQueue.VALIDATION,
                "the mechanism changed status because the rules changed, not because the "
                "evidence did — validate the rule before claiming anything from it",
                rule_change_caveat=True,
            )
        return SurfacingDecision(
            Urgency.LOW,
            ReviewQueue.VALIDATION,
            "rule behaviour changed without changing the therapeutic direction",
            rule_change_caveat=True,
        )

    def _clinical(self, attribution: Attribution) -> SurfacingDecision:
        """Route a change the world caused.

        A table rather than a ladder: which kinds reach a clinician, and at what
        urgency, is the policy this module exists to state, and it should be
        readable at a glance rather than reconstructed from control flow.
        """
        delta = attribution.delta
        caveat = attribution.is_rule_change_implicated
        cause = attribution.explain()
        withdrawn = ", ".join(m.value for m in delta.withdrawn_modalities)
        opened = ", ".join(m.value for m in delta.opened_modalities)

        routings: dict[DeltaKind, tuple[Urgency, ReviewQueue, str]] = {
            DeltaKind.MECHANISM_INVERTED: (
                Urgency.CRITICAL,
                ReviewQueue.CLINICAL_SIGNOUT,
                f"the mechanism inverted from {delta.before.mechanism} to "
                f"{delta.after.mechanism} — every modality below it changes with it",
            ),
            DeltaKind.MODALITY_WITHDRAWN: (
                Urgency.CRITICAL,
                ReviewQueue.CLINICAL_SIGNOUT,
                f"a route previously offered is now ruled out ({withdrawn})",
            ),
            DeltaKind.MECHANISM_RESOLVED: (
                Urgency.HIGH,
                ReviewQueue.CLINICAL_SIGNOUT,
                "a mechanism was established where none was before",
            ),
            DeltaKind.MECHANISM_LOST: (
                Urgency.HIGH,
                ReviewQueue.CLINICAL_SIGNOUT,
                "the mechanism no longer resolves; anything downstream of it is "
                "unsupported until it does",
            ),
            DeltaKind.MODALITY_OPENED: (
                Urgency.ROUTINE,
                ReviewQueue.CLINICAL_SIGNOUT,
                f"a route not previously available is now open ({opened})",
            ),
            DeltaKind.CONFIDENCE_CHANGED: (
                Urgency.LOW,
                ReviewQueue.WATCHLIST,
                f"the mechanism holds but its footing moved "
                f"({delta.before.confidence} to {delta.after.confidence})",
            ),
        }

        routing = routings.get(delta.kind)
        if routing is None:
            return SurfacingDecision(Urgency.SILENT, ReviewQueue.NONE, "nothing changed")
        urgency, queue, reason = routing
        return SurfacingDecision(urgency, queue, f"{reason} — {cause}", rule_change_caveat=caveat)
