"""Rendering a mechanism call for a human.

The layout follows one principle: the reasoning comes before the conclusion is
usable. A reader who stops after the first line has the answer; a reader who
needs to disagree with it has, in the next few lines, everything they would need
to say where the disagreement is — which rule, on what grounds, from which
version of the rule file.
"""

from __future__ import annotations

from repairbench.modality import ModalitySelection
from repairbench.model import Confidence, MechanismCall

_INDENT = "  "


def render(call: MechanismCall) -> str:
    """Render a call as plain text."""
    lines: list[str] = []

    lines.append(f"{call.gene}  {call.transcript}")
    lines.append(f"{_INDENT}mechanism   {call.mechanism}  ({call.confidence})")
    if call.confidence is Confidence.ESTABLISHED:
        lines.append(
            f"{_INDENT}            an expert curation settled this; the rules below still ran"
        )

    if call.evidence:
        lines.append(f"{_INDENT}because")
        for evidence in call.evidence:
            lines.append(f"{_INDENT * 2}{evidence.rule_id} [{evidence.strength}]")
            lines.append(f"{_INDENT * 3}{evidence.because}")
            if evidence.citation:
                lines.append(f"{_INDENT * 3}— {evidence.citation}")
    else:
        lines.append(f"{_INDENT}because     no rule fired for this variant")

    if call.conflicts:
        lines.append(f"{_INDENT}but")
        for conflict in call.conflicts:
            argues = (
                "argues the evidence does not settle it"
                if not conflict.supports.is_determined
                else f"argues for {conflict.supports}"
            )
            lines.append(f"{_INDENT * 2}{conflict.rule_id} [{conflict.strength}] {argues}")
            lines.append(f"{_INDENT * 3}{conflict.because}")

    feasibility = call.feasibility
    lines.append(f"{_INDENT}not ruled out on mechanistic grounds")
    for label, value in (
        ("gene addition addresses the mechanism", feasibility.gene_addition_coherent),
        ("coding sequence fits a viral payload", feasibility.fits_viral_payload),
        ("silenced allele available to reactivate", feasibility.silenced_allele_available),
        ("allele-specific silencing indicated", feasibility.allele_specific_silencing_indicated),
    ):
        lines.append(f"{_INDENT * 2}{'yes' if value else 'no ':<4} {label}")
    if feasibility.exon_skipping_preserves_frame is not None:
        answer = "yes" if feasibility.exon_skipping_preserves_frame else "no "
        lines.append(f"{_INDENT * 2}{answer:<4} skipping the affected exon preserves the frame")

    if feasibility.notes:
        lines.append(f"{_INDENT}notes")
        for note in feasibility.notes:
            lines.append(f"{_INDENT * 2}· {note}")

    lines.append(f"{_INDENT}ruleset     {call.ruleset_version}")
    if call.needs_review:
        lines.append(f"{_INDENT}review      required before any downstream use")

    return "\n".join(lines)


def render_selection(selection: ModalitySelection) -> str:
    """Render a modality assessment.

    Contraindications are printed before indications, in defiance of how a
    reader would like to see it. That ordering is the module's opinion: the
    output that changes what somebody does is the one that rules a route out,
    and burying it under a list of possibilities would invert the priority the
    engine works to.
    """
    lines: list[str] = [f"{selection.gene}  mechanism {selection.mechanism}"]

    if selection.is_blocked:
        lines.append(f"{_INDENT}nothing assessed")
        lines.append(f"{_INDENT * 2}{selection.blocked_reason}")
        lines.append(f"{_INDENT}ruleset     {selection.ruleset_version}")
        return "\n".join(lines)

    contraindicated = selection.contraindicated
    if contraindicated:
        lines.append(f"{_INDENT}ruled out")
        for assessment in contraindicated:
            lines.append(f"{_INDENT * 2}{assessment.modality}")
            for evidence in assessment.contraindications:
                lines.append(f"{_INDENT * 3}{evidence.rule_id}: {evidence.because}")

    indicated = selection.indicated
    if indicated:
        lines.append(f"{_INDENT}not ruled out, in order of accumulated support")
        for rank, assessment in enumerate(indicated, start=1):
            lines.append(f"{_INDENT * 2}{rank}. {assessment.modality}  ({assessment.points} pts)")
            for evidence in assessment.indications:
                lines.append(f"{_INDENT * 3}{evidence.rule_id} [{evidence.strength}]")
                lines.append(f"{_INDENT * 4}{evidence.because}")
                if evidence.citation:
                    lines.append(f"{_INDENT * 4}— {evidence.citation}")
    else:
        lines.append(f"{_INDENT}no modality is coherent with this mechanism")

    if selection.caveats:
        lines.append(f"{_INDENT}read with caution")
        for caveat in selection.caveats:
            lines.append(f"{_INDENT * 2}· {caveat}")

    lines.append(f"{_INDENT}ruleset     {selection.ruleset_version}")
    lines.append(
        f"{_INDENT}note        'not ruled out' is not a recommendation. Delivery, tissue, "
        "dose and safety are all downstream of this and none of them are assessed here."
    )
    return "\n".join(lines)
