"""Rendering a design so that what is missing is as visible as what is there.

The ordering of this report is the argument. Refusals print before candidates,
bystanders print under the candidate that causes them rather than in a summary,
and the line saying nothing here is ranked by efficiency prints whether or not
anybody reads it. A design report that leads with a 20-mer invites somebody to
order the 20-mer.
"""

from __future__ import annotations

from repairbench.design.aso import AsoOutcome
from repairbench.design.candidate import DesignOutcome
from repairbench.design.prime import PrimeOutcome
from repairbench.design.risk import OffTargetAssessment, RiskTier
from repairbench.modality import Verdict
from repairbench.plan import Plan

_INDENT = "  "


def render_design(outcome: DesignOutcome) -> str:
    lines = [
        f"{outcome.gene}  {outcome.chromosome}:{outcome.position}  "
        f"{outcome.patient_base}→{outcome.wild_type_base}",
    ]

    for note in outcome.notes:
        lines.append(f"{_INDENT}note        {_wrap(note)}")

    if outcome.refusals:
        lines.append(f"{_INDENT}not designed")
        lines.extend(f"{_INDENT}  · {_wrap(refusal)}" for refusal in outcome.refusals)

    if outcome.candidates:
        clean = len(outcome.clean)
        lines.append(
            f"{_INDENT}candidates  {len(outcome.candidates)} "
            f"({clean} with no other editable base in the window)"
        )
        for candidate in outcome.candidates:
            lines.append(_INDENT + "  " + candidate.describe().replace("\n", f"\n{_INDENT}  "))

    if outcome.considered:
        lines.append(f"{_INDENT}considered  {', '.join(outcome.considered)}")
    if outcome.ranking:
        lines.append(f"{_INDENT}ranking     {_wrap(outcome.ranking)}")
    lines.append(f"{_INDENT}catalogue   {outcome.catalogue_pin}")
    lines.append(
        f"{_INDENT}note        a candidate here is a placement, not a therapy. Nothing above "
        f"\n{_INDENT}            accounts for off-target activity, delivery, or what a bystander "
        f"\n{_INDENT}            edit does to the protein."
    )
    return "\n".join(lines)


def render_offtarget(assessment: OffTargetAssessment) -> str:
    """The ranked hit list, worst first, with the reason each one ranks where it does."""
    lines = [
        f"off-target review  {assessment.guide}",
        f"{_INDENT}hits        {len(assessment.hits)} read from {assessment.source_pin}",
    ]
    if assessment.scoring:
        lines.append(f"{_INDENT}scoring     {_wrap(assessment.scoring)}")

    for tier in RiskTier.worst_first():
        in_tier = [hit for hit in assessment.hits if hit.tier is tier]
        if not in_tier:
            continue
        lines.append(f"{_INDENT}{tier.value}  ({len(in_tier)})")
        for hit in in_tier:
            count = hit.hit.mismatches
            mismatches = "1 mismatch" if count == 1 else f"{count} mismatches"
            lines.append(
                f"{_INDENT}  {hit.hit.chromosome}:{hit.hit.position} ({hit.hit.strand})  "
                f"{mismatches}  {hit.where}"
            )
            for evidence in hit.evidence:
                lines.append(f"{_INDENT}      {evidence.rule_id}")
                lines.append(f"{_INDENT}        {_wrap(evidence.because, indent=' ' * 8)}")

    if assessment.unranked:
        lines.append(f"{_INDENT}unranked    {_wrap(assessment.unranked)}")
    lines.append(f"{_INDENT}rules       {assessment.ruleset_pin}")
    return "\n".join(lines)


def _wrap(text: str, width: int = 76, indent: str = " " * 14) -> str:
    words, lines, current = text.split(), [], ""
    for word in words:
        if len(current) + len(word) + 1 > width:
            lines.append(current)
            current = word
        else:
            current = f"{current} {word}".strip()
    lines.append(current)
    return f"\n{indent}".join(lines)


def render_pegrnas(outcome: PrimeOutcome, *, limit: int = 5) -> str:
    """The pegRNA report, truncated on purpose and saying by how much.

    A single edit admits hundreds of pegRNAs, most differing only in primer
    length. Printing them all buries the choices that matter; printing a few
    without saying how many there were would imply the list is the design space
    rather than a slice of it.
    """
    lines = [
        f"{outcome.gene}  {outcome.chromosome}:{outcome.position}  {outcome.edit}",
    ]
    for note in outcome.notes:
        lines.append(f"{_INDENT}note        {_wrap(note)}")

    if outcome.refusals:
        lines.append(f"{_INDENT}not designed")
        lines.extend(f"{_INDENT}  · {_wrap(refusal)}" for refusal in outcome.refusals)

    if outcome.candidates:
        blocked = len(outcome.candidates) - len(outcome.usable)
        lines.append(
            f"{_INDENT}pegRNAs     {len(outcome.candidates)} across {outcome.protospacers} "
            f"protospacer(s); {blocked} blocked by a rule"
        )
        for candidate in outcome.usable[:limit]:
            lines.append(_INDENT + "  " + candidate.describe().replace("\n", f"\n{_INDENT}  "))
        remaining = len(outcome.usable) - limit
        if remaining > 0:
            lines.append(f"{_INDENT}  … and {remaining} more, not shown")

    if outcome.ranking:
        lines.append(f"{_INDENT}ranking     {_wrap(outcome.ranking)}")
    lines.append(f"{_INDENT}rules       {outcome.ruleset_pin}")
    return "\n".join(lines)


def render_asos(outcome: AsoOutcome, *, limit: int = 8) -> str:
    """The tiling report.

    Leads with how many windows were tiled and how many survived, because the
    ratio is the only honest summary of a tiling run: a list of two hundred
    oligonucleotides reads as an answer, and it is a starting point.
    """
    lines = [
        f"{outcome.gene}  {outcome.chromosome}:{outcome.span[0]}-{outcome.span[1]}  "
        f"{outcome.chemistry.id} ({outcome.chemistry.action}, {outcome.chemistry.length} nt)",
        f"{_INDENT}tiled       {outcome.tiled} windows; {len(outcome.usable)} with nothing "
        f"blocking, {len(outcome.candidates) - len(outcome.usable)} blocked",
    ]
    for note in outcome.notes:
        lines.append(f"{_INDENT}note        {_wrap(note)}")

    for candidate in outcome.usable[:limit]:
        lines.append(_INDENT + "  " + candidate.describe().replace("\n", f"\n{_INDENT}  "))
    remaining = len(outcome.usable) - limit
    if remaining > 0:
        lines.append(f"{_INDENT}  … and {remaining} more, not shown")

    blocked = [candidate for candidate in outcome.candidates if candidate.is_blocked]
    if blocked:
        lines.append(f"{_INDENT}blocked     {len(blocked)}, worst first")
        for candidate in blocked[:2]:
            lines.append(_INDENT + "  " + candidate.describe().replace("\n", f"\n{_INDENT}  "))

    if outcome.ranking:
        lines.append(f"{_INDENT}ranking     {_wrap(outcome.ranking)}")
    lines.append(f"{_INDENT}rules       {outcome.ruleset_pin}")
    return "\n".join(lines)


def render_plan(plan: Plan, *, limit: int = 3) -> str:
    """The whole case in one document: why, what class, which molecule.

    Ordered the way a sign-out is read rather than the way the pipeline runs.
    The mechanism first, because everything below is conditional on it. Then
    what was ruled out — before what was designed, for the same reason M6 prints
    its refusals first: the output that changes what somebody does is the one
    that closes a door. Then the molecules. Then the pins, because a plan that
    cannot name the files it was made under cannot be compared with a later one.
    """
    lines = [
        f"{plan.gene}  {plan.call.transcript}",
        f"{_INDENT}mechanism   {plan.call.mechanism} ({plan.call.confidence})",
    ]
    for evidence in plan.call.evidence:
        lines.append(f"{_INDENT * 2}{evidence.rule_id} [{evidence.strength}]")
    if plan.call.needs_review:
        lines.append(f"{_INDENT * 2}· needs review: the rules that fired do not agree")

    ruled_out = plan.ruled_out
    if ruled_out:
        lines.append(f"{_INDENT}not designed, ruled out by the modality rules")
        for design in ruled_out:
            lines.append(f"{_INDENT * 2}{design.modality}")

    designed = plan.designed
    if designed:
        lines.append(f"{_INDENT}designed")
        for design in designed:
            lines.append(f"{_INDENT * 2}{design.modality} → {design.designer}")
            body = _render_outcome(design.outcome, limit=limit)
            lines.append("\n".join(f"{_INDENT * 3}{line}" for line in body.splitlines()))

    considered = [
        design
        for design in plan.designs
        if not design.designed and design.verdict is not Verdict.CONTRAINDICATED
    ]
    if considered:
        lines.append(f"{_INDENT}considered, nothing designed")
        for design in considered:
            lines.append(f"{_INDENT * 2}{design.modality}  ({design.verdict})")
            lines.append(f"{_INDENT * 3}{_wrap(design.refusal, indent=' ' * 8)}")

    for note in plan.notes:
        lines.append(f"{_INDENT}note        {_wrap(note)}")

    lines.append(f"{_INDENT}made under")
    lines.extend(f"{_INDENT * 2}{pin}" for pin in plan.pins)
    return "\n".join(lines)


def _render_outcome(outcome: object, *, limit: int) -> str:
    """Whichever designer produced this, rendered by its own report."""
    if isinstance(outcome, DesignOutcome):
        return render_design(outcome)
    if isinstance(outcome, PrimeOutcome):
        return render_pegrnas(outcome, limit=limit)
    if isinstance(outcome, AsoOutcome):
        return render_asos(outcome, limit=limit)
    return "nothing was designed"
