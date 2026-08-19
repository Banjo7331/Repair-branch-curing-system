"""Running the reanalysis reference set: episodes rather than mechanisms.

The mechanism and modality reference sets ask whether the rules reach the answer
the literature reached. This one cannot ask that, because reanalysis has no
answers of that kind — it has *episodes*. A release lands, an assessment moves or
does not, and what is under test is the causal claim the system makes and who it
tells.

Shared between the test suite and the ``reference --reanalysis`` command for the
same reason the other two runners are: a reference set that only runs under
pytest is a reference set nobody looks at, and one that only runs from the
command line is one CI does not enforce.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from repairbench.model import Consequence, MissenseDistribution, RepairbenchError, Zygosity
from repairbench.reanalysis.attribution import Attributor
from repairbench.reanalysis.catalogue import SourceCatalogue
from repairbench.reanalysis.drift import AttributionPattern, DeltaKind, compare
from repairbench.reanalysis.engine import RepairbenchEngine, WatchedVariant
from repairbench.reanalysis.routing import ReviewQueue, SurfacingDecision, Urgency
from repairbench.reanalysis.surfacing import SurfacingPolicy
from repairbench.reanalysis.world import DriftAxis, Pin, World


class ReferenceSetError(RepairbenchError):
    """The reference set or the deployment it names is malformed."""


@dataclass(frozen=True, slots=True)
class Episode:
    """One week in a laboratory, as the reference set describes it."""

    name: str
    note: str
    moves: dict[str, str]
    expect: dict[str, Any]


@dataclass(frozen=True, slots=True)
class EpisodeResult:
    """What the system did, and where that differs from what was expected."""

    episode: Episode
    kind: DeltaKind
    pattern: AttributionPattern | None
    decision: SurfacingDecision
    causal_axes: tuple[DriftAxis, ...]
    contributing_axes: tuple[DriftAxis, ...]
    mismatches: tuple[str, ...] = ()

    @property
    def reproduced(self) -> bool:
        return not self.mismatches

    def summarise(self) -> str:
        axes = ", ".join(axis.value for axis in self.causal_axes) or "nothing"
        return f"{self.kind} → {axes} → {self.decision.queue} ({self.decision.urgency})"


def load_episodes(path: str | Path) -> tuple[dict[str, Any], dict[str, str], list[Episode]]:
    """Read the set: the variant it watches, the baseline world, the episodes."""
    document = yaml.safe_load(Path(path).read_text())
    for required in ("variant", "baseline", "episodes"):
        if required not in document:
            raise ReferenceSetError(f"{path}: reference set has no {required!r}")

    episodes = [
        Episode(
            name=str(entry["name"]),
            note=" ".join(str(entry.get("note", "")).split()),
            moves=dict(entry.get("moves") or {}),
            expect=dict(entry.get("expect") or {}),
        )
        for entry in document["episodes"]
    ]
    return document["variant"], dict(document["baseline"]), episodes


def watched_from(spec: dict[str, Any]) -> WatchedVariant:
    return WatchedVariant(
        key=str(spec["key"]),
        gene=str(spec["gene"]),
        consequence=Consequence(spec["consequence"]),
        zygosity=Zygosity(spec.get("zygosity", "unknown")),
        cds_position=int(spec["cds_position"]),
        distribution=MissenseDistribution(**(spec.get("distribution") or {})),
    )


def world_from(
    catalogue: SourceCatalogue, baseline: dict[str, str], moves: dict[str, str] | None = None
) -> World:
    """Assemble a world from named versions, with the phenotype pinned by hand.

    The phenotype is the one axis no catalogue holds: it moves on a ward round
    rather than on a release schedule, so it is pinned to the case rather than
    read from a file.
    """
    versions = {**baseline, **(moves or {})}
    pins = [
        catalogue.pin_for(axis, versions[axis.value])
        for axis in DriftAxis
        if not axis.is_case_scoped and axis.value in versions
    ]
    return World.of([*pins, Pin(axis=DriftAxis.PHENOTYPE, version="day-1", digest="hpo-day-1")])


def run_episode(
    episode: Episode,
    engine: RepairbenchEngine,
    catalogue: SourceCatalogue,
    baseline: dict[str, str],
    variant_key: str,
) -> EpisodeResult:
    """Run one episode end to end and check it against what was expected.

    Nothing here is mocked: the counterfactual probes reload the older files off
    disk, which is the only way "we re-ran it with January's curation" can be a
    true sentence rather than a plausible one.
    """
    before = engine.assess(variant_key, world_from(catalogue, baseline))
    after = engine.assess(variant_key, world_from(catalogue, baseline, episode.moves))
    delta = compare(before, after)

    pattern: AttributionPattern | None = None
    causal: tuple[DriftAxis, ...] = ()
    contributing: tuple[DriftAxis, ...] = ()

    if delta.is_material:
        attribution = Attributor().attribute(delta, engine.assess)
        pattern = attribution.pattern
        causal = tuple(
            axis for axis, role in sorted(attribution.roles.items()) if role.is_causal
        )
        contributing = tuple(
            axis for axis, role in sorted(attribution.roles.items()) if not role.is_causal
        )
        decision = SurfacingPolicy().decide(attribution)
    else:
        decision = SurfacingDecision(Urgency.SILENT, ReviewQueue.NONE, "nothing changed")

    return EpisodeResult(
        episode=episode,
        kind=delta.kind,
        pattern=pattern,
        decision=decision,
        causal_axes=causal,
        contributing_axes=contributing,
        mismatches=_check(episode.expect, delta.kind, pattern, decision, causal, contributing),
    )


def _check(
    expect: dict[str, Any],
    kind: DeltaKind,
    pattern: AttributionPattern | None,
    decision: SurfacingDecision,
    causal: tuple[DriftAxis, ...],
    contributing: tuple[DriftAxis, ...],
) -> tuple[str, ...]:
    """Compare what happened with what the set says should have.

    Every expectation is optional and only the stated ones are checked, so an
    episode can pin the causal claim without also pinning an urgency nobody has
    an opinion about. What it may not do is pass by saying nothing: an episode
    with an empty ``expect`` is a case that tests nothing, and is refused.
    """
    if not expect:
        return ("the episode states no expectations, so it cannot fail",)

    problems: list[str] = []
    if "kind" in expect and kind.value != expect["kind"]:
        problems.append(f"kind is {kind.value}, expected {expect['kind']}")
    if "pattern" in expect and (pattern.value if pattern else None) != expect["pattern"]:
        problems.append(f"pattern is {pattern}, expected {expect['pattern']}")
    if "queue" in expect and decision.queue.value != expect["queue"]:
        problems.append(f"queue is {decision.queue.value}, expected {expect['queue']}")
    if "urgency" in expect and decision.urgency.value != expect["urgency"]:
        problems.append(f"urgency is {decision.urgency.value}, expected {expect['urgency']}")
    if "rule_change_caveat" in expect and decision.rule_change_caveat != expect[
        "rule_change_caveat"
    ]:
        problems.append(
            f"rule_change_caveat is {decision.rule_change_caveat}, "
            f"expected {expect['rule_change_caveat']}"
        )

    for field, actual in (("causal_axes", causal), ("contributing_axes", contributing)):
        if field in expect:
            names = sorted(axis.value for axis in actual)
            if names != sorted(expect[field]):
                problems.append(f"{field} are {names}, expected {sorted(expect[field])}")

    return tuple(problems)


def run_reference_set(
    set_path: str | Path, deployment: str | Path
) -> tuple[EpisodeResult, ...]:
    """Every episode, against a deployment on disk."""
    variant_spec, baseline, episodes = load_episodes(set_path)
    catalogue = SourceCatalogue.load(Path(deployment) / "catalogue.yaml")
    watched = watched_from(variant_spec)
    engine = RepairbenchEngine(catalogue, {watched.key: watched})

    return tuple(
        run_episode(episode, engine, catalogue, baseline, watched.key) for episode in episodes
    )
