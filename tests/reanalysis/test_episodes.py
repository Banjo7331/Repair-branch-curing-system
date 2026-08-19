"""The reanalysis reference set, run against files on disk.

The unit tests around it establish that each piece behaves — that a
counterfactual loads the older file, that a rule change routes to validation,
that an acknowledged transition does not resurface. What they cannot establish
is that the *episodes* come out right: a week in a laboratory is several of
those pieces at once, and the interesting failures live in the combinations.

Each case here is one such week. The set is data, so an episode can be added
without touching Python, and the same runner backs ``repairbench reference
--reanalysis`` so the file is not a thing only pytest ever reads.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from repairbench.reanalysis.catalogue import SourceCatalogue
from repairbench.reanalysis.drift import DeltaKind
from repairbench.reanalysis.engine import RepairbenchEngine
from repairbench.reanalysis.reference import (
    Episode,
    EpisodeResult,
    ReferenceSetError,
    load_episodes,
    run_episode,
    run_reference_set,
    watched_from,
)
from repairbench.reanalysis.routing import ReviewQueue, Urgency
from repairbench.reanalysis.world import DriftAxis

SET = Path(__file__).parents[1] / "reference" / "reanalysis.yaml"
DEPLOYMENT = Path(__file__).parents[1] / "data" / "deployment"


def results() -> tuple[EpisodeResult, ...]:
    return run_reference_set(SET, DEPLOYMENT)


def result_named(fragment: str) -> EpisodeResult:
    return next(result for result in results() if fragment in result.episode.name)


@pytest.mark.parametrize("result", results(), ids=lambda r: r.episode.name)
def test_reference_episode(result: EpisodeResult) -> None:
    assert result.reproduced, (
        f"{result.episode.name}: {'; '.join(result.mismatches)}\n"
        f"  got: {result.summarise()}"
    )


# --------------------------------------------------------------------------
# The two claims the whole module exists to keep apart
# --------------------------------------------------------------------------


def test_the_same_outcome_from_two_causes_goes_to_two_different_queues():
    """Both episodes end with a settled mechanism unsettled. One is the field
    learning something and reaches a clinician; the other is us editing a
    threshold and reaches validation. If these ever converge, the system has
    started counting its own corrections as discoveries."""
    recuration = result_named("A curation removes")
    rule_edit = result_named("Our own rule edit")

    assert recuration.kind is rule_edit.kind is DeltaKind.MECHANISM_LOST
    assert recuration.decision.queue is ReviewQueue.CLINICAL_SIGNOUT
    assert rule_edit.decision.queue is ReviewQueue.VALIDATION
    assert rule_edit.decision.rule_change_caveat
    assert not recuration.decision.rule_change_caveat


def test_our_own_change_is_less_urgent_than_the_fields():
    """Lower urgency for the same outcome, and the reason is in the reference
    file: nothing has been learned about the patient, so nothing is
    time-critical — what is needed is somebody confirming the new rule."""
    recuration = result_named("A curation removes")
    rule_edit = result_named("Our own rule edit")

    assert rule_edit.decision.urgency.rank > recuration.decision.urgency.rank


def test_an_axis_that_moved_without_mattering_is_not_blamed():
    """The episode that makes attribution worth its cost. Two releases in one
    window; the constraint moved without crossing anything a rule reads, and
    guessing would have named it."""
    result = result_named("one of them matters")

    assert result.causal_axes == (DriftAxis.GENE_CURATION,)
    assert DriftAxis.POPULATION_FREQUENCY in result.contributing_axes


def test_a_quiet_release_reaches_nobody():
    """The commonest week, and the one that decides whether the output is still
    being read a year in."""
    result = result_named("nothing changes")

    assert result.kind is DeltaKind.NONE
    assert result.decision.queue is ReviewQueue.NONE
    assert result.decision.urgency is Urgency.SILENT


# --------------------------------------------------------------------------
# The set itself
# --------------------------------------------------------------------------


def test_the_set_covers_more_than_one_causal_pattern():
    """A set where every episode has one cause would pass with a much dumber
    attributor than this one."""
    patterns = {result.pattern for result in results() if result.pattern}

    assert len(patterns) > 1


def test_an_episode_that_expects_nothing_is_refused():
    """A case with no expectations passes whatever happens, which is worse than
    having no case: it makes the count look like coverage."""
    variant_spec, baseline, _ = load_episodes(SET)
    catalogue = SourceCatalogue.load(DEPLOYMENT / "catalogue.yaml")
    watched = watched_from(variant_spec)
    engine = RepairbenchEngine(catalogue, {watched.key: watched})

    silent = Episode(name="says nothing", note="", moves={"rules": "v2"}, expect={})
    result = run_episode(silent, engine, catalogue, baseline, watched.key)

    assert not result.reproduced
    assert "cannot fail" in result.mismatches[0]


def test_every_shipped_episode_states_something():
    assert all(result.episode.expect for result in results())


def test_a_set_without_a_baseline_is_refused(tmp_path: Path):
    path = tmp_path / "set.yaml"
    path.write_text("variant: {}\nepisodes: []\n")

    with pytest.raises(ReferenceSetError, match="baseline"):
        load_episodes(path)


def test_every_episode_carries_the_story_it_stands_for():
    """The note is the part a reviewer reads. An episode without one is a
    coordinate change nobody can evaluate."""
    for result in results():
        assert len(result.episode.note) > 80
