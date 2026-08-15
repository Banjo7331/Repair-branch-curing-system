"""How well would this edit work? This package does not know, and says so.

Editing efficiency at a given site is what BE-Hive, DeepBaseEditor and their
successors predict, from large screens of measured outcomes. None of them is
attached here: no weights, no inference, no way to run one.

The interesting design question is what to do about that. Three answers were
available and two of them are worse than nothing:

* Invent a heuristic — GC content, position in the window, a scoring formula
  that looks like the real thing. It would produce a ranked list, and a ranked
  list is read as knowledge. Somebody would order the top guide.
* Say nothing and return candidates in whatever order they were found. The
  order still reads as a ranking, because lists do.
* Make the absence a first-class object that every report prints.

The third is what this module is. ``NoModelAttached`` is not a stub to be filled
in later — it is the correct implementation of "nobody has scored these", and it
refuses in the same shape a real model would answer, so attaching one is a
constructor argument rather than a rewrite.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from repairbench.design.candidate import EditCandidate


@dataclass(frozen=True, slots=True)
class EfficiencyScore:
    """One model's estimate, and the model that made it.

    The model name is not decoration. Two candidates scored by different models
    are not comparable, for the same reason two mechanism calls made under
    different rule files are not.
    """

    value: float
    model: str

    def __post_init__(self) -> None:
        if not 0.0 <= self.value <= 1.0:
            raise ValueError(f"{self.model} returned {self.value}, which is not a fraction")


class EfficiencyModel(Protocol):
    """Whatever can estimate how often an edit is made.

    ``score`` returns ``None`` when the model declines this candidate — an
    editor it was not trained on, a context outside its domain — which is a
    different answer from a low score and is kept different.
    """

    @property
    def name(self) -> str: ...

    @property
    def availability(self) -> str:
        """One line a report prints about what ranking is, or is not, on offer."""
        ...

    def score(self, candidate: EditCandidate) -> EfficiencyScore | None: ...


@dataclass(frozen=True, slots=True)
class NoModelAttached:
    """The default, and an honest answer rather than a placeholder."""

    @property
    def name(self) -> str:
        return "none"

    @property
    def availability(self) -> str:
        return (
            "no efficiency model is attached, so the candidates below are not ranked by how "
            "well they would work — only grouped by whether the editor can reach another base "
            "in the same window. Which of these actually edits, and how often, is what a "
            "screen-trained model answers and what this package does not"
        )

    def score(self, candidate: EditCandidate) -> EfficiencyScore | None:
        return None


def ordered(
    candidates: tuple[EditCandidate, ...], model: EfficiencyModel
) -> tuple[EditCandidate, ...]:
    """Order candidates by whatever the model can say, then by stated criteria.

    With no model, the criteria are declared rather than implied: fewest
    bystanders first, then genomic position, which is a stable order and not a
    claim about quality. Sorting silently by something quality-shaped — window
    centrality, GC content — would be the heuristic this module exists to refuse.
    """
    scored = [(candidate, model.score(candidate)) for candidate in candidates]
    return tuple(
        candidate
        for candidate, _ in sorted(
            scored,
            key=lambda pair: (
                -(pair[1].value if pair[1] else 0.0),
                len(pair[0].bystanders),
                pair[0].span[0],
                pair[0].editor.id,
            ),
        )
    )
