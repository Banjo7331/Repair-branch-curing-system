"""The composition root: running the whole project as of a given world.

Every other module has been a piece. This is where they meet, and the shape of
it is the reason the two halves of this project were merged into one package.

To assess a variant in a world, the engine loads *the files that world names* —
the rule files at their pinned versions, the gene context at its pinned release,
the annotation at its pinned version — and runs the mechanism and modality rules
against them. Honouring a world is not a promise here, it is a file path.

Which is what makes the counterfactuals in ``attribution.py`` mean something. A
probe that advances one axis loads one different file and leaves the rest alone,
so "we re-ran it with January's curation and the change did not happen" is
literally what happened.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

from repairbench.annotation.gff import parse_gff3
from repairbench.annotation.store import TranscriptStore
from repairbench.context.expression import Tissue
from repairbench.context.registry import GeneContextRegistry
from repairbench.engine import resolve
from repairbench.features import MechanismQuery, Variant
from repairbench.modality_rules import ModalityRuleset, load_modality_ruleset
from repairbench.model import Consequence, MissenseDistribution, Zygosity
from repairbench.reanalysis.catalogue import SourceCatalogue
from repairbench.reanalysis.drift import Assessment
from repairbench.reanalysis.world import DriftAxis, Pin, World
from repairbench.ruleset import Ruleset, load_ruleset
from repairbench.selector import select
from repairbench.transcript import Transcript


@dataclass(frozen=True, slots=True)
class WatchedVariant:
    """A variant a case is watching, as the case file describes it.

    These facts do not vary with the world — the patient's genotype is not a
    release — so they live on the case rather than being looked up per run.
    """

    key: str
    gene: str
    consequence: Consequence
    zygosity: Zygosity
    #: Set when the transcript is resolved from annotation rather than asserted.
    chromosome: str | None = None
    position: int | None = None
    cds_position: int | None = None
    distribution: MissenseDistribution = field(default_factory=MissenseDistribution)


class RepairbenchEngine:
    """Assesses a variant as of a world, by loading what that world names."""

    def __init__(
        self,
        catalogue: SourceCatalogue,
        variants: dict[str, WatchedVariant],
        tissue: Tissue | None = None,
    ) -> None:
        self._catalogue = catalogue
        self._variants = variants
        #: The tissue the case is about. Case-scoped like the phenotype, so it
        #: sits on the engine for this case rather than on any release.
        self._tissue = tissue

    def latest_global_pins(self) -> list[Pin]:
        """Satisfies ``SnapshotCatalog``: what a run should be assessed against."""
        return self._catalogue.latest_pins()

    def assess(self, variant_key: str, world: World) -> Assessment:
        watched = self._variants.get(variant_key)
        if watched is None:
            raise KeyError(
                f"{variant_key} is not among the variants this case is watching: "
                f"{', '.join(sorted(self._variants)) or 'none'}"
            )

        rules_dir = self._path(world, DriftAxis.RULES)
        mechanism_rules = _rules(rules_dir)
        modality_rules = _modality_rules(rules_dir)
        registry = _context(
            self._path(world, DriftAxis.GENE_CURATION),
            self._path(world, DriftAxis.POPULATION_FREQUENCY),
            rules_dir,
            self._path(world, DriftAxis.EXPRESSION),
        )
        transcript, cds_position = self._place(watched, world)
        sourced = registry.gene(watched.gene, distribution=watched.distribution)

        query = MechanismQuery(
            variant=Variant(
                gene=watched.gene,
                consequence=watched.consequence,
                cds_position=cds_position,
                zygosity=watched.zygosity,
            ),
            transcript=transcript,
            gene=sourced.gene,
            tissue=self._tissue,
            expression=sourced.expression,
        )
        call = resolve(query, mechanism_rules)
        return Assessment.of(variant_key, world, call, select(call, query, modality_rules))

    def _place(self, watched: WatchedVariant, world: World) -> tuple[Transcript, int]:
        """Resolve the transcript from the annotation release the world names."""
        annotation = _annotation(self._path(world, DriftAxis.ANNOTATION), watched.gene)
        if watched.chromosome is not None and watched.position is not None:
            resolved = annotation.resolve(watched.gene, watched.chromosome, watched.position)
            return resolved.transcript, resolved.cds_position
        if watched.cds_position is None:
            raise ValueError(
                f"{watched.key}: needs either a genomic coordinate or a CDS position"
            )
        record, _ = annotation.preferred_for(watched.gene)
        return record.to_transcript(), watched.cds_position

    def _path(self, world: World, axis: DriftAxis) -> Path:
        """Where this axis's release lives, at the version this world names.

        The ``rules`` axis is a *directory* holding three files that version
        together — the two rule files and our own gene curation. They travel as
        one because the same people edit them in the same review, and pinning
        them separately would let a report cite a combination nobody ever
        reviewed as a whole.
        """
        return self._catalogue.path_for(axis, world.pin_for(axis).version)


@lru_cache(maxsize=32)
def _rules(path: Path) -> Ruleset:
    return load_ruleset(_single(path, "mechanism-*.yaml"))


@lru_cache(maxsize=32)
def _modality_rules(path: Path) -> ModalityRuleset:
    return load_modality_ruleset(_single(path, "modality-*.yaml"))


@lru_cache(maxsize=32)
def _context(
    dosage: Path, constraint: Path, curation_dir: Path, expression: Path
) -> GeneContextRegistry:
    return GeneContextRegistry.load(
        dosage=dosage,
        constraint=constraint,
        local=_single(curation_dir, "*curation*.yaml"),
        expression_matrix=expression,
    )


@lru_cache(maxsize=32)
def _annotation(path: Path, gene: str) -> TranscriptStore:
    return TranscriptStore(parse_gff3(path, genes={gene}))


def _single(directory: Path, pattern: str) -> Path:
    """Exactly one file, or an error naming what was found.

    A rule directory holding two mechanism files is ambiguous, and picking the
    first would make the answer depend on filesystem ordering.
    """
    if directory.is_file():
        return directory
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise ValueError(
            f"{directory} should hold exactly one {pattern}, found "
            f"{[p.name for p in matches] or 'none'}"
        )
    return matches[0]
