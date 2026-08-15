"""Which releases exist on disk, and where.

This is what lets the engine honour a world rather than merely accept one. To
re-assess a variant as it would have been assessed last March, the March files
have to still be here — so the catalogue is a map from ``(axis, version)`` to a
path, and every entry is digested when it is read.

The refusal is the point. Asked for a version it does not hold, the catalogue
raises. "I cannot reproduce that world" is the honest answer, and it is what
makes an attribution trustworthy: a counterfactual that silently fell back to
today's files would produce a confident causal claim about an experiment that
was never run.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from repairbench.context.source import ContextError, Source
from repairbench.reanalysis.world import DriftAxis, Pin


@dataclass(frozen=True, slots=True)
class Release:
    """One version of one axis, as a file on disk."""

    axis: DriftAxis
    version: str
    path: Path
    #: Newer releases sort later. The catalogue keeps the order the file gives
    #: rather than guessing from version strings, because release labels are not
    #: comparable across sources — "2026-07" and "v4.1" do not share a scheme.
    ordinal: int


class SourceCatalogue:
    """Every release this deployment can still reproduce."""

    def __init__(self, releases: list[Release]) -> None:
        self._by_axis: dict[DriftAxis, list[Release]] = {}
        for release in releases:
            self._by_axis.setdefault(release.axis, []).append(release)
        for axis_releases in self._by_axis.values():
            axis_releases.sort(key=lambda release: release.ordinal)

        missing = [
            axis.value
            for axis in DriftAxis
            if not axis.is_case_scoped and axis not in self._by_axis
        ]
        if missing:
            raise ContextError(
                f"the catalogue has no releases for: {', '.join(missing)}. "
                "A world cannot be assembled without every axis, and a run that "
                "assembled one anyway would be unattributable by construction."
            )

    @classmethod
    def load(cls, path: str | Path) -> SourceCatalogue:
        """Read a catalogue file.

        The format is a list per axis, oldest first::

            gene_curation:
              - {version: "2025-10", path: data/clingen-2025-10.tsv}
              - {version: "2026-01", path: data/clingen-2026-01.tsv}
        """
        path = Path(path)
        document = yaml.safe_load(path.read_text()) or {}
        base = path.parent

        releases: list[Release] = []
        for axis_name, entries in document.items():
            try:
                axis = DriftAxis(axis_name)
            except ValueError as error:
                raise ContextError(f"{path.name}: {axis_name!r} is not a drift axis") from error
            if axis.is_case_scoped:
                raise ContextError(
                    f"{path.name}: {axis_name} moves per patient, not per release, and "
                    "belongs on the case rather than in the catalogue"
                )
            for ordinal, entry in enumerate(entries or []):
                file_path = (base / entry["path"]).resolve()
                if not file_path.exists():
                    raise ContextError(
                        f"{path.name}: {axis_name}@{entry['version']} points at "
                        f"{file_path}, which does not exist"
                    )
                releases.append(
                    Release(
                        axis=axis,
                        version=str(entry["version"]),
                        path=file_path,
                        ordinal=ordinal,
                    )
                )
        return cls(releases)

    def path_for(self, axis: DriftAxis, version: str) -> Path:
        for release in self._by_axis.get(axis, []):
            if release.version == version:
                return release.path
        available = ", ".join(r.version for r in self._by_axis.get(axis, [])) or "none"
        raise ContextError(
            f"cannot reproduce {axis.value}@{version}: this deployment holds {available}. "
            "Falling back to the current release would make every counterfactual below "
            "a claim about an experiment that was never run."
        )

    def pin_for(self, axis: DriftAxis, version: str) -> Pin:
        path = self.path_for(axis, version)
        source = Source.of(axis.value, path, version)
        return Pin(axis=axis, version=version, digest=source.digest)

    def latest_pins(self) -> list[Pin]:
        """The newest release of every non-case-scoped axis, digested."""
        return [
            self.pin_for(axis, self._by_axis[axis][-1].version)
            for axis in DriftAxis
            if not axis.is_case_scoped
        ]

    def versions(self, axis: DriftAxis) -> list[str]:
        return [release.version for release in self._by_axis.get(axis, [])]
