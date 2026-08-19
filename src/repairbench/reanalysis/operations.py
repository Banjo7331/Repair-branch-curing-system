"""Registering a case and running one comparison, for whatever asks.

These two things used to live in ``cli.py``, which was fine while the command
line was the only way in. It stopped being fine the moment a browser could do
them too: two copies of "what registering a case means" is a guarantee that the
terminal and the web page drift apart, and the one that drifts is always the
one nobody is testing.

So the operations live here and both entry points call them. What stays in
``cli.py`` is argument parsing and printing; what stays in ``webapp.py`` is
routing and HTML. Neither owns a rule about what a case *is*.

The one invariant worth stating, because a web button makes it easy to lose:
**neither of these edits a conclusion.** ``register`` records which variants a
case is watching. ``run`` performs exactly the comparison a scheduled process
would, loading the files the world names, and writes what it found. There is no
path here that lets a caller assert a mechanism, an urgency or a queue — those
come out of the rule files or they do not exist.
"""

from __future__ import annotations

import hashlib
import json
import logging
import time
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from repairbench.context.expression import Tissue
from repairbench.model import Consequence, RepairbenchError, Zygosity
from repairbench.observability import Metrics
from repairbench.reanalysis.catalogue import SourceCatalogue
from repairbench.reanalysis.engine import RepairbenchEngine, WatchedVariant
from repairbench.reanalysis.ledger import ReanalysisReport
from repairbench.reanalysis.store import JsonAssessmentStore, JsonCaseRepository, StoreError
from repairbench.reanalysis.usecase import ReanalyseCase
from repairbench.reanalysis.world import DriftAxis, Pin


class OperationError(RepairbenchError):
    """A case cannot be registered or run as asked, and the message says why."""


def parse_variant(spec: str) -> WatchedVariant:
    """``gene:consequence:cds_position[:zygosity]``.

    Zygosity is optional and defaults to unknown rather than to heterozygous,
    because guessing it would silently offer every modality that needs an intact
    allele. Left out, it produces a caveat instead.
    """
    parts = [part.strip() for part in spec.strip().split(":")]
    if len(parts) < 3:
        raise OperationError(f"{spec!r}: expected gene:consequence:cds_position[:zygosity]")
    gene, consequence, position = parts[0], parts[1], parts[2]
    if not position.lstrip("-").isdigit():
        raise OperationError(f"{spec!r}: {position!r} is not a CDS position")
    try:
        parsed_consequence = Consequence(consequence)
    except ValueError:
        known = ", ".join(sorted(item.value for item in Consequence))
        raise OperationError(
            f"{spec!r}: unknown consequence {consequence!r}. Known: {known}"
        ) from None
    try:
        zygosity = Zygosity(parts[3]) if len(parts) > 3 and parts[3] else Zygosity.UNKNOWN
    except ValueError:
        known = ", ".join(sorted(item.value for item in Zygosity))
        raise OperationError(f"{spec!r}: unknown zygosity {parts[3]!r}. Known: {known}") from None
    return WatchedVariant(
        key=f"{gene}-c{position}",
        gene=gene,
        consequence=parsed_consequence,
        zygosity=zygosity,
        cds_position=int(position),
    )


@dataclass(frozen=True, slots=True)
class CaseFiles:
    """Where the parts of a case that are not the ledger are kept.

    The ledger holds what the case is watching; these two sidecars hold what the
    watched variants *are* and which tissue the disease affects. Both are facts
    about the patient rather than about any release, which is why they sit
    beside the ledger rather than inside a world.
    """

    state: Path

    @property
    def directory(self) -> Path:
        return self.state / "cases"

    def variants_path(self, case_id: str) -> Path:
        return self.directory / f"{case_id}.variants.json"

    def tissue_path(self, case_id: str) -> Path:
        return self.directory / f"{case_id}.tissue"

    def write_variants(self, case_id: str, variants: list[WatchedVariant]) -> None:
        self.variants_path(case_id).write_text(
            json.dumps(
                [
                    {
                        "key": variant.key,
                        "gene": variant.gene,
                        "consequence": variant.consequence.value,
                        "cds_position": variant.cds_position,
                        "zygosity": variant.zygosity.value,
                    }
                    for variant in variants
                ],
                indent=2,
            )
        )

    def read_variants(self, case_id: str) -> dict[str, WatchedVariant]:
        path = self.variants_path(case_id)
        if not path.exists():
            raise OperationError(
                f"{case_id} has a ledger but no variant file. A run will not invent the "
                "list of variants it is meant to be watching — register the case again."
            )
        return {
            entry["key"]: WatchedVariant(
                key=entry["key"],
                gene=entry["gene"],
                consequence=Consequence(entry["consequence"]),
                zygosity=Zygosity(entry["zygosity"]),
                cds_position=entry["cds_position"],
            )
            for entry in json.loads(path.read_text())
        }

    def write_tissue(self, case_id: str, tissue: str) -> None:
        self.tissue_path(case_id).write_text(tissue)

    def read_tissue(self, case_id: str) -> str:
        path = self.tissue_path(case_id)
        return path.read_text().strip() if path.exists() else ""


@dataclass(frozen=True, slots=True)
class Registration:
    """What registering produced, for whoever wants to report it."""

    case_id: str
    variants: tuple[WatchedVariant, ...]
    tissue: str
    #: The world the case was assessed against at registration, when a
    #: catalogue was supplied — or why it could not be.
    baseline: str = ""

    @property
    def caveat(self) -> str:
        """The one thing a caller should say out loud if it applies."""
        if self.tissue:
            return ""
        return (
            "no tissue given: nothing will be checked against where these genes are "
            "switched on, and the modality lists will say so"
        )


def register(
    state: Path,
    case_id: str,
    variants: list[WatchedVariant],
    *,
    phenotype: str = "unrecorded",
    tissue: str = "",
    overwrite: bool = False,
    catalogue_path: Path | None = None,
    logger: logging.Logger | None = None,
) -> Registration:
    """Record which variants a case is watching, and what we conclude today.

    Refuses an existing case unless told otherwise. Re-registering silently
    would replace the watched list while keeping the ledger's history, so a case
    could end up carrying events about variants it is no longer watching — which
    reads as a system that lost track rather than as one that was told to change
    its mind.

    Given a catalogue, registration also **lays down the baseline**: it assesses
    every watched variant against today's releases and stores the result. That
    is not an optimisation, it is what "watching" has to mean. Without it the
    ledger has nothing to compare against, so the *first* run after registration
    can never report drift however much moved in between — register a case in
    January, run it in April, and April's answer becomes the baseline while
    February's curation change goes unnoticed. It took a browser to make that
    visible: from a command line the first run is one of many, and from a page
    it is the button somebody presses expecting an answer.
    """
    case_id = case_id.strip()
    if not case_id:
        raise OperationError("a case needs an identifier")
    if "/" in case_id or "\\" in case_id or case_id.startswith("."):
        # The identifier becomes a file name. A traversal here would write a
        # ledger outside the state directory, which for a browser-facing server
        # is the difference between a form and a way in.
        raise OperationError(f"{case_id!r}: a case identifier cannot contain a path")
    if not variants:
        raise OperationError("a case with no variants would be watching nothing")

    cases = JsonCaseRepository(state)
    if not overwrite and case_id in cases.case_ids():
        raise OperationError(
            f"{case_id} is already registered. Registering it again would replace the "
            "variants it watches while keeping its history, so the case could carry "
            "events about variants it no longer watches."
        )

    pin = Pin(
        axis=DriftAxis.PHENOTYPE,
        version=phenotype or "unrecorded",
        digest=hashlib.sha256((phenotype or "unrecorded").encode()).hexdigest(),
    )
    cases.register(case_id, [variant.key for variant in variants], pin)
    files = CaseFiles(Path(state))
    files.write_variants(case_id, variants)
    files.write_tissue(case_id, tissue)

    baseline = ""
    if catalogue_path is not None:
        # One ordinary run. It reports nothing by construction — every variant
        # is being seen for the first time — and what it leaves behind is the
        # assessment the next run will compare against.
        try:
            report = run(state, catalogue_path, case_id, logger or logging.getLogger("repairbench"))
            baseline = report.candidate.describe()
        except (RepairbenchError, OSError, ValueError) as failure:
            # The case stays registered. A baseline that could not be taken is
            # worth reporting and is not worth discarding the registration over:
            # the variants are recorded, and the next run establishes it.
            baseline = f"no baseline taken: {failure}"

    return Registration(
        case_id=case_id, variants=tuple(variants), tissue=tissue, baseline=baseline
    )


def run(
    state: Path,
    catalogue_path: Path,
    case_id: str,
    logger: logging.Logger,
    *,
    tissue: str = "",
) -> ReanalysisReport:
    """One comparison: today against whatever this case was last compared with.

    Exactly what a scheduled process does, with the same inputs — which is the
    property that makes a "run now" button safe to offer. It cannot produce a
    conclusion a cron run would not have produced, because it *is* the cron run,
    started by a person instead of by a clock.
    """
    catalogue = SourceCatalogue.load(catalogue_path)
    files = CaseFiles(Path(state))
    cases = JsonCaseRepository(state)
    if case_id not in cases.case_ids():
        raise OperationError(f"no case {case_id!r} in {state}. Register it first.")

    name = tissue or files.read_tissue(case_id)
    engine = RepairbenchEngine(
        catalogue, files.read_variants(case_id), Tissue(name) if name else None
    )
    metrics = Metrics()

    started = time.monotonic()
    report = ReanalyseCase(
        engine,
        engine,
        cases,
        JsonAssessmentStore(state),
        SystemClock(),
        SequentialIds(Path(state)),
        LoggingNotifier(logger),
    ).execute(case_id)
    metrics.run_completed(report, time.monotonic() - started)
    (Path(state) / "metrics.prom").write_text(metrics.expose())
    return report


def run_all(
    state: Path,
    catalogue_path: Path,
    logger: logging.Logger,
) -> tuple[list[ReanalysisReport], dict[str, str]]:
    """Every registered case, and the ones that could not be run.

    Failures are collected rather than raised. One case with a missing variant
    file must not stop the other forty from being examined — and the reason is
    returned rather than logged away, because a caller reporting "39 of 40 ran"
    without saying which one did not is reporting a number nobody can act on.
    """
    reports: list[ReanalysisReport] = []
    failures: dict[str, str] = {}
    for case_id in JsonCaseRepository(state).case_ids():
        try:
            reports.append(run(state, catalogue_path, case_id, logger))
        except (OperationError, StoreError, RepairbenchError, OSError, ValueError) as failure:
            failures[case_id] = str(failure)
    return reports, failures


class SystemClock:
    def now(self) -> datetime:
        return datetime.now(UTC)


class SequentialIds:
    """Event ids that keep counting across invocations, because each run may be
    a separate process and a counter starting at one every time would collide."""

    def __init__(self, state: Path) -> None:
        self._path = Path(state) / "next-event-id"

    def next_id(self) -> str:
        current = int(self._path.read_text()) if self._path.exists() else 0
        current += 1
        self._path.write_text(str(current))
        return f"evt-{current:06d}"


class LoggingNotifier:
    def __init__(self, logger: logging.Logger) -> None:
        self._logger = logger

    def publish(self, report: object) -> None:
        self._logger.warning(
            "reanalysis needs a human",
            extra={"fields": {"case": getattr(report, "case_id", "?")}},
        )


#: The synthetic case the demo seeds. Named so that nobody mistakes it for a
#: patient, and pinned to the fixture gene whose two releases disagree.
DEMO_CASE = "DEMO-NICU-014"
DEMO_VARIANT = "PLUSG:stop_gained:158:heterozygous"


def seed_demo(state: Path, older_catalogue: Path, logger: logging.Logger) -> str:
    """Put a state directory into the one situation worth demonstrating.

    Drift needs two things that a fresh install cannot have at once: a case
    assessed against *last* month's releases, and *this* month's releases to
    compare with. Producing that by hand meant registering under one catalogue,
    restarting the server under another, and running again — three steps in a
    terminal to see something the browser is supposed to be for.

    So this registers the case against the older releases and takes its
    baseline there. Start the server against the current catalogue afterwards
    and the very first button press has something real to report.

    The case identifier says DEMO, the gene is the fixture's, and both are
    synthetic. This seeds a demonstration; it does not simulate a patient.
    """
    if DEMO_CASE in JsonCaseRepository(state).case_ids():
        return DEMO_CASE
    register(
        state,
        DEMO_CASE,
        [parse_variant(DEMO_VARIANT)],
        phenotype="hpo-day-1",
        catalogue_path=older_catalogue,
        logger=logger,
    )
    return DEMO_CASE
