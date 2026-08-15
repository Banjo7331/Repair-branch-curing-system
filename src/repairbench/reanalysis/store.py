"""Persisting a watched case between runs.

A scheduled reanalysis is a process that starts, compares today with last month,
and exits. That only works if last month survived, so this is the smallest thing
that makes the module runnable from cron: JSON files under a state directory,
one per case and one per assessment.

Two decisions worth defending.

**Assessments are written append-only, and the newest is an index rather than an
overwrite.** "What did we know, and when" has to survive being asked years
later, and a store that overwrites cannot answer it. This is a directory of
files rather than a database, so the property is enforced by never rewriting a
file — a real deployment would move to Postgres and enforce it with a schema.

**A stored assessment does not carry the reasoning that produced it.** It has
the mechanism, the confidence and the two modality sets, which is exactly what a
comparison reads. The evidence trail belongs to the run that produced it and is
in that run's report; duplicating it here would make the store the authority on
something it cannot re-derive.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, replace
from datetime import UTC, datetime
from pathlib import Path

from repairbench.modality import Modality
from repairbench.model import Confidence, Mechanism, RepairbenchError
from repairbench.reanalysis.drift import Assessment
from repairbench.reanalysis.ledger import CaseLedger, DriftEvent, EventStatus, LedgerEntry
from repairbench.reanalysis.world import DriftAxis, Pin, World


class StoreError(RepairbenchError):
    """The state directory is missing something it needs, or holds something it
    should not."""


def _world_to_json(world: World) -> list[dict[str, str]]:
    return [
        {"axis": pin.axis.value, "version": pin.version, "digest": pin.digest}
        for pin in world.pins
    ]


def _world_from_json(raw: list[dict[str, str]]) -> World:
    return World.of(
        Pin(axis=DriftAxis(entry["axis"]), version=entry["version"], digest=entry["digest"])
        for entry in raw
    )


class JsonAssessmentStore:
    """Assessments on disk, newest indexed, nothing overwritten."""

    def __init__(self, root: Path) -> None:
        self._root = Path(root)
        self._root.mkdir(parents=True, exist_ok=True)

    def _case_dir(self, case_id: str) -> Path:
        path = self._root / "assessments" / case_id
        path.mkdir(parents=True, exist_ok=True)
        return path

    @staticmethod
    def _slug(variant_key: str) -> str:
        return variant_key.replace(":", "_").replace("/", "_")

    def latest_for(self, case_id: str, variant_key: str) -> Assessment | None:
        index = self._case_dir(case_id) / f"{self._slug(variant_key)}.latest.json"
        if not index.exists():
            return None
        raw = json.loads(index.read_text())
        return Assessment(
            variant_key=raw["variant_key"],
            world=_world_from_json(raw["world"]),
            mechanism=Mechanism(raw["mechanism"]),
            confidence=Confidence(raw["confidence"]),
            indicated=frozenset(Modality(name) for name in raw["indicated"]),
            contraindicated=frozenset(Modality(name) for name in raw["contraindicated"]),
        )

    def record(self, case_id: str, assessment: Assessment) -> None:
        payload = {
            "variant_key": assessment.variant_key,
            "world": _world_to_json(assessment.world),
            "world_digest": assessment.world.digest,
            "mechanism": assessment.mechanism.value,
            "confidence": assessment.confidence.value,
            "indicated": sorted(m.value for m in assessment.indicated),
            "contraindicated": sorted(m.value for m in assessment.contraindicated),
            "recorded_at": datetime.now(UTC).isoformat(),
        }
        directory = self._case_dir(case_id)
        slug = self._slug(assessment.variant_key)

        # The history file is named for the world, so recording the same world
        # twice is idempotent and recording a different one never overwrites.
        history = directory / f"{slug}.{assessment.world.short_digest}.json"
        history.write_text(json.dumps(payload, indent=2, sort_keys=True))
        serialised = json.dumps(payload, indent=2, sort_keys=True)
        (directory / f"{slug}.latest.json").write_text(serialised)

    def history_for(self, case_id: str, variant_key: str) -> list[dict[str, object]]:
        slug = self._slug(variant_key)
        files = sorted(
            path
            for path in self._case_dir(case_id).glob(f"{slug}.*.json")
            if not path.name.endswith(".latest.json")
        )
        return [json.loads(path.read_text()) for path in files]


class JsonCaseRepository:
    """Case ledgers on disk.

    Events are stored as summaries rather than as full attributions. The
    attribution is a claim about a run — which counterfactuals were tried, what
    they showed — and re-reading it from a file a year later would present it as
    a standing fact. What the ledger needs to keep is the fingerprint, so an
    acknowledged transition stays acknowledged, and enough to show a reviewer
    what they signed.
    """

    def __init__(self, root: Path) -> None:
        self._root = Path(root)
        (self._root / "cases").mkdir(parents=True, exist_ok=True)

    def _path(self, case_id: str) -> Path:
        return self._root / "cases" / f"{case_id}.json"

    def get(self, case_id: str) -> CaseLedger:
        path = self._path(case_id)
        if not path.exists():
            raise StoreError(
                f"no case {case_id!r} in {self._root}. Register it first — a reanalysis "
                "run will not invent the list of variants it is meant to be watching."
            )
        raw = json.loads(path.read_text())
        ledger = CaseLedger(
            case_id=raw["case_id"],
            variant_keys=tuple(raw["variant_keys"]),
            phenotype=Pin(
                axis=DriftAxis.PHENOTYPE,
                version=raw["phenotype"]["version"],
                digest=raw["phenotype"]["digest"],
            ),
            last_world=_world_from_json(raw["last_world"]) if raw.get("last_world") else None,
        )
        ledger.events = [_event_from_json(entry) for entry in raw.get("events", [])]
        return ledger

    def save(self, ledger: CaseLedger) -> None:
        payload = {
            "case_id": ledger.case_id,
            "variant_keys": list(ledger.variant_keys),
            "phenotype": {
                "version": ledger.phenotype.version,
                "digest": ledger.phenotype.digest,
            },
            "last_world": _world_to_json(ledger.last_world) if ledger.last_world else None,
            "events": [_event_to_json(event) for event in ledger.events],
        }
        self._path(ledger.case_id).write_text(json.dumps(payload, indent=2, sort_keys=True))

    def register(self, case_id: str, variant_keys: list[str], phenotype: Pin) -> CaseLedger:
        ledger = CaseLedger(
            case_id=case_id, variant_keys=tuple(variant_keys), phenotype=phenotype
        )
        self.save(ledger)
        return ledger

    def case_ids(self) -> list[str]:
        return sorted(path.stem for path in (self._root / "cases").glob("*.json"))


@dataclass(frozen=True, slots=True)
class StoredEvent:
    """An event read back from disk.

    It satisfies ``LedgerEntry`` — the fingerprint and the status a later run
    consults — and carries the sentence a reviewer read. It does not carry an
    attribution, and that absence is the design: the causal claim belonged to
    the run that made it, and a store that handed one back would be asserting
    something it never established.
    """

    event_id: str
    variant_key: str
    fingerprint: str
    status: EventStatus
    summary: str
    superseded_by: str | None = None

    @property
    def attribution(self) -> object:
        """Deliberately unavailable, with the reason attached.

        An ``AttributeError`` here would read as an oversight. This is not one:
        the causal claim belonged to the run that established it, and a store
        that handed one back a year later would be asserting an experiment it
        never ran. The summary the reviewer read is kept; the claim is not.
        """
        raise StoreError(
            f"the attribution behind {self.event_id} was not reconstructed from disk. "
            f"The run that raised it recorded: {self.summary}"
        )

    def acknowledged(self) -> StoredEvent:
        return replace(self, status=EventStatus.ACKNOWLEDGED)

    def superseded(self, by_event_id: str) -> StoredEvent:
        return replace(self, status=EventStatus.SUPERSEDED, superseded_by=by_event_id)


def _event_to_json(event: LedgerEntry) -> dict[str, object]:
    """Serialise either kind of entry.

    A ledger loaded from disk and saved again contains entries that were never
    live in this process. Round-tripping them has to preserve what they carry
    without inventing what they do not — so a stored entry writes back its
    summary, and only a live event writes the transition and routing behind it.

    Both branches are explicit rather than one being the fallback: a third kind
    of entry appearing later should fail here loudly, not be written out as
    whichever shape the `else` happened to be.
    """
    if isinstance(event, DriftEvent):
        delta = event.attribution.delta
        return {
            "event_id": event.event_id,
            "variant_key": event.variant_key,
            "fingerprint": event.fingerprint,
            "kind": delta.kind.value,
            "mechanism_before": delta.before.mechanism.value,
            "mechanism_after": delta.after.mechanism.value,
            "urgency": event.decision.urgency.value,
            "queue": event.decision.queue.value,
            "reason": event.decision.reason,
            "status": event.status.value,
            "raised_at": event.raised_at.isoformat() if event.raised_at else None,
            "superseded_by": event.superseded_by,
        }

    if isinstance(event, StoredEvent):
        return {
            "event_id": event.event_id,
            "variant_key": event.variant_key,
            "fingerprint": event.fingerprint,
            "kind": "restored",
            "mechanism_before": "-",
            "mechanism_after": "-",
            "urgency": "-",
            "queue": "-",
            "reason": event.summary,
            "status": event.status.value,
            "raised_at": None,
            "superseded_by": event.superseded_by,
        }

    raise StoreError(f"cannot serialise a ledger entry of type {type(event).__name__}")


def _event_from_json(raw: dict[str, object]) -> StoredEvent:
    """Rebuild the parts of an event a later run acts on.

    A stored event is a record, not a replayable object: the attribution that
    produced it is deliberately not reconstructed, because re-presenting a
    year-old causal claim as if it had just been established would misrepresent
    its footing. What survives is what the surfacing policy consults — the
    fingerprint and the status — plus what a reviewer reads.
    """
    return StoredEvent(
        event_id=str(raw["event_id"]),
        variant_key=str(raw["variant_key"]),
        fingerprint=str(raw["fingerprint"]),
        status=EventStatus(str(raw["status"])),
        summary=f"{raw['mechanism_before']} → {raw['mechanism_after']}: {raw['reason']}",
        superseded_by=raw.get("superseded_by"),  # type: ignore[arg-type]
    )


