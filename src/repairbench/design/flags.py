"""Design rules that flag a candidate rather than choose between candidates.

M5 and M6 ask "which answer is right"; a designer asks "what is wrong with this
one". The shape differs — a flag attaches to one candidate and carries a
severity rather than arguing for a conclusion — but the language does not, and
that is the point of this module: pegRNA rules and ASO rules are written in the
same predicate dialect as the mechanism rules, loaded by one parser, and pinned
by the digest of their file.

Why these live in files at all, when they look like engineering constants:
almost none of them is. "Avoid four consecutive Ts" is a fact about Pol III
termination. "Keep the PBS near 30 °C" is a fitted recommendation from one
group's screen. "Nick 40 to 90 nucleotides away for PE3" is a range measured on
a handful of loci. Each will be revised, each is a claim somebody should be able
to check against its citation, and none of them belongs in a Python constant
where it is invisible.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from enum import StrEnum
from pathlib import Path
from typing import Any

import yaml

from repairbench.ruleset import RulesetError, evaluate_predicate


class Severity(StrEnum):
    """How much a flag matters, worst first."""

    BLOCKING = "blocking"
    """The candidate should not be ordered. Not a preference — a reason it
    cannot work as written."""

    CAUTION = "caution"
    """A published reason to expect this one to underperform, or a hazard worth
    a bench check before committing to it."""

    NOTE = "note"
    """Worth knowing, and not a reason to avoid the candidate."""

    @property
    def rank(self) -> int:
        return list(Severity).index(self)


@dataclass(frozen=True, slots=True)
class FlatFeatures:
    """A candidate flattened into the names its rule file uses.

    Same contract as the mechanism feature record and for the same reason: an
    unknown name raises rather than evaluating false, so a typo in a rule file
    fails the run instead of producing a rule that silently never fires.
    """

    values: dict[str, Any]

    def get(self, name: str) -> Any:
        if name not in self.values:
            raise KeyError(name)
        return self.values[name]

    def known_names(self) -> list[str]:
        return sorted(self.values)


@dataclass(frozen=True, slots=True)
class FlagRule:
    id: str
    severity: Severity
    when: dict[str, Any]
    because: str
    citation: str = ""

    def fires(self, features: FlatFeatures) -> bool:
        return evaluate_predicate(self.when, features, self.id)


@dataclass(frozen=True, slots=True)
class Flag:
    """One rule's verdict on one candidate."""

    rule_id: str
    severity: Severity
    because: str
    citation: str = ""

    def describe(self) -> str:
        return f"[{self.severity}] {self.rule_id}: {self.because}"


@dataclass(frozen=True, slots=True)
class FlagRuleset:
    """Rules, thresholds and the digest of the file both came from."""

    version: str
    description: str
    rules: tuple[FlagRule, ...]
    digest: str
    thresholds: dict[str, Any] = field(default_factory=dict)
    #: Whatever else the file declares — chemistries, nucleases — left as read.
    #: The caller knows what its own file carries; this loader does not need to.
    extra: dict[str, Any] = field(default_factory=dict)

    @property
    def pin(self) -> str:
        return f"{self.version}@{self.digest[:12]}"

    def threshold(self, name: str, default: Any = None) -> Any:
        """A named number from the file, with the default declared at the call site.

        Reading it by name rather than into a dataclass keeps each designer's
        thresholds in its own file without a parallel class to maintain, and an
        absent one falls back to a default the caller states out loud.
        """
        return self.thresholds.get(name, default)

    def raise_flags(self, features: FlatFeatures) -> tuple[Flag, ...]:
        return tuple(
            Flag(rule.id, rule.severity, rule.because, rule.citation)
            for rule in self.rules
            if rule.fires(features)
        )


def sort_weight(severity: Severity | None) -> int:
    """Where a candidate sorts: unflagged first, blocking last.

    A separate function from ``Severity.rank`` because the two orders are
    opposite and conflating them is an easy, silent mistake — ``rank`` counts
    down from the worst so that "take the minimum" means "take the worst",
    while a report wants the cleanest candidate at the top. A candidate with no
    flags is not a severity at all, and it sorts ahead of every one that is.
    """
    if severity is None:
        return 0
    return len(Severity) - severity.rank


def worst_of(flags: tuple[Flag, ...]) -> Severity | None:
    """The severity a candidate carries overall.

    The worst one, never an average. Reasons a design will not work do not
    cancel out against reasons it might.
    """
    if not flags:
        return None
    return min((flag.severity for flag in flags), key=lambda severity: severity.rank)


def load_flag_rules(path: str | Path) -> FlagRuleset:
    """Read and validate a design rule file."""
    raw = Path(path).read_bytes()
    document = yaml.safe_load(raw)

    if not isinstance(document, dict):
        raise RulesetError(f"{path}: design rule file must be a mapping")
    for required in ("version", "rules"):
        if required not in document:
            raise RulesetError(f"{path}: design rule file has no {required!r}")

    rules: list[FlagRule] = []
    seen: set[str] = set()
    for index, entry in enumerate(document["rules"], start=1):
        if not isinstance(entry, dict):
            raise RulesetError(f"{path}: rule {index} is not a mapping")
        for required in ("id", "severity", "when", "because"):
            if required not in entry:
                raise RulesetError(f"{path}: rule {index} has no {required!r}")
        try:
            severity = Severity(entry["severity"])
        except ValueError as error:
            raise RulesetError(
                f"{path}: rule {entry['id']!r} has unknown severity {entry['severity']!r}; "
                f"known: {', '.join(level.value for level in Severity)}"
            ) from error
        if entry["id"] in seen:
            raise RulesetError(f"{path}: duplicate rule id {entry['id']!r}")
        seen.add(str(entry["id"]))
        rules.append(
            FlagRule(
                id=str(entry["id"]),
                severity=severity,
                when=entry["when"],
                because=" ".join(str(entry["because"]).split()),
                citation=str(entry.get("citation", "")),
            )
        )

    if not rules:
        raise RulesetError(f"{path}: design rule file declares no rules")

    return FlagRuleset(
        version=str(document["version"]),
        description=str(document.get("description", "")),
        rules=tuple(rules),
        digest=hashlib.sha256(raw).hexdigest(),
        thresholds=dict(document.get("thresholds") or {}),
        extra={
            key: value
            for key, value in document.items()
            if key not in {"version", "description", "rules", "thresholds"}
        },
    )
