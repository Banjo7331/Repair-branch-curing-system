"""What an operator can see: metrics, health, and structured logs.

A scheduled reanalysis is invisible by default. It runs at three in the morning,
finishes, and exits — and the failure that matters is not an error, it is the job
quietly not running at all. Nothing in an error rate detects that. The age of the
last successful run does.

So the most useful thing in this module is one gauge, and the alert that reads it:

    time() - repairbench_last_run_timestamp_seconds > 172800

The Prometheus text exposition format is emitted directly rather than through a
client library. The format is stable and documented, the arithmetic is a dozen
lines, and a package whose whole argument is reproducibility is better off
without a dependency it does not need. Swapping in ``prometheus_client`` behind
the same interface is a one-file change if histograms ever justify it.
"""

from __future__ import annotations

import json
import logging
import sys
import threading
import time
from dataclasses import dataclass, field
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any

from repairbench.reanalysis.ledger import ReanalysisReport


@dataclass
class Metrics:
    """A tiny registry, exposing what an operator actually asks about a run:
    did it happen, what did it cost, and is anything urgent going unread."""

    runs_total: dict[str, float] = field(default_factory=dict)
    events_total: dict[tuple[str, str, str], float] = field(default_factory=dict)
    rule_evaluations_total: float = 0.0
    run_seconds: dict[str, float] = field(default_factory=dict)
    last_run_unix: float = 0.0
    _lock: threading.Lock = field(default_factory=threading.Lock, repr=False)

    def run_completed(self, report: ReanalysisReport, elapsed_seconds: float) -> None:
        with self._lock:
            outcome = "raised" if report.needs_a_human else "quiet"
            self.runs_total[outcome] = self.runs_total.get(outcome, 0.0) + 1
            self.rule_evaluations_total += report.rule_evaluations
            self.run_seconds[report.case_id] = elapsed_seconds
            self.last_run_unix = time.time()
            for event in report.events:
                key = (
                    event.attribution.delta.kind.value,
                    event.queue.value,
                    event.urgency.value,
                )
                self.events_total[key] = self.events_total.get(key, 0.0) + 1

    def expose(self) -> str:
        """Render the Prometheus text exposition format."""
        with self._lock:
            lines: list[str] = []

            _write(lines, "repairbench_runs_total", "counter",
                   "Reanalysis runs completed, by whether anything needed a human.")
            for outcome in sorted(self.runs_total):
                value = self.runs_total[outcome]
                lines.append(f'repairbench_runs_total{{outcome="{outcome}"}} {value:g}')

            _write(lines, "repairbench_rule_evaluations_total", "counter",
                   "Rule evaluations, including counterfactual probes. The cost of attribution.")
            lines.append(f"repairbench_rule_evaluations_total {self.rule_evaluations_total:g}")

            _write(lines, "repairbench_events_total", "counter",
                   "Drift events surfaced, by transition, queue and urgency.")
            for kind, queue, urgency in sorted(self.events_total):
                value = self.events_total[(kind, queue, urgency)]
                lines.append(
                    f'repairbench_events_total{{kind="{kind}",queue="{queue}",'
                    f'urgency="{urgency}"}} {value:g}'
                )

            _write(lines, "repairbench_run_duration_seconds", "gauge",
                   "Wall time of the most recent run, by case.")
            for case_id in sorted(self.run_seconds):
                lines.append(
                    f'repairbench_run_duration_seconds{{case="{case_id}"}} '
                    f"{self.run_seconds[case_id]:g}"
                )

            # The single most useful alerting signal this system has. A scheduled
            # job that stops running is invisible unless something measures its
            # absence, so alert on this gauge's age rather than on failures.
            _write(lines, "repairbench_last_run_timestamp_seconds", "gauge",
                   "Unix time of the last completed run. Alert on its age, not on errors.")
            lines.append(f"repairbench_last_run_timestamp_seconds {self.last_run_unix:g}")

            return "\n".join(lines) + "\n"


class NoMetrics:
    """The default. Keeps the use case free of a hard dependency on being watched."""

    def run_completed(self, report: ReanalysisReport, elapsed_seconds: float) -> None:
        """Deliberately empty."""


def _write(lines: list[str], name: str, kind: str, help_text: str) -> None:
    lines.append(f"# HELP {name} {help_text}")
    lines.append(f"# TYPE {name} {kind}")


def configure_logging(*, json_output: bool = True, level: int = logging.INFO) -> logging.Logger:
    """Structured logs, JSON by default because a log pipeline has to index them."""
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(_JsonFormatter() if json_output else logging.Formatter("%(message)s"))
    logger = logging.getLogger("repairbench")
    logger.handlers = [handler]
    logger.setLevel(level)
    logger.propagate = False
    return logger


class _JsonFormatter(logging.Formatter):
    def format(self, record: logging.LogRecord) -> str:
        payload: dict[str, Any] = {
            "time": self.formatTime(record, "%Y-%m-%dT%H:%M:%S%z"),
            "level": record.levelname.lower(),
            "message": record.getMessage(),
        }
        payload.update(getattr(record, "fields", {}))
        return json.dumps(payload, sort_keys=True)


@dataclass(frozen=True, slots=True)
class Health:
    """What ``/health`` reports.

    ``last_run`` is included on purpose: a service answering "ok" while never
    having run is exactly the failure the endpoint should expose rather than
    hide.
    """

    version: str
    metrics: Metrics

    def as_json(self) -> str:
        return json.dumps(
            {
                "status": "ok",
                "version": self.version,
                "last_run_unix": self.metrics.last_run_unix,
            },
            sort_keys=True,
        )


def serve(addr: str, version: str, metrics: Metrics, logger: logging.Logger) -> None:
    """Expose /health and /metrics until interrupted.

    There is no API here beyond those two, and that is not an omission. The
    reanalysis itself is a one-shot command that exits — cron owns the schedule,
    the process owns one run — and this server exists only so the periods between
    runs are observable.
    """
    health = Health(version=version, metrics=metrics)

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self) -> None:
            if self.path.startswith("/health"):
                self._respond(200, "application/json", health.as_json())
            elif self.path.startswith("/metrics"):
                self._respond(200, "text/plain; version=0.0.4; charset=utf-8", metrics.expose())
            else:
                self._respond(404, "text/plain", "not found\n")

        def _respond(self, status: int, content_type: str, body: str) -> None:
            encoded = body.encode()
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(encoded)))
            self.end_headers()
            self.wfile.write(encoded)

        def log_message(self, format: str, *args: Any) -> None:
            """Route the server's own logging through ours rather than to stderr raw."""
            logger.debug("http", extra={"fields": {"request": format % args}})

    host, _, port = addr.rpartition(":")
    server = ThreadingHTTPServer((host or "0.0.0.0", int(port)), Handler)
    logger.info("listening", extra={"fields": {"addr": addr, "version": version}})
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        logger.info("shutting down")
    finally:
        server.server_close()
