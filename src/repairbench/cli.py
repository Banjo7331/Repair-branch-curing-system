"""Command-line entry point.

Three commands, and the shape of them is the argument of the module. ``rules``
prints what the system believes, ``explain`` applies it to one case, and
``reference`` re-runs the published cases the rules are validated against.

There is no command that returns a mechanism without its reasoning, and there is
no command that proposes a therapy. The second is M6's job and this module ends
before it.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import sys
import time
from collections.abc import Callable
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as package_version
from pathlib import Path
from typing import Any

import yaml

from repairbench.annotation.fasta import IndexedFasta
from repairbench.annotation.gff import parse_gff3
from repairbench.annotation.normalise import left_align, verify_reference
from repairbench.annotation.store import TranscriptStore
from repairbench.context.expression import Tissue, median_tpm
from repairbench.context.expression import ingest as expression_ingest
from repairbench.context.genelists import load_gene_lists
from repairbench.context.registry import GeneContextRegistry
from repairbench.context.source import Provenance, Source
from repairbench.design.aso import Exon, tile
from repairbench.design.designer import CorrectionRequest, design
from repairbench.design.editors import load_editors
from repairbench.design.flags import load_flag_rules
from repairbench.design.offtarget import read_casoffinder
from repairbench.design.prime import EditRequest, design_pegrnas
from repairbench.design.report import (
    render_asos,
    render_design,
    render_offtarget,
    render_pegrnas,
    render_plan,
)
from repairbench.design.risk import RiskTier, assess, load_risk_rules
from repairbench.engine import resolve
from repairbench.features import MechanismQuery, SplicePrediction, Variant
from repairbench.modality_rules import load_modality_ruleset
from repairbench.model import (
    Consequence,
    DosageScore,
    Gene,
    Imprinting,
    MissenseDistribution,
    RepairbenchError,
    Zygosity,
)
from repairbench.observability import Metrics, configure_logging, serve
from repairbench.plan import Designers, Locus, load_routing, plan
from repairbench.reanalysis.catalogue import SourceCatalogue
from repairbench.reanalysis.engine import RepairbenchEngine, WatchedVariant
from repairbench.reanalysis.store import JsonAssessmentStore, JsonCaseRepository
from repairbench.reanalysis.usecase import ReanalyseCase
from repairbench.reanalysis.world import DriftAxis, Pin, World
from repairbench.report import render, render_selection
from repairbench.ruleset import Ruleset, load_ruleset
from repairbench.selector import select
from repairbench.transcript import Transcript
from repairbench.vcf import VcfReader

DEFAULT_RULES = Path(__file__).parents[2] / "rules" / "mechanism-v1.yaml"
DEFAULT_MODALITY_RULES = Path(__file__).parents[2] / "rules" / "modality-v1.yaml"
DEFAULT_EDITORS = Path(__file__).parents[2] / "rules" / "editors-v1.yaml"
DEFAULT_RISK_RULES = Path(__file__).parents[2] / "rules" / "offtarget-v1.yaml"
DEFAULT_PRIME_RULES = Path(__file__).parents[2] / "rules" / "prime-v1.yaml"
DEFAULT_ASO_RULES = Path(__file__).parents[2] / "rules" / "aso-v1.yaml"
DEFAULT_ROUTING = Path(__file__).parents[2] / "rules" / "routing-v1.yaml"
REFERENCE_DIR = Path(__file__).parents[2] / "tests" / "reference"
CONTEXT_DIR = Path(__file__).parents[2] / "tests" / "data" / "context"


def version(package: str) -> str:
    try:
        return package_version(package)
    except PackageNotFoundError:
        return "dev"


def _add_design_commands(subcommands: Any) -> None:
    """The M7 half of the command line, split out because one function that
    declares every command is unreadable long before it is unmaintainable."""
    design_command = subcommands.add_parser(
        "design", help="place a base-editing protospacer over one substitution"
    )
    design_command.add_argument("--gene", required=True)
    design_command.add_argument("--at", required=True, help="chromosome:position, 1-based")
    design_command.add_argument("--patient", required=True, help="the base the patient has")
    design_command.add_argument("--wild-type", required=True, help="the base it should read")
    design_command.add_argument("--fasta", type=Path, required=True, help="indexed reference")
    design_command.add_argument("--editors", type=Path, default=DEFAULT_EDITORS)
    design_command.add_argument("--annotation", type=Path, help="GFF3, to place bystanders")
    design_command.add_argument("--zygosity", default="unknown")

    pegrna = subcommands.add_parser(
        "pegrna", help="design pegRNAs for one edit, including insertions and deletions"
    )
    pegrna.add_argument("--gene", required=True)
    pegrna.add_argument("--at", required=True, help="chromosome:position, 1-based")
    pegrna.add_argument("--patient", required=True, help="the patient's allele, VCF-style")
    pegrna.add_argument("--wild-type", required=True, help="the allele it should read")
    pegrna.add_argument("--fasta", type=Path, required=True, help="indexed reference")
    pegrna.add_argument("--prime-rules", type=Path, default=DEFAULT_PRIME_RULES)
    pegrna.add_argument("--limit", type=int, default=5, help="how many to print")

    aso_command = subcommands.add_parser(
        "aso", help="tile antisense oligonucleotides along a target region"
    )
    aso_command.add_argument("--gene", required=True)
    aso_command.add_argument("--at", required=True, help="chromosome:start-end, 1-based inclusive")
    aso_command.add_argument("--fasta", type=Path, required=True, help="indexed reference")
    aso_command.add_argument("--chemistry", required=True, help="an id from the ASO rule file")
    aso_command.add_argument("--aso-rules", type=Path, default=DEFAULT_ASO_RULES)
    aso_command.add_argument(
        "--exon", help="exon bounds as start-end, to place each window against the splice sites"
    )
    aso_command.add_argument("--limit", type=int, default=8, help="how many to print")

    plan_command = subcommands.add_parser(
        "plan", help="one case end to end: mechanism, then modalities, then molecules"
    )
    plan_command.add_argument("case", type=Path, help="YAML case file with a genomic block")
    plan_command.add_argument("--fasta", type=Path, help="indexed reference, to design against")
    plan_command.add_argument("--annotation", type=Path, help="GFF3, for the affected exon")
    plan_command.add_argument("--modality-rules", type=Path, default=DEFAULT_MODALITY_RULES)
    plan_command.add_argument("--editors", type=Path, default=DEFAULT_EDITORS)
    plan_command.add_argument("--prime-rules", type=Path, default=DEFAULT_PRIME_RULES)
    plan_command.add_argument("--aso-rules", type=Path, default=DEFAULT_ASO_RULES)
    plan_command.add_argument("--routing", type=Path, default=DEFAULT_ROUTING)
    plan_command.add_argument("--limit", type=int, default=3, help="molecules to print each")

    offtarget_command = subcommands.add_parser(
        "offtarget", help="rank a Cas-OFFinder hit list by where the hits land"
    )
    offtarget_command.add_argument("hits", type=Path, help="Cas-OFFinder output")
    offtarget_command.add_argument("--guide", help="restrict to one crRNA")
    offtarget_command.add_argument("--risk-rules", type=Path, default=DEFAULT_RISK_RULES)
    offtarget_command.add_argument("--annotation", type=Path, help="GFF3, to place each hit")
    offtarget_command.add_argument("--gene-lists", type=Path, help="symbol/list TSV")
    offtarget_command.add_argument("--expression", type=Path, help="GTEx median TPM matrix")
    offtarget_command.add_argument("--tissue", help="the tissue the therapy is aimed at")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="repairbench", description=__doc__)
    parser.add_argument("--rules", type=Path, default=DEFAULT_RULES, help="path to the rule file")
    subcommands = parser.add_subparsers(dest="command", required=True)

    subcommands.add_parser("rules", help="print the loaded rules and what each one claims")

    explain = subcommands.add_parser("explain", help="determine the mechanism of one case")
    explain.add_argument("case", type=Path, help="YAML case file")
    explain.add_argument("--annotation", type=Path, help="GFF3 to resolve the transcript from")
    explain.add_argument("--fasta", type=Path, help="indexed reference FASTA for normalisation")

    context = subcommands.add_parser(
        "context", help="show gene-level facts and where each one came from"
    )
    context.add_argument("gene", nargs="*", help="restrict to these symbols")
    context.add_argument("--dosage", type=Path, default=CONTEXT_DIR / "clingen_dosage.tsv")
    context.add_argument("--constraint", type=Path, default=CONTEXT_DIR / "gnomad_constraint.tsv")
    context.add_argument("--local", type=Path, default=CONTEXT_DIR / "local_curation.yaml")
    context.add_argument(
        "--expression",
        type=Path,
        default=CONTEXT_DIR / "gtex_median_tpm.tsv",
        help="GTEx-shaped median TPM matrix (GCT preamble stripped)",
    )
    context.add_argument("--tissue", help="report expression in this tissue only")

    annotation = subcommands.add_parser(
        "annotation", help="summarise a GFF3: which transcripts, which are MANE Select"
    )
    annotation.add_argument("gff", type=Path)
    annotation.add_argument("--gene", action="append", help="restrict to these symbols")

    assess = subcommands.add_parser(
        "assess", help="determine the mechanism, then which modalities it admits"
    )
    assess.add_argument("case", type=Path, help="YAML case file")
    assess.add_argument("--modality-rules", type=Path, default=DEFAULT_MODALITY_RULES)

    _add_design_commands(subcommands)

    watch = subcommands.add_parser(
        "watch", help="register a case so scheduled runs know what to re-examine"
    )
    watch.add_argument("case_id")
    watch.add_argument("--state", type=Path, required=True)
    watch.add_argument("--catalogue", type=Path, required=True)
    watch.add_argument("--variant", action="append",
                       help="gene:consequence:cds_position[:zygosity], repeatable")
    watch.add_argument("--vcf", type=Path, help="read the variants from a VCF instead")
    watch.add_argument("--sample", help="which sample in the VCF")
    watch.add_argument("--fasta", type=Path, help="indexed reference, to left-align indels")
    watch.add_argument("--phenotype", default="unrecorded", help="HPO snapshot label")
    watch.add_argument("--tissue", help="the tissue the disease affects, as GTEx names it")

    reanalyse = subcommands.add_parser(
        "reanalyse", help="re-examine a registered case against the current releases"
    )
    reanalyse.add_argument("case_id")
    reanalyse.add_argument("--state", type=Path, required=True)
    reanalyse.add_argument("--catalogue", type=Path, required=True)
    reanalyse.add_argument("--json", action="store_true", help="structured logs")
    reanalyse.add_argument("--tissue", help="override the tissue recorded at registration")

    serve_command = subcommands.add_parser(
        "serve", help="expose /health and /metrics between scheduled runs"
    )
    serve_command.add_argument("--addr", default=":9090")

    reference = subcommands.add_parser(
        "reference", help="re-run the published cases the rules are validated against"
    )
    reference.add_argument("--set", type=Path, default=REFERENCE_DIR / "mechanisms.yaml")
    reference.add_argument(
        "--modalities",
        action="store_true",
        help="run the modality reference set instead of the mechanism one",
    )
    reference.add_argument("--modality-rules", type=Path, default=DEFAULT_MODALITY_RULES)

    args = parser.parse_args(argv)

    try:
        return _dispatch(args, load_ruleset(args.rules))
    except (RepairbenchError, OSError) as error:
        print(f"repairbench: {error}", file=sys.stderr)
        return 1


def _dispatch(args: argparse.Namespace, ruleset: Ruleset) -> int:
    """Route a parsed command line.

    A table rather than a ladder, for the same reason the rule files are tables:
    the set of commands should be readable in one place rather than
    reconstructed from control flow.
    """
    handlers: dict[str, Callable[[], int]] = {
        "rules": lambda: _print_rules(ruleset),
        "context": lambda: _context(args, ruleset.thresholds.expressed_above_tpm),
        "annotation": lambda: _annotation(args.gff, set(args.gene) if args.gene else None),
        "design": lambda: _design(args),
        "pegrna": lambda: _pegrna(args),
        "aso": lambda: _aso(args),
        "plan": lambda: _plan(ruleset, args),
        "offtarget": lambda: _offtarget(args),
        "explain": lambda: _explain(ruleset, args.case, args.annotation, args.fasta),
        "assess": lambda: _assess(ruleset, args.modality_rules, args.case),
        "watch": lambda: _watch(args),
        "reanalyse": lambda: _reanalyse(args),
        "serve": lambda: _serve(args),
        "reference": lambda: (
            _modality_reference(ruleset, args.modality_rules, REFERENCE_DIR / "modalities.yaml")
            if args.modalities
            else _reference(ruleset, args.set)
        ),
    }
    return handlers[args.command]()


def _print_rules(ruleset: Ruleset) -> int:
    print(f"{ruleset.version}  {ruleset.short_digest}")
    print(ruleset.description.strip())
    print()
    for rule in ruleset.rules:
        print(f"{rule.id}  →  {rule.supports}  [{rule.strength}]")
        print(f"    {rule.because}")
        if rule.citation:
            print(f"    — {rule.citation}")
        print()
    print(f"{len(ruleset.rules)} rules")
    return 0


def _explain(
    ruleset: Ruleset,
    case_path: Path,
    annotation_path: Path | None = None,
    fasta_path: Path | None = None,
) -> int:
    case = yaml.safe_load(case_path.read_text())
    query, provenance = build_query_with_provenance(case, annotation_path, fasta_path)
    for line in provenance:
        print(f"# {line}")
    if provenance:
        print()
    print(render(resolve(query, ruleset)))
    return 0


def _context(args: argparse.Namespace, expressed_above_tpm: float) -> int:
    """Print gene facts with a citation on every line.

    The point of the output shape: a fact from ClinGen says ClinGen, and a fact
    we decided says we decided. Two of the four facts the rules read have no
    public table, and this is where that stops being invisible.
    """
    wanted = set(args.gene) if args.gene else None
    registry = GeneContextRegistry.load(
        dosage=args.dosage,
        constraint=args.constraint,
        local=args.local,
        expression_matrix=args.expression,
        genes=wanted,
    )
    tissue = Tissue(args.tissue) if args.tissue else None
    print("sources")
    for source in registry.pins:
        print(f"  {source.pin}  ({source.path.name})")
    print()
    for symbol in registry.genes:
        print(registry.provenance_for(symbol).describe())
        if tissue is not None:
            _print_expression(registry.gene(symbol).expression, tissue, expressed_above_tpm)
        print()
    print(f"{len(registry.genes)} genes")
    return 0


def _split_coordinate(text: str) -> tuple[str, int]:
    chromosome, _, position = text.partition(":")
    if not position.isdigit():
        raise RepairbenchError(f"{text!r}: expected chromosome:position, 1-based")
    return chromosome, int(position)


def _design(args: argparse.Namespace) -> int:
    """Design base-editing candidates for one substitution."""
    chromosome, position = _split_coordinate(args.at)
    catalogue = load_editors(args.editors)

    in_coding_sequence: Callable[[int], bool] | None = None
    if args.annotation:
        store = TranscriptStore(parse_gff3(args.annotation))
        print(f"# annotation  {store.pin}")
        in_coding_sequence = lambda genomic: store.locate(  # noqa: E731
            chromosome, genomic
        ).in_coding_sequence

    with IndexedFasta(args.fasta) as genome:
        outcome = design(
            CorrectionRequest(
                gene=args.gene,
                chromosome=chromosome,
                position=position,
                patient_base=args.patient,
                wild_type_base=args.wild_type,
                zygosity=Zygosity(args.zygosity),
            ),
            genome,
            catalogue,
            coding=in_coding_sequence,
        )
    print(render_design(outcome))
    return 0 if outcome.has_candidates else 2


def _split_span(text: str) -> tuple[str, int, int]:
    chromosome, _, span = text.partition(":")
    start, _, end = span.partition("-")
    if not (start.isdigit() and end.isdigit()):
        raise RepairbenchError(f"{text!r}: expected chromosome:start-end, 1-based inclusive")
    return chromosome, int(start), int(end)


def _pegrna(args: argparse.Namespace) -> int:
    """Design pegRNAs for one edit."""
    chromosome, position = _split_coordinate(args.at)
    rules = load_flag_rules(args.prime_rules)

    with IndexedFasta(args.fasta) as genome:
        outcome = design_pegrnas(
            EditRequest(
                gene=args.gene,
                chromosome=chromosome,
                position=position,
                patient_allele=args.patient,
                wild_type_allele=args.wild_type,
            ),
            genome,
            rules,
        )
    print(render_pegrnas(outcome, limit=args.limit))
    return 0 if outcome.usable else 2


def _aso(args: argparse.Namespace) -> int:
    """Tile antisense oligonucleotides along a region."""
    chromosome, start, end = _split_span(args.at)
    rules = load_flag_rules(args.aso_rules)

    exon = None
    if args.exon:
        exon_start, _, exon_end = args.exon.partition("-")
        exon = Exon(start=int(exon_start), end=int(exon_end))

    with IndexedFasta(args.fasta) as genome:
        target = genome.fetch(chromosome, start, end)

    outcome = tile(
        args.gene,
        chromosome,
        start,
        target,
        rules,
        chemistry=args.chemistry,
        exon=exon,
    )
    print(render_asos(outcome, limit=args.limit))
    return 0 if outcome.usable else 2


def _plan(ruleset: Ruleset, args: argparse.Namespace) -> int:
    """One case from why to what, in a single document.

    The exit code is about the assessment rather than the molecules: a case
    where every modality is ruled out is a finished, useful answer, and a case
    where the mechanism could not be resolved is not.
    """
    case = yaml.safe_load(args.case.read_text())
    store = TranscriptStore(parse_gff3(args.annotation)) if args.annotation else None
    query, provenance = build_query_with_provenance(case, args.annotation, args.fasta)
    for line in provenance:
        print(f"# {line}")
    if provenance:
        print()

    modality_rules = load_modality_ruleset(args.modality_rules)
    call = resolve(query, ruleset)
    selection = select(call, query, modality_rules)

    designers = Designers(
        editors=load_editors(args.editors),
        prime=load_flag_rules(args.prime_rules),
        aso=load_flag_rules(args.aso_rules),
        routing=load_routing(args.routing),
    )
    locus = Locus.from_case(case, store)

    if args.fasta and locus is not None:
        with IndexedFasta(args.fasta) as genome:
            assembled = plan(
                query,
                call,
                selection,
                designers,
                locus=locus,
                sequences=genome,
                mechanism_rules=ruleset,
                modality_rules=modality_rules,
            )
    else:
        assembled = plan(
            query,
            call,
            selection,
            designers,
            locus=locus,
            mechanism_rules=ruleset,
            modality_rules=modality_rules,
        )

    print(render_plan(assembled, limit=args.limit))
    return 2 if selection.is_blocked else 0


def _offtarget(args: argparse.Namespace) -> int:
    """Rank a hit list somebody else's search produced."""
    hits = read_casoffinder(args.hits, guide=args.guide)
    rules = load_risk_rules(args.risk_rules)

    locate = None
    if args.annotation:
        store = TranscriptStore(parse_gff3(args.annotation))
        locate = store.locate

    lists = load_gene_lists(args.gene_lists) if args.gene_lists else None

    expression: dict[str, dict[str, float]] | None = None
    if args.expression:
        collected: dict[str, Provenance] = {}
        expression_ingest(Source.of("expression", args.expression, "unversioned"), collected)
        expression = {
            gene: dict(provenance.facts["expression"].value)
            for gene, provenance in collected.items()
            if "expression" in provenance.facts
        }

    assessment = assess(
        hits,
        rules,
        source_pin=str(args.hits),
        locate=locate,
        lists=lists,
        expression=expression,
        tissue=Tissue(args.tissue) if args.tissue else None,
    )
    print(render_offtarget(assessment))
    return 1 if assessment.worst is RiskTier.PROHIBITIVE else 0


def _print_expression(
    profile: dict[str, float] | None, tissue: Tissue, expressed_above_tpm: float
) -> None:
    """One tissue's number, with the three answers kept apart.

    No profile means no expression release covered this gene; a profile without
    this tissue means the release did not measure it; a number means it did.
    Only the third is evidence, and only the third should read as one.
    """
    value = median_tpm(profile, tissue)
    if value is None and profile is None:
        print(f"      {tissue.name}: no expression source covers this gene")
    elif value is None:
        print(f"      {tissue.name}: not measured in this release")
    else:
        state = "expressed" if value >= expressed_above_tpm else "essentially silent"
        print(f"      {tissue.name} ({tissue.system}): {value:g} TPM — {state}")


def _parse_variant(spec: str) -> WatchedVariant:
    """``gene:consequence:cds_position[:zygosity]``.

    Zygosity is optional and defaults to unknown rather than to heterozygous,
    because guessing it would silently offer every modality that needs an intact
    allele. Left out, it produces a caveat instead.
    """
    parts = spec.split(":")
    if len(parts) < 3:
        raise RepairbenchError(
            f"{spec!r}: expected gene:consequence:cds_position[:zygosity]"
        )
    gene, consequence, position = parts[0], parts[1], int(parts[2])
    zygosity = Zygosity(parts[3]) if len(parts) > 3 else Zygosity.UNKNOWN
    return WatchedVariant(
        key=f"{gene}-c{position}",
        gene=gene,
        consequence=Consequence(consequence),
        zygosity=zygosity,
        cds_position=position,
    )


def _watch(args: argparse.Namespace) -> int:
    """Register a case. A run will not invent the variants it is watching."""
    catalogue = SourceCatalogue.load(args.catalogue)
    if args.vcf:
        variants = _variants_from_vcf(args, catalogue)
    elif args.variant:
        variants = [_parse_variant(spec) for spec in args.variant]
    else:
        raise RepairbenchError("give either --vcf or at least one --variant")
    cases = JsonCaseRepository(args.state)
    phenotype = Pin(
        axis=DriftAxis.PHENOTYPE,
        version=args.phenotype,
        digest=hashlib.sha256(args.phenotype.encode()).hexdigest(),
    )
    cases.register(args.case_id, [v.key for v in variants], phenotype)
    (Path(args.state) / "cases" / f"{args.case_id}.variants.json").write_text(
        json.dumps(
            [
                {
                    "key": v.key,
                    "gene": v.gene,
                    "consequence": v.consequence.value,
                    "cds_position": v.cds_position,
                    "zygosity": v.zygosity.value,
                }
                for v in variants
            ],
            indent=2,
        )
    )
    (Path(args.state) / "cases" / f"{args.case_id}.tissue").write_text(args.tissue or "")
    if not args.tissue:
        print(
            "# no --tissue given: nothing will be checked against where these genes are "
            "switched on, and the modality lists will say so"
        )
    print(f"watching {args.case_id}: {', '.join(v.key for v in variants)}")
    print(f"current world: {World.of([*catalogue.latest_pins(), phenotype]).describe()}")
    return 0


def _variants_from_vcf(
    args: argparse.Namespace, catalogue: SourceCatalogue
) -> list[WatchedVariant]:
    """Turn a patient's VCF into the variants a case will watch.

    Three filters, and each one reports what it dropped rather than doing it
    quietly. A run that silently watches four of a patient's forty variants is
    worse than one that refuses, because nobody would know to ask.
    """
    sequences = None
    if args.fasta:
        sequences = IndexedFasta(args.fasta)
    try:
        observed = VcfReader(
            path=args.vcf, sample=args.sample, expected_assembly="GRCh38"
        ).read(sequences)
    finally:
        if sequences is not None:
            sequences.close()

    if not args.fasta:
        print(
            "# no --fasta given: indels are recorded as the caller wrote them, which may "
            "not match how the same change is written elsewhere"
        )

    annotation_path = catalogue.path_for(
        DriftAxis.ANNOTATION, catalogue.versions(DriftAxis.ANNOTATION)[-1]
    )

    watched: list[WatchedVariant] = []
    uninterpretable = 0
    unplaceable: list[str] = []
    for variant in observed:
        gene, consequence = variant.gene, variant.consequence
        if gene is None or consequence is None:
            uninterpretable += 1
            continue
        store = TranscriptStore(parse_gff3(annotation_path, genes={gene}))
        try:
            resolved = store.resolve(gene, variant.chromosome, variant.position)
        except RepairbenchError as error:
            unplaceable.append(f"{variant.key} ({error})")
            continue
        watched.append(
            WatchedVariant(
                key=f"{gene}-c{resolved.cds_position}",
                gene=gene,
                consequence=consequence,
                zygosity=variant.zygosity,
                cds_position=resolved.cds_position,
            )
        )

    print(f"# {len(observed)} carried alleles read from {args.vcf.name}")
    if uninterpretable:
        print(f"# {uninterpretable} without a gene or consequence annotation — run VEP first")
    for note in unplaceable:
        print(f"# not placed on a transcript: {note}")
    if not watched:
        raise RepairbenchError(
            "nothing to watch: no variant in this VCF could be placed on a coding transcript"
        )
    return watched


def _reanalyse(args: argparse.Namespace) -> int:
    """One run, then exit. Cron owns the schedule; this owns one comparison."""
    logger = configure_logging(json_output=args.json)
    catalogue = SourceCatalogue.load(args.catalogue)
    cases = JsonCaseRepository(args.state)
    variants = {
        entry["key"]: WatchedVariant(
            key=entry["key"],
            gene=entry["gene"],
            consequence=Consequence(entry["consequence"]),
            zygosity=Zygosity(entry["zygosity"]),
            cds_position=entry["cds_position"],
        )
        for entry in json.loads(
            (Path(args.state) / "cases" / f"{args.case_id}.variants.json").read_text()
        )
    }
    tissue_file = Path(args.state) / "cases" / f"{args.case_id}.tissue"
    recorded = tissue_file.read_text().strip() if tissue_file.exists() else ""
    name = args.tissue or recorded
    engine = RepairbenchEngine(catalogue, variants, Tissue(name) if name else None)
    metrics = Metrics()

    started = time.monotonic()
    report = ReanalyseCase(
        engine,
        engine,
        cases,
        JsonAssessmentStore(args.state),
        _SystemClock(),
        _SequentialIds(args.state),
        _LoggingNotifier(logger),
    ).execute(args.case_id)
    metrics.run_completed(report, time.monotonic() - started)

    print(report.headline())
    for event in report.events:
        print(f"  → {event.summary()}")
    (Path(args.state) / "metrics.prom").write_text(metrics.expose())
    return 0


def _serve(args: argparse.Namespace) -> int:
    serve(args.addr, version("repairbench"), Metrics(), configure_logging())
    return 0


class _SystemClock:
    def now(self) -> datetime:
        return datetime.now(UTC)


class _SequentialIds:
    """Event ids that keep counting across invocations, because each run is a
    separate process and a counter starting at one every time would collide."""

    def __init__(self, state: Path) -> None:
        self._path = Path(state) / "next-event-id"

    def next_id(self) -> str:
        current = int(self._path.read_text()) if self._path.exists() else 0
        current += 1
        self._path.write_text(str(current))
        return f"evt-{current:06d}"


class _LoggingNotifier:
    def __init__(self, logger: logging.Logger) -> None:
        self._logger = logger

    def publish(self, report: object) -> None:
        self._logger.warning(
            "reanalysis needs a human",
            extra={"fields": {"case": getattr(report, "case_id", "?")}},
        )


def _annotation(gff_path: Path, genes: set[str] | None) -> int:
    """Show what an annotation file actually contains before trusting it."""
    parsed = parse_gff3(gff_path, genes=genes)
    store = TranscriptStore(parsed)
    print(f"{parsed.pin}  {len(parsed)} coding transcripts across {len(store.genes)} genes")
    print()
    for gene in store.genes:
        record, reason = store.preferred_for(gene)
        print(
            f"{gene:<12} {record.accession:<16} {record.strand}  "
            f"{len(record.cds_blocks)} coding exons, {record.coding_length} nt   [{reason}]"
        )
    return 0


def resolve_from_annotation(
    case: dict[str, Any],
    annotation_path: Path,
    fasta_path: Path | None,
) -> tuple[dict[str, Any], list[str]]:
    """Turn a genomic coordinate into the transcript facts the rules need.

    This is the whole point of the annotation package: until it existed, a case
    file asserted its own exon structure and CDS offset, and the NMD arithmetic
    — which decides between two opposite therapies — ran on numbers nobody had
    checked. Here the structure comes from the annotation, the coordinate is
    verified against the reference, and both sources are named in the output.
    """
    genomic = case["genomic"]
    chromosome, position = str(genomic["chromosome"]), int(genomic["position"])
    reference, alternate = genomic.get("reference"), genomic.get("alternate")
    gene = case["variant"]["gene"]

    notes: list[str] = []
    parsed = parse_gff3(annotation_path, genes={gene})
    store = TranscriptStore(parsed)
    notes.append(f"annotation  {store.pin}")

    if fasta_path and reference and alternate:
        with IndexedFasta(fasta_path) as genome:
            verify_reference(chromosome, position, reference, genome)
            normalised = left_align(chromosome, position, reference, alternate, genome)
        if normalised.shifted_by:
            notes.append(
                f"normalised  {chromosome}:{position} {reference}>{alternate} moved "
                f"{normalised.shifted_by} bases left to {normalised.key}"
            )
        position = normalised.position

    resolved = store.resolve(gene, chromosome, position, accession=case.get("transcript_accession"))
    notes.append(
        f"transcript  {resolved.record.accession} ({resolved.selection_reason}), "
        f"{len(resolved.record.cds_blocks)} coding exons, "
        f"{resolved.record.coding_length} nt; variant at c.{resolved.cds_position}"
    )

    enriched = dict(case)
    enriched["variant"] = {**case["variant"], "cds_position": resolved.cds_position}
    enriched["transcript"] = {
        "accession": resolved.record.accession,
        "exon_lengths": list(resolved.record.coding_exon_lengths),
        "mane_select": resolved.record.mane_select,
    }
    return enriched, notes


def build_query_with_provenance(
    case: dict[str, Any],
    annotation_path: Path | None = None,
    fasta_path: Path | None = None,
) -> tuple[MechanismQuery, list[str]]:
    if annotation_path and "genomic" in case:
        enriched, notes = resolve_from_annotation(case, annotation_path, fasta_path)
        return build_query(enriched), notes
    if "genomic" in case:
        return build_query(case), [
            "no annotation supplied — the transcript structure below is asserted by the "
            "case file rather than read from a reference"
        ]
    return build_query(case), []


def _assess(ruleset: Ruleset, modality_rules: Path, case_path: Path) -> int:
    case = yaml.safe_load(case_path.read_text())
    query = build_query(case)
    call = resolve(query, ruleset)
    print(render(call))
    print()
    print(render_selection(select(call, query, load_modality_ruleset(modality_rules))))
    return 0


def _modality_reference(ruleset: Ruleset, modality_rules: Path, set_path: Path) -> int:
    """Re-run the modality reference set: diseases where a route was in fact taken."""
    modality_ruleset = load_modality_ruleset(modality_rules)
    cases = yaml.safe_load(set_path.read_text())["cases"]
    misses = 0

    for case in cases:
        query = build_query(case)
        call = resolve(query, ruleset)
        selection = select(call, query, modality_ruleset)
        expected = case["expect"]

        if expected.get("blocked"):
            agrees = selection.is_blocked
            detail = "everything blocked" if agrees else "expected everything to be blocked"
        else:
            indicated = {a.modality.value for a in selection.indicated}
            ruled_out = {a.modality.value for a in selection.contraindicated}
            agrees = set(expected.get("indicated", [])) <= indicated and set(
                expected.get("contraindicated", [])
            ) <= ruled_out
            detail = ", ".join(a.modality.value for a in selection.indicated) or "nothing indicated"

        if not agrees:
            misses += 1
        print(f"{'ok  ' if agrees else 'MISS'}  {case['name']}")
        print(f"        {call.mechanism} → {detail}")

    print()
    print(f"{len(cases) - misses}/{len(cases)} modality reference cases reproduced "
          f"under {modality_ruleset.pin}")
    return 1 if misses else 0


def _reference(ruleset: Ruleset, set_path: Path) -> int:
    """Re-run every published case and report agreement.

    This is the same check the test suite makes, exposed as a command because
    the people who most need to run it — the ones editing the rule file — are
    not necessarily running pytest.
    """
    cases = yaml.safe_load(set_path.read_text())["cases"]
    disagreements = 0

    for case in cases:
        call = resolve(build_query(case), ruleset)
        expected = case["expect"]["mechanism"]
        agrees = call.mechanism.value == expected
        marker = "ok  " if agrees else "MISS"
        if not agrees:
            disagreements += 1
        print(f"{marker}  {case['name']}")
        print(f"        expected {expected}, got {call.mechanism} ({call.confidence})")
        if not agrees:
            print(f"        fired: {[e.rule_id for e in call.evidence] or 'nothing'}")

    print()
    print(f"{len(cases) - disagreements}/{len(cases)} reference cases reproduced "
          f"under {ruleset.pin}")
    return 1 if disagreements else 0


def build_query(case: dict[str, Any]) -> MechanismQuery:
    """Build a query from the case-file format used by the reference set."""
    variant_spec = case["variant"]
    transcript_spec = case["transcript"]
    gene_spec = dict(case.get("gene", {}))

    exon_lengths = transcript_spec.get("exon_lengths")
    if exon_lengths is None and "uniform_exon_length" in transcript_spec:
        # For cases where what matters is the frame of the affected exon rather
        # than the true structure — a fixture convenience, and the reference
        # files say so where they use it.
        exon_lengths = [transcript_spec["uniform_exon_length"]] * transcript_spec["exon_count"]
    if exon_lengths is None:
        count = transcript_spec["exon_count"]
        total = transcript_spec["coding_length"]
        base = total // count
        exon_lengths = [base] * count
        exon_lengths[-1] += total - base * count

    distribution = MissenseDistribution(**gene_spec.pop("distribution", {}))
    gene = Gene(
        symbol=variant_spec["gene"],
        haploinsufficiency=DosageScore(gene_spec.pop("haploinsufficiency", "no_evidence")),
        triplosensitivity=DosageScore(gene_spec.pop("triplosensitivity", "no_evidence")),
        loeuf=gene_spec.pop("loeuf", None),
        forms_multimer=gene_spec.pop("forms_multimer", False),
        truncating_variants_are_milder=gene_spec.pop("truncating_variants_are_milder", False),
        imprinting=Imprinting(gene_spec.pop("imprinting", "not_imprinted")),
        silenced_allele_intact=gene_spec.pop("silenced_allele_intact", False),
        distribution=distribution,
        curated_mechanism=gene_spec.pop("curated_mechanism", None),
        curated_mechanism_source=gene_spec.pop("curated_mechanism_source", None),
    )

    # Both optional, and both absent means "not checked" rather than "checked and
    # found nothing" — a case file that says nothing about tissue gets a caveat,
    # not a silent pass.
    tissue_name = case.get("tissue")
    expression = case.get("expression")

    return MechanismQuery(
        variant=Variant(
            gene=variant_spec["gene"],
            consequence=Consequence(variant_spec["consequence"]),
            cds_position=variant_spec["cds_position"],
            protein_change=variant_spec.get("protein_change", ""),
            hgvs_c=variant_spec.get("hgvs_c", ""),
            zygosity=Zygosity(variant_spec.get("zygosity", "unknown")),
        ),
        transcript=Transcript(
            accession=transcript_spec["accession"],
            gene=variant_spec["gene"],
            coding_exon_lengths=tuple(exon_lengths),
            mane_select=transcript_spec.get("mane_select", True),
        ),
        gene=gene,
        splice=SplicePrediction(max_delta=case.get("splice_max_delta")),
        tissue=Tissue(tissue_name) if tissue_name else None,
        expression=dict(expression) if expression else None,
    )


if __name__ == "__main__":
    raise SystemExit(main())
