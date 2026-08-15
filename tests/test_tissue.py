"""Where the gene is switched on, and what that is allowed to change.

Tissue is the dimension where it is easiest to over-claim, so these tests are
mostly about restraint. Three questions run through them:

* Does the parse keep *measured zero* and *never measured* apart? They look the
  same in a report and mean opposite things, and only one of them is evidence.
* Does a silent gene caution the mechanism rather than refute it? Bulk GTEx
  hides a cell type that is two percent of the tissue, and it is post-mortem
  adult tissue besides.
* Does silence rule out only the modalities that work through the native locus?
  Gene addition delivers its own promoter, and contraindicating it here would be
  a wrong answer with a plausible reason attached.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from repairbench.cli import build_query, main
from repairbench.context.expression import (
    Tissue,
    TissueSystem,
    ingest,
    median_tpm,
    system_for,
)
from repairbench.context.registry import GeneContextRegistry
from repairbench.context.source import ContextError, Provenance, Source
from repairbench.engine import resolve
from repairbench.features import build_features
from repairbench.modality import Modality, Verdict
from repairbench.modality_rules import load_modality_ruleset
from repairbench.model import Confidence, Mechanism
from repairbench.ruleset import load_ruleset
from repairbench.selector import select

DATA = Path(__file__).parent / "data" / "context"
MATRIX = DATA / "gtex_median_tpm.tsv"
MECHANISM_RULES = Path(__file__).parents[1] / "rules" / "mechanism-v1.yaml"
MODALITY_RULES = Path(__file__).parents[1] / "rules" / "modality-v1.yaml"

BRAIN = Tissue("Brain - Cortex")
MUSCLE = Tissue("Muscle - Skeletal")


@pytest.fixture
def profiles() -> dict[str, Provenance]:
    collected: dict[str, Provenance] = {}
    ingest(Source.of("expression", MATRIX, "gtex-test"), collected)
    return collected


def profile_of(collected: dict[str, Provenance], gene: str) -> dict[str, float]:
    return dict(collected[gene].facts["expression"].value)


# --------------------------------------------------------------------------
# Reading the matrix
# --------------------------------------------------------------------------


def test_the_whole_row_is_kept_rather_than_one_number(profiles):
    """Which tissue matters is a property of the case, and the case is not known
    at ingest time. Reducing here would mean re-reading the file per case."""
    plusg = profile_of(profiles, "PLUSG")

    assert plusg["Brain - Cortex"] == 42.5
    assert plusg["Muscle - Skeletal"] == 0.1
    assert len(plusg) == 5


def test_a_gene_the_matrix_does_not_list_has_no_profile(profiles):
    assert "NOT_A_GENE" not in profiles


def test_a_value_that_is_not_a_number_stops_the_run(tmp_path: Path):
    """Skipping it would leave the gene looking measured and silent, which is
    the one wrong answer this dimension can produce that reads as a finding."""
    path = tmp_path / "broken.tsv"
    path.write_text("Name\tDescription\tLiver\nENSG1\tGENEX\tNA\n")

    with pytest.raises(ContextError, match="'NA', which is not a number"):
        ingest(Source.of("expression", path, "broken"), {})


def test_the_expression_source_is_pinned_like_any_other(profiles):
    source = profiles["PLUSG"].source_for("expression")

    assert source is not None
    assert source.digest  # earned from the file's contents, not declared


def test_only_the_genes_asked_for_are_read(tmp_path: Path):
    collected: dict[str, Provenance] = {}
    ingest(Source.of("expression", MATRIX, "gtex-test"), collected, genes={"DMD"})

    assert list(collected) == ["DMD"]


def test_expression_reaches_the_registry_beside_the_other_sources():
    registry = GeneContextRegistry.load(
        dosage=DATA / "clingen_dosage.tsv",
        expression_matrix=MATRIX,
        expression_version="gtex-test",
    )

    assert registry.gene("DMD").expression["Muscle - Skeletal"] == 120.0


# --------------------------------------------------------------------------
# Measured zero is not the same as unmeasured
# --------------------------------------------------------------------------


def test_a_measured_zero_is_evidence(profiles):
    """GTEx looked at PLUSG in liver and found nothing. That is a fact about the
    gene, and the rules are entitled to act on it."""
    assert median_tpm(profile_of(profiles, "PLUSG"), Tissue("Liver")) == 0.0


def test_a_tissue_nobody_measured_is_none_rather_than_zero(profiles):
    assert median_tpm(profile_of(profiles, "PLUSG"), Tissue("Pancreas")) is None


def test_without_an_expression_release_nothing_is_claimed(profiles):
    assert median_tpm(None, BRAIN) is None
    assert median_tpm(profile_of(profiles, "PLUSG"), None) is None


# --------------------------------------------------------------------------
# The coarse system
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("tissue", "system"),
    [
        ("Brain - Cortex", TissueSystem.CNS),
        ("Brain - Nucleus accumbens (basal ganglia)", TissueSystem.CNS),
        ("Nerve - Tibial", TissueSystem.CNS),
        ("Muscle - Skeletal", TissueSystem.MUSCLE),
        ("Heart - Left Ventricle", TissueSystem.HEART),
        ("Whole Blood", TissueSystem.BLOOD),
        ("Skin - Sun Exposed (Lower leg)", TissueSystem.SKIN),
    ],
)
def test_a_subregion_maps_onto_its_system(tissue: str, system: TissueSystem):
    """Delivery rules ask "could anything get to the brain", never "could
    anything get to the nucleus accumbens". The two vocabularies are kept apart
    so a rule cannot accidentally key on the fine one."""
    assert system_for(tissue) is system
    assert Tissue(tissue).system is system


def test_a_tissue_this_release_did_not_have_is_other_not_an_error():
    """A new GTEx release adding a tissue should produce an unclassified answer,
    not stop a clinical run."""
    assert system_for("Fallopian Tube") is TissueSystem.OTHER


# --------------------------------------------------------------------------
# Features
# --------------------------------------------------------------------------


def case(**overrides) -> dict:
    """A predicted-null variant in a haploinsufficient gene: loss of function on
    the arithmetic alone, with no curation to make the call decisive."""
    spec = {
        "variant": {
            "gene": "TESTG",
            "consequence": "stop_gained",
            "cds_position": 400,
            "zygosity": "heterozygous",
        },
        "transcript": {"accession": "NM_000001.1", "exon_count": 8, "coding_length": 2400},
        "gene": {"haploinsufficiency": "sufficient_evidence", "loeuf": 0.2},
    }
    return spec | overrides


def features_for(**overrides):
    return build_features(build_query(case(**overrides)))


def test_a_case_with_no_tissue_reports_no_expression_rather_than_silence():
    """The absence of a tissue must not read as "measured at zero", or every
    case run without one would be cautioned for the wrong reason."""
    features = features_for()

    assert features.get("tissue.known") is False
    assert features.get("expression.measured") is False
    assert features.get("expression.tpm_in_affected_tissue") is None


def test_a_tissue_the_gene_was_not_measured_in_is_also_not_silence():
    features = features_for(tissue="Pancreas", expression={"Liver": 40.0})

    assert features.get("tissue.known") is True
    assert features.get("expression.measured") is False


def test_the_measured_value_and_the_system_both_reach_the_rules():
    features = features_for(tissue="Brain - Cortex", expression={"Brain - Cortex": 42.5})

    assert features.get("expression.tpm_in_affected_tissue") == 42.5
    assert features.get("tissue.system") == "central_nervous_system"


# --------------------------------------------------------------------------
# What silence does to the mechanism
# --------------------------------------------------------------------------


def mechanism_of(**overrides):
    return resolve(build_query(case(**overrides)), load_ruleset(MECHANISM_RULES))


def test_without_the_tissue_the_call_stands_on_its_own_evidence():
    call = mechanism_of()

    assert call.mechanism is Mechanism.LOSS_OF_FUNCTION
    assert call.confidence is Confidence.PROBABLE


def test_a_gene_expressed_in_the_affected_tissue_changes_nothing():
    call = mechanism_of(tissue="Brain - Cortex", expression={"Brain - Cortex": 42.5})

    assert call.confidence is Confidence.PROBABLE
    assert "GENE_SILENT_IN_THE_AFFECTED_TISSUE" not in {e.rule_id for e in call.conflicts}


def test_a_silent_gene_is_a_caution_not_a_refutation():
    """The mechanism survives — the arithmetic that produced it has not changed.
    What drops is the confidence, and a reviewer is told why."""
    call = mechanism_of(tissue="Brain - Cortex", expression={"Brain - Cortex": 0.02})

    assert call.mechanism is Mechanism.LOSS_OF_FUNCTION
    assert call.confidence is Confidence.POSSIBLE
    assert "GENE_SILENT_IN_THE_AFFECTED_TISSUE" in {e.rule_id for e in call.conflicts}
    assert call.needs_review


def test_the_threshold_is_a_convention_and_lives_in_the_rule_file():
    """1 TPM is where the field draws the line, not something measured. A
    laboratory that wants a different floor edits the rule, not the code."""
    ruleset = load_ruleset(MECHANISM_RULES)

    assert ruleset.thresholds.expressed_above_tpm == 1.0


# --------------------------------------------------------------------------
# What silence does to the modalities
# --------------------------------------------------------------------------


def selection_for(**overrides):
    query = build_query(case(**overrides))
    call = resolve(query, load_ruleset(MECHANISM_RULES))
    return call, select(call, query, load_modality_ruleset(MODALITY_RULES))


SILENT = {"tissue": "Brain - Cortex", "expression": {"Brain - Cortex": 0.02}}
EXPRESSED = {"tissue": "Brain - Cortex", "expression": {"Brain - Cortex": 42.5}}


@pytest.mark.parametrize(
    "modality",
    [
        Modality.WILD_TYPE_UPREGULATION,
        Modality.SILENCED_ALLELE_REACTIVATION,
        Modality.ALLELE_SPECIFIC_SILENCING,
    ],
)
def test_the_modalities_that_work_through_the_native_locus_are_ruled_out(modality: Modality):
    """All three need the gene to be transcribed where the disease is. Raising,
    releasing or removing a transcript that is not there does nothing."""
    _, selection = selection_for(**SILENT)

    assert selection.verdict_for(modality) is Verdict.CONTRAINDICATED


def test_gene_addition_is_not_ruled_out_by_a_silent_native_locus():
    """The case this block exists to get right. A delivered transgene carries
    its own promoter, so the native gene being off says nothing about whether a
    supplied copy would be expressed. Contraindicating it here would be a wrong
    answer with a plausible reason attached."""
    _, selection = selection_for(**SILENT)

    contraindications = [
        evidence.rule_id
        for assessment in selection.assessments
        if assessment.modality is Modality.GENE_ADDITION
        for evidence in assessment.contraindications
    ]

    assert not any("SILENT" in rule_id for rule_id in contraindications)


def test_an_expressed_gene_leaves_those_modalities_alone():
    _, selection = selection_for(**EXPRESSED)

    ruled_out = {a.modality for a in selection.contraindicated}
    assert Modality.WILD_TYPE_UPREGULATION not in ruled_out


# --------------------------------------------------------------------------
# Caveats
# --------------------------------------------------------------------------


def test_a_case_without_a_tissue_says_the_check_was_not_made():
    _, selection = selection_for()

    assert any("no affected tissue was given" in caveat for caveat in selection.caveats)


def test_a_case_with_a_tissue_says_delivery_there_is_still_unanswered():
    """The honest limit of the whole package. Knowing a gene is on in the cortex
    is not knowing that anything can be got into the cortex, and the second
    question is the one that decides most of these in practice."""
    _, selection = selection_for(**EXPRESSED)

    delivery = [c for c in selection.caveats if "delivery to" in c]
    assert delivery
    assert "central_nervous_system" in delivery[0]
    assert "Brain - Cortex" in delivery[0]


# --------------------------------------------------------------------------
# What the operator sees
# --------------------------------------------------------------------------


def context_output(*argv: str, capsys) -> str:
    assert main(["context", *argv]) == 0
    return capsys.readouterr().out


def test_the_report_says_expressed_or_silent_against_the_rule_file_threshold(capsys):
    out = context_output("SCN1A", "--tissue", "Brain - Cortex", capsys=capsys)

    assert "Brain - Cortex (central_nervous_system): 35 TPM — expressed" in out


def test_the_report_will_not_call_an_unmeasured_tissue_silent(capsys):
    """The distinction the whole module rests on, at the point an operator reads
    it: nobody looked is not the same as nothing there."""
    out = context_output("SCN1A", "--tissue", "Pancreas", capsys=capsys)

    assert "Pancreas: not measured in this release" in out
    assert "silent" not in out


def test_the_profile_does_not_bury_the_facts_beside_it(capsys):
    """Fifty numbers printed in full would push the four single-valued facts off
    the screen, so the profile is summarised rather than dumped."""
    out = context_output("SCN1A", capsys=capsys)

    assert "expression = 5 tissues measured" in out


def test_the_two_tissue_caveats_are_never_both_raised():
    """They contradict each other: one says no tissue was given, the other names
    it."""
    for overrides in ({}, EXPRESSED, SILENT):
        _, selection = selection_for(**overrides)
        tissue_caveats = [
            c
            for c in selection.caveats
            if "no affected tissue was given" in c or "delivery to" in c
        ]
        assert len(tissue_caveats) <= 1
