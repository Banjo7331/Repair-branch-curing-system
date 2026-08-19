"""Reading gene-level facts from files instead of from a fixture.

Two things are under test and they are not the same thing.

Most of what follows tests the **ingest**: that ClinGen's non-scale is not
flattened, that a non-MANE constraint row is skipped, that a local claim without
a citation is refused, that two sources disagreeing is an error rather than a
silent preference.

The last test is the **bridge**: it runs the whole mechanism reference set with
the gene context assembled from these files rather than typed into the case, and
asserts the mechanisms come out the same. That proves the wiring end to end. It
does *not* prove the biology, because the fixture files below were generated
from the values that were already inline — every number in them was originally
typed by hand, and pointing the loader at the real ClinGen and gnomAD downloads
is what would make "reproduced" mean something stronger.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest
import yaml

from repairbench.context import gnomad
from repairbench.context.registry import GeneContextRegistry
from repairbench.context.source import ContextError, Fact, Provenance, Source, read_tsv
from repairbench.engine import resolve
from repairbench.features import MechanismQuery, Variant
from repairbench.model import (
    Consequence,
    DosageScore,
    Gene,
    Mechanism,
    MissenseDistribution,
    Zygosity,
)
from repairbench.ruleset import load_ruleset
from repairbench.transcript import Transcript

DATA = Path(__file__).parent / "data" / "context"
DOSAGE = DATA / "clingen_dosage.tsv"
CONSTRAINT = DATA / "gnomad_constraint.tsv"
LOCAL = DATA / "local_curation.yaml"
REFERENCE = Path(__file__).parent / "reference" / "mechanisms.yaml"
RULES = Path(__file__).parents[1] / "rules" / "mechanism-v1.yaml"


@pytest.fixture
def registry() -> GeneContextRegistry:
    return GeneContextRegistry.load(
        dosage=DOSAGE,
        constraint=CONSTRAINT,
        local=LOCAL,
        dosage_version="2026-01",
        constraint_version="v4.1",
        local_version="rev3",
    )


# --------------------------------------------------------------------------
# Pins that are earned rather than declared
# --------------------------------------------------------------------------


def test_each_source_is_pinned_by_the_digest_of_its_bytes(registry: GeneContextRegistry):
    """The whole reason for this package: before it, a world could claim
    'gene_curation@2026-01' with nothing behind the label."""
    pins = {source.name: source for source in registry.pins}

    assert set(pins) == {"clingen_dosage", "gnomad_constraint", "local_curation"}
    for source in pins.values():
        assert len(source.digest) == 64
        assert source.pin.count("@") == 1


def test_the_digest_changes_when_the_file_does(tmp_path: Path):
    edited = tmp_path / "dosage.tsv"
    edited.write_text(DOSAGE.read_text() + "# a curator's note\n")

    assert Source.of("d", edited, "x").digest != Source.of("d", DOSAGE, "x").digest


def test_provenance_is_recorded_per_fact_not_per_gene(registry: GeneContextRegistry):
    """A report saying 'dominant negative, because null alleles are milder' is
    only reviewable if a reader can tell that claim came from our own curation
    while the dosage score came from ClinGen."""
    sourced = registry.gene("COL1A1")

    assert "clingen_dosage" in sourced.cite("haploinsufficiency")
    assert "gnomad_constraint" in sourced.cite("loeuf")
    assert "local_curation" in sourced.cite("truncating_variants_are_milder")


# --------------------------------------------------------------------------
# ClinGen's non-scale
# --------------------------------------------------------------------------


def test_a_refutation_is_not_flattened_into_a_low_score(registry: GeneContextRegistry):
    """30 and 40 are not points on the 0-3 scale. 'No evidence yet' is a gap and
    'refuted' is a finding, and a predicted null means something different
    under each."""
    refuted = registry.gene("SYNTHETIC_HI_REFUTED").gene

    assert refuted.haploinsufficiency is DosageScore.UNLIKELY
    assert refuted.haploinsufficiency.refutes_haploinsufficiency
    assert not refuted.haploinsufficiency.supports_haploinsufficiency


def test_an_unevaluated_gene_contributes_nothing_rather_than_a_default(tmp_path: Path):
    """Recording 'no evidence' where ClinGen has not looked would put a claim in
    the provenance that ClinGen never made."""
    sparse = tmp_path / "dosage.tsv"
    sparse.write_text(
        "#Gene Symbol\tHaploinsufficiency Score\tTriplosensitivity Score\n"
        "UNSEEN\tNot yet evaluated\tNot yet evaluated\n"
    )

    registry = GeneContextRegistry.load(dosage=sparse)

    with pytest.raises(ContextError, match="not present in any loaded source"):
        registry.provenance_for("UNSEEN")


def test_an_unrecognised_dosage_code_is_refused(tmp_path: Path):
    bad = tmp_path / "dosage.tsv"
    bad.write_text(
        "#Gene Symbol\tHaploinsufficiency Score\tTriplosensitivity Score\nG\t7\t0\n"
    )

    with pytest.raises(ContextError, match="not a ClinGen dosage code"):
        GeneContextRegistry.load(dosage=bad)


# --------------------------------------------------------------------------
# Constraint
# --------------------------------------------------------------------------


def test_only_the_mane_select_row_is_read(registry: GeneContextRegistry):
    """The file has one row per transcript and the numbers differ between them;
    taking whichever came first would make the value depend on file order."""
    assert registry.gene("COL1A1").gene.loeuf == 0.25


def test_missing_constraint_stays_missing(tmp_path: Path):
    """Constraint is undefined for short genes. The rules already refuse to fire
    on a missing value, so absent is the honest answer."""
    partial = tmp_path / "constraint.tsv"
    partial.write_text("gene\tmane_select\tlof.oe_ci.upper\nSHORTG\ttrue\tNA\n")

    registry = GeneContextRegistry.load(dosage=DOSAGE, constraint=partial)

    assert registry.gene("SCN1A").gene.loeuf is None


def test_a_renamed_column_fails_the_run_rather_than_reading_the_wrong_one(tmp_path: Path):
    """Column order changes between releases; a positional parser would keep
    working and silently read the wrong column."""
    renamed = tmp_path / "constraint.tsv"
    renamed.write_text("gene\tmane_select\tloeuf\nG\ttrue\t0.2\n")

    with pytest.raises(ContextError, match=r"lof\.oe_ci\.upper"):
        GeneContextRegistry.load(constraint=renamed)


# --------------------------------------------------------------------------
# The facts nobody publishes
# --------------------------------------------------------------------------


def test_a_local_claim_without_a_citation_is_refused(tmp_path: Path):
    """A clinical assertion with no pointer to where it was established is
    folklore with a version number."""
    uncited = tmp_path / "local.yaml"
    uncited.write_text(
        yaml.safe_dump({"genes": {"G": {"forms_multimer": {"value": True}}}})
    )

    with pytest.raises(ContextError, match="no citation"):
        GeneContextRegistry.load(local=uncited)


def test_local_curation_may_not_override_a_published_fact(tmp_path: Path):
    """Only the two fields with no public table belong here. A local override of
    something ClinGen publishes is a way to be quietly wrong."""
    overreaching = tmp_path / "local.yaml"
    published = {"value": "sufficient_evidence", "citation": "me"}
    overreaching.write_text(yaml.safe_dump({"genes": {"G": {"haploinsufficiency": published}}}))

    with pytest.raises(ContextError, match="may not supply"):
        GeneContextRegistry.load(local=overreaching)


def test_a_curated_mechanism_carries_its_citation_into_the_gene(registry: GeneContextRegistry):
    """The model already refuses a curated mechanism with no source; this is
    where the source comes from."""
    mecp2 = registry.gene("MECP2").gene

    assert mecp2.curated_mechanism == "loss_of_function"
    assert mecp2.curated_mechanism_source


# --------------------------------------------------------------------------
# Disagreement
# --------------------------------------------------------------------------


def test_two_sources_disagreeing_is_an_error_not_a_preference():
    """Picking one silently is how a pipeline develops opinions nobody chose."""
    provenance = Provenance(gene="G")
    first = Source(name="a", path=DOSAGE, digest="x" * 64, version="1")
    second = Source(name="b", path=CONSTRAINT, digest="y" * 64, version="1")

    provenance.record(Fact(field="loeuf", value=0.2, source=first))

    with pytest.raises(ContextError, match="Two sources disagree"):
        provenance.record(Fact(field="loeuf", value=0.9, source=second))


def test_a_registry_with_no_sources_is_refused():
    with pytest.raises(ContextError, match="no sources"):
        GeneContextRegistry.load()


# --------------------------------------------------------------------------
# The bridge
# --------------------------------------------------------------------------


def test_the_reference_set_reproduces_under_ingested_context(registry: GeneContextRegistry):
    """Every mechanism reference case, with the gene context read from files.

    This proves the wiring, not the biology: the fixture files were generated
    from the values that were already inline in the reference set. What it rules
    out is a whole class of quiet failure — a field ingested under the wrong
    name, a score mapped to the wrong enum, a MANE row picked wrongly — any of
    which would move a mechanism here and nowhere else.
    """
    rules = load_ruleset(RULES)
    cases = yaml.safe_load(REFERENCE.read_text())["cases"]

    for case in cases:
        symbol = case["variant"]["gene"]
        inline = case.get("gene", {})
        distribution = MissenseDistribution(**inline.get("distribution", {}))

        transcript_spec = case["transcript"]
        exon_lengths = transcript_spec.get("exon_lengths")
        if exon_lengths is None:
            count, total = transcript_spec["exon_count"], transcript_spec["coding_length"]
            base = total // count
            exon_lengths = [base] * count
            exon_lengths[-1] += total - base * count

        query = MechanismQuery(
            variant=Variant(
                gene=symbol,
                consequence=Consequence(case["variant"]["consequence"]),
                cds_position=case["variant"]["cds_position"],
                zygosity=Zygosity(case["variant"].get("zygosity", "unknown")),
            ),
            transcript=Transcript(
                accession=transcript_spec["accession"],
                gene=symbol,
                coding_exon_lengths=tuple(exon_lengths),
                mane_select=True,
            ),
            gene=registry.gene(symbol, distribution=distribution).gene,
        )

        call = resolve(query, rules)

        assert call.mechanism is Mechanism(case["expect"]["mechanism"]), (
            f"{case['name']}: under ingested context the mechanism came out as "
            f"{call.mechanism}, not {case['expect']['mechanism']} — "
            f"dosage cited as {registry.gene(symbol).cite('haploinsufficiency')}"
        )


# --------------------------------------------------------------------------
# The shapes real files actually arrive in
# --------------------------------------------------------------------------


def test_a_commented_header_below_a_preamble_is_found(tmp_path: Path):
    """ClinGen's list opens with five lines of provenance and then a header
    that is itself commented. Skipping every '#' line — which this used to do —
    made the first line of provenance into the column names, and the failure
    surfaced as "no Gene Symbol column" about a file whose column is Gene
    Symbol."""
    path = tmp_path / "clingen.tsv"
    path.write_text(
        "#ClinGen Gene Curation Results\n"
        "#17 Aug,2026\n"
        "#Genomic Locations are reported on GRCh38\n"
        "#Gene Symbol\tHaploinsufficiency Score\tTriplosensitivity Score\n"
        "GENEX\t3\t0\n"
    )

    rows = list(read_tsv(Source.of("clingen", path, "test"), required={"Gene Symbol"}))

    assert rows == [
        {"Gene Symbol": "GENEX", "Haploinsufficiency Score": "3", "Triplosensitivity Score": "0"}
    ]


def test_a_gct_preamble_is_read_as_the_format_declares_it(tmp_path: Path):
    """A GCT says what it is on line one and gives its dimensions on line two.
    That second line is tab-separated and uncommented, so anything hunting for
    "the first line with tabs" takes it for the header and reads the row count
    as a gene name."""
    path = tmp_path / "expression.gct"
    path.write_text("#1.2\n2\t1\nName\tDescription\tLiver\nENSG1\tGENEX\t42.5\n")

    rows = list(read_tsv(Source.of("expression", path, "test"), required={"Description"}))

    assert rows == [{"Name": "ENSG1", "Description": "GENEX", "Liver": "42.5"}]


def test_a_gzipped_source_is_read_without_being_unpacked_first(tmp_path: Path):
    """Every one of these releases ships compressed and several ship only
    compressed."""
    path = tmp_path / "constraint.tsv.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("gene\tlof.oe_ci.upper\nGENEX\t0.25\n")

    rows = list(read_tsv(Source.of("constraint", path, "test"), required={"gene"}))

    assert rows == [{"gene": "GENEX", "lof.oe_ci.upper": "0.25"}]


def test_an_empty_file_yields_nothing_rather_than_raising(tmp_path: Path):
    path = tmp_path / "empty.tsv"
    path.write_text("")

    assert list(read_tsv(Source.of("empty", path, "test"), required=set())) == []


# --------------------------------------------------------------------------
# gnomAD ships its columns under different names in different releases
# --------------------------------------------------------------------------


def v2_constraint(path: Path) -> Path:
    """A miniature file in the v2.1.1 shape: bgzipped, Ensembl canonical."""
    with gzip.open(path, "wt") as handle:
        handle.write("gene\ttranscript\tcanonical\tobs_lof\toe_lof_upper\tpLI\n")
        handle.write("GENEX\tENST00000303395\ttrue\t3\t0.14\t1.0\n")
        handle.write("GENEX\tENST00000000000\tfalse\t9\t0.99\t0.0\n")
    return path


def test_the_v2_schema_is_detected_from_the_header(tmp_path: Path):
    """The release calls LOEUF `oe_lof_upper` and marks the row with
    `canonical`; v4 calls them `lof.oe_ci.upper` and `mane_select`. This module
    was written against v4 names, and the file the downloads page links is v2."""
    source = Source.of("gnomad", v2_constraint(tmp_path / "v2.txt.bgz"), "v2.1.1")

    assert gnomad.detect_schema(source).release == "v2.1.1"


def test_the_row_the_release_prefers_is_the_one_read(tmp_path: Path):
    into: dict[str, Provenance] = {}
    gnomad.ingest(Source.of("gnomad", v2_constraint(tmp_path / "v2.txt.bgz"), "v2.1.1"), into)

    assert into["GENEX"].facts["loeuf"].value == 0.14


def test_which_transcript_the_value_came_through_is_recorded(tmp_path: Path):
    """MANE Select and Ensembl canonical usually agree and are not the same
    claim, so a reader is told which one this number is about."""
    into: dict[str, Provenance] = {}
    gnomad.ingest(Source.of("gnomad", v2_constraint(tmp_path / "v2.txt.bgz"), "v2.1.1"), into)

    assert "Ensembl canonical" in into["GENEX"].facts["loeuf"].citation


def test_the_v4_schema_still_reads(tmp_path: Path):
    into: dict[str, Provenance] = {}
    gnomad.ingest(Source.of("gnomad", DATA / "gnomad_constraint.tsv", "v4.1"), into)

    assert into["COL1A1"].facts["loeuf"].value == 0.25
    assert "MANE Select" in into["COL1A1"].facts["loeuf"].citation


def test_a_file_in_neither_schema_names_every_spelling_it_looked_for(tmp_path: Path):
    """A constraint value read out of the wrong column is a number that looks
    entirely reasonable, so nothing here matches by position or resemblance."""
    path = tmp_path / "other.tsv"
    path.write_text("gene\tsome_other_metric\nGENEX\t0.5\n")

    with pytest.raises(ContextError, match="oe_lof_upper"):
        gnomad.ingest(Source.of("gnomad", path, "unknown"), {})


# --------------------------------------------------------------------------
# "Nobody looked" is not "somebody looked and found nothing"
# --------------------------------------------------------------------------


def test_a_gene_clingen_never_evaluated_is_not_scored(tmp_path: Path):
    """The default used to be `no_evidence`, which reads as a ClinGen finding.
    ClinGen has curated a few thousand genes, not twenty thousand, and for the
    rest it has said nothing at all."""
    registry = GeneContextRegistry.load(
        dosage=DATA / "clingen_dosage.tsv", local=DATA / "local_curation.yaml"
    )

    gene = registry.gene("COL1A1").gene

    assert gene.haploinsufficiency is DosageScore.SUFFICIENT_EVIDENCE
    assert Gene(symbol="NEVER_CURATED").haploinsufficiency is DosageScore.NOT_EVALUATED


def test_not_evaluated_neither_supports_nor_refutes():
    """Both directions must be false, or an absence of curation would argue for
    something."""
    absent = DosageScore.NOT_EVALUATED

    assert not absent.supports_haploinsufficiency
    assert not absent.refutes_haploinsufficiency
    assert not absent.is_curated


def test_no_evidence_is_a_curation_and_says_so():
    """ClinGen scoring a gene 0 means it looked. That is weaker than a
    refutation and stronger than silence, and it is the middle of the three."""
    scored = DosageScore.NO_EVIDENCE

    assert scored.is_curated
    assert not scored.supports_haploinsufficiency
    assert not scored.refutes_haploinsufficiency


# --------------------------------------------------------------------------
# ClinVar, joined in like any other source
# --------------------------------------------------------------------------

VARIANT_SUMMARY = DATA / "clinvar_variant_summary.txt"


def test_the_distribution_is_now_an_ingested_fact_with_a_source():
    """It used to be a parameter and nothing else, which meant the most
    inferential rule in the package read a number with no provenance at all."""
    registry = GeneContextRegistry.load(
        dosage=DOSAGE, variant_summary=VARIANT_SUMMARY, genes={"COL1A1"}
    )

    sourced = registry.gene("COL1A1")
    assert sourced.gene.distribution.pathogenic_missense_total == 3
    assert sourced.cite("distribution").startswith("clinvar@")


def test_a_passed_distribution_still_wins_over_the_counted_one():
    """A curated count from a disease database is better evidence than a
    submission tally, and the reference set has to be able to hold a gene still
    while a rule is examined."""
    registry = GeneContextRegistry.load(
        dosage=DOSAGE, variant_summary=VARIANT_SUMMARY, genes={"COL1A1"}
    )

    given = MissenseDistribution(
        pathogenic_missense_total=700,
        pathogenic_missense_in_hotspot=600,
        pathogenic_truncating_total=200,
    )
    assert registry.gene("COL1A1", distribution=given).gene.distribution == given


def test_loading_clinvar_without_naming_genes_is_refused():
    """The file is millions of submissions. Reading all of them to build
    context for a handful would take minutes and look like it was working."""
    with pytest.raises(ContextError, match="needs the genes to look for"):
        GeneContextRegistry.load(dosage=DOSAGE, variant_summary=VARIANT_SUMMARY)


def test_a_gene_clinvar_has_nothing_for_gets_no_distribution_rather_than_zero():
    """Zero truncating variants is half the gain-of-function rule. A gene that
    was simply not counted must not supply it."""
    registry = GeneContextRegistry.load(
        dosage=DOSAGE, variant_summary=VARIANT_SUMMARY, genes={"COL1A1", "SCN1A"}
    )

    assert "distribution" not in registry.provenance_for("SCN1A").facts
    assert registry.gene("SCN1A").cite("distribution") == "not supplied by any source"


def test_an_uncounted_distribution_prints_as_uncounted_not_as_zeroes():
    assert str(MissenseDistribution()) == "nothing counted"
