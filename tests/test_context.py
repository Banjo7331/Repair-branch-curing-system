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

from pathlib import Path

import pytest
import yaml

from repairbench.context.registry import GeneContextRegistry
from repairbench.context.source import ContextError, Fact, Provenance, Source
from repairbench.engine import resolve
from repairbench.features import MechanismQuery, Variant
from repairbench.model import (
    Consequence,
    DosageScore,
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
