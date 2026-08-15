"""The seam: a contraindicated modality is never handed to a designer.

Every other test file in this package checks one layer. This one checks that the
layers are wired together in the one direction that matters. M6 can rule gene
addition out perfectly and it buys nothing if M7 will still design for it when
asked separately — which is exactly what the package did before this module
existed, because the two commands did not know about each other.

So the assertions here are mostly negative: no molecule for a ruled-out
modality, no molecule under an unresolved mechanism, no molecule from a designer
that was never routed to. A negative assertion is a weak test in general and the
right one here, because the failure it guards against is output appearing where
there should be none.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
import yaml

from repairbench.annotation.fasta import IndexedFasta
from repairbench.annotation.gff import parse_gff3
from repairbench.annotation.store import TranscriptStore
from repairbench.cli import build_query, main
from repairbench.design.aso import AsoOutcome
from repairbench.design.candidate import DesignOutcome
from repairbench.design.editors import load_editors
from repairbench.design.flags import load_flag_rules
from repairbench.design.prime import PrimeOutcome
from repairbench.design.sequence import reverse_complement
from repairbench.engine import resolve
from repairbench.modality import Modality, Verdict
from repairbench.modality_rules import load_modality_ruleset
from repairbench.plan import Designers, Locus, PlanError, load_routing, plan
from repairbench.ruleset import RulesetError, load_ruleset
from repairbench.selector import select

ROOT = Path(__file__).parents[1]
DATA = Path(__file__).parent / "data" / "design"
CASE = DATA / "case.yaml"
FASTA = DATA / "target.fa"
GFF = DATA / "target.gff3"

MECHANISM_RULES = ROOT / "rules" / "mechanism-v1.yaml"
MODALITY_RULES = ROOT / "rules" / "modality-v1.yaml"
ROUTING = ROOT / "rules" / "routing-v1.yaml"

#: The fixture's variant: a missense at 17:301, reference G, patient A.
POSITION, REFERENCE, ALTERNATE = 301, "G", "A"


@pytest.fixture
def designers():
    return Designers(
        editors=load_editors(ROOT / "rules" / "editors-v1.yaml"),
        prime=load_flag_rules(ROOT / "rules" / "prime-v1.yaml"),
        aso=load_flag_rules(ROOT / "rules" / "aso-v1.yaml"),
        routing=load_routing(ROOT / "rules" / "routing-v1.yaml"),
    )


@pytest.fixture
def case() -> dict[str, Any]:
    return yaml.safe_load(CASE.read_text())


@pytest.fixture
def genome():
    with IndexedFasta(FASTA) as fasta:
        yield fasta


def assemble(case, designers, genome=None, **kwargs):
    query = build_query(case)
    mechanism_rules = load_ruleset(MECHANISM_RULES)
    modality_rules = load_modality_ruleset(MODALITY_RULES)
    call = resolve(query, mechanism_rules)
    selection = select(call, query, modality_rules)
    store = TranscriptStore(parse_gff3(GFF))
    return plan(
        query,
        call,
        selection,
        designers,
        locus=kwargs.pop("locus", Locus.from_case(case, store)),
        sequences=genome,
        mechanism_rules=mechanism_rules,
        modality_rules=modality_rules,
        **kwargs,
    )


def design_for(result, modality: Modality):
    return next(design for design in result.designs if design.modality is modality)


# --------------------------------------------------------------------------
# The rule this module exists for
# --------------------------------------------------------------------------


def test_a_contraindicated_modality_is_never_designed(case, designers, genome):
    """The variant is a dominant negative, so supplying a normal copy leaves
    every poisoning subunit where it was. M6 rules gene addition out; nothing
    below M6 may produce a sequence for it."""
    result = assemble(case, designers, genome)
    gene_addition = design_for(result, Modality.GENE_ADDITION)

    assert gene_addition.verdict is Verdict.CONTRAINDICATED
    assert gene_addition.outcome is None
    assert "ruled out by the modality rules" in gene_addition.refusal


def test_no_ruled_out_modality_anywhere_carries_a_molecule(case, designers, genome):
    """Stated over the whole set rather than one modality, because the failure
    this guards against is a route added later that skips the check."""
    result = assemble(case, designers, genome)

    for design in result.designs:
        if design.verdict is Verdict.CONTRAINDICATED:
            assert design.outcome is None
            assert design.refusal


def test_the_refusal_is_a_refusal_and_not_a_warning(case, designers, genome):
    """A sequence is something a reader can order; a caveat is something a
    reader can skim. The distinction is the whole design of this seam."""
    result = assemble(case, designers, genome)
    ruled_out = result.ruled_out

    assert ruled_out
    assert all(design.outcome is None for design in ruled_out)


def test_an_unresolved_mechanism_designs_nothing(case, designers, genome):
    """Every modality depends on the mechanism, so a molecule designed under an
    undetermined one is designed against nothing."""
    weakened = {
        **case,
        "gene": {},
        "variant": {**case["variant"], "consequence": "synonymous_variant"},
    }
    result = assemble(weakened, designers, genome)

    assert not result.has_designs
    assert all(design.outcome is None for design in result.designs)
    assert any("unresolved" in note for note in result.notes)


# --------------------------------------------------------------------------
# Routing
# --------------------------------------------------------------------------


def test_each_indicated_modality_reaches_the_designer_that_belongs_to_it(
    case, designers, genome
):
    result = assemble(case, designers, genome)

    assert isinstance(design_for(result, Modality.BASE_EDITING).outcome, DesignOutcome)
    assert isinstance(design_for(result, Modality.PRIME_EDITING).outcome, PrimeOutcome)
    assert isinstance(
        design_for(result, Modality.ALLELE_SPECIFIC_SILENCING).outcome, AsoOutcome
    )


def test_silencing_gets_a_cleaving_chemistry_and_skipping_would_not(designers):
    """The one confusion that inverts a therapy, decided in the routing file
    rather than in code: destroying the transcript is the point of silencing and
    the opposite of the point of skipping."""
    routing = designers.routing

    assert routing.route_for(Modality.ALLELE_SPECIFIC_SILENCING).chemistry == "gapmer-2MOE"
    assert routing.route_for(Modality.EXON_SKIPPING).chemistry.startswith("steric")
    assert routing.route_for(Modality.SPLICE_CORRECTION).chemistry.startswith("steric")


def test_a_modality_with_no_designer_says_what_would_have_to_be_designed(designers):
    """"No candidates" and "the thing to design here is a vector and a promoter"
    look identical in an empty list and mean entirely different things."""
    route = designers.routing.route_for(Modality.GENE_ADDITION)

    assert not route.has_designer
    assert "vector" in route.because


def test_a_routing_table_missing_a_modality_is_refused(tmp_path: Path):
    """An unlisted modality would produce no output and look exactly like a
    modality nothing was found for."""
    path = tmp_path / "routing.yaml"
    path.write_text(
        "version: partial\ndesigners:\n  - {modality: base_editing, designer: base_editor}\n"
    )

    with pytest.raises(RulesetError, match="no routing entry for"):
        load_routing(path)


def test_every_modality_in_the_shipped_table_has_a_reason(designers):
    for modality in Modality:
        assert designers.routing.route_for(modality).because


# --------------------------------------------------------------------------
# What the designers are given
# --------------------------------------------------------------------------


def test_the_patient_carries_the_alternate_allele(case, designers, genome):
    """Swapping the two produces a design that installs the disease, and it is
    the kind of swap that never raises."""
    result = assemble(case, designers, genome)
    outcome = design_for(result, Modality.BASE_EDITING).outcome

    assert outcome.patient_base == ALTERNATE
    assert outcome.wild_type_base == REFERENCE


def test_a_discriminating_oligonucleotide_must_cover_the_variant(case, designers, genome):
    """The variant is the only thing telling the two transcripts apart. A gapmer
    that misses it knocks down the healthy allele as well — which for a
    dominant-negative variant removes the good product along with the bad."""
    result = assemble(case, designers, genome)
    outcome = design_for(result, Modality.ALLELE_SPECIFIC_SILENCING).outcome

    assert outcome.candidates
    for candidate in outcome.candidates:
        assert candidate.span[0] <= POSITION <= candidate.span[1]
    assert any("telling the two transcripts apart" in note for note in outcome.notes)


def test_the_oligonucleotide_is_complementary_to_the_patients_allele(case, designers, genome):
    """Tiled against the patient's sequence rather than the reference. One
    complementary to the reference base is complementary to the transcript it
    was supposed to spare."""
    result = assemble(case, designers, genome)
    outcome = design_for(result, Modality.ALLELE_SPECIFIC_SILENCING).outcome
    candidate = outcome.candidates[0]

    offset = POSITION - candidate.span[0]
    assert candidate.target[offset] == ALTERNATE
    assert candidate.sequence == reverse_complement(candidate.target)


def test_a_case_with_no_genomic_coordinate_designs_nothing_and_says_why(case, designers, genome):
    """M5 and M6 reason in CDS offsets; a protospacer is placed on the genome,
    and without an annotation the two are not interconvertible."""
    result = assemble(case, designers, genome, locus=None)

    assert not result.has_designs
    assert any("no genomic coordinate" in note for note in result.notes)


def test_without_a_reference_sequence_nothing_is_designed(case, designers):
    result = assemble(case, designers, genome=None)

    assert not result.has_designs
    assert any("no reference sequence" in note for note in result.notes)


def test_a_genomic_block_without_alleles_is_refused(case):
    partial = {**case, "genomic": {"chromosome": "17", "position": 301}}

    with pytest.raises(PlanError, match="reference and alternate alleles"):
        Locus.from_case(partial, None)


# --------------------------------------------------------------------------
# Reproducibility
# --------------------------------------------------------------------------


def test_a_plan_names_every_rule_file_it_was_made_under(case, designers, genome):
    """Six files decide a plan. One that cannot name them cannot be reproduced,
    and cannot be compared with the plan the same case gets next year."""
    result = assemble(case, designers, genome)
    joined = " ".join(result.pins)

    for version in ("mechanism-v1@", "modality-v1@", "editors-v1@", "prime-v1@", "aso-v1@",
                    "routing-v1@"):
        assert version in joined


# --------------------------------------------------------------------------
# What the operator sees
# --------------------------------------------------------------------------


def test_the_report_puts_the_ruled_out_list_before_the_molecules(capsys):
    """Same ordering M6 uses, for the same reason: the output that changes what
    somebody does is the one that closes a door."""
    assert (
        main(
            [
                "plan",
                str(CASE),
                "--fasta", str(FASTA),
                "--annotation", str(GFF),
                "--limit", "1",
            ]
        )
        == 0
    )
    out = capsys.readouterr().out

    assert out.index("not designed, ruled out") < out.index("\n  designed\n")
    assert "gene_addition" in out
    assert "routing-v1@" in out


def test_the_command_fails_when_the_mechanism_is_unresolved(capsys, tmp_path: Path):
    case = yaml.safe_load(CASE.read_text())
    case["gene"] = {}
    case["variant"]["consequence"] = "synonymous_variant"
    path = tmp_path / "weak.yaml"
    path.write_text(yaml.safe_dump(case))

    assert main(["plan", str(path), "--fasta", str(FASTA), "--annotation", str(GFF)]) == 2
    assert "unresolved" in capsys.readouterr().out
