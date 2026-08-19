"""The reference set, run against counts nobody here typed.

``tests/reference/mechanisms.yaml`` carries a ``distribution:`` block per gene:
how many pathogenic missense variants, how many of them fall in one stretch of
protein, how many truncating variants. It is the input to the least direct
inference in the package — the rule that reads clustering-without-truncation as
gain of function — and every one of those numbers was typed in from memory of
the literature.

This file replaces them with counts read out of a real ClinVar release and asks
whether the reference set still reproduces. That is a different question from
the one ``test_reference_set.py`` asks. That file tests whether the rules are
consistent with what they were fed; this one tests whether they survive being
fed something they were not tuned against.

There are three ways it can come out, and all three are worth having:

* **Same mechanism.** The rule was reading a real pattern and the invented
  numbers happened to describe it. That is the result the reference set claims.
* **Different mechanism.** Either the invented numbers were flattering the rule,
  or the threshold is in the wrong place. The disagreement gets written down in
  the README rather than tuned away — that is the standing rule for this set.
* **Nothing to count.** ClinVar files the gene under another symbol, or every
  submission is below the review threshold. Reported as a skip for that gene,
  never as zero: zero truncating variants is half of the gain-of-function rule,
  and inventing it would argue for a mechanism.

Skips without ``refdata/variant_summary.txt.gz`` — 250 MB that is not ours to
redistribute, fetched by ``scripts/fetch-reference-data.sh clinvar`` — so CI
runs the rest of the suite and this runs where the data is.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from repairbench.context.clinvar import distribution_for, group_by_gene, read_variants
from repairbench.context.source import Source
from repairbench.engine import resolve
from repairbench.model import Mechanism
from repairbench.ruleset import load_ruleset
from test_reference_set import RULES, build_query, load_cases

ROOT = Path(__file__).parents[1]
REFDATA = ROOT / "refdata"
#: The whole release, or an extract of it filtered to these genes. The extract
#: is 3 MB against 440, produced by ``fetch-reference-data.sh clinvar`` from the
#: file it just downloaded, and preferred because a test that takes a minute to
#: read a file it uses a thousandth of is a test people stop running. Both are
#: parsed identically — the extract keeps ClinVar's own header, and the gene
#: filter that produced it is a superset of the one the parser applies.
CANDIDATES = (REFDATA / "clinvar_refset.tsv.gz", REFDATA / "variant_summary.txt.gz")
SUMMARY = next((path for path in CANDIDATES if path.exists()), CANDIDATES[-1])

pytestmark = pytest.mark.skipif(
    not SUMMARY.exists(),
    reason=f"{SUMMARY} absent — run scripts/fetch-reference-data.sh clinvar",
)

#: Cases whose gene is invented, and which therefore have nothing to look up.
#: Named rather than detected, so that a real gene ClinVar happens to be silent
#: about is reported as a gap instead of quietly joining this list.
SYNTHETIC = {"SYNTHETIC_HI_REFUTED", "SYNTHETIC_SPARSE"}


@pytest.fixture(scope="module")
def counted() -> dict[str, object]:
    """Every reference-set gene's real distribution, in one pass of the file.

    One pass because the file is millions of rows: a fixture per case would read
    it once per case, and a suite that takes an hour is a suite nobody runs.
    """
    wanted = {
        case["variant"]["gene"] for case in load_cases() if case["variant"]["gene"] not in SYNTHETIC
    }
    source = Source.of("clinvar", SUMMARY, "refdata")
    grouped = group_by_gene(read_variants(source, genes=wanted))
    return {gene: distribution_for(found) for gene, found in grouped.items()}


@pytest.mark.parametrize("case", load_cases(), ids=lambda c: c["name"])
def test_the_reference_case_survives_real_counts(case: dict, counted: dict) -> None:
    gene = case["variant"]["gene"]
    if gene in SYNTHETIC:
        pytest.skip(f"{gene} is invented, and the case says so")
    if gene not in counted:
        pytest.skip(f"ClinVar has no pathogenic submissions at ≥1★ under the symbol {gene}")

    query = build_query(case)
    real = counted[gene]
    with_real_counts = replace(query, gene=replace(query.gene, distribution=real))

    call = resolve(with_real_counts, load_ruleset(RULES))

    assert call.mechanism is Mechanism(case["expect"]["mechanism"]), (
        f"{case['name']}\n"
        f"  typed counts:   {query.gene.distribution}\n"
        f"  ClinVar counts: {real}\n"
        f"  expected {case['expect']['mechanism']}, got {call.mechanism}\n"
        f"  evidence: {[e.rule_id for e in call.evidence]}"
    )


def test_the_gain_of_function_gene_really_does_lack_truncating_variants(counted: dict) -> None:
    """The claim underneath the clustering rule, checked directly.

    *PIK3CA* causes disease by the product doing something new. If the rule is
    reading a real signal rather than a tuned one, the real file should show its
    pathogenic variation to be overwhelmingly missense — truncating variants in
    this gene are not a disease mechanism, they are a different gene's biology.
    """
    if "PIK3CA" not in counted:
        pytest.skip("no PIK3CA submissions counted")

    pik3ca = counted["PIK3CA"]
    assert pik3ca.pathogenic_missense_total > 0
    assert pik3ca.pathogenic_truncating_total < pik3ca.pathogenic_missense_total


def test_the_haploinsufficiency_genes_really_do_have_truncating_variants(counted: dict) -> None:
    """The mirror image, and the reason the rule can tell them apart.

    *SCN1A* and *DMD* cause disease by there being less product. Truncating
    variants are the commonest way that happens, and a release where they were
    absent would mean the count is not measuring what this package thinks.
    """
    for gene in ("SCN1A", "DMD"):
        if gene not in counted:
            continue
        assert counted[gene].pathogenic_truncating_total > 0, (
            f"{gene}: no truncating variants counted, which contradicts the mechanism "
            "the reference set asserts for it"
        )
