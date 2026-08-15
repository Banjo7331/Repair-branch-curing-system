"""Where the gene-level facts come from.

Until this package existed, every fact the rules read about a gene — the dosage
curation, the constraint, whether the product multimerises — was typed into a
fixture by hand. "Twenty-one reference cases reproduced" then meant only that
the rules were consistent with what somebody had already fed them.

This closes half of that, and is explicit about which half.

**Two of the facts are published as tables and are ingested here.** ClinGen
curates dosage sensitivity; gnomAD publishes constraint. Both are files with
release identities, both get parsed and digested, and the digest becomes the pin
that a mechanism call cites.

**Two of them are not published as tables at all.** Whether a gene product
assembles into a complex, and whether its null alleles are milder than its
missense ones, are judgements read out of the literature. Nobody ships them as a
TSV, and pretending otherwise would be the dishonest move — so they live in
``curation.py``, a local file that demands a citation per entry and is pinned
exactly like the downloaded sources. The difference is visible in every report:
a fact from ClinGen says ClinGen, and a fact we decided says we decided.
"""
