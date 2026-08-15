"""Watching a mechanism call drift as the world underneath it changes.

The genome does not change. What changes is everything the rules consulted while
reading it: the ClinVar release, the population frequencies, the gene curation,
the transcript annotation, the patient's phenotype — and the rule files
themselves.

So a mechanism call is not a fact about a variant. It is a fact about a variant
*and a world*, and this package is what makes that difference usable. Name the
world and two things follow. A call becomes reproducible: re-pin, re-run, expect
the same answer. And the difference between two calls becomes attributable,
because the worlds differ in a small, enumerable number of coordinates.

The distinction the whole package is built to protect: "ClinVar learned
something" and "we changed our own rules" produce identical-looking
reclassifications, and must never be reported alike.
"""
