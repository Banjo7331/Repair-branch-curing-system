# Pointing this at real data

Everything in `tests/data/` is synthetic and labelled as such. It tests the
arithmetic, which does not care whether the coordinates correspond to a real
locus. What it cannot test is whether the rules survive contact with real
transcript structures, real curations and real constraint values — and that is
the validation this project's README has been promising in three places.

```bash
./scripts/fetch-reference-data.sh
```

About 500 MB of downloads, roughly 1.5 GB on disk once unpacked, into
`refdata/` — which is gitignored, because none of these files is ours to
redistribute.

## What arrives, and what each one is for

| file | ~size | what it settles |
| --- | --- | --- |
| `GRCh38_latest_genomic.gff` | 1.2 GB | Real exon structures, which the reference set now carries — they replaced uniform invented ones, and doing so exposed a case asserting the opposite of what its own structure said. |
| `chr{2,3,5,15,17,X}.fa` | 900 MB | Real sequence for the genes in the reference set. Every designer reads bases; so far it has read invented ones. |
| `ClinGen_gene_curation_list_GRCh38.tsv` | 2 MB | Dosage sensitivity as ClinGen actually scores it, including the 30 and 40 codes the rules treat as findings rather than absences. |
| `gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz` | 13 MB | Real LOEUF, per transcript, with the canonical row that has to be picked correctly. The downloads page links v2.1.1, whose column names are not v4's — which is why the parser declares both schemas and detects one from the header. |
| `GTEx_..._gene_median_tpm.gct.gz` | 15 MB | Real expression. Note the GCT preamble — two lines before the header, which the parser refuses to guess at. |
| `variant_summary.txt.gz` | 250 MB | Where pathogenic variation sits in each gene: the last input the reference set was still inventing, and the one feeding the least direct rule in the package. |

## The one wrinkle worth knowing before you start

**The annotation and the sequence use different chromosome names.** NCBI's GFF3
says `NC_000017.11`; UCSC's FASTA says `chr17`. They are the same sequence, and
nothing in this package silently reconciles them — a coordinate looked up under
the wrong name either fails loudly or, worse, succeeds against the wrong
contig.

Both conventions are downloaded on purpose rather than by accident. The
alternative is NCBI's own 3 GB genome FASTA, which matches its GFF3 and costs an
hour of downloading for sequence that is identical base for base.

## What to expect

Some of the reference set will probably not survive this. That is the point of
running it: a case that passes against uniform 141-nucleotide fixture exons and
fails against the real structure was passing for the wrong reason, and finding
that out is worth more than a green badge.

That has already happened five times — a transcript chooser returning a computed
prediction, two parsers reading a preamble line as a header, a schema that was
not the one the downloads page links, and one gene whose real biology took down
the whole annotation file. The README's "What happened when it was pointed at
real data" section is the record. None of them could have been found against a
fixture, because a fixture is built by the same person who is about to read it.
