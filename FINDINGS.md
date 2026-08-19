# What broke, and what it changed

This project's rules, parsers and designers were written against fixtures, and
fixtures are written by whoever is about to read them. Everything below is what
happened when the same code met data and molecules produced by other people:
public releases, a real genome, and three published therapeutic oligonucleotides.

Not one of the defects below was found by the test suite. Every one of them had
passing tests over it, and two of the three designers were returning molecules
that could not work — one identical to its own target, one encoding a deletion
nobody asked for.

That is not a complaint about the tests. It is what tests cannot do: a fixture's
expected value is written by the same person who is about to assert it, so a
suite can only ever confirm that the code agrees with its author. Disagreement
has to come from outside.

They are kept because the corrections are the argument. A green badge says the
code agrees with itself; this says where it did not agree with the world, and
what was changed as a result.

---

## The public releases

The reference set used to assert its own transcript structures, with exon lengths made uniform. It now carries the **real** ones, read from RefSeq's GRCh38 annotation (release 2024-08) for the transcript this package resolves for each gene. Here is what that move turned up, because the finding is the point and a green badge afterwards is not.

**Every mechanism call survived.** All nine real-gene cases resolve to the same mechanism at the same confidence against real exon structures. The NMD predictions are unchanged too — every predicted-null case lands on the same side of the junction boundary it did before.

**But three of them had the variant in the wrong exon.** *SCN1A* c.2000 was exon 9 under uniform exons and is exon 11 under real ones; *COL1A1* c.1000 moved from exon 12 to 15; *UBE3A* from exon 5 of 12 to exon 4 of 11. The NMD answer happened to survive because the distance to the last junction stayed on the same side of a 50-nucleotide line. It was right for a reason the fixture could not have guaranteed.

**One case was asserting the opposite of what its own structure said.** The *DMD* case claims the affected exon is not a multiple of three, which is why skipping it alone would not restore the frame — and under real lengths c.4000 sits in exon 29, which is 150 nucleotides and *does* preserve the frame. The variant is now at c.7400 in exon 51: 233 nucleotides, out of frame, and the exon the first skipping oligonucleotide was aimed at. The case now tests what it says it tests.

**Four of ten modality cases moved, all on the same rule.** Exon skipping is the one modality that reads exon-level detail, and it is the one that changed: *UBE3A* lost the indication (real exon 4 is 1247 nucleotides, so skipping it shifts the frame), while *SMN1* and *SCN1A* gained it. Nothing else moved — which is the reassuring half of the result: the rules that read gene-level evidence do not care what the exons look like, and the one rule that does was the one getting it wrong.

**And the transcript chooser was quietly returning a prediction.** This is the finding worth the download on its own. *SMN1* — the flagship case of the whole modality set — has no MANE Select tag in RefSeq; its curated `NM_000344.4` sits on an unplaced scaffold, because the SMA locus is duplicated; and the assembled chromosome 5 carries only computed `XM_` models. The old rule was "take the longest coding sequence", so it returned `XM_054329962.1`: a transcript no human ever curated, with nothing in the output saying so.

The fix is an explicit ladder, and every rung is reported in the selection reason:

```
SMN1  NM_000344.4  8 exons, 885 nt
      no MANE Select transcript is annotated; curated over 33 modelled
      transcript(s); on NT_187651.1, which is not an assembled chromosome;
      longest coding sequence (885 nt of 3)
```

Curated over modelled first, because `NM_`/`NR_` were read by a human and `XM_`/`XR_` are a pipeline's guess. Then the assembled chromosome over a scaffold — but only *within* that choice, because for a duplicated locus the curated transcript on a scaffold is still the one the literature is about. Then length, to break what is left. A synthetic fixture cannot produce this situation, which is exactly why it survived until real data arrived.

Two smaller corrections came with it: *UBE3A*'s MANE Select transcript is `NM_130839.5` and not `NM_130838.4`, and *MECP2*'s is `NM_001110792.2` — the e1 isoform, three coding exons — rather than the e2 `NM_004992.4` the fixture named.

### The context files

**Two parsers were reading the wrong line.** ClinGen's gene curation list opens with five lines of provenance and then a header that is *itself* commented — `#Gene Symbol`, tab-separated. The reader skipped every `#` line, so the first line of provenance became the column names, and the failure surfaced as *"this file has no Gene Symbol column"* about a file whose first column is Gene Symbol. GTEx's GCT has the mirror-image problem: it declares itself on line one and gives its dimensions on line two, and that dimensions line is tab-separated and uncommented, so anything hunting for "the first line with tabs" reads `56200` as a gene name.

The rule that reads both is now explicit: while lines begin with `#`, keep the last one that has tabs in it, because a header has columns and a sentence does not; and a file that opens with the GCT magic has exactly two lines skipped, which is reading the format rather than guessing at it. Gzip is handled too — every one of these releases ships compressed and several ship only compressed. The shipped ClinGen fixture now carries a real preamble, so the shape is exercised by every test that touches it rather than only when somebody downloads 250 MB.

**Every dosage value the reference set asserts is what ClinGen actually says.** Seven of the nine genes are curated, and all seven match: *SCN1A*, *UBE3A*, *DMD*, *COL1A1*, *SCN2A* and *MECP2* at sufficient evidence for haploinsufficiency, *PIK3CA* at none — correct, it is an oncogene and its variants are gain of function. Those numbers were typed in by hand months before the file was downloaded, and they were right.

Two genes are absent from ClinGen's list entirely: *KRT14* and *SMN1*. The cases assert `no_evidence` for them, and that is now the one claim in the set that puts a word in ClinGen's mouth — "nobody has curated this" and "curated, and found nothing" are different, the model already distinguishes them, and the reference set does not yet.

**The expression numbers are real and they behave.** *KRT14* at 7304 TPM in skin, *DMD* at 23 in skeletal muscle, *SCN1A* at 8.9 in cortex and **0.006 in liver**. Pointed at a liver-directed therapy for the *SCN1A* case, the tissue rules do what they were written to do:

```
Brain - Cortex       8.866 TPM  ->  loss_of_function (probable)
                                    ruled out: allele_specific_silencing, gene_addition, ...
Liver                0.006 TPM  ->  loss_of_function (possible)
                                    caution: GENE_SILENT_IN_THE_AFFECTED_TISSUE
                                    ruled out: ..., wild_type_upregulation
```

Confidence drops, the caution names itself, and raising output from the intact allele is withdrawn — because there is nothing to raise where the gene is not transcribed. That path had only ever been exercised against numbers this project invented.

### The constraint file, and the schema that was not there

gnomAD's downloads page links the **v2.1.1** per-transcript file, not the v4 one this package's parser was written against — and the two releases spell the same two columns differently. v2.1.1 calls LOEUF `oe_lof_upper` and marks the row to use with `canonical`, Ensembl's own pick; v4 calls them `lof.oe_ci.upper` and `mane_select`. Pointed at the file the site actually serves, the parser failed with "no such column", which is the correct failure and a useless one to somebody holding the right file.

Both schemas are now declared and detected from the header — not from the file name, because a file renamed on the way to disk is ordinary and its columns are not — and the refusal names every spelling it looked for. Nothing matches a column by position or by resemblance: a constraint value read out of the wrong column is a number that looks entirely reasonable. `.bgz` is read too, which is what genomics ships and what `.gz`-only handling misses.

One difference is recorded rather than smoothed over. MANE Select is an agreement between RefSeq and Ensembl about the clinical reference transcript; Ensembl canonical is Ensembl's own. They usually agree and are not the same claim, so every value says which one it came through:

```
loeuf = 0.071  [gnomad_constraint@v2.1.1/d54d60c8f87e]  gnomAD v2.1.1, Ensembl canonical transcript
```

**The reference set now carries the real numbers, and every mechanism holds.** Some of the hand-typed values were well off — *PIK3CA* was asserted at 0.6 and is really 0.117, *MECP2* at 0.12 and is really 0.407, *SMN1* at 0.6 and is really **1.929** — and none of the calls moved.

*PIK3CA* is the one worth dwelling on. It is among the most loss-of-function-constrained genes in the genome, and its disease mechanism is gain of function; a system that read constraint as evidence for haploinsufficiency would have flipped it the moment the real number arrived. It did not, because the rule that reads LOEUF is scoped to predicted-null variants and marked *supporting* — constraint is a claim about a gene across a population, not about the variant in front of you. That scoping was written months before there was a real number to test it against, and this is the first evidence that it was the right call rather than a plausible one.

*SMN1* at 1.929 is the mirror image: the gene behind spinal muscular atrophy is, by this measure, unconstrained. Which is correct and easy to misread — SMA is recessive, carriers are common, and the duplicated locus makes calling hard. A pipeline treating high LOEUF as "this gene tolerates loss" would draw exactly the wrong conclusion about the best-known recessive disease in the list.

### One real variant, all the way through

Everything above validates the *inputs*. The designers had still never seen a real chromosome: base editing, pegRNAs and antisense tiling had only ever run against an 800-nucleotide repeated trimer with PAMs placed by hand at known distances — a fixture built so the answer would be knowable, which is exactly why it cannot test this.

So: *COL1A1* p.(Gly821Ser), c.2461G>A, on real GRCh38. The case is chosen for three properties. The gene is on the **minus strand**, so every coordinate, complement and exon ordering runs in the direction where mistakes hide. The change is G>A on the coding strand and therefore **C>T on the plus strand**, which is where a designer reasoning in genomic coordinates without thinking about strand produces a confident wrong molecule. And it is the textbook dominant negative.

It runs, and it verifies from the outside:

```
$ repairbench plan col1a1-gly821ser.yaml --fasta chr17.fa --annotation GRCh38_latest_genomic.gff.gz

# annotation  GRCh38_latest_genomic.gff.gz@b13607bdad3a
# transcript  NM_000088.4 (MANE Select), 51 coding exons, 4395 nt; variant at c.2461

COL1A1  NM_000088.4
  mechanism   dominant_negative (probable)
  not designed, ruled out by the modality rules
    gene_addition, truncated_construct, wild_type_upregulation, ...
  designed
    allele_specific_silencing → aso     AGCCAACCTGGTGCTAAAGG  chr17:50190080-50190099
    base_editing → base_editor          GACAGCCAACCTGGTGCTAA agg  chr17:50190083-50190102 (-)
```

Four things checked by hand against the reference, none of them by the code that produced them:

* The genomic coordinate maps to **c.2461**, which is what the literature calls this variant. The case file carried both numbers, written independently; on a minus-strand gene they agree only if the exon reversal and the CDS arithmetic are both right.
* The codon at that coordinate reads **GGC** on the coding strand — a glycine. The reference set has been calling this a glycine substitution since before there was a genome to check it against.
* The base editor chose the **minus strand**, because C→T on the plus strand is A→G on the minus, and only an adenine editor can make it. Read the span out of the reference, substitute the patient's allele, complement it, and the guide comes out identical to the reported one.
* The antisense oligonucleotide is the reverse complement of the patient's window, not the reference's — the distinction that decides whether it silences the affected allele or the healthy one.

The chromosome-naming collision turned out to be the least of it: the case says `chr17`, RefSeq says `NC_000017.11`, and the aliases read out of the file's own region records reconciled them without anybody noticing.

**What did break was the parser, on a gene nobody would think to invent.** *PEG10* is a retrotransposon-derived gene translated through a programmed ribosomal frameshift: the ribosome slips back one base and reads on in another frame, so two of its CDS blocks legitimately share a coordinate. RefSeq annotates that correctly. This package refuses it — a position inside the overlap has two CDS offsets and every calculation downstream assumes one — and *the refusal took down the entire file*. Three transcripts of one gene were denying access to the other **136,249**.

They are now dropped individually, with the reason kept and readable, because "this transcript is unusable" and "this file is broken" are different sentences and only one of them is true.

The whole thing is a test — `tests/test_real_locus.py`, with the verified values written down — that skips when `refdata/` is absent, so CI runs the other 368 and this one runs where the data is.

### A default that was quietly making a claim

The last thing the real ClinGen file turned up was in the model rather than in a parser. `DosageScore` had six values and the *default* was `no_evidence` — so every gene nobody had curated came out saying, in the output and in the provenance, that ClinGen had looked and found nothing.

ClinGen has evaluated a few thousand genes. There are twenty thousand. `KRT14` and `SMN1` are two of the ones it has not, and both were carrying a finding it never made.

There are now three absences, and they are different claims:

| | meaning |
| --- | --- |
| `not_evaluated` | nobody has looked — and it is the default, because for most of the genome that is true |
| `no_evidence` | ClinGen curated the gene and found nothing to support dosage sensitivity yet |
| `dosage_sensitivity_unlikely` | ClinGen looked and concluded it does not apply — a finding, not a gap |

Nothing in the rules changed and no reference case moved, which is the point: the correction is about what the system *says it knows*, not about what it concludes. A feature — `gene.haploinsufficiency.curated` — is now available to any rule that wants to hedge on an uncurated gene, and none does yet; putting that choice in the rule file rather than in code is where it belongs.

The *SMN1* case had the same problem in sharper form. It asserted `autosomal_recessive`, which is ClinGen's dosage code 30, for a gene ClinGen has no dosage row for at all. SMA *is* recessive — that is not in doubt — but the claim was hanging on the wrong hook. It now hangs on the two that hold it honestly: the gene-disease validity curation, which does exist for *SMN1*, and the patient's zygosity, which is what actually rules the modalities out.

### The last invented input, which cost two reference cases

Everything above was a parser meeting a file that was shaped differently than expected. This one is different: the file was read correctly and said the rules were wrong.

`distribution:` — how many pathogenic missense variants a gene has, how many cluster, how many truncating — was the last input the reference set was still inventing. It is also the input to the least direct thing the package attempts: telling gain of function from loss of function by noticing that missense variants pile into one stretch of protein while truncating variants stay away. Nine genes, counted from ClinVar at ≥1★:

| gene | missense | truncating | in densest 20 aa | truncating share | mechanism in the set |
| --- | ---: | ---: | ---: | ---: | --- |
| PIK3CA | 76 | 4 | 18% | **5%** | gain of function |
| KRT14 | 42 | 21 | **52%** | 33% | dominant negative |
| SCN2A | 369 | 208 | 7% | 36% | gain of function |
| SCN1A | 994 | 862 | 4% | 46% | loss of function |
| COL1A1 | 446 | 657 | 6% | 60% | loss of function / dominant negative |
| SMN1 | 18 | 41 | 44% | 69% | loss of function |
| UBE3A | 40 | 159 | 12% | 80% | loss of function |
| MECP2 | 136 | 612 | 25% | 82% | loss of function |
| DMD | 24 | 1806 | 12% | 99% | loss of function |

**Read the clustering column downward.** The most tightly clustered gene in the set is *KRT14*, a textbook dominant negative. Then *SMN1* and *MECP2*, both loss of function. Both gain-of-function genes are near the bottom. The rule that demanded 60% clustering to call gain of function was, on real data, a rule that would have preferred the dominant-negative gene.

Two things make the measurement blind, and both are properties of the file rather than of the biology.

**ClinVar counts alleles, not patients.** "Hotspot" is a claim about recurrence — the same residue, over and over, in unrelated people. The variant summary has one row per distinct allele, so p.His1047Arg contributes exactly as much as a variant seen once in one family, and a long tail of private missense outvotes the hotspot every time. The signal is not entirely gone — weighting by the number of submitting laboratories puts *PIK3CA* residues 1047, 545, 546 and 542 at the top, which are precisely the canonical hotspots — but counting rows deletes it.

**A hotspot is rarely one window.** The rule's own citation named "the voltage sensors *and* the pore" — four domains, hundreds of residues apart. The densest single window catches one of them. Widening to the five densest non-overlapping windows and weighting by submitters still only reaches 31% for *SCN2A*.

**What actually separates the genes is the last column.** *PIK3CA* sits at 5% truncating; the next gene up sits at 33%. That is the only clean separation in the table, and it is the claim the rule was really making all along — truncating variants are not a disease mechanism in this gene. But the rule that read it demanded *zero* truncating variants, and real *PIK3CA* has four. Any database of a hundred submissions holds a few rows that disagree with the rest, and a boolean hands each of them a veto.

So `mechanism-v1` changed in two ways. `NO_PATHOGENIC_TRUNCATING_REPORTED` became `TRUNCATING_VARIATION_DEPLETED`, reading a share against a threshold in the rule file rather than a yes/no, and carrying the *strong* weight. `PATHOGENIC_MISSENSE_CLUSTERING` was demoted from *strong* to *supporting*: it still fires, it still argues for gain of function, and it can no longer reach that conclusion on its own.

Two reference cases moved, and neither was tuned back:

* ***PIK3CA*** dropped from `probable` to `possible`. Probable used to come from two rules firing — clustering and no-truncating — and **those two rules read the same dataset**. Two views of one distribution are not two pieces of evidence; scoring them as though they were is how a package talks itself into certainty. One strong rule now fires, and possible is what one rule buys.
* ***SCN2A*** stopped resolving. The case always said the gene carries two mechanisms; what the real numbers show is that a gene-level distribution is the *sum* of both, and reads as neither — 36% truncating, 7% clustering. The variant in the case, p.Arg853Gln, genuinely is gain of function, but nothing the package can see says so. Answering the question would need variant-level curation or a functional assay. `undetermined` is now the expected result, and the case explains why.

The reference set carries real CDS positions too, on the transcript each case pins: *COL1A1* c.2461G>A (p.Gly821Ser), *MECP2* c.916C>T at 3★, *DMD* c.4088del inside exon 30 — 162 nucleotides, still a multiple of three, so the exon-skipping case keeps the property it was built to test. Every input in the reference set now came from somewhere else.

`tests/test_real_clinvar.py` runs the same cases against the file itself where it exists, so the copies in the YAML cannot drift from the release without something going red.


## The published molecules

Every other check in this repository compares the package with a *statement* — the literature's conclusion, a curated release, its own reference set. This one compares it with an *object*: a molecule that exists, whose sequence is printed on an FDA label, which was manufactured and given to patients.

**Eteplirsen** (EXONDYS 51) is a thirty-nucleotide morpholino that makes the spliceosome skip *DMD* exon 51. Its thirty bases are quoted verbatim in the label. The question is exact and has no room for interpretation: point this package at *DMD* exon 51 in GRCh38 and does the approved molecule come out?

**The first run said no, and the reason was worse than a near miss.**

`tile` reverse-complemented every window unconditionally. For a plus-strand gene that is right — the messenger is the forward sequence, so an oligonucleotide complementary to it is the reverse complement. For a **minus-strand** gene it is exactly backwards: the messenger *is* the reverse complement, so the oligonucleotide must carry the forward sequence. What the module returned for such a gene was the sequence of the transcript itself — a molecule identical to its own target, which hybridises with nothing.

*DMD*, *COL1A1*, *MECP2*: minus strand, most of the reference set. Every antisense oligonucleotide this package had ever printed for them was inert by construction.

Nothing caught it, and the two reasons are worth naming:

* **The synthetic fixtures have no orientation.** A made-up sequence on a made-up contig has no strand to get wrong, so the arithmetic was tested and the biology was not.
* **The one value that would have caught it was written down and never asserted.** `tests/test_real_locus.py` had `OLIGONUCLEOTIDE = "AGCCAACCTGGTGCTAAAGG"` sitting in a block headed *verified by hand* — and no test referenced it. It was the reverse complement, copied out of this function's own output, so even the hand verification had checked the code against itself.

The fix gives `tile` a required `strand` with no default, because the defect was precisely a silent assumption of one. A plan without an annotation now refuses to design an antisense oligonucleotide at all rather than guessing the orientation.

With that, the molecule comes out base for base:

```
$ repairbench aso --gene DMD --at chrX:31773960-31774192 --chemistry steric-PMO-30 \
      --strand - --fasta refdata/chrX.fa

CTCCAACATCAAGGAAGATGGCATTTCTAG   chrX:31774098-31774127   ← eteplirsen, FDA label 206488
```

Two smaller things fell out of the same exercise, and both are corrections to the rule file rather than to code:

* **PMO length is not a property of the chemistry.** One entry declared 25 nucleotides; the approved exon-skipping morpholinos are 30, 25, 22 and 21. A catalogue that could not express eteplirsen's length could not have reproduced eteplirsen whatever the strand logic did.
* **`tm_max_c: 75` was a threshold no rule read.** Writing the missing ceiling would have made things worse: the Wallace approximation returns about 86 °C for a thirty-mer, which is not a melting temperature, so the ceiling would have flagged an approved drug for a property this file cannot measure. The threshold is gone and a rule now says when the number is outside the approximation's range.

### The base editor, where the disagreement was the interesting part

The same exercise for M7's other designer, on a real disease variant: ***FAH* c.1062+5G>A**, hereditary tyrosinemia type 1, the fifth base of a splice donor. A published worked example gives two guides — one that installs the variant in a cell line with a cytosine editor, one that corrects it with ABE7.10 — and prints both sequences.

Handed nothing but the patient's allele and chromosome 15, this package independently produced the correction: **ABE7.10, PAM AGG, plus strand, target at protospacer position 5** — every structural property of the published guide, none of it supplied.

And the protospacer differed from the printed one at **one base out of twenty**.

That single base is the result worth the whole exercise. Position 3 of our guide reads G, where the paper reads A. Neither is a mistake. The authors' correction guide was written against their *cell line*, in which the variant had been installed by a cytosine editor — and a cytosine editor's window covers more than the base it was aimed at. Their allele carries the variant **plus a bystander**; a patient carries only the variant.

The check that turns that from a story into a finding: run the *installing* guide through this package and ask what bystanders it predicts.

```
BE4max-SpCas9  GATACTCACCGGCCCGCTGA tgg  strand -  position 5
bystanders:
   position 7  g.80180228  →  T
```

`g.80180228` is c.1062+3 — exactly the base by which the published correction guide differs from ours. The package predicted the discrepancy before anyone knew what it was, and both guides land where the paper says they land: position 5, PAM TGG on the minus strand for the installer, PAM AGG on the plus strand for the corrector, all confirmed against GRCh38 rather than taken on trust.

### Prime editing, and the worst defect of the three

The third designer, against the paper that introduced the method: Anzalone et al. 2019, and its pegRNA for the Ashkenazi Tay-Sachs allele ***HEXA* c.1274_1277dupTATC**. The supplementary tables give the spacer, the 3′ extension, and the PBS and template lengths separately.

Handed a patient allele and chromosome 15, this package produced:

```
spacer ATCCTTCCAGTCAGGGCCAT   PAM AGG   strand +   nick g.72346574, 5 nt from the edit
PBS  10 nt  GCCCTGACTG
RTT  21 nt  ACCTGAACCGTATATCCTATG
3' extension  ACCTGAACCGTATATCCTATGGCCCTGACTG      ← identical to the published pegRNA
```

— and offered the paper's own second nick, `TACCTGAACCGTATATCCTA`, classified **PE3b**, which is correct and not obvious: that guide is the *installing* pegRNA's spacer, so it matches the corrected allele and cannot fire until the edit is made.

Getting there took two corrections, and the first is the most serious defect this project has had.

**The patient's sequence was built by deleting the wrong number of bases.** Both the patient and the edited sequence spliced `len(patient_allele)` bases out of the reference, where they had to splice out `len(wild_type_allele)` — what the reference actually carries there. For a substitution the two are equal and nothing shows. For an insertion or a deletion — *which is the entire reason this module exists rather than the base editor* — both came out wrong: a four-base insertion silently consumed four reference bases, so the "patient" sequence was the reference untouched, and the "edited" sequence carried a four-base deletion nobody asked for. Every pegRNA the module produced for an indel encoded a template that writes the wrong product.

It survived a full module of tests because the fixtures exercise substitutions, where the bug is invisible by construction. The published molecule found it in one run, and found it *precisely*: the primer binding site matched the paper exactly — it is read upstream of the nick, where the bug does not reach — while the template encoded something else. A partial match is a much sharper diagnostic than a total mismatch.

**And the design space could not express the molecule.** The module scanned primer binding site length across its whole range but emitted exactly one template length: the shortest that carries the edit plus the minimum homology arm. The published pegRNA uses a 16-nucleotide arm, and the three HEXA pegRNAs in that table alone use templates of 14, 21 and 27 nucleotides. Template length is one of the two parameters people vary at the bench, and fixing it to a constant meant the molecule somebody made was not in our design space at all. Both parameters are scanned now, and the candidate count grew from 40 to 600 — which the module accepts, because it ranks none of them and says so.

**What this does and does not establish.** All three designers now agree with molecules that exist: strand handling, coordinate mapping, PAM scanning, nick placement, window arithmetic, bystander prediction, primer binding site and template derivation, PE3b classification and the chemistry catalogue all reproduce published objects at real loci.

Two of the three were **wrong** when first asked, and neither defect was subtle once seen: the antisense module returned molecules identical to their own targets for every minus-strand gene, and the prime module built the patient's sequence by deleting the wrong bases for every insertion and deletion. Both had full test modules that passed. What the tests could not do is disagree with the package about what the answer *is* — a fixture's expected value is written by whoever writes the fixture, and against real molecules that stops being true.

It establishes nothing about efficiency. The rule file cautions about eteplirsen, correctly, for binding an exon interior rather than a splice site, and the drug works anyway; that is why the caution is a caution and why nothing here ranks candidates.


## Before any of that: zygosity, and the answer that was wrong

The first version of this package had no concept of zygosity, and its absence produced a *wrong* answer rather than a missing one: `wild_type_upregulation` came back indicated for the Duchenne case. Duchenne is X-linked, the boy has one X, and it carries the variant. There is no intact allele to raise.

That is the shape of mistake the model now guards against. A mechanism is a property of a variant and a gene; how much functional product the patient has left is a property of the **patient**, and roughly half the interventions in M6 depend on the second rather than the first. Three contraindications follow from it — no upregulation without an intact allele, no allele-specific silencing of the only copy, and no X reactivation for a hemizygous male — and a test asserts across the whole reference set that no modality requiring a wild-type allele is ever offered to a patient without one.

Unknown zygosity is deliberately not treated as either answer. Reading it as "no wild-type allele" would rule out real options on missing data; reading it as "yes" would offer options that may not exist. The verdicts stand and the selection carries a caveat naming exactly which of them depend on the gap:

```
read with caution
  · zygosity was not supplied, and these modalities only work if the patient has an
    unaffected copy of the gene: wild_type_upregulation. Their verdicts above assume
    one exists and must not be read as established until it is confirmed.
```

Which modalities need an intact allele is declared in the rule file, not in code, because it is a clinical claim like any other.



## The review surface

Building the reanalysis queue as a page turned up three defects in the code
beneath it, which is the usual result of making something visible.

**Urgency was being lost on every save.** A ledger written back to disk recorded
its events' urgency and queue as `"-"`, so a case loaded and saved twice forgot
how loudly it had been asking. Nothing read those fields until a queue needed to
be listed, and by then the information had been quietly discarded for months.

**`case_ids()` returned things that were not cases.** The command line writes
`<case>.variants.json` beside each ledger, and a `*.json` glob read that back as
a case named `NICU-014.variants`. The first thing the dashboard ever did was
crash on it. Cases are now identified by shape — a ledger is an object with a
`case_id` — rather than by file name.

**The ledger could not tell "quiet" from "never ran".** It recorded when a case
last *changed*, never when it was last *examined*, so the first version of the
page reported every healthy silent case as dead. That distinction is the whole
reason the page exists: an empty queue is the commonest correct output this
system produces, and it is also exactly what a scheduler that died in March
looks like. `last_examined_at` is now written on every run, especially the ones
that find nothing.

**And registration was not taking a baseline.** It recorded which variants a
case was watching and nothing about what we currently concluded, so the first
run after registration was always a no-op — every variant was being seen for the
first time, and that run's answer silently became the reference. Register a case
in January, run it in April, and February's curation change goes unreported.
From a terminal the first run is one of many and the gap is invisible; from a
page it is the button somebody presses expecting an answer.
