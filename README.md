# repairbench

**Two nonsense variants in the same gene can have opposite mechanisms and need opposite therapies. This works out which one you have, which interventions that admits — and then keeps watching, because the answer can change without anybody learning anything.**

Every clinical judgement lives in a rule file — `mechanism-v1`, `modality-v1`, `editors-v1`, `prime-v1`, `aso-v1`, `offtarget-v1`, `routing-v1`, all under `rules/` — and not in the code. If a conclusion here is wrong, the error is in one of those files, in a sentence a geneticist can read. That is the point.

The layers, one package:

| | |
| --- | --- |
| `vcf.py` | the patient's own file in; carried alleles, zygosity and read consequences out |
| `context/` | ClinGen dosage, gnomAD constraint and GTEx expression in; gene-level facts out, with a citation on every one |
| `annotation/` | GFF3 and reference FASTA in; real transcripts, verified coordinates, left-aligned indels out |
| mechanism (M5) | why this variant causes disease — loss of function, gain of function, dominant negative — with its evidence |
| modality (M6) | which classes of intervention that mechanism admits, and which it rules out |
| design (M7) | the molecule: a base editor's protospacer, a pegRNA, or a tiling of antisense oligonucleotides — each with what is wrong with it |
| `plan.py` | the seam: one case from mechanism to molecules, where a ruled-out modality is never designed |
| `reanalysis/` | watching all of the above drift as releases land, and proving which release moved it |

![Python](https://img.shields.io/badge/python-3.11%2B-3776AB)
![Rules](https://img.shields.io/badge/rules-as%20data-6f42c1)
![Tests](https://img.shields.io/badge/tests-333-brightgreen)
![Reference set](https://img.shields.io/badge/reference%20cases-21%2F21-brightgreen)

> ⚠️ **Not a medical device.** Research and educational use. A candidate here is a *design* — a protospacer, a pegRNA, an oligonucleotide, and what is wrong with each — not a therapy. No model of editing efficiency, prime-editing yield or target accessibility is attached to any of them, and "not ruled out" is not a recommendation.

---

## Why this and not something more exciting

The obvious first module to build is the one that designs sequences. It is the wrong one to build first, because a sequence designed against the wrong mechanism is worse than no sequence at all.

Take a nonsense variant in *COL1A1*. If the premature stop sits far enough upstream of the last exon-exon junction, the transcript is destroyed, half the normal collagen is made, and the child has the mild form of osteogenesis imperfecta. The therapeutic question is how to supply more product.

Take a glycine substitution in the same gene. The protein is made, it is incorporated into the triple helix, and it ruins it. The phenotype is *worse* than having no protein at all — and supplying more normal collagen does not remove the defective chains. The therapeutic question inverts: how to get rid of the mutant allele.

Same gene. Opposite intervention. A module that designs correction sequences without settling this first will, half the time, confidently design the one that makes things worse.

## What it does

```
$ repairbench explain col1a1-gly821ser.yaml

COL1A1  NM_000088.4
  mechanism   dominant_negative  (probable)
  because
    NULL_ALLELES_ARE_MILDER [strong]
      Carriers of null alleles in this gene are milder than carriers of in-frame
      variants. If having none of the protein is better than having a broken one,
      the broken one is doing harm rather than being absent.
      — Classic contrast in COL1A1 — haploinsufficiency gives mild osteogenesis
        imperfecta type I, glycine substitutions in the triple helix give the severe forms.
    MISSENSE_IN_MULTIMER [supporting]
      A missense product that still assembles into its complex can interfere with it.
  not ruled out on mechanistic grounds
    no   gene addition addresses the mechanism
    yes  coding sequence fits a viral payload
    no   silenced allele available to reactivate
    yes  allele-specific silencing indicated
  notes
    · the mutant allele is doing harm rather than being absent, so supplying a
      normal copy does not address the mechanism and may worsen it
  ruleset     mechanism-v1@7dd6bd311e2b
```

## The architecture is the rule file

Everything clinical is declarative:

```yaml
- id: NMD_ESCAPE_IN_MULTIMER
  supports: dominant_negative
  strength: strong
  when:
    all:
      - {feature: consequence.is_predicted_null, is: true}
      - {feature: nmd.outcome, eq: escape}
      - {feature: gene.forms_multimer, is: true}
  because: >
    The stop codon escapes decay, so a truncated protein is made, and the product
    assembles into a complex where a defective subunit can poison the whole.
    Supplying a normal copy does not remove the poisoning subunit.
  citation: ClinGen SVI PVS1 decision tree; Veitia 2007, Trends Genet
```

Three properties follow, and each was the reason for a design decision:

**A geneticist can review it without reading Python.** The predicate language is *interpreted, not executed* — no `eval`, no expression strings, no way for a rule to reach anything the feature record does not name. A rule file is a document you can hand to a clinician without also handing them code they are running.

**A typo fails the run.** A rule naming a feature that does not exist raises, listing the features that do. The alternative — evaluating to false — produces a rule that silently never fires, which is indistinguishable from a rule that correctly never fires, and only one of those is a bug.

**Every call names the rules that produced it.** `mechanism-v1@7dd6bd311e2b` is the version plus the digest of the file's bytes. Two calls made under different rules are not comparable, and a system that cannot tell them apart will eventually report a rule change as a biological finding. This is the same discipline the reanalysis service keeps about knowledge snapshots, applied to our own judgement.

Thresholds are rule data too — the NMD junction boundary is quoted as 50–55 nt across the literature, and a laboratory preferring 55 should be able to say so in the file rather than in a patch.

## Three things it refuses to do

**It will not return a mechanism without evidence.** There is no constructor for one. A determined call always carries the rules that produced it, and a test asserts that across the whole reference set.

**It will not hide a disagreement.** A curated mechanism from an expert panel outranks anything inferred — but the other rules still run, and any that disagree are recorded as conflicts. *MECP2* is the case that forces this: most of its coding sequence sits in the last exon, so common truncating variants escape decay and a protein *is* made, yet the curated mechanism is loss of function. The call follows the curation and reports the escape, because a reviewer needs both facts.

**It will say it does not know.** `undetermined` is a first-class answer with its own rules — "the dosage curation was refuted, not merely absent", "a truncated protein is produced and nothing here says whether it is inert or harmful", "four reported missense variants is not a distribution". A wrong mechanism points at a wrong therapy; an honest hedge costs nothing.

## Validation

The reference set is the specification; the rule file is the implementation.

```
$ repairbench reference

ok    SCN1A nonsense, Dravet syndrome                  loss_of_function (probable)
ok    UBE3A nonsense, Angelman syndrome                loss_of_function (probable)
ok    DMD frameshift, Duchenne muscular dystrophy      loss_of_function (probable)
ok    COL1A1 nonsense, OI type I                       loss_of_function (probable)
ok    COL1A1 glycine substitution, OI type II-IV       dominant_negative (probable)
ok    KRT14 missense, epidermolysis bullosa simplex    dominant_negative (probable)
ok    PIK3CA hotspot missense                          gain_of_function (probable)
ok    SCN2A missense, both mechanisms in one gene      gain_of_function (possible)
ok    MECP2 nonsense in the final exon, Rett           loss_of_function (established)
ok    dosage sensitivity refuted                       undetermined (none)
ok    too few missense variants to see a pattern       undetermined (none)

11/11 reference cases reproduced under mechanism-v1@7dd6bd311e2b
```

Two of these deserve a note.

**SCN2A resolves to *possible*, and that is the correct answer.** Gain-of-function variants cluster in the channel's voltage sensors and cause early epilepsy; truncating variants cause a later, different phenotype. A resolver that reported one confident mechanism for this gene would be wrong about half the patients.

**The reference set already caught a rule that over-claimed.** The clustering rule originally read concentration of pathogenic missense variants as gain of function — and so called the *COL1A1* glycine substitutions, the textbook dominant negative, gain of function. Clustering says *this domain is critical*, which is not the same as *this variant confers a new activity*. The rule is now gated on there being no dominant-negative explanation available, and the gate is commented with the case that forced it. This is what the validation is for, and it is worth more than a passing test suite.

## Running it

```bash
repairbench plan case.yaml --fasta ref.fa --annotation ref.gff3   # the whole thing
make test        # 333 tests
make reference   # re-run both reference sets
make rules       # print what the system believes, with citations
repairbench assess case.yaml   # mechanism, then the modalities it admits
```

## M6 — which interventions the mechanism admits

M5 answers *why*. M6 answers *what class of intervention that mechanism allows*, and stops there.

**A contraindication outranks any number of indications.** Reasons to try something accumulate; one reason not to ends the matter. The report prints the ruled-out list first, in defiance of what a reader would rather see, because the output that changes what somebody does is the one that closes a door.

**An unresolved mechanism blocks everything.** If M5 returned `undetermined`, no modality is assessed at all — the alternative is rules firing on transcript facts alone and producing a plausible-looking list of interventions resting on nothing.

The modality reference set is diseases where a route was actually taken:

```
$ repairbench reference --modalities

ok    SMN1 loss of function, spinal muscular atrophy
        loss_of_function → gene_addition, base_editing, prime_editing
ok    DMD frameshift in a skippable exon, Duchenne muscular dystrophy
        loss_of_function → truncated_construct, exon_skipping
ok    UBE3A nonsense, Angelman syndrome
        loss_of_function → gene_addition, silenced_allele_reactivation, base_editing
ok    SCN1A nonsense, Dravet syndrome
        loss_of_function → wild_type_upregulation, truncated_construct
ok    MECP2 nonsense, Rett syndrome (heterozygous female)
        loss_of_function → silenced_allele_reactivation, base_editing, prime_editing
ok    MECP2 nonsense in a hemizygous male
        loss_of_function → base_editing, prime_editing
ok    COL1A1 glycine substitution, severe osteogenesis imperfecta
        dominant_negative → allele_specific_silencing, base_editing, prime_editing
ok    PIK3CA hotspot missense
        gain_of_function → allele_specific_silencing, base_editing, prime_editing
ok    Unresolved mechanism blocks every modality
        undetermined → everything blocked

10/10 modality reference cases reproduced under modality-v1@6e96b18194a3
```

Three of those are worth reading closely.

**SCN1A does not get gene addition, and that is right.** Dravet syndrome is haploinsufficiency, so replacement looks obvious — until you notice the coding sequence is six kilobases and does not fit in a vector. What the clinical programme actually does is raise output from the intact allele instead, and the module gets there by arithmetic rather than by having been told.

**MECP2 gets gene addition *contraindicated*, not merely unranked.** Loss of MECP2 causes Rett syndrome and duplication of MECP2 causes a syndrome of its own. Dosage runs in both directions, so an untitrated working copy trades one disease for another. A module that only asked "is this loss of function" would propose it.

**Writing the MECP2 case exposed a missing modality.** Affected girls are heterozygous and mosaic for X inactivation, so an intact wild-type copy is present and silent in about half their cells. That is not imprinting — but it is the same therapeutic opportunity, and it was the only route that case had. The model now covers both under `silenced_allele_reactivation`, and the enum says plainly why a taxonomic distinction was set aside.

## Context: which facts can be ingested, and which nobody publishes

Until `context/` existed, every fact the rules read about a gene was typed into a fixture by hand, and "twenty-one reference cases reproduced" meant only that the rules were consistent with what somebody had already fed them. This closes half of that, and the interesting part is being explicit about which half.

**Two of the four facts are published as tables.** ClinGen curates dosage sensitivity; gnomAD publishes constraint. Both are parsed, digested, and the digest becomes the pin a call cites — which is what turns `gene_curation@2026-01` from a label into something checkable.

**Two of them are not published at all.** Whether a gene product assembles into a complex, and whether its null alleles are milder than its missense ones, are judgements read out of the literature. Nobody ships them as a TSV. So they live in a local curation file that *demands a citation per entry* and is pinned like everything else — and the difference shows up in every line of output:

```
$ repairbench context COL1A1

COL1A1
  forms_multimer = True  [local_curation@rev3/bbf8f9fd50de]
      Type I collagen is a heterotrimer; a mutant pro-alpha1 chain is incorporated into the helix
  haploinsufficiency = sufficient_evidence  [clingen_dosage@2026-01/9ba328295043]
  loeuf = 0.25  [gnomad_constraint@v4.1/118d41eba7e8]
  truncating_variants_are_milder = True  [local_curation@rev3/bbf8f9fd50de]
      Marini et al. 2007, Hum Mutat — glycine substitutions in the triple helix versus
      haploinsufficiency in OI type I
```

A fact from ClinGen says ClinGen. A fact we decided says we decided.

Four refusals are worth naming, and one of them was found by its own test:

* **ClinGen's scale is not a scale.** 0–3 rank evidence, 30 means recessive, 40 means dosage sensitivity was *actively refuted*. Flattening 40 into "low" would erase the distinction the mechanism rules most need — a predicted null variant means something different under "nobody looked" than under "somebody looked and said no".
* **A gene ClinGen has not evaluated contributes nothing, not a default.** The first version created the provenance entry before parsing, so an unevaluated gene came out *present with nothing in it* — which reads downstream as "we have context for this gene". The test caught it.
* **Only the MANE Select constraint row is read.** The file has a row per transcript with different numbers; taking whichever came first would make the value depend on file ordering.
* **Local curation may not override a published fact.** Only the two fields with no public table are allowed there. A local override of something ClinGen publishes is a way to be quietly wrong.

### What this validates, and what it does not

The whole mechanism reference set runs a second time with the gene context read from files rather than typed into each case, and the mechanisms come out identical. That proves the **wiring** — it rules out a field ingested under the wrong name, a score mapped to the wrong enum, a MANE row picked wrongly, any of which would move a mechanism here and nowhere else.

It does **not** prove the biology. The fixture files were generated from the values that were already inline, so every number in them was originally typed by hand. Pointing the loader at the real ClinGen gene curation list and the real gnomAD constraint file is a flag away, and it is what would make "reproduced" mean something stronger than "self-consistent". Until then, the pins are earned but the values are still ours.

## Tissue: where the gene is switched on

A mechanism is a claim about a gene product. Where that product is not made, the claim has nothing to attach to — and until `context/expression.py` existed, the package had no way to say so. GTEx median TPM joins ClinGen and gnomAD as a third published table, pinned by content digest like the rest, and the affected tissue is supplied per case (`--tissue "Brain - Cortex"`) because it is a fact about the patient, not about a release.

It buys two things and nothing more.

**A gene silent where the disease is becomes a caution, not a refutation.** The mechanism survives — the arithmetic that produced it has not changed — and the confidence drops with the reason attached:

```
SCN1A  NM_001165963.4
  mechanism   loss_of_function  (possible)
  caution     GENE_SILENT_IN_THE_AFFECTED_TISSUE
              The gene is measured as essentially silent in the tissue this disease
              affects. A variant in a gene that is not transcribed there is unlikely to
              be what is causing the phenotype ... This is a caution rather than a
              refutation: bulk tissue hides a cell type that is two percent of it, and
              adult measurements say nothing about the developmental window that mattered.
```

That hedging is the honest reading of what GTEx is. It is bulk, adult, post-mortem tissue: a gene switched off at fifty may have been essential at four weeks, and averaging over a whole cortex hides the interneuron the disease is actually about. Treating a low TPM as proof of nothing-to-see would be a confident wrong answer, so the rule caps confidence and says why.

**Silence rules out the modalities that work through the native locus, and only those.** Upregulation needs a locus that is being transcribed; reactivation lifts a brake on a locus the cell would otherwise read; allele-specific silencing needs a transcript to remove. Where the gene is off, all three have nothing to act on. **Gene addition is deliberately absent from that block** — a delivered transgene carries its own promoter, so the native gene being silent says nothing about whether a supplied copy would be expressed. Contraindicating it here would be a wrong answer with a plausible reason attached, and a test exists to keep it out.

Three refusals:

* **Measured zero and never measured are different answers.** A gene GTEx assayed in liver and found at 0.0 TPM is evidence; a tissue nobody assayed is `None`, and no rule fires on it. Collapsing the two would let an absent column silence a gene.
* **The 1 TPM floor is a convention and lives in the rule file.** It is where the field draws the line, not something anyone measured. A laboratory that wants 0.5 edits `expressed_above_tpm`, not the code — and the pin changes with it, so the reanalysis layer attributes any reclassification to `rules` rather than to a discovery.
* **The fine vocabulary and the coarse one are kept apart.** Rules about the gene ask "is it on in *Brain - Nucleus accumbens (basal ganglia)*"; rules about delivery ask "could anything reach the central nervous system". Conflating them would let a delivery rule key off a basal ganglia subregion.

Expression also earns a drift axis, which is the expensive choice: `World` refuses to build without every axis pinned, so adding it broke every existing world until each was re-pinned. That cost is the point — a new GTEx release can move an answer without anything about the patient changing, and an input nobody has to acknowledge is an input the system reads without naming.

The limit is stated in every selection that gets this far, because it is the limit of the whole package:

```
read with caution
  · delivery to central_nervous_system is not assessed anywhere in this package.
    Whether a vector, an oligonucleotide or an editor reaches Brain - Cortex at a
    useful dose is the question that decides most of these in practice, and nothing
    here answers it
```

Knowing a gene is on in the cortex is not knowing that anything can be got into the cortex. A case run with no tissue at all gets the opposite caveat — that the check was never made — rather than a silent pass.

## Reading the patient's file

Everything in `annotation/` and `context/` is reference material. `vcf.py` reads the one input that describes the patient — and it is the only file the project consumes that nobody else versions.

**Genotype is where zygosity comes from**, which is why this module and the zygosity rules were worth having in the same package. The case worth knowing about is `1/2`: two *different* alternate alleles at one site means the sample carries no reference allele at all. Every intervention that works by raising output from an intact copy is off the table — the same conclusion as a homozygote, reached by a different route, and the record keeps the distinction even though the consequence is identical.

There is a case the reader cannot see, and it is written down rather than hidden: many callers emit `1/1` for a male X chromosome instead of a haploid `1`. That reads here as homozygous rather than hemizygous. The two agree on the only question the rules ask of zygosity, so the misreading is inert — but distinguishing them properly needs the sample's sex, which a VCF does not carry.

Four more refusals:

* **Decomposition without normalisation is wrong.** Splitting a multi-allelic record gives representations that are not canonical, and a variant written non-canonically fails to match the same variant written properly — including the patient's own earlier report. Handed a reference, the reader left-aligns on the way out; without one it sets `normalised = False` rather than pretending.
* **Consequence is read, never predicted.** If the VCF carries a CSQ or ANN field the reader takes the most severe term it recognises; if it does not, the record is simply not interpretable. Guessing "missense because the alleles are the same length" would be a prediction dressed as a parse.
* **A multi-sample VCF must name its sample.** Taking the first column would silently interpret a parent's genotype as the child's, and every answer downstream would be about the wrong person.
* **The assembly is checked against the header.** A VCF called against a different build gives coordinates that are all plausible and all about the wrong part of the genome.

Registering a case from a VCF reports what it dropped — uninterpretable records, variants that could not be placed on a coding transcript — because a run that silently watches four of a patient's forty variants is worse than one that refuses. Nobody would know to ask.

## Reanalysis: the answer can change without anybody learning anything

The genome does not change. What changes is everything the rules consulted while reading it — and `reanalysis/` is what makes that difference usable.

A mechanism call is not a fact about a variant. It is a fact about a variant **and a world**: seven pinned coordinates covering ClinVar, population frequencies, gene curation, the panel, the patient's phenotype, the transcript annotation, and our own rule files. Name the world and two things follow. A call becomes reproducible — re-pin, re-run, expect the same answer. And the difference between two calls becomes **attributable**, because the worlds differ in a small, enumerable number of coordinates.

**Attribution is established by experiment, not by correlation.** When several axes move in one week, the module re-runs the rules on counterfactual worlds: the old world with exactly one axis advanced, the new world with exactly one held back. Each moved axis comes back *decisive*, *sufficient*, *necessary* or *contributing* — the last meaning it moved and changed nothing. The cost is at most two re-evaluations per moved axis, and when a single axis moved it is free, because both endpoints are already in hand.

**Our own corrections are not discoveries.** `rules` is the one non-clinical axis. A change attributable only to a rule-file edit goes to a validation queue, not to clinical sign-out — it may well *become* a finding once someone confirms the new rule is right, but the confirmation comes first. Crucially this is not "ignore anything where the rule file changed": a rule edit landing in the same week as a real recuration comes back *contributing*, and the counterfactuals are what tell the two apart.

**`annotation` is an axis, and the test that justifies it is the sharpest case in the project.** A new transcript release moves the exon boundaries. The same stop codon that sat comfortably upstream of the last junction now sits within fifty nucleotides of it, so the transcript escapes nonsense-mediated decay, a truncated protein is made, and a collagen-like gene goes from haploinsufficiency to dominant-negative. Nobody learned anything about the patient. Gene addition goes from indicated to contraindicated anyway — and the run reports it as `mechanism_inverted`, `critical`, `caused by annotation`.

What is watched is deliberately the *mechanism* rather than an ACMG tier. A verdict moving from likely-pathogenic to pathogenic is a tier shift that changes nothing anybody does; a mechanism moving from loss-of-function to dominant-negative inverts the therapy. The vocabulary reflects that: `mechanism_inverted` and `modality_withdrawn` are the two critical kinds, and a withdrawn route outranks a newly opened one because something may already have been planned around it.

## Running it as a scheduled job

A reanalysis run is a process that starts, compares today with last month, and exits. Cron owns the schedule; the process owns one comparison.

```bash
repairbench watch NICU-014 --state /var/lib/repairbench --catalogue catalogue.yaml \
    --vcf child.vcf --sample CHILD --fasta GRCh38.fa --phenotype hpo-day-1
# 5 carried alleles read from child.vcf
# not placed on a transcript: 17-999-G-GA (intronic, and this module reasons in CDS coordinates)

repairbench reanalyse NICU-014 --state /var/lib/repairbench --catalogue catalogue.yaml
# NICU-014: 1 change(s) across 1 variants, 1 needing prompt review (axes moved: gene_curation)
#   → PLUSG-c158  loss_of_function → undetermined  [mechanism_lost, high]
#     the mechanism no longer resolves; anything downstream of it is unsupported
#     until it does — caused by gene_curation

repairbench serve --addr :9090      # /health and /metrics between runs
```

**The catalogue is what lets a counterfactual be real.** It maps `(axis, version)` to a file, so re-assessing a variant "as of January" loads January's file. Asked for a release the deployment no longer holds, it *refuses* — because falling back to the current one would produce a confident causal claim about an experiment that was never run. That refusal is the difference between attribution and a plausible story.

**A case must be registered before it can be run.** A scheduled run will not invent the list of variants it is meant to be watching.

**State is plain JSON on disk, and assessments are never overwritten.** An operator debugging a run at three in the morning should not need the library to find out what it last concluded. A stored event keeps its fingerprint — so an acknowledged transition stays acknowledged across processes — but *refuses to hand back its attribution*: the causal claim belonged to the run that established it, and reconstituting one a year later would misrepresent its footing.

**The one alert worth setting** is not an error rate. A scheduled job that stops running is invisible unless something measures its absence:

```promql
time() - repairbench_last_run_timestamp_seconds > 172800
```

### The sharpest test in the suite

Two changes produce an identical outcome — a settled mechanism becomes unsettled — and the system must not report them alike:

| | outcome | routed to |
| --- | --- | --- |
| ClinGen refutes dosage sensitivity | `mechanism_lost` | **clinical sign-out**, high — "caused by gene_curation" |
| We move the NMD boundary from 50 nt to 55 | `mechanism_lost` | **validation**, routine — "not because the evidence did" |

Both are defensible readings; only one is the field learning something. The test runs against real files in `tests/data/deployment/`, with two rule revisions and two curation releases on disk, because that claim is not testable against a stubbed engine.

## Scope and its edges

The feasibility flags are **necessary conditions, never sufficient ones** — `fits_viral_payload` means the coding sequence is under 4.4 kb, not that a vector exists, reaches the tissue, or is safe. The field names are chosen so a reader is not tempted to hear the stronger claim.

Designing candidates is M7 and is not here. The roadmap builds in this order for the reason both modules exist: each depends on the one before being right.

### Zygosity, and the answer that was wrong

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

## Annotation: turning asserted inputs into earned ones

Until this layer existed, a case file asserted its own exon structure and its own CDS offset. The NMD calculation — the one that decides between two opposite therapies — ran on numbers nobody had checked. `src/repairbench/annotation/` closes that: GFF3 in, a real transcript out, with the coordinate verified against the reference genome on the way.

```
$ repairbench explain case.yaml --annotation refseq.gff3 --fasta GRCh38.fa

# annotation  mini.gff3@103d7eefd327
# normalised  17:3305 T>TT moved 6 bases left to 17-3299-G-GT
# transcript  NM_000006.1 (MANE Select), 2 coding exons, 900 nt; variant at c.300

DEMOG  NM_000006.1
  mechanism   loss_of_function  (probable)
  ...
```

Four things it does that a naive version gets wrong, each of which corrupts an answer rather than raising one:

**Strand.** On the minus strand the first coding exon has the *highest* genomic coordinate. Sorting CDS blocks by start position and calling that transcript order inverts the exon numbering for half the genome — and puts the last junction at the wrong end of the gene, which inverts the NMD prediction with it.

**Coding, not exonic.** UTRs are absent from CDS records, which is what makes these offsets CDS offsets. Mixing exon and CDS coordinates shifts every position by the length of the 5′ UTR and lands the variant in the wrong exon.

**Left-alignment, finally.** M5 refused to left-align and said why: trimming shared affixes is safe with no external input, shifting is not. With a reference in hand the refusal is lifted, and the algorithm is the one bcftools and vt implement — the version where an allele may be *empty* internally, because that is what makes the loop's two moves alternate cleanly. Every one of the ten ways to write the same insertion into a run of ten A's now normalises to one key. Before, they were ten different variants that would each have failed to join against the patient's own earlier report.

**Reference verification.** The cheapest detection of the most damaging error there is: a VCF called against a different assembly than the annotation. Every coordinate stays plausible, every lookup succeeds, and every answer is about the wrong part of the genome. One comparison catches it.

The FASTA is read through its `.fai` index rather than loaded — left-alignment needs single bases upstream of a variant and a human genome is three gigabytes. A missing index is reported rather than silently built, because building one means reading the whole file, which is what the index exists to avoid.

Every source is pinned by content digest, so a call can name the annotation release it was made against — the same discipline the rule files already keep.

### What this is and is not validated against

The fixtures in `tests/data/` are **synthetic and labelled as such**: a 5 kb reference and a hand-written GFF3 with a plus-strand gene, a minus-strand gene, a gene with no MANE Select transcript, a non-coding transcript that must be dropped, and two homopolymer runs — one intronic, one inside a coding exon. What they test is the parser, the strand handling, the coordinate arithmetic and the shift algorithm, none of which care whether the coordinates correspond to a real locus.

They are not a test against real RefSeq annotation, and the README will not pretend otherwise. Pointing this at a real GRCh38 GFF3 is a `--annotation` flag away and is the obvious next validation: parse the genes in the reference set, and check that the mechanism calls survive the move from uniform fixture exons to real ones. Where they do not, the reference set was passing for the wrong reason — which is exactly what this layer exists to find out.

## The seam: one case, end to end

Every layer above was a separate command, and that was a hole rather than a design. You could run `assess`, learn that gene addition is contraindicated, then run `design` on the same variant and get a page of protospacers — because the two commands did not know about each other. M6's whole safety property stopped at M6's edge.

`repairbench plan` closes it:

```
$ repairbench plan case.yaml --fasta target.fa --annotation target.gff3

# annotation  target.gff3@44a62720b05a
# transcript  NM_000099.1 (MANE Select), 1 coding exons, 101 nt; variant at c.52

TARG  NM_000099.1
  mechanism   dominant_negative (probable)
    NULL_ALLELES_ARE_MILDER [strong]
    MISSENSE_IN_MULTIMER [supporting]
  not designed, ruled out by the modality rules
    gene_addition
    truncated_construct
    wild_type_upregulation
    ...
  designed
    allele_specific_silencing → aso
      TARG  17:241-361  gapmer-2MOE (cleaves, 20 nt)
        tiled       102 windows; 20 with nothing blocking, 0 blocked
        note        only windows covering 17:301 were kept. That base is the one thing
                    telling the two transcripts apart ...
    base_editing → base_editor
      ...
  made under
    mechanism  mechanism-v1@0adb1c56d9d9
    modality   modality-v1@1e2a1e66c493
    editors    editors-v1@f6a7a69f266d
    prime      prime-v1@f44edb7f0284
    aso        aso-v1@5e926df550ae
    routing    routing-v1@fd18a39e4191
```

**A modality M6 ruled out is never designed.** Not designed and marked; not designed with a warning at the bottom — not designed. A sequence is something a reader can order and a caveat is something a reader can skim, and that asymmetry is the entire argument for making this a refusal. The check runs before any other, so no combination of missing inputs or later failures can route around it. A test asserts it across every modality rather than one, because the failure it guards against is a route added later that skips the check.

Two weaker rules follow. An unresolved mechanism designs nothing at all. And a modality that was merely *not indicated* is also not designed, but is listed with its verdict — the reader should see that it was considered and why nothing came of it.

**The routing table is data too.** `rules/routing-v1.yaml` says which modality goes to which designer, and the entries with *no* designer are the ones that earn their place: gene addition returns no molecule, and a package that quietly returns nothing there looks like a package that found nothing. What it says instead is that the sequence to deliver is already known and what has to be designed is a vector, a promoter and a dose — none of which is sequence design. The table is refused at load time if it does not mention every modality, because an unlisted one produces silence that reads as "none found".

Two corrections the seam forced, both of which the separate commands could never have caught:

* **An allele-discriminating oligonucleotide is tiled against the patient's sequence, not the reference.** One complementary to the reference base is complementary to the *healthy* transcript — the allele it was supposed to spare.
* **And every window must cover the variant.** That base is the only thing telling the two transcripts apart, so a gapmer that misses it knocks down both alleles, which for a dominant-negative variant removes the good product along with the bad. Whether covering it discriminates *enough* — a single mismatch in a 20-mer often does not — is stated as unanswered.

## M7 — designing the edit, and the two things nobody can design

M6 ends at a class of intervention. M7 takes one of those classes — base editing — and produces the actual placement: which editor, which protospacer, on which strand, and every unintended base that comes with it.

```
$ repairbench design --gene TARG --at 17:301 --patient A --wild-type G \
    --fasta target.fa --annotation target.gff3

TARG  17:301  A→G
  candidates  6 (1 with no other editable base in the window)
    ABE7.10-SpCas9  GCTGATACTGCTGCTGCTGC agg  17:297-316 (+)
        target at position 5 (g.301), A>G
        bystanders (1):
          · position 7 (g.303) → G, in coding sequence
        ! a bystander above is in coding sequence, where a silent change and a missense
          change look identical from here
  considered  ABE8e-SpCas9 on the + strand, ABE7.10-SpCas9 on the + strand, ...
  ranking     no efficiency model is attached, so the candidates below are not ranked by
              how well they would work — only grouped by whether the editor can reach
              another base in the same window
  catalogue   editors-v1@f6a7a69f266d
```

**The editors are data.** `rules/editors-v1.yaml` carries the conversion, the PAM and the editing window for each one, with the paper each number came from. All three are measurements that get revised — ABE8e's window is what a group measured, not a property of arithmetic — so widening a window is an edit to a pinned file, visible in the digest of every design made after it, rather than a constant somebody changed.

Four things it does that a first version gets wrong, each of which produces plausible output rather than an error:

**It designs against the patient's sequence, not the reference.** The base being corrected is, by definition, the one the reference does not have. A scan over reference bases looks for a PAM around a base that is not there — and where the variant creates or destroys a PAM, it invents guides this patient does not have, or misses the ones they do.

**It decides the strand from the conversion.** A deaminase makes A→G and C→T, full stop. Correcting a patient's T to a C is not a cytosine editor's job; it is an adenine editor's job on the *other* strand, where the same base pair reads A and needs to read G. Roughly half of all correctable variants are only reachable that way, and a single-strand scanner reports them as "no candidates" — a wrong answer wearing the clothes of a missing one.

**It numbers positions from the PAM-distal end**, which is how every window in the literature is quoted. The other direction puts every window at the wrong end of every guide and still lands the target inside one often enough to look right.

**It lists bystanders individually, with coordinates.** Every other editable base in the window is an unintended change to a person's genome. They are never summarised as a count, and where an annotation is supplied each one is marked as coding or not — because a silent change and a missense change look identical from here, and the report says so rather than implying it worked out.

Refusals are first-class and carry their reason: a transversion is refused with the note that base editing does not make one and prime editing is the route; an indel is refused before any sequence is read; no PAM at a usable distance is a different message from either.

### Efficiency: the model that is deliberately absent

How often an editor actually converts its target is what BE-Hive and its successors predict, from screens of measured outcomes. None is attached here — no weights, no way to run one — and there were three options:

Invent a heuristic that looks like the real thing. It would produce a ranked list, a ranked list reads as knowledge, and somebody would order the top guide. Return candidates in found-order and say nothing; the order still reads as a ranking, because lists do. Or make the absence an object every report prints.

`NoModelAttached` is the third, and it is the correct implementation rather than a stub: it refuses in the same shape a real model would answer, so attaching PRIDICT or BE-Hive later is a constructor argument, not a rewrite. What ordering exists is declared — fewest bystanders first, then position — which is the difference between a criterion and a heuristic.

### Off-target: the part the search tools leave out

Finding every site in a genome within six mismatches is Cas-OFFinder's job, on a GPU, over three gigabytes. This package reads its output rather than reimplementing it badly. What it adds is the ranking, and the claim is that a hit list sorted by mismatch count is sorted by the wrong thing:

```
$ repairbench offtarget hits.txt --annotation refseq.gff3 \
    --gene-lists lists.tsv --expression gtex.tsv --tissue "Muscle - Skeletal"

  prohibitive  (3)
    17:1550 (+)  4 mismatches  PLUSG (NM_000001.1), coding sequence
        CODING_HIT_IN_AN_ESSENTIAL_GENE
    17:545 (+)   0 mismatches  PLUSG (NM_000001.1), coding sequence
        A_HIT_WITH_NO_MISMATCHES
  moderate  (1)
    17:700 (+)   1 mismatch   PLUSG (NM_000001.1), intron or untranslated region
        HIT_IN_A_GENE_SILENT_IN_THE_TARGET_TISSUE
  low  (1)
    11:5000 (+)  5 mismatches  not in any transcript in this annotation
```

The four-mismatch hit outranks the two-mismatch one, and that inversion is the whole point. Mismatch count says how likely the nuclease is to bind; it says nothing about what happens if it does, and the second question is the one being asked. The rules live in `rules/offtarget-v1.yaml` in the same predicate language as the mechanism and modality rules, because which off-target hits are unacceptable is a clinical judgement like any other.

This is also where the tissue dimension pays for itself twice: a non-coding hit in a gene that is silent in the target tissue is downgraded — not dismissed, because an edit is permanent and a cell's expression programme is not.

Three refusals worth naming:

* **No CFD score is invented.** Weighting a mismatch by its position and identity needs a published table this package does not carry. The mismatch count is reported as the search gave it, and the report says that is what you are reading.
* **Nothing is ever marked safe.** The lowest tier means no rule fired — a statement about the rule file, not about the genome. And a hit list run without an annotation comes back `unassessed` rather than `low`, so "we did not look" cannot be read as "we looked and it was fine".
* **A missing gene list is not clearance.** Absence from DepMap or COSMIC is absence of evidence about a gene. A coding hit still ranks on the generic rule, and the rule file says why: DepMap measures essentiality in cell lines, and a gene dispensable in culture can be indispensable in a neuron.

The command exits non-zero on a prohibitive hit, so a pipeline stops rather than logging it somewhere nobody reads.

### Prime editing: writing the edit into the guide

Base editing covers four transitions. Prime editing covers everything — all twelve substitutions, and small insertions and deletions besides — and pays for it with a design space large enough that picking badly out of it is the normal outcome.

```
$ repairbench pegrna --gene TARG --at 17:500 --patient A --wild-type C --fasta target.fa

TARG  17:500  A→C (substitution)
  pegRNAs     10 across 1 protospacer(s); 0 blocked by a rule
    TGCTGCTGCTGCTGCTGCCC agg  17:470-489 (+)
        nick at g.486, 14 nt from the edit; the PAM survives the edit
        PBS  9 nt  CAGCAGCAG
        RTT 19 nt  GCAGCGGCAGCAGCCTGGG  (5 nt of homology past the edit)
        3' extension  GCAGCGGCAGCAGCCTGGGCAGCAGCAG
          · PE3b  AGCAGCAGCGGCAGCAGCCT ggg (-)  nick at g.493, 7 nt away
        ! [caution] THE_PAM_SURVIVES_THE_EDIT: ... the editor re-nicks what it has
          just corrected, and the second pass can install an indel instead
        ! [note] A_PE3B_NICK_IS_AVAILABLE: ... it cannot fire on the unedited allele,
          so the two nicks are never open at once
```

**One inequality carries the module: the nick must fall before the edit.** Reverse transcription runs from the nick forwards, so a protospacer whose nick lands past the edit writes over sequence the polymerase never reaches. It looks perfect — good PAM, right distance, clean spacer — and it does nothing. That constraint eliminates roughly half the PAMs near any variant, and the refusal names it rather than reporting "no candidates", because the fix is a PAM-relaxed nickase and not a longer template.

**Everything else is four reverse complements, and each is a silent failure.** The primer binding site is the reverse complement of the bases *upstream* of the nick; the template is the reverse complement of the genome *downstream* of it, with the edit written in; a minus-strand protospacer complements both again. Get one wrong and the pegRNA has exactly the right length and installs nothing. So the test suite verifies the template from the outside: reverse-complement it back, and it must equal the genome as it should read, base for base.

**PE3b is found and named.** A second nick raises efficiency several-fold and is also where PE3's indels come from, because for a while both strands are cut at once. A PE3b guide is one whose own spacer matches the *edited* allele — so it cannot fire until the edit is installed, and the two nicks are never open together. It has to overlap the edit, which puts its nick far closer than PE3's 40–90 nt window, so the rule file gives it its own search range and says why.

The nuclease, the window, the PBS melting target, the nick-to-edit ceiling all live in `rules/prime-v1.yaml`, because almost none of them is a fact about chemistry — they are fitted recommendations from screens, and they will be revised.

### Antisense oligonucleotides: the easiest thing here to design badly

```
$ repairbench aso --gene TARG --at 17:240-360 --fasta target.fa \
    --chemistry gapmer-2MOE --exon 250-350

TARG  17:240-360  gapmer-2MOE (cleaves, 20 nt)
  tiled       102 windows; 76 with nothing blocking, 26 blocked
    AGCAGCAGCAGCAGCAGCAG  17:253-272  gapmer-2MOE, exon interior
        20 nt  GC 65%  Tm ~66 °C  CpG 0
  blocked     26, worst first
    AGCAGCAGCAGCAGCAGCAG  17:241-260  gapmer-2MOE, acceptor site
        ! [blocking] A_SPLICE_SITE_TARGET_NEEDS_A_STERIC_BLOCKER: A gapmer recruits
          RNase H, which cuts the transcript it binds. Aimed at a splice site to
          redirect splicing, it destroys the very transcript the redirection was
          meant to rescue.
```

Tiling is trivial: slide a window, take the reverse complement, count the GC. What makes it worth writing down is the one confusion that inverts the therapy rather than degrading it. **A gapmer destroys what it binds; a steric blocker occupies it.** Aim a gapmer at a splice site to skip an exon and it degrades the transcript the skip was meant to rescue — every composition rule passes, and the therapy is backwards. So the chemistry declares an action in the rule file, and a cleaving chemistry over a splice boundary is blocked outright.

The mirror-image flag matters too: a steric blocker in the middle of an exon binds and changes nothing, unless it happens to cover a silencer or an enhancer — and this package reads CDS records only, so that one is a caution rather than a block, with the reason stated.

The rest is composition, from `rules/aso-v1.yaml`: four consecutive guanines block a window (G-quadruplex), CpG dinucleotides are counted and flagged (TLR9), GC bounds and melting temperature are cautions, self-complementarity is measured crudely and labelled as crude.

**And the honest part:** what decides an antisense oligonucleotide is whether the site is *accessible* — whether it is paired up inside the transcript's own fold. Two windows identical on every rule above can differ tenfold for that reason alone. Folding the target is RNAfold's job, no folding model is attached, and so a tiling run here produces starting points rather than candidates. The report says exactly that, next to the number.

## What is missing

* **Real transcript structures.** The reference set uses simplified exon lengths — real counts and coding lengths, uniform exons — because what it tests is the rule layer. Real annotation belongs with the ingest module.
* **More cautions than confidences.** The rule file currently has four rules that argue for uncertainty and nine that argue for a mechanism. That ratio should probably grow, not shrink.
* **Splice rules are thin.** They classify the mechanism but do nothing with the consequence; deciding whether a skipped exon rescues the frame or moves the problem downstream needs the exon-level detail that ingest will provide.
* **Tissue is bulk, adult and post-mortem.** GTEx answers "is this gene on in this tissue in an adult who died of something else". It does not answer which cell type within it, or what was happening in the developmental window that produced the phenotype — which, for most of the diseases this package is aimed at, is the window that mattered. Single-cell and developmental atlases exist; nothing here reads them.
* **The seam is one-way.** A plan runs mechanism → modality → molecule and stops. What it does not do is feed the design back: a modality with no workable candidate — every pegRNA blocked, every window refused — should arguably weaken the modality's verdict, and today it does not. The plan reports both facts and leaves the reader to connect them.
* **Delivery is not assessed at all.** Whether a vector, an oligonucleotide or an editor reaches the affected tissue at a useful dose is the question that decides most of these in practice. Every selection says so; none of them answers it.
* **Three designers, three absent models.** Base editing has no BE-Hive, prime editing no PRIDICT, antisense no folding model. All three are one class against a Protocol that already exists; obtaining and validating the weights is the actual work, and until it is done every list here is grouped rather than ranked.
* **No silent PAM edits.** When the PAM survives a prime edit the report says so and stops there. Designing the silent change that destroys it — checking the reading frame, picking a synonymous codon — is the obvious next thing the module should do for itself.
* **Splice-regulatory elements are unread.** The ASO module places a window against CDS boundaries and nothing else. Silencers and enhancers — ISS-N1, the element nusinersen was built around — are exactly what a steric blocker aims at, and no annotation here contains them.
* **The design fixture is synthetic.** `tests/data/design/target.fa` is a repeating trimer with PAMs placed by hand at known distances, which tests the arithmetic and nothing else. A real locus with a real published guide — one somebody actually made and measured — is the validation that would mean something.
