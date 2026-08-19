# repairbench

**Two nonsense variants in the same gene can have opposite mechanisms and need opposite therapies. This works out which one you have, which interventions that admits — and then keeps watching, because the answer can change without anybody learning anything.**

Every clinical judgement lives in a rule file — `mechanism-v1`, `modality-v1`, `editors-v1`, `prime-v1`, `aso-v1`, `offtarget-v1`, `routing-v1`, all under `rules/` — and not in the code. If a conclusion here is wrong, the error is in one of those files, in a sentence a geneticist can read. That is the point.

The layers, one package:

| | |
| --- | --- |
| `vcf.py` | the patient's own file in; carried alleles, zygosity and read consequences out |
| `context/` | ClinGen dosage, gnomAD constraint, GTEx expression and ClinVar submissions in; gene-level facts out, with a citation on every one |
| `annotation/` | GFF3 and reference FASTA in; real transcripts, verified coordinates, left-aligned indels out |
| mechanism (M5) | why this variant causes disease — loss of function, gain of function, dominant negative — with its evidence |
| modality (M6) | which classes of intervention that mechanism admits, and which it rules out |
| design (M7) | the molecule: a base editor's protospacer, a pegRNA, or a tiling of antisense oligonucleotides — each with what is wrong with it |
| `plan.py` | the seam: one case from mechanism to molecules, where a ruled-out modality is never designed |
| `reanalysis/` | watching all of the above drift as releases land, and proving which release moved it |

![Python](https://img.shields.io/badge/python-3.11%2B-3776AB)
![Rules](https://img.shields.io/badge/rules-as%20data-6f42c1)
![Tests](https://img.shields.io/badge/tests-495-brightgreen)
![Reference set](https://img.shields.io/badge/reference%20cases-26%2F26-brightgreen)

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
make check       # lint, types, 424 tests, all three reference sets
make reference   # re-run the reference sets alone, case by case
make rules       # print what the system believes, with citations
repairbench assess case.yaml   # mechanism, then the modalities it admits
```

[`RUNNING.md`](RUNNING.md) is the longer version: setup, one command per layer,
and what each one's output is claiming.

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

### ClinVar: the number the most inferential rule was reading

The four facts above are the ones a gene *has*. There is a fifth the rules read, and for most of this project's life it was the weakest thing in the package: **where pathogenic variation sits in the gene**.

It matters because of one rule. Nothing in ClinGen or gnomAD distinguishes gain of function from loss of function — both describe how badly the gene is needed, not what a variant does to it. What separates them is a pattern: pathogenic missense variants piling into one stretch of protein while truncating variants are *absent*, which is what a gene that causes disease by doing something new looks like in a variant database. That inference is the least direct thing the rule file attempts, and it was resting on `distribution:` blocks typed into a YAML file from memory of the literature. Plausible numbers, invented — the exact failure this package exists to be about.

`context/clinvar.py` counts them. Three refusals shape it, and each one has a way of producing a confident wrong number:

* **"Conflicting classifications of pathogenicity" contains the word "Pathogenic".** A parser matching the substring counts submitter disagreement as support. The accepted classifications are written out rather than matched.
* **Review status is not decoration.** One laboratory with no assertion criteria and a reviewed expert panel are not the same evidence, and a count that averages them is how a tally becomes confident nonsense. The star rating is a threshold, kept per variant, and every citation says what was counted at what level: `6 pathogenic submissions at ≥1★ (1×3★, 2×2★, 3×1★)`.
* **What kind of variant this is comes from the protein, not from ClinVar's `Type` column.** `Type` says how the *sequence* changed — deletion, duplication, single nucleotide variant. The rules ask what happened to the *product*, and one single nucleotide variant can be missense or nonsense, which argue for opposite mechanisms.

And one thing it is careful to *name* rather than claim. A hotspot here is the densest window of N residues, N being a rule-file threshold, computed by a two-pointer sweep over sorted positions rather than by binning — binning would make the answer depend on where the bin edges happened to fall, which for a gene whose cluster straddles one understates exactly the clustering that matters. That is a measurement of tightness. It is *not* a curated functional domain, this file has no access to one, and the output says "densest 20-residue window" everywhere it could be mistaken for one.

The count is now an ingested fact like the other four, with a pin and a citation:

```
$ repairbench context COL1A1 --clinvar refdata/variant_summary.txt.gz

COL1A1
  distribution = 3/3 missense in the densest window (100%), 2 truncating  [clinvar@2026-08/5d866ed862b2]
      6 pathogenic submissions at ≥1★ (1×3★, 2×2★, 3×1★); hotspot = densest 20-residue window
```

Loading it without naming genes is refused. The file is millions of submissions; reading all of them to build context for nine would take minutes and look like it was working.

The reference set still holds its counts inline rather than reading ClinVar at test time — CI has to run without a 250 MB download — so `repairbench clinvar` prints them in the shape the YAML wants, which makes updating after a release transcription rather than judgement. It prints real c. and p. positions alongside, on the transcript the submitters used, because a c. position is only meaningful against the transcript it was written on.

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

### The queue, as a page somebody opens

Everything above renders one answer to one question asked now. Reanalysis is not that — it runs at three in the morning and exits — so `repairbench dashboard --state /var/lib/repairbench --out queue.html` writes the work list as a self-contained page: which cases moved, what the change was, which queue it went to, and the world each case was last compared against.

The design question was not what to show. It was **what an empty page means**.

A queue with nothing in it is the commonest correct output this system produces. It is also exactly what a scheduler that died in March looks like, and the two must not render alike — so the loudest element on the page is not an event, it is a case nobody has examined recently. Under a healthy run the page says so in words: *an empty queue below therefore means nothing moved, rather than nothing ran*.

Three refusals hold it up. **Nothing is recomputed** — every urgency and queue is reproduced exactly as the run recorded it, so the page cannot promote an event or quietly demote one. **Quiet cases still get a row**, because three rows have to mean "three cases exist" rather than "three cases need attention". And **no JavaScript, no network, no build step**: one file that opens from disk, because a dashboard that needs a server running is a dashboard that is down exactly when the pipeline is.

Building it turned up three defects underneath it, which is the usual result of making something visible:

* **Urgency was being lost on every save.** A ledger written back to disk recorded its events' urgency and queue as `"-"`, so a case loaded and saved twice forgot how loudly it had been asking. Nothing read those fields until a queue needed to be listed.
* **`case_ids()` returned things that were not cases.** The command line writes `<case>.variants.json` beside each ledger, and a `*.json` glob read that back as a case named `NICU-014.variants`. The first thing the dashboard ever did was crash on it. Cases are now identified by shape — a ledger is an object with a `case_id` — rather than by file name.
* **The ledger could not tell "quiet" from "never ran".** It recorded when a case last *changed*, never when it was last *examined*, so the first version of the page reported every healthy silent case as dead. `last_examined_at` is now written on every run, especially the ones that find nothing.

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

### The third reference set: episodes rather than mechanisms

The other two reference sets ask whether the rules reach the answer the literature reached. Reanalysis has no answers of that kind — it has *episodes*: a release lands, an assessment moves or does not, and what is under test is the causal claim and who hears about it.

```
$ repairbench reference --reanalysis

ok    A curation removes a settled answer
        mechanism_lost → gene_curation → clinical_signout (high)
ok    Our own rule edit reaches the same outcome
        mechanism_lost → rules → validation (routine)
ok    Two releases land, one of them matters
        mechanism_lost → gene_curation → clinical_signout (high)
ok    Two causes, either would have done it
        mechanism_lost → gene_curation, rules → clinical_signout (high)
ok    A release lands and nothing changes
        none → nothing → none (silent)

5/5 reanalysis episodes reproduced
```

Five shapes, and the last one is the one a laboratory meets most weeks: an annotation release re-issues the gene with identical structure, the digest moves, and nothing a rule reads differs. It must reach nobody. A reanalysis system that reports that week is a system somebody stops opening, and then it misses the week that mattered.

The set is data, in `tests/reference/reanalysis.yaml`, with the story of each episode written next to it — and it did its job on the first run. The rule-edit episode was written expecting *high* urgency, and the system returned *routine*. The system was right: nothing had been learned about the patient, so nothing was time-critical; what was needed was somebody confirming the new rule was the better reading. The expectation was wrong and is now corrected, with the correction recorded in the file.


### The browser as the operating surface

`repairbench review --state /var/lib/repairbench --catalogue catalogue.yaml` serves the queue with the daily loop on it: register a case, re-examine one or all of them, read what moved, sign it off. After that command the terminal is not in the loop.

`repairbench demo --state /tmp/try` is the same thing with a synthetic case already assessed against last month's releases, so the first button press has something real to report.

Building that turned up a defect the command line had been hiding. **Registration did not take a baseline** — it recorded which variants a case was watching and nothing about what we currently concluded, so the first run after registration was always a no-op: every variant was being seen for the first time, and that run's answer silently became the reference. Register a case in January, run it in April, and February's curation change goes unreported. From a terminal the first run is one of many and the gap is invisible; from a page it is the button somebody presses expecting an answer. Registration now assesses the case as of today, and a baseline that could not be taken is reported without discarding the registration — the variants are the part somebody typed.

The line that keeps this honest is **what the server is allowed to change**. It can start a run — which is *exactly* the comparison a scheduled process performs, from the same pinned files, so a button press cannot reach a conclusion cron would not have reached; it is the same run, started by a person instead of by a clock. It can record that a named person read a change. It cannot assert a mechanism, an urgency or a queue, because those come out of the rule files or they do not exist. A button that started a run is a person doing what the clock does; a button that changed a verdict would be a different product.

Two consequences of that line, both visible in the code. Registering and running moved out of `cli.py` into `reanalysis/operations.py`, called by both entry points — two copies of "what registering a case means" is a guarantee the terminal and the page drift apart, and the one that drifts is the one nobody tests. And the case identifier is checked for path separators before it becomes a file name, because a form facing a browser is a different threat model from an argument typed by the person who owns the machine.

### The one thing a reviewer needs to *do*

The dashboard writes a page, and a page cannot do anything. `repairbench review --state /var/lib/repairbench` serves the same queue with one button on it, and the button closes a hole that had been in this package from the beginning.

**Acknowledging was unreachable.** The surfacing policy suppresses any change whose fingerprint has already been acknowledged — that is what stops a queue from being a list that only grows. `CaseLedger.acknowledge()` existed, the policy read it, and *nothing exposed it*. A reviewer who read a change, decided it needed nothing, and closed the tab would be shown the same change at the next release, and the release after that.

What closing it forced was a question about evidence rather than about UI. Acknowledging is the only write in this package that makes the system **quieter**, and an anonymous switch that suppresses future alerts is precisely what an incident review cannot reconstruct. So `acknowledge` now requires a name, refuses a blank one, and records the note alongside it — because *why* a change needed nothing is the part a later reviewer wants and the part a boolean throws away. The ledger enforces that, not the form, so `repairbench acknowledge CASE EVENT --by … --note …` is held to the same standard as the button.

And the server says what it is worth. It binds to loopback, and the line under the button reads: *the name is recorded as attribution — this server does not authenticate it.* That is true, it is a real limitation, and writing it on the page is better than a password box that would imply otherwise. A deployment where the distinction matters needs an identity provider in front of this.

It is `http.server` from the standard library, the same thing the metrics endpoint already runs — no framework, no template engine, no database, no JavaScript. Every page is a full document and every action is a form post followed by a redirect, so a refresh cannot replay it. The server can change exactly one thing: whether an event is marked as read, and by whom. It cannot re-run an analysis, edit a rule, or alter an urgency, so the worst a misclick can do is mark one event read — and the ledger keeps who did it.

## Reproducing a drug, and the defect it found

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

### And the base editor, where the disagreement was the interesting part

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

### And prime editing, where it found the worst defect of the three

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

`./scripts/fetch-reference-data.sh` downloads the real files — RefSeq annotation, six chromosomes of GRCh38, ClinGen dosage, gnomAD constraint, GTEx expression — and `scripts/README.md` says what each one settles.

## What happened when it was pointed at real data

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

### And then the context files

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

### And then the last invented input, which cost two reference cases

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

### The name every chromosome has twice

Pointing this at real files produced a collision within a minute of the download finishing. NCBI's annotation calls chromosome 17 `NC_000017.11`; UCSC's FASTA calls it `chr17`; a clinical VCF calls it `17` or `chr17` depending on who wrote the pipeline. All four are the same sequence, and string equality says they are four.

The fix is not a table of accessions. That would go stale with every assembly patch, and would silently invent an answer for anything it had not heard of. NCBI opens each chromosome with a `region` record carrying `chromosome=17` in its attributes — so **the mapping is in the annotation, written by the people who assigned both names**, and `parse_gff3` reads it while it is streaming past anyway.

Two refusals hold it in place. An unplaced scaffold carries a `chromosome` attribute too, so only records marked `genome=chromosome` contribute an alias — without that check a patch contig would alias itself to the chromosome it patches, and coordinates would resolve against the wrong sequence. And where the annotation says nothing, exactly one substitution is attempted: adding or removing the `chr` prefix, which carries no information. Anything more ambitious would be the package guessing which contig somebody meant.

## What is missing

* **Recurrence is unread.** The clustering measure counts alleles, and what a hotspot means is that many patients carry the same one. Number of submitting laboratories is a weak proxy and is not read by any rule; the honest input would be per-variant observation counts, which the variant summary does not carry.
* **One window, and hotspots come in fours.** The clustering statistic takes the densest single window. A multi-domain protein has several, and the measure sees one of them — which is half of why the rule had to be demoted rather than repaired.
* **More cautions than confidences.** The rule file currently has four rules that argue for uncertainty and nine that argue for a mechanism. That ratio should probably grow, not shrink.
* **Splice rules are thin.** They classify the mechanism but do nothing with the consequence; deciding whether a skipped exon rescues the frame or moves the problem downstream needs the exon-level detail that ingest will provide.
* **Tissue is bulk, adult and post-mortem.** GTEx answers "is this gene on in this tissue in an adult who died of something else". It does not answer which cell type within it, or what was happening in the developmental window that produced the phenotype — which, for most of the diseases this package is aimed at, is the window that mattered. Single-cell and developmental atlases exist; nothing here reads them.
* **The reanalysis episodes are synthetic.** The shapes are real and the decisions they pin are the ones a laboratory needs, but the four files each episode runs against are ours. Reproducing a real week would mean holding a real ClinGen release, a real gnomAD release and a real annotation release from two different dates, and no laboratory publishes that set.
* **The seam is one-way.** A plan runs mechanism → modality → molecule and stops. What it does not do is feed the design back: a modality with no workable candidate — every pegRNA blocked, every window refused — should arguably weaken the modality's verdict, and today it does not. The plan reports both facts and leaves the reader to connect them.
* **Delivery is not assessed at all.** Whether a vector, an oligonucleotide or an editor reaches the affected tissue at a useful dose is the question that decides most of these in practice. Every selection says so; none of them answers it.
* **Three designers, three absent models.** Base editing has no BE-Hive, prime editing no PRIDICT, antisense no folding model. All three are one class against a Protocol that already exists; obtaining and validating the weights is the actual work, and until it is done every list here is grouped rather than ranked.
* **No silent PAM edits.** When the PAM survives a prime edit the report says so and stops there. Designing the silent change that destroys it — checking the reading frame, picking a synonymous codon — is the obvious next thing the module should do for itself.
* **Splice-regulatory elements are unread.** The ASO module places a window against CDS boundaries and nothing else. Silencers and enhancers — ISS-N1, the element nusinersen was built around — are exactly what a steric blocker aims at, and no annotation here contains them.
* **Three molecules is three, not a benchmark.** Each designer has been checked against one published object, and each check was worth it — two of the three found real defects. What has not been done is the systematic version: a set of published guides per modality, run as a suite, with the disagreements counted rather than investigated one at a time.
