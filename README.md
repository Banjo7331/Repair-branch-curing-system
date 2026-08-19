# repairbench

**Two nonsense variants in the same gene can have opposite mechanisms and need opposite therapies. This works out which one you have, which interventions that admits, which molecule to design — and then keeps watching, because the answer can change without anybody learning anything.**

Every clinical judgement lives in a rule file — `mechanism-v1`, `modality-v1`, `editors-v1`, `prime-v1`, `aso-v1`, `offtarget-v1`, `routing-v1`, all under `rules/` — and not in the code. If a conclusion here is wrong, the error is in one of those files, in a sentence a geneticist can read and change.

![Python](https://img.shields.io/badge/python-3.11%2B-3776AB)
![Rules](https://img.shields.io/badge/rules-as%20data-6f42c1)
![Tests](https://img.shields.io/badge/tests-495-brightgreen)
![Reference set](https://img.shields.io/badge/reference%20cases-26%2F26-brightgreen)
![Published molecules](https://img.shields.io/badge/published%20molecules-3%2F3-brightgreen)

> ⚠️ **Not a medical device.** Research and educational use. A candidate here is a *design* — a protospacer, a pegRNA, an oligonucleotide, and what is wrong with each — not a therapy. No model of editing efficiency, prime-editing yield or target accessibility is attached to any of them, and "not ruled out" is not a recommendation.

---

## Why mechanism first

The obvious first module to build is the one that designs sequences. It is the wrong one, because a sequence designed against the wrong mechanism is worse than no sequence at all.

Take a nonsense variant in *COL1A1*. If the premature stop sits far enough upstream of the last exon–exon junction, the transcript is destroyed, half the normal collagen is made, and the child has the mild form of osteogenesis imperfecta. The therapeutic question is how to supply more product.

Take a glycine substitution in the same gene. The protein is made, it is incorporated into the triple helix, and it ruins it. The phenotype is *worse* than having no protein at all — and supplying more normal collagen does not remove the defective chains. The question inverts: how to get rid of the mutant allele.

Same gene. Opposite intervention. A tool that designs correction sequences without settling this first will, half the time, confidently design the one that makes things worse.

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
        imperfecta type I, glycine substitutions give the severe forms.
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

Three properties follow, and each was the reason for a design decision.

**A geneticist can review it without reading Python.** The predicate language is *interpreted, not executed* — no `eval`, no expression strings, no way for a rule to reach anything the feature record does not name. A rule file is a document you can hand to a clinician without also handing them code they are running.

**A typo fails the run.** A rule naming a feature that does not exist raises, listing the features that do. The alternative — evaluating to false — produces a rule that silently never fires, which is indistinguishable from a rule that correctly never fires, and only one of those is a bug.

**Every call names the rules that produced it.** `mechanism-v1@7dd6bd311e2b` is the version plus the digest of the file's bytes. Two calls made under different rules are not comparable, and a system that cannot tell them apart will eventually report a rule change as a biological finding.

Thresholds are rule data too — the NMD junction boundary is quoted as 50–55 nt across the literature, and a laboratory preferring 55 says so in the file rather than in a patch.

## The pipeline

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

### Gene context, and which facts nobody publishes

Five facts feed the mechanism rules. Three are published as tables — ClinGen dosage sensitivity, gnomAD constraint, and where pathogenic variation sits in the gene, counted from ClinVar. Two are not: whether a gene product assembles into a complex, and whether its null alleles are milder than its missense ones, are judgements read out of the literature that nobody ships as a TSV.

Those two live in a local curation file that **demands a citation per entry**, and the difference shows in every line of output:

```
$ repairbench context COL1A1

COL1A1
  distribution = 25/446 missense in the densest window (6%), 657 truncating
      [clinvar@2026-08/5d866ed862b2]  1386 submissions at ≥1★ (353×2★, 1033×1★)
  forms_multimer = True  [local_curation@rev3/bbf8f9fd50de]
      Type I collagen is a heterotrimer; a mutant pro-alpha1 chain is incorporated
  haploinsufficiency = sufficient_evidence  [clingen_dosage@2026-01/9ba328295043]
  loeuf = 0.105  [gnomad_constraint@v2.1.1/118d41eba7e8]
```

A fact from ClinGen says ClinGen. A fact we decided says we decided. Every source is pinned by the digest of its bytes, which is what turns `gene_curation@2026-01` from a label into something checkable.

### Which interventions the mechanism admits

M5 answers *why*. M6 answers *what class of intervention that mechanism allows*, and stops there.

**A contraindication outranks any number of indications.** Reasons to try something accumulate; one reason not to ends the matter. The report prints the ruled-out list first, in defiance of what a reader would rather see, because the output that changes what somebody does is the one that closes a door.

**An unresolved mechanism blocks everything.** If M5 returned `undetermined`, no modality is assessed — the alternative is rules firing on transcript facts alone and producing a plausible list of interventions resting on nothing.

**Half of it depends on the patient, not the variant.** A mechanism is a property of a variant and a gene; how much functional product the patient has left is a property of the patient. No upregulation without an intact allele, no allele-specific silencing of the only copy, no X reactivation for a hemizygous male. Unknown zygosity is treated as neither answer — the verdicts stand and carry a caveat naming exactly which of them depend on the gap.

### Designing the molecule

Three designers, each producing a placement rather than a recommendation.

**Base editing** — which editor, which protospacer, on which strand, and every unintended base that comes with it. Designed against the *patient's* sequence, because a guide matching the reference is a guide matching the healthy allele.

**Prime editing** — pegRNA geometry: the nick three bases from the PAM, the primer binding site read upstream of it, the reverse-transcription template written downstream, and the second nick classified PE3 or PE3b.

**Antisense oligonucleotides** — a tiling of the target with composition, motif and placement rules: G-quadruplex risk, CpG content in a phosphorothioate backbone, and whether the window covers anything worth covering.

Two absences are objects rather than omissions. **No efficiency model is attached** — no BE-Hive, no PRIDICT, no folding model — and `NoModelAttached` refuses in the same shape a real model would answer, so attaching one later is a constructor argument rather than a rewrite. What ordering exists is declared: fewest bystanders first, then position, which is the difference between a criterion and a heuristic. **No CFD table is attached** to the off-target reader either; hits are ranked by *where they land* — coding sequence in an essential gene, an oncogene, a gene expressed in the target tissue — and the mismatch count is reported as the search gave it.

### The answer changes without anybody learning anything

The genome does not change. What changes is everything the rules consulted while reading it.

A mechanism call is not a fact about a variant. It is a fact about a variant **and a world**: eight pinned coordinates covering ClinVar, population frequencies, gene curation, the panel, the patient's phenotype, the transcript annotation, tissue expression and our own rule files. Name the world and two things follow — a call becomes reproducible, and the difference between two calls becomes attributable.

**Attribution is established by experiment, not by correlation.** When several axes move in one week, the module re-runs the rules on counterfactual worlds: the old world with exactly one axis advanced, the new world with exactly one held back. Each moved axis comes back *decisive*, *sufficient*, *necessary* or *contributing*.

**Our own corrections are not discoveries.** A change attributable only to a rule-file edit goes to a validation queue, not to clinical sign-out. It may become a finding once someone confirms the new rule is right — the confirmation comes first.

**Annotation is an axis, and it is the sharpest case here.** A new transcript release moves the exon boundaries. The same stop codon that sat upstream of the last junction now sits within fifty nucleotides of it, the transcript escapes decay, a truncated protein is made, and a collagen-like gene goes from haploinsufficiency to dominant negative. Nobody learned anything about the patient. Gene addition goes from indicated to contraindicated anyway.

### The queue, as a page somebody opens

A scheduled reanalysis runs at three in the morning and exits, so `repairbench review` serves the work list: which cases moved, what the change was, which queue it went to, and one button to sign it off.

The design question is not what to show — it is **what an empty page means**. An empty queue is the commonest correct output this system produces, and it is exactly what a scheduler that died in March looks like. So the loudest element on the page is not an event; it is a case nobody has examined recently.

Nothing on the page is recomputed. Every urgency and queue is reproduced as the run recorded it, so the page cannot promote an event or quietly demote one. The server can change exactly one thing — whether an event is marked read, and by whom — and acknowledging requires a name, because it suppresses every future alert with that fingerprint and an anonymous suppression is one nobody can be asked about later.

## Three things it refuses to do

**It will not return a mechanism without evidence.** There is no constructor for one. A determined call always carries the rules that produced it, and a test asserts that across the whole reference set.

**It will not hide a disagreement.** A curated mechanism from an expert panel outranks anything inferred — but the other rules still run, and any that disagree are recorded as conflicts. *MECP2* forces this: most of its coding sequence sits in the last exon, so common truncating variants escape decay and a protein *is* made, yet the curated mechanism is loss of function. The call follows the curation and reports the escape, because a reviewer needs both facts.

**It will say it does not know.** `undetermined` is a first-class answer with its own rules — "the dosage curation was refuted, not merely absent", "a truncated protein is produced and nothing here says whether it is inert or harmful", "four reported missense variants is not a distribution". A wrong mechanism points at a wrong therapy; an honest hedge costs nothing.

## Validation

Three kinds, in increasing order of how hard they are to fake.

**Reference sets — 26 cases.** Genes whose mechanism the field has settled, diseases where a route was actually taken, and reanalysis episodes of the shape a laboratory meets. The reference set is the specification and the rule file is the implementation: when they disagree, the first question is which rule is wrong.

```
$ repairbench reference

ok    SCN1A nonsense, Dravet syndrome                  loss_of_function (probable)
ok    COL1A1 nonsense, OI type I                       loss_of_function (probable)
ok    COL1A1 glycine substitution, OI type II-IV       dominant_negative (probable)
ok    PIK3CA hotspot missense                          gain_of_function (possible)
ok    MECP2 nonsense in the final exon, Rett           loss_of_function (established)
ok    dosage sensitivity refuted                       undetermined (none)
...
11/11 reference cases reproduced under mechanism-v1@1594adb029e2
```

**Real releases.** Every input is a public file, pinned by content digest: RefSeq GRCh38 annotation, UCSC sequence, ClinGen dosage curation, gnomAD constraint, GTEx expression, ClinVar submissions. Nothing in the reference set is invented any more — the exon structures, the constraint values, the variant counts and the CDS positions all came from somewhere else.

**Published molecules — 3 of 3.** The strongest check, because a molecule somebody manufactured cannot be argued with:

| | reproduced from GRCh38 |
| --- | --- |
| **eteplirsen** (EXONDYS 51) | 30-mer against *DMD* exon 51, base for base as the FDA label prints it |
| ***FAH* c.1062+5G>A** | the published ABE7.10 correction guide — editor, PAM, strand, position 5; the one differing base turned out to be a bystander this package predicts |
| ***HEXA* c.1274_1277dup** | the pegRNA from the paper that introduced prime editing: spacer, nick, PBS, RT template and PE3b nicking guide, identical |

Two of those three designers were **wrong** when first asked. Both had full test modules that passed. What a test cannot do is disagree with the package about what the answer is — a fixture's expected value is written by whoever writes the fixture.

[`FINDINGS.md`](FINDINGS.md) is the record: the defects real data and real molecules exposed, what each one was, and what changed.

## Running it

```bash
pip install -e ".[dev]"

make check                            # lint, types, 495 tests, all three reference sets
repairbench rules                     # what the system believes, with citations
repairbench plan case.yaml --fasta ref.fa --annotation ref.gff3
repairbench demo --state /tmp/try     # the reanalysis loop, in a browser
```

[`RUNNING.md`](RUNNING.md) is the tour: setup, one command per layer, and what each output is claiming. [`scripts/README.md`](scripts/README.md) covers the reference data — about 500 MB, none of it ours to redistribute, and `refdata/` is gitignored.

## What is missing

* **Recurrence is unread.** The clustering measure counts alleles, and what a hotspot means is that many patients carry the same one. The variant summary does not carry observation counts.
* **One window, and hotspots come in fours.** The clustering statistic takes the densest single window; a multi-domain protein has several.
* **Splice rules are thin.** They classify the mechanism but do nothing with the consequence.
* **Tissue is bulk, adult and post-mortem.** GTEx does not answer which cell type, or what was happening in the developmental window that produced the phenotype.
* **The reanalysis episodes are synthetic.** The shapes are real; the four files each runs against are ours, because no laboratory publishes that set.
* **The seam is one-way.** A modality with no workable candidate should arguably weaken its own verdict, and today it does not.
* **Delivery is not assessed at all.** Whether a vector, an oligonucleotide or an editor reaches the affected tissue at a useful dose decides most of these in practice.
* **Three designers, three absent models.** Obtaining and validating the weights is the actual work; until then every list is grouped rather than ranked.
* **No silent PAM edits.** When the PAM survives a prime edit the report says so and stops there.
* **Splice-regulatory elements are unread.** ISS-N1, the element nusinersen was built around, is exactly what a steric blocker aims at, and no annotation here contains one.
* **Three molecules is three, not a benchmark.** Each designer has been checked against one published object. The systematic version — a set per modality, run as a suite, disagreements counted — has not been done.

## Sources

RefSeq (O'Leary et al. 2016) · UCSC hg38 (Kent et al. 2002) · ClinGen (Rehm et al. 2015) · gnomAD (Karczewski et al. 2020) · GTEx (GTEx Consortium 2020) · ClinVar (Landrum et al. 2018) · prime editing (Anzalone et al. 2019) · base editing (Komor et al. 2016; Gaudelli et al. 2017)

MIT licensed. Not affiliated with any of the above.
