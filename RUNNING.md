# Running and testing repairbench

Two ways to read this file. **Setup** and **The one command** are what you need
to convince yourself the thing works. Everything after that is a tour: one
section per layer, each with a command that produces output you can check
against what the section claims.

Nothing here needs the reference data except the sections marked **needs
`refdata/`**. The suite runs without it and says which tests it skipped.

---

## Setup

Python 3.11 or newer. The package itself depends on one library, PyYAML —
everything heavier belongs behind an interface and is not installed.

```bash
cd ~/bioinformatics/repairbench
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

The editable install puts a `repairbench` command on the path. If you would
rather not install anything, every command below also works as
`PYTHONPATH=src python -m repairbench.cli …`, and `make` sets that itself.

A note on the Python you use. The `python3` on macOS may be a version that has
nothing installed into it — that is what produced `ModuleNotFoundError: No
module named 'repairbench'` the first time. Inside the virtual environment above
this cannot happen, because `pip install -e` puts the package where that
interpreter looks.

---

## The one command

```bash
make check
```

Four things in sequence, and it fails on the first that does not hold:

| step | what it proves |
| --- | --- |
| `ruff check src tests` | style, imports, a few real bug classes |
| `mypy` (strict) | every function annotated, no `Any` leaking |
| `pytest tests` | 424 tests |
| three reference sets | 11 mechanism cases, 10 modality cases, 5 reanalysis episodes |

Expect **422 passed, 2 skipped in ~55 s** with the reference data present, and
**414 passed, 10 skipped in ~9 s** without it. The skips are the tests that read
a real genome or a real ClinVar release; they name what is missing rather than
passing quietly.

Individually, when one of them is what you want:

```bash
make test        # just the suite
make lint        # ruff + mypy, exactly as CI runs them
make reference   # the three reference sets, printed case by case
make rules       # every rule, what it claims, and what it cites
```

`make reference` is the one worth reading rather than watching. It prints each
case, the mechanism the rules produced, and the modalities that follow — so a
disagreement shows up as a line, not as a stack trace.

---

## What the rules believe

```bash
repairbench rules
```

Every rule in `rules/mechanism-v1.yaml`, with its strength, its reasoning in a
sentence, and its citation. This is the file to argue with. If a conclusion the
package reaches is wrong, the error is here rather than in the code, and this
command is how you find which line to change.

---

## One case, from the outside in

There is a case file in the repository that exercises the whole seam against
synthetic sequence, so it needs no downloads:

```bash
repairbench explain tests/data/design/case.yaml   # why this variant causes disease
repairbench assess  tests/data/design/case.yaml   # …and which interventions that admits
repairbench plan    tests/data/design/case.yaml \
    --fasta tests/data/design/target.fa \
    --annotation tests/data/design/target.gff3    # …and the actual molecules
```

Read `plan`'s output in three parts. **`mechanism`** with the rules that fired
underneath it. **`not designed, ruled out by the modality rules`** — the list
that matters most, because nothing on it gets a molecule designed for it however
well the sequence would work. **`designed`**, where each candidate carries what
is wrong with it rather than a score.

The case is a dominant negative, so gene addition is ruled out and
allele-specific silencing is indicated. That is the property to check: change
`truncating_variants_are_milder` to `false` in the case file and watch the
mechanism, and then the whole plan, move.

---

## The designers on their own

```bash
repairbench design  --gene TARG --at 17:301 --patient A --wild-type G \
    --fasta tests/data/design/target.fa                       # base editing
repairbench pegrna  --gene TARG --at 17:301 --patient A --wild-type G \
    --fasta tests/data/design/target.fa                       # prime editing
repairbench aso     --gene TARG --at 17:280-330 --chemistry gapmer-2MOE \
    --fasta tests/data/design/target.fa                       # antisense
repairbench offtarget tests/data/design/hits.txt              # rank a hit list
```

Each of the three refuses in a different way, and the refusals are the point:
`design` will tell you when no editor's window reaches the base, `pegrna` when
the PAM survives the edit, `aso` when a window has nothing to discriminate on.

---

## Gene facts, and where each one came from

```bash
repairbench context COL1A1
```

Every fact the rules read about a gene, each on its own line with the pin of the
file it came from. Two of the five have no public table and live in a local
curation file that demands a citation per entry — those lines say so.

**needs `refdata/`** — with the real releases instead of the fixtures:

```bash
repairbench context COL1A1 PIK3CA \
    --dosage refdata/ClinGen_gene_curation_list_GRCh38.tsv \
    --constraint refdata/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz \
    --expression refdata/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz \
    --clinvar refdata/clinvar_refset.tsv.gz
```

---

## Counting pathogenic variation

**needs `refdata/`**

```bash
repairbench clinvar refdata/clinvar_refset.tsv.gz \
    COL1A1 DMD KRT14 MECP2 PIK3CA SCN1A SCN2A SMN1 UBE3A --release 2026-08
```

Prints each gene's `distribution:` block in the shape `tests/reference/*.yaml`
wants, the submission count broken down by review stars, and real c./p.
positions on the transcript the submitters used. This is how the reference set
gets updated after a ClinVar release — transcription rather than judgement.

Point it at `variant_summary.txt.gz` instead and it reads the whole release;
that takes about a minute rather than a second and gives the same answer, since
the extract is a superset of what the parser keeps.

---

## The real locus, end to end

**needs `refdata/`** — the point of the whole download.

```bash
repairbench plan tests/data/real/col1a1-gly821ser.yaml \
    --fasta refdata/chr17.fa \
    --annotation refdata/GRCh38_latest_genomic.gff.gz \
    --limit 3
```

About 30 seconds, most of it parsing 1.5 GB of annotation down to one gene.
*COL1A1* p.(Gly821Ser) is the case worth spending a real genome on: the gene is
on the minus strand, so every coordinate and complement is exercised in the
direction where mistakes hide, and the change is G>A on the coding strand and
therefore C>T on the plus strand — which is where a designer that ignores strand
produces a confident, wrong molecule.

Two values in the output are hand-verified against GRCh38 and asserted in
`tests/test_real_locus.py`: the protospacer `GACAGCCAACCTGGTGCTAA` and the
oligonucleotide `AGCCAACCTGGTGCTAAAGG`. If a refactor moves either, that test
goes red.

```bash
repairbench annotation refdata/GRCh38_latest_genomic.gff.gz --gene COL1A1
```

Which transcripts exist for a gene, which is MANE Select, and — the part worth
looking at when a gene comes back odd — which transcripts the file describes
that this package cannot use, with the reason for each. Drop `--gene` and it
reads the whole annotation: about four minutes, 136,000 transcripts, and three
rejections that are real biology rather than a broken file.

---

## Reproducing an approved drug

**needs `refdata/`** — chrX and the annotation.

```bash
repairbench aso --gene DMD --at chrX:31773960-31774192 \
    --chemistry steric-PMO-30 --strand - --fasta refdata/chrX.fa --limit 204
```

Somewhere in those 204 windows is `CTCCAACATCAAGGAAGATGGCATTTCTAG` at
chrX:31774098-31774127 — eteplirsen, base for base as the FDA label prints it.
`tests/test_published_molecules.py` asserts it, and the same file records the
strand defect that reproducing it exposed.

Note `--strand`, which is required and has no default. An antisense
oligonucleotide is complementary to the messenger, so which genomic strand it
copies depends on the gene's orientation — and assuming one silently is what
made every oligonucleotide this package printed for a minus-strand gene a copy
of its own target.

Two more reproductions run in the same file — the *FAH* base-editing guide and
the *HEXA* pegRNA from the paper that introduced prime editing. All three live
in `tests/test_published_molecules.py`; it is the slowest file in the suite
because it reads the real annotation, and it is the only one that can disagree
with the package about what the answer is.

## Reanalysis, which is the part that runs unattended

A reanalysis run starts, compares today with last month, and exits. Cron owns
the schedule; the process owns one comparison.

```bash
repairbench watch DEMO-1 \
    --state /tmp/rbstate \
    --catalogue tests/data/deployment/catalogue.yaml \
    --vcf tests/data/deployment/case.vcf --sample CHILD

repairbench reanalyse DEMO-1 \
    --state /tmp/rbstate \
    --catalogue tests/data/deployment/catalogue.yaml
```

The queue across every watched case, as a page:

```bash
repairbench dashboard --state /tmp/rbstate               # short digest, for a terminal or an email
repairbench dashboard --state /tmp/rbstate --out queue.html   # self-contained page, opens from disk
```

Read the page top down. The cards are the counts; the red one is cases nobody
has examined recently, and it is red because an empty queue under a dead
scheduler looks exactly like an empty queue under a healthy one. Then the two
review queues, then every case including the quiet ones.

To see it with something in it, register a case under the older releases and
then let the newer ones land:

```bash
D=tests/data/deployment
repairbench watch DEMO-2 --state /tmp/rbdemo --catalogue $D/catalogue-old.yaml \
    --variant "PLUSG:stop_gained:158:heterozygous" --phenotype hpo-day-1
repairbench reanalyse DEMO-2 --state /tmp/rbdemo --catalogue $D/catalogue-old.yaml
repairbench reanalyse DEMO-2 --state /tmp/rbdemo --catalogue $D/catalogue.yaml
repairbench dashboard --state /tmp/rbdemo --out /tmp/queue.html
```

The second reanalysis is the one that moves: ClinGen refutes dosage sensitivity
for the gene between the two releases, a settled loss-of-function call stops
resolving, and the change routes to clinical sign-out at high urgency.

`watch` reads the patient's VCF, says how many alleles it carried and which ones
it could not place on a transcript and why, and pins the world it registered
against. `reanalyse` re-examines the case against the current world and reports
what moved, along which of the eight drift axes, and whether it needs a person.
Run twice in a row it says *nothing moved*, which is the answer it gives most
often and the one it has to give honestly.

To see it actually move something, point the catalogue at the second release of
a source — `tests/data/deployment/` holds two ClinGen releases, two gnomAD
releases, two annotations and two rule versions for exactly this.

```bash
repairbench serve --addr :9090     # /health and /metrics between runs
```

### Everything from the browser

To try it, one command — and then nothing else in the terminal:

```bash
repairbench demo --state ~/rb-demo
```

It seeds a synthetic case assessed against last month's releases, loads this
month's, and serves the page. Open `http://127.0.0.1:8080`, press **Re-examine
every case**, and a HIGH change appears in the clinical sign-out queue: ClinGen
refuted dosage sensitivity for the gene, and a settled loss-of-function call
stopped resolving. Click the case, type your name and a note, press Acknowledge,
go back — the queue is empty and the page says an empty queue means nothing
moved rather than nothing ran.

Why the demo has to seed anything: drift needs a case assessed against *older*
releases and *newer* ones to compare with, and a fresh directory has neither.
Everything the seed does is available from the page — it is a shortcut past the
month of waiting, not past the interface.

For real use, point it at your own state and catalogue:

```bash
repairbench review --state ~/rb-state --catalogue ~/releases/catalogue.yaml
# http://127.0.0.1:8080
```

The catalogue is deployment configuration, like a database URL: whatever fetches
releases appends to it, and the server re-reads it on every run. A reviewer
never restarts anything.

The page carries the queue, a form to start watching a new case, a button to
re-examine every case, and a button per case. Register a patient, press
**Re-examine every case**, read what moved, open the case, sign it off. The
banner at the top says what the last action did — including *nothing moved*,
which is the commonest answer and the one worth stating.

Started without `--catalogue` the server still shows and signs off a queue that
cron fills; it says on the page that it cannot start a run, rather than hiding
the button.

Open the queue, click a case, read the change, type your name and a note, press
Acknowledge. The same fingerprint will not be raised again — which is why the
name is required and a blank one is refused.

Two things to know before you use it in anger. The name is **attribution, not
authentication**: it records who *said* they reviewed a change, which is what
the ledger needs and is not the same as proving it, so the server binds to
loopback and nothing more. And the server can change exactly one thing —
whether an event is marked as read — so a misclick cannot move a mechanism,
edit a rule or alter an urgency.

Without a browser, for a script:

```bash
repairbench acknowledge NICU-014 evt-000001 --state /tmp/rbdemo \
    --by "B. Kowalski" --note "confirmed against ClinGen, plan paused"
```

---

## Getting the reference data

```bash
./scripts/fetch-reference-data.sh            # everything, ~500 MB
./scripts/fetch-reference-data.sh clinvar    # or one section at a time
```

Sections: `annotation`, `sequence`, `curation`, `expression`, `clinvar`, `all`.
Everything lands in `refdata/`, which is gitignored — none of these files is
ours to redistribute. `scripts/README.md` says what each one settles and what
went wrong the first time each was read.

The ClinVar section writes a 3 MB extract beside the 440 MB release, filtered to
the reference-set genes, so the real-data tests run in a second.

---

## When something fails

**A reference case fails.** The first question is not "which assertion do I
relax" but "which rule is wrong". The reference set is the specification and the
rule file is the implementation, not the other way round. Two cases have already
moved this way — *PIK3CA* and *SCN2A*, both when real ClinVar counts replaced
invented ones — and both times the rule changed and the disagreement was written
into the case note.

**A test skips.** It says what is missing. `refdata/` skips are expected without
the downloads.

**`mypy` fails after an edit to a rule file.** It will not — rule files are data
and are not type-checked. What checks them is `load_ruleset`, which refuses an
unknown feature name, an unknown strength, and a rule with no reasoning. Run
`repairbench rules` after editing one; a malformed file fails there first, with
the line.
