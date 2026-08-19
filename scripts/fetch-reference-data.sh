#!/usr/bin/env bash
#
# Fetch the public reference data this package can be pointed at.
#
# Everything in tests/data/ is synthetic and says so. This script downloads the
# real thing, so that "the reference set reproduces the literature" can be
# checked against real transcript structures and real curations rather than
# against numbers that were originally typed in by hand.
#
# Nothing here is redistributable, which is why it is a script rather than a
# directory of files. Each source is public, free, and asks to be cited rather
# than mirrored.
#
#   ./scripts/fetch-reference-data.sh            # everything, ~500 MB
#   ./scripts/fetch-reference-data.sh annotation # one section at a time
#
# Downloads land in refdata/, which is gitignored. Total on disk after
# unpacking: roughly 1.5 GB, almost all of it the chromosome sequences.
#
# If a URL has moved — and one of them will, these files are re-issued
# constantly — the directory listing above each file is the place to look, and
# the parsers in this package read by column name rather than by position for
# exactly that reason.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${REPAIRBENCH_REFDATA:-$HERE/refdata}"
mkdir -p "$OUT"

# The chromosomes carrying the genes in tests/reference/mechanisms.yaml.
# Downloading six is about 150 MB; downloading the whole genome is 3 GB and
# buys nothing until the reference set grows.
CHROMOSOMES=(chr2 chr3 chr5 chr15 chr17 chrX)

PYTHON="${PYTHON:-python3}"

say() { printf '\n\033[1m%s\033[0m\n' "$*"; }

index_fasta() {
  # Write the .fai index without requiring the package to be installed, and
  # without requiring samtools either.
  #
  # It runs the indexer directly rather than through the command line, and the
  # reason is worth stating: the CLI module imports the whole package including
  # PyYAML, so going through it would make indexing a downloaded chromosome
  # depend on a dependency that has nothing to do with indexing. The indexer
  # itself imports nothing outside the standard library.
  local fasta="$1"
  if PYTHONPATH="$HERE/src" "$PYTHON" -c 'import sys
from repairbench.annotation.fasta import write_index
print(f"  index {write_index(sys.argv[1]).name}")' "$fasta" 2>/dev/null; then
    return 0
  fi
  if command -v samtools >/dev/null 2>&1; then
    printf '  index %s (samtools)\n' "$(basename "$fasta").fai"
    samtools faidx "$fasta"
    return 0
  fi
  cat >&2 <<MISSING

  Could not index $(basename "$fasta").

  The indexer needs Python 3.11 or newer and this package on the path. Either:

      PYTHON=/path/to/python3.11 ./scripts/fetch-reference-data.sh sequence

  or install the package and re-run:

      python3 -m venv .venv && source .venv/bin/activate
      pip install -e ".[dev]"
      ./scripts/fetch-reference-data.sh sequence

  The downloads above are already on disk, so a re-run only indexes.

MISSING
  return 1
}

get() {
  local url="$1" name="${2:-$(basename "$1")}"
  if [ -s "$OUT/$name" ]; then
    printf '  have  %s\n' "$name"
    return 0
  fi
  printf '  get   %s\n' "$name"
  # --location follows the redirects NCBI and UCSC both use; --fail turns an
  # HTML error page into a non-zero exit rather than a file full of HTML that
  # a parser would later choke on in a much more confusing place.
  curl --fail --location --progress-bar --output "$OUT/$name.part" "$url"
  mv "$OUT/$name.part" "$OUT/$name"
}

fetch_annotation() {
  say "RefSeq annotation for GRCh38  (~60 MB compressed)"
  # Listing: https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/
  #
  # Note which naming this file uses: RefSeq accessions (NC_000017.11), not
  # UCSC names (chr17). The sequence below uses the other convention, and
  # reconciling them is a real step rather than an oversight — see
  # scripts/README.md.
  get "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz"
  say "Unpacking the annotation"
  [ -s "$OUT/GRCh38_latest_genomic.gff" ] || gunzip --keep --force "$OUT/GRCh38_latest_genomic.gff.gz"
}

fetch_sequence() {
  say "Reference sequence, ${#CHROMOSOMES[@]} chromosomes  (~150 MB compressed)"
  # Listing: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/
  for chromosome in "${CHROMOSOMES[@]}"; do
    get "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/${chromosome}.fa.gz"
  done

  say "Unpacking and indexing"
  for chromosome in "${CHROMOSOMES[@]}"; do
    [ -s "$OUT/${chromosome}.fa" ] || gunzip --keep --force "$OUT/${chromosome}.fa.gz"
    # The reader will not build an index implicitly — that would mean reading
    # gigabytes at a moment nobody asked for it — so it is built here, once.
    [ -s "$OUT/${chromosome}.fa.fai" ] || index_fasta "$OUT/${chromosome}.fa"
  done
}

fetch_curation() {
  say "ClinGen dosage sensitivity  (~2 MB)"
  # Landing page, if this moves: https://search.clinicalgenome.org/kb/gene-dosage
  get "https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"

  say "gnomAD constraint  (~30 MB)"
  # Landing page: https://gnomad.broadinstitute.org/downloads#v2-constraint
  #
  # This is the v2.1.1 per-transcript file, which is what the downloads page
  # links. Two things follow from the release, and both are recorded in the
  # provenance rather than smoothed over: it calls LOEUF `oe_lof_upper` and
  # marks the row to use with `canonical` (Ensembl's pick) rather than
  # `mane_select`, and its coordinates are GRCh37 — which does not matter here,
  # because this package joins constraint to a gene by symbol and takes only the
  # ratio, never a position. The parser detects either schema from the header.
  get "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_transcript.txt.bgz"
}

fetch_clinvar() {
  say "ClinVar variant summary  (~440 MB compressed)"
  # Listing: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/
  #
  # The largest download here, and the one that removed the last invented input:
  # the variant positions in the reference set, and the counts of pathogenic
  # missense and truncating variants per gene, used to be ours. This file has
  # both, for every gene, with the review status attached — which matters,
  # because a single submitter's assertion and a reviewed expert-panel one are
  # not the same evidence and this package will not average them.
  get "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
  extract_refset_genes
}

# The reference set asks about nine genes; the file holds every gene there is.
# Reading 440 MB to answer a question about a thousandth of it takes the better
# part of a minute, and a test that costs a minute is a test that gets skipped —
# so the release is kept and an extract is written beside it.
#
# The filter is loose on purpose: it matches a symbol anywhere in the line
# rather than in the gene column, which is what lets one pass of awk do it. That
# makes the extract a *superset* of what the parser wants, and the parser
# applies the exact filter itself — so a row that slips through cannot change a
# count, it is simply read and dropped. ClinVar's own header line is kept
# verbatim, leading `#AlleleID` and all, because the reader expects it.
extract_refset_genes() {
  whole="$OUT/variant_summary.txt.gz"
  part="$OUT/clinvar_refset.tsv.gz"
  genes="COL1A1|DMD|KRT14|MECP2|PIK3CA|SCN1A|SCN2A|SMN1|UBE3A"
  [ -f "$whole" ] || return 0
  say "extracting the reference-set genes"
  gzip -dc "$whole" | awk -v pattern="$genes" 'NR == 1 || $0 ~ pattern' | gzip -c > "$part"
}

fetch_expression() {
  say "GTEx median expression  (~15 MB)"
  # Landing page: https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression
  # The file name carries the release, so this one dates fastest of the four.
  get "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
}

case "${1:-all}" in
  annotation) fetch_annotation ;;
  sequence)   fetch_sequence ;;
  curation)   fetch_curation ;;
  expression) fetch_expression ;;
  clinvar)    fetch_clinvar ;;
  all)        fetch_annotation; fetch_sequence; fetch_curation; fetch_expression; fetch_clinvar ;;
  *)          echo "usage: $0 [all|annotation|sequence|curation|expression|clinvar]" >&2; exit 2 ;;
esac

say "In $OUT"
ls -lh "$OUT" | tail -n +2 | awk '{printf "  %-64s %s\n", $9, $5}'

cat <<'NOTE'

Cite what you use. None of these is ours:

  RefSeq       O'Leary et al. 2016, Nucleic Acids Res
  UCSC hg38    Kent et al. 2002, Genome Res
  ClinGen      Rehm et al. 2015, N Engl J Med — dosage sensitivity curation
  gnomAD       Karczewski et al. 2020, Nature 581:434 — the constraint spectrum,
               which is the paper behind the v2.1.1 file this fetches
  GTEx         GTEx Consortium 2020, Science
  ClinVar      Landrum et al. 2018, Nucleic Acids Res

NOTE
