#!/bin/bash

set -euo pipefail
shopt -s nullglob


### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
OUTDIR="${ROOT}/01_variant_calling/supplementary/010_supplementary_read_trimming"

### external: raw paired-end reads per sequencing batch and of the outgroup (copy manually into these folders)
INPUT_DIRS=(
  "${ROOT}/01_variant_calling/data/sequences/cg_igat_batch1/raw_reads"
  "${ROOT}/01_variant_calling/data/sequences/cg_igat_batch2/raw_reads"
  "${ROOT}/01_variant_calling/data/sequences/cg_igat_batch3/raw_reads"
  "${ROOT}/01_variant_calling/data/sequences/magic_igat_batch1/raw_reads"
  "${ROOT}/01_variant_calling/data/sequences/magic_igat_batch2/raw_reads"
  "${ROOT}/01_variant_calling/data/outgroups"
)

### made by hand: raw file stems to exclude, one per line (into ${ROOT}/01_variant_calling/supplementary/010_supplementary_read_trimming)
EXCLUDE_FILE="${ROOT}/01_variant_calling/supplementary/010_supplementary_read_trimming/exclude.txt"


### parameters
### prefix shared by the four lists
PREFIX="ceno_genomics_igat"

### R1 file pattern and the tags mapping an R1 path to its R2 mate
R1_GLOB="*_R1_*.fastq.gz"
R1_TAG="_R1_"
R2_TAG="_R2_"

### exclusions by extended regex on the sample name, e.g. EXCLUDE=( '^M4-' 'FI-KOX-E2E' ) (EXCLUDE_FILE: fixed substrings of the raw file name, full stems)
EXCLUDE=()

### sample name from the R1 file name (check whenever the naming changes)
###   189176_ID3825_36-FI-KOX-E2E-DNA-P1-D05_S1_L001_R1_001.fastq.gz -> FI-KOX-E2E
###   114388_ID2907_215-M4-01-P04-G03_S1_L001_R1_001.fastq.gz        -> M4-01
sample_name() {
  basename "$1" \
    | sed -E 's/_S[0-9]+_L[0-9]+_R[12]_[0-9]+\.fastq\.gz$//;
              s/^[0-9]+_ID[0-9]+_[0-9]+-//;
              s/-(DNA-)?P[0-9]+-[A-H][0-9]+$//'
}


### output lists (line N is the same sample in all four, the 010 array index)
raw_r1="${OUTDIR}/${PREFIX}_raw_reads_list_R1.txt"
raw_r2="${OUTDIR}/${PREFIX}_raw_reads_list_R2.txt"
trim_r1="${OUTDIR}/${PREFIX}_trimmed_reads_list_R1.txt"
trim_r2="${OUTDIR}/${PREFIX}_trimmed_reads_list_R2.txt"
: > "$raw_r1"; : > "$raw_r2"; : > "$trim_r1"; : > "$trim_r2"

### exclusion file (blank lines stripped so -F cannot match everything)
exclude_tmp=""
if [[ -n "$EXCLUDE_FILE" ]]; then
  ## check if the exclusion file exists
  [[ -f "$EXCLUDE_FILE" ]] || { echo "ERROR: EXCLUDE_FILE not found: $EXCLUDE_FILE" >&2; exit 1; }
  exclude_tmp=$(mktemp)
  trap 'rm -f "$exclude_tmp"' EXIT
  grep -v '^[[:space:]]*$' "$EXCLUDE_FILE" > "$exclude_tmp" || true
fi

### R1 files (sorted within each folder, folders in listed order)
mapfile -t r1_files < <(
  for d in "${INPUT_DIRS[@]}"; do
    for f in "$d"/$R1_GLOB; do printf '%s\n' "$f"; done | sort
  done
)

declare -A seen
n=0
skipped=0
for r1 in "${r1_files[@]}"; do
  base=$(basename "$r1")
  sample=$(sample_name "$r1")

  ### exclusion by regex on the sample name
  excluded=0
  for pat in "${EXCLUDE[@]:-}"; do
    [[ -n "$pat" ]] || continue
    if grep -qE "$pat" <<<"$sample"; then excluded=1; break; fi
  done

  ### exclusion by fixed substring of the raw file name
  if [[ "$excluded" -eq 0 && -s "${exclude_tmp:-/dev/null}" ]] \
     && grep -Fqf "$exclude_tmp" <<<"$base"; then
    excluded=1
  fi

  if [[ "$excluded" -eq 1 ]]; then
    echo "excluded: $sample" >&2
    skipped=$((skipped + 1))
    continue
  fi

  r2="${r1/$R1_TAG/$R2_TAG}"
  ## check if the R2 mate exists
  if [[ ! -f "$r2" ]]; then
    echo "WARNING: no R2 mate for ${base} -- skipping" >&2
    skipped=$((skipped + 1))
    continue
  fi

  ## check if the sample name is taken already (two raw files, one output name)
  if [[ -n "${seen[$sample]:-}" ]]; then
    echo "WARNING: duplicate sample name '${sample}' (${seen[$sample]} and ${base})" >&2
  fi
  seen[$sample]="$base"

  printf '%s\n' "$r1"          >> "$raw_r1"
  printf '%s\n' "$r2"          >> "$raw_r2"
  printf '%s\n' "${sample}_R1" >> "$trim_r1"
  printf '%s\n' "${sample}_R2" >> "$trim_r2"
  n=$((n + 1))
done


### output files
echo "Wrote ${n} sample(s):"
printf '  %s\n' "$raw_r1" "$raw_r2" "$trim_r1" "$trim_r2"
if [[ "$skipped" -gt 0 ]]; then echo "Skipped/excluded: ${skipped}"; fi
