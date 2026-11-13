#!/bin/bash
#SBATCH --job-name=phi_test    ### job name
#SBATCH --cpus-per-task=16     ### cores
#SBATCH --mem-per-cpu=1G       ### memory per core
#SBATCH --time=08:00:00        ### walltime
#SBATCH --output /dev/null     ### stdout (logged by tee)
#SBATCH --error /dev/null      ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
PHYLO="${ROOT}/02_structure_phylogeny"
RESULTS="${ROOT}/03_delimitation/results/021_results_phi_test"
LOGS="${ROOT}/03_delimitation/logs/021_logs_phi_test"
### from local/03_delimitation/001_lineage_labels.R (copy local/03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv and lineage_labels.tsv.md5 to ${ROOT}/03_delimitation/supplementary)
LABELS_DIR="${ROOT}/03_delimitation/supplementary"
### frozen labels (one representative per MLG: clustered = TRUE)
LABELS="${LABELS_DIR}/lineage_labels.tsv"
### from local/01_variant_calling/103_mating_type.R (copy local/01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv to ${ROOT}/03_delimitation/supplementary)
META="${LABELS_DIR}/sample_metadata_ALL_FIXED.tsv"
### from hpc/02_structure_phylogeny/042_rename_headers.sh (per-gene alignments)
ALN_DIR="${PHYLO}/results/042_results_rename_headers/per_og_renamed"

UNIT_DIR="${RESULTS}/units"
JOBS="${RESULTS}/jobs.tsv"
UNIT_LEVELS="${RESULTS}/unit_levels.tsv"
PER_GENE="${RESULTS}/phi_per_gene.tsv"
SUMMARY="${RESULTS}/phi_summary.tsv"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-16}"
### levels whose units are tested
K_LEVELS="k4 k5 k6"
### mating type calls of possible mixed cultures (sample set no_mixed)
MAT_EXCLUDE="MAT12,ambiguous"
### random draws at the size of the smallest unit (sample sets equal_n_1 ...)
N_DRAWS=5
SEED=42
### PHI window in bp (PhiPack and SplitsTree default)
WINDOW=100
### permutations
N_PERM=1000
### C locale (no locale warnings from R, fixed sort order)
export LC_ALL=C

### create results folders
mkdir -p "${RESULTS}" "${UNIT_DIR}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_phi_test.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_R

START_TIME=$(date +%s)

### run info
echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"

## check if the labels are the frozen version
( cd "${LABELS_DIR}" && md5sum -c lineage_labels.tsv.md5 ) \
  || { echo "ERROR: lineage_labels.tsv does not match its md5"; exit 1; }

## check if the inputs are there
for F in "${LABELS}" "${META}"; do
    if [[ ! -s "${F}" ]]; then
        echo "ERROR: input not found: ${F}"
        exit 1
    fi
done

## check if the alignments are there
N_GENES=$(ls "${ALN_DIR}"/*_final_alignment.fasta | wc -l)
if [[ "${N_GENES}" -eq 0 ]]; then
    echo "ERROR: no alignments in ${ALN_DIR}"
    exit 1
fi

### temporary folder
TMP="${RESULTS}/temp_${SLURM_JOB_ID:-local}"
mkdir -p "${TMP}"


### =====================================================================
### step 1: sample sets per unit
### =====================================================================

echo ""
echo "=== Step 1: Sample sets per unit ==="

### one representative per MLG; sets all, no_mixed and equal_n draws; units of each level
Rscript - "${LABELS}" "${META}" "${K_LEVELS}" "${MAT_EXCLUDE}" "${N_DRAWS}" "${SEED}" \
    "${UNIT_DIR}" "${UNIT_LEVELS}" <<'RSCRIPT'
a <- as.list(commandArgs(trailingOnly = TRUE))
names(a) <- c("labels", "meta", "k_levels", "mat_exclude", "n_draws", "seed", "out_dir",
              "out_levels")

k_levels <- strsplit(a$k_levels, " ")[[1]]

lab  <- read.delim(a$labels, colClasses = "character")
meta <- read.delim(a$meta, colClasses = "character")
rep  <- lab[lab$clustered == "TRUE", ]

mixed <- meta$sample[meta$mat_genome %in% strsplit(a$mat_exclude, ",")[[1]]]

### members of every unit (a unit name stands for the same samples at every level)
units   <- sort(unique(unlist(rep[k_levels])))
members <- lapply(units, function(u) {
  m <- unique(Filter(length, lapply(k_levels, function(k) sort(rep$sample[rep[[k]] == u]))))
  if (length(m) != 1) stop("unit ", u, " differs between levels")
  m[[1]]
})
names(members) <- units

### sample sets (draws at the size of the smallest unit)
n_sub <- min(lengths(members))
sets  <- list(all = members, no_mixed = lapply(members, setdiff, mixed))
set.seed(as.integer(a$seed))
for (d in seq_len(as.integer(a$n_draws))) {
  sets[[paste0("equal_n_", d)]] <- lapply(members, function(m) sort(sample(m, n_sub)))
}

for (s in names(sets)) {
  dir.create(file.path(a$out_dir, s), showWarnings = FALSE)
  for (u in units) writeLines(sets[[s]][[u]], file.path(a$out_dir, s, paste0(u, ".txt")))
}

### level and unit (a unit of several levels is tested once)
lev <- do.call(rbind, lapply(k_levels, function(k)
  data.frame(k = k, unit = sort(unique(rep[[k]])))))
write.table(lev, a$out_levels, sep = "\t", quote = FALSE, row.names = FALSE)

cat("representatives:       ", nrow(rep), "\n")
cat("mixed cultures:        ", sum(rep$sample %in% mixed), "\n")
cat("equal_n draw size:     ", n_sub, "\n")
cat("sample sets:           ", paste(names(sets), collapse = " "), "\n\n")
print(sapply(sets[c("all", "no_mixed")], lengths))
RSCRIPT

SETS=$(ls "${UNIT_DIR}")
UNITS=$(ls "${UNIT_DIR}/all" | sed 's/\.txt$//')


### =====================================================================
### step 2: PHI per sample set, unit and gene
### =====================================================================

echo ""
echo "=== Step 2: PHI per sample set, unit and gene ==="

conda deactivate
conda activate cg_bioinfo_v1_seqkit
echo "seqkit version: $(seqkit version)"
SEQKIT=$(command -v seqkit)

conda deactivate
conda activate cg_bioinfo_v1_phipack
echo "PhiPack: $(conda list phipack | tail -1)"
PHI=$(command -v Phi)

## check if every sample of a unit is in the alignments
FIRST=$(ls "${ALN_DIR}"/*_final_alignment.fasta | head -1)
"${SEQKIT}" seq -n -i "${FIRST}" | sort > "${TMP}/aln_samples.txt"
N_MISS=$(cat "${UNIT_DIR}"/all/*.txt | sort -u | comm -23 - "${TMP}/aln_samples.txt" | wc -l)
if [[ "${N_MISS}" -gt 0 ]]; then
    echo "ERROR: ${N_MISS} unit samples are not in the alignments"
    exit 1
fi

echo "Genes:        ${N_GENES}"
echo "Window:       ${WINDOW} bp"
echo "Permutations: ${N_PERM}"

### one job per sample set, unit and gene
for S in ${SETS}; do
    for U in ${UNITS}; do
        ls "${ALN_DIR}"/*_final_alignment.fasta | sed "s|^|${S} ${U} |"
    done
done > "${JOBS}"

echo "PHI runs:     $(wc -l < "${JOBS}")"

### function: run_phi (unit subset of one gene, PHI with permutation test; PhiPack writes into the working directory)
run_phi () {
    local set="$1" unit="$2" aln="$3"
    local gene d n_seq out n_inf p_perm p_norm
    gene=$(basename "${aln}" _final_alignment.fasta)
    d=$(mktemp -d -p "${TMP}")
    "${SEQKIT}" grep --quiet -f "${UNIT_DIR}/${set}/${unit}.txt" "${aln}" > "${d}/aln.fasta"
    n_seq=$(grep -c ">" "${d}/aln.fasta")
    out=$(cd "${d}" && "${PHI}" -f aln.fasta -w "${WINDOW}" -p "${N_PERM}" 2>&1)
    n_inf=$(echo "${out}" | grep -m1 "informative sites" | grep -o "[0-9]\+")
    p_perm=$(echo "${out}" | awk '/PHI \(Permutation\)/ {print $3}')
    p_norm=$(echo "${out}" | awk '/PHI \(Normal\)/ {print $3}')
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "${set}" "${unit}" "${gene}" "${n_seq}" \
        "${n_inf:-NA}" "${p_perm:-NA}" "${p_norm:-NA}"
    rm -rf "${d}"
}
export -f run_phi
export SEQKIT PHI TMP UNIT_DIR WINDOW N_PERM

{
    printf "set\tunit\tgene\tn_seq\tn_informative\tp_perm\tp_normal\n"
    xargs -P "${THREADS}" -n 3 bash -c 'run_phi "$1" "$2" "$3"' _ < "${JOBS}" \
        | sort -k1,1 -k2,2 -k3,3
} > "${PER_GENE}"


### =====================================================================
### checks
### =====================================================================

echo ""
echo "=== Checks ==="

## check if every sample set, unit and gene has a row
N_ROWS=$(( $(wc -l < "${PER_GENE}") - 1 ))
N_JOBS=$(wc -l < "${JOBS}")
echo "Rows: ${N_ROWS}   jobs: ${N_JOBS}"
[[ "${N_ROWS}" -eq "${N_JOBS}" ]] || echo "WARNING: rows and jobs differ"

## check if every gene subset holds all samples of its unit
for S in ${SETS}; do
    for U in ${UNITS}; do
        N_U=$(wc -l < "${UNIT_DIR}/${S}/${U}.txt")
        N_SHORT=$(grep -P "^${S}\t${U}\t" "${PER_GENE}" | cut -f4 | grep -vx "${N_U}" | wc -l)
        [[ "${N_SHORT}" -eq 0 ]] || echo "WARNING: ${S} ${U}: ${N_SHORT} genes hold fewer than ${N_U} samples"
    done
done

### temporary folder removed
rm -rf "${TMP}"


### =====================================================================
### step 3: summary
### =====================================================================

echo ""
echo "=== Genes with PHI p < 0.05 per level and unit (5% expected without recombination) ==="

conda deactivate
conda activate cg_bioinfo_v1_R

Rscript - "${PER_GENE}" "${UNIT_LEVELS}" "${SUMMARY}" <<'RSCRIPT'
a <- as.list(commandArgs(trailingOnly = TRUE))
names(a) <- c("per_gene", "levels", "out")

d <- read.delim(a$per_gene, na.strings = c("NA", "--"))

### per sample set and unit
s <- do.call(rbind, lapply(split(d, list(d$set, d$unit), drop = TRUE), function(x) {
  t <- x[!is.na(x$p_perm), ]
  data.frame(set = x$set[1], unit = x$unit[1], n_seq = x$n_seq[1], genes = nrow(x),
             tested = nrow(t), median_informative = median(t$n_informative),
             pct_p05 = round(100 * mean(t$p_perm < 0.05), 1),
             n_fdr05 = sum(p.adjust(t$p_perm, "BH") < 0.05))
}))
s <- merge(read.delim(a$levels), s, by = "unit")
s <- s[order(s$k, s$unit, s$set), c("k", "unit", setdiff(names(s), c("k", "unit")))]
write.table(s, a$out, sep = "\t", quote = FALSE, row.names = FALSE)

### one table per level: all, no_mixed and the mean of the equal_n draws
s$group <- sub("_[0-9]+$", "", s$set)
w <- aggregate(pct_p05 ~ k + unit + group, s, mean)
w <- reshape(w, idvar = c("k", "unit"), timevar = "group", direction = "wide")
names(w) <- sub("pct_p05\\.", "pct_", names(w))
for (k in unique(w$k)) {
  cat("\n", k, "\n", sep = "")
  print(w[w$k == k, -1], row.names = FALSE, digits = 3)
}
RSCRIPT

### output files (copy phi_per_gene.tsv and unit_levels.tsv to local/03_delimitation/data/021_data_phi_test for local/03_delimitation/021_phi_test.R; units/ and unit_levels.tsv to local/03_delimitation/data/022_data_index_association for local/03_delimitation/022_index_association.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  units/<set>/<unit>.txt  samples per sample set and unit (K4, K5 and K6 units)"
echo "  unit_levels.tsv         units of each level"
echo "  jobs.tsv                sample set, unit and alignment per PHI run"
echo "  phi_per_gene.tsv        informative sites and PHI p values per sample set, unit and gene"
echo "  phi_summary.tsv         genes tested and genes with PHI p < 0.05 per level, unit and sample set"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "PHI tests completed in ${RUNTIME} seconds ($((RUNTIME / 60)) min)."
