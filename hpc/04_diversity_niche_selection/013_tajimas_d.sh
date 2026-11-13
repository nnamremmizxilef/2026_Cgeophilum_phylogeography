#!/bin/bash
#SBATCH --job-name=tajimas_d    ### job name
#SBATCH --cpus-per-task=2       ### cores
#SBATCH --mem-per-cpu=4G        ### memory per core
#SBATCH --time=06:00:00         ### walltime
#SBATCH --output /dev/null      ### stdout (logged by tee)
#SBATCH --error /dev/null       ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
PHYLO="${ROOT}/02_structure_phylogeny"
RESULTS="${ROOT}/04_diversity_niche_selection/results/013_results_tajimas_d"
LOGS="${ROOT}/04_diversity_niche_selection/logs/013_logs_tajimas_d"
### from hpc/04_diversity_niche_selection/010_sequence_distances.sh (representatives, units and mixed-MAT flag)
CLASSES="${ROOT}/04_diversity_niche_selection/results/010_results_sequence_distances/classes.tsv"
### from hpc/02_structure_phylogeny/042_rename_headers.sh (per-gene alignments)
ALN_DIR="${PHYLO}/results/042_results_rename_headers/per_og_renamed"

GROUPS_TSV="${RESULTS}/groups.tsv"
SUMS="${RESULTS}/locus_group_sums.tsv.gz"

### parameters
### levels
K_LEVELS="k4 k5 k6"
### sample size every site is projected to (sites with fewer called sequences in a unit are dropped for that unit; the smallest unit without mixed-MAT samples)
N_PROJ=16
### C locale (fixed sort order)
export LC_ALL=C

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_tajimas_d.log"
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

## check if the inputs are there
if [[ ! -s "${CLASSES}" ]]; then
    echo "ERROR: input not found: ${CLASSES} (run 010_sequence_distances.sh first)"
    exit 1
fi

## check if the alignments are there
N_GENES=$(ls "${ALN_DIR}"/*_final_alignment.fasta | wc -l)
if [[ "${N_GENES}" -eq 0 ]]; then
    echo "ERROR: no alignments in ${ALN_DIR}"
    exit 1
fi
echo "Alignments: ${N_GENES}"


### =====================================================================
### step 1: groups (unit of a level and sample set)
### =====================================================================

echo ""
echo "=== Step 1: Groups ==="

### every unit of every level, for all representatives and without mixed-MAT samples (identical member sets, e.g. L1 at K = 4, 5 and 6, computed once)
Rscript - "${CLASSES}" "${K_LEVELS}" "${N_PROJ}" "${GROUPS_TSV}" <<'RSCRIPT'
a <- as.list(commandArgs(trailingOnly = TRUE))
names(a) <- c("classes", "k_levels", "n_proj", "out")

k_levels <- strsplit(a$k_levels, " ")[[1]]
cls <- read.delim(a$classes, colClasses = "character")
cls$mixed <- as.logical(cls$mixed)

sets <- list(all = rep(TRUE, nrow(cls)), no_mixed = !cls$mixed)

grp <- do.call(rbind, lapply(k_levels, function(k) do.call(rbind, lapply(names(sets), function(s) {
  u <- sort(unique(cls[[k]][sets[[s]]]))
  data.frame(k = k, set = s, unit = u,
             members = sapply(u, function(x) paste(sort(cls$sample[sets[[s]] & cls[[k]] == x]), collapse = ",")),
             n = sapply(u, function(x) sum(sets[[s]] & cls[[k]] == x)))
}))))
grp$group <- paste0("g", match(grp$members, unique(grp$members)))

## check if every group is at least as large as the projection
if (any(grp$n < as.integer(a$n_proj))) stop("groups smaller than the projected sample size: ",
                                            paste(unique(grp$unit[grp$n < as.integer(a$n_proj)]), collapse = ", "))

write.table(grp[, c("k", "set", "unit", "n", "group", "members")], a$out, sep = "\t", quote = FALSE, row.names = FALSE)

cat("groups (level, set, unit):", nrow(grp), "   distinct member sets:", length(unique(grp$group)), "\n")
print(grp[, c("k", "set", "unit", "n", "group")], row.names = FALSE)
RSCRIPT

## check if the groups were written
if [[ ! -s "${GROUPS_TSV}" ]]; then
    echo "ERROR: no groups written"
    exit 1
fi


### =====================================================================
### step 2: per locus and group, sites, pairwise diversity and segregating sites
### =====================================================================

echo ""
echo "=== Step 2: Diversity and segregating sites per locus and group ==="

### per site and group: called sequences n (A, C, G, T), pairwise diversity from all n, and the probability that a random subsample of N_PROJ sequences is polymorphic (hypergeometric projection; sites with fewer than N_PROJ called sequences in a group left out for that group)
Rscript - "${ALN_DIR}" "${CLASSES}" "${GROUPS_TSV}" "${N_PROJ}" "${SUMS}" <<'RSCRIPT'
suppressMessages(library(ape))

a <- as.list(commandArgs(trailingOnly = TRUE))
names(a) <- c("aln_dir", "classes", "groups", "n_proj", "out")

m   <- as.integer(a$n_proj)
cls <- read.delim(a$classes, colClasses = "character")
grp <- read.delim(a$groups, colClasses = "character")
grp <- grp[!duplicated(grp$group), c("group", "members")]

### membership of the representatives in the distinct groups
memb <- sapply(strsplit(grp$members, ","), function(s) as.numeric(cls$sample %in% s))
colnames(memb) <- grp$group

files <- list.files(a$aln_dir, pattern = "_final_alignment\\.fasta$", full.names = TRUE)

con <- gzfile(a$out, "w")
writeLines("locus\tgroup\tsites\tpi_sum\tseg_proj\tseg_raw", con)

t0 <- Sys.time()
for (i in seq_along(files)) {
  x <- read.FASTA(files[i])
  if (!all(cls$sample %in% names(x))) stop("representatives missing from ", basename(files[i]))
  ch <- as.character(as.matrix(x[cls$sample]))

  ### base counts per group and site
  cnt <- lapply(c("a", "c", "g", "t"), function(b) crossprod(memb, (ch == b) * 1))
  n   <- Reduce(`+`, cnt)
  ok  <- n >= m

  ### pairwise diversity per site (unbiased, from all called sequences)
  hom <- Reduce(`+`, lapply(cnt, function(cb) cb * (cb - 1)))
  pi_site <- 1 - hom / (n * (n - 1))

  ### probability that a subsample of m is monomorphic
  mono <- Reduce(`+`, lapply(cnt, function(cb) exp(lchoose(cb, m) - lchoose(n, m))))
  raw  <- Reduce(`+`, lapply(cnt, function(cb) (cb > 0) * 1)) > 1

  pi_site[!ok] <- 0
  mono[!ok]    <- 1
  raw[!ok]     <- FALSE

  locus <- sub("_final_alignment\\.fasta$", "", basename(files[i]))
  writeLines(paste(locus, grp$group, rowSums(ok), rowSums(pi_site), rowSums(1 - mono), rowSums(raw), sep = "\t"), con)

  if (i %% 500 == 0) message(i, " loci, ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
}
close(con)

cat("loci:", length(files), "   groups:", nrow(grp), "\n")
RSCRIPT

## check if every locus was written for every group
N_GROUPS=$(cut -f5 "${GROUPS_TSV}" | tail -n +2 | sort -u | wc -l)
N_ROWS=$(( $(zcat "${SUMS}" | wc -l) - 1 ))
echo "Rows: ${N_ROWS} (expected $(( N_GENES * N_GROUPS )))"
[[ "${N_ROWS}" -eq $(( N_GENES * N_GROUPS )) ]] || { echo "ERROR: rows missing"; exit 1; }

### output files (copy groups.tsv and locus_group_sums.tsv.gz to local/04_diversity_niche_selection/data/013_data_tajimas_d for local/04_diversity_niche_selection/013_tajimas_d.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  groups.tsv                 level, sample set, unit, size and group of every unit"
echo "  locus_group_sums.tsv.gz    per locus and group: sites, summed pairwise diversity,"
echo "                             expected segregating sites at ${N_PROJ} sequences, raw segregating sites"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Tajima's D input completed in ${RUNTIME} seconds ($((RUNTIME / 3600))h $(((RUNTIME % 3600) / 60))m)."
