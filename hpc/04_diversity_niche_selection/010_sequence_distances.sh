#!/bin/bash
#SBATCH --job-name=sequence_distances    ### job name
#SBATCH --cpus-per-task=4                ### cores
#SBATCH --mem-per-cpu=4G                 ### memory per core
#SBATCH --time=12:00:00                  ### walltime
#SBATCH --output /dev/null               ### stdout (logged by tee)
#SBATCH --error /dev/null                ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
PHYLO="${ROOT}/02_structure_phylogeny"
RESULTS="${ROOT}/04_diversity_niche_selection/results/010_results_sequence_distances"
LOGS="${ROOT}/04_diversity_niche_selection/logs/010_logs_sequence_distances"
### from local/03_delimitation/001_lineage_labels.R (copy local/03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv and lineage_labels.tsv.md5 to ${ROOT}/04_diversity_niche_selection/supplementary)
LABELS_DIR="${ROOT}/04_diversity_niche_selection/supplementary"
### frozen labels
LABELS="${LABELS_DIR}/lineage_labels.tsv"
### from local/01_variant_calling/103_mating_type.R (copy local/01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv to ${ROOT}/04_diversity_niche_selection/supplementary)
META="${LABELS_DIR}/sample_metadata_ALL_FIXED.tsv"
### from hpc/02_structure_phylogeny/042_rename_headers.sh (per-gene alignments)
ALN_DIR="${PHYLO}/results/042_results_rename_headers/per_og_renamed"

CLASSES="${RESULTS}/classes.tsv"
EXCLUDED="${RESULTS}/samples_excluded_no_coordinates.txt"
CLASS_SUMS="${RESULTS}/class_pair_sums.tsv.gz"
PAIR_DIFFS="${RESULTS}/pair_differences.tsv.gz"
PAIR_SITES="${RESULTS}/pair_sites.tsv.gz"
LOCI="${RESULTS}/loci.tsv"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"
### levels whose units are carried along (classes are units of the finest level)
K_LEVELS="k4 k5 k6"
K_FINE=k6
### mating type calls of possible mixed cultures (flagged, not removed)
MAT_EXCLUDE="MAT12,ambiguous"
### cores for the matrix products (OpenBLAS)
export OPENBLAS_NUM_THREADS="${THREADS}"
### C locale (fixed sort order)
export LC_ALL=C

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_sequence_distances.log"
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
echo "Alignments: ${N_GENES}"


### =====================================================================
### step 1: representatives and classes
### =====================================================================

echo ""
echo "=== Step 1: Representatives and classes ==="

### one representative per MLG with sampling coordinates, its class (unit of the finest level and mixed-MAT flag), and the labelled samples without coordinates (kept in the phylogeny and delimitation, not here)
Rscript - "${LABELS}" "${META}" "${K_LEVELS}" "${K_FINE}" "${MAT_EXCLUDE}" "${CLASSES}" "${EXCLUDED}" <<'RSCRIPT'
a <- as.list(commandArgs(trailingOnly = TRUE))
names(a) <- c("labels", "meta", "k_levels", "k_fine", "mat_exclude", "out", "excluded")

k_levels <- strsplit(a$k_levels, " ")[[1]]

lab  <- read.delim(a$labels, colClasses = "character")
meta <- read.delim(a$meta, colClasses = "character")

### labelled samples without sampling coordinates
exc <- merge(lab[, c("sample", "clustered")], meta[, c("sample", "project", "site", "latitude", "longitude")], by = "sample")
exc <- exc[is.na(as.numeric(exc$latitude)) | is.na(as.numeric(exc$longitude)), c("sample", "project", "site", "clustered")]
exc$reason <- rep("no sampling coordinates", nrow(exc))

write.table(exc, a$excluded, sep = "\t", quote = FALSE, row.names = FALSE)

rep <- merge(lab[lab$clustered == "TRUE", c("sample", k_levels)],
             meta[, c("sample", "site", "mat_genome")], by = "sample")

## check if every representative is in the metadata
if (nrow(rep) != sum(lab$clustered == "TRUE")) stop("representatives missing from the metadata")

rep <- rep[!rep$sample %in% exc$sample, ]
rep$mixed <- rep$mat_genome %in% strsplit(a$mat_exclude, ",")[[1]]
rep$class <- paste(rep[[a$k_fine]], ifelse(rep$mixed, "mixed", "single"), sep = "_")
rep <- rep[order(rep$class, rep$sample), c("sample", "site", k_levels, "mixed", "class")]

write.table(rep, a$out, sep = "\t", quote = FALSE, row.names = FALSE)

cat("excluded (no coordinates):", nrow(exc), "labelled samples,", sum(exc$clustered == "TRUE"), "representatives\n")
if (nrow(exc) > 0) print(exc, row.names = FALSE)
cat("representatives: ", nrow(rep), "\n")
cat("mixed-MAT:       ", sum(rep$mixed), "\n\n")
print(table(rep$class))
RSCRIPT

## check if the classes and the exclusion list were written
if [[ ! -s "${CLASSES}" || ! -f "${EXCLUDED}" ]]; then
    echo "ERROR: no classes written"
    exit 1
fi


### =====================================================================
### step 2: pairwise differences and compared sites
### =====================================================================

echo ""
echo "=== Step 2: Pairwise differences and compared sites ==="

### differences and compared sites per locus and class pair (within a class each pair once) and per pair of representatives over all loci (a site counts for a pair when both carry A, C, G or T)
Rscript - "${ALN_DIR}" "${CLASSES}" "${CLASS_SUMS}" "${PAIR_DIFFS}" "${PAIR_SITES}" "${LOCI}" <<'RSCRIPT'
suppressMessages(library(ape))

a <- as.list(commandArgs(trailingOnly = TRUE))
names(a) <- c("aln_dir", "classes", "out_sums", "out_diffs", "out_sites", "out_loci")

cls <- read.delim(a$classes, colClasses = "character")
n   <- nrow(cls)

### membership of the representatives in the classes
classes <- sort(unique(cls$class))
memb    <- sapply(classes, function(k) as.numeric(cls$class == k))

### class pairs (upper triangle with the diagonal)
up <- which(upper.tri(diag(length(classes)), diag = TRUE), arr.ind = TRUE)

files <- list.files(a$aln_dir, pattern = "_final_alignment\\.fasta$", full.names = TRUE)

diffs_all <- matrix(0, n, n, dimnames = list(cls$sample, cls$sample))
sites_all <- diffs_all
loci      <- vector("list", length(files))

con <- gzfile(a$out_sums, "w")
writeLines("locus\tclass_a\tclass_b\tdifferences\tsites", con)

t0 <- Sys.time()
for (i in seq_along(files)) {
  x <- read.FASTA(files[i])
  if (!all(cls$sample %in% names(x))) stop("representatives missing from ", basename(files[i]))
  ch <- as.character(as.matrix(x[cls$sample]))

  ### indicator per base, and called sites (A, C, G, T)
  base  <- lapply(c("a", "c", "g", "t"), function(b) (ch == b) * 1)
  valid <- Reduce(`+`, base)

  ### per pair: compared sites and differences
  s <- tcrossprod(valid)
  d <- s - Reduce(`+`, lapply(base, tcrossprod))
  diag(s) <- 0

  diffs_all <- diffs_all + d
  sites_all <- sites_all + s

  ### per class pair (within a class each pair counted once)
  cd <- crossprod(memb, d %*% memb)
  cs <- crossprod(memb, s %*% memb)
  diag(cd) <- diag(cd) / 2
  diag(cs) <- diag(cs) / 2

  locus <- sub("_final_alignment\\.fasta$", "", basename(files[i]))
  writeLines(paste(locus, classes[up[, 1]], classes[up[, 2]], cd[up], cs[up], sep = "\t"), con)

  present   <- sapply(base, function(b) colSums(b) > 0)
  loci[[i]] <- data.frame(locus = locus, length = ncol(ch), called = round(mean(valid), 4),
                          variable = sum(rowSums(present) > 1))

  if (i %% 500 == 0) message(i, " loci, ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")
}
close(con)

loci <- do.call(rbind, loci)
write.table(loci, a$out_loci, sep = "\t", quote = FALSE, row.names = FALSE)

write.table(diffs_all, gzfile(a$out_diffs), sep = "\t", quote = FALSE, col.names = NA)
write.table(sites_all, gzfile(a$out_sites), sep = "\t", quote = FALSE, col.names = NA)

## check if every pair of representatives was compared at some site
n_zero <- sum(sites_all[upper.tri(sites_all)] == 0)
if (n_zero > 0) warning(n_zero, " pairs of representatives without a compared site")

cat("loci:                        ", nrow(loci), "\n")
cat("sites per locus (median):    ", median(loci$length), "\n")
cat("variable sites (median):     ", median(loci$variable), "\n")
cat("called share (median):       ", median(loci$called), "\n")
cat("compared sites per pair (median):", median(sites_all[upper.tri(sites_all)]), "\n")
cat("p-distance per pair (median):    ",
    signif(median(diffs_all[upper.tri(diffs_all)] / sites_all[upper.tri(sites_all)]), 3), "\n")
RSCRIPT

## check if every alignment was processed
N_LOCI=$(( $(wc -l < "${LOCI}") - 1 ))
echo "Loci processed: ${N_LOCI} of ${N_GENES}"
[[ "${N_LOCI}" -eq "${N_GENES}" ]] || { echo "ERROR: not every alignment was processed"; exit 1; }

### output files (copy classes.tsv, class_pair_sums.tsv.gz, pair_differences.tsv.gz, pair_sites.tsv.gz, loci.tsv and samples_excluded_no_coordinates.txt to local/04_diversity_niche_selection/data/010_data_sequence_distances for local/04_diversity_niche_selection/011, 012, 030, 031 and 032)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  samples_excluded_no_coordinates.txt   labelled samples without sampling coordinates (not used here)"
echo "  classes.tsv                 representatives, site, units, mixed-MAT flag, class"
echo "  class_pair_sums.tsv.gz      per locus and class pair: differences and compared sites"
echo "  pair_differences.tsv.gz     differences per pair of representatives, all loci"
echo "  pair_sites.tsv.gz           compared sites per pair of representatives, all loci"
echo "  loci.tsv                    per locus: length, called share, variable sites"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Sequence distances completed in ${RUNTIME} seconds ($((RUNTIME / 3600))h $(((RUNTIME % 3600) / 60))m)."
