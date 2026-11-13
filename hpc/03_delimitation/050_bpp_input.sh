#!/bin/bash
#SBATCH --job-name=bpp_input    ### job name
#SBATCH --cpus-per-task=1       ### cores
#SBATCH --mem-per-cpu=8G        ### memory per core
#SBATCH --time=06:00:00         ### walltime
#SBATCH --output /dev/null      ### stdout (logged by tee)
#SBATCH --error /dev/null       ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
PHYLO="${ROOT}/02_structure_phylogeny"
RESULTS="${ROOT}/03_delimitation/results/050_results_bpp_input"
LOGS="${ROOT}/03_delimitation/logs/050_logs_bpp_input"
### from local/03_delimitation/001_lineage_labels.R (copy local/03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv and lineage_labels.tsv.md5 to ${ROOT}/03_delimitation/supplementary)
LABELS_DIR="${ROOT}/03_delimitation/supplementary"
### frozen labels (one representative per MLG: clustered = TRUE)
LABELS="${LABELS_DIR}/lineage_labels.tsv"
### from local/01_variant_calling/103_mating_type.R (copy local/01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv to ${ROOT}/03_delimitation/supplementary)
META="${LABELS_DIR}/sample_metadata_ALL_FIXED.tsv"
### from hpc/02_structure_phylogeny/042_rename_headers.sh (per-gene alignments)
ALN_DIR="${PHYLO}/results/042_results_rename_headers/per_og_renamed"
### from hpc/03_delimitation/010_aster_lineages.sh (unit trees)
TREE_DIR="${ROOT}/03_delimitation/results/010_results_aster_lineages"

INDIVIDUALS="${RESULTS}/individuals.tsv"
LOCI="${RESULTS}/loci.tsv"
DIVERSITY="${RESULTS}/diversity.tsv"
PRIORS="${RESULTS}/priors.tsv"
RUNS="${RESULTS}/runs.tsv"
RUN_DIR="${RESULTS}/runs"

### parameters
### levels analysed (A00 at every level, A11 from the finest)
K_LEVELS="k4 k5 k6"
K_FINE=k6
### mating type calls of possible mixed cultures (left out)
MAT_EXCLUDE="MAT12,ambiguous"
### loci per set, disjoint random sets, most missing data allowed per sequence, seed (individuals per unit: the smallest unit of the finest level without mixed-MAT samples)
N_LOCI=1000
N_SETS=3
MAX_MISSING=0.5
SEED=42
### loci per set for A11 (the first of each set; joint delimitation and species tree are far slower than A00)
N_LOCI_A11=500
### inverse-gamma shape for theta and root tau (diffuse, mean from the data)
PRIOR_ALPHA=3
### theta estimated in A00 and A11 (integrated theta runs serial code and drifts in BPP 4.8.7)
THETA_MODE=e
### locus rates: gamma with mean 1 fixed and shape 2 for coding loci (BPP manual)
LOCUSRATE="1 0 0 2 iid"
### rjMCMC for A11 (algorithm 1, a = 2, m = 1) and uniform prior on rooted trees
DELIMITATION="1 1 2 1"
MODEL_PRIOR=1
### MCMC (500,000 iterations), checkpoint at every 100,000 and at the end, independent chains, BPP threads per run (written into the control files)
BURNIN=50000
SAMPFREQ=10
NSAMPLE=45000
CHECKPOINT="100000 100000"
N_CHAINS=2
THREADS=20
### C locale (fixed sort order)
export LC_ALL=C

### create results folders
mkdir -p "${RESULTS}" "${RUN_DIR}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_bpp_input.log"
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

## check if the unit trees are there
for K in ${K_LEVELS}; do
    if [[ ! -s "${TREE_DIR}/cluster_tree_${K}_rooted.tre" ]]; then
        echo "ERROR: unit tree not found: ${TREE_DIR}/cluster_tree_${K}_rooted.tre"
        exit 1
    fi
done

## check if the alignments are there
N_GENES=$(ls "${ALN_DIR}"/*_final_alignment.fasta | wc -l)
if [[ "${N_GENES}" -eq 0 ]]; then
    echo "ERROR: no alignments in ${ALN_DIR}"
    exit 1
fi


### =====================================================================
### step 1: individuals, Imap files and species trees
### =====================================================================

echo ""
echo "=== Step 1: Individuals, Imap files and species trees ==="

### the same individuals at every level (a coarser unit holds the individuals of its parts), reused if present (A00 runs read these files)
if [[ -s "${INDIVIDUALS}" ]]; then
    echo "Individuals, Imap files and species trees found, reused."
else
Rscript - "${LABELS}" "${META}" "${TREE_DIR}" "${K_LEVELS}" "${K_FINE}" "${MAT_EXCLUDE}" \
    "${SEED}" "${RESULTS}" <<'RSCRIPT'
suppressMessages(library(ape))

a <- as.list(commandArgs(trailingOnly = TRUE))
names(a) <- c("labels", "meta", "tree_dir", "k_levels", "k_fine", "mat_exclude", "seed", "out_dir")

k_levels <- strsplit(a$k_levels, " ")[[1]]

lab  <- read.delim(a$labels, colClasses = "character")
meta <- read.delim(a$meta, colClasses = "character")

### one representative per MLG, mixed-MAT samples left out
rep  <- merge(lab[lab$clustered == "TRUE", c("sample", k_levels)],
              meta[, c("sample", "site", "mat_genome", "pct_scos_present")], by = "sample")
pool <- rep[!rep$mat_genome %in% strsplit(a$mat_exclude, ",")[[1]], ]

## check if every representative is in the metadata
if (nrow(rep) != sum(lab$clustered == "TRUE")) stop("representatives missing from the metadata")

### individuals per unit: the smallest unit of the finest level without mixed-MAT samples
n_ind <- min(table(pool[[a$k_fine]]))

### per unit: sites in random order, the most complete sample of a site first, one per site in turn
set.seed(as.integer(a$seed))
pick <- do.call(rbind, lapply(split(pool, pool[[a$k_fine]]), function(u) {
  u$site_rank <- match(u$site, sample(unique(u$site)))
  u <- u[order(u$site_rank, -as.numeric(u$pct_scos_present)), ]
  u$turn <- ave(u$site_rank, u$site_rank, FUN = seq_along)
  head(u[order(u$turn, u$site_rank), ], n_ind)
}))
pick <- pick[order(pick[[a$k_fine]], pick$sample), c("sample", "site", "pct_scos_present", k_levels)]

write.table(pick, file.path(a$out_dir, "individuals.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

### Imap (individual and unit) and species tree without the outgroup, per level
for (k in k_levels) {
  write.table(pick[, c("sample", k)], file.path(a$out_dir, sprintf("imap_%s.txt", k)),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  tr <- read.tree(file.path(a$tree_dir, sprintf("cluster_tree_%s_rooted.tre", k)))
  tr <- drop.tip(tr, "outgroup")
  tr$edge.length <- NULL
  tr$node.label  <- NULL

  ## check if the tree holds exactly the units of the level
  if (!setequal(tr$tip.label, unique(pick[[k]]))) stop("tree tips and units differ at ", k)

  write.tree(tr, file.path(a$out_dir, sprintf("tree_%s.nwk", k)))
  cat(k, write.tree(tr), "\n")
}

cat("\nrepresentatives:           ", nrow(rep), "\n")
cat("mixed-MAT left out:        ", nrow(rep) - nrow(pool), "\n")
cat("individuals per unit (", a$k_fine, "): ", n_ind, "\n", sep = "")
cat("individuals in total:      ", nrow(pick), "\n")
cat("sites per unit:\n")
print(tapply(pick$site, pick[[a$k_fine]], function(s) length(unique(s))))
RSCRIPT
fi

## check if individuals were selected
if [[ ! -s "${INDIVIDUALS}" ]]; then
    echo "ERROR: no individuals selected"
    exit 1
fi


### =====================================================================
### step 2: locus sets and priors
### =====================================================================

echo ""
echo "=== Step 2: Locus sets and priors ==="

### loci in random order with at most MAX_MISSING missing per sequence (all-missing columns dropped) until N_SETS x N_LOCI, prior means from the diversity within and between units of the finest level, reused if present (A00 runs read these files)
if [[ -s "${LOCI}" && -s "${PRIORS}" ]]; then
    echo "Locus sets and priors found, reused."
else
Rscript - "${ALN_DIR}" "${INDIVIDUALS}" "${K_FINE}" "${N_LOCI}" "${N_SETS}" "${MAX_MISSING}" \
    "${SEED}" "${PRIOR_ALPHA}" "${RESULTS}" <<'RSCRIPT'
suppressMessages(library(ape))

a <- as.list(commandArgs(trailingOnly = TRUE))
names(a) <- c("aln_dir", "individuals", "k_fine", "n_loci", "n_sets", "max_missing", "seed",
              "prior_alpha", "out_dir")

n_loci <- as.integer(a$n_loci)
n_sets <- as.integer(a$n_sets)
alpha  <- as.numeric(a$prior_alpha)

ind   <- read.delim(a$individuals, colClasses = "character")
unit  <- setNames(ind[[a$k_fine]], ind$sample)
units <- sort(unique(unit))

files <- list.files(a$aln_dir, pattern = "_final_alignment\\.fasta$", full.names = TRUE)
set.seed(as.integer(a$seed))
files <- sample(files)

con <- lapply(seq_len(n_sets), function(s) file(file.path(a$out_dir, sprintf("loci_set%d.txt", s)), "w"))

### mean pairwise p-distance within and between units, summed over loci
pair_sum <- matrix(0, length(units), length(units), dimnames = list(units, units))
loci     <- list()
n_read   <- 0

for (f in files) {
  if (length(loci) == n_loci * n_sets) break
  n_read <- n_read + 1

  x <- read.FASTA(f)
  if (!all(ind$sample %in% names(x))) stop("selected samples missing from ", basename(f))
  x <- as.matrix(x[ind$sample])

  ### missing share per sequence; skip the locus if any sequence lacks too much
  ch   <- as.character(x)
  miss <- ch %in% c("n", "-", "?")
  dim(miss) <- dim(ch)
  if (max(rowMeans(miss)) > as.numeric(a$max_missing)) next

  keep <- colSums(!miss) > 0
  x    <- x[, keep]
  ch   <- toupper(ch[, keep])

  s <- length(loci) %/% n_loci + 1
  cat("\n", nrow(ch), " ", ncol(ch), "\n\n", sep = "", file = con[[s]])
  writeLines(paste0(ind$sample, "^", ind$sample, "  ", apply(ch, 1, paste, collapse = "")), con[[s]])

  d <- dist.dna(x, model = "raw", pairwise.deletion = TRUE, as.matrix = TRUE)
  diag(d) <- NA
  pair_sum <- pair_sum + tapply(d, list(unit[row(d)], unit[col(d)]), mean, na.rm = TRUE)[units, units]

  loci[[length(loci) + 1]] <- data.frame(set = s, locus = sub("_final_alignment\\.fasta$", "", basename(f)),
                                         sites = ncol(ch), variable = sum(apply(ch, 2, function(v)
                                           length(setdiff(unique(v), c("N", "-", "?"))) > 1)),
                                         max_missing = round(max(rowMeans(miss)), 4))
}

invisible(lapply(con, close))

## check if enough loci passed
if (length(loci) < n_loci * n_sets) stop("only ", length(loci), " loci passed the missing-data filter")

loci <- do.call(rbind, loci)
write.table(loci, file.path(a$out_dir, "loci.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

### diversity per locus averaged over loci (within-unit theta and between-unit distance)
div <- pair_sum / nrow(loci)
div_long <- data.frame(unit_a = rep(units, times = length(units)), unit_b = rep(units, each = length(units)),
                       p_distance = signif(as.vector(div), 4))
write.table(div_long[div_long$unit_a <= div_long$unit_b, ], file.path(a$out_dir, "diversity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

### prior means: mean theta over units; root tau from the most distant pair, (d - theta) / 2
theta_mean <- mean(diag(div))
tau_mean   <- (max(div[upper.tri(div)]) - theta_mean) / 2
priors <- data.frame(parameter = c("theta", "tau_root"), alpha = alpha,
                     beta = signif(c(theta_mean, tau_mean) * (alpha - 1), 2))
priors$mean <- priors$beta / (priors$alpha - 1)
write.table(priors, file.path(a$out_dir, "priors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("alignments read:           ", n_read, "\n")
cat("loci kept:                 ", nrow(loci), "(", n_sets, "sets of", n_loci, ")\n")
cat("sites per locus (median):  ", median(loci$sites), "\n")
cat("variable sites (median):   ", median(loci$variable), "\n\n")
cat("within-unit p-distance (theta):\n")
print(signif(diag(div), 3))
cat("\nbetween-unit p-distance:\n")
print(signif(div, 3))
cat("\npriors (inverse gamma):\n")
print(priors)
RSCRIPT
fi

## check if the locus sets and priors were written
if [[ ! -s "${PRIORS}" || ! -s "${LOCI}" ]]; then
    echo "ERROR: locus sets or priors not written"
    exit 1
fi

### A11 locus sets: the first N_LOCI_A11 loci of each set
for S in $(seq 1 "${N_SETS}"); do
    awk -v n="${N_LOCI_A11}" '/^[0-9]+ [0-9]+$/ {i++} i <= n {print}' \
        "${RESULTS}/loci_set${S}.txt" > "${RESULTS}/loci_set${S}_a11.txt"
done


### =====================================================================
### step 3: control files
### =====================================================================

echo ""
echo "=== Step 3: Control files ==="

### A00 for every level, locus set and chain; A11 from the finest level for every set and chain
Rscript - "${RESULTS}" "${K_LEVELS}" "${K_FINE}" "${N_SETS}" "${N_LOCI}" "${N_LOCI_A11}" "${N_CHAINS}" \
    "${SEED}" "${THETA_MODE}" "${LOCUSRATE}" "${DELIMITATION}" "${MODEL_PRIOR}" "${BURNIN}" \
    "${SAMPFREQ}" "${NSAMPLE}" "${CHECKPOINT}" "${THREADS}" <<'RSCRIPT'
a <- as.list(commandArgs(trailingOnly = TRUE))
names(a) <- c("dir", "k_levels", "k_fine", "n_sets", "n_loci", "n_loci_a11", "n_chains", "seed",
              "theta_mode", "locusrate", "delimitation", "model_prior", "burnin", "sampfreq", "nsample",
              "checkpoint", "threads")

k_levels <- strsplit(a$k_levels, " ")[[1]]
priors   <- read.delim(file.path(a$dir, "priors.tsv"))
pr       <- setNames(sprintf("invgamma %g %g", priors$alpha, priors$beta), priors$parameter)

runs <- rbind(expand.grid(chain = seq_len(as.integer(a$n_chains)), set = seq_len(as.integer(a$n_sets)),
                          k = k_levels, analysis = "A00", stringsAsFactors = FALSE),
              expand.grid(chain = seq_len(as.integer(a$n_chains)), set = seq_len(as.integer(a$n_sets)),
                          k = a$k_fine, analysis = "A11", stringsAsFactors = FALSE))
runs <- runs[, c("analysis", "k", "set", "chain")]
### A11 runs named by their locus count
runs$run  <- ifelse(runs$analysis == "A00",
                    sprintf("a00_%s_set%d_chain%d", runs$k, runs$set, runs$chain),
                    sprintf("a11_%s_n%s_set%d_chain%d", runs$k, a$n_loci_a11, runs$set, runs$chain))
runs$seed <- as.integer(a$seed) * 1000 + seq_len(nrow(runs))

for (i in seq_len(nrow(runs))) {
  r    <- runs[i, ]
  imap <- read.delim(file.path(a$dir, sprintf("imap_%s.txt", r$k)), header = FALSE, col.names = c("ind", "unit"))
  n    <- table(imap$unit)
  tree <- readLines(file.path(a$dir, sprintf("tree_%s.nwk", r$k)))
  a00  <- r$analysis == "A00"

  ctl <- c(
    sprintf("          seed = %d", r$seed),
    sprintf("       seqfile = %s", file.path(a$dir, sprintf(if (a00) "loci_set%d.txt" else "loci_set%d_a11.txt", r$set))),
    sprintf("      Imapfile = %s", file.path(a$dir, sprintf("imap_%s.txt", r$k))),
    sprintf("       jobname = %s", r$run),
    "",
    sprintf("speciesdelimitation = %s", if (a00) "0" else a$delimitation),
    sprintf("   speciestree = %s", if (a00) "0" else "1"),
    if (!a00) sprintf("speciesmodelprior = %s", a$model_prior),
    "",
    sprintf("  species&tree = %d  %s", length(n), paste(names(n), collapse = " ")),
    sprintf("                    %s", paste(as.integer(n), collapse = " ")),
    sprintf("                    %s", tree),
    sprintf("         phase = %s", paste(rep(0, length(n)), collapse = " ")),
    "",
    "       usedata = 1",
    sprintf("         nloci = %s", if (a00) a$n_loci else a$n_loci_a11),
    "     cleandata = 0",
    "",
    sprintf("    thetaprior = %s %s", pr[["theta"]], a$theta_mode),
    sprintf("      tauprior = %s", pr[["tau_root"]]),
    sprintf("     locusrate = %s", a$locusrate),
    "",
    "      finetune = 1",
    "         print = 1 0 0 0 0",
    sprintf("        burnin = %s", a$burnin),
    sprintf("      sampfreq = %s", a$sampfreq),
    sprintf("       nsample = %s", a$nsample),
    sprintf("    checkpoint = %s", a$checkpoint),
    sprintf("       threads = %s", a$threads))

  dir.create(file.path(a$dir, "runs", r$run), showWarnings = FALSE)
  writeLines(ctl, file.path(a$dir, "runs", r$run, paste0(r$run, ".ctl")))
}

write.table(runs, file.path(a$dir, "runs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("runs:", nrow(runs), "\n")
print(table(runs$analysis, runs$k))
for (i in c(1, nrow(runs))) {
  cat("\nexample control file (", runs$run[i], "):\n", sep = "")
  writeLines(readLines(file.path(a$dir, "runs", runs$run[i], paste0(runs$run[i], ".ctl"))))
}
RSCRIPT

N_RUNS=$(( $(wc -l < "${RUNS}") - 1 ))

## check if every run has a control file
N_CTL=$(awk -F'\t' -v d="${RUN_DIR}" 'NR > 1 {print d "/" $5 "/" $5 ".ctl"}' "${RUNS}" \
    | while read -r C; do [[ -s "${C}" ]] && echo "${C}"; done | wc -l)
echo ""
echo "Runs: ${N_RUNS}   control files: ${N_CTL}"
[[ "${N_RUNS}" -eq "${N_CTL}" ]] || { echo "ERROR: runs and control files differ"; exit 1; }

## check if every locus set holds N_LOCI loci
for S in $(seq 1 "${N_SETS}"); do
    N=$(grep -cE '^[0-9]+ [0-9]+$' "${RESULTS}/loci_set${S}.txt")
    echo "loci_set${S}.txt: ${N} loci"
    [[ "${N}" == "${N_LOCI}" ]] || { echo "ERROR: loci_set${S}.txt does not hold ${N_LOCI} loci"; exit 1; }

    ## check if the A11 set holds N_LOCI_A11 loci
    N=$(grep -cE '^[0-9]+ [0-9]+$' "${RESULTS}/loci_set${S}_a11.txt")
    echo "loci_set${S}_a11.txt: ${N} loci"
    [[ "${N}" == "${N_LOCI_A11}" ]] || { echo "ERROR: loci_set${S}_a11.txt does not hold ${N_LOCI_A11} loci"; exit 1; }
done

### output files (copy runs.tsv, individuals.tsv, priors.tsv, loci.tsv and tree_k*.nwk to local/03_delimitation/data/050_data_bpp_input for local/03_delimitation/051_bpp.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  individuals.tsv          selected individuals, site and unit at every level"
echo "  imap_kX.txt              individual and unit per level"
echo "  tree_kX.nwk              unit tree of the level without the outgroup (from 010)"
echo "  loci_setS.txt            ${N_SETS} disjoint sets of ${N_LOCI} loci in BPP format"
echo "  loci_setS_a11.txt        the first ${N_LOCI_A11} loci of each set (A11)"
echo "  loci.tsv                 set, locus, sites, variable sites, most missing per locus"
echo "  diversity.tsv            mean pairwise p-distance within and between units"
echo "  priors.tsv               inverse-gamma priors for theta and root tau"
echo "  runs.tsv                 one row per BPP run (input to 051_bpp.sh)"
echo "  runs/<run>/<run>.ctl     control file per run"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "BPP input completed in ${RUNTIME} seconds ($((RUNTIME / 3600))h $(((RUNTIME % 3600) / 60))m)."
