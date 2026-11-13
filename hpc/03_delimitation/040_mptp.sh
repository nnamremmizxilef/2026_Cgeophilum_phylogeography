#!/bin/bash
#SBATCH --job-name=mptp       ### job name
#SBATCH --cpus-per-task=1     ### cores
#SBATCH --mem-per-cpu=8G      ### memory per core
#SBATCH --time=1-00:00:00     ### walltime
#SBATCH --output /dev/null    ### stdout (logged by tee)
#SBATCH --error /dev/null     ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
PHYLO="${ROOT}/02_structure_phylogeny"
RESULTS="${ROOT}/03_delimitation/results/040_results_mptp"
LOGS="${ROOT}/03_delimitation/logs/040_logs_mptp"
### from local/03_delimitation/001_lineage_labels.R (copy local/03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv and lineage_labels.tsv.md5 to ${ROOT}/03_delimitation/supplementary)
LABELS_DIR="${ROOT}/03_delimitation/supplementary"
### frozen labels (one representative per MLG: clustered = TRUE)
LABELS="${LABELS_DIR}/lineage_labels.tsv"
### from hpc/02_structure_phylogeny/050_species_tree.sh (concatenation tree)
CONCAT_TREE="${PHYLO}/results/050_results_species_tree/species_tree.treefile"
### from hpc/02_structure_phylogeny/050_species_tree.sh (IQ-TREE log, read for the alignment length)
CONCAT_LOG="${PHYLO}/results/050_results_species_tree/species_tree.log"

TREE="${RESULTS}/species_tree_representatives.nwk"
OUTGROUP_FILE="${RESULTS}/outgroup.txt"
SENSITIVITY="${RESULTS}/mptp_ml_sensitivity.tsv"

### parameters
### mPTP binary
MPTP=mptp
### outgroup (samples by prefix, assemblies by name)
OUTGROUP_PREFIX="M4-"
OUTGROUP_GENOMES="PsefloM402_1,PsefloM405_1,PsefloM406_1"
### minimum branch length: one substitution over the alignment times these factors (the first the main run, the others the sensitivity check)
MINBR_FACTORS="1 10 100"
### MCMC: iterations, sampling interval, burn-in, independent runs from random starts, seed
MCMC_N=100000000
MCMC_SAMPLE=10000
MCMC_BURNIN=10000000
MCMC_RUNS=4
SEED=42
export LC_ALL=C

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_mptp.log"
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
for F in "${LABELS}" "${CONCAT_TREE}" "${CONCAT_LOG}"; do
    if [[ ! -s "${F}" ]]; then
        echo "ERROR: input not found: ${F}"
        exit 1
    fi
done


### =====================================================================
### step 1: representatives tree and minimum branch length
### =====================================================================

echo ""
echo "=== Step 1: Tree and minimum branch length ==="

### alignment length from the IQ-TREE log
N_COLS=$(grep -m 1 '^Alignment has' "${CONCAT_LOG}" | sed 's/.* with \([0-9]*\) columns.*/\1/')

## check if the alignment length was read
if ! [[ "${N_COLS}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: alignment length not found in ${CONCAT_LOG}"
    exit 1
fi
echo "Alignment columns: ${N_COLS}"

### one representative per MLG plus the outgroup (clonemates, Cenococcum assemblies and reference dropped)
Rscript - "${CONCAT_TREE}" "${LABELS}" "${OUTGROUP_PREFIX}" "${OUTGROUP_GENOMES}" "${TREE}" \
    "${OUTGROUP_FILE}" <<'RSCRIPT'
suppressMessages(library(ape))

a <- as.list(commandArgs(trailingOnly = TRUE))
names(a) <- c("tree", "labels", "og_prefix", "og_genomes", "out_tree", "out_og")

tr  <- read.tree(a$tree)
lab <- read.delim(a$labels, colClasses = "character")

reps <- lab$sample[lab$clustered == "TRUE"]
og   <- tr$tip.label[startsWith(tr$tip.label, a$og_prefix) | tr$tip.label %in% strsplit(a$og_genomes, ",")[[1]]]

## check if every representative is in the tree
miss <- setdiff(reps, tr$tip.label)
if (length(miss) > 0) stop(length(miss), " representatives not in the tree: ", paste(head(miss), collapse = ", "))

sub <- keep.tip(tr, c(reps, og))
sub$node.label <- NULL

## check if the outgroup is one clade (unrooted bipartition)
if (!is.monophyletic(sub, og, reroot = TRUE)) stop("outgroup is not monophyletic")

write.tree(sub, a$out_tree)
writeLines(paste(og, collapse = ","), a$out_og)

cat("tips:", Ntip(sub), "(representatives", length(reps), "+ outgroup", length(og), ")\n")
cat("branch length: min", signif(min(sub$edge.length), 3), " median", signif(median(sub$edge.length), 3),
    " branches below 1e-6:", sum(sub$edge.length < 1e-6), "\n")
RSCRIPT

## check if the representatives tree was written
if [[ ! -s "${TREE}" || ! -s "${OUTGROUP_FILE}" ]]; then
    echo "ERROR: no representatives tree written"
    exit 1
fi

OUTGROUP=$(cat "${OUTGROUP_FILE}")

### mPTP
conda deactivate
conda activate cg_bioinfo_v1_mptp

## check if the mPTP binary is there
if ! command -v "${MPTP}" > /dev/null 2>&1; then
    echo "ERROR: ${MPTP} not found in cg_bioinfo_v1_mptp (mPTP 0.2.5 built from source)"
    exit 1
fi
echo "$("${MPTP}" --version 2>&1 | head -1)"


### =====================================================================
### step 2: maximum-likelihood delimitation and minimum branch length sensitivity
### =====================================================================

echo ""
echo "=== Step 2: ML delimitation ==="

printf "factor\tminbr\tspecies\tlrt_p\n" > "${SENSITIVITY}"

for F in ${MINBR_FACTORS}; do

    MINBR=$(awk -v f="${F}" -v n="${N_COLS}" 'BEGIN {printf "%.6e", f / n}')
    PREFIX="${RESULTS}/mptp_ml_minbr${F}"

    "${MPTP}" --ml --multi \
        --tree_file "${TREE}" \
        --outgroup "${OUTGROUP}" --outgroup_crop \
        --minbr "${MINBR}" \
        --output_file "${PREFIX}" > "${PREFIX}.log" 2>&1

    ## check if a delimitation was written
    if [[ ! -s "${PREFIX}.txt" ]]; then
        echo "ERROR: mPTP wrote no delimitation for minbr ${MINBR}, see ${PREFIX}.log"
        exit 1
    fi

    N_SP=$(grep -m 1 'Number of delimited species' "${PREFIX}.txt" | awk '{print $NF}')
    P_LRT=$(grep -m 1 'LRT computed p-value' "${PREFIX}.txt" | awk '{print $NF}')
    printf "%s\t%s\t%s\t%s\n" "${F}" "${MINBR}" "${N_SP}" "${P_LRT}" >> "${SENSITIVITY}"
done

awk -F'\t' '{printf "%-8s %-14s %-8s %s\n", $1, $2, $3, $4}' "${SENSITIVITY}"


### =====================================================================
### step 3: MCMC support at the main minimum branch length
### =====================================================================

echo ""
echo "=== Step 3: MCMC ==="

F_MAIN=$(echo "${MINBR_FACTORS}" | awk '{print $1}')
MINBR=$(awk -v f="${F_MAIN}" -v n="${N_COLS}" 'BEGIN {printf "%.6e", f / n}')
PREFIX="${RESULTS}/mptp_mcmc"

"${MPTP}" --mcmc "${MCMC_N}" --multi \
    --mcmc_sample "${MCMC_SAMPLE}" \
    --mcmc_burnin "${MCMC_BURNIN}" \
    --mcmc_runs "${MCMC_RUNS}" \
    --mcmc_startrandom \
    --mcmc_log \
    --seed "${SEED}" \
    --tree_file "${TREE}" \
    --outgroup "${OUTGROUP}" --outgroup_crop \
    --minbr "${MINBR}" \
    --output_file "${PREFIX}" > "${PREFIX}.log" 2>&1

## check if the combined support tree was written
if [[ ! -s "${PREFIX}.${SEED}.combined.tree" ]]; then
    echo "ERROR: no combined support tree, see ${PREFIX}.log"
    exit 1
fi

### credible interval of the species number per run, agreement of the runs
grep -E '^CCI|^HPD|ML average support|Average standard deviation' "${PREFIX}.log"

### output files (copy mptp_ml_minbr*.txt, mptp_ml_sensitivity.tsv, mptp_mcmc.*.stats and mptp_mcmc.42.combined.tree to local/03_delimitation/data/040_data_mptp for local/03_delimitation/040_mptp.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  species_tree_representatives.nwk   050 tree on the representatives and the outgroup"
echo "  mptp_ml_minbrF.txt                 ML delimitation per minimum branch length factor"
echo "  mptp_ml_sensitivity.tsv            species number per minimum branch length"
echo "  mptp_mcmc.${SEED}.combined.tree          support per node, runs combined"
echo "  mptp_mcmc.<seed>.stats             posterior of the species number per run"
echo "  mptp_mcmc.log                      MCMC log (credible intervals, agreement of the runs)"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "mPTP completed in ${RUNTIME} seconds ($((RUNTIME / 3600))h $(((RUNTIME % 3600) / 60))m)."
