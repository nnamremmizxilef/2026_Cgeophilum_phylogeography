#!/bin/bash
#SBATCH --job-name=species_tree_decloned    ### job name
#SBATCH --cpus-per-task=32                  ### cores
#SBATCH --mem-per-cpu=20G                   ### memory per core
#SBATCH --time=21-00:00:00                  ### walltime
#SBATCH --output /dev/null                  ### stdout (logged by tee)
#SBATCH --error /dev/null                   ### stderr (logged by tee)

### no core dumps
ulimit -c 0

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/02_structure_phylogeny/results/052_results_species_tree_decloned"
LOGS="${ROOT}/02_structure_phylogeny/logs/052_logs_species_tree_decloned"
### from hpc/02_structure_phylogeny/050_species_tree.sh (concatenation)
CONCATENATED="${ROOT}/02_structure_phylogeny/results/050_results_species_tree/cenococcum_orthologs_concatenated.fasta"
### from hpc/02_structure_phylogeny/050_species_tree.sh (best tree)
FULL_TREE="${ROOT}/02_structure_phylogeny/results/050_results_species_tree/species_tree.treefile"
### from local/03_delimitation/001_lineage_labels.R (copy local/03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv and lineage_labels.tsv.md5 to ${ROOT}/03_delimitation/supplementary)
LABELS_DIR="${ROOT}/03_delimitation/supplementary"
### frozen labels (one representative per MLG: clustered = TRUE)
LABELS="${LABELS_DIR}/lineage_labels.tsv"

KEEP="${RESULTS}/tips_kept.txt"
DECLONED="${RESULTS}/cenococcum_orthologs_concatenated_decloned.fasta"
START_TREE="${RESULTS}/start_tree.nwk"
PREFIX="${RESULTS}/species_tree_decloned"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-32}"
### outgroup (as in 050) and tips kept besides the MLG representatives (assemblies and the reference)
OUTGROUP="PsefloM402_1,PsefloM405_1,PsefloM406_1,M4-01,M4-02,M4-03,M4-04,M4-05,M4-06"
KEEP_PATTERN="^(Cenge|Pseflo|M4-|consensus$)"
### IQ-TREE memory (RAM 050 reports as required, scaled by the share of tips kept, with a margin)
IQ_MEM="600G"
### C locale (fixed sort order)
export LC_ALL=C

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_species_tree_decloned.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_seqkit

START_TIME=$(date +%s)

### run info
echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"

## check if the labels are the frozen version
( cd "${LABELS_DIR}" && md5sum -c lineage_labels.tsv.md5 ) \
  || { echo "ERROR: lineage_labels.tsv does not match its md5"; exit 1; }

## check if the inputs are there
for F in "${CONCATENATED}" "${LABELS}"; do
    if [[ ! -s "${F}" ]]; then
        echo "ERROR: input not found: ${F}"
        exit 1
    fi
done


### =====================================================================
### step 1: tips kept
### =====================================================================

echo ""
echo "=== Step 1: Tips kept ==="

seqkit seq -n -i "${CONCATENATED}" | sort > "${RESULTS}/tips_all.txt"

### MLG representatives (column clustered of the labels)
awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "clustered") c = i; next} $c == "TRUE" {print $1}' \
    "${LABELS}" | sort > "${RESULTS}/representatives.txt"

## check if every representative is in the concatenation
N_MISS=$(comm -23 "${RESULTS}/representatives.txt" "${RESULTS}/tips_all.txt" | wc -l)
if [[ "${N_MISS}" -gt 0 ]]; then
    echo "ERROR: ${N_MISS} representatives are not in the concatenation"
    exit 1
fi

{ cat "${RESULTS}/representatives.txt"; grep -E "${KEEP_PATTERN}" "${RESULTS}/tips_all.txt"; } \
    | sort -u > "${KEEP}"

echo "Tips in the concatenation: $(wc -l < "${RESULTS}/tips_all.txt")"
echo "MLG representatives:       $(wc -l < "${RESULTS}/representatives.txt")"
echo "Tips kept:                 $(wc -l < "${KEEP}")"

## check if every outgroup tip is kept
for T in ${OUTGROUP//,/ }; do
    grep -qx "${T}" "${KEEP}" || { echo "ERROR: outgroup tip ${T} not kept"; exit 1; }
done


### =====================================================================
### step 2: decloned alignment
### =====================================================================

echo ""
echo "=== Step 2: Decloned alignment ==="

### reused when newer than the concatenation and the tip list (resubmission goes straight to the tree)
if [[ -s "${DECLONED}" && "${DECLONED}" -nt "${CONCATENATED}" && "${DECLONED}" -nt "${KEEP}" ]]; then
    echo "Decloned alignment found, reused."
else
    seqkit grep -f "${KEEP}" "${CONCATENATED}" > "${DECLONED}"
fi

## check if the alignment holds every kept tip
N_SEQ=$(grep -c '^>' "${DECLONED}")
if [[ "${N_SEQ}" -ne "$(wc -l < "${KEEP}")" ]]; then
    echo "ERROR: ${N_SEQ} sequences in the decloned alignment, $(wc -l < "${KEEP}") tips kept"
    exit 1
fi
echo "Sequences: ${N_SEQ}"


### =====================================================================
### step 3: starting tree
### =====================================================================

echo ""
echo "=== Step 3: Starting tree ==="

### best tree of 050 pruned to the kept tips (saves the initial search)
conda deactivate
conda activate cg_bioinfo_v1_R

Rscript -e '
  a  <- commandArgs(trailingOnly = TRUE)
  tr <- ape::read.tree(a[1])
  tr$node.label <- NULL
  keep <- readLines(a[2])
  tr <- ape::drop.tip(tr, setdiff(tr$tip.label, keep))
  if (!setequal(tr$tip.label, keep)) stop("kept tips missing from the full tree")
  ape::write.tree(tr, a[3])
  cat("starting tree tips:", ape::Ntip(tr), "\n")' "${FULL_TREE}" "${KEEP}" "${START_TREE}"

START="${START_TREE}"
echo "Starting tree: pruned best tree of 050"


### =====================================================================
### step 4: IQ-TREE
### =====================================================================

echo ""
echo "=== Step 4: IQ-TREE ==="

conda deactivate
conda activate cg_bioinfo_v1_iqtree
echo "Running in environment: $CONDA_DEFAULT_ENV"
echo "iqtree version: $(iqtree2 --version | head -1)"

if [[ -s "${PREFIX}.ckp.gz" ]]; then
    echo "Checkpoint found, the run continues from it."
fi

### GTR+G as in 050, 1000 ultrafast bootstraps optimised by NNI (bnni), rooted on the Pseudocenococcum assemblies
iqtree2 -s "${DECLONED}" \
    -B 1000 \
    --bnni \
    -nm 1000 \
    -t "${START}" \
    -T "${THREADS}" \
    -mem "${IQ_MEM}" \
    -m GTR+G \
    -o "${OUTGROUP}" \
    --prefix "${PREFIX}"

## check if IQ-TREE wrote a tree
if [[ ! -s "${PREFIX}.treefile" ]]; then
    echo "ERROR: IQ-TREE wrote no tree"
    exit 1
fi

echo "Tips in the tree: $(tr ',' '\n' < "${PREFIX}.treefile" | grep -c '[A-Za-z]')"

### output files (copy species_tree_decloned.treefile to local/02_structure_phylogeny/data/052_data_species_tree_decloned for local/02_structure_phylogeny/052_species_tree_decloned.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  representatives.txt                MLG representatives (clustered = TRUE)"
echo "  tips_kept.txt                      representatives, assemblies, reference and outgroup"
echo "  $(basename "${DECLONED}")"
echo "  start_tree.nwk                     best tree of 050 pruned to the kept tips"
echo "  species_tree_decloned.treefile     ML tree with UFBoot support"
echo "  species_tree_decloned.contree      consensus tree"
echo "  species_tree_decloned.iqtree       report, model fit and convergence"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Decloned species tree completed in ${RUNTIME} seconds ($((RUNTIME / 3600))h $(((RUNTIME % 3600) / 60))m)."
