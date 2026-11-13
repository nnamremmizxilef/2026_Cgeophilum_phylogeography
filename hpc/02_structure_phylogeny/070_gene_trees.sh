#!/bin/bash
#SBATCH --job-name=gene_trees    ### job name
#SBATCH --array=1-5616%100       ### array (one task per orthogroup in og_list.txt)
#SBATCH --cpus-per-task=1        ### cores
#SBATCH --mem-per-cpu=8G         ### memory per core
#SBATCH --time=08:00:00          ### walltime
#SBATCH --output /dev/null       ### stdout (logged by tee)
#SBATCH --error /dev/null        ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
VARCALL="${ROOT}/01_variant_calling"
RESULTS="${ROOT}/02_structure_phylogeny/results/070_results_gene_trees"
LOGS="${ROOT}/02_structure_phylogeny/logs/070_logs_gene_trees/logs_per_og"
### from hpc/02_structure_phylogeny/042_rename_headers.sh
ALN_DIR="${ROOT}/02_structure_phylogeny/results/042_results_rename_headers/per_og_renamed"
### from hpc/01_variant_calling/050_cds_extraction.sh
OG_LIST="${VARCALL}/results/050_results_sco_consensus/og_list.txt"

TREE_DIR="${RESULTS}/per_og_tree"
WORK_DIR="${RESULTS}/per_og_work"
EMPTY_DIR="${RESULTS}/per_og_empty_tips"

### parameters
### array task (first argument outside SLURM)
TASK_ID="${SLURM_ARRAY_TASK_ID:-$1}"

### create results folders
mkdir -p "${TREE_DIR}" "${WORK_DIR}" "${EMPTY_DIR}" "${LOGS}"

### orthogroup of this array task
OG=$(sed -n "${TASK_ID}p" "${OG_LIST}")

## check if the array task is within the orthogroup list
if [[ -z "${OG}" ]]; then
    echo "ERROR: array task ${TASK_ID} is beyond the end of ${OG_LIST}" \
      > "${LOGS}/${SLURM_JOB_ID:-local}_task${TASK_ID}.log"
    exit 1
fi

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_${OG}.log"
exec > >(tee "${LOGFILE}") 2>&1

### run info
echo "Date: $(date)"
echo "Orthogroup: ${OG}"

ALN="${ALN_DIR}/${OG}_final_alignment.fasta"

## check if the alignment is there
if [[ ! -s "${ALN}" ]]; then
    echo "ERROR: alignment not found: ${ALN}"
    exit 1
fi

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_seqkit


### =====================================================================
### step 1: empty tips
### =====================================================================

### tips without a called base (no information; rejected by iqtree2 as "contains only gaps or missing data")
EMPTY="${EMPTY_DIR}/${OG}_empty_tips.txt"

seqkit fx2tab -n -i -B ACGTacgt "${ALN}" \
    | awk '$2 <= 0.001 {print $1}' > "${EMPTY}"

echo "Empty tips: $(wc -l < "${EMPTY}")"

### tips dropped from this gene (only the empty ones; clonemates kept, same taxa as the concatenation)
DROP="${WORK_DIR}/${OG}_dropped.txt"
WORK="${WORK_DIR}/${OG}_work.fasta"

sort -u "${EMPTY}" > "${DROP}"

### genes that lose no tip copied through untouched (no empty pattern file for seqkit grep)
if [[ -s "${DROP}" ]]; then
    seqkit grep -v -n -f "${DROP}" "${ALN}" > "${WORK}"
else
    cp "${ALN}" "${WORK}"
fi

echo "Tips: $(grep -c '^>' "${WORK}")"

### genes with fewer than four tips recorded and skipped
if [[ "$(grep -c '^>' "${WORK}")" -lt 4 ]]; then
    echo "SKIP: fewer than four tips left for ${OG}"
    rm -f "${DROP}" "${WORK}"
    exit 0
fi


### =====================================================================
### step 2: gene tree
### =====================================================================

conda deactivate
conda activate cg_bioinfo_v1_iqtree
echo "iqtree version: $(iqtree2 --version | head -1)"

### GTR as in the species tree, rate model per gene by BIC (no bootstrap)
iqtree2 -s "${WORK}" \
    -m MFP -mset GTR -mrate E,G,I,I+G \
    -T 1 \
    --prefix "${TREE_DIR}/${OG}"

## check if IQ-TREE wrote a gene tree
if [[ ! -s "${TREE_DIR}/${OG}.treefile" ]]; then
    echo "ERROR: no gene tree written for ${OG}"
    exit 1
fi

echo "Model: $(grep -m1 'Best-fit model' "${TREE_DIR}/${OG}.iqtree" | sed 's/.*: //')"

### remove intermediate files
rm -f "${DROP}" "${WORK}" "${WORK}.ckp.gz" \
      "${TREE_DIR}/${OG}.uniqueseq.phy" "${TREE_DIR}/${OG}.ckp.gz"

echo "${OG}: gene tree complete."

### output files (copy per_og_tree (*.treefile, *.iqtree) and per_og_empty_tips to local/02_structure_phylogeny/data/070_data_gene_trees for local/02_structure_phylogeny/070_gene_trees.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  per_og_tree/${OG}.treefile"
echo "  per_og_tree/${OG}.iqtree"
echo "  per_og_empty_tips/${OG}_empty_tips.txt"
