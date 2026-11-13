#!/bin/bash
#SBATCH --job-name=species_tree     ### job name
#SBATCH --cpus-per-task=32          ### cores
#SBATCH --mem-per-cpu=32G           ### memory per core
#SBATCH --time=21-00:00:00          ### walltime
#SBATCH --output /dev/null          ### stdout (logged by tee)
#SBATCH --error /dev/null           ### stderr (logged by tee)

### no core dumps
ulimit -c 0

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/02_structure_phylogeny/results/050_results_species_tree"
LOGS="${ROOT}/02_structure_phylogeny/logs/050_logs_species_tree"
### from hpc/02_structure_phylogeny/042_rename_headers.sh
RENAMED_DIR="${ROOT}/02_structure_phylogeny/results/042_results_rename_headers/per_og_renamed"
CONCATENATED="${RESULTS}/cenococcum_orthologs_concatenated.fasta"
PARTITIONS="${RESULTS}/cenococcum_orthologs_partitions.txt"
PREFIX="${RESULTS}/species_tree"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-32}"
### outgroup tips (assemblies and short-read samples)
OUTGROUP="PsefloM402_1,PsefloM405_1,PsefloM406_1,M4-01,M4-02,M4-03,M4-04,M4-05,M4-06"

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_species_tree.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_amas

START_TIME=$(date +%s)

### run info
echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"

N_FILES=$(ls -1 "${RENAMED_DIR}"/*_final_alignment.fasta 2>/dev/null | wc -l)

## check if the renamed alignments are there
if [[ "${N_FILES}" -eq 0 ]]; then
    echo "ERROR: no renamed alignments in ${RENAMED_DIR}, run 042 first"
    exit 1
fi

echo "Alignments: ${N_FILES}"


### =====================================================================
### step 1: concatenate
### =====================================================================

echo ""
echo "=== Step 1: Concatenation ==="

### concatenation reused if present (a resubmission goes straight to the tree)
if [[ -s "${CONCATENATED}" ]]; then
    echo "Concatenation found, reused."
else

    cd "${RENAMED_DIR}"

    AMAS.py concat \
        -i *_final_alignment.fasta \
        -f fasta \
        -d dna \
        -t "${CONCATENATED}" \
        -p "${PARTITIONS}" \
        -u fasta

    cd "${ROOT}/02_structure_phylogeny"

fi

## check if the concatenation was written
if [[ ! -s "${CONCATENATED}" ]]; then
    echo "ERROR: concatenation produced no output"
    exit 1
fi

echo "Concatenated: $(grep -c '^>' "${CONCATENATED}") sequences"
echo "Columns: $(awk '/^>/ {if (seen) exit; seen = 1; next} {n += length($0)} END {print n}' "${CONCATENATED}")"


### =====================================================================
### step 2: IQ-TREE
### =====================================================================

echo ""
echo "=== Step 2: IQ-TREE ==="

conda deactivate
conda activate cg_bioinfo_v1_iqtree
echo "Running in environment: $CONDA_DEFAULT_ENV"
echo "iqtree version: $(iqtree2 --version | head -1)"

if [[ -s "${PREFIX}.ckp.gz" ]]; then
    echo "Checkpoint found, the run continues from it."
fi

### GTR with gamma rate heterogeneity, 1000 ultrafast bootstrap replicates, at most 1000 search iterations, parsimony starting tree, rooted on the outgroup
iqtree2 -s "${CONCATENATED}" \
    -B 1000 \
    -nm 1000 \
    -t PARS \
    -T "${THREADS}" \
    -mem 900G \
    -m GTR+G \
    -o "${OUTGROUP}" \
    --prefix "${PREFIX}"

## check if IQ-TREE wrote a tree
if [[ ! -s "${PREFIX}.treefile" ]]; then
    echo "ERROR: IQ-TREE wrote no tree"
    exit 1
fi

echo "Tips in the tree: $(tr ',' '\n' < "${PREFIX}.treefile" | grep -c '[A-Za-z]')"

### output files (copy species_tree.treefile to local/02_structure_phylogeny/data/050_data_species_tree for local/02_structure_phylogeny/050_species_tree.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  species_tree.treefile     ML tree with UFBoot support"
echo "  species_tree.contree      consensus tree"
echo "  species_tree.iqtree       report, model fit and convergence"
echo "  $(basename "${CONCATENATED}")"
echo "  $(basename "${PARTITIONS}")"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Species tree completed in ${RUNTIME} seconds ($((RUNTIME / 3600))h $(((RUNTIME % 3600) / 60))m)."
