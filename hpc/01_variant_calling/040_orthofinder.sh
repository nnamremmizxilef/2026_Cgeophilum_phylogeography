#!/bin/bash
#SBATCH --job-name=orthofinder      ### job name
#SBATCH --cpus-per-task=16          ### cores
#SBATCH --mem-per-cpu=8G            ### memory per core
#SBATCH --time=2-00:00:00           ### walltime
#SBATCH --output /dev/null          ### stdout (logged by tee)
#SBATCH --error /dev/null           ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/040_results_orthofinder"
LOGS="${ROOT}/01_variant_calling/logs/040_logs_orthofinder"
### external: predicted proteomes of the reference genomes, one FASTA per genome (copy manually into ${ROOT}/01_variant_calling/data/references/protein_sequences)
PROTEOMES="${ROOT}/01_variant_calling/data/references/protein_sequences"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-16}"

### create results folders (OrthoFinder creates RESULTS)
mkdir -p "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_orthofinder.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_orthofinder

START_TIME=$(date +%s)

echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"
echo "Number of proteomes: $(ls -1 ${PROTEOMES}/*.fasta | wc -l)"


### =====================================================================
### step 1: OrthoFinder
### =====================================================================

echo "Starting OrthoFinder analysis..."
orthofinder \
    -f "${PROTEOMES}" \
    -o "${RESULTS}" \
    -t "${THREADS}" \
    -a "${THREADS}" \
    -M msa \
    -S diamond

echo "OrthoFinder analysis complete!"


### output files
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  Results_<run date>/"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo "OrthoFinder completed in ${RUNTIME} seconds."
