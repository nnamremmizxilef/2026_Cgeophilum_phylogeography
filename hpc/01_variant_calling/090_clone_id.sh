#!/bin/bash
#SBATCH --job-name=clone_id         ### job name
#SBATCH --cpus-per-task=4           ### cores
#SBATCH --mem-per-cpu=8G            ### memory per core
#SBATCH --time=04:00:00             ### walltime
#SBATCH --output /dev/null          ### stdout (logged by tee)
#SBATCH --error /dev/null           ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/090_results_clone_id"
LOGS="${ROOT}/01_variant_calling/logs/090_logs_clone_id"
### from hpc/01_variant_calling/080_filtering.sh
FILTERED_DIR="${ROOT}/01_variant_calling/results/080_results_filtering"
INPUT_VCF="${FILTERED_DIR}/consensus_ingroup_noMAF.vcf.gz"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_clone_id.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_plink

START_TIME=$(date +%s)

### run info
echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"
echo "plink version: $(plink --version 2>&1 | head -1)"
echo ""

echo "Input VCF: ${INPUT_VCF}"


### =====================================================================
### step 1: pairwise counts and per-sample missingness
### =====================================================================

echo ""
echo "=== Step 1: plink ==="

### --allow-extra-chr (orthologue contigs), --double-id (sample name intact), --set-missing-var-ids (no variant IDs)
### haploid calls read as homozygous (IBS1 = 0, IBS0 = differing sites, IBS0 + IBS2 = sites called in both)
plink --vcf "${INPUT_VCF}" \
    --allow-extra-chr \
    --double-id \
    --set-missing-var-ids @:# \
    --genome full \
    --missing \
    --out "${RESULTS}/plink" \
    --threads "${THREADS}"

echo "Done."

N_PAIRS=$(( $(wc -l < "${RESULTS}/plink.genome") - 1 ))
N_SAMPLES=$(( $(wc -l < "${RESULTS}/plink.imiss") - 1 ))

echo ""
echo "Samples: ${N_SAMPLES}"
echo "Pairs:   ${N_PAIRS} (expected $((N_SAMPLES * (N_SAMPLES - 1) / 2)))"

### output files (copy plink.genome and plink.imiss to local/01_variant_calling/data/090_data_clone_id for local/01_variant_calling/090_clone_id.R and 091_genet_extent.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  plink.genome"
echo "  plink.imiss"
echo ""

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Completed in ${RUNTIME} seconds ($((RUNTIME / 60)) min)."
