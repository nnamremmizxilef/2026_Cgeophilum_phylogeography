#!/bin/bash
#SBATCH --job-name=consensus_align_trim    ### job name
#SBATCH --array=1-5807%100                 ### array (one task per orthogroup in og_list.txt)
#SBATCH --cpus-per-task=1                  ### cores
#SBATCH --mem-per-cpu=4G                   ### memory per core
#SBATCH --time=01:00:00                    ### walltime
#SBATCH --output /dev/null                 ### stdout (logged by tee)
#SBATCH --error /dev/null                  ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/054_results_sco_consensus_core"
LOGS="${ROOT}/01_variant_calling/logs/054_logs_sco_consensus_core/logs_align_trim"
### from hpc/01_variant_calling/054_cds_extraction_core.sh
INPUT_DIR="${RESULTS}/per_og_fasta"
### from hpc/01_variant_calling/054_cds_extraction_core.sh
OG_LIST="${RESULTS}/og_list.txt"
ALN_DIR="${RESULTS}/per_og_aln"
TRIM_DIR="${RESULTS}/per_og_trimmed"

### parameters
### array task (first argument outside SLURM)
TASK_ID="${SLURM_ARRAY_TASK_ID:-$1}"
### trimAl gap threshold (columns with more than half gaps removed)
GAP_THRESHOLD=0.5

### create results folders
mkdir -p "${ALN_DIR}" "${TRIM_DIR}" "${LOGS}"

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

echo "Date: $(date)"
echo "Orthogroup: ${OG}"

INPUT="${INPUT_DIR}/${OG}.fa"
ALN="${ALN_DIR}/${OG}.aln.fa"
TRIM="${TRIM_DIR}/${OG}.trimmed.fa"

## check if the orthogroup FASTA exists
if [[ ! -s "${INPUT}" ]]; then
    echo "ERROR: missing input ${INPUT}"
    exit 1
fi

echo "Sequences: $(grep -c '^>' "${INPUT}")"


### =====================================================================
### step 1: alignment
### =====================================================================

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_mafft

echo "MAFFT version: $(mafft --version 2>&1)"

mafft --auto --quiet "${INPUT}" > "${ALN}"

## check if MAFFT wrote an alignment
if [[ ! -s "${ALN}" ]]; then
    echo "ERROR: MAFFT produced no output for ${OG}"
    exit 1
fi


### =====================================================================
### step 2: trimming
### =====================================================================

conda deactivate
conda activate cg_bioinfo_v1_trimal

trimal -in "${ALN}" -out "${TRIM}" -gt "${GAP_THRESHOLD}"

## check if trimAl wrote an alignment
if [[ ! -s "${TRIM}" ]]; then
    echo "ERROR: trimAl produced no output for ${OG}"
    exit 1
fi

echo "${OG}: alignment and trimming complete."

### output files
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  per_og_aln/${OG}.aln.fa"
echo "  per_og_trimmed/${OG}.trimmed.fa"
