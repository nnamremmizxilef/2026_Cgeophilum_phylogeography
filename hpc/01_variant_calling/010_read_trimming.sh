#!/bin/bash
#SBATCH --job-name=read_trimming    ### job name
#SBATCH --array=1-608%32            ### array (one task per sample)
#SBATCH --cpus-per-task=4           ### cores
#SBATCH --mem-per-cpu=4G            ### memory per core
#SBATCH --time=12:00:00             ### walltime
#SBATCH --output /dev/null          ### stdout (logged by tee)
#SBATCH --error /dev/null           ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/010_results_read_trimming"
LOGS="${ROOT}/01_variant_calling/logs/010_logs_read_trimming"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"
### array task (first argument outside SLURM)
TASK_ID="${SLURM_ARRAY_TASK_ID:-$1}"

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### define lists
### from hpc/01_variant_calling/000_prep.sh
DATA_RAW_R1=$(sed -n "${TASK_ID}p" "${ROOT}/01_variant_calling/supplementary/010_supplementary_read_trimming/ceno_genomics_igat_raw_reads_list_R1.txt")  ### list of raw R1 fastq files
DATA_RAW_R2=$(sed -n "${TASK_ID}p" "${ROOT}/01_variant_calling/supplementary/010_supplementary_read_trimming/ceno_genomics_igat_raw_reads_list_R2.txt")  ### list of raw R2 fastq files
DATA_TRIMMED_R1=$(sed -n "${TASK_ID}p" "${ROOT}/01_variant_calling/supplementary/010_supplementary_read_trimming/ceno_genomics_igat_trimmed_reads_list_R1.txt") ### list of R1 output names
DATA_TRIMMED_R2=$(sed -n "${TASK_ID}p" "${ROOT}/01_variant_calling/supplementary/010_supplementary_read_trimming/ceno_genomics_igat_trimmed_reads_list_R2.txt") ### list of R2 output names

### sample name (_R1/_R2 stripped)
SAMPLE_NAME=$(basename "${DATA_TRIMMED_R1}" | sed 's/_R[12]$//')

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_${SAMPLE_NAME}_read_trimming.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_fastp

START_TIME=$(date +%s)

echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"
echo "fastp version: $(fastp --version 2>&1)"
echo "Starting analysis..."
echo "Processing: ${DATA_RAW_R1} + ${DATA_RAW_R2} -> ${DATA_TRIMMED_R1}_trimmed.fastq.gz / ${DATA_TRIMMED_R2}_trimmed.fastq.gz"


### =====================================================================
### step 1: read trimming
### =====================================================================

fastp \
  --in1 "${DATA_RAW_R1}" \
  --in2 "${DATA_RAW_R2}" \
  --out1 "${RESULTS}/${DATA_TRIMMED_R1}_trimmed.fastq.gz" \
  --out2 "${RESULTS}/${DATA_TRIMMED_R2}_trimmed.fastq.gz" \
  --html "${RESULTS}/${DATA_TRIMMED_R1}_trimmed_fastp.html" \
  --json "${RESULTS}/${DATA_TRIMMED_R1}_trimmed_fastp.json" \
  --thread "${THREADS}" \
  --verbose \
  --n_base_limit 3 \
  --length_required 50


### output files
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  ${DATA_TRIMMED_R1}_trimmed.fastq.gz"
echo "  ${DATA_TRIMMED_R2}_trimmed.fastq.gz"
echo "  ${DATA_TRIMMED_R1}_trimmed_fastp.html"
echo "  ${DATA_TRIMMED_R1}_trimmed_fastp.json"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo "Array task ${TASK_ID} / sample ${SAMPLE_NAME} completed in ${RUNTIME} seconds."
