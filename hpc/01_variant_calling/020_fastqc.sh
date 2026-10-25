#!/bin/bash
#SBATCH --job-name=fastqc           ### job name
#SBATCH --array=1-1216%32           ### array (one task per read file)
#SBATCH --cpus-per-task=4           ### cores
#SBATCH --mem-per-cpu=8G            ### memory per core
#SBATCH --time=04:00:00             ### walltime
#SBATCH --output /dev/null          ### stdout (logged by tee)
#SBATCH --error /dev/null           ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS_RAW="${ROOT}/01_variant_calling/results/020_results_fastqc/fastqc_results_raw_reads"
RESULTS_TRIMMED="${ROOT}/01_variant_calling/results/020_results_fastqc/fastqc_results_trimmed_reads"
LOGS="${ROOT}/01_variant_calling/logs/020_logs_fastqc"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"
### array task (first argument outside SLURM)
TASK_ID="${SLURM_ARRAY_TASK_ID:-$1}"

### create results folders
mkdir -p "${RESULTS_RAW}" "${RESULTS_TRIMMED}" "${LOGS}"

### define lists
### made by hand: raw read files, trimmed read files (hpc/01_variant_calling/010_read_trimming.sh) and sample names, one read file per line (into ${ROOT}/01_variant_calling/supplementary/020_supplementary_fastqc)
RAW_FILE=$(sed -n "${TASK_ID}p" "${ROOT}/01_variant_calling/supplementary/020_supplementary_fastqc/list1_raw_reads.txt")        ### list of raw fastq files
TRIMMED_FILE=$(sed -n "${TASK_ID}p" "${ROOT}/01_variant_calling/supplementary/020_supplementary_fastqc/list2_trimmed_reads.txt") ### list of trimmed fastq files

### sample name (read from list)
SAMPLE_NAME=$(sed -n "${TASK_ID}p" "${ROOT}/01_variant_calling/supplementary/020_supplementary_fastqc/list3_sample_names.txt")

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_${SAMPLE_NAME}_fastqc.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### locale
export LANG=C
export LC_ALL=C
export LANGUAGE=C

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_fastqc

START_TIME=$(date +%s)

echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"
echo "FastQC version: $(fastqc --version 2>&1)"


### =====================================================================
### step 1: FastQC on the raw reads
### =====================================================================

echo "Starting analysis on raw reads..."
echo "Processing: ${RAW_FILE}"

fastqc \
  -t "${THREADS}" \
  -o "${RESULTS_RAW}" \
  "${RAW_FILE}"

echo "FastQC analysis on raw reads complete!"


### =====================================================================
### step 2: FastQC on the trimmed reads
### =====================================================================

echo "Starting analysis on trimmed reads..."
echo "Processing: ${TRIMMED_FILE}"

fastqc \
  -t "${THREADS}" \
  -o "${RESULTS_TRIMMED}" \
  "${TRIMMED_FILE}"

echo "FastQC analysis on trimmed reads complete!"


### output files
echo ""
echo "=== Output files ==="
echo "  FastQC report (html, zip) of $(basename "${RAW_FILE}") in ${RESULTS_RAW}"
echo "  FastQC report (html, zip) of $(basename "${TRIMMED_FILE}") in ${RESULTS_TRIMMED}"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo "Array task ${TASK_ID} / sample ${SAMPLE_NAME} completed in ${RUNTIME} seconds."
