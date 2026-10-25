#!/bin/bash
#SBATCH --job-name=multiqc          ### job name
#SBATCH --cpus-per-task=1           ### cores
#SBATCH --mem-per-cpu=16G           ### memory per core
#SBATCH --time=01:00:00             ### walltime
#SBATCH --output /dev/null          ### stdout (logged by tee)
#SBATCH --error /dev/null           ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/030_results_multiqc"
LOGS="${ROOT}/01_variant_calling/logs/030_logs_multiqc"
### from hpc/01_variant_calling/020_fastqc.sh
FASTQC_RAW="${ROOT}/01_variant_calling/results/020_results_fastqc/fastqc_results_raw_reads"        ### fastqc reports of raw reads
### from hpc/01_variant_calling/020_fastqc.sh
FASTQC_TRIMMED="${ROOT}/01_variant_calling/results/020_results_fastqc/fastqc_results_trimmed_reads" ### fastqc reports of trimmed reads
### from hpc/01_variant_calling/010_read_trimming.sh
FASTP="${ROOT}/01_variant_calling/results/010_results_read_trimming"                               ### fastp reports

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_multiqc.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_multiqc

START_TIME=$(date +%s)

echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"
echo "MultiQC version: $(multiqc --version 2>&1)"


### =====================================================================
### step 1: MultiQC on the FastQC reports of the raw reads
### =====================================================================

echo "Starting analysis on raw reads..."
echo "Processing: ${FASTQC_RAW}/"

multiqc \
  "${FASTQC_RAW}/" \
  -o "${RESULTS}" \
  --filename multiqc_report_raw_reads \
  --force \
  --interactive \
  --export \
  --fullnames \
  --verbose

echo "MultiQC analysis on raw reads complete!"


### =====================================================================
### step 2: MultiQC on the FastQC reports of the trimmed reads
### =====================================================================

echo "Starting analysis on trimmed reads..."
echo "Processing: ${FASTQC_TRIMMED}/"

multiqc \
  "${FASTQC_TRIMMED}/" \
  -o "${RESULTS}" \
  --filename multiqc_report_trimmed_reads \
  --force \
  --interactive \
  --export \
  --fullnames \
  --verbose

echo "MultiQC analysis on trimmed reads complete!"


### =====================================================================
### step 3: MultiQC on the fastp reports
### =====================================================================

echo "Starting analysis on fastp reports..."
echo "Processing: ${FASTP}/"

multiqc \
  "${FASTP}/" \
  -o "${RESULTS}" \
  --filename multiqc_report_fastp \
  --force \
  --interactive \
  --export \
  --verbose

echo "MultiQC analysis on fastp reports complete!"


### output files (copy the three multiqc_report_*_data folders to local/01_variant_calling/data/030_data_multiqc,
### renamed to start with pre (raw reads), post (trimmed reads) and fastp, for local/01_variant_calling/030_multiqc.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  multiqc_report_raw_reads.html      multiqc_report_raw_reads_data/"
echo "  multiqc_report_trimmed_reads.html  multiqc_report_trimmed_reads_data/"
echo "  multiqc_report_fastp.html          multiqc_report_fastp_data/"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo "MultiQC analysis completed in ${RUNTIME} seconds."
