#!/bin/bash
#SBATCH --job-name=mapping          ### job name
#SBATCH --array=1-608%32            ### array (one task per sample)
#SBATCH --cpus-per-task=4           ### cores
#SBATCH --mem-per-cpu=8G            ### memory per core
#SBATCH --time=7-00:00:00           ### walltime
#SBATCH --output /dev/null          ### stdout (logged by tee)
#SBATCH --error /dev/null           ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/060_results_mapping/consensus"
LOGS="${ROOT}/01_variant_calling/logs/060_logs_mapping"
### from hpc/01_variant_calling/052_consensus_build.sh (indexed by hpc/01_variant_calling/053_indexing.sh)
REFERENCE="${ROOT}/01_variant_calling/results/050_results_sco_consensus/consensus_sco_reference.fasta" ### single-copy orthologue consensus

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"
### array task (first argument outside SLURM)
TASK_ID="${SLURM_ARRAY_TASK_ID:-$1}"

### create results folders (QC files in their own subfolders)
mkdir -p "${RESULTS}/flagstat" "${RESULTS}/coverage" "${LOGS}"

### define lists
### made by hand: trimmed R1 and R2 files (hpc/01_variant_calling/010_read_trimming.sh) and sample names, one sample per line (into ${ROOT}/01_variant_calling/supplementary/060_supplementary_mapping)
DATA_TRIMMED_R1=$(sed -n "${TASK_ID}p" "${ROOT}/01_variant_calling/supplementary/060_supplementary_mapping/trimmed_reads_R1_list.txt") ### list of trimmed R1 fastq files
DATA_TRIMMED_R2=$(sed -n "${TASK_ID}p" "${ROOT}/01_variant_calling/supplementary/060_supplementary_mapping/trimmed_reads_R2_list.txt") ### list of trimmed R2 fastq files
SAMPLE_NAME=$(sed -n "${TASK_ID}p" "${ROOT}/01_variant_calling/supplementary/060_supplementary_mapping/sample_names_list.txt") ### list of sample names

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_${SAMPLE_NAME}_mapping.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_bwa

START_TIME=$(date +%s)

echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"
echo "BWA version: $(bwa 2>&1 | head -3 | tail -1)"
echo "Starting analysis..."
echo "Processing sample: ${SAMPLE_NAME}"
echo "  R1: ${DATA_TRIMMED_R1}"
echo "  R2: ${DATA_TRIMMED_R2}"


### =====================================================================
### step 1: mapping to the consensus reference
### =====================================================================

echo ""
echo "########## Mapping to consensus reference ##########"

### map (read group carries the sample name)
bwa mem -t "${THREADS}" \
  -R "@RG\tID:${SAMPLE_NAME}\tSM:${SAMPLE_NAME}\tPL:ILLUMINA\tLB:${SAMPLE_NAME}" \
  "${REFERENCE}" \
  "${DATA_TRIMMED_R1}" \
  "${DATA_TRIMMED_R2}" \
  > "${RESULTS}/${SAMPLE_NAME}.sam"
echo "BWA MEM complete."


### =====================================================================
### step 2: duplicate marking and indexing
### =====================================================================

conda deactivate
conda activate cg_bioinfo_v1_samtools

### mate-score tags (needed by markdup)
samtools fixmate -m "${RESULTS}/${SAMPLE_NAME}.sam" "${RESULTS}/${SAMPLE_NAME}_fixmate.bam"

### sort by coordinate
samtools sort -@ "${THREADS}" "${RESULTS}/${SAMPLE_NAME}_fixmate.bam" -o "${RESULTS}/${SAMPLE_NAME}_sorted.bam"
rm "${RESULTS}/${SAMPLE_NAME}.sam" "${RESULTS}/${SAMPLE_NAME}_fixmate.bam"

### PCR and optical duplicates (marked, not removed)
samtools markdup -@ "${THREADS}" "${RESULTS}/${SAMPLE_NAME}_sorted.bam" "${RESULTS}/${SAMPLE_NAME}_sorted_markdup.bam"
rm "${RESULTS}/${SAMPLE_NAME}_sorted.bam"

### index
samtools index -@ "${THREADS}" "${RESULTS}/${SAMPLE_NAME}_sorted_markdup.bam"
echo "Duplicate marking and indexing complete."


### =====================================================================
### step 3: QC files
### =====================================================================

### flagstat and per-contig coverage (one contig per orthologue; duplicates and secondary alignments excluded by default)
samtools flagstat -@ "${THREADS}" "${RESULTS}/${SAMPLE_NAME}_sorted_markdup.bam" \
  > "${RESULTS}/flagstat/${SAMPLE_NAME}_flagstat.txt"
samtools coverage "${RESULTS}/${SAMPLE_NAME}_sorted_markdup.bam" \
  > "${RESULTS}/coverage/${SAMPLE_NAME}_coverage.txt"
cat "${RESULTS}/flagstat/${SAMPLE_NAME}_flagstat.txt"


### output files (copy ${RESULTS}/flagstat and ${RESULTS}/coverage to local/01_variant_calling/data/060_data_mapping for local/01_variant_calling/060_mapping.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  ${SAMPLE_NAME}_sorted_markdup.bam (.bai)"
echo "  flagstat/${SAMPLE_NAME}_flagstat.txt"
echo "  coverage/${SAMPLE_NAME}_coverage.txt"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Array task ${TASK_ID} / sample ${SAMPLE_NAME} completed in ${RUNTIME} seconds."
