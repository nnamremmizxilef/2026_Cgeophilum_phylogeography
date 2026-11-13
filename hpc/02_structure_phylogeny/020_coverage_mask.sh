#!/bin/bash
#SBATCH --job-name=coverage_mask    ### job name
#SBATCH --array=1-608%100           ### array (one task per sample)
#SBATCH --cpus-per-task=1           ### cores
#SBATCH --mem-per-cpu=4G            ### memory per core
#SBATCH --time=02:00:00             ### walltime
#SBATCH --output /dev/null          ### stdout (logged by tee)
#SBATCH --error /dev/null           ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
VARCALL="${ROOT}/01_variant_calling"
RESULTS="${ROOT}/02_structure_phylogeny/results/020_results_coverage_mask"
LOGS="${ROOT}/02_structure_phylogeny/logs/020_logs_coverage_mask/logs_per_sample"
### from hpc/01_variant_calling/060_mapping.sh
BAM_DIR="${VARCALL}/results/060_results_mapping/consensus"
### made by hand: sample names, one per line (the list of hpc/01_variant_calling/060_mapping.sh, in ${ROOT}/01_variant_calling/supplementary/060_supplementary_mapping)
SAMPLE_LIST="${VARCALL}/supplementary/060_supplementary_mapping/sample_names_list.txt"
BED_DIR="${RESULTS}/per_sample_bed"

### parameters
### array task (first argument outside SLURM)
TASK_ID="${SLURM_ARRAY_TASK_ID:-$1}"
### minimum depth of a called genotype (mirrors DP in hpc/02_structure_phylogeny/010_filtering_phylogeny.sh)
MIN_DP=3
### minimum mapping quality (mirrors mpileup -q in hpc/01_variant_calling/070_variant_calling_consensus.sh)
MIN_MQ=5

### create results folders
mkdir -p "${BED_DIR}" "${LOGS}"

### define lists
### sample of this array task
SAMPLE=$(sed -n "${TASK_ID}p" "${SAMPLE_LIST}")

## check if the array task is within the sample list
if [[ -z "${SAMPLE}" ]]; then
    echo "ERROR: array task ${TASK_ID} is beyond the end of ${SAMPLE_LIST}" \
      > "${LOGS}/${SLURM_JOB_ID:-local}_task${TASK_ID}.log"
    exit 1
fi

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_${SAMPLE}.log"
exec > >(tee "${LOGFILE}") 2>&1

### run info
echo "Date: $(date)"
echo "Sample: ${SAMPLE}"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_bedtools

BAM="${BAM_DIR}/${SAMPLE}_sorted_markdup.bam"
BED="${BED_DIR}/${SAMPLE}_lowcov.bed"

## check if the BAM is there
if [[ ! -s "${BAM}" ]]; then
    echo "ERROR: BAM not found: ${BAM}"
    exit 1
fi


### =====================================================================
### step 1: low-coverage mask
### =====================================================================

### positions with fewer than MIN_DP reads (written as N by 030; read filters as in the pileup)
samtools view -b -q "${MIN_MQ}" -F 0x704 "${BAM}" \
    | bedtools genomecov -ibam - -bga \
    | awk -v d="${MIN_DP}" 'BEGIN {OFS = "\t"} $4 < d {print $1, $2, $3}' \
    > "${BED}"

if [[ ! -s "${BED}" ]]; then
    echo "No positions below ${MIN_DP}x, empty mask written."
fi

### masked share of the reference (reference index from hpc/01_variant_calling/053_indexing.sh)
MASKED=$(awk '{n += $3 - $2} END {print n + 0}' "${BED}")
REF_SIZE=$(awk '{n += $2} END {print n}' \
    "${VARCALL}/results/050_results_sco_consensus/consensus_sco_reference.fasta.fai")

echo "Intervals: $(wc -l < "${BED}")"
echo "Masked:    ${MASKED} of ${REF_SIZE} bp ($(awk -v a="${MASKED}" -v b="${REF_SIZE}" 'BEGIN {printf "%.2f", a / b * 100}')%)"

echo "${SAMPLE}: coverage mask complete."

### output files
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  per_sample_bed/${SAMPLE}_lowcov.bed"
