#!/bin/bash
#SBATCH --job-name=mating_type      ### job name
#SBATCH --array=1-608%60            ### array (one task per sample)
#SBATCH --cpus-per-task=4           ### cores
#SBATCH --mem-per-cpu=4G            ### memory per core
#SBATCH --time=04:00:00             ### walltime
#SBATCH --output /dev/null          ### stdout (logged by tee)
#SBATCH --error /dev/null           ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/102_results_mating_type"
LOGS="${ROOT}/01_variant_calling/logs/102_logs_mating_type"

### from hpc/01_variant_calling/010_read_trimming.sh
TRIM_DIR="${ROOT}/01_variant_calling/results/010_results_read_trimming"
### from hpc/01_variant_calling/000_prep.sh
SUPP="${ROOT}/01_variant_calling/supplementary/010_supplementary_read_trimming"

### idiomorph reference and its core, from hpc/01_variant_calling/101_mat_idiomorphs.sh
REF="${ROOT}/01_variant_calling/results/101_results_mat_idiomorphs/mat_idiomorphs.fasta"
CORE="${ROOT}/01_variant_calling/results/101_results_mat_idiomorphs/mat_core.bed"

PER_SAMPLE="${RESULTS}/per_sample"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"

### array task (first argument outside SLURM)
TASK_ID="${SLURM_ARRAY_TASK_ID:-$1}"

### minimum mapping quality
MAPQ=20

### breadth for presence
BREADTH_PRESENT=0.75

### breadth for absence
BREADTH_ABSENT=0.40

### mean depth for presence
DEPTH_MIN=3

### share of depth below which an idiomorph is absent
MINOR_MAX=0.10

### create results folders
mkdir -p "${RESULTS}" "${PER_SAMPLE}" "${LOGS}"

### define list(s)
### trimmed read names, from hpc/01_variant_calling/000_prep.sh
SAMPLE_LIST="${SUPP}/ceno_genomics_igat_trimmed_reads_list_R1.txt"

TRIM_R1=$(sed -n "${TASK_ID}p" "${SAMPLE_LIST}")
TRIM_R2=$(sed -n "${TASK_ID}p" \
          "${SUPP}/ceno_genomics_igat_trimmed_reads_list_R2.txt")

## check if the task is within the read lists
if [[ -z "${TRIM_R1}" || -z "${TRIM_R2}" ]]; then
    echo "ERROR: array task ${TASK_ID} is beyond the end of the read lists" \
      > "${LOGS}/${SLURM_JOB_ID:-local}_task${TASK_ID}.log"
    exit 1
fi

SAMPLE=$(basename "${TRIM_R1}" | sed 's/_R[12]$//')

R1="${TRIM_DIR}/${TRIM_R1}_trimmed.fastq.gz"
R2="${TRIM_DIR}/${TRIM_R2}_trimmed.fastq.gz"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_${SAMPLE}_mating_type.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate

START_TIME=$(date +%s)

echo "Date: $(date)"
echo "Sample: ${SAMPLE}"

## check if the indexed idiomorph reference and its core exist
if [[ ! -s "${REF}" || ! -s "${REF}.fai" || ! -s "${REF}.bwt" || ! -s "${CORE}" ]]; then
    echo "ERROR: no indexed idiomorph reference or core, run 101_mat_idiomorphs.sh first"
    echo "       ${REF}"
    exit 1
fi

## check if both trimmed files exist
if [[ ! -s "${R1}" || ! -s "${R2}" ]]; then
    echo "ERROR: trimmed reads missing for ${SAMPLE}"
    echo "       ${R1}"
    echo "       ${R2}"
    exit 1
fi

echo "R1: $(basename "${R1}")"
echo "R2: $(basename "${R2}")"


### =====================================================================
### step 1: mapping
### =====================================================================

echo ""
echo "=== Step 1: mapping ==="

conda deactivate
conda activate cg_bioinfo_v1_bwa
echo "Running in environment: $CONDA_DEFAULT_ENV"

BAM="${PER_SAMPLE}/${SAMPLE}.bam"

### map, filter (MAPQ, secondary, supplementary) and sort
bwa mem -t "${THREADS}" -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" "${REF}" "${R1}" "${R2}" 2>/dev/null \
    | samtools view -b -q "${MAPQ}" -F 0x900 - \
    | samtools sort -@ 2 -o "${BAM}" -

samtools index "${BAM}"

echo "Reads kept: $(samtools view -c "${BAM}")"


### =====================================================================
### step 2: coverage
### =====================================================================

echo ""
echo "=== Step 2: coverage ==="

conda deactivate
conda activate cg_bioinfo_v1_bedtools
echo "Running in environment: $CONDA_DEFAULT_ENV"

### depth at every core position
bedtools coverage -a "${CORE}" -b "${BAM}" -d > "${PER_SAMPLE}/${SAMPLE}.depth"


### =====================================================================
### step 3: the call
### =====================================================================

echo ""
echo "=== Step 3: calling ==="

OUT="${PER_SAMPLE}/${SAMPLE}_mating_type.tsv"

### breadth, mean depth and call per sample
awk -v s="${SAMPLE}" -v bp="${BREADTH_PRESENT}" -v ba="${BREADTH_ABSENT}" \
    -v dm="${DEPTH_MIN}" -v mm="${MINOR_MAX}" '
{
    n[$1]++
    sum[$1] += $5
    if ($5 > 0) cov[$1]++
}
END {
    b1 = n["MAT1_1"] ? cov["MAT1_1"] / n["MAT1_1"] : 0
    b2 = n["MAT1_2"] ? cov["MAT1_2"] / n["MAT1_2"] : 0
    d1 = n["MAT1_1"] ? sum["MAT1_1"] / n["MAT1_1"] : 0
    d2 = n["MAT1_2"] ? sum["MAT1_2"] / n["MAT1_2"] : 0

    tot    = d1 + d2
    minor  = (tot > 0) ? ((d1 < d2 ? d1 : d2) / tot) : 0
    minor1 = (tot > 0 && d1 / tot < mm)
    minor2 = (tot > 0 && d2 / tot < mm)

    ### absent on either criterion, present on all
    a1 = (b1 <= ba || minor1)
    a2 = (b2 <= ba || minor2)
    p1 = (b1 >= bp && d1 >= dm && !minor1)
    p2 = (b2 >= bp && d2 >= dm && !minor2)

    if      (p1 && a2) call = "MAT1"
    else if (p2 && a1) call = "MAT2"
    else if (p1 && p2) call = "MAT12"
    else if (a1 && a2) call = "none"
    else               call = "ambiguous"

    printf "%s\t%.2f\t%.3f\t%.2f\t%.3f\t%.3f\t%s\n",
           s, d1, b1, d2, b2, minor, call
}' "${PER_SAMPLE}/${SAMPLE}.depth" > "${OUT}"

cat "${OUT}"


### output files
echo ""
echo "=== Output files in ${PER_SAMPLE} ==="
echo "  ${SAMPLE}.bam"
echo "  ${SAMPLE}.bam.bai"
echo "  ${SAMPLE}.depth"
echo "  ${SAMPLE}_mating_type.tsv"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Array task ${TASK_ID} / sample ${SAMPLE} completed in ${RUNTIME} seconds."
