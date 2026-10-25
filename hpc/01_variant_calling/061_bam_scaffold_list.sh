#!/bin/bash
#SBATCH --job-name=bam_scaffold_list    ### job name
#SBATCH --cpus-per-task=1               ### cores
#SBATCH --mem-per-cpu=4G                ### memory per core
#SBATCH --time=00:30:00                 ### walltime
#SBATCH --output /dev/null              ### stdout (logged by tee)
#SBATCH --error /dev/null               ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/061_results_bam_scaffold_list"
LOGS="${ROOT}/01_variant_calling/logs/061_logs_bam_scaffold_list"
### from hpc/01_variant_calling/052_consensus_build.sh (indexed by hpc/01_variant_calling/053_indexing.sh)
REFERENCE="${ROOT}/01_variant_calling/results/050_results_sco_consensus/consensus_sco_reference.fasta" ### single-copy orthologue consensus
### from hpc/01_variant_calling/060_mapping.sh
BAMS="${ROOT}/01_variant_calling/results/060_results_mapping/consensus"

### define lists
### made by hand: sample names, one per line, the 060 array list (into ${ROOT}/01_variant_calling/supplementary/060_supplementary_mapping)
SAMPLE_LIST="${ROOT}/01_variant_calling/supplementary/060_supplementary_mapping/sample_names_list.txt"

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_bam_scaffold_list.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_samtools

START_TIME=$(date +%s)

echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"
echo "samtools version: $(samtools --version | head -1)"
echo "Sample list: ${SAMPLE_LIST}"


### =====================================================================
### step 1: contig list (one contig per orthologue, for reference only)
### =====================================================================

echo ""
echo "########## Contigs ##########"

### reference index (if missing)
if [[ ! -f "${REFERENCE}.fai" ]]; then
    echo "Creating reference index..."
    samtools faidx "${REFERENCE}"
fi

CONTIGS="${RESULTS}/consensus_scaffold_list.txt"
awk '{print $1}' "${REFERENCE}.fai" > "${CONTIGS}"

N_CONTIGS=$(wc -l < "${CONTIGS}")
REF_SIZE=$(awk '{sum += $2} END {print sum}' "${REFERENCE}.fai")

echo "Contigs:        ${N_CONTIGS}"
echo "Reference size: ${REF_SIZE} bp"


### =====================================================================
### step 2: BAM list (in sample list order)
### =====================================================================

echo ""
echo "########## BAMs ##########"

## check if the sample list exists
if [[ ! -s "${SAMPLE_LIST}" ]]; then
    echo "ERROR: sample list not found: ${SAMPLE_LIST}"
    exit 1
fi

BAM_LIST="${RESULTS}/consensus_bam_list.txt"
MISSING="${RESULTS}/consensus_missing_bams.txt"
: > "${BAM_LIST}"
: > "${MISSING}"

while read -r SAMPLE; do
    [[ -z "${SAMPLE}" ]] && continue
    BAM="${BAMS}/${SAMPLE}_sorted_markdup.bam"
    if [[ -s "${BAM}" ]]; then
        echo "${BAM}" >> "${BAM_LIST}"
    else
        echo "${SAMPLE}" >> "${MISSING}"
    fi
done < "${SAMPLE_LIST}"

N_SAMPLES=$(grep -cve '^[[:space:]]*$' "${SAMPLE_LIST}")
N_BAMS=$(wc -l < "${BAM_LIST}")
N_MISSING=$(wc -l < "${MISSING}")

echo "Samples in the list: ${N_SAMPLES}"
echo "BAMs found:          ${N_BAMS}"
echo "BAMs missing:        ${N_MISSING}"

if [[ "${N_MISSING}" -gt 0 ]]; then
    echo "Samples without a BAM:"
    sed 's/^/  /' "${MISSING}"
fi


### =====================================================================
### step 3: checks
### =====================================================================

echo ""
echo "########## Checks ##########"

## check if every sample has a BAM
if [[ "${N_BAMS}" -ne "${N_SAMPLES}" ]]; then
    echo "WARNING: ${N_BAMS} of ${N_SAMPLES} samples have a BAM."
    echo "         re-run the failed 060 array tasks before calling."
fi

## check if every BAM carries its own read-group sample name (the VCF sample names)
N_RG=$(while read -r BAM; do
           samtools view -H "${BAM}" | awk '/^@RG/ {for (i=1;i<=NF;i++) if ($i ~ /^SM:/) print substr($i,4)}'
       done < "${BAM_LIST}" | sort -u | wc -l)

if [[ "${N_RG}" -ne "${N_BAMS}" ]]; then
    echo "WARNING: ${N_RG} distinct read-group sample names for ${N_BAMS} BAMs."
    echo "         two BAMs may carry the same SM tag."
fi


### output files
echo ""
echo "########## Output ##########"
echo "Contig list: ${CONTIGS}   (${N_CONTIGS} contigs, for reference only)"
echo "BAM list:    ${BAM_LIST}   (${N_BAMS} samples)"
echo "Missing:     ${MISSING}   (${N_MISSING} samples)"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Prep completed in ${RUNTIME} seconds."
