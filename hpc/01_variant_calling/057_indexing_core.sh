#!/bin/bash
#SBATCH --job-name=ref_indexing     ### job name
#SBATCH --cpus-per-task=4           ### cores
#SBATCH --mem-per-cpu=4G            ### memory per core
#SBATCH --time=02:00:00             ### walltime
#SBATCH --output /dev/null          ### stdout (logged by tee)
#SBATCH --error /dev/null           ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
LOGS="${ROOT}/01_variant_calling/logs/057_logs_indexing_core"
### from hpc/01_variant_calling/056_consensus_build_core.sh
REF_CONS="${ROOT}/01_variant_calling/results/054_results_sco_consensus_core/consensus_sco_reference_core.fasta"

### create results folders
mkdir -p "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_ref_indexing.log"
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


### =====================================================================
### function: index a reference (bwa index, samtools faidx, samtools dict)
### =====================================================================

index_reference() {

    local REF=$1
    local LABEL=$2

    echo ""
    echo "================================================================="
    echo "Indexing: ${LABEL}"
    echo "File:    ${REF}"
    echo "================================================================="

    ## check if the reference exists
    if [[ ! -f "${REF}" ]]; then
        echo "ERROR: Reference file not found: ${REF}"
        return 1
    fi

    ### bwa index (skipped if present)
    if [[ -f "${REF}.bwt" ]]; then
        echo "  BWA index already exists, skipping."
    else
        echo "  Running BWA index..."
        bwa index "${REF}"
        echo "  BWA index complete."
    fi

    ### samtools faidx (skipped if present)
    if [[ -f "${REF}.fai" ]]; then
        echo "  Samtools faidx already exists, skipping."
    else
        echo "  Running samtools faidx..."
        samtools faidx "${REF}"
        echo "  Samtools faidx complete."
    fi

    ### samtools dict (skipped if present)
    local DICT="${REF%.fasta}.dict"
    if [[ -f "${DICT}" ]]; then
        echo "  Sequence dictionary already exists, skipping."
    else
        echo "  Creating sequence dictionary..."
        samtools dict "${REF}" > "${DICT}"
        echo "  Sequence dictionary created."
    fi

    ### summary
    local N_SEQS=$(wc -l < "${REF}.fai")
    local REF_SIZE=$(awk '{sum += $2} END {print sum}' "${REF}.fai")
    echo ""
    echo "  ${LABEL}: ${REF_SIZE} bp across ${N_SEQS} sequences"
    echo "  Index files:"
    ls -lh "${REF}".{bwt,pac,ann,amb,sa} "${REF}.fai" "${DICT}" 2>/dev/null | awk '{print "    " $NF " (" $5 ")"}'
}


### =====================================================================
### step 1: index the core reference
### =====================================================================

index_reference "${REF_CONS}" "Consensus SCO panel"


### output files
echo ""
echo "=== Output files ==="
echo "  bwa index, .fai and .dict next to ${REF_CONS}"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Indexing completed in ${RUNTIME} seconds."
