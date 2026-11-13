#!/bin/bash
#SBATCH --job-name=merge_alignments    ### job name
#SBATCH --array=1-5616%100             ### array (one task per orthogroup in og_list.txt)
#SBATCH --cpus-per-task=4              ### cores
#SBATCH --mem-per-cpu=4G               ### memory per core
#SBATCH --time=04:00:00                ### walltime
#SBATCH --output /dev/null             ### stdout (logged by tee)
#SBATCH --error /dev/null              ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
VARCALL="${ROOT}/01_variant_calling"
RESULTS="${ROOT}/02_structure_phylogeny/results/040_results_merge_alignments"
LOGS="${ROOT}/02_structure_phylogeny/logs/040_logs_merge_alignments/logs_per_og"
### from hpc/02_structure_phylogeny/030_consensus_seq.sh
CONSENSUS_DIR="${ROOT}/02_structure_phylogeny/results/030_results_consensus_seq/per_og_consensus"
### from hpc/01_variant_calling/050_cds_extraction.sh
OG_LIST="${VARCALL}/results/050_results_sco_consensus/og_list.txt"
### from hpc/01_variant_calling/051_consensus_align_trim.sh (coding sequences of the 26 genomes, aligned per orthogroup)
REFERENCE_DIR="${VARCALL}/results/050_results_sco_consensus/per_og_aln"
FINAL_DIR="${RESULTS}/per_og_final"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"
### array task (first argument outside SLURM)
TASK_ID="${SLURM_ARRAY_TASK_ID:-$1}"

### create results folders
mkdir -p "${FINAL_DIR}" "${LOGS}"

### define lists
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

### run info
echo "Date: $(date)"
echo "Orthogroup: ${OG}"

CONSENSUS_FASTA="${CONSENSUS_DIR}/${OG}_consensus_and_all_samples.fasta"
REFERENCE_ALIGNMENT="${REFERENCE_DIR}/${OG}.aln.fa"
MERGED="${RESULTS}/${OG}_merged.fasta"
RAW_ALIGNMENT="${RESULTS}/${OG}_raw_alignment.fasta"
FINAL_ALIGNMENT="${FINAL_DIR}/${OG}_final_alignment.fasta"

## check if the consensus file is there
if [[ ! -s "${CONSENSUS_FASTA}" ]]; then
    echo "ERROR: consensus file not found: ${CONSENSUS_FASTA}"
    exit 1
fi

## check if the reference alignment is there
if [[ ! -s "${REFERENCE_ALIGNMENT}" ]]; then
    echo "ERROR: reference alignment not found: ${REFERENCE_ALIGNMENT}"
    exit 1
fi


### =====================================================================
### step 1: merge
### =====================================================================

### genome sequences added to the samples (assemblies as taxa of the phylogeny)
cat "${CONSENSUS_FASTA}" "${REFERENCE_ALIGNMENT}" > "${MERGED}"

echo "Sequences: $(grep -c '^>' "${MERGED}")"


### =====================================================================
### step 2: alignment
### =====================================================================

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_mafft

echo "MAFFT version: $(mafft --version 2>&1)"

mafft --auto --thread "${THREADS}" "${MERGED}" > "${RAW_ALIGNMENT}"

## check if MAFFT wrote an alignment
if [[ ! -s "${RAW_ALIGNMENT}" ]]; then
    echo "ERROR: MAFFT produced no output for ${OG}"
    exit 1
fi


### =====================================================================
### step 3: trimming
### =====================================================================

conda deactivate
conda activate cg_bioinfo_v1_trimal

trimal -in "${RAW_ALIGNMENT}" -out "${FINAL_ALIGNMENT}" -automated1

## check if trimAl wrote an alignment
if [[ ! -s "${FINAL_ALIGNMENT}" ]]; then
    echo "ERROR: trimAl produced no output for ${OG}"
    exit 1
fi

### remove intermediate files
rm -f "${MERGED}" "${RAW_ALIGNMENT}"

echo "${OG}: merging, alignment and trimming complete."

### output files
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  per_og_final/${OG}_final_alignment.fasta"
