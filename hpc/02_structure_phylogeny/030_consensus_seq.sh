#!/bin/bash
#SBATCH --job-name=consensus_seq    ### job name
#SBATCH --array=1-5616%100          ### array (one task per orthogroup in og_list.txt)
#SBATCH --cpus-per-task=1           ### cores
#SBATCH --mem-per-cpu=4G            ### memory per core
#SBATCH --time=02:00:00             ### walltime
#SBATCH --output /dev/null          ### stdout (logged by tee)
#SBATCH --error /dev/null           ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
VARCALL="${ROOT}/01_variant_calling"
RESULTS="${ROOT}/02_structure_phylogeny/results/030_results_consensus_seq"
LOGS="${ROOT}/02_structure_phylogeny/logs/030_logs_consensus_seq/logs_per_og"
### from hpc/01_variant_calling/052_consensus_build.sh (indexed by hpc/01_variant_calling/053_indexing.sh)
REFERENCE="${VARCALL}/results/050_results_sco_consensus/consensus_sco_reference.fasta"
### from hpc/01_variant_calling/050_cds_extraction.sh
OG_LIST="${VARCALL}/results/050_results_sco_consensus/og_list.txt"
### from hpc/02_structure_phylogeny/020_coverage_mask.sh
MASK_DIR="${ROOT}/02_structure_phylogeny/results/020_results_coverage_mask/per_sample_bed"
### from hpc/02_structure_phylogeny/010_filtering_phylogeny.sh (all raw sites, failed sites and genotypes set to missing)
VCF_PHYLO="${ROOT}/02_structure_phylogeny/results/010_results_filtering_phylogeny/consensus_phylogeny.vcf.gz"
### calls and sample set from the same file
VCF_INPUT="${VCF_PHYLO}"
SAMPLE_SOURCE="${VCF_PHYLO}"
FASTA_DIR="${RESULTS}/per_og_consensus"

### parameters
### array task (first argument outside SLURM)
TASK_ID="${SLURM_ARRAY_TASK_ID:-$1}"

### create results folders
mkdir -p "${FASTA_DIR}" "${LOGS}"

### define lists
### orthogroup of this array task
OG=$(sed -n "${TASK_ID}p" "${OG_LIST}")

## check if the array task is within the orthogroup list
if [[ -z "${OG}" ]]; then
    echo "ERROR: array task ${TASK_ID} is beyond the end of ${OG_LIST}" \
      > "${LOGS}/${SLURM_JOB_ID:-local}_task${TASK_ID}.log"
    exit 1
fi

SCG="${OG}_consensus" ### contig name in the reference

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_${OG}.log"
exec > >(tee "${LOGFILE}") 2>&1

### run info
echo "Date: $(date)"
echo "Orthogroup: ${OG}"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_samtools

## check if the phylogeny VCF is there
if [[ ! -s "${VCF_INPUT}" ]]; then
    echo "ERROR: phylogeny VCF not found: ${VCF_INPUT}"
    exit 1
fi


### =====================================================================
### step 1: single-gene reference and records
### =====================================================================

### single-gene reference
SCG_SEQ="${RESULTS}/${SCG}_ref.fasta"
samtools faidx "${REFERENCE}" "${SCG}" > "${SCG_SEQ}"

## check if the contig is in the reference
if [[ ! -s "${SCG_SEQ}" ]]; then
    echo "ERROR: ${SCG} not found in ${REFERENCE}"
    exit 1
fi

conda deactivate
conda activate cg_bioinfo_v1_bcftools

### records of this gene (small file for the per-sample loop)
SCG_VCF="${RESULTS}/${SCG}.vcf.gz"
bcftools view -r "${SCG}" "${VCF_INPUT}" -Oz -o "${SCG_VCF}"
bcftools index "${SCG_VCF}"

echo "Variants: $(bcftools index -n "${SCG_VCF}")"


### =====================================================================
### step 2: sample set and coverage masks
### =====================================================================

### sample list (read from a file, not a pipe, so an exit in the loop stops the task)
SAMPLES="${RESULTS}/${SCG}_samples.txt"
bcftools query -l "${SAMPLE_SOURCE}" > "${SAMPLES}"

N_SAMPLES=$(wc -l < "${SAMPLES}")
echo "Samples:  ${N_SAMPLES}"

## check if every sample has a coverage mask (before the loop, so a missing mask stops the task)
N_MISSING=0
while read -r SAMPLE; do
    if [[ ! -f "${MASK_DIR}/${SAMPLE}_lowcov.bed" ]]; then
        echo "ERROR: coverage mask not found: ${MASK_DIR}/${SAMPLE}_lowcov.bed"
        N_MISSING=$((N_MISSING + 1))
    fi
done < "${SAMPLES}"

if [[ "${N_MISSING}" -gt 0 ]]; then
    echo "ERROR: ${N_MISSING} samples have no coverage mask, run 020 for them"
    exit 1
fi


### =====================================================================
### step 3: consensus per sample
### =====================================================================

### reference sequence first
OUTPUT="${FASTA_DIR}/${SCG}_and_all_samples.fasta"
cat "${SCG_SEQ}" > "${OUTPUT}"

### one consensus per sample (N for missing genotypes and failed sites: -M N; N below the 020 depth: -m --mask-with N)
while read -r SAMPLE; do

    bcftools consensus -s "${SAMPLE}" -M N \
        -m "${MASK_DIR}/${SAMPLE}_lowcov.bed" --mask-with N \
        -f "${SCG_SEQ}" "${SCG_VCF}" \
        | sed "s/>${SCG}/>${OG}_${SAMPLE}/" >> "${OUTPUT}"

done < "${SAMPLES}"

### sequences written
N_WRITTEN=$(grep -c '^>' "${OUTPUT}")
N_EXPECTED=$((N_SAMPLES + 1))

echo "Sequences: ${N_WRITTEN} (expected ${N_EXPECTED}, the reference and one per sample)"

## check if every sample is in the file (a failed consensus writes nothing)
if [[ "${N_WRITTEN}" -ne "${N_EXPECTED}" ]]; then
    echo "ERROR: ${N_WRITTEN} sequences written, ${N_EXPECTED} expected"
    exit 1
fi

### remove intermediate files
rm -f "${SCG_SEQ}" "${SCG_SEQ}.fai" "${SCG_VCF}" "${SCG_VCF}.csi" "${SAMPLES}"

echo "${OG}: consensus sequences complete."

### output files
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  per_og_consensus/${SCG}_and_all_samples.fasta"
