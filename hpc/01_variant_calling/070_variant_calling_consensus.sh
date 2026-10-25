#!/bin/bash
#SBATCH --job-name=variant_calling_consensus    ### job name
#SBATCH --cpus-per-task=4                       ### cores
#SBATCH --mem-per-cpu=8G                        ### memory per core
#SBATCH --time=7-00:00:00                       ### walltime
#SBATCH --output /dev/null                      ### stdout (logged by tee)
#SBATCH --error /dev/null                       ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/070_results_variant_calling_consensus"
LOGS="${ROOT}/01_variant_calling/logs/070_logs_variant_calling_consensus"
### from hpc/01_variant_calling/052_consensus_build.sh (indexed by hpc/01_variant_calling/053_indexing.sh)
REFERENCE="${ROOT}/01_variant_calling/results/050_results_sco_consensus/consensus_sco_reference.fasta"
### from hpc/01_variant_calling/061_bam_scaffold_list.sh
LISTS="${ROOT}/01_variant_calling/results/061_results_bam_scaffold_list"
### BAM list (all mapped samples)
BAM_LIST="${LISTS}/consensus_bam_list.txt"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"
### rows in the random subset of the metrics table and its seed
MAX_ROWS=500000
SEED=42

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_variant_calling_consensus.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_bcftools

START_TIME=$(date +%s)

### run info
echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"
echo "bcftools version: $(bcftools --version | head -1)"
echo "Reference: ${REFERENCE}"
echo "Number of samples: $(wc -l < "${BAM_LIST}")"


### =====================================================================
### step 1: variant calling over all consensus contigs
### =====================================================================

echo ""
echo "########## Variant calling ##########"

### one pass over all contigs (no -r, no array)
### mpileup (-q 5 minimum mapping quality, -I no indels, FMT/DP as the only per-sample field)
### call (--ploidy 1 haploid, -m multiallelic caller, -v variant sites only, -f GQ genotype quality)
### fill-tags (population-level INFO tags for filtering)
bcftools mpileup -Ou \
    --threads "${THREADS}" \
    -f "${REFERENCE}" \
    --bam-list "${BAM_LIST}" \
    -q 5 \
    -I \
    -a FMT/DP,INFO/AD,INFO/ADF,INFO/ADR,INFO/MQ0F,INFO/SCR,INFO/SGB \
| bcftools call \
    --threads "${THREADS}" \
    --ploidy 1 \
    -f GQ -mv \
    -Ou \
| bcftools +fill-tags \
    -Oz -o "${RESULTS}/consensus_raw.vcf.gz" \
    -- \
    -t F_MISSING,MAF,NS,AC,AN,AF

echo "mpileup, call and fill-tags complete."

### index
bcftools index "${RESULTS}/consensus_raw.vcf.gz"


### =====================================================================
### step 2: per-variant metrics
### =====================================================================

echo ""
echo "########## Per-variant metrics ##########"

### one row per variant with the INFO fields
METRICS="${RESULTS}/consensus_raw_all_variants_metrics.txt"

printf "CHROM\tPOS\tQUAL\tDP\tMQ\tMAF\tF_MISSING\tNS\tAC\tAN\tAF\tMQ0F\tSCR\tSGB\tFILE\n" > "${METRICS}"

bcftools query \
    -f '%CHROM\t%POS\t%QUAL\t%INFO/DP\t%INFO/MQ\t%INFO/MAF\t%INFO/F_MISSING\t%INFO/NS\t%INFO/AC\t%INFO/AN\t%INFO/AF\t%INFO/MQ0F\t%INFO/SCR\t%INFO/SGB\n' \
    "${RESULTS}/consensus_raw.vcf.gz" \
| awk 'BEGIN{OFS="\t"} {print $0, "consensus_raw"}' >> "${METRICS}"

TOTAL=$(( $(wc -l < "${METRICS}") - 1 ))
echo "Variants written to metrics file: ${TOTAL}"

### random subset (table small enough to plot on a laptop)
if [[ "${TOTAL}" -gt "${MAX_ROWS}" ]]; then
    SUBSET="${RESULTS}/consensus_raw_metrics_subsample.txt"
    head -1 "${METRICS}" > "${SUBSET}"
    ### seeded random source (GNU coreutils manual, "Sources of random data")
    tail -n +2 "${METRICS}" \
      | shuf -n "${MAX_ROWS}" \
          --random-source=<(openssl enc -aes-256-ctr -pass pass:"${SEED}" -nosalt </dev/zero 2>/dev/null) \
      >> "${SUBSET}"
    echo "Subsampled to ${MAX_ROWS} rows: ${SUBSET}"
else
    echo "Below ${MAX_ROWS} rows, no subsample needed."
fi


### =====================================================================
### step 3: summary
### =====================================================================

echo ""
echo "########## Result ##########"

N_VARIANTS=$(bcftools index -n "${RESULTS}/consensus_raw.vcf.gz")
N_CONTIGS=$(bcftools query -f '%CHROM\n' "${RESULTS}/consensus_raw.vcf.gz" | sort -u | wc -l)
TOTAL_CONTIGS=$(wc -l < "${REFERENCE}.fai")
REF_SIZE=$(awk '{sum += $2} END {print sum}' "${REFERENCE}.fai")

echo "Raw variants:       ${N_VARIANTS}"
echo "Contigs with vars:  ${N_CONTIGS} / ${TOTAL_CONTIGS}"
echo "Reference size:     ${REF_SIZE} bp"
echo "Variant density:    $(awk -v n="${N_VARIANTS}" -v s="${REF_SIZE}" 'BEGIN{printf "%.2f", n/s*1000}') variants/kb"
echo "Samples in the VCF: $(bcftools query -l "${RESULTS}/consensus_raw.vcf.gz" | wc -l)"
echo "FORMAT fields kept: $(bcftools view -h "${RESULTS}/consensus_raw.vcf.gz" | grep '^##FORMAT' | grep -o 'ID=[A-Za-z]*' | sed 's/ID=//' | tr '\n' ' ')"

### output files (copy consensus_raw_all_variants_metrics.txt to local/01_variant_calling/data/070_data_variant_calling for local/01_variant_calling/070_variant_calling.R and 080_filtering.R)
echo ""
echo "Output:"
echo "  VCF:     ${RESULTS}/consensus_raw.vcf.gz"
echo "  metrics: ${METRICS}"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
HOURS=$((RUNTIME / 3600))
MINUTES=$(((RUNTIME % 3600) / 60))
echo ""
echo "Variant calling on the consensus completed in ${RUNTIME} seconds (${HOURS}h ${MINUTES}m)."
