#!/bin/bash
#SBATCH --job-name=filtering_phylogeny    ### job name
#SBATCH --cpus-per-task=4                 ### cores
#SBATCH --mem-per-cpu=8G                  ### memory per core
#SBATCH --time=1-00:00:00                 ### walltime
#SBATCH --output /dev/null                ### stdout (logged by tee)
#SBATCH --error /dev/null                 ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
VARCALL="${ROOT}/01_variant_calling"
RESULTS="${ROOT}/02_structure_phylogeny/results/010_results_filtering_phylogeny"
LOGS="${ROOT}/02_structure_phylogeny/logs/010_logs_filtering_phylogeny"
### from hpc/01_variant_calling/070_variant_calling_consensus.sh (unfiltered call set, filtered here for the phylogeny)
VCF_RAW="${VARCALL}/results/070_results_variant_calling_consensus/consensus_raw.vcf.gz"

OUT_VCF="${RESULTS}/consensus_phylogeny.vcf.gz"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"
### thresholds of hpc/01_variant_calling/080_filtering.sh (failing sites and genotypes masked, not removed)
### site-level, masked
QUAL=30           ### minimum variant quality score
MQ=50             ### minimum root-mean-square mapping quality
### genotype-level, masked
GQ=20             ### minimum genotype quality
DP=3              ### minimum read depth
### sample-level (the only filter that drops samples, outgroup exempt)
SAMPLE_MISS=0.50
OUTGROUP_PREFIX="M4-"

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_filtering_phylogeny.log"
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

### parameters to the log
echo ""
echo "=== Filtering parameters ==="
echo "Site-level:    QUAL >= ${QUAL}, MQ >= ${MQ}, failing sites masked"
echo "Genotype:      GQ >= ${GQ}, DP >= ${DP}, failing genotypes masked"
echo "Sample-level:  ingroup samples above ${SAMPLE_MISS} missing are removed,"
echo "               outgroup samples (${OUTGROUP_PREFIX}*) never fail this filter"
echo "Not applied:   biallelic restriction, site missingness, MAF"

## check if the input is there
if [[ ! -s "${VCF_RAW}" ]]; then
    echo "ERROR: input VCF not found: ${VCF_RAW}"
    exit 1
fi

TEMP_DIR="${RESULTS}/temp_${SLURM_JOB_ID:-local}"
mkdir -p "${TEMP_DIR}"

N_INPUT=$(bcftools index -n "${VCF_RAW}")
echo ""
echo "Input:   ${VCF_RAW}"
echo "Sites:   ${N_INPUT}"
echo "Samples: $(bcftools query -l "${VCF_RAW}" | wc -l)"


### =====================================================================
### what the site rules will mask
### =====================================================================

### counted from the INFO fields (masking keeps the record count, so the index cannot count the steps)
echo ""
echo "########## Sites the site rules will mask ##########"

bcftools query -f '%QUAL\t%INFO/MQ\t%ALT\n' "${VCF_RAW}" \
| awk -F'\t' -v q="${QUAL}" -v m="${MQ}" '
{
    n++
    if (split($3, a, ",") > 1) ma++
    if ($1 + 0 <  q) lq++
    if ($2 == "." )  nomq++
    else if ($2 + 0 < m) lm++
    if ($1 + 0 < q || ($2 != "." && $2 + 0 < m)) masked++
}
END {
    printf "  sites                       %d\n", n
    printf "  QUAL < %-5s                %d (%.1f%%)\n", q, lq, 100 * lq / n
    printf "  MQ < %-5s                  %d (%.1f%%)\n", m, lm, 100 * lm / n
    printf "  MQ absent                   %d\n", nomq
    printf "  masked by either            %d (%.1f%%)\n", masked, 100 * masked / n
    printf "  multiallelic, kept here     %d (%.1f%%)\n", ma, 100 * ma / n
}' | tee "${RESULTS}/site_rule_counts.txt"


### =====================================================================
### site and genotype masking
### =====================================================================

echo ""
echo "########## Masking ##########"

### step 1: site-level rules (-S . sets every genotype at a failing site to missing and keeps the record; sites without INFO/MQ not masked)
bcftools filter -e "QUAL < ${QUAL} || MQ < ${MQ}" -S . \
    --threads "${THREADS}" \
    "${VCF_RAW}" -Oz -o "${TEMP_DIR}/step1.vcf.gz"
bcftools index --threads "${THREADS}" "${TEMP_DIR}/step1.vcf.gz"
echo "Step 1 - QUAL and MQ:              failing sites masked, records kept"

### step 2: genotype-level rules (single & so both conditions are evaluated per sample)
bcftools filter -i "FMT/GQ >= ${GQ} & FMT/DP >= ${DP}" -S . \
    --threads "${THREADS}" \
    "${TEMP_DIR}/step1.vcf.gz" -Oz -o "${TEMP_DIR}/step2.vcf.gz"
bcftools index --threads "${THREADS}" "${TEMP_DIR}/step2.vcf.gz"
echo "Step 2 - GQ and DP:                failing genotypes masked"

### step 3: count-based tags recalculated after masking
bcftools +fill-tags "${TEMP_DIR}/step2.vcf.gz" \
    -Oz -o "${TEMP_DIR}/step3.vcf.gz" \
    -- -t F_MISSING,MAF,NS,AC,AN,AF
bcftools index --threads "${THREADS}" "${TEMP_DIR}/step3.vcf.gz"
echo "Step 3 - tags recalculated:        F_MISSING, MAF, NS, AC, AN, AF"


### =====================================================================
### per-sample missingness
### =====================================================================

### after masking, over the sites with a call for at least one sample (NS > 0)
echo ""
echo "########## Per-sample missingness ##########"

SMISS="${RESULTS}/phylogeny_sample_missingness.tsv"

bcftools query -l "${TEMP_DIR}/step3.vcf.gz" > "${TEMP_DIR}/sample_names.txt"

bcftools query -i 'NS > 0' -f '[\t%GT]\n' "${TEMP_DIR}/step3.vcf.gz" \
| awk '{
    for (i = 1; i <= NF; i++) {
        total[i]++
        if ($i == "." || $i == "./.") miss[i]++
    }
}
END {
    for (i = 1; i <= NF; i++) print total[i], miss[i] + 0, (miss[i] + 0) / total[i]
}' > "${TEMP_DIR}/sample_miss_counts.txt"

printf "sample\tn_genotypes\tn_missing\tf_missing\n" > "${SMISS}"
paste "${TEMP_DIR}/sample_names.txt" "${TEMP_DIR}/sample_miss_counts.txt" \
| awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4}' >> "${SMISS}"

echo "Written: $(basename "${SMISS}")  ($(( $(wc -l < "${SMISS}") - 1 )) samples)"
echo "Counted over $(bcftools view -H -i 'NS > 0' "${TEMP_DIR}/step3.vcf.gz" | wc -l) sites with a call for at least one sample"

### ingroup samples above the threshold (outgroup exempt)
awk -F'\t' -v thresh="${SAMPLE_MISS}" -v opfx="${OUTGROUP_PREFIX}" \
    'NR>1 && index($1, opfx) != 1 && $4 + 0 > thresh + 0 {print $1}' \
    "${SMISS}" > "${TEMP_DIR}/exclude.txt"

N_EXCLUDE=$(wc -l < "${TEMP_DIR}/exclude.txt")
echo "Ingroup above ${SAMPLE_MISS} missing:      ${N_EXCLUDE}"

if [[ "${N_EXCLUDE}" -gt 0 ]]; then
    echo "Samples removed:"
    while read -r S; do
        awk -F'\t' -v s="${S}" '$1 == s {printf "  %-14s %.4f\n", $1, $4}' "${SMISS}"
    done < "${TEMP_DIR}/exclude.txt"
fi

### step 4: sample selection (removes samples, never sites)
if [[ "${N_EXCLUDE}" -gt 0 ]]; then
    bcftools view -S ^"${TEMP_DIR}/exclude.txt" --force-samples \
        --threads "${THREADS}" \
        "${TEMP_DIR}/step3.vcf.gz" -Oz -o "${TEMP_DIR}/step4.vcf.gz"
else
    cp "${TEMP_DIR}/step3.vcf.gz" "${TEMP_DIR}/step4.vcf.gz"
fi
bcftools index --threads "${THREADS}" "${TEMP_DIR}/step4.vcf.gz"
echo "Step 4 - sample selection:         $(bcftools query -l "${TEMP_DIR}/step4.vcf.gz" | wc -l) samples"

### step 5: tags recalculated after sample removal
bcftools +fill-tags "${TEMP_DIR}/step4.vcf.gz" \
    -Oz -o "${OUT_VCF}" \
    -- -t F_MISSING,MAF,NS,AC,AN,AF
bcftools index --threads "${THREADS}" "${OUT_VCF}"
echo "Step 5 - tags recalculated:        F_MISSING, MAF, NS, AC, AN, AF"


### =====================================================================
### checks
### =====================================================================

echo ""
echo "########## Checks ##########"

N_OUT=$(bcftools index -n "${OUT_VCF}")

## check if any site was lost (sites out must equal sites in)
if [[ "${N_OUT}" -ne "${N_INPUT}" ]]; then
    echo "ERROR: ${N_INPUT} sites in, ${N_OUT} sites out. a step deleted records"
    exit 1
fi

### genotypes missing, masked sites included
bcftools query -f '[%GT\t]\n' "${OUT_VCF}" \
| awk '{ for (i = 1; i <= NF; i++) { n++; if ($i == ".") m++ } }
       END { printf "Genotypes missing, masked sites included: %d of %d (%.1f%%)\n", m, n, 100 * m / n }'

### outgroup samples in the output
echo "Outgroup samples in the output:    $(bcftools query -l "${OUT_VCF}" | grep -c "^${OUTGROUP_PREFIX}")"


### =====================================================================
### per-variant metrics
### =====================================================================

echo ""
echo "########## Per-variant metrics ##########"

METRICS="${RESULTS}/consensus_phylogeny_metrics.txt.gz"

{
    printf "CHROM\tPOS\tQUAL\tDP\tMQ\tMAF\tF_MISSING\tNS\tAC\tAN\tAF\tFILE\n"
    bcftools query \
        -f '%CHROM\t%POS\t%QUAL\t%INFO/DP\t%INFO/MQ\t%INFO/MAF\t%INFO/F_MISSING\t%INFO/NS\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' \
        "${OUT_VCF}" \
    | awk 'BEGIN{OFS="\t"} {print $0, "consensus_phylogeny"}'
} | gzip > "${METRICS}"

echo "Written: $(basename "${METRICS}")"


### =====================================================================
### clean up
### =====================================================================

rm -rf "${TEMP_DIR}"


### output files
echo ""
echo "########## Output ##########"
echo ""
echo "  VCF:     ${OUT_VCF}"
echo "  sites:   ${N_OUT}"
echo "  samples: $(bcftools query -l "${OUT_VCF}" | wc -l)"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
HOURS=$((RUNTIME / 3600))
MINUTES=$(((RUNTIME % 3600) / 60))
echo ""
echo "Filtering for the phylogeny completed in ${RUNTIME} seconds (${HOURS}h ${MINUTES}m)."
