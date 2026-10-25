#!/bin/bash
#SBATCH --job-name=filtering        ### job name
#SBATCH --cpus-per-task=4           ### cores
#SBATCH --mem-per-cpu=8G            ### memory per core
#SBATCH --time=3-00:00:00           ### walltime
#SBATCH --output /dev/null          ### stdout (logged by tee)
#SBATCH --error /dev/null           ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/080_results_filtering"
LOGS="${ROOT}/01_variant_calling/logs/080_logs_filtering"
### from hpc/01_variant_calling/070_variant_calling_consensus.sh (all mapped samples)
VCF_CONSENSUS="${ROOT}/01_variant_calling/results/070_results_variant_calling_consensus/consensus_raw.vcf.gz"
### from hpc/01_variant_calling/060_mapping.sh (per-orthologue coverage, one row per orthologue)
COVERAGE_DIR="${ROOT}/01_variant_calling/results/060_results_mapping/consensus/coverage"
### from hpc/01_variant_calling/052_consensus_build.sh (indexed by hpc/01_variant_calling/053_indexing.sh)
SCO_REFERENCE="${ROOT}/01_variant_calling/results/050_results_sco_consensus/consensus_sco_reference.fasta"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"
### outgroup sample prefix
OUTGROUP_PREFIX="M4-"
### site-level filters
QUAL=30           ### minimum variant quality score
MQ=50             ### minimum root-mean-square mapping quality
SITE_MISS=0.10    ### maximum fraction of missing genotypes per site
### genotype-level filters (failing genotypes set to missing)
GQ=20             ### minimum genotype quality
DP=3              ### minimum read depth
### sample-level filter (maximum fraction of missing genotypes, outgroup exempt)
SAMPLE_MISS=0.50
### breadth (%) at which an orthologue counts as present
SCO_PRESENT=50
### MAF levels
MAF_V1=0.01       ### permissive
MAF_V2=0.05       ### stringent

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_filtering.log"
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
echo "Site-level:    biallelic SNPs, QUAL >= ${QUAL}, MQ >= ${MQ},"
echo "               site missingness <= ${SITE_MISS}"
echo "Genotype:      GQ >= ${GQ}, DP >= ${DP}, failing genotypes set to missing"
echo "Sample-level:  ingroup samples above ${SAMPLE_MISS} missing are removed,"
echo "               outgroup samples (${OUTGROUP_PREFIX}*) never fail this filter"
echo "SCO presence:  a SCO counts as present at >= ${SCO_PRESENT}% breadth"
echo "MAF levels:    none, >= ${MAF_V1}, >= ${MAF_V2}"


### =====================================================================
### single copy orthologues present per sample
### =====================================================================

### from the 060 coverage tables (column 6 = percentage of the contig covered; diagnostic, not a filter)
echo ""
echo "########## Single copy orthologues per sample ##########"

SCO_TABLE="${RESULTS}/sco_presence.tsv"

### orthologues in the reference
N_SCO_REF=$(wc -l < "${SCO_REFERENCE}.fai")
echo "Orthologues in the consensus reference: ${N_SCO_REF}"

printf "sample\tn_scos_total\tn_scos_present\tpct_scos_present\n" > "${SCO_TABLE}"

for F in "${COVERAGE_DIR}"/*_coverage.txt; do
    SAMPLE=$(basename "${F}" _coverage.txt)
    awk -v s="${SAMPLE}" -v thr="${SCO_PRESENT}" 'BEGIN{OFS="\t"}
    /^#/ { next }
    {
        total++
        if ($6 + 0 >= thr) present++
    }
    END {
        print s, total, present + 0, (total > 0 ? present / total * 100 : 0)
    }' "${F}"
done >> "${SCO_TABLE}"

echo "Samples processed: $(( $(wc -l < "${SCO_TABLE}") - 1 ))"
echo "Written: $(basename "${SCO_TABLE}")"

echo "Distinct orthologue totals across samples:"
awk -F'\t' 'NR>1 {print $2}' "${SCO_TABLE}" | sort -u | sed 's/^/  /'

## check if every sample lists every orthologue of the reference (a truncated coverage file lowers the count)
N_BAD=$(awk -F'\t' -v n="${N_SCO_REF}" 'NR>1 && $2 + 0 != n' "${SCO_TABLE}" | wc -l)
if [[ "${N_BAD}" -gt 0 ]]; then
    echo "WARNING: ${N_BAD} samples do not list all ${N_SCO_REF} orthologues:"
    awk -F'\t' -v n="${N_SCO_REF}" 'NR>1 && $2 + 0 != n {printf "  %-14s %d rows\n", $1, $2}' \
        "${SCO_TABLE}"
    echo "         their coverage files are truncated, re-run 060 for them"
fi

echo ""
echo "Ten samples with the fewest orthologues present:"
awk -F'\t' 'NR>1' "${SCO_TABLE}" | sort -t$'\t' -k4,4g | head -10 \
| awk -F'\t' '{printf "  %-14s %6d of %d  (%.1f%%)\n", $1, $3, $2, $4}'


### =====================================================================
### site and genotype filters
### =====================================================================

echo ""
echo "########## Site and genotype filters ##########"
echo "Input: ${VCF_CONSENSUS}"

## check if the input VCF exists
if [[ ! -s "${VCF_CONSENSUS}" ]]; then
    echo "ERROR: input VCF not found: ${VCF_CONSENSUS}"
    exit 1
fi

TEMP_DIR="${RESULTS}/temp_${SLURM_JOB_ID:-local}"
mkdir -p "${TEMP_DIR}"

### site-filtered base of the sample sets (kept in the results folder)
BASE_VCF="${RESULTS}/consensus_sitefiltered_base.vcf.gz"

### record counts from the index
N_INPUT=$(bcftools index -n "${VCF_CONSENSUS}")
echo ""
echo "Samples in the input: $(bcftools query -l "${VCF_CONSENSUS}" | wc -l)"
echo "Step 0 - input:                    ${N_INPUT} variants"

### step 1: biallelic SNPs only
bcftools view -m2 -M2 -v snps --threads "${THREADS}" \
    "${VCF_CONSENSUS}" -Oz -o "${TEMP_DIR}/step1.vcf.gz"
bcftools index --threads "${THREADS}" "${TEMP_DIR}/step1.vcf.gz"
N1=$(bcftools index -n "${TEMP_DIR}/step1.vcf.gz")
echo "Step 1 - biallelic SNPs:           ${N1} variants (removed $((N_INPUT - N1)))"

### step 2: variant quality
bcftools filter -i "QUAL >= ${QUAL}" --threads "${THREADS}" \
    "${TEMP_DIR}/step1.vcf.gz" -Oz -o "${TEMP_DIR}/step2.vcf.gz"
bcftools index --threads "${THREADS}" "${TEMP_DIR}/step2.vcf.gz"
N2=$(bcftools index -n "${TEMP_DIR}/step2.vcf.gz")
echo "Step 2 - QUAL >= ${QUAL}:                ${N2} variants (removed $((N1 - N2)))"

### step 3: mapping quality of the reads at the site
bcftools filter -i "MQ >= ${MQ}" --threads "${THREADS}" \
    "${TEMP_DIR}/step2.vcf.gz" -Oz -o "${TEMP_DIR}/step3.vcf.gz"
bcftools index --threads "${THREADS}" "${TEMP_DIR}/step3.vcf.gz"
N3=$(bcftools index -n "${TEMP_DIR}/step3.vcf.gz")
echo "Step 3 - MQ >= ${MQ}:                ${N3} variants (removed $((N2 - N3)))"

### step 4: genotype masking (single & so both conditions are evaluated per sample)
bcftools filter -i "FMT/GQ >= ${GQ} & FMT/DP >= ${DP}" -S . \
    --threads "${THREADS}" \
    "${TEMP_DIR}/step3.vcf.gz" -Oz -o "${TEMP_DIR}/step4.vcf.gz"
bcftools index --threads "${THREADS}" "${TEMP_DIR}/step4.vcf.gz"
echo "Step 4 - GQ >= ${GQ} and DP >= ${DP}:      genotypes below either set to missing"

### step 5: count-based tags recalculated after masking
bcftools +fill-tags "${TEMP_DIR}/step4.vcf.gz" \
    -Oz -o "${BASE_VCF}" \
    -- -t F_MISSING,MAF,NS,AC,AN,AF
bcftools index --threads "${THREADS}" "${BASE_VCF}"
echo "Step 5 - tags recalculated:        F_MISSING, MAF, NS, AC, AN, AF"


### =====================================================================
### per-sample missingness
### =====================================================================

### on the base VCF (after genotype masking, before sample removal)
echo ""
echo "########## Per-sample missingness ##########"

SMISS="${RESULTS}/consensus_sample_missingness.tsv"

bcftools query -l "${BASE_VCF}" > "${TEMP_DIR}/sample_names.txt"

bcftools query -f '[\t%GT]\n' "${BASE_VCF}" \
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

### outgroup samples (removed from the ingroup set)
awk -F'\t' -v opfx="${OUTGROUP_PREFIX}" 'NR>1 && index($1, opfx) == 1 {print $1}' \
    "${SMISS}" > "${TEMP_DIR}/list_outgroup.txt"

### ingroup samples above the missingness threshold
awk -F'\t' -v thresh="${SAMPLE_MISS}" -v opfx="${OUTGROUP_PREFIX}" \
    'NR>1 && index($1, opfx) != 1 && $4 + 0 > thresh + 0 {print $1}' \
    "${SMISS}" > "${TEMP_DIR}/list_failing.txt"

echo ""
echo "Outgroup samples:                 $(wc -l < "${TEMP_DIR}/list_outgroup.txt")"
echo "Ingroup above ${SAMPLE_MISS} missing:      $(wc -l < "${TEMP_DIR}/list_failing.txt")"

if [[ -s "${TEMP_DIR}/list_failing.txt" ]]; then
    echo "Samples that fail the missingness threshold:"
    while read -r S; do
        awk -F'\t' -v s="${S}" '$1 == s {printf "  %-14s %.4f\n", $1, $4}' "${SMISS}"
    done < "${TEMP_DIR}/list_failing.txt"
fi

### exclusion lists, one per sample set (named after the set suffix)
cat "${TEMP_DIR}/list_outgroup.txt" "${TEMP_DIR}/list_failing.txt" \
    > "${TEMP_DIR}/exclude_ingroup.txt"
cp  "${TEMP_DIR}/list_failing.txt" \
    "${TEMP_DIR}/exclude_plusoutgroup.txt"


### =====================================================================
### variant count table
### =====================================================================

### filled as each VCF is written
COUNTS="${RESULTS}/variant_counts.tsv"
printf "file\treference\tsample_set\tmaf_level\tn_variants\tn_samples\n" > "${COUNTS}"

record_count() {
    local VCF=$1
    local SET=$2
    local MAF_LEVEL=$3
    printf "%s\tconsensus\t%s\t%s\t%d\t%d\n" \
        "$(basename "${VCF}" .vcf.gz)" "${SET}" "${MAF_LEVEL}" \
        "$(bcftools index -n "${VCF}")" \
        "$(bcftools query -l "${VCF}" | wc -l)" >> "${COUNTS}"
}


### =====================================================================
### function: per-variant metrics for one VCF
### =====================================================================

### full table, written through gzip, row count from the index
extract_metrics() {

    local VCF=$1
    local OUT_PREFIX=$2
    local METRICS="${RESULTS}/${OUT_PREFIX}_metrics.txt.gz"

    {
        printf "CHROM\tPOS\tQUAL\tDP\tMQ\tMAF\tF_MISSING\tNS\tAC\tAN\tAF\tMQ0F\tSCR\tSGB\tFILE\n"

        bcftools query \
            -f '%CHROM\t%POS\t%QUAL\t%INFO/DP\t%INFO/MQ\t%INFO/MAF\t%INFO/F_MISSING\t%INFO/NS\t%INFO/AC\t%INFO/AN\t%INFO/AF\t%INFO/MQ0F\t%INFO/SCR\t%INFO/SGB\n' \
            "${VCF}" \
        | awk -v file="${OUT_PREFIX}" 'BEGIN{OFS="\t"} {print $0, file}'

    } | gzip > "${METRICS}"

    echo "  metrics: $(bcftools index -n "${VCF}") rows -> $(basename "${METRICS}")"
}


### =====================================================================
### function: one sample set and its three MAF levels
### =====================================================================

make_sample_set() {

    local SUFFIX=$1        ### _ingroup or _plusoutgroup
    local LABEL=$2

    local EXCLUDE="${TEMP_DIR}/exclude${SUFFIX}.txt"
    local WORK="${TEMP_DIR}/set${SUFFIX}"
    local OUT_NOMAF="${RESULTS}/consensus${SUFFIX}_noMAF.vcf.gz"
    local N_EXCLUDE N_SAMPLES N_BEFORE N_AFTER N

    echo ""
    echo "================================================================="
    echo "Sample set: ${LABEL}"
    echo "================================================================="

    ## check if the exclusion list exists (a missing list would keep every sample)
    if [[ ! -f "${EXCLUDE}" ]]; then
        echo "ERROR: exclusion list not found: ${EXCLUDE}"
        return 1
    fi

    N_EXCLUDE=$(wc -l < "${EXCLUDE}")

    ### step 6: remove the samples this set does not use
    if [[ "${N_EXCLUDE}" -gt 0 ]]; then
        bcftools view -S ^"${EXCLUDE}" --force-samples \
            --threads "${THREADS}" \
            "${BASE_VCF}" -Oz -o "${WORK}_samples.vcf.gz"
    else
        cp "${BASE_VCF}" "${WORK}_samples.vcf.gz"
    fi
    bcftools index --threads "${THREADS}" "${WORK}_samples.vcf.gz"

    N_SAMPLES=$(bcftools query -l "${WORK}_samples.vcf.gz" | wc -l)
    echo "Step 6 - sample selection:         ${N_SAMPLES} samples (removed ${N_EXCLUDE})"

    ### step 7: tags recalculated (sample removal changes allele counts and missingness)
    bcftools +fill-tags "${WORK}_samples.vcf.gz" \
        -Oz -o "${WORK}_tags.vcf.gz" \
        -- -t F_MISSING,MAF,NS,AC,AN,AF
    bcftools index --threads "${THREADS}" "${WORK}_tags.vcf.gz"
    echo "Step 7 - tags recalculated:        F_MISSING, MAF, NS, AC, AN, AF"

    ### step 8: site missingness (genotypes masked, sample set final)
    bcftools filter -i "F_MISSING <= ${SITE_MISS}" --threads "${THREADS}" \
        "${WORK}_tags.vcf.gz" -Oz -o "${OUT_NOMAF}"
    bcftools index --threads "${THREADS}" "${OUT_NOMAF}"

    N_BEFORE=$(bcftools index -n "${WORK}_tags.vcf.gz")
    N_AFTER=$(bcftools index -n "${OUT_NOMAF}")
    echo "Step 8 - site missing <= ${SITE_MISS}:     ${N_AFTER} variants (removed $((N_BEFORE - N_AFTER)))"

    record_count "${OUT_NOMAF}" "${SUFFIX#_}" "noMAF"
    extract_metrics "${OUT_NOMAF}" "consensus${SUFFIX}_noMAF"

    ### step 9: the two MAF levels, both from the noMAF set
    for LEVEL in "maf01:${MAF_V1}" "maf05:${MAF_V2}"; do
        local NAME="${LEVEL%%:*}"
        local VALUE="${LEVEL##*:}"
        local OUT="${RESULTS}/consensus${SUFFIX}_${NAME}.vcf.gz"

        bcftools filter -i "MAF >= ${VALUE}" --threads "${THREADS}" \
            "${OUT_NOMAF}" -Oz -o "${OUT}"
        bcftools index --threads "${THREADS}" "${OUT}"

        N=$(bcftools index -n "${OUT}")
        echo ""
        echo "MAF >= ${VALUE}:                       ${N} variants"

        record_count "${OUT}" "${SUFFIX#_}" "${NAME}"
        extract_metrics "${OUT}" "consensus${SUFFIX}_${NAME}"
    done

    ### remove the intermediates of this set
    rm -f "${WORK}_samples.vcf.gz" "${WORK}_samples.vcf.gz.csi" \
          "${WORK}_tags.vcf.gz" "${WORK}_tags.vcf.gz.csi"
}


### =====================================================================
### the two sample sets
### =====================================================================

make_sample_set "_ingroup" \
    "ingroup only, samples above the missingness threshold removed"

make_sample_set "_plusoutgroup" \
    "ingroup and outgroup, samples above the missingness threshold removed"


### =====================================================================
### clean up
### =====================================================================

rm -f "${TEMP_DIR}"/step*.vcf.gz "${TEMP_DIR}"/step*.vcf.gz.csi
rm -rf "${TEMP_DIR}"


### output files (copy consensus_ingroup_noMAF_metrics.txt.gz, variant_counts.tsv, consensus_sample_missingness.tsv and sco_presence.tsv to local/01_variant_calling/data/080_data_filtering for local/01_variant_calling/080_filtering.R)
echo ""
echo "########## Output VCFs ##########"
echo ""
printf "%-46s %12s %8s\n" "file" "variants" "samples"
printf "%-46s %12s %8s\n" "----------------------------------------------" "------------" "--------"
awk -F'\t' 'NR>1 {printf "%-46s %12d %8d\n", $1 ".vcf.gz", $5, $6}' "${COUNTS}"

echo ""
echo "$(( $(wc -l < "${COUNTS}") - 1 )) VCFs written (expected 6)"

echo ""
echo "Tables for the R script:"
echo "  $(basename "${COUNTS}")"
echo "  $(basename "${SMISS}")"
echo "  $(basename "${SCO_TABLE}")"
echo "  consensus_*_metrics.txt.gz"

echo ""
echo "The site-filtered base is kept at:"
echo "  ${BASE_VCF}"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
HOURS=$((RUNTIME / 3600))
MINUTES=$(((RUNTIME % 3600) / 60))
echo ""
echo "Filtering completed in ${RUNTIME} seconds (${HOURS}h ${MINUTES}m)."
