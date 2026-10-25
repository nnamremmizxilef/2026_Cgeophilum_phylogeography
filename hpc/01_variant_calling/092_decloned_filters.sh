#!/bin/bash
#SBATCH --job-name=decloned_filters    ### job name
#SBATCH --cpus-per-task=4              ### cores
#SBATCH --mem-per-cpu=4G               ### memory per core
#SBATCH --time=04:00:00                ### walltime
#SBATCH --output /dev/null             ### stdout (logged by tee)
#SBATCH --error /dev/null              ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/092_results_decloned_filters"
LOGS="${ROOT}/01_variant_calling/logs/092_logs_decloned_filters"
### from hpc/01_variant_calling/091_decolned.sh
DECLONED_DIR="${ROOT}/01_variant_calling/results/091_results_decloned"
### decloned sets with all sites (noMAF sets: all three levels; site-filtered base: polymorphic level only)
INPUTS=(
    "consensus_ingroup_noMAF_decloned"
    "consensus_plusoutgroup_noMAF_decloned"
    "consensus_sitefiltered_base_decloned"
)

### parameters
### thresholds (same as hpc/01_variant_calling/080_filtering.sh)
MAF_V1=0.01
MAF_V2=0.05
SITE_MISS=0.10
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_decloned_filters.log"
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
echo "MAF levels: > 0, >= ${MAF_V1}, >= ${MAF_V2}"
echo "Site missingness: <= ${SITE_MISS}"
echo ""


### =====================================================================
### step 1: site filters on the decloned sets
### =====================================================================

echo "=== Filtering after clone correction ==="
echo ""

### MAF > 0 = polymorphic subset (drops sites without alternate or without reference carriers and all-missing sites)
### F_MISSING as recalculated in 091 for the retained samples; no fill-tags (removing sites leaves per-site INFO valid)
for BASENAME in "${INPUTS[@]}"; do

    INPUT="${DECLONED_DIR}/${BASENAME}.vcf.gz"

    ### noMAF sets: site missingness filter and three levels; base set: polymorphic level only
    if [[ "${BASENAME}" == *_noMAF_decloned ]]; then
        STEM="${BASENAME%_noMAF_decloned}"
        LEVELS=("polymorphic:MAF > 0" "maf01:MAF >= ${MAF_V1}" "maf05:MAF >= ${MAF_V2}")
        MISS_CLAUSE=" && F_MISSING <= ${SITE_MISS}"
    else
        STEM="${BASENAME%_decloned}"
        LEVELS=("polymorphic:MAF > 0")
        MISS_CLAUSE=""
    fi

    echo "  ${BASENAME}: $(bcftools index -n "${INPUT}") sites, $(bcftools query -l "${INPUT}" | wc -l) samples"

    for LEVEL in "${LEVELS[@]}"; do

        NAME="${LEVEL%%:*}"
        EXPR="${LEVEL#*:}${MISS_CLAUSE}"
        OUTPUT="${RESULTS}/${STEM}_decloned_${NAME}.vcf.gz"

        bcftools filter --threads "${THREADS}" -i "${EXPR}" \
            "${INPUT}" -Oz -o "${OUTPUT}"

        bcftools index --threads "${THREADS}" "${OUTPUT}"

        printf "    %-12s %-40s %10d sites (removed %d)\n" \
            "${NAME}" "${EXPR}" \
            "$(bcftools index -n "${OUTPUT}")" \
            "$(( $(bcftools index -n "${INPUT}") - $(bcftools index -n "${OUTPUT}") ))"

    done
    echo ""

done

### output files
echo ""
echo "=== Output files in ${RESULTS} ==="
printf "%-58s  %10s  %8s\n" "File" "Variants" "Samples"
for F in "${RESULTS}"/*.vcf.gz; do
    printf "%-58s  %10d  %8d\n" "$(basename "${F}")" \
        "$(bcftools index -n "${F}")" "$(bcftools query -l "${F}" | wc -l)"
done

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo "Completed in ${RUNTIME} seconds ($((RUNTIME / 60)) min)."
