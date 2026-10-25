#!/bin/bash
#SBATCH --job-name=declone          ### job name
#SBATCH --cpus-per-task=4           ### cores
#SBATCH --mem-per-cpu=4G            ### memory per core
#SBATCH --time=04:00:00             ### walltime
#SBATCH --output /dev/null          ### stdout (logged by tee)
#SBATCH --error /dev/null           ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/091_results_decloned"
LOGS="${ROOT}/01_variant_calling/logs/091_logs_decloned"
### from hpc/01_variant_calling/080_filtering.sh
FILTERED_DIR="${ROOT}/01_variant_calling/results/080_results_filtering"
SUPP="${ROOT}/01_variant_calling/supplementary/091_supplementary_decloned"
### from local/01_variant_calling/090_clone_id.R (copy local/01_variant_calling/results/090_results_clone_id/per_sample_mlg.tsv to ${ROOT}/01_variant_calling/supplementary/091_supplementary_decloned)
MLG_TABLE="${SUPP}/per_sample_mlg.tsv"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_declone.log"
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
echo ""


### =====================================================================
### step 1: clonemates to drop
### =====================================================================

echo "=== Step 1: Clonemates ==="

### samples that are not the representative of their multilocus genotype (columns by name)
awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) col[$i] = i; next}
            $(col["is_representative"]) == "FALSE" {print $(col["sample"])}
' "${MLG_TABLE}" > "${RESULTS}/samples_excluded_clones.txt"

N_EXCLUDED=$(wc -l < "${RESULTS}/samples_excluded_clones.txt")

echo "MLG table:  ${MLG_TABLE}"
echo "Clonemates: ${N_EXCLUDED}"
echo ""


### =====================================================================
### step 2: drop them from every filtered VCF
### =====================================================================

echo "=== Step 2: Decloning ==="
echo ""

for VCF_FILE in "${FILTERED_DIR}"/*.vcf.gz; do

    BASENAME=$(basename "${VCF_FILE}" .vcf.gz)
    OUTPUT="${RESULTS}/${BASENAME}_decloned.vcf.gz"

    ### every sample except the clonemates (outgroup kept; -x whole-line matches, 1-42 must not match 1-422)
    bcftools query -l "${VCF_FILE}" > "${RESULTS}/tmp_samples.txt"
    grep -vxF -f "${RESULTS}/samples_excluded_clones.txt" \
        "${RESULTS}/tmp_samples.txt" > "${RESULTS}/tmp_keep.txt"

    N_IN=$(wc -l < "${RESULTS}/tmp_samples.txt")
    N_KEEP=$(wc -l < "${RESULTS}/tmp_keep.txt")

    ### sample columns only, sites unchanged (shared coordinates with the uncorrected set, monomorphic sites counted, not removed)
    ### fill-tags recalculates INFO for the retained samples
    bcftools view --threads "${THREADS}" -S "${RESULTS}/tmp_keep.txt" \
        "${VCF_FILE}" -Ou | \
    bcftools +fill-tags --threads "${THREADS}" -Oz -o "${OUTPUT}" \
        -- -t AC,AN,AF,MAF,F_MISSING,NS

    bcftools index --threads "${THREADS}" "${OUTPUT}"

    N_V_IN=$(bcftools index -n "${VCF_FILE}")
    N_V_OUT=$(bcftools index -n "${OUTPUT}")
    N_MONO=$(bcftools view -H -i 'COUNT(GT="alt")=0 || COUNT(GT="ref")=0' \
        "${OUTPUT}" | wc -l)

    echo "  ${BASENAME}"
    echo "    samples:  ${N_IN} -> ${N_KEEP} ($((N_IN - N_KEEP)) removed)"
    echo "    variants: ${N_V_IN} -> ${N_V_OUT} (must be equal)"
    echo "    of those ${N_MONO} are monomorphic in the retained samples"
    echo ""

done

rm -f "${RESULTS}/tmp_samples.txt" "${RESULTS}/tmp_keep.txt"

### output files
echo ""
echo "=== Output files in ${RESULTS} ==="
printf "%-55s  %10s  %8s\n" "File" "Variants" "Samples"
for F in "${RESULTS}"/*_decloned.vcf.gz; do
    printf "%-55s  %10d  %8d\n" "$(basename "${F}")" \
        "$(bcftools index -n "${F}")" "$(bcftools query -l "${F}" | wc -l)"
done
echo ""
echo "  samples_excluded_clones.txt  ${N_EXCLUDED} clonemates"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo "Completed in ${RUNTIME} seconds ($((RUNTIME / 60)) min)."
