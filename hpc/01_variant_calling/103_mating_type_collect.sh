#!/bin/bash
#SBATCH --job-name=mating_type_collect   ### job name
#SBATCH --cpus-per-task=4                ### cores
#SBATCH --mem-per-cpu=4G                 ### memory per core
#SBATCH --time=04:00:00                  ### walltime
#SBATCH --output /dev/null               ### stdout (logged by tee)
#SBATCH --error /dev/null                ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/103_results_mating_type_collect"
LOGS="${ROOT}/01_variant_calling/logs/103_logs_mating_type_collect"

### from hpc/01_variant_calling/000_prep.sh
SUPP="${ROOT}/01_variant_calling/supplementary/010_supplementary_read_trimming"

### idiomorph reference, core and primer sites, from hpc/01_variant_calling/101_mat_idiomorphs.sh
REF="${ROOT}/01_variant_calling/results/101_results_mat_idiomorphs/mat_idiomorphs.fasta"
CORE="${ROOT}/01_variant_calling/results/101_results_mat_idiomorphs/mat_core.bed"
PRIMERS_BED="${ROOT}/01_variant_calling/results/101_results_mat_idiomorphs/mat_primer_sites.bed"

### per-sample calls and alignments, from hpc/01_variant_calling/102_mating_type.sh
PER_SAMPLE="${ROOT}/01_variant_calling/results/102_results_mating_type/per_sample"

### output files
BAM_LIST="${RESULTS}/bam_list.txt"
VCF="${RESULTS}/mat_joint.vcf.gz"

SUMMARY="${RESULTS}/mating_type.tsv"
VARIANTS="${RESULTS}/mat_variants.tsv"
CALLABLE="${RESULTS}/mat_callable.tsv"
PRIMER_GT="${RESULTS}/mat_primer_genotypes.tsv"

### parameters
### minimum mapping quality
MAPQ=20

### minimum depth of a callable site
DEPTH_MIN=3

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### define list(s)
### trimmed read names, from hpc/01_variant_calling/000_prep.sh
SAMPLE_LIST="${SUPP}/ceno_genomics_igat_trimmed_reads_list_R1.txt"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_mating_type_collect.log"
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

## check if the reference, core and primer sites exist
if [[ ! -s "${REF}" || ! -s "${CORE}" || ! -s "${PRIMERS_BED}" ]]; then
    echo "ERROR: reference, core or primer sites missing, run 101_mat_idiomorphs.sh first"
    exit 1
fi


### =====================================================================
### step 1: calls
### =====================================================================

echo ""
echo "=== Step 1: calls ==="

N_EXPECTED=$(wc -l < "${SAMPLE_LIST}")
N_DONE=$(ls "${PER_SAMPLE}" | grep -c '_mating_type.tsv$')

## check if every sample has a call
if [[ "${N_DONE}" -ne "${N_EXPECTED}" ]]; then
    echo "WARNING: ${N_DONE} of ${N_EXPECTED} samples called, missing:"
    comm -23 <(xargs -n1 basename < "${SAMPLE_LIST}" | sed 's/_R1$//' | sort) \
             <(ls "${PER_SAMPLE}" | grep '_mating_type.tsv$' | sed 's/_mating_type.tsv$//' | sort)
fi

### calls of all samples
{
    printf "sample\tdepth_mat1\tbreadth_mat1\tdepth_mat2\tbreadth_mat2\tminor_fraction\tmat_called\n"
    cat "${PER_SAMPLE}"/*_mating_type.tsv | sort
} > "${SUMMARY}"

echo "calls:"
tail -n +2 "${SUMMARY}" | cut -f7 | sort | uniq -c


### =====================================================================
### step 2: joint call
### =====================================================================

echo ""
echo "=== Step 2: joint call ==="

conda deactivate
conda activate cg_bioinfo_v1_bcftools
echo "Running in environment: $CONDA_DEFAULT_ENV"

ls "${PER_SAMPLE}"/*.bam > "${BAM_LIST}"

### all samples over the core (invariant sites and reads without a proper pair kept)
bcftools mpileup -f "${REF}" -b "${BAM_LIST}" -R "${CORE}" -a AD,DP \
    -q "${MAPQ}" -d 2000 -A -Ou 2>/dev/null \
  | bcftools call -m --ploidy 1 -Oz -o "${VCF}" 2>/dev/null

bcftools index -f "${VCF}"

echo "Samples:        $(bcftools query -l "${VCF}" | wc -l)"
echo "Sites:          $(bcftools view -H "${VCF}" | wc -l)"
echo "Variable sites: $(bcftools view -H -m2 "${VCF}" | wc -l)"


### =====================================================================
### step 3: tables
### =====================================================================

echo ""
echo "=== Step 3: tables ==="

### variable sites per sample with reads
{
    printf "sample\tidiomorph\tpos\tref\talt\tqual\tgt\tdp\tad\n"
    bcftools view -m2 -Ou "${VCF}" \
      | bcftools query -i 'FMT/DP > 0' \
          -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%DP\t%AD\n]'
} > "${VARIANTS}"

### callable sites per sample and idiomorph
{
    printf "sample\tidiomorph\tn_callable\n"
    bcftools query -i "FMT/DP >= ${DEPTH_MIN}" -f '[%SAMPLE\t%CHROM\n]' "${VCF}" \
      | sort | uniq -c | sed -E 's/^ *([0-9]+) (.*)$/\2\t\1/'
} > "${CALLABLE}"

### genotypes at the primer sites
{
    printf "sample\tidiomorph\tpos\tref\talt\tgt\tdp\n"
    bcftools query -R "${PRIMERS_BED}" \
        -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%GT\t%DP\n]' "${VCF}"
} > "${PRIMER_GT}"

echo "Variant rows:  $(tail -n +2 "${VARIANTS}" | wc -l)"
echo "Callable rows: $(tail -n +2 "${CALLABLE}" | wc -l)"
echo "Primer rows:   $(tail -n +2 "${PRIMER_GT}" | wc -l)"


### output files (copy mating_type.tsv, mat_variants.tsv, mat_callable.tsv and mat_primer_genotypes.tsv to local/01_variant_calling/data/103_data_mating_type_collect for local/01_variant_calling/103_mating_type.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  bam_list.txt"
echo "  mat_joint.vcf.gz"
echo "  mating_type.tsv"
echo "  mat_variants.tsv"
echo "  mat_callable.tsv"
echo "  mat_primer_genotypes.tsv"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Mating type calls collected in ${RUNTIME} seconds."
