#!/bin/bash
#SBATCH --job-name=pruning_sco_consensus    ### job name
#SBATCH --cpus-per-task=4                   ### cores
#SBATCH --mem-per-cpu=4G                    ### memory per core
#SBATCH --time=08:00:00                     ### walltime
#SBATCH --output /dev/null                  ### stdout (logged by tee)
#SBATCH --error /dev/null                   ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/110_results_pruning_sco_consensus"
LOGS="${ROOT}/01_variant_calling/logs/110_logs_pruning_sco_consensus"
### from hpc/01_variant_calling/080_filtering.sh, 091_decolned.sh and 092_decloned_filters.sh
INPUT_DIRS=(
    "${ROOT}/01_variant_calling/results/080_results_filtering"
    "${ROOT}/01_variant_calling/results/091_results_decloned"
    "${ROOT}/01_variant_calling/results/092_results_decloned_filters"
)

### parameters
### random seed (same seed and site list give the same draw)
SEED=42
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_pruning.log"
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
echo "gawk version: $(gawk --version | head -1)"
echo "Seed: ${SEED}"
echo ""


### =====================================================================
### step 1: one random site per contig and subset to it
### =====================================================================

echo "=== Step 1: One random SNP per SCO ==="
echo ""

for INPUT_DIR in "${INPUT_DIRS[@]}"; do
  for VCF_FILE in "${INPUT_DIR}"/*.vcf.gz; do

      BASENAME=$(basename "${VCF_FILE}" .vcf.gz)
      POSITIONS="${RESULTS}/${BASENAME}_selected_sites.tsv"
      OUTPUT="${RESULTS}/${BASENAME}_1snp_per_sco.vcf.gz"

      ### one draw per contig (n-th site replaces the held one with probability 1/n; contigs in VCF order, as bcftools view -T expects)
      ### gawk (mawk ignores the srand seed)
      bcftools query -f '%CHROM\t%POS\n' "${VCF_FILE}" | \
      gawk -v seed="${SEED}" 'BEGIN {srand(seed)}
      {
          if (!($1 in n)) order[++k] = $1
          n[$1]++
          if (rand() < 1 / n[$1]) pos[$1] = $2
      }
      END {for (i = 1; i <= k; i++) print order[i] "\t" pos[order[i]]}
      ' > "${POSITIONS}"

      bcftools view --threads "${THREADS}" -T "${POSITIONS}" \
          "${VCF_FILE}" -Oz -o "${OUTPUT}"

      bcftools index --threads "${THREADS}" "${OUTPUT}"

      N_V_IN=$(bcftools index -n "${VCF_FILE}")
      N_V_OUT=$(bcftools index -n "${OUTPUT}")
      N_SITES=$(wc -l < "${POSITIONS}")

      echo "  ${BASENAME}"
      echo "    contigs with a variant: ${N_SITES}"
      echo "    variants: ${N_V_IN} -> ${N_V_OUT} (must equal the contig count)"
      echo ""

  done
done

### output files
echo ""
echo "=== Output files in ${RESULTS} ==="
printf "%-70s  %8s  %8s\n" "File" "Variants" "Samples"
for F in "${RESULTS}"/*_1snp_per_sco.vcf.gz; do
    printf "%-70s  %8d  %8d\n" "$(basename "${F}")" \
        "$(bcftools index -n "${F}")" "$(bcftools query -l "${F}" | wc -l)"
done
echo ""
echo "  *_selected_sites.tsv  the drawn sites"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo "Completed in ${RUNTIME} seconds ($((RUNTIME / 60)) min)."
