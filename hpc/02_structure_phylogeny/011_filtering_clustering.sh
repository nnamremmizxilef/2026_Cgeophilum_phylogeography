#!/bin/bash
#SBATCH --job-name=filtering_clustering    ### job name
#SBATCH --cpus-per-task=4                  ### cores
#SBATCH --mem-per-cpu=8G                   ### memory per core
#SBATCH --time=08:00:00                    ### walltime
#SBATCH --output /dev/null                 ### stdout (logged by tee)
#SBATCH --error /dev/null                  ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
VARCALL="${ROOT}/01_variant_calling"
RESULTS="${ROOT}/02_structure_phylogeny/results/011_results_filtering_clustering"
REP_DIR="${RESULTS}/replicates"
LOGS="${ROOT}/02_structure_phylogeny/logs/011_logs_filtering_clustering"

### from hpc/02_structure_phylogeny/010_filtering_phylogeny.sh (masked call set)
VCF_IN="${ROOT}/02_structure_phylogeny/results/010_results_filtering_phylogeny/consensus_phylogeny.vcf.gz"

### from hpc/01_variant_calling/100_metadata.sh (read for in_ingroup, mlg and clone_status)
META="${VARCALL}/results/100_results_metadata/sample_metadata.tsv"

### output names (the 010 name, then one suffix per step)
BASE="consensus_phylogeny_ingroup_decloned_maf01"
POOL="${RESULTS}/${BASE}.vcf.gz"
OUT_VCF="${RESULTS}/${BASE}_1snp_per_sco.vcf.gz"
OUT_SITES="${RESULTS}/${BASE}_1snp_per_sco_selected_sites.tsv"
KEEP="${RESULTS}/clustering_samples.txt"
DRAWS="${RESULTS}/draws.tsv"

### parameters
### site missingness (on the retained samples)
F_MISSING=0.10
### minor allele frequency
MAF=0.01
### outgroup
OUTGROUP_PREFIX="M4-"
### seed of the primary draw (as in hpc/01_variant_calling/110_pruning_sco_consensus.sh)
SEED=42
### replicate draws (seeds 1 to N_REP)
N_REP=20
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"

### create results folders
mkdir -p "${RESULTS}" "${REP_DIR}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_filtering_clustering.log"
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

### parameters to the log
echo ""
echo "=== Filtering parameters ==="
echo "Samples:   ingroup, one per MLG (clone_status Unique or Representative)"
echo "Sites:     sites masked at 010 removed, biallelic SNPs,"
echo "           F_MISSING <= ${F_MISSING}, MAF >= ${MAF}"
echo "Thinning:  one random SNP per SCO, primary seed ${SEED},"
echo "           ${N_REP} replicate draws with seeds 1 to ${N_REP}"

## check if the inputs are there
if [[ ! -s "${VCF_IN}" ]]; then
    echo "ERROR: input VCF not found, run 010_filtering_phylogeny.sh first"
    echo "       ${VCF_IN}"
    exit 1
fi

if [[ ! -s "${META}" ]]; then
    echo "ERROR: sample metadata not found"
    echo "       ${META}"
    exit 1
fi

TEMP_DIR="${RESULTS}/temp_${SLURM_JOB_ID:-local}"
mkdir -p "${TEMP_DIR}"

N_INPUT=$(bcftools index -n "${VCF_IN}")

echo ""
echo "Input:   $(basename "${VCF_IN}")"
echo "Sites:   ${N_INPUT}"
echo "Samples: $(bcftools query -l "${VCF_IN}" | wc -l)"


### =====================================================================
### step 1: the sample set
### =====================================================================

echo ""
echo "########## Step 1: sample selection ##########"

### one sample per MLG (the representatives of every decloned analysis; carriage returns stripped)
gawk -F'\t' '
{ sub(/\r$/, "") }
NR == 1 {
    for (i = 1; i <= NF; i++) h[$i] = i
    split("sample in_ingroup mlg clone_status", need, " ")
    for (i in need)
        if (!(need[i] in h)) { print "MISSING_COLUMN " need[i] > "/dev/stderr"; bad = 1 }
    if (bad) exit 1
    next
}
$(h["in_ingroup"]) == "TRUE" &&
($(h["clone_status"]) == "Unique" || $(h["clone_status"]) == "Representative") {
    print $(h["sample"])
}' "${META}" | LC_ALL=C sort > "${TEMP_DIR}/meta_keep.txt"

## check if any sample was selected
if [[ ! -s "${TEMP_DIR}/meta_keep.txt" ]]; then
    echo "ERROR: no samples selected, check the column names in ${META}"
    exit 1
fi

### ingroup samples in the VCF (010 drops those above its missingness threshold)
bcftools query -l "${VCF_IN}" | grep -v "^${OUTGROUP_PREFIX}" | LC_ALL=C sort \
  > "${TEMP_DIR}/vcf_samples.txt"

LC_ALL=C comm -12 "${TEMP_DIR}/meta_keep.txt" "${TEMP_DIR}/vcf_samples.txt" > "${KEEP}"

N_META=$(wc -l < "${TEMP_DIR}/meta_keep.txt")
N_VCF=$(wc -l < "${TEMP_DIR}/vcf_samples.txt")
N_KEEP=$(wc -l < "${KEEP}")

echo "ingroup samples in the VCF:        ${N_VCF}"
echo "MLG representatives in the table:  ${N_META}"
echo "kept, present in both:             ${N_KEEP}"

## check if at least two samples are kept
if [[ "${N_KEEP}" -lt 2 ]]; then
    echo "ERROR: ${N_KEEP} samples kept, the two name sets do not match"
    exit 1
fi

### representatives the VCF does not hold
LC_ALL=C comm -23 "${TEMP_DIR}/meta_keep.txt" "${TEMP_DIR}/vcf_samples.txt" \
  > "${RESULTS}/representatives_absent_from_vcf.txt"

## check if every MLG representative is in the VCF
N_ABSENT=$(wc -l < "${RESULTS}/representatives_absent_from_vcf.txt")
[[ "${N_ABSENT}" -gt 0 ]] && \
  echo "WARNING: ${N_ABSENT} MLG representatives are not in the VCF, see representatives_absent_from_vcf.txt"


### =====================================================================
### step 2: site filters
### =====================================================================

echo ""
echo "########## Step 2: site filters ##########"

### sites masked at 010 (no call left for any sample)
bcftools view --threads "${THREADS}" -i 'NS > 0' \
    "${VCF_IN}" -Oz -o "${TEMP_DIR}/step1.vcf.gz"
bcftools index --threads "${THREADS}" "${TEMP_DIR}/step1.vcf.gz"
N_CALLED=$(bcftools index -n "${TEMP_DIR}/step1.vcf.gz")
echo "sites masked at 010, removed:      $(( N_INPUT - N_CALLED ))"
echo "sites with a call:                 ${N_CALLED}"

### biallelic SNPs on the kept samples
bcftools view --threads "${THREADS}" -S "${KEEP}" -m2 -M2 -v snps \
    "${TEMP_DIR}/step1.vcf.gz" -Oz -o "${TEMP_DIR}/step2.vcf.gz"
bcftools index --threads "${THREADS}" "${TEMP_DIR}/step2.vcf.gz"
echo "biallelic SNPs, kept samples:      $(bcftools index -n "${TEMP_DIR}/step2.vcf.gz")"

### tags on the kept samples
bcftools +fill-tags "${TEMP_DIR}/step2.vcf.gz" \
    -Oz -o "${TEMP_DIR}/step3.vcf.gz" -- -t F_MISSING,MAF,NS,AC,AN,AF
bcftools index --threads "${THREADS}" "${TEMP_DIR}/step3.vcf.gz"

### site missingness
bcftools view --threads "${THREADS}" -i "F_MISSING <= ${F_MISSING}" \
    "${TEMP_DIR}/step3.vcf.gz" -Oz -o "${TEMP_DIR}/step4.vcf.gz"
bcftools index --threads "${THREADS}" "${TEMP_DIR}/step4.vcf.gz"
echo "F_MISSING <= ${F_MISSING}:                 $(bcftools index -n "${TEMP_DIR}/step4.vcf.gz")"

### MAF (the pool every draw is taken from)
bcftools view --threads "${THREADS}" -i "MAF >= ${MAF}" \
    "${TEMP_DIR}/step4.vcf.gz" -Oz -o "${POOL}"
bcftools index --threads "${THREADS}" "${POOL}"
N_POOL=$(bcftools index -n "${POOL}")
echo "MAF >= ${MAF}:                       ${N_POOL}"


### =====================================================================
### step 3: one random SNP per SCO
### =====================================================================

echo ""
echo "########## Step 3: thinning ##########"

### positions of the pool (read once for every draw)
bcftools query -f '%CHROM\t%POS\n' "${POOL}" > "${TEMP_DIR}/pool_positions.tsv"

### function: draw (reservoir sampling as in hpc/01_variant_calling/110_pruning_sco_consensus.sh; gawk as mawk ignores the seed)
draw () {
    local seed="$1" vcf="$2" sites="$3"

    gawk -v seed="${seed}" 'BEGIN {srand(seed)}
    {
        if (!($1 in n)) order[++k] = $1
        n[$1]++
        if (rand() < 1 / n[$1]) pos[$1] = $2
    }
    END {for (i = 1; i <= k; i++) print order[i] "\t" pos[order[i]]}
    ' "${TEMP_DIR}/pool_positions.tsv" > "${sites}"

    bcftools view --threads "${THREADS}" -T "${sites}" "${POOL}" -Oz -o "${vcf}"
    bcftools index --threads "${THREADS}" "${vcf}"
}

### primary draw
draw "${SEED}" "${OUT_VCF}" "${OUT_SITES}"

N_OUT=$(bcftools index -n "${OUT_VCF}")
N_CONTIG=$(wc -l < "${OUT_SITES}")

echo "contigs with a variant:            ${N_CONTIG}"
echo "primary, seed ${SEED}:                ${N_OUT} sites"

LC_ALL=C sort "${OUT_SITES}" > "${TEMP_DIR}/primary_sorted.tsv"

printf "file\tdraw\tseed\tn_sites\tn_shared_with_primary\n" > "${DRAWS}"
printf "%s\tprimary\t%s\t%s\t%s\n" \
       "$(basename "${OUT_VCF}")" "${SEED}" "${N_OUT}" "${N_OUT}" >> "${DRAWS}"

### replicate draws
echo ""
for (( i = 1; i <= N_REP; i++ )); do

    R=$(printf "rep%02d" "${i}")
    V="${REP_DIR}/${BASE}_1snp_per_sco_${R}.vcf.gz"
    S="${REP_DIR}/${BASE}_1snp_per_sco_${R}_selected_sites.tsv"

    draw "${i}" "${V}" "${S}"

    N_S=$(bcftools index -n "${V}")
    N_SH=$(LC_ALL=C comm -12 "${TEMP_DIR}/primary_sorted.tsv" <(LC_ALL=C sort "${S}") | wc -l)

    printf "%s\t%s\t%s\t%s\t%s\n" "$(basename "${V}")" "${R}" "${i}" "${N_S}" "${N_SH}" >> "${DRAWS}"
    echo "  ${R}  seed ${i}  ${N_S} sites, ${N_SH} shared with the primary"
done

### the metadata beside the VCFs
cp "${META}" "${RESULTS}/sample_metadata.tsv"


### =====================================================================
### checks
### =====================================================================

echo ""
echo "########## Checks ##########"

## check if the primary draw kept exactly one site per contig
if [[ "${N_OUT}" -ne "${N_CONTIG}" ]]; then
    echo "ERROR: ${N_CONTIG} contigs drawn, ${N_OUT} sites in the primary draw"
    exit 1
fi

## check if every replicate holds the same number of sites as the primary
N_DIFF=$(awk -F'\t' -v n="${N_OUT}" 'NR > 1 && $4 != n' "${DRAWS}" | wc -l)
if [[ "${N_DIFF}" -gt 0 ]]; then
    echo "ERROR: ${N_DIFF} draws differ in site count from the primary, see ${DRAWS}"
    exit 1
fi

### mean overlap of the replicates with the primary
awk -F'\t' 'NR > 2 {s += $5 / $4; n++} END {printf "mean overlap with the primary:     %.2f%% of sites\n", 100 * s / n}' \
    "${DRAWS}"

## check if any site in the pool is invariant among the retained samples
N_INVAR=$(bcftools view -H -i 'AC == 0 || AC == AN' "${POOL}" | wc -l)
if [[ "${N_INVAR}" -gt 0 ]]; then
    echo "ERROR: ${N_INVAR} invariant sites survived the MAF filter"
    exit 1
fi

### outgroup samples in the output
echo "outgroup samples in the output:    $(bcftools query -l "${OUT_VCF}" | grep -c "^${OUTGROUP_PREFIX}")"

### worst site missingness in the pool
bcftools query -f '%INFO/F_MISSING\n' "${POOL}" \
| awk '{if ($1 + 0 > m) m = $1 + 0} END {printf "max site missingness:              %.4f\n", m}'

## check if no two retained samples share an MLG
gawk -F'\t' '
{ sub(/\r$/, "") }
NR == FNR { keep[$1] = 1; next }
FNR == 1 { for (i = 1; i <= NF; i++) h[$i] = i; next }
($(h["sample"]) in keep) {
    m = $(h["mlg"])
    if (m != "NA" && m != "") { seen[m]++; if (seen[m] > 1) dup++ }
}
END { printf "MLGs represented more than once:   %d\n", dup + 0; if (dup) exit 1 }
' "${KEEP}" "${META}" || { echo "ERROR: the retained set still holds clonemates"; exit 1; }


### =====================================================================
### clean up
### =====================================================================

rm -rf "${TEMP_DIR}"


### output files (copy consensus_phylogeny_ingroup_decloned_maf01_1snp_per_sco.vcf.gz and replicates/ to local/02_structure_phylogeny/data/014_data_clustering for local/02_structure_phylogeny/014_clustering.R and 015_clustering_publication.R)
echo ""
echo "########## Output ##########"
echo ""
echo "  pool:     ${N_POOL} sites"
echo "  draws:    ${N_OUT} sites each, primary plus ${N_REP} replicates"
echo "  samples:  $(bcftools query -l "${OUT_VCF}" | wc -l)"
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  ${BASE}.vcf.gz"
echo "      all SNPs passing the filters, the pool every draw comes from"
echo "  ${BASE}_1snp_per_sco.vcf.gz"
echo "      primary draw, seed ${SEED}"
echo "  replicates/${BASE}_1snp_per_sco_rep01 ... rep$(printf "%02d" "${N_REP}").vcf.gz"
echo "      replicate draws, seeds 1 to ${N_REP}"
echo "  *_selected_sites.tsv                the drawn sites of each draw"
echo "  draws.tsv                           seed, site count and overlap per draw"
echo "  sample_metadata.tsv                 copied from 100"
echo "  clustering_samples.txt              the retained samples"
echo "  representatives_absent_from_vcf.txt MLGs lost at 010"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Filtering for the clustering completed in ${RUNTIME} seconds ($((RUNTIME / 60)) min)."
