#!/bin/bash
#SBATCH --job-name=conversion_clustering    ### job name
#SBATCH --cpus-per-task=4                   ### cores
#SBATCH --mem-per-cpu=8G                    ### memory per core
#SBATCH --time=12:00:00                     ### walltime
#SBATCH --output /dev/null                  ### stdout (logged by tee)
#SBATCH --error /dev/null                   ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/02_structure_phylogeny/results/013_results_conversion_clustering"
LOGS="${ROOT}/02_structure_phylogeny/logs/013_logs_conversion_clustering"

### from hpc/02_structure_phylogeny/011_filtering_clustering.sh (SNP pool of the clustering, 412 MLG representatives, all SNPs passing 011)
POOL="${ROOT}/02_structure_phylogeny/results/011_results_filtering_clustering/consensus_phylogeny_ingroup_decloned_maf01.vcf.gz"

### output names
OUT="${RESULTS}/fineradstructure_input.txt"
KEY="${RESULTS}/sample_key.tsv"
LOCI="${RESULTS}/loci.tsv"

### parameters
### site missingness allowed (0 = sites called in every sample; RADpainter scores an unseen donor base as a near-certain difference)
SITE_MISS=0
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_conversion_clustering.log"
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

## check if the input is there
if [[ ! -s "${POOL}" ]]; then
    echo "ERROR: SNP pool not found, run 011_filtering_clustering.sh first"
    echo "       ${POOL}"
    exit 1
fi

TEMP_DIR="${RESULTS}/temp_${SLURM_JOB_ID:-local}"
mkdir -p "${TEMP_DIR}"

N_POOL=$(bcftools index -n "${POOL}")

echo ""
echo "Input:   $(basename "${POOL}")"
echo "Sites:   ${N_POOL}"
echo "Samples: $(bcftools query -l "${POOL}" | wc -l)"


### =====================================================================
### step 1: sample key
### =====================================================================

echo ""
echo "########## Step 1: sample key ##########"

### short ids in VCF order (sampleLD.R reads the header with read.table, which rewrites names such as 1-093)
bcftools query -l "${POOL}" \
  | gawk 'BEGIN {OFS = "\t"; print "id", "sample"} {printf "S%03d\t%s\n", NR, $1}' \
  > "${KEY}"

N_SAMPLES=$(( $(wc -l < "${KEY}") - 1 ))
echo "samples:                           ${N_SAMPLES}"


### =====================================================================
### step 2: sites
### =====================================================================

echo ""
echo "########## Step 2: sites ##########"

### missingness from the genotypes (bare F_MISSING, not the INFO tag)
bcftools view --threads "${THREADS}" -i "F_MISSING <= ${SITE_MISS}" \
    "${POOL}" -Oz -o "${TEMP_DIR}/sites.vcf.gz"
bcftools index --threads "${THREADS}" "${TEMP_DIR}/sites.vcf.gz"

N_SITES=$(bcftools index -n "${TEMP_DIR}/sites.vcf.gz")
echo "pool:                              ${N_POOL}"
echo "F_MISSING <= ${SITE_MISS}:                    ${N_SITES}"


### =====================================================================
### step 3: one haplotype per sample and SCO
### =====================================================================

echo ""
echo "########## Step 3: haplotypes ##########"

### header (sample ids only, no Chr column, so RADpainter reads a simple matrix)
gawk 'NR > 1 {printf "%s%s", (NR > 2 ? "\t" : ""), $1} END {print ""}' "${KEY}" > "${OUT}"

printf "row\tlocus\tn_snps\n" > "${LOCI}"

### bases of each sample at the SNPs of one SCO in position order (uncalled base N, no call at the SCO an empty cell)
bcftools query -f '%CHROM[\t%TGT]\n' "${TEMP_DIR}/sites.vcf.gz" \
| gawk -F'\t' -v loci="${LOCI}" '
function flush(   i, h, row) {
    row = ""
    for (i = 1; i <= n; i++) {
        h = hap[i]
        if (h ~ /^N+$/) h = ""
        row = row (i > 1 ? "\t" : "") h
        hap[i] = ""
    }
    print row
    rows++
    print rows "\t" cur "\t" nsnp >> loci
}
{
    if ($1 != cur) {
        if (cur != "") flush()
        if ($1 in seen) { print "ERROR: " $1 " is not contiguous" > "/dev/stderr"; failed = 1; exit 1 }
        seen[$1] = 1
        cur = $1
        nsnp = 0
        n = NF - 1
    }
    nsnp++
    for (i = 2; i <= NF; i++) {
        b = $i
        if (b !~ /^[ACGT]$/) b = "N"
        hap[i - 1] = hap[i - 1] b
    }
}
END { if (cur != "" && !failed) flush() }
' >> "${OUT}"

## check if the conversion ran through
if [[ "${PIPESTATUS[1]}" -ne 0 ]]; then
    echo "ERROR: conversion stopped, see the message above"
    exit 1
fi

N_LOCI=$(( $(wc -l < "${OUT}") - 1 ))
echo "SCOs with a site:                  ${N_LOCI}"

gawk -F'\t' 'NR > 1 {print $3}' "${LOCI}" | sort -n \
  | gawk '{v[NR] = $1; s += $1}
    END {printf "SNPs per SCO:                      min %d  median %d  mean %.1f  max %d\n",
         v[1], v[int((NR + 1) / 2)], s / NR, v[NR]}'


### =====================================================================
### checks
### =====================================================================

echo ""
echo "########## Checks ##########"

## check if every row holds one cell per sample
N_BAD=$(gawk -F'\t' -v n="${N_SAMPLES}" 'NF != n' "${OUT}" | wc -l)
if [[ "${N_BAD}" -gt 0 ]]; then
    echo "ERROR: ${N_BAD} rows without ${N_SAMPLES} cells"
    exit 1
fi

## check if every haplotype has the length of its SCO
N_LEN=$(gawk -F'\t' 'NR == FNR {if (FNR > 1) len[$1] = $3; next}
    FNR > 1 {for (i = 1; i <= NF; i++) if ($i != "" && length($i) != len[FNR - 1]) bad++}
    END {print bad + 0}' "${LOCI}" "${OUT}")
if [[ "${N_LEN}" -gt 0 ]]; then
    echo "ERROR: ${N_LEN} haplotypes differ in length from their SCO"
    exit 1
fi

## check if the SNP count matches the sites kept
N_SUM=$(gawk -F'\t' 'NR > 1 {s += $3} END {print s + 0}' "${LOCI}")
if [[ "${N_SUM}" -ne "${N_SITES}" ]]; then
    echo "ERROR: ${N_SUM} SNPs in the haplotypes, ${N_SITES} sites kept"
    exit 1
fi

### empty cells and uncalled bases (none at SITE_MISS = 0)
gawk -F'\t' 'FNR > 1 {for (i = 1; i <= NF; i++) {if ($i == "") e++; else n += gsub(/N/, "", $i)}}
    END {printf "empty cells:                       %d\nuncalled bases:                    %d\n", e + 0, n + 0}' "${OUT}"

## check if any sample id collides with the RADpainter format keywords
if head -1 "${OUT}" | cut -f1,3 | grep -qw "Chr"; then
    echo "ERROR: header read as Stacks or Matrix format"
    exit 1
fi


### =====================================================================
### clean up
### =====================================================================

rm -rf "${TEMP_DIR}"


### output files
echo ""
echo "########## Output ##########"
echo ""
echo "  samples:  ${N_SAMPLES}"
echo "  SCOs:     ${N_LOCI}"
echo "  SNPs:     ${N_SITES} of ${N_POOL} in the pool"
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  fineradstructure_input.txt  one row per SCO, one haplotype per sample"
echo "  sample_key.tsv              id in the header and the sample it stands for"
echo "  loci.tsv                    SCO and SNP count of each row"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Conversion completed in ${RUNTIME} seconds ($((RUNTIME / 60)) min)."
