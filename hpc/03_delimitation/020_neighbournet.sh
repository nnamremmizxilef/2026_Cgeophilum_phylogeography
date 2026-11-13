#!/bin/bash
#SBATCH --job-name=neighbournet    ### job name
#SBATCH --cpus-per-task=4          ### cores
#SBATCH --mem-per-cpu=4G           ### memory per core
#SBATCH --time=04:00:00            ### walltime
#SBATCH --output /dev/null         ### stdout (logged by tee)
#SBATCH --error /dev/null          ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
PHYLO="${ROOT}/02_structure_phylogeny"
RESULTS="${ROOT}/03_delimitation/results/020_results_neighbournet"
LOGS="${ROOT}/03_delimitation/logs/020_logs_neighbournet"
### from local/03_delimitation/001_lineage_labels.R (copy local/03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv and lineage_labels.tsv.md5 to ${ROOT}/03_delimitation/supplementary)
LABELS_DIR="${ROOT}/03_delimitation/supplementary"
### frozen labels (one representative per MLG: clustered = TRUE)
LABELS="${LABELS_DIR}/lineage_labels.tsv"
### from hpc/02_structure_phylogeny/011_filtering_clustering.sh (decloned SNP pool: biallelic, F_MISSING <= 0.10, MAF >= 0.01)
POOL="${PHYLO}/results/011_results_filtering_clustering/consensus_phylogeny_ingroup_decloned_maf01.vcf.gz"

KEEP="${RESULTS}/network_samples.txt"
PREFIX="${RESULTS}/network_ibs"
NET_PREFIX=neighbournet

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"
### split weights kept in the network (fast-nnt default)
CUTOFF=1e-4

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_neighbournet.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_R

START_TIME=$(date +%s)

### run info
echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"

## check if the labels are the frozen version
( cd "${LABELS_DIR}" && md5sum -c lineage_labels.tsv.md5 ) \
  || { echo "ERROR: lineage_labels.tsv does not match its md5"; exit 1; }

## check if the inputs are there
for F in "${LABELS}" "${POOL}"; do
    if [[ ! -s "${F}" ]]; then
        echo "ERROR: input not found: ${F}"
        exit 1
    fi
done


### =====================================================================
### step 1: samples
### =====================================================================

echo ""
echo "=== Step 1: Samples ==="

### one representative per MLG
Rscript -e "
lab  <- read.delim('${LABELS}', colClasses = 'character')
keep <- lab\$sample[lab\$clustered == 'TRUE']
writeLines(keep, '${KEEP}')
cat('representatives:', length(keep), '\n')
"

conda deactivate
conda activate cg_bioinfo_v1_bcftools
echo "bcftools version: $(bcftools --version | head -1)"

## check if the pool holds exactly the representatives
if ! diff -q <(sort "${KEEP}") <(bcftools query -l "${POOL}" | sort) > /dev/null; then
    echo "ERROR: samples in the pool and representatives in the labels differ"
    exit 1
fi

echo "Samples: $(bcftools query -l "${POOL}" | wc -l)   SNPs: $(bcftools index -n "${POOL}")"


### =====================================================================
### step 2: pairwise distances
### =====================================================================

echo ""
echo "=== Step 2: Distances ==="

conda deactivate
conda activate cg_bioinfo_v1_plink
echo "plink version: $(plink --version)"

### 1 - IBS, square matrix (haploid calls read as homozygous, so 1 - IBS is the share of compared SNPs at which two samples differ)
plink --vcf "${POOL}" \
    --allow-extra-chr \
    --double-id \
    --set-missing-var-ids @:# \
    --distance square 1-ibs flat-missing \
    --threads "${THREADS}" \
    --out "${PREFIX}"

## check if plink wrote the distance matrix
if [[ ! -s "${PREFIX}.mdist" ]]; then
    echo "ERROR: plink wrote no distance matrix"
    exit 1
fi


### =====================================================================
### checks
### =====================================================================

echo ""
echo "=== Checks ==="

## check if the matrix holds every kept sample
N_ROWS=$(wc -l < "${PREFIX}.mdist")
N_KEEP=$(wc -l < "${KEEP}")
echo "Matrix rows: ${N_ROWS}   samples kept: ${N_KEEP}"
[[ "${N_ROWS}" -eq "${N_KEEP}" ]] || echo "WARNING: matrix rows and kept samples differ"

## check if any distance could not be computed
N_NAN=$(grep -oi 'nan' "${PREFIX}.mdist" | wc -l)
echo "Missing distances: ${N_NAN}"
[[ "${N_NAN}" -eq 0 ]] || echo "WARNING: ${N_NAN} distances are nan"

### distance range
tr '\t' '\n' < "${PREFIX}.mdist" | sort -g | awk 'NR == 1 {min = $1} {max = $1} END {print "Distance range: " min " to " max}'


### =====================================================================
### step 3: Neighbour-Net
### =====================================================================

echo ""
echo "=== Step 3: Neighbour-Net ==="

conda deactivate
conda activate cg_bioinfo_v1_fastnnt
echo "Running in environment: $CONDA_DEFAULT_ENV"
echo "fast-nnt version: $(fast_nnt --version)"

### matrix as csv with sample names as header and first column
{
    printf ","
    cut -f2 "${PREFIX}.mdist.id" | paste -sd,
    paste -d, <(cut -f2 "${PREFIX}.mdist.id") <(tr '\t' ',' < "${PREFIX}.mdist")
} > "${PREFIX}.csv"

### splits by NNLS, equal-angle layout, written as nexus
fast_nnt neighbour_net \
    -i "${PREFIX}.csv" \
    -d "${RESULTS}" \
    -o "${NET_PREFIX}" \
    --cutoff "${CUTOFF}" \
    -t "${THREADS}"

## check if fast-nnt wrote the network
if [[ ! -s "${RESULTS}/${NET_PREFIX}.nex" ]]; then
    echo "ERROR: fast-nnt wrote no network"
    exit 1
fi

### least-squares fit and split count
echo "Least-squares fit: $(grep -o '"fit_percent": [0-9.]*' "${RESULTS}/run_log.json" | cut -d' ' -f2) %"
grep -m1 "nsplits" "${RESULTS}/${NET_PREFIX}.nex"

### output files (copy neighbournet.nex, run_log.json, network_ibs.mdist and network_ibs.mdist.id to local/03_delimitation/data/020_data_neighbournet for local/03_delimitation/020_neighbournet.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  network_samples.txt       one representative per MLG"
echo "  network_ibs.mdist         1 - IBS distance matrix (read by 020_neighbournet.R)"
echo "  network_ibs.mdist.id      sample order of the matrix"
echo "  network_ibs.csv           the same matrix as fast-nnt reads it"
echo "  neighbournet.nex          splits and equal-angle layout (also opens in SplitsTree)"
echo "  run_log.json              fast-nnt settings, split counts and least-squares fit"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Neighbour-Net completed in ${RUNTIME} seconds ($((RUNTIME / 60)) min)."
