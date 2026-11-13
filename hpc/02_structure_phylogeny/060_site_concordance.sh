#!/bin/bash
#SBATCH --job-name=site_concordance    ### job name
#SBATCH --cpus-per-task=32             ### cores
#SBATCH --mem-per-cpu=30G              ### memory per core
#SBATCH --time=7-00:00:00              ### walltime
#SBATCH --output /dev/null             ### stdout (logged by tee)
#SBATCH --error /dev/null              ### stderr (logged by tee)

### no core dumps
ulimit -c 0

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/02_structure_phylogeny/results/060_results_site_concordance"
LOGS="${ROOT}/02_structure_phylogeny/logs/060_logs_site_concordance"
### from hpc/02_structure_phylogeny/052_species_tree_decloned.sh (tree and decloned alignment)
TREE_DIR="${ROOT}/02_structure_phylogeny/results/052_results_species_tree_decloned"
TREEFILE="${TREE_DIR}/species_tree_decloned.treefile"
CONCATENATED="${TREE_DIR}/cenococcum_orthologs_concatenated_decloned.fasta"
PREFIX="${RESULTS}/site_concordance"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-32}"
### quartets drawn per branch
N_QUARTETS=100

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_site_concordance.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_iqtree

START_TIME=$(date +%s)

### run info
echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"
echo "iqtree version: $(iqtree2 --version | head -1)"

## check if the tree is there
if [[ ! -s "${TREEFILE}" ]]; then
    echo "ERROR: tree not found: ${TREEFILE}"
    exit 1
fi

## check if the alignment is there
if [[ ! -s "${CONCATENATED}" ]]; then
    echo "ERROR: alignment not found: ${CONCATENATED}"
    exit 1
fi


### =====================================================================
### step 1: site concordance factors
### =====================================================================

echo ""
echo "=== Site concordance ==="

### sCF per branch (decisive sites against the two alternative quartet resolutions)
iqtree2 -te "${TREEFILE}" \
    -s "${CONCATENATED}" \
    --scfl "${N_QUARTETS}" \
    -m GTR+G \
    -T "${THREADS}" \
    -mem 900G \
    --prefix "${PREFIX}" \
    -redo

## check if IQ-TREE wrote the concordance tree
if [[ ! -s "${PREFIX}.cf.tree" ]]; then
    echo "ERROR: no concordance tree written"
    exit 1
fi


### =====================================================================
### step 2: summary
### =====================================================================

echo ""
echo "=== Distribution of sCF ==="

### sCF distribution (floor about 33 %, three resolutions of a quartet: values near it mean no signal)
awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
    $(c["sCF"]) != "NA" {print $(c["sCF"])}' "${PREFIX}.cf.stat" \
| sort -g \
| awk '{v[NR] = $1; s += $1; if ($1 + 0 < 40) lo++}
    END {printf "branches   %d\nmin        %.1f\nmedian     %.1f\nmean       %.1f\nmax        %.1f\nbelow 40%%  %d\n",
         NR, v[1], v[int((NR + 1) / 2)], s / NR, v[NR], lo + 0}'

### output files (copy site_concordance.cf.stat, site_concordance.cf.tree and site_concordance.cf.branch to local/02_structure_phylogeny/data/060_data_site_concordance for local/02_structure_phylogeny/060_site_concordance.R and 100_cluster_monophylly.R, and site_concordance.cf.stat to local/02_structure_phylogeny/data/090_data_gene_concordance for local/02_structure_phylogeny/090_gene_concordance.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  site_concordance.cf.tree    tree with sCF on the branches"
echo "  site_concordance.cf.stat    sCF and decisive site counts per branch"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Site concordance completed in ${RUNTIME} seconds ($((RUNTIME / 3600))h $(((RUNTIME % 3600) / 60))m)."
