#!/bin/bash
#SBATCH --job-name=tree_comparison    ### job name
#SBATCH --cpus-per-task=32            ### cores
#SBATCH --mem-per-cpu=32G             ### memory per core
#SBATCH --time=3-00:00:00             ### walltime
#SBATCH --output /dev/null            ### stdout (logged by tee)
#SBATCH --error /dev/null             ### stderr (logged by tee)

### no core dumps
ulimit -c 0

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/02_structure_phylogeny/results/091_results_tree_comparison"
LOGS="${ROOT}/02_structure_phylogeny/logs/091_logs_tree_comparison"
### from hpc/02_structure_phylogeny/050_species_tree.sh (concatenation)
CONCATENATED="${ROOT}/02_structure_phylogeny/results/050_results_species_tree/cenococcum_orthologs_concatenated.fasta"
### the four topologies (all four needed for the test)
### from hpc/02_structure_phylogeny/052_species_tree_decloned.sh (concatenation tree on the MLG representatives)
CONCAT_TREE="${ROOT}/02_structure_phylogeny/results/052_results_species_tree_decloned/species_tree_decloned.treefile"
### from hpc/02_structure_phylogeny/080_aster_collapsed.sh
COLLAPSED_TREE="${ROOT}/02_structure_phylogeny/results/080_results_aster_collapsed/aster_species_tree_collapsed_plain.tre"
### from hpc/02_structure_phylogeny/081_aster.sh
ASTER_TREE="${ROOT}/02_structure_phylogeny/results/081_results_aster/aster_species_tree_plain.tre"
### from hpc/02_structure_phylogeny/082_caster_site.sh
CASTER_TREE="${ROOT}/02_structure_phylogeny/results/082_results_caster_site/caster_site_species_tree_unrooted.tre"

COMMON_TIPS="${RESULTS}/common_tips.txt"
ALN_SUB="${RESULTS}/concatenated_common_tips.fasta"
TREE_SET="${RESULTS}/candidate_topologies.treefile"
PREFIX="${RESULTS}/topology_test"
RF_PREFIX="${RESULTS}/topology_distances"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-32}"
### replicates for the resampling behind the AU test
N_RELL=10000

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_tree_comparison.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_R

### C locale
export LC_ALL=C

START_TIME=$(date +%s)

### run info
echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"

## check if every topology is there (a test over a subset answers a different question)
for T in "${CONCAT_TREE}" "${COLLAPSED_TREE}" "${ASTER_TREE}" "${CASTER_TREE}"; do
    if [[ ! -s "${T}" ]]; then
        echo "ERROR: topology not found: ${T}"
        exit 1
    fi
done

## check if the concatenation is there
if [[ ! -s "${CONCATENATED}" ]]; then
    echo "ERROR: concatenation not found: ${CONCATENATED}"
    exit 1
fi


### =====================================================================
### step 1: a common tip set
### =====================================================================

echo ""
echo "=== Step 1: Common tip set ==="

### all four trees pruned to their shared tips (clonemates in the coalescent trees, not in the 052 tree)
Rscript -e "
library(ape)

tr <- list(concat          = read.tree('${CONCAT_TREE}'),
           aster_collapsed = read.tree('${COLLAPSED_TREE}'),
           aster           = read.tree('${ASTER_TREE}'),
           caster          = read.tree('${CASTER_TREE}'))

for (n in names(tr)) cat(n, 'tips:', length(tr[[n]]\$tip.label), '\n')

common <- Reduce(intersect, lapply(tr, function(x) x\$tip.label))
cat('shared tips:', length(common), '\n')

if (length(common) < 4) stop('fewer than four shared tips')

writeLines(common, '${COMMON_TIPS}')

### unrooted topologies only (annotation and branch lengths dropped, iqtree2 refits the lengths)
pruned <- lapply(tr, function(x) {
  y <- unroot(drop.tip(x, setdiff(x\$tip.label, common)))
  y\$node.label  <- NULL
  y\$edge.length <- NULL
  y
})

### one topology per line, in the order concat, aster_collapsed, aster, caster
write.tree(pruned[['concat']],          '${TREE_SET}')
write.tree(pruned[['aster_collapsed']], '${TREE_SET}', append = TRUE)
write.tree(pruned[['aster']],           '${TREE_SET}', append = TRUE)
write.tree(pruned[['caster']],          '${TREE_SET}', append = TRUE)
"

## check if the pruning wrote topologies
if [[ ! -s "${TREE_SET}" ]]; then
    echo "ERROR: pruning produced no topologies"
    exit 1
fi

echo "Topologies written: $(wc -l < "${TREE_SET}")"
echo "Shared tips:        $(wc -l < "${COMMON_TIPS}")"


### =====================================================================
### step 2: the alignment on the same tips
### =====================================================================

echo ""
echo "=== Step 2: Alignment subset ==="

conda deactivate
conda activate cg_bioinfo_v1_seqkit

### alignment on the taxa of the tree set (subset reused on resubmission)
if [[ -s "${ALN_SUB}" ]]; then
    echo "Subset alignment found, reused."
else
    seqkit grep -n -f "${COMMON_TIPS}" "${CONCATENATED}" > "${ALN_SUB}"
fi

echo "Sequences: $(grep -c '^>' "${ALN_SUB}")"


### =====================================================================
### step 3: topological distances
### =====================================================================

echo ""
echo "=== Step 3: Robinson-Foulds distances ==="

conda deactivate
conda activate cg_bioinfo_v1_iqtree
echo "iqtree version: $(iqtree2 --version | head -1)"

### pairwise Robinson-Foulds distances between the four topologies (independent of the alignment)
iqtree2 -t "${TREE_SET}" --rf-all --prefix "${RF_PREFIX}" -redo

if [[ -s "${RF_PREFIX}.rfdist" ]]; then
    echo "--- rfdist ---"
    cat "${RF_PREFIX}.rfdist"
fi


### =====================================================================
### step 4: topology test
### =====================================================================

echo ""
echo "=== Step 4: AU test ==="

### AU test of the four topologies on the concatenated alignment (model as in 050: only the topology differs)
iqtree2 -s "${ALN_SUB}" \
    -m GTR+G \
    -z "${TREE_SET}" \
    -n 0 \
    -zb "${N_RELL}" \
    -zw \
    -au \
    -T "${THREADS}" \
    -mem 900G \
    --prefix "${PREFIX}" \
    -redo

## check if the topology test wrote a report
if [[ ! -s "${PREFIX}.iqtree" ]]; then
    echo "ERROR: the topology test wrote no report"
    exit 1
fi

echo ""
echo "--- tree topology test ---"
sed -n '/USER TREES/,/^$/p' "${PREFIX}.iqtree" | head -40

### output files (copy candidate_topologies.treefile, topology_distances.rfdist and topology_test.iqtree to local/02_structure_phylogeny/data/091_data_tree_comparison for local/02_structure_phylogeny/091_tree_comparison.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  common_tips.txt                  the tips all four topologies share"
echo "  candidate_topologies.treefile    concat, aster_collapsed, aster and caster, in that order"
echo "  concatenated_common_tips.fasta   the alignment on those tips"
echo "  topology_distances.rfdist        pairwise Robinson-Foulds distances"
echo "  topology_test.iqtree             AU, KH, SH and c-ELW per topology"
echo "  topology_test.sitelh             per site log likelihood per topology"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Tree comparison completed in ${RUNTIME} seconds ($((RUNTIME / 3600))h $(((RUNTIME % 3600) / 60))m)."
