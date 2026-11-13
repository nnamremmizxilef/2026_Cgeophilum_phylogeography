#!/bin/bash
#SBATCH --job-name=aster      ### job name
#SBATCH --cpus-per-task=8     ### cores
#SBATCH --mem-per-cpu=7G      ### memory per core
#SBATCH --time=2-00:00:00     ### walltime
#SBATCH --output /dev/null    ### stdout (logged by tee)
#SBATCH --error /dev/null     ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
VARCALL="${ROOT}/01_variant_calling"
RESULTS="${ROOT}/02_structure_phylogeny/results/081_results_aster"
LOGS="${ROOT}/02_structure_phylogeny/logs/081_logs_aster"
### from hpc/02_structure_phylogeny/070_gene_trees.sh
TREE_DIR="${ROOT}/02_structure_phylogeny/results/070_results_gene_trees/per_og_tree"
### from hpc/01_variant_calling/050_cds_extraction.sh
OG_LIST="${VARCALL}/results/050_results_sco_consensus/og_list.txt"

GENE_TREES="${RESULTS}/gene_trees.treefile"
ASTER_TREE="${RESULTS}/aster_species_tree.tre"
ASTER_ROOTED="${RESULTS}/aster_species_tree_rooted.tre"
ASTER_PLAIN="${RESULTS}/aster_species_tree_plain.tre"
POLY_TREE="${RESULTS}/aster_polytomy_test.tre"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-8}"
### ASTER binary
ASTER="astral4"
### outgroup
OUTGROUP="PsefloM402_1,PsefloM405_1,PsefloM406_1,M4-01,M4-02,M4-03,M4-04,M4-05,M4-06"
### java heap for the ASTRAL-III scoring step, below the allocation
export _JAVA_OPTIONS="-Xmx48g"
export LC_ALL=C

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_aster.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_astral

START_TIME=$(date +%s)

### run info
echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"


### =====================================================================
### step 1: collect the gene trees
### =====================================================================

echo ""
echo "=== Step 1: Gene trees ==="

cat "${TREE_DIR}"/*.treefile > "${GENE_TREES}"

N_EXPECTED=$(wc -l < "${OG_LIST}")
N_TREES=$(wc -l < "${GENE_TREES}")

echo "Expected:  ${N_EXPECTED}"
echo "Collected: ${N_TREES}"

## check if gene trees were collected
if [[ "${N_TREES}" -eq 0 ]]; then
    echo "ERROR: no gene trees collected from ${TREE_DIR}"
    exit 1
fi


### =====================================================================
### step 2: coalescent species tree
### =====================================================================

echo ""
echo "=== Step 2: ASTER ==="

### help in the log (documents the flags)
echo "Binary: ${ASTER}"
echo "--- help ---"
"${ASTER}" --help 2>&1 | head -40
echo "--- end of help ---"

### quartet criterion of ASTRAL-III (-u 2: q1, q2, q3 and the posterior per branch; reused on resubmission)
if [[ -s "${ASTER_TREE}" ]]; then
    echo "Species tree found, reused."
else
    "${ASTER}" -i "${GENE_TREES}" -o "${ASTER_TREE}" -t "${THREADS}" -u 2
fi

## check if ASTER wrote a tree
if [[ ! -s "${ASTER_TREE}" ]]; then
    echo "ERROR: ASTER wrote no tree"
    exit 1
fi

echo "Tips in the tree: $(( $(tr -cd ',' < "${ASTER_TREE}" | wc -c) + 1 ))"


### =====================================================================
### step 3: polytomy test
### =====================================================================

echo ""
echo "=== Step 3: Polytomy test ==="

### ASTRAL-III jar of the active environment (the astral wrapper is shadowed by ASTER)
JAR=$(ls "${CONDA_PREFIX}"/share/astral-tree-*/astral.*.jar 2>/dev/null | head -1)

### topology only (annotation blocks in the node labels stripped)
conda deactivate
conda activate cg_bioinfo_v1_R

Rscript -e "
library(ape)
tr <- read.tree('${ASTER_TREE}')
tr\$node.label <- NULL
write.tree(tr, '${ASTER_PLAIN}')
cat('tips:', length(tr\$tip.label), '\n')
"

conda deactivate
conda activate cg_bioinfo_v1_astral

echo "jar: ${JAR}"

### polytomy test of ASTRAL-III on the fixed topology (not implemented in ASTER)
java -Xmx48g -jar "${JAR}" \
    -q "${ASTER_PLAIN}" \
    -i "${GENE_TREES}" \
    -o "${POLY_TREE}" \
    -t 10

## check if the polytomy test wrote a tree
if [[ ! -s "${POLY_TREE}" ]]; then
    echo "WARNING: polytomy test wrote no tree"
fi


### =====================================================================
### step 4: root
### =====================================================================

echo ""
echo "=== Step 4: Rooting ==="

conda deactivate
conda activate cg_bioinfo_v1_R

Rscript -e "
library(ape)
tr <- read.tree('${ASTER_TREE}')
og <- strsplit('${OUTGROUP}', ',')[[1]]
og <- intersect(og, tr\$tip.label)
if (length(og) == 0) stop('no outgroup tips in the tree')
cat('outgroup tips found:', length(og), '\n')
tr <- root(tr, outgroup = og, resolve.root = TRUE, edgelabel = TRUE)
write.tree(tr, '${ASTER_ROOTED}')
cat('tips:', length(tr\$tip.label), '\n')
"

### output files (copy aster_species_tree.tre and aster_polytomy_test.tre to local/02_structure_phylogeny/data/081_data_aster for local/02_structure_phylogeny/081_aster.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  gene_trees.treefile               the ${N_TREES} gene trees"
echo "  aster_species_tree.tre          coalescent tree from ASTER"
echo "  aster_species_tree_rooted.tre   the same tree rooted on the outgroup"
echo "  aster_species_tree_plain.tre    the same topology without annotations"
echo "  aster_polytomy_test.tre         per branch p value against a polytomy"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "ASTER completed in ${RUNTIME} seconds ($((RUNTIME / 3600))h $(((RUNTIME % 3600) / 60))m)."
