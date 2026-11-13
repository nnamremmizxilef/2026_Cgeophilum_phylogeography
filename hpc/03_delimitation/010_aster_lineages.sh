#!/bin/bash
#SBATCH --job-name=aster_lineages    ### job name
#SBATCH --cpus-per-task=4            ### cores
#SBATCH --mem-per-cpu=8G             ### memory per core
#SBATCH --time=1-00:00:00            ### walltime
#SBATCH --output /dev/null           ### stdout (logged by tee)
#SBATCH --error /dev/null            ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
PHYLO="${ROOT}/02_structure_phylogeny"
RESULTS="${ROOT}/03_delimitation/results/010_results_aster_lineages"
LOGS="${ROOT}/03_delimitation/logs/010_logs_aster_lineages"
### from local/03_delimitation/001_lineage_labels.R (copy local/03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv and lineage_labels.tsv.md5 to ${ROOT}/03_delimitation/supplementary)
LABELS_DIR="${ROOT}/03_delimitation/supplementary"
### frozen labels (one representative per MLG: clustered = TRUE)
LABELS="${LABELS_DIR}/lineage_labels.tsv"
### from hpc/02_structure_phylogeny/080_aster_collapsed.sh (gene trees, internal branches below 1e-5 collapsed)
TREES_ALL="${PHYLO}/results/080_results_aster_collapsed/gene_trees_collapsed.treefile"

GENE_TREES="${RESULTS}/gene_trees_pruned.treefile"
INDIVIDUALS="${RESULTS}/individuals.tsv"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"
### ASTER binary (the one 080 and 081 used)
ASTER=astral4
### levels mapped to units, each a column of the labels table
K_LEVELS="k4 k5 k6"
### outgroup (samples by prefix, assemblies by name)
OUTGROUP_PREFIX="M4-"
OUTGROUP_GENOMES="PsefloM402_1,PsefloM405_1,PsefloM406_1"
### reference tip (not a taxon)
REF_TIP=consensus
### Cenococcum assemblies as units of their own (no: pruned)
KEEP_GENOMES=no
### fewest tips for a gene tree to be kept
MIN_TIPS=4
### java heap for the ASTRAL-III polytomy test
export _JAVA_OPTIONS="-Xmx28g"
export LC_ALL=C

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_aster_lineages.log"
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
for F in "${LABELS}" "${TREES_ALL}"; do
    if [[ ! -s "${F}" ]]; then
        echo "ERROR: input not found: ${F}"
        exit 1
    fi
done


### =====================================================================
### step 1: gene trees on one representative per MLG, mapping per level
### =====================================================================

echo ""
echo "=== Step 1: Gene trees and mapping ==="

echo "Gene trees from 080: $(wc -l < "${TREES_ALL}")"

### individuals and their unit at every level, gene trees pruned to them, mapping per level (ASTER: individual unit; ASTRAL-III: unit:individuals)
Rscript - "${LABELS}" "${TREES_ALL}" "${K_LEVELS}" "${OUTGROUP_PREFIX}" \
    "${OUTGROUP_GENOMES}" "${REF_TIP}" "${KEEP_GENOMES}" "${MIN_TIPS}" "${RESULTS}" <<'RSCRIPT'
suppressMessages(library(ape))

a <- as.list(commandArgs(trailingOnly = TRUE))
names(a) <- c("labels", "trees", "k_levels", "og_prefix", "og_genomes", "ref_tip",
              "keep_genomes", "min_tips", "out_dir")

k_levels <- strsplit(a$k_levels, " ")[[1]]

### labels
lab <- read.delim(a$labels, colClasses = "character")

### one representative per MLG
ind <- lab[lab$clustered == "TRUE", c("sample", k_levels)]
names(ind)[1] <- "tip"

### tips across all gene trees
trees <- read.tree(a$trees)
tips  <- unique(unlist(lapply(trees, `[[`, "tip.label")))

### outgroup samples and assemblies, Cenococcum assemblies
og    <- c(tips[startsWith(tips, a$og_prefix)], intersect(tips, strsplit(a$og_genomes, ",")[[1]]))
cenge <- setdiff(tips, c(lab$sample, og, a$ref_tip))

extra <- data.frame(tip = og, unit = "outgroup")
if (a$keep_genomes == "yes") extra <- rbind(extra, data.frame(tip = cenge, unit = cenge))
for (k in k_levels) extra[[k]] <- extra$unit

ind <- rbind(ind, extra[, c("tip", k_levels)])
ind <- ind[ind$tip %in% tips, ]

### gene trees reduced to the mapped tips
n_in <- vapply(trees, function(t) sum(t$tip.label %in% ind$tip), integer(1))
keep <- n_in >= as.integer(a$min_tips)

pruned <- lapply(trees[keep], function(t) keep.tip(t, intersect(t$tip.label, ind$tip)))
class(pruned) <- "multiPhylo"
write.tree(pruned, file.path(a$out_dir, "gene_trees_pruned.treefile"))

write.table(ind, file.path(a$out_dir, "individuals.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

### mapping per level
for (k in k_levels) {
  write.table(ind[, c("tip", k)], file.path(a$out_dir, sprintf("mapping_%s_aster.txt", k)),
              sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
  grp <- tapply(ind$tip, ind[[k]], paste, collapse = ",")
  writeLines(paste0(names(grp), ":", grp),
             file.path(a$out_dir, sprintf("mapping_%s_astral.txt", k)))
}

cat("representatives:             ", sum(lab$clustered == "TRUE"), "\n")
cat("outgroup tips:               ", length(og), "\n")
cat("Cenococcum assemblies:       ", length(cenge), if (a$keep_genomes == "yes") "(kept)" else "(pruned)", "\n")
cat("gene trees kept:             ", sum(keep), "of", length(trees), "\n")
cat("tips per gene tree (median): ", median(n_in[keep]), "\n")
for (k in k_levels) {
  cat("individuals per unit,", k, "\n")
  print(table(ind[[k]][!ind[[k]] %in% cenge]))
}
RSCRIPT

## check if the pruning wrote the gene trees and the individuals table
if [[ ! -s "${GENE_TREES}" || ! -s "${INDIVIDUALS}" ]]; then
    echo "ERROR: pruning wrote no gene trees or no individuals table"
    exit 1
fi


### =====================================================================
### step 2: cluster-level tree, polytomy test and rooting per level
### =====================================================================

for K in ${K_LEVELS}; do

    echo ""
    echo "=== Step 2: ${K} ==="

    MAP_ASTER="${RESULTS}/mapping_${K}_aster.txt"
    MAP_ASTRAL="${RESULTS}/mapping_${K}_astral.txt"
    TREE="${RESULTS}/aster_cluster_tree_${K}.tre"
    PLAIN="${RESULTS}/aster_cluster_tree_${K}_plain.tre"
    ROOTED="${RESULTS}/aster_cluster_tree_${K}_rooted.tre"
    POLY="${RESULTS}/aster_cluster_tree_${K}_polytomy.tre"
    UNITS_TREE="${RESULTS}/cluster_tree_${K}_rooted.tre"

    ### ASTER with individuals mapped to units (reused on resubmission)
    conda deactivate
    conda activate cg_bioinfo_v1_astral

    ## check if the ASTER binary is there
    if ! command -v "${ASTER}" > /dev/null 2>&1; then
        echo "ERROR: ${ASTER} not found in cg_bioinfo_v1_astral"
        exit 1
    fi

    if [[ -s "${TREE}" ]]; then
        echo "Cluster tree found, reused."
    else
        "${ASTER}" -i "${GENE_TREES}" -a "${MAP_ASTER}" -o "${TREE}" \
            -t "${THREADS}" -u 2
    fi

    ## check if ASTER wrote a tree
    if [[ ! -s "${TREE}" ]]; then
        echo "ERROR: ASTER wrote no tree for ${K}"
        exit 1
    fi

    ### topology without annotation, rooted tree, and the units of the level with the outgroup
    conda deactivate
    conda activate cg_bioinfo_v1_R

    Rscript - "${TREE}" "${LABELS}" "${K}" "${PLAIN}" "${ROOTED}" "${UNITS_TREE}" <<'RSCRIPT'
suppressMessages(library(ape))

a <- as.list(commandArgs(trailingOnly = TRUE))
names(a) <- c("tree", "labels", "k", "out_plain", "out_rooted", "out_units")

tr <- read.tree(a$tree)

plain <- tr
plain$node.label <- NULL
write.tree(plain, a$out_plain)

tr <- root(tr, outgroup = "outgroup", resolve.root = TRUE, edgelabel = TRUE)
write.tree(tr, a$out_rooted)

units <- c(unique(read.delim(a$labels, colClasses = "character")[[a$k]]), "outgroup")
u <- drop.tip(tr, setdiff(tr$tip.label, units))
u$node.label <- NULL
write.tree(u, a$out_units)

cat("units in the tree:", length(tr$tip.label), "\n")
cat("cluster tree:", write.tree(u), "\n")
RSCRIPT

    ### polytomy test on the fixed topology
    conda deactivate
    conda activate cg_bioinfo_v1_astral

    ### ASTRAL-III jar of the active environment
    JAR=$(ls "${CONDA_PREFIX}"/share/astral-tree-*/astral.*.jar 2>/dev/null | head -1)

    java -Xmx28g -jar "${JAR}" \
        -q "${PLAIN}" \
        -i "${GENE_TREES}" \
        -a "${MAP_ASTRAL}" \
        -o "${POLY}" \
        -t 10

    ## check if the polytomy test wrote a tree
    [[ -s "${POLY}" ]] || echo "WARNING: polytomy test wrote no tree for ${K}"

    ## check if every mapped unit is in the tree
    N_UNITS=$(cut -d' ' -f2 "${MAP_ASTER}" | sort -u | wc -l)
    N_TIPS=$(( $(tr -cd ',' < "${TREE}" | wc -c) + 1 ))
    echo "Units mapped: ${N_UNITS}   units in the tree: ${N_TIPS}"
    [[ "${N_UNITS}" -eq "${N_TIPS}" ]] || echo "WARNING: mapped units and tree tips differ for ${K}"

done

### output files (copy aster_cluster_tree_k*.tre, aster_cluster_tree_k*_polytomy.tre and individuals.tsv to local/03_delimitation/data/010_data_aster_lineages for local/03_delimitation/010_aster_lineages.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  individuals.tsv                        tip and its unit at every level"
echo "  gene_trees_pruned.treefile             gene trees on the mapped individuals"
echo "  mapping_kX_aster.txt, _astral.txt      mapping per level for ASTER and ASTRAL-III"
echo "  aster_cluster_tree_kX.tre              cluster tree, q1/q2/q3 and LPP per branch"
echo "  aster_cluster_tree_kX_plain.tre        the same topology without annotation"
echo "  aster_cluster_tree_kX_rooted.tre       rooted on the outgroup"
echo "  aster_cluster_tree_kX_polytomy.tre     per branch p value against a polytomy"
echo "  cluster_tree_kX_rooted.tre             the units of the level and the outgroup only"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Cluster trees completed in ${RUNTIME} seconds ($((RUNTIME / 3600))h $(((RUNTIME % 3600) / 60))m)."
