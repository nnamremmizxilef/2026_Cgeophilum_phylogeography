#!/bin/bash
#SBATCH --job-name=aster_collapsed    ### job name
#SBATCH --cpus-per-task=8             ### cores
#SBATCH --mem-per-cpu=7G              ### memory per core
#SBATCH --time=2-00:00:00             ### walltime
#SBATCH --output /dev/null            ### stdout (logged by tee)
#SBATCH --error /dev/null             ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
VARCALL="${ROOT}/01_variant_calling"
RESULTS="${ROOT}/02_structure_phylogeny/results/080_results_aster_collapsed"
LOGS="${ROOT}/02_structure_phylogeny/logs/080_logs_aster_collapsed"
### from hpc/02_structure_phylogeny/070_gene_trees.sh
TREE_DIR="${ROOT}/02_structure_phylogeny/results/070_results_gene_trees/per_og_tree"
### from hpc/01_variant_calling/050_cds_extraction.sh
OG_LIST="${VARCALL}/results/050_results_sco_consensus/og_list.txt"
### from hpc/02_structure_phylogeny/081_aster.sh (same gene trees, uncollapsed; run 081 before 080)
ASTER_RAW="${ROOT}/02_structure_phylogeny/results/081_results_aster/aster_species_tree.tre"
### from hpc/02_structure_phylogeny/081_aster.sh (polytomy test of the uncollapsed tree)
POLY_RAW="${ROOT}/02_structure_phylogeny/results/081_results_aster/aster_polytomy_test.tre"

GENE_TREES="${RESULTS}/gene_trees.treefile"
GENE_COLLAPSED="${RESULTS}/gene_trees_collapsed.treefile"
ASTER_TREE="${RESULTS}/aster_species_tree_collapsed.tre"
ASTER_PLAIN="${RESULTS}/aster_species_tree_collapsed_plain.tre"
POLY_TREE="${RESULTS}/aster_polytomy_test_collapsed.tre"
COMPARISON="${RESULTS}/aster_collapsed_comparison.tsv"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-8}"
### ASTER binary
ASTER="astral4"
### internal branches below this length collapsed (IQ-TREE minimum 1e-6, one substitution in a median gene about 7e-4)
BL_ZERO=1e-5
### polytomy test level
ALPHA=0.05
### java heap for the ASTRAL-III scoring step, below the allocation
export _JAVA_OPTIONS="-Xmx48g"
export LC_ALL=C

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_aster_collapsed.log"
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


### =====================================================================
### step 1: collapse the gene trees
### =====================================================================

echo ""
echo "=== Step 1: Collapse gene trees ==="

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

## check if the gene tree count matches the orthogroup list
if [[ "${N_TREES}" -ne "${N_EXPECTED}" ]]; then
    echo "WARNING: gene tree count differs from the orthogroup list"
fi

Rscript -e "
library(ape)
tr <- read.tree('${GENE_TREES}')
tr <- lapply(seq_along(tr), function(i) tr[[i]])
n_int <- function(t) sum(t\$edge[, 2] > length(t\$tip.label))
tc <- lapply(tr, di2multi, tol = ${BL_ZERO})
## check if collapsing kept every tip
if (!identical(lapply(tr, function(t) sort(t\$tip.label)),
               lapply(tc, function(t) sort(t\$tip.label)))) stop('tips changed')
before <- vapply(tr, n_int, numeric(1))
after  <- vapply(tc, n_int, numeric(1))
class(tc) <- 'multiPhylo'
write.tree(tc, '${GENE_COLLAPSED}')
cat('gene trees:', length(tc), '\n')
cat('internal branches before:', sum(before), '\n')
cat('internal branches after: ', sum(after), '\n')
cat(sprintf('collapsed: %d (%.1f%%)\n', sum(before - after), 100 * sum(before - after) / sum(before)))
cat('median internal branches per gene tree:', median(after), 'of', median(before), '\n')
"

## check if the collapsed file holds every gene tree
if [[ "$(wc -l < "${GENE_COLLAPSED}")" -ne "${N_TREES}" ]]; then
    echo "ERROR: collapsed file does not hold ${N_TREES} trees"
    exit 1
fi


### =====================================================================
### step 2: coalescent species tree
### =====================================================================

echo ""
echo "=== Step 2: ASTER ==="

conda deactivate
conda activate cg_bioinfo_v1_astral

echo "Binary: ${ASTER}"

### search reused on resubmission
if [[ -s "${ASTER_TREE}" ]]; then
    echo "Species tree found, reused."
else
    "${ASTER}" -i "${GENE_COLLAPSED}" -o "${ASTER_TREE}" -t "${THREADS}" -u 2
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

### ASTRAL-III jar of the active environment
JAR=$(ls "${CONDA_PREFIX}"/share/astral-tree-*/astral.*.jar 2>/dev/null | head -1)

### topology only
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

java -Xmx48g -jar "${JAR}" \
    -q "${ASTER_PLAIN}" \
    -i "${GENE_COLLAPSED}" \
    -o "${POLY_TREE}" \
    -t 10

## check if the polytomy test wrote a tree
if [[ ! -s "${POLY_TREE}" ]]; then
    echo "WARNING: polytomy test wrote no tree"
fi


### =====================================================================
### step 4: comparison with 081
### =====================================================================

echo ""
echo "=== Step 4: Comparison with the uncollapsed tree (081) ==="

conda deactivate
conda activate cg_bioinfo_v1_R

Rscript -e "
library(ape)

### ASTER fields by name
field <- function(x, key) {
  x[is.na(x)] <- ''
  pat <- paste0('(^|[^A-Za-z])', key, '=(-?[0-9.]+([eE][+-]?[0-9]+)?)')
  v <- rep(NA_real_, length(x))
  hit <- grepl(pat, x)
  v[hit] <- as.numeric(sub(paste0('.*', pat, '.*'), '\\\\2', x[hit]))
  v
}

### one row per internal branch (side without the first tip as key)
bip <- function(f) {
  tr <- read.tree(f)
  n <- length(tr\$tip.label)
  lab <- tr\$tip.label
  ref <- sort(lab)[1]
  nl <- if (is.null(tr\$node.label)) rep(NA_character_, tr\$Nnode) else tr\$node.label
  desc <- vector('list', n + tr\$Nnode)
  for (i in seq_len(n)) desc[[i]] <- i
  po <- reorder(tr, 'postorder')\$edge
  for (k in seq_len(nrow(po))) desc[[po[k, 1]]] <- c(desc[[po[k, 1]]], desc[[po[k, 2]]])
  rows <- lapply(seq_len(tr\$Nnode), function(i) {
    s <- lab[desc[[n + i]]]
    if (length(s) < 2 || length(s) > n - 2) return(NULL)
    if (ref %in% s) s <- setdiff(lab, s)
    data.frame(key = paste(sort(s), collapse = ','),
               n_side = min(length(s), n - length(s)), label = nl[i])
  })
  d <- do.call(rbind, rows)
  d <- d[order(is.na(d\$label) | d\$label == ''), ]
  d[!duplicated(d\$key), ]
}

raw <- bip('${ASTER_RAW}')
col <- bip('${ASTER_TREE}')

## check if both trees carry the same tips
if (!identical(sort(read.tree('${ASTER_RAW}')\$tip.label),
               sort(read.tree('${ASTER_TREE}')\$tip.label))) stop('tip sets differ')

raw <- data.frame(key = raw\$key, n_side = raw\$n_side, in_raw = TRUE,
                  pp1_raw = field(raw\$label, 'pp1'), q1_raw = field(raw\$label, 'q1'))
col <- data.frame(key = col\$key, n_side = col\$n_side, in_collapsed = TRUE,
                  pp1_collapsed = field(col\$label, 'pp1'), q1_collapsed = field(col\$label, 'q1'))
m <- merge(raw, col, by = c('key', 'n_side'), all = TRUE)
m\$in_raw[is.na(m\$in_raw)] <- FALSE
m\$in_collapsed[is.na(m\$in_collapsed)] <- FALSE

### polytomy test p values
add_p <- function(m, f, name) {
  p <- bip(f)
  p <- data.frame(key = p\$key, v = suppressWarnings(as.numeric(p\$label)))
  names(p)[2] <- name
  merge(m, p, by = 'key', all.x = TRUE)
}
m <- add_p(m, '${POLY_RAW}', 'p_raw')
m <- add_p(m, '${POLY_TREE}', 'p_collapsed')

m <- m[order(-m\$n_side), c('n_side', 'in_raw', 'in_collapsed', 'pp1_raw', 'pp1_collapsed',
                           'q1_raw', 'q1_collapsed', 'p_raw', 'p_collapsed', 'key')]
names(m)[names(m) == 'key'] <- 'tips'
write.table(m, '${COMPARISON}', sep = '\t', quote = FALSE, row.names = FALSE)

shared <- m\$in_raw & m\$in_collapsed
lost   <- m[m\$in_raw & !m\$in_collapsed, ]
hi     <- m\$in_raw & !is.na(m\$pp1_raw) & m\$pp1_raw >= 0.95
n_tip  <- length(read.tree('${ASTER_RAW}')\$tip.label)

cat('tips:', n_tip, '\n')
cat('internal branches, raw:', sum(m\$in_raw), ' collapsed:', sum(m\$in_collapsed),
    ' shared:', sum(shared), '\n')
cat(sprintf('normalised RF distance: %.4f\n',
            (sum(m\$in_raw & !m\$in_collapsed) + sum(!m\$in_raw & m\$in_collapsed)) /
            (2 * (n_tip - 3))))
cat('raw branches with pp1 >= 0.95:', sum(hi), ' kept:', sum(hi & m\$in_collapsed), '\n')
cat(sprintf('shared branches, median |delta pp1|: %.3f, |delta pp1| > 0.1: %d\n',
            median(abs(m\$pp1_raw[shared] - m\$pp1_collapsed[shared]), na.rm = TRUE),
            sum(abs(m\$pp1_raw[shared] - m\$pp1_collapsed[shared]) > 0.1, na.rm = TRUE)))
if (any(!is.na(m\$p_raw[shared])) && any(!is.na(m\$p_collapsed[shared]))) {
  cat('shared branches not separable from a polytomy (p > ${ALPHA}), raw:',
      sum(m\$p_raw[shared] > ${ALPHA}, na.rm = TRUE), ' collapsed:',
      sum(m\$p_collapsed[shared] > ${ALPHA}, na.rm = TRUE), '\n')
}
cat('\nlargest raw branches missing from the collapsed tree:\n')
print(head(lost[, c('n_side', 'pp1_raw', 'q1_raw', 'p_raw')], 15), row.names = FALSE)
"

### output files (copy aster_species_tree_collapsed.tre and aster_polytomy_test_collapsed.tre to local/02_structure_phylogeny/data/080_data_aster_collapsed for local/02_structure_phylogeny/080_aster_collapsed.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  gene_trees.treefile                      the ${N_TREES} gene trees"
echo "  gene_trees_collapsed.treefile            the same, branches < ${BL_ZERO} collapsed"
echo "  aster_species_tree_collapsed.tre         coalescent tree from ASTER"
echo "  aster_species_tree_collapsed_plain.tre   the same topology without annotations"
echo "  aster_polytomy_test_collapsed.tre        per branch p value against a polytomy"
echo "  aster_collapsed_comparison.tsv           per branch, 081 against 080"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "ASTER collapsed completed in ${RUNTIME} seconds ($((RUNTIME / 3600))h $(((RUNTIME % 3600) / 60))m)."
