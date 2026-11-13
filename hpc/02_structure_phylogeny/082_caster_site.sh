#!/bin/bash
#SBATCH --job-name=caster_site    ### job name
#SBATCH --cpus-per-task=20        ### cores
#SBATCH --mem-per-cpu=2G          ### memory per core
#SBATCH --time=3-00:00:00         ### walltime
#SBATCH --output /dev/null        ### stdout (logged by tee)
#SBATCH --error /dev/null         ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/02_structure_phylogeny/results/082_results_caster_site"
LOGS="${ROOT}/02_structure_phylogeny/logs/082_logs_caster_site"
### from hpc/02_structure_phylogeny/050_species_tree.sh (concatenation)
ALIGNMENT="${ROOT}/02_structure_phylogeny/results/050_results_species_tree/cenococcum_orthologs_concatenated.fasta"
### from hpc/02_structure_phylogeny/050_species_tree.sh (concatenation tree)
CONCAT_TREE="${ROOT}/02_structure_phylogeny/results/050_results_species_tree/species_tree.treefile"

TIPS="${RESULTS}/alignment_tips.txt"
TREE="${RESULTS}/caster_site_species_tree.tre"
UNROOTED="${RESULTS}/caster_site_species_tree_unrooted.tre"
SUPPORT="${RESULTS}/caster_site_support.tsv"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-20}"
### CASTER binary (ASTER package, same environment as 080/081)
CASTER=caster-site
### tip the output is rooted on (the R script roots on the whole outgroup)
ROOT_TIP=PsefloM402_1
### outgroup (samples by prefix, assemblies by name)
OUTGROUP_PREFIX="M4-"
OUTGROUP_GENOMES="PsefloM402_1,PsefloM405_1,PsefloM406_1"
### thorough search (-R = 16 rounds of placement and 16 of subsampling)
SEARCH="-R"
### seed
SEED=42
### support threshold (local bootstrap)
MIN_SUPPORT=95
export LC_ALL=C

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_caster_site.log"
RUNLOG="${LOGS}/${SLURM_JOB_ID:-local}_caster_site_run.log"
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

## check if the CASTER binary is there
if ! command -v "${CASTER}" > /dev/null 2>&1; then
    echo "ERROR: ${CASTER} not found in cg_bioinfo_v1_astral (bioconda package aster)"
    exit 1
fi
echo "Binary: $(command -v "${CASTER}")"
echo "CASTER-site $("${CASTER}" -h 2>&1 | grep -m 1 '^Version')"

## check if the alignment is there
if [[ ! -s "${ALIGNMENT}" ]]; then
    echo "ERROR: missing input ${ALIGNMENT}"
    exit 1
fi


### =====================================================================
### step 1: alignment of the 050 concatenation
### =====================================================================

echo ""
echo "=== Step 1: Alignment ==="

### tip names and columns
FORMAT=fasta
grep '^>' "${ALIGNMENT}" | sed 's/^>//; s/[[:space:]].*//' > "${TIPS}"
N_COLS=$(awk '/^>/ {if (n++) exit; next} {l += length($0)} END {print l}' "${ALIGNMENT}")

echo "Alignment: ${ALIGNMENT}"
echo "Format:    ${FORMAT}"
echo "Tips:      $(wc -l < "${TIPS}")"
echo "Columns:   ${N_COLS}"

## check if the root tip is in the alignment
if ! grep -qx "${ROOT_TIP}" "${TIPS}"; then
    echo "ERROR: root tip ${ROOT_TIP} not in the alignment"
    exit 1
fi

## check if the tips are unique
if [[ $(sort "${TIPS}" | uniq -d | wc -l) -gt 0 ]]; then
    echo "ERROR: duplicated tip names in the alignment"
    exit 1
fi


### =====================================================================
### step 2: species tree from site patterns
### =====================================================================

echo ""
echo "=== Step 2: CASTER-site ==="

### reused if newer than the alignment
if [[ -s "${TREE}" && "${TREE}" -nt "${ALIGNMENT}" ]]; then
    echo "CASTER-site tree found, reused."
else
    "${CASTER}" -t "${THREADS}" ${SEARCH} \
        --root "${ROOT_TIP}" \
        --seed "${SEED}" \
        -i "${ALIGNMENT}" \
        -o "${TREE}" 2> "${RUNLOG}"
    echo "Exit status: $?"
    echo "Run log: ${RUNLOG}"
    head -3 "${RUNLOG}"
    echo "..."
    tail -5 "${RUNLOG}"
fi

## check if CASTER wrote a tree
if [[ ! -s "${TREE}" ]]; then
    echo "ERROR: CASTER-site wrote no tree, see ${RUNLOG}"
    exit 1
fi


### =====================================================================
### step 3: tree checks and branch support
### =====================================================================

echo ""
echo "=== Step 3: Tree checks and support ==="

conda deactivate
conda activate cg_bioinfo_v1_R

Rscript - "${TREE}" "${TIPS}" "${CONCAT_TREE}" "${OUTGROUP_PREFIX}" "${OUTGROUP_GENOMES}" \
    "${MIN_SUPPORT}" "${UNROOTED}" "${SUPPORT}" <<'RSCRIPT'
suppressMessages(library(ape))

a <- as.list(commandArgs(trailingOnly = TRUE))
names(a) <- c("tree", "tips", "concat", "og_prefix", "og_genomes", "min_support",
              "out_unrooted", "out_support")

tr   <- read.tree(a$tree)
tips <- readLines(a$tips)

### outgroup tips
og <- tr$tip.label[startsWith(tr$tip.label, a$og_prefix) |
                   tr$tip.label %in% strsplit(a$og_genomes, ",")[[1]]]

## check if the tree carries exactly the alignment tips
if (!setequal(tr$tip.label, tips)) {
  warning("tree and alignment tips differ: ",
          length(setdiff(tips, tr$tip.label)), " missing, ",
          length(setdiff(tr$tip.label, tips)), " extra")
}

## check if the tree is fully resolved
if (!is.binary(tr)) warning("CASTER-site tree has polytomies")

## check if the outgroup is one clade (unrooted bipartition)
if (!is.monophyletic(tr, og, reroot = TRUE)) warning("outgroup is not monophyletic")

### unrooted topology with support (no branch lengths: terminal lengths are undefined)
unr <- unroot(tr)
unr$edge.length <- NULL
write.tree(unr, a$out_unrooted)

### local bootstrap per internal branch (root node has none)
n_tip <- Ntip(tr)
sup   <- suppressWarnings(as.numeric(tr$node.label))
desc  <- vapply(seq_len(tr$Nnode) + n_tip, function(nd) length(extract.clade(tr, nd)$tip.label),
                integer(1))
tab   <- data.frame(node = seq_len(tr$Nnode) + n_tip, tips_below = desc, support = sup)
tab   <- tab[!is.na(tab$support), ]
write.table(tab, a$out_support, sep = "\t", quote = FALSE, row.names = FALSE)

### distance to the concatenation topology (Robinson-Foulds, unrooted)
ct     <- read.tree(a$concat)
shared <- intersect(tr$tip.label, ct$tip.label)
t1     <- unroot(keep.tip(tr, shared))
t2     <- unroot(keep.tip(ct, shared))
rf     <- dist.topo(t1, t2)
rf_max <- 2 * (length(shared) - 3)
cat(sprintf("%-28s %d of %d (%.3f)\n", "RF distance to 050 tree", as.integer(rf), rf_max, rf / rf_max))

min_s <- as.numeric(a$min_support)
cat(sprintf("%-28s %d\n", "tips", n_tip))
cat(sprintf("%-28s %d\n", "outgroup tips", length(og)))
cat(sprintf("%-28s %d\n", "branches with support", nrow(tab)))
if (nrow(tab) > 0) {
  cat(sprintf("%-28s %.1f [%.1f, %.1f]\n", "support median [IQR]", median(tab$support),
              quantile(tab$support, 0.25), quantile(tab$support, 0.75)))
  cat(sprintf("%-28s %d\n", paste("support >=", min_s), sum(tab$support >= min_s)))
  cat(sprintf("%-28s %d\n", paste("support <", min_s), sum(tab$support < min_s)))
} else {
  cat("no numeric support labels on the tree\n")
}
RSCRIPT

## check if the unrooted tree was written
if [[ ! -s "${UNROOTED}" ]]; then
    echo "ERROR: no unrooted tree written"
    exit 1
fi

### output files (copy caster_site_species_tree.tre to local/02_structure_phylogeny/data/082_data_caster_site for local/02_structure_phylogeny/082_caster_site.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  caster_site_species_tree.tre            species tree, local bootstrap per branch, rooted on ${ROOT_TIP}"
echo "  caster_site_species_tree_unrooted.tre   the same topology unrooted, support kept (input to 090 and 091)"
echo "  caster_site_support.tsv                 local bootstrap and tips below per internal branch"
echo "  alignment_tips.txt                      tips of the 050 alignment"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "CASTER-site completed in ${RUNTIME} seconds ($((RUNTIME / 3600))h $(((RUNTIME % 3600) / 60))m)."
