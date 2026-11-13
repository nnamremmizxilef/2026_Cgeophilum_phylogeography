#!/bin/bash
#SBATCH --job-name=gene_concordance    ### job name
#SBATCH --cpus-per-task=4              ### cores
#SBATCH --mem-per-cpu=8G               ### memory per core
#SBATCH --time=04:00:00                ### walltime
#SBATCH --output /dev/null             ### stdout (logged by tee)
#SBATCH --error /dev/null              ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
VARCALL="${ROOT}/01_variant_calling"
RESULTS="${ROOT}/02_structure_phylogeny/results/090_results_gene_concordance"
LOGS="${ROOT}/02_structure_phylogeny/logs/090_logs_gene_concordance"
### from hpc/02_structure_phylogeny/070_gene_trees.sh (gene trees and IQ-TREE reports)
TREE_DIR="${ROOT}/02_structure_phylogeny/results/070_results_gene_trees/per_og_tree"
### from hpc/01_variant_calling/050_cds_extraction.sh
OG_LIST="${VARCALL}/results/050_results_sco_consensus/og_list.txt"
### from hpc/02_structure_phylogeny/080_aster_collapsed.sh (gene trees with internal branches < 1e-5 collapsed)
GENE_TREES_COLLAPSED="${ROOT}/02_structure_phylogeny/results/080_results_aster_collapsed/gene_trees_collapsed.treefile"
### from hpc/02_structure_phylogeny/052_species_tree_decloned.sh (concatenation tree on the MLG representatives)
CONCAT_TREE="${ROOT}/02_structure_phylogeny/results/052_results_species_tree_decloned/species_tree_decloned.treefile"
### from hpc/02_structure_phylogeny/080_aster_collapsed.sh
COLLAPSED_TREE="${ROOT}/02_structure_phylogeny/results/080_results_aster_collapsed/aster_species_tree_collapsed_plain.tre"
### from hpc/02_structure_phylogeny/081_aster.sh
ASTER_TREE="${ROOT}/02_structure_phylogeny/results/081_results_aster/aster_species_tree_plain.tre"
### from hpc/02_structure_phylogeny/082_caster_site.sh
CASTER_TREE="${ROOT}/02_structure_phylogeny/results/082_results_caster_site/caster_site_species_tree_unrooted.tre"

GENE_TREES="${RESULTS}/gene_trees.treefile"
GENE_TREES_DECLONED="${RESULTS}/gene_trees_decloned.treefile"
MISSING="${RESULTS}/missing_gene_trees.txt"
MODELS="${RESULTS}/gene_models.tsv"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_gene_concordance.log"
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

## check if the orthogroup list is there
if [[ ! -s "${OG_LIST}" ]]; then
    echo "ERROR: missing input ${OG_LIST}"
    exit 1
fi

## check if the collapsed gene trees are there
if [[ ! -s "${GENE_TREES_COLLAPSED}" ]]; then
    echo "ERROR: missing input ${GENE_TREES_COLLAPSED}"
    exit 1
fi

### reference trees
REF_NAMES=("concat" "aster_collapsed" "aster" "caster")
REF_PATHS=("${CONCAT_TREE}" "${COLLAPSED_TREE}" "${ASTER_TREE}" "${CASTER_TREE}")


### =====================================================================
### step 1: collect the gene trees
### =====================================================================

echo ""
echo "=== Step 1: Gene trees ==="

### the collapsed gene trees of 080 (same trees and order as 070)
cp "${GENE_TREES_COLLAPSED}" "${GENE_TREES}"

N_EXPECTED=$(wc -l < "${OG_LIST}")
N_TREES=$(wc -l < "${GENE_TREES}")

echo "Expected:  ${N_EXPECTED}"
echo "Collected: ${N_TREES}"

## check if the file holds gene trees
if [[ "${N_TREES}" -eq 0 ]]; then
    echo "ERROR: no gene trees in ${GENE_TREES_COLLAPSED}"
    exit 1
fi

### orthogroups without a gene tree
while read -r OG; do
    [[ -s "${TREE_DIR}/${OG}.treefile" ]] || echo "${OG}"
done < "${OG_LIST}" > "${MISSING}"

echo "Missing:   $(wc -l < "${MISSING}")"

## check if every orthogroup has a gene tree
if [[ -s "${MISSING}" ]]; then
    echo "WARNING: gCF is computed over an incomplete gene set, see ${MISSING}"
fi


### =====================================================================
### step 2: the gene tree panel
### =====================================================================

echo ""
echo "=== Step 2: Gene tree panel ==="

### tips per tree, one comma per tip after the first
echo "Tips per tree:"
awk '{print gsub(/,/, ",") + 1}' "${GENE_TREES}" | sort -n | uniq -c | sort -k2 -n \
    | awk '{printf "  %6d trees with %d tips\n", $1, $2}' | tail -5

### total branch length per tree (spread among the genes)
awk '{
        s = $0
        n = 0
        while (match(s, /:[0-9.eE+-]+/)) {
            n += substr(s, RSTART + 1, RLENGTH - 1)
            s = substr(s, RSTART + RLENGTH)
        }
        print n
    }' "${GENE_TREES}" | sort -g \
| awk '{v[NR] = $1; s += $1}
    END {printf "Tree length: min %.3f  median %.3f  mean %.3f  max %.3f\n",
         v[1], v[int((NR + 1) / 2)], s / NR, v[NR]}'

### rate model chosen per gene (substitution matrix fixed to GTR)
grep -h 'Best-fit model' "${TREE_DIR}"/*.iqtree 2>/dev/null \
    | sed 's/.*: //; s/ chosen.*//' | sort | uniq -c | sort -rn \
    | awk 'BEGIN {print "model\tgenes"} {print $2 "\t" $1}' > "${MODELS}"

echo "Model chosen per gene:"
awk -F'\t' 'NR > 1 {t += $2; printf "  %-14s %d\n", $1, $2}
    END {printf "  %-14s %d\n", "total", t}' "${MODELS}"


### =====================================================================
### step 3: reference trees
### =====================================================================

echo ""
echo "=== Step 3: Reference trees ==="

### reference trees as they are (the 052 tree carries the MLG representatives only)
for i in "${!REF_NAMES[@]}"; do
    echo "${REF_NAMES[$i]}: ${REF_PATHS[$i]} ($(awk '{print gsub(/,/, ",") + 1}' "${REF_PATHS[$i]}") tips)"
done

### gene trees pruned to the tips of the concatenation tree and unrooted (--gcf needs every gene tree tip in the reference)
conda deactivate
conda activate cg_bioinfo_v1_R

Rscript -e '
  a    <- commandArgs(trailingOnly = TRUE)
  keep <- ape::read.tree(a[1])$tip.label
  gt   <- ape::read.tree(a[2])
  gt   <- lapply(gt, function(x) ape::unroot(ape::drop.tip(x, setdiff(x$tip.label, keep))))
  class(gt) <- "multiPhylo"
  ape::write.tree(gt, a[3])
  cat("gene trees pruned to the concatenation tips:", length(gt), "\n")' \
  "${CONCAT_TREE}" "${GENE_TREES}" "${GENE_TREES_DECLONED}"

conda deactivate
conda activate cg_bioinfo_v1_iqtree

## check if every gene tree was pruned
if [[ "$(wc -l < "${GENE_TREES_DECLONED}")" -ne "${N_TREES}" ]]; then
    echo "ERROR: $(wc -l < "${GENE_TREES_DECLONED}") pruned gene trees, ${N_TREES} expected"
    exit 1
fi


### =====================================================================
### step 4: gene concordance factors
### =====================================================================

### function: report_cf (gCF summary per reference, the same for each)
report_cf () {

    local NAME="$1"
    local STAT="$2"

    echo ""
    echo "--- ${NAME}: distribution of gCF ---"

    awk -F'\t' '/^#/ {next}
        h == 0 {for (i = 1; i <= NF; i++) c[$i] = i; h = 1; next}
        $(c["gCF"]) != "NA" {print $(c["gCF"])}' "${STAT}" \
    | sort -g \
    | awk '{v[NR] = $1; s += $1; if ($1 >= 50) hi++; if ($1 < 10) lo++}
        END {printf "branches   %d\nmin        %.1f\nQ1         %.1f\nmedian     %.1f\nmean       %.1f\nQ3         %.1f\nmax        %.1f\nabove 50   %d\nbelow 10   %d\n",
             NR, v[1], v[int(NR * 0.25)], v[int((NR + 1) / 2)], s / NR,
             v[int(NR * 0.75)], v[NR], hi + 0, lo + 0}'

    echo ""
    echo "--- ${NAME}: decisive genes per branch ---"

    ### decisive genes per branch (low gCF over few genes: no information; over many: conflict)
    awk -F'\t' '/^#/ {next}
        h == 0 {for (i = 1; i <= NF; i++) c[$i] = i; h = 1; next}
        {print $(c["gN"])}' "${STAT}" \
    | sort -n \
    | awk '{v[NR] = $1} END {printf "min %d  Q1 %d  median %d  Q3 %d  max %d\n",
         v[1], v[int(NR * 0.25)], v[int((NR + 1) / 2)], v[int(NR * 0.75)], v[NR]}'

    echo ""
    echo "--- ${NAME}: support against concordance ---"

    ### branches with high support and low gene concordance
    awk -F'\t' '/^#/ {next}
        h == 0 {for (i = 1; i <= NF; i++) c[$i] = i; h = 1; next}
        c["Label"] && $(c["Label"]) != "" && $(c["gCF"]) != "NA" {
            n++
            b = $(c["Label"]) + 0
            g = $(c["gCF"]) + 0
            if (b >= 95 && g >= 50) a++
            if (b >= 95 && g <  50 && g >= 20) m++
            if (b >= 95 && g <  20) x++
            if (b <  95) w++
        }
        END {
            if (n == 0) {print "no branch labels on the tree, skipped"; exit}
            printf "branches with a label      %d\n", n
            printf "support >= 95, gCF >= 50    %d\n", a + 0
            printf "support >= 95, gCF 20 to 50 %d\n", m + 0
            printf "support >= 95, gCF < 20     %d\n", x + 0
            printf "support < 95                %d\n", w + 0
        }' "${STAT}"
}

for i in "${!REF_NAMES[@]}"; do

    NAME="${REF_NAMES[$i]}"
    REF="${REF_PATHS[$i]}"
    PREFIX="${RESULTS}/gene_concordance_${NAME}"

    ### the pruned gene trees for the concatenation tree
    GT="${GENE_TREES}"
    [[ "${NAME}" == "concat" ]] && GT="${GENE_TREES_DECLONED}"

    echo ""
    echo "=== Step 4: Gene concordance against the ${NAME} tree ==="

    ### gCF (share of gene trees recovering each branch of the reference tree)
    iqtree2 -t "${REF}" \
        --gcf "${GT}" \
        -T "${THREADS}" \
        --prefix "${PREFIX}" \
        -redo

    ## check if IQ-TREE wrote the concordance statistics
    if [[ ! -s "${PREFIX}.cf.stat" ]]; then
        echo "ERROR: no concordance statistics written for ${NAME}"
        continue
    fi

    report_cf "${NAME}" "${PREFIX}.cf.stat"

done

### output files (copy gene_concordance_<tree>.cf.stat, gene_concordance_<tree>.cf.branch and gene_models.tsv to local/02_structure_phylogeny/data/090_data_gene_concordance for local/02_structure_phylogeny/090_gene_concordance.R and 100_cluster_monophylly.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
for NAME in "${REF_NAMES[@]}"; do
    echo "  gene_concordance_${NAME}.cf.tree   tree with gCF on the branches"
    echo "  gene_concordance_${NAME}.cf.stat   gCF and decisive gene counts per branch"
done
echo "  gene_trees.treefile                the ${N_TREES} gene trees, branches < 1e-5 collapsed"
echo "  gene_trees_decloned.treefile       the same pruned to the tips of the concatenation tree"
echo "  gene_models.tsv                    rate model chosen per gene, counted"
echo "  missing_gene_trees.txt             orthogroups without a tree"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Gene concordance completed in ${RUNTIME} seconds ($((RUNTIME / 3600))h $(((RUNTIME % 3600) / 60))m)."
