#!/bin/bash
#SBATCH --job-name=dsuite     ### job name
#SBATCH --cpus-per-task=6     ### cores
#SBATCH --mem-per-cpu=4G      ### memory per core
#SBATCH --time=24:00:00       ### walltime
#SBATCH --output /dev/null    ### stdout (logged by tee)
#SBATCH --error /dev/null     ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
PHYLO="${ROOT}/02_structure_phylogeny"
RESULTS="${ROOT}/03_delimitation/results/030_results_dsuite"
LOGS="${ROOT}/03_delimitation/logs/030_logs_dsuite"
### from local/03_delimitation/001_lineage_labels.R (copy local/03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv and lineage_labels.tsv.md5 to ${ROOT}/03_delimitation/supplementary)
LABELS_DIR="${ROOT}/03_delimitation/supplementary"
### frozen labels (one representative per MLG: clustered = TRUE)
LABELS="${LABELS_DIR}/lineage_labels.tsv"
### from local/01_variant_calling/103_mating_type.R (copy local/01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv to ${ROOT}/03_delimitation/supplementary)
META="${LABELS_DIR}/sample_metadata_ALL_FIXED.tsv"
### from hpc/02_structure_phylogeny/010_filtering_phylogeny.sh (masked call set, holds the outgroup)
VCF_ALL="${PHYLO}/results/010_results_filtering_phylogeny/consensus_phylogeny.vcf.gz"
### from hpc/02_structure_phylogeny/011_filtering_clustering.sh (decloned SNP pool, its sites)
POOL="${PHYLO}/results/011_results_filtering_clustering/consensus_phylogeny_ingroup_decloned_maf01.vcf.gz"
### from hpc/03_delimitation/010_aster_lineages.sh (unit trees)
TREE_DIR="${ROOT}/03_delimitation/results/010_results_aster_lineages"

KEEP="${RESULTS}/dsuite_samples.txt"
SITES="${RESULTS}/dsuite_sites.tsv"
VCF="${RESULTS}/dsuite_snps.vcf.gz"
SET_DIR="${RESULTS}/sets"
OUT_DIR="${RESULTS}/dtrios"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-6}"
### levels whose units are tested
K_LEVELS="k4 k5 k6"
### outgroup samples (Pseudocenococcum)
OUTGROUP_PREFIX="M4-"
### mating type calls of possible mixed cultures (sample set no_mixed)
MAT_EXCLUDE="MAT12,ambiguous"
### jackknife blocks (Dsuite default)
JK_BLOCKS=20
### C locale (fixed sort order)
export LC_ALL=C

### create results folders
mkdir -p "${RESULTS}" "${SET_DIR}" "${OUT_DIR}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_dsuite.log"
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

## check if the labels are the frozen version
( cd "${LABELS_DIR}" && md5sum -c lineage_labels.tsv.md5 ) \
  || { echo "ERROR: lineage_labels.tsv does not match its md5"; exit 1; }

## check if the inputs are there
for F in "${LABELS}" "${META}" "${VCF_ALL}" "${POOL}"; do
    if [[ ! -s "${F}" ]]; then
        echo "ERROR: input not found: ${F}"
        exit 1
    fi
done

## check if the unit trees are there
for K in ${K_LEVELS}; do
    if [[ ! -s "${TREE_DIR}/cluster_tree_${K}_rooted.tre" ]]; then
        echo "ERROR: unit tree not found: ${TREE_DIR}/cluster_tree_${K}_rooted.tre"
        exit 1
    fi
done


### =====================================================================
### step 1: samples, sets files and trees
### =====================================================================

echo ""
echo "=== Step 1: Samples, sets and trees ==="

bcftools query -l "${VCF_ALL}" > "${RESULTS}/vcf_samples.txt"

conda deactivate
conda activate cg_bioinfo_v1_R

### one representative per MLG plus the outgroup; one sets file per level and sample set (all, no_mixed with mixed cultures as xxx); unit trees with the outgroup named Outgroup
Rscript - "${LABELS}" "${META}" "${RESULTS}/vcf_samples.txt" "${K_LEVELS}" "${OUTGROUP_PREFIX}" \
    "${MAT_EXCLUDE}" "${TREE_DIR}" "${SET_DIR}" "${KEEP}" <<'RSCRIPT'
suppressMessages(library(ape))

a <- as.list(commandArgs(trailingOnly = TRUE))
names(a) <- c("labels", "meta", "vcf_samples", "k_levels", "og_prefix", "mat_exclude", "tree_dir",
              "set_dir", "keep")

k_levels <- strsplit(a$k_levels, " ")[[1]]

lab  <- read.delim(a$labels, colClasses = "character")
meta <- read.delim(a$meta, colClasses = "character")
rep  <- lab[lab$clustered == "TRUE", ]

og    <- grep(paste0("^", a$og_prefix), readLines(a$vcf_samples), value = TRUE)
mixed <- meta$sample[meta$mat_genome %in% strsplit(a$mat_exclude, ",")[[1]]]

writeLines(c(rep$sample, og), a$keep)

for (k in k_levels) {
  for (s in c("all", "no_mixed")) {
    unit <- rep[[k]]
    if (s == "no_mixed") unit[rep$sample %in% mixed] <- "xxx"
    sets <- data.frame(sample = c(rep$sample, og), unit = c(unit, rep("Outgroup", length(og))))
    write.table(sets, file.path(a$set_dir, sprintf("sets_%s_%s.txt", k, s)), sep = "\t",
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }

  tr <- read.tree(file.path(a$tree_dir, sprintf("cluster_tree_%s_rooted.tre", k)))
  tr$tip.label[tr$tip.label == "outgroup"] <- "Outgroup"
  tr$edge.length <- NULL
  write.tree(tr, file.path(a$set_dir, sprintf("tree_%s.nwk", k)))

  ## check if the tree holds exactly the units of the level and the outgroup
  if (!setequal(tr$tip.label, c(unique(rep[[k]]), "Outgroup"))) stop("tree and units differ for ", k)
}

cat("representatives:", nrow(rep), "  outgroup samples:", length(og),
    "  mixed cultures:", sum(rep$sample %in% mixed), "\n")
RSCRIPT

## check if the sets files were written
if [[ $(ls "${SET_DIR}"/sets_*.txt 2>/dev/null | wc -l) -eq 0 ]]; then
    echo "ERROR: no sets files written"
    exit 1
fi

### unit tree per level
for K in ${K_LEVELS}; do
    echo "${K}: $(cat "${SET_DIR}/tree_${K}.nwk")"
done


### =====================================================================
### step 2: SNPs of the pool in the representatives and the outgroup
### =====================================================================

echo ""
echo "=== Step 2: SNPs ==="

conda deactivate
conda activate cg_bioinfo_v1_bcftools

## check if every kept sample is in the 010 call set
N_MISS=$(sort "${KEEP}" | comm -23 - <(sort "${RESULTS}/vcf_samples.txt") | wc -l)
if [[ "${N_MISS}" -gt 0 ]]; then
    echo "ERROR: ${N_MISS} kept samples are not in the 010 call set"
    exit 1
fi

### the sites of the pool, genotyped in the representatives and the outgroup (streamed with -T; -R jumps the index per site)
bcftools query -f '%CHROM\t%POS\n' "${POOL}" > "${SITES}"
bcftools view --threads "${THREADS}" -S "${KEEP}" -T "${SITES}" -m2 -M2 -v snps \
    "${VCF_ALL}" -Oz -o "${VCF}"
bcftools index --threads "${THREADS}" "${VCF}"

N_POOL=$(bcftools index -n "${POOL}")
N_VCF=$(bcftools index -n "${VCF}")
echo "SNPs in the pool: ${N_POOL}   in the Dsuite set: ${N_VCF}"
echo "Samples: $(bcftools query -l "${VCF}" | wc -l)"

## check if the Dsuite set holds exactly the pool sites
[[ "${N_VCF}" -eq "${N_POOL}" ]] || echo "WARNING: the Dsuite set and the pool differ in SNPs"

### SNPs with at least one called outgroup sample
grep "^${OUTGROUP_PREFIX}" "${KEEP}" > "${RESULTS}/outgroup_samples.txt"
N_OG=$(bcftools view -S "${RESULTS}/outgroup_samples.txt" "${VCF}" -Ou \
    | bcftools +fill-tags -Ou -- -t NS \
    | bcftools view -H -i 'INFO/NS > 0' | wc -l)
echo "SNPs with an outgroup call: ${N_OG} of ${N_VCF}"


### =====================================================================
### step 3: D and f4-ratio per trio, f-branch
### =====================================================================

echo ""
echo "=== Step 3: Dtrios and Fbranch ==="

conda deactivate
conda activate cg_bioinfo_v1_dsuite
echo "Running in environment: $CONDA_DEFAULT_ENV"
echo "Dsuite: $(Dsuite 2>&1 | grep -m1 Version)"
echo "Jackknife blocks: ${JK_BLOCKS}"

### one Dtrios run per level and sample set, in parallel (waiting on these runs only, not on the log tee)
PIDS=()
for K in ${K_LEVELS}; do
    for S in all no_mixed; do
        Dsuite Dtrios -k "${JK_BLOCKS}" -c -t "${SET_DIR}/tree_${K}.nwk" -o "${OUT_DIR}/${K}_${S}" \
            "${VCF}" "${SET_DIR}/sets_${K}_${S}.txt" > "${OUT_DIR}/${K}_${S}_dtrios.log" 2>&1 &
        PIDS+=($!)
    done
done
wait "${PIDS[@]}"

### f-branch per level and sample set
for K in ${K_LEVELS}; do
    for S in all no_mixed; do
        Dsuite Fbranch "${SET_DIR}/tree_${K}.nwk" "${OUT_DIR}/${K}_${S}_tree.txt" \
            > "${OUT_DIR}/${K}_${S}_fbranch.txt"
    done
done


### =====================================================================
### checks
### =====================================================================

echo ""
echo "=== Checks ==="

## check if every run wrote all trios of its level (n units choose 3)
for K in ${K_LEVELS}; do
    N_UNITS=$(cut -f2 "${SET_DIR}/sets_${K}_all.txt" | grep -v -x "Outgroup" | sort -u | wc -l)
    N_TRIOS=$(( N_UNITS * (N_UNITS - 1) * (N_UNITS - 2) / 6 ))
    for S in all no_mixed; do
        N_OUT=$(( $(wc -l < "${OUT_DIR}/${K}_${S}_tree.txt") - 1 ))
        echo "${K} ${S}: ${N_OUT} trios (expected ${N_TRIOS})"
        [[ "${N_OUT}" -eq "${N_TRIOS}" ]] || echo "WARNING: ${K} ${S} trios differ from the expected number"
        [[ -s "${OUT_DIR}/${K}_${S}_fbranch.txt" ]] || echo "WARNING: ${K} ${S} wrote no f-branch matrix"
    done
done


### =====================================================================
### step 4: summary
### =====================================================================

echo ""
echo "=== D and f4-ratio per trio, all representatives (ordered by the unit tree) ==="

for K in ${K_LEVELS}; do
    echo ""
    echo "${K}"
    cut -f1-7 "${OUT_DIR}/${K}_all_tree.txt"
done

### output files (copy dtrios/k*_<set>_tree.txt and dtrios/k*_<set>_fbranch.txt to local/03_delimitation/data/030_data_dsuite/dtrios for local/03_delimitation/030_dsuite.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  dsuite_samples.txt           representatives and outgroup samples"
echo "  dsuite_snps.vcf.gz           the pool sites in these samples"
echo "  sets/sets_kX_<set>.txt       sample and unit per level (all; no_mixed with mixed cultures as xxx)"
echo "  sets/tree_kX.nwk             unit tree per level with the outgroup as Outgroup"
echo "  dtrios/kX_<set>_BBAA.txt     D and f4-ratio per trio, P1 and P2 as the pair sharing most derived alleles"
echo "  dtrios/kX_<set>_Dmin.txt     D and f4-ratio per trio at the arrangement with the lowest D"
echo "  dtrios/kX_<set>_tree.txt     D and f4-ratio per trio arranged by the unit tree"
echo "  dtrios/kX_<set>_fbranch.txt  f-branch matrix"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Dsuite completed in ${RUNTIME} seconds ($((RUNTIME / 60)) min)."
