#!/bin/bash
#SBATCH --job-name=clustering       ### job name
#SBATCH --cpus-per-task=2           ### cores
#SBATCH --mem-per-cpu=16G           ### memory per core
#SBATCH --time=4-00:00:00           ### walltime
#SBATCH --output /dev/null          ### stdout (logged by tee)
#SBATCH --error /dev/null           ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/02_structure_phylogeny/results/014_results_clustering"
LOGS="${ROOT}/02_structure_phylogeny/logs/014_logs_clustering"

### from hpc/02_structure_phylogeny/013_conversion_clustering.sh (haplotype matrix and sample key)
CONV="${ROOT}/02_structure_phylogeny/results/013_results_conversion_clustering"
INPUT="${CONV}/fineradstructure_input.txt"
KEY="${CONV}/sample_key.tsv"

### external: sampleLD.R of fineRADstructure v0.3.2r109 (bioconda installs the binaries only; copy manually into ${ROOT}/02_structure_phylogeny/supplementary/014_supplementary_clustering)
SAMPLELD="${ROOT}/02_structure_phylogeny/supplementary/014_supplementary_clustering/sampleLD.R"

### output names (RADpainter names its output after its input)
BASE="${RESULTS}/fineradstructure_ldsorted"
SORTED="${BASE}.txt"
CHUNKS="${BASE}_chunks.out"
POPS="${RESULTS}/fineradstructure_populations.tsv"

### parameters
### seed of the LD sampling and the tree
SEED=42
### draws per SCO for the LD estimate (sampleLD.R default)
LD_SAMPLES=500
### MCMC burn-in, sampling and thinning (100 samples per chain)
BURNIN=1000000
SAMPLES=1000000
THIN=10000
### tree: maximisation attempts
TREE_ITER=10000
### two chains (seeds), compared for convergence, the first gives the tree
CHAIN_SEEDS=(1 2)

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_clustering.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

START_TIME=$(date +%s)

### run info
echo "Date: $(date)"

## check if the inputs are there
for F in "${INPUT}" "${KEY}" "${SAMPLELD}"; do
    if [[ ! -s "${F}" ]]; then
        echo "ERROR: missing ${F}"
        exit 1
    fi
done

N_SAMPLES=$(( $(wc -l < "${KEY}") - 1 ))
N_LOCI=$(( $(wc -l < "${INPUT}") - 1 ))

echo ""
echo "Input:   $(basename "${INPUT}")"
echo "Samples: ${N_SAMPLES}"
echo "SCOs:    ${N_LOCI}"


### =====================================================================
### step 1: SCOs ordered by LD
### =====================================================================

echo ""
echo "########## Step 1: LD ordering ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_R
echo "R version: $(R --version | head -1)"

### SCOs in LD next to each other (the RADpainter jackknife blocks see their dependence; the SCO order carries no genome position)
echo "sampleLD.R md5: $(md5sum "${SAMPLELD}" | cut -d' ' -f1)"
Rscript "${SAMPLELD}" -s "${SEED}" -n "${LD_SAMPLES}" "${INPUT}" "${SORTED}" \
  | grep -v "Performing sample"

conda deactivate

## check if the reordering kept every row and the header
if [[ "$(head -1 "${INPUT}")" != "$(head -1 "${SORTED}")" ]]; then
    echo "ERROR: header changed in ${SORTED}"
    exit 1
fi
if ! cmp -s <(tail -n +2 "${INPUT}" | LC_ALL=C sort) <(tail -n +2 "${SORTED}" | LC_ALL=C sort); then
    echo "ERROR: rows of ${SORTED} differ from ${INPUT}"
    exit 1
fi


### =====================================================================
### step 2: coancestry matrix
### =====================================================================

echo ""
echo "########## Step 2: RADpainter ##########"

conda activate cg_bioinfo_v1_fineradstructure
echo "Running in environment: $CONDA_DEFAULT_ENV"
echo "fineradstructure version: $(conda list fineradstructure | gawk '$1 == "fineradstructure" {print $2}')"

## check if both programs are there
for P in RADpainter finestructure; do
    if ! command -v "${P}" > /dev/null; then
        echo "ERROR: ${P} not found in ${CONDA_DEFAULT_ENV}"
        exit 1
    fi
done

### haploid (-p 1)
RADpainter paint -p 1 "${SORTED}" 2>&1 | grep -v "^Processed:"

## check if the coancestry matrix was written
if [[ ! -s "${CHUNKS}" ]]; then
    echo "ERROR: no coancestry matrix written"
    exit 1
fi

## check if the jackknife c is a positive number (finestructure silently uses c = 1 otherwise, e.g. -nan from pairs that never share coancestry)
CFACTOR=$(head -1 "${CHUNKS}" | gawk '$1 == "#Cfactor" {print $2}')
if ! [[ "${CFACTOR}" =~ ^[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$ ]] \
   || ! gawk -v c="${CFACTOR}" 'BEGIN {exit !(c + 0 > 0)}'; then
    echo "ERROR: jackknife c is '${CFACTOR}', not a positive number"
    exit 1
fi
echo "jackknife c:                       ${CFACTOR}"

## check if the matrix holds every sample
N_ROWS=$(gawk 'NR > 2' "${CHUNKS}" | wc -l)
if [[ "${N_ROWS}" -ne "${N_SAMPLES}" ]]; then
    echo "ERROR: ${N_ROWS} rows in the coancestry matrix, ${N_SAMPLES} samples"
    exit 1
fi


### =====================================================================
### step 3: MCMC
### =====================================================================

echo ""
echo "########## Step 3: MCMC ##########"

### chains in parallel (progress bars to their own logs; waited for by pid, a bare wait would also wait for the tee behind the log)
PIDS=()
for S in "${CHAIN_SEEDS[@]}"; do
    finestructure -s "${S}" -x "${BURNIN}" -y "${SAMPLES}" -z "${THIN}" \
        "${CHUNKS}" "${BASE}_chain${S}.mcmc.xml" \
        > "${LOGS}/${SLURM_JOB_ID:-local}_chain${S}.log" 2>&1 &
    PIDS+=($!)
done
wait "${PIDS[@]}"

N_EXP=$(( SAMPLES / THIN ))

for S in "${CHAIN_SEEDS[@]}"; do

    MCMC="${BASE}_chain${S}.mcmc.xml"

    ## check if the chain finished
    N_IT=$(grep -c "<Iteration>" "${MCMC}" 2>/dev/null)
    if [[ "${N_IT:-0}" -lt "${N_EXP}" ]]; then
        echo "ERROR: chain ${S} holds ${N_IT} samples, ${N_EXP} expected"
        exit 1
    fi

    ### number of populations and posterior, first and second half of the samples
    grep -o "<K>[0-9]*</K>\|<Posterior>[-0-9.e]*</Posterior>" "${MCMC}" \
      | sed 's/<[^>]*>//g' | paste - - \
      | gawk -v s="${S}" '{post[NR] = $1; k[NR] = $2}
        END {
            h = int(NR / 2)
            for (i = 1; i <= NR; i++) {
                if (i <= h) {p1 += post[i]; n1++} else {p2 += post[i]; n2++}
                if (min == "" || k[i] < min) min = k[i]
                if (k[i] > max) max = k[i]
                sk += k[i]
            }
            printf "chain %s: %d samples, K %d to %d (mean %.1f), posterior %.1f / %.1f (halves)\n",
                   s, NR, min, max, sk / NR, p1 / n1, p2 / n2
        }'

    ### pairwise coincidence (share of samples in which two individuals cluster)
    finestructure -e meancoincidence "${CHUNKS}" "${MCMC}" \
        "${BASE}_chain${S}_meancoincidence.csv" > /dev/null 2>&1
done

### agreement between the chains (mean and largest difference in pairwise coincidence)
C1="${BASE}_chain${CHAIN_SEEDS[0]}_meancoincidence.csv"
C2="${BASE}_chain${CHAIN_SEEDS[1]}_meancoincidence.csv"
gawk -F',' 'NR == FNR {if (FNR > 1) for (j = 2; j <= NF; j++) a[FNR, j] = $j + 0; next}
    FNR > 1 {for (j = 2; j <= NF; j++) if (FNR != j) {d = $j - a[FNR, j]; if (d < 0) d = -d; s += d; n++; if (d > m) m = d}}
    END {printf "coincidence between chains:        mean |diff| %.4f, max %.4f\n", s / n, m}' "${C1}" "${C2}"


### =====================================================================
### step 4: tree
### =====================================================================

echo ""
echo "########## Step 4: tree ##########"

MCMC1="${BASE}_chain${CHAIN_SEEDS[0]}.mcmc.xml"
TREE="${BASE}_chain${CHAIN_SEEDS[0]}.mcmcTree.xml"

finestructure -s "${SEED}" -m T -x "${TREE_ITER}" "${CHUNKS}" "${MCMC1}" "${TREE}" \
    > "${LOGS}/${SLURM_JOB_ID:-local}_tree.log" 2>&1

conda deactivate

## check if the tree was written
if ! grep -q "<Tree>" "${TREE}" 2>/dev/null; then
    echo "ERROR: no tree in ${TREE}, see ${LOGS}/${SLURM_JOB_ID:-local}_tree.log"
    exit 1
fi

### MAP populations of the tree run, ids back to sample names
grep -m1 -o "<Pop>.*</Pop>" "${TREE}" | sed 's/<\/*Pop>//g; s/)(/\n/g; s/[()]//g' \
  | gawk -v key="${KEY}" 'BEGIN {OFS = "\t"; while ((getline l < key) > 0) {split(l, f, "\t"); name[f[1]] = f[2]}
        print "sample", "id", "population", "population_size"}
    {n = split($0, ids, ","); for (i = 1; i <= n; i++) print name[ids[i]], ids[i], NR, n}' \
  > "${POPS}"

N_POP=$(gawk -F'\t' 'NR > 1 {print $3}' "${POPS}" | sort -u | wc -l)
N_ASSIGNED=$(( $(wc -l < "${POPS}") - 1 ))

## check if every sample is in exactly one population
if [[ "${N_ASSIGNED}" -ne "${N_SAMPLES}" ]] \
   || [[ "$(gawk -F'\t' 'NR > 1 {print $1}' "${POPS}" | sort | uniq -d | wc -l)" -gt 0 ]]; then
    echo "ERROR: ${N_ASSIGNED} assignments for ${N_SAMPLES} samples"
    exit 1
fi

echo "MAP populations:                   ${N_POP}"
echo "population sizes:                  $(gawk -F'\t' 'NR > 1 {print $4 "\t" $3}' "${POPS}" | sort -u | sort -k1,1nr | cut -f1 | paste -sd' ')"


### output files (copy fineradstructure_populations.tsv to local/02_structure_phylogeny/data/014_data_clustering for local/02_structure_phylogeny/014_clustering.R)
echo ""
echo "########## Output ##########"
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  fineradstructure_ldsorted.txt                   input, SCOs ordered by LD"
echo "  fineradstructure_ldsorted_chunks.out            coancestry matrix (with c)"
echo "  fineradstructure_ldsorted_missingness*.out      RADpainter missingness"
echo "  fineradstructure_ldsorted_chain*.mcmc.xml       MCMC samples per chain"
echo "  fineradstructure_ldsorted_chain*_meancoincidence.csv"
echo "  fineradstructure_ldsorted_chain${CHAIN_SEEDS[0]}.mcmcTree.xml   tree and MAP state"
echo "  fineradstructure_populations.tsv                MAP population per sample"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "fineRADstructure completed in ${RUNTIME} seconds ($((RUNTIME / 3600))h $(((RUNTIME % 3600) / 60))m)."
