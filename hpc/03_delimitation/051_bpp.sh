#!/bin/bash
#SBATCH --job-name=bpp        ### job name
#SBATCH --array=1-24          ### array (one task per row of runs.tsv: 18 A00, 6 A11)
#SBATCH --cpus-per-task=20    ### cores
#SBATCH --mem-per-cpu=1G      ### memory per core
#SBATCH --time=16-00:00:00    ### walltime
#SBATCH --output /dev/null    ### stdout (logged by tee)
#SBATCH --error /dev/null     ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/03_delimitation/results/051_results_bpp"
LOGS="${ROOT}/03_delimitation/logs/051_logs_bpp"
### from hpc/03_delimitation/050_bpp_input.sh (run table and control files)
INPUT="${ROOT}/03_delimitation/results/050_results_bpp_input"

RUNS="${INPUT}/runs.tsv"

### parameters
### array task (first argument outside SLURM)
TASK_ID="${SLURM_ARRAY_TASK_ID:-$1}"
### oldest BPP version with every option used (int, checkpoint, --extend)
BPP_MIN=4.8.0
### samples added to a finished run from its last checkpoint (0 = none; set at submission: sbatch --array=<rows> --export=ALL,EXTEND=<samples> 051_bpp.sh)
EXTEND=${EXTEND:-0}
### C locale (fixed sort order)
export LC_ALL=C

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_ARRAY_JOB_ID:-local}_${TASK_ID}_bpp.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_bpp

START_TIME=$(date +%s)

### run info
echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"

### BPP version
BPP_VERSION=$(bpp --version 2>&1 | grep -o 'v[0-9]\+\.[0-9]\+\.[0-9]\+' | head -1 | tr -d v)
echo "BPP version: ${BPP_VERSION}"

## check if BPP is recent enough
if [[ "$(printf '%s\n' "${BPP_MIN}" "${BPP_VERSION}" | sort -V | head -1)" != "${BPP_MIN}" ]]; then
    echo "ERROR: BPP ${BPP_VERSION} is older than ${BPP_MIN}"
    exit 1
fi

## check if the run table is there
if [[ ! -s "${RUNS}" ]]; then
    echo "ERROR: input not found: ${RUNS}"
    exit 1
fi


### =====================================================================
### step 1: the run of this task
### =====================================================================

### row of runs.tsv (header first): analysis, k, set, chain, run, seed
read -r ANALYSIS K SET CHAIN RUN SEED <<< "$(sed -n "$((TASK_ID + 1))p" "${RUNS}")"

## check if the row holds a run
if [[ -z "${RUN}" ]]; then
    echo "ERROR: no run in row ${TASK_ID} of ${RUNS}"
    exit 1
fi

DIR="${RESULTS}/${RUN}"
CTL="${DIR}/${RUN}.ctl"

mkdir -p "${DIR}"
cd "${DIR}" || exit 1

echo ""
echo "=== Step 1: ${RUN} (${ANALYSIS}, ${K}, locus set ${SET}, chain ${CHAIN}, seed ${SEED}) ==="

### a finished run is skipped, or extended once by EXTEND samples
EXTEND_NOW=0
if [[ -f "${DIR}/done" ]]; then
    if [[ "${EXTEND}" -gt 0 && ! -f "${DIR}/extended_${EXTEND}" ]]; then
        EXTEND_NOW=${EXTEND}
        rm -f "${DIR}/done"
        echo "Run finished before, extended by ${EXTEND_NOW} samples."
    else
        echo "Run finished before, nothing to do."
        exit 0
    fi
fi

### control file of the run
cp "${INPUT}/runs/${RUN}/${RUN}.ctl" "${CTL}"

## check if the threads match the allocation
THREADS=$(awk -F= '$1 ~ /threads/ {print $2 + 0}' "${CTL}")
[[ "${THREADS}" -eq "${SLURM_CPUS_PER_TASK:-20}" ]] \
  || echo "WARNING: ${THREADS} threads in the control file, ${SLURM_CPUS_PER_TASK:-20} cores allocated"


### =====================================================================
### step 2: MCMC (resumed from the last checkpoint if one exists)
### =====================================================================

echo ""
echo "=== Step 2: MCMC ==="

### last checkpoint, if any
LAST_CHK=$(ls -v "${RUN}".*.chk 2>/dev/null | tail -1)

### threads not pinned to fixed cores (other jobs may share the node); an interrupted extension resumes without --extend (the checkpoint holds the extended length)
if [[ "${EXTEND_NOW}" -gt 0 ]]; then
    ## check if the final checkpoint is there
    if [[ -z "${LAST_CHK}" ]]; then
        echo "ERROR: no checkpoint to extend from"
        touch "${DIR}/done"
        exit 1
    fi
    echo "Extending from ${LAST_CHK}"
    bpp --resume "${LAST_CHK}" --extend "${EXTEND_NOW}" --no-pin
elif [[ -n "${LAST_CHK}" ]]; then
    echo "Resuming from ${LAST_CHK}"
    bpp --resume "${LAST_CHK}" --no-pin
else
    bpp --cfile "${CTL}" --no-pin
fi

## check if BPP stopped without error
STATUS=$?
if [[ "${STATUS}" -ne 0 ]]; then
    echo "ERROR: BPP stopped with status ${STATUS}"
    exit 1
fi


### =====================================================================
### step 3: summary and checks
### =====================================================================

echo ""
echo "=== Step 3: Summary ==="

### posterior summary of the full sample (A00: parameters; A11: delimitations and trees)
bpp --summary "${CTL}" --no-pin > "${RUN}.summary.txt" 2>&1

## check if the sample is complete (A00 has a header line)
NSAMPLE=$(awk -F= '$1 ~ /nsample/ {print $2 + 0}' "${CTL}")
N_LINES=$(wc -l < "${RUN}.mcmc.txt")
[[ "${ANALYSIS}" == "A00" ]] && N_LINES=$((N_LINES - 1))
echo "Samples: ${N_LINES} (nsample in the control file: ${NSAMPLE})"
if [[ "${N_LINES}" -lt "${NSAMPLE}" ]]; then
    echo "ERROR: incomplete MCMC sample"
    exit 1
fi

## check if the summary was written
if [[ "${ANALYSIS}" == "A00" ]]; then
    grep -q "^ *param *mean" "${RUN}.summary.txt" || { echo "ERROR: no parameter summary"; exit 1; }
else
    grep -q "Posterior probability for # of species" "${RUN}.summary.txt" || { echo "ERROR: no delimitation summary"; exit 1; }
fi

### the last checkpoint is kept (to extend the run with --resume and --extend)
ls -v "${RUN}".*.chk 2>/dev/null | head -n -1 | xargs -r rm -f

### markers of a finished (and extended) run
touch "${DIR}/done"
[[ "${EXTEND}" -gt 0 ]] && touch "${DIR}/extended_${EXTEND}"

### output files (copy every run folder (<run>.ctl, <run>.mcmc.txt, <run>.summary.txt, extended_* markers) to local/03_delimitation/data/051_data_bpp for local/03_delimitation/051_bpp.R)
echo ""
echo "=== Output files in ${DIR} ==="
echo "  ${RUN}.ctl           control file (copied from 050)"
echo "  ${RUN}.mcmc.txt      MCMC sample"
echo "  ${RUN}.summary.txt   posterior summary"
echo "  ${RUN}.txt           BPP output file"
echo "  ${RUN}.N.chk         last checkpoint"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "BPP run completed in ${RUNTIME} seconds ($((RUNTIME / 3600))h $(((RUNTIME % 3600) / 60))m)."
