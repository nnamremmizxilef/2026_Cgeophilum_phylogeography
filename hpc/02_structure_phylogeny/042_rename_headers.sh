#!/bin/bash
#SBATCH --job-name=rename_headers    ### job name
#SBATCH --cpus-per-task=16           ### cores
#SBATCH --mem-per-cpu=2G             ### memory per core
#SBATCH --time=02:00:00              ### walltime
#SBATCH --output /dev/null           ### stdout (logged by tee)
#SBATCH --error /dev/null            ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/02_structure_phylogeny/results/042_results_rename_headers"
LOGS="${ROOT}/02_structure_phylogeny/logs/042_logs_rename_headers"
### from hpc/02_structure_phylogeny/040_merge_alignments.sh
INPUT_DIR="${ROOT}/02_structure_phylogeny/results/040_results_merge_alignments/per_og_final"
RENAMED_DIR="${RESULTS}/per_og_renamed"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-16}"

### create results folders
mkdir -p "${RENAMED_DIR}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_rename_headers.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_seqkit

START_TIME=$(date +%s)

### run info
echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"
echo "seqkit version: $(seqkit version)"


### =====================================================================
### step 1: one tip name per taxon
### =====================================================================

### orthogroup prefix stripped (one tip per sample across all genes, reference as "consensus", genomes without prefix; one seqkit call per file, in parallel)
ls -1 "${INPUT_DIR}"/*_final_alignment.fasta \
    | xargs -P "${THREADS}" -I{} bash -c \
      'seqkit replace -p "^OG[0-9]+_" -r "" "$1" -o "'"${RENAMED_DIR}"'/$(basename "$1")"' _ {}

N_IN=$(ls -1 "${INPUT_DIR}"/*_final_alignment.fasta | wc -l)
N_OUT=$(ls -1 "${RENAMED_DIR}"/*_final_alignment.fasta | wc -l)
N_TIPS=$(seqkit seq -n "$(ls "${RENAMED_DIR}"/*_final_alignment.fasta | head -1)" | wc -l)

echo "Files in:  ${N_IN}"
echo "Files out: ${N_OUT}"
echo "Tips in the first file: ${N_TIPS}"

## check if every alignment was renamed
if [[ "${N_OUT}" -ne "${N_IN}" ]]; then
    echo "ERROR: ${N_IN} alignments in, ${N_OUT} out"
    exit 1
fi

## check if a prefix is left (it would make that sample a separate taxon)
if grep -l '^>OG[0-9]*_' "${RENAMED_DIR}"/*_final_alignment.fasta | head -1 | grep -q .; then
    echo "ERROR: orthogroup prefix left in at least one renamed file"
    exit 1
fi

### files per taxon (the same for every taxon when all genes carry the same taxa)
grep -h '^>' "${RENAMED_DIR}"/*_final_alignment.fasta | sed 's/^>//' \
    | sort | uniq -c | awk '{print $1}' | sort -u > "${RESULTS}/tip_counts.txt"

echo "Distinct per-taxon file counts: $(wc -l < "${RESULTS}/tip_counts.txt")"
echo "Files per taxon: $(head -1 "${RESULTS}/tip_counts.txt")"

## check if all genes carry the same taxa
if [[ $(wc -l < "${RESULTS}/tip_counts.txt") -ne 1 ]]; then
    echo "ERROR: the genes do not all carry the same taxa"
    exit 1
fi

### output files
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  per_og_renamed/  ${N_OUT} alignments with one tip name per taxon"
echo "  tip_counts.txt"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Renaming completed in ${RUNTIME} seconds ($((RUNTIME / 60)) min)."
