#!/bin/bash
#SBATCH --job-name=cds_extraction       ### job name
#SBATCH --cpus-per-task=1               ### cores
#SBATCH --mem-per-cpu=8G                ### memory per core
#SBATCH --time=02:00:00                 ### walltime
#SBATCH --output /dev/null              ### stdout (logged by tee)
#SBATCH --error /dev/null               ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/050_results_sco_consensus"
LOGS="${ROOT}/01_variant_calling/logs/050_logs_sco_consensus"
### from hpc/01_variant_calling/040_orthofinder.sh (replace Results_Date with the run-date folder OrthoFinder wrote, e.g. Results_Feb18)
ORTHOFINDER="${ROOT}/01_variant_calling/results/040_results_orthofinder/Results_Date"
### external: coding sequences of the reference genomes, one FASTA per genome (copy manually into ${ROOT}/01_variant_calling/data/references/coordinates)
CDS_DIR="${ROOT}/01_variant_calling/data/references/coordinates"

FASTA_DIR="${RESULTS}/per_og_fasta"
SCO_TABLE="${RESULTS}/sco_orthogroups.tsv"
OG_LIST="${RESULTS}/og_list.txt"
TMP="${RESULTS}/tmp_maps"

### create results folders
mkdir -p "${FASTA_DIR}" "${TMP}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_cds_extraction.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_samtools

START_TIME=$(date +%s)

echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"


### =====================================================================
### step 1: input line endings
### =====================================================================

### Windows line endings stripped (they break the column parsing)
echo ""
echo "=== Sanitising input files ==="

for FILE in "${CDS_DIR}"/*_cds_*.fasta \
            "${ORTHOFINDER}/Orthogroups/Orthogroups.tsv" \
            "${ORTHOFINDER}/Orthogroups/Orthogroups_SingleCopyOrthologues.txt"; do
    if file "${FILE}" | grep -q CRLF; then
        echo "  converting $(basename "${FILE}")"
        sed -i 's/\r$//' "${FILE}"
    fi
done


### =====================================================================
### step 2: single-copy orthogroups
### =====================================================================

### single-copy rows of the orthogroup table (matched on the first column, so a blank line cannot match everything)
echo ""
echo "=== Selecting single-copy orthogroups ==="

awk -F'\t' '
FNR == NR { if ($1 != "") sco[$1]; next }
FNR == 1 || ($1 in sco)
' "${ORTHOFINDER}/Orthogroups/Orthogroups_SingleCopyOrthologues.txt" \
  "${ORTHOFINDER}/Orthogroups/Orthogroups.tsv" > "${SCO_TABLE}"

N_SCO=$(tail -n +2 "${SCO_TABLE}" | wc -l)
echo "  single-copy orthogroups: ${N_SCO}"

### genome names (column headers of the orthogroup table)
mapfile -t GENOMES < <(head -1 "${SCO_TABLE}" | cut -f2- | tr '\t' '\n')
N_GENOMES=${#GENOMES[@]}

echo "  genomes (including outgroup): ${N_GENOMES}"


### =====================================================================
### step 3: CDS file per genome
### =====================================================================

### species and version prefix of the genome name (e.g. Cenge1058_1_FilteredModels1_... -> Cenge1058_1)
echo ""
echo "=== Locating CDS files ==="

declare -A CDS_FILE PREFIX

for G in "${GENOMES[@]}"; do
    P=$(echo "${G}" | sed -E 's/^([A-Za-z]+[0-9]+_[0-9]+).*/\1/')
    F=$(ls "${CDS_DIR}/${P}_"*_cds_*.fasta 2>/dev/null | head -1)

    ## check if the genome has a CDS file
    if [[ -z "${F}" ]]; then
        echo "ERROR: no CDS file for ${G} (prefix ${P})"
        exit 1
    fi

    PREFIX["${G}"]="${P}"
    CDS_FILE["${G}"]="${F}"
    echo "  ${P} -> $(basename "${F}")"
done


### =====================================================================
### step 4: CDS per orthogroup
### =====================================================================

### one pass per genome (each sequence appended to the FASTA of its orthogroup)
echo ""
echo "=== Extracting CDS per orthogroup ==="

rm -f "${FASTA_DIR}"/*.fa

for I in "${!GENOMES[@]}"; do
    G="${GENOMES[$I]}"
    COL=$((I + 2))   ### column 1 holds the orthogroup name

    ### protein identifier -> orthogroup, this genome only
    awk -F'\t' -v col="${COL}" 'NR > 1 && $col != "" {print $col "\t" $1}' \
        "${SCO_TABLE}" > "${TMP}/map.tsv"

    awk -v label="${PREFIX[$G]}" -v dir="${FASTA_DIR}" '
    FNR == NR { og[$1] = $2; next }
    /^>/ {
        if (out) close(out)
        id  = substr($1, 2)
        out = (id in og) ? dir "/" og[id] ".fa" : ""
        if (out) print ">" label >> out
        next
    }
    out { print >> out }
    END { if (out) close(out) }
    ' "${TMP}/map.tsv" "${CDS_FILE[$G]}"

    echo "  ${PREFIX[$G]} done ($((I + 1)) / ${N_GENOMES})"
done

rm -rf "${TMP}"


### =====================================================================
### step 5: complete orthogroups
### =====================================================================

### orthogroups with a sequence from every genome
echo ""
echo "=== Checking completeness ==="

N_DROPPED=0

for FA in "${FASTA_DIR}"/*.fa; do
    if [[ $(grep -c '^>' "${FA}") -ne ${N_GENOMES} ]]; then
        echo "  dropping $(basename "${FA}" .fa): incomplete"
        rm -f "${FA}"
        N_DROPPED=$((N_DROPPED + 1))
    fi
done

ls "${FASTA_DIR}"/*.fa | xargs -n1 basename | sed 's/\.fa$//' | sort > "${OG_LIST}"
N_OK=$(wc -l < "${OG_LIST}")

echo "  complete:  ${N_OK}"
echo "  dropped:   ${N_DROPPED}"
echo "  missing:   $((N_SCO - N_OK - N_DROPPED))"

### output files
echo ""
echo "  orthogroup FASTAs: ${FASTA_DIR}/"
echo "  orthogroup list:   ${OG_LIST}"

### runtime
END_TIME=$(date +%s)
echo ""
echo "CDS extraction completed in $((END_TIME - START_TIME)) seconds."
