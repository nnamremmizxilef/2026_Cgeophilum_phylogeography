#!/bin/bash
#SBATCH --job-name=consensus_build      ### job name
#SBATCH --cpus-per-task=4               ### cores
#SBATCH --mem-per-cpu=8G                ### memory per core
#SBATCH --time=04:00:00                 ### walltime
#SBATCH --output /dev/null              ### stdout (logged by tee)
#SBATCH --error /dev/null               ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/054_results_sco_consensus_core"
LOGS="${ROOT}/01_variant_calling/logs/054_logs_sco_consensus_core"
### from hpc/01_variant_calling/055_consensus_align_trim_core.sh
ALN_DIR="${RESULTS}/per_og_aln"
### from hpc/01_variant_calling/055_consensus_align_trim_core.sh
TRIM_DIR="${RESULTS}/per_og_trimmed"
### from hpc/01_variant_calling/054_cds_extraction_core.sh
OG_LIST="${RESULTS}/og_list.txt"
CONS_DIR="${RESULTS}/per_og_consensus"
STATS="${RESULTS}/stats"
REFERENCE="${RESULTS}/consensus_sco_reference_core.fasta"

### create results folders
mkdir -p "${CONS_DIR}" "${STATS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_consensus_build.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda
eval "$(conda shell.bash hook)"
conda deactivate
conda activate HTT_emboss

START_TIME=$(date +%s)

echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"


### =====================================================================
### step 1: trimmed alignments
### =====================================================================

echo ""
echo "=== Checking trimmed alignments ==="

MISSING=0
while read -r OG; do
    [[ -s "${TRIM_DIR}/${OG}.trimmed.fa" ]] || { echo "  missing: ${OG}"; MISSING=$((MISSING + 1)); }
done < "${OG_LIST}"

## check if every orthogroup has a trimmed alignment
if [[ ${MISSING} -gt 0 ]]; then
    echo "ERROR: ${MISSING} trimmed alignments missing, resubmit the failed 051 tasks."
    exit 1
fi

TOTAL_OGS=$(wc -l < "${OG_LIST}")

### consensus threshold (strict majority of the genomes per alignment)
FIRST_OG=$(head -1 "${OG_LIST}")
N_GENOMES=$(grep -c '^>' "${TRIM_DIR}/${FIRST_OG}.trimmed.fa")
THRESHOLD=$(( N_GENOMES / 2 + 1 ))

echo "  genomes per alignment: ${N_GENOMES}"
echo "  consensus threshold:   ${THRESHOLD} of ${N_GENOMES} (plurality, strict majority)"


### =====================================================================
### step 2: consensus sequences
### =====================================================================

### consensus per orthogroup (plurality alone decides, N below the threshold, lower case at exactly the threshold)
echo ""
echo "=== Building consensus sequences ==="

FAILED=0
while read -r OG; do
    cons -sequence "${TRIM_DIR}/${OG}.trimmed.fa" \
         -outseq "${CONS_DIR}/${OG}_consensus.fasta" \
         -plurality "${THRESHOLD}" \
         -identity 0 \
         -setcase "${THRESHOLD}" \
         -name "${OG}_consensus" \
         -auto 2>/dev/null

    [[ -s "${CONS_DIR}/${OG}_consensus.fasta" ]] || { echo "  failed: ${OG}"; FAILED=$((FAILED + 1)); }
done < "${OG_LIST}"

echo "  built:  $((TOTAL_OGS - FAILED))"
echo "  failed: ${FAILED}"


### =====================================================================
### step 3: mapping reference
### =====================================================================

### consensus sequences concatenated into the mapping reference (case kept for the base composition)
echo ""
echo "=== Building the mapping reference ==="

: > "${REFERENCE}"
while read -r OG; do
    CONS="${CONS_DIR}/${OG}_consensus.fasta"
    [[ -s "${CONS}" ]] && cat "${CONS}" >> "${REFERENCE}"
done < "${OG_LIST}"

N_SEQS=$(grep -c '^>' "${REFERENCE}")
REF_SIZE=$(awk '!/^>/ {n += length($0)} END {print n}' "${REFERENCE}")

echo "  ${N_SEQS} sequences, ${REF_SIZE} bp"


### =====================================================================
### step 4: alignment lengths
### =====================================================================

echo ""
echo "=== Measuring alignments ==="

### function: alignment length per file (whole first record, wrapped lines summed, one gawk call per directory)
measure() {
    local DIR=$1 SUFFIX=$2 OUT=$3
    ls "${DIR}"/*"${SUFFIX}" | xargs gawk -v suffix="${SUFFIX}" '
    FNR == 1 { rec = 0; len = 0 }
    /^>/     { rec++; if (rec == 2) nextfile; next }
    rec == 1 { len += length($0) }
    ENDFILE  { og = FILENAME; sub(/.*\//, "", og); sub(suffix "$", "", og); print og "\t" len }
    ' | sort > "${OUT}"
}

measure "${ALN_DIR}"  ".aln.fa"     "${STATS}/lengths_aln.tsv"
measure "${TRIM_DIR}" ".trimmed.fa" "${STATS}/lengths_trim.tsv"

echo "  written: ${STATS}/lengths_aln.tsv and ${STATS}/lengths_trim.tsv"


### =====================================================================
### step 5: per-orthogroup statistics
### =====================================================================

### agreement between the genomes per orthogroup
echo ""
echo "=== Per-orthogroup statistics ==="

cat > "${RESULTS}/compute_og_stats.py" <<'PYEOF'
#!/usr/bin/env python3
"""How well do the genomes agree, orthogroup by orthogroup?
 
Sites are classified the same way EMBOSS cons decides a base, by an absolute
count out of all genomes: unanimous if every genome carries the same base,
majority if the commonest base reaches the strict-majority threshold, below
threshold otherwise. Sites below the threshold are exactly the ones cons
writes as N, and the summary checks that the two counts agree.
 
Occupancy and agreement are reported separately so that gaps and genuine
disagreement are not confused: occupancy is the fraction of genomes carrying
any base at a site, agreement the fraction of those that carry the commonest.
 
Usage: compute_og_stats.py trim_dir cons_dir og_list threshold output.tsv
"""
import os
import sys
from collections import Counter
 
HEADER = ["orthogroup", "n_seqs", "trimmed_length", "consensus_length",
          "n_gap_only_cols", "n_unanimous", "n_majority", "n_below_threshold",
          "mean_agreement", "mean_occupancy", "cons_resolved", "cons_n"]
 
 
def read_fasta(path):
    name, seq = None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.upper())
    if name is not None:
        yield name, "".join(seq)
 
 
def stats(trim_file, cons_file, og, threshold):
    aln = [s for _, s in read_fasta(trim_file)]
    cons = next((s for _, s in read_fasta(cons_file)), "")
    n_genomes = len(aln)
 
    unanimous = majority = below = gap_only = 0
    agreement, occupancy = [], []
 
    for i in range(len(aln[0]) if aln else 0):
        counts = Counter(s[i] for s in aln if s[i] in "ACGT")
        if not counts:
            gap_only += 1
            continue
 
        top = counts.most_common(1)[0][1]
        present = sum(counts.values())
 
        agreement.append(top / present)
        occupancy.append(present / n_genomes)
 
        ## the same rule cons applies: an absolute count, gaps do not help
        if top == n_genomes:
            unanimous += 1
        elif top >= threshold:
            majority += 1
        else:
            below += 1
 
    return {
        "orthogroup": og,
        "n_seqs": n_genomes,
        "trimmed_length": len(aln[0]) if aln else 0,
        "consensus_length": len(cons),
        "n_gap_only_cols": gap_only,
        "n_unanimous": unanimous,
        "n_majority": majority,
        "n_below_threshold": below,
        "mean_agreement": sum(agreement) / len(agreement) if agreement else 0,
        "mean_occupancy": sum(occupancy) / len(occupancy) if occupancy else 0,
        ## count both cases, EMBOSS lowercases bases carried by exactly the
        ## threshold number of genomes
        "cons_resolved": sum(1 for b in cons if b in "ACGT"),
        "cons_n": cons.count("N"),
    }
 
 
def main():
    trim_dir, cons_dir, og_list, threshold, out_tsv = sys.argv[1:6]
    threshold = int(threshold)
 
    with open(og_list) as handle:
        ogs = [line.strip() for line in handle if line.strip()]
 
    with open(out_tsv, "w") as out:
        out.write("\t".join(HEADER) + "\n")
 
        for i, og in enumerate(ogs, start=1):
            if i % 500 == 0:
                print(f"  {i} / {len(ogs)}...")
 
            trim_file = os.path.join(trim_dir, f"{og}.trimmed.fa")
            cons_file = os.path.join(cons_dir, f"{og}_consensus.fasta")
 
            if not (os.path.exists(trim_file) and os.path.exists(cons_file)):
                print(f"  skipping {og}, files missing")
                continue
 
            row = stats(trim_file, cons_file, og, threshold)
            out.write("\t".join(
                f"{row[k]:.4f}" if isinstance(row[k], float) else str(row[k])
                for k in HEADER) + "\n")
 
    print(f"  statistics for {len(ogs)} orthogroups written.")
 
 
if __name__ == "__main__":
    main()
PYEOF

PER_OG="${STATS}/per_orthogroup_stats.tsv"
python3 "${RESULTS}/compute_og_stats.py" "${TRIM_DIR}" "${CONS_DIR}" "${OG_LIST}" "${THRESHOLD}" "${PER_OG}"


### =====================================================================
### step 6: base composition
### =====================================================================

### base composition of the reference (upper case above the threshold, lower case at exactly the threshold, both resolved)
echo ""
echo "=== Base composition ==="

### header outside the pipe (sort would move it)
BASES="${STATS}/base_composition.tsv"
{
    printf 'base\tcount\n'
    awk '!/^>/ { for (i = 1; i <= length($0); i++) c[substr($0, i, 1)]++ }
    END { for (b in c) print b "\t" c[b] }' "${REFERENCE}" | sort
} > "${BASES}"

read -r N_UPPER N_LOWER N_AMBIG N_OTHER < <(awk -F'\t' '
NR > 1 {
    if      ($1 ~ /^[ACGT]$/) up += $2
    else if ($1 ~ /^[acgt]$/) lo += $2
    else if ($1 ~ /^[Nn]$/)   am += $2
    else                      ot += $2
}
END { print up + 0, lo + 0, am + 0, ot + 0 }' "${BASES}")

N_RESOLVED=$((N_UPPER + N_LOWER))

echo "  resolved, more than ${THRESHOLD} genomes:  ${N_UPPER}"
echo "  resolved, exactly ${THRESHOLD} genomes:    ${N_LOWER}"
echo "  ambiguous (N):                     ${N_AMBIG}"
echo "  unexpected characters:             ${N_OTHER}"
echo "  resolved in total:                 ${N_RESOLVED}"


### =====================================================================
### step 7: reference in upper case
### =====================================================================

### upper case only (lower case reads as soft masking in some tools)
echo ""
echo "=== Normalising the reference ==="

awk '/^>/ {print; next} {print toupper($0)}' "${REFERENCE}" > "${REFERENCE}.tmp" \
    && mv "${REFERENCE}.tmp" "${REFERENCE}"

echo "  reference written in upper case."


### =====================================================================
### step 8: summary
### =====================================================================

echo ""
echo "=== Summary ==="

MEAN_ALN=$(awk '{n += $2} END {printf "%.0f", n / NR}' "${STATS}/lengths_aln.tsv")
MEAN_TRIM=$(awk '{n += $2} END {printf "%.0f", n / NR}' "${STATS}/lengths_trim.tsv")
PCT_TRIMMED=$(awk 'NR == FNR {a[$1] = $2; next} $1 in a && a[$1] > 0 {
    s += (a[$1] - $2) * 100 / a[$1]; n++} END {printf "%.2f", s / n}' \
    "${STATS}/lengths_aln.tsv" "${STATS}/lengths_trim.tsv")
MEAN_CONS=$(awk -F'\t' 'NR > 1 {n += $4; c++} END {printf "%.0f", n / c}' "${PER_OG}")
MEDIAN_CONS=$(awk -F'\t' 'NR > 1 {print $4}' "${PER_OG}" | sort -n | \
    awk '{v[NR] = $1} END {print (NR % 2) ? v[(NR + 1) / 2] : int((v[NR / 2] + v[NR / 2 + 1]) / 2)}')

awk -F'\t' -v mean_aln="${MEAN_ALN}" -v mean_trim="${MEAN_TRIM}" \
    -v pct="${PCT_TRIMMED}" -v median="${MEDIAN_CONS}" -v thr="${THRESHOLD}" \
    -v ngen="${N_GENOMES}" -v n_amb="${N_AMBIG}" '
NR > 1 { n++; cons += $4; unan += $6; maj += $7; below += $8; agree += $9; occ += $10 }
END {
    sites = unan + maj + below
    printf "Orthogroups:                    %d\n", n
    printf "Mean alignment length:          %s bp\n", mean_aln
    printf "Mean length after trimming:     %s bp\n", mean_trim
    printf "Mean columns removed:           %s%%\n", pct
    printf "Median consensus length:        %s bp\n", median
    printf "Total consensus bp:             %d\n", cons
    printf "\n"
    printf "Mean site occupancy:            %.4f (genomes carrying a base)\n", occ / n
    printf "Mean per-site agreement:        %.4f (among those genomes)\n", agree / n
    printf "\n"
    printf "Across all alignment sites, by number of genomes sharing a base:\n"
    printf "  Unanimous (all %d):            %d (%.2f%%)\n", ngen, unan, unan / sites * 100
    printf "  Majority (%d to %d):           %d (%.2f%%)\n", thr, ngen - 1, maj, maj / sites * 100
    printf "  Below threshold (< %d):        %d (%.2f%%)\n", thr, below, below / sites * 100
    printf "\n"
    printf "Sites below threshold vs N in the reference: %d vs %d%s\n", \
           below, n_amb, (below == n_amb ? "  (consistent)" : "  MISMATCH")
}' "${PER_OG}" | tee "${STATS}/summary_stats.txt"

### distributions for plotting
awk -F'\t' 'NR > 1 {print $9}'  "${PER_OG}" | sort -g > "${STATS}/agreement_distribution.txt"
awk -F'\t' 'NR > 1 {print $10}' "${PER_OG}" | sort -g > "${STATS}/occupancy_distribution.txt"
awk -F'\t' 'NR > 1 {print $4}'  "${PER_OG}" | sort -g > "${STATS}/consensus_length_distribution.txt"


### =====================================================================
### step 9: publication table
### =====================================================================

echo ""
echo "=== Publication table ==="

PUB="${STATS}/publication_table.tsv"
{
    printf "Parameter\tValue\n"

    printf "Input and parameters\t\n"
    printf "  Genomes in the consensus\t%s\n" "${N_GENOMES}"
    printf "  Single-copy orthogroups\t%s\n" "${TOTAL_OGS}"
    printf "  trimAl gap threshold\t0.5\n"
    printf "  EMBOSS cons plurality threshold (strict majority)\t%s of %s\n" \
           "${THRESHOLD}" "${N_GENOMES}"

    printf "Alignments\t\n"
    printf "  Mean alignment length (bp)\t%s\n" "${MEAN_ALN}"
    printf "  Mean length after trimming (bp)\t%s\n" "${MEAN_TRIM}"
    printf "  Mean columns removed (%%)\t%s\n" "${PCT_TRIMMED}"
    printf "  Total columns removed (bp)\t%s\n" \
           "$(awk 'NR == FNR {a[$1] = $2; next} $1 in a {s += a[$1] - $2} END {print s}' \
              "${STATS}/lengths_aln.tsv" "${STATS}/lengths_trim.tsv")"

    printf "Consensus panel\t\n"
    printf "  Consensus sequences\t%s\n" "${N_SEQS}"
    printf "  Total length (bp)\t%s\n" "${REF_SIZE}"
    printf "  Mean consensus length (bp)\t%s\n" "${MEAN_CONS}"
    printf "  Median consensus length (bp)\t%s\n" "${MEDIAN_CONS}"
    printf "  Resolved bases (A/C/G/T)\t%s (%.2f%%)\n" "${N_RESOLVED}" \
           "$(awk -v a="${N_RESOLVED}" -v b="${REF_SIZE}" 'BEGIN {print a / b * 100}')"
    printf "    carried by more than %s genomes\t%s (%.2f%%)\n" "${THRESHOLD}" "${N_UPPER}" \
           "$(awk -v a="${N_UPPER}" -v b="${REF_SIZE}" 'BEGIN {print a / b * 100}')"
    printf "    carried by exactly %s genomes\t%s (%.2f%%)\n" "${THRESHOLD}" "${N_LOWER}" \
           "$(awk -v a="${N_LOWER}" -v b="${REF_SIZE}" 'BEGIN {print a / b * 100}')"
    printf "  Ambiguous bases (N)\t%s (%.2f%%)\n" "${N_AMBIG}" \
           "$(awk -v a="${N_AMBIG}" -v b="${REF_SIZE}" 'BEGIN {print a / b * 100}')"

    printf "Agreement among genomes\t\n"
    awk -F'\t' -v thr="${THRESHOLD}" -v ngen="${N_GENOMES}" '
    NR > 1 { n++; agree += $9; occ += $10; unan += $6; maj += $7; below += $8 }
    END {
        sites = unan + maj + below
        printf "  Mean site occupancy\t%.4f\n", occ / n
        printf "  Mean per-site agreement\t%.4f\n", agree / n
        printf "  Sites shared by all %d genomes (%%)\t%.2f\n", ngen, unan / sites * 100
        printf "  Sites shared by %d to %d genomes (%%)\t%.2f\n", thr, ngen - 1, maj / sites * 100
        printf "  Sites below the threshold (%%)\t%.2f\n", below / sites * 100
    }' "${PER_OG}"
} > "${PUB}"

cat "${PUB}"

### output files
echo ""
echo "  mapping reference: ${REFERENCE}"
echo "  per-orthogroup:    ${PER_OG}"
echo "  summary:           ${STATS}/summary_stats.txt"
echo "  publication table: ${PUB}"

### runtime
END_TIME=$(date +%s)
echo ""
echo "Consensus build completed in $((END_TIME - START_TIME)) seconds."
