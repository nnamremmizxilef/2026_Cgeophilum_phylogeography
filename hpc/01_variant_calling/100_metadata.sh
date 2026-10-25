#!/bin/bash
#SBATCH --job-name=metadata         ### job name
#SBATCH --cpus-per-task=1           ### cores
#SBATCH --mem-per-cpu=4G            ### memory per core
#SBATCH --time=00:30:00             ### walltime
#SBATCH --output /dev/null          ### stdout (logged by tee)
#SBATCH --error /dev/null           ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/100_results_metadata"
LOGS="${ROOT}/01_variant_calling/logs/100_logs_metadata"
### from hpc/01_variant_calling/080_filtering.sh
FILTERED_DIR="${ROOT}/01_variant_calling/results/080_results_filtering"
### from hpc/01_variant_calling/060_mapping.sh
MAPPING_DIR="${ROOT}/01_variant_calling/results/060_results_mapping/consensus"
### external: raw paired-end reads (copy manually into ${ROOT}/01_variant_calling/data/sequences/<batch>/raw_reads and ${ROOT}/01_variant_calling/data/outgroups)
SEQ_DIR="${ROOT}/01_variant_calling/data/sequences"
OUTGROUP_DIR="${ROOT}/01_variant_calling/data/outgroups"
SUPP="${ROOT}/01_variant_calling/supplementary/091_supplementary_decloned"

### sample list (608 samples including the outgroup)
SRC_VCF="${FILTERED_DIR}/consensus_sitefiltered_base.vcf.gz"
### ingroup set (samples passing the missingness filter)
INGROUP_VCF="${FILTERED_DIR}/consensus_ingroup_noMAF.vcf.gz"
### orthologues present per sample
SCO_TABLE="${FILTERED_DIR}/sco_presence.tsv"
### per-sample mapping output, one file per sample
FLAGSTAT_DIR="${MAPPING_DIR}/flagstat"
COVERAGE_DIR="${MAPPING_DIR}/coverage"
### made by hand: raw read paths, one read file per line (into ${ROOT}/01_variant_calling/supplementary/020_supplementary_fastqc)
RAW_LIST="${ROOT}/01_variant_calling/supplementary/020_supplementary_fastqc/list1_raw_reads.txt"
### from local/01_variant_calling/090_clone_id.R (copy local/01_variant_calling/results/090_results_clone_id/per_sample_mlg.tsv to ${ROOT}/01_variant_calling/supplementary/091_supplementary_decloned)
MLG_TABLE="${SUPP}/per_sample_mlg.tsv"

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_metadata.log"
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
echo ""
echo "Sample list:      ${SRC_VCF}"
echo "Ingroup set:      ${INGROUP_VCF}"
echo "Orthologues:      ${SCO_TABLE}"
echo "Flagstat:         ${FLAGSTAT_DIR}"
echo "Coverage:         ${COVERAGE_DIR}"
echo "Raw read list:    ${RAW_LIST}"
echo "MLG table:        ${MLG_TABLE}"
echo ""


### =====================================================================
### step 1: lookups
### =====================================================================

echo "=== Step 1: Lookups ==="

TMP="${RESULTS}/tmp"
mkdir -p "${TMP}"

### run (the *_igat_batchN folder in the raw read path; outgroups in data/outgroups)
awk '{
    n = split($0, p, "/"); fn = p[n];
    run = "";
    for (i = 1; i <= n; i++) {
        if (p[i] ~ /_igat_batch[0-9]/) run = p[i]
        if (p[i] == "outgroups")       run = "outgroup"
    }
    s = fn;
    sub(/^[0-9]+_ID[0-9]+_[0-9]+-/, "", s)   ### seqid_IDnnnn_plate-
    sub(/-DNA-.*/, "", s)                    ### cg trailing
    sub(/-P[0-9]+-.*/, "", s)                ### magic and M4 trailing
    sub(/_S[0-9]+_L[0-9]+_R[12].*/, "", s)   ### safety net
    if (s != "" && run != "") print s "\t" run
}' "${RAW_LIST}" | sort -u > "${TMP}/run.tsv"

### reads and primary mapped share from the flagstat files (primary mapped excludes secondary and supplementary; gawk for ENDFILE)
gawk '
function pct(line) {i = index(line, "("); r = substr(line, i + 1); j = index(r, "%"); return substr(r, 1, j - 1)}
FNR==1 {s = FILENAME; sub(/.*\//, "", s); sub(/_flagstat.txt$/, "", s); tot = ""; p = ""; m = ""}
/in total/ {tot = $1}
/primary mapped \(/ {p = pct($0)}
/ mapped \(/ && !/primary/ {m = pct($0)}
ENDFILE {print s "\t" tot "\t" (p != "" ? p : m)}
' "${FLAGSTAT_DIR}"/*_flagstat.txt | sort -u > "${TMP}/flagstat.tsv"

### depth and breadth over the consensus from the coverage files (mean depth weighted by contig length)
gawk -F'\t' '
FNR==1 {s = FILENAME; sub(/.*\//, "", s); sub(/_coverage.txt$/, "", s); wd = 0; cb = 0; ln = 0}
!/^#/ {len = $3 - $2 + 1; ln += len; wd += $7 * len; cb += $5}
ENDFILE {printf "%s\t%.2f\t%.2f\n", s, (ln > 0 ? wd / ln : 0), (ln > 0 ? cb / ln * 100 : 0)}
' "${COVERAGE_DIR}"/*_coverage.txt | sort -u > "${TMP}/coverage.tsv"

### orthologues recovered
awk -F'\t' '
NR==1 {for (i = 1; i <= NF; i++) col[$i] = i; next}
{printf "%s\t%s\t%.1f\n", $(col["sample"]), $(col["n_scos_present"]), $(col["pct_scos_present"])}
' "${SCO_TABLE}" | sort -u > "${TMP}/sco.tsv"

### multilocus genotype, group size and role
awk -F'\t' '
NR==1 {for (i = 1; i <= NF; i++) col[$i] = i; next}
{print $(col["sample"]) "\t" $(col["mlg"]) "\t" $(col["n_members"]) "\t" $(col["status"])}
' "${MLG_TABLE}" | sort -u > "${TMP}/mlg.tsv"

### samples passing the per-sample missingness filter
bcftools query -l "${INGROUP_VCF}" > "${TMP}/ingroup.txt"

echo "  run:          $(wc -l < "${TMP}/run.tsv") samples"
echo "  flagstat:     $(wc -l < "${TMP}/flagstat.tsv") samples"
echo "  coverage:     $(wc -l < "${TMP}/coverage.tsv") samples"
echo "  orthologues:  $(wc -l < "${TMP}/sco.tsv") samples"
echo "  genotype:     $(wc -l < "${TMP}/mlg.tsv") samples"
echo "  ingroup:      $(wc -l < "${TMP}/ingroup.txt") samples"
echo ""

### lookups into memory (one lookup per sample)
declare -A RUN READS MAPPED DEPTH BREADTH SCO_N SCO_PCT MLG MLG_N STATUS INGROUP
while IFS=$'\t' read -r k a;     do RUN["${k}"]="${a}"; done < "${TMP}/run.tsv"
while IFS=$'\t' read -r k a b;   do READS["${k}"]="${a}"; MAPPED["${k}"]="${b}"; done < "${TMP}/flagstat.tsv"
while IFS=$'\t' read -r k a b;   do DEPTH["${k}"]="${a}"; BREADTH["${k}"]="${b}"; done < "${TMP}/coverage.tsv"
while IFS=$'\t' read -r k a b;   do SCO_N["${k}"]="${a}"; SCO_PCT["${k}"]="${b}"; done < "${TMP}/sco.tsv"
while IFS=$'\t' read -r k a b c; do MLG["${k}"]="${a}"; MLG_N["${k}"]="${b}"; STATUS["${k}"]="${c}"; done < "${TMP}/mlg.tsv"
while read -r k;                 do INGROUP["${k}"]=1; done < "${TMP}/ingroup.txt"

### run from the sequence folders (samples missing from the raw read list)
lookup_run_fs () {
    local hit
    hit=$(find "${SEQ_DIR}" "${OUTGROUP_DIR}" -type f -name "*${1}*_R1_*.fastq.gz" 2>/dev/null | head -1)
    if [[ -z "${hit}" ]]; then
        echo ""
    elif [[ "${hit}" == */outgroups/* ]]; then
        echo "outgroup"
    else
        sed -E 's#.*/([a-z]+_igat_batch[0-9]+)/.*#\1#' <<<"${hit}"
    fi
}


### =====================================================================
### step 2: one row per sample
### =====================================================================

echo "=== Step 2: Metadata ==="

### project and sampling design from the name (magic 1-XXX, country CH; outgroup M4-XX; phytoakmeter CC-LLL-DNR: country CC, site CC-LLL, pit D, depth N, replicate R)
OUT="${RESULTS}/sample_metadata.tsv"

{
printf "sample\tproject\tcountry\tsite\tlocation\tdirection\tsampling_depth\treplicate"
printf "\trun_id\tn_reads\tprimary_mapped_pct\tmean_seq_depth\tbreadth_pct"
printf "\tn_scos_present\tpct_scos_present\tin_ingroup"
printf "\tmlg\tclone_group_size\tclone_status\n"

bcftools query -l "${SRC_VCF}" | while read -r s; do

    ### name without the trailing disambiguation underscore
    c="${s%_}"

    project="" ; country="" ; site="" ; location=""
    direction="" ; sdepth="" ; replicate=""

    if [[ "${c}" == 1-* ]]; then
        project="magic" ; country="CH"
    elif [[ "${c}" == M4-* ]]; then
        project="outgroup"
    else
        project="phytoakmeter"
        country=$(cut -d- -f1 <<<"${c}")
        location=$(cut -d- -f2 <<<"${c}")
        site="${country}-${location}"
        field3=$(cut -d- -f3 <<<"${c}")
        direction="${field3:0:1}"
        sdepth=$(tr -cd '0-9' <<<"${field3}")
        replicate="${field3: -1}"
    fi

    run="${RUN[${s}]:-}"
    [[ -z "${run}" ]] && run="${RUN[${c}]:-}"
    [[ -z "${run}" ]] && run=$(lookup_run_fs "${c}")

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${s}" "${project}" "${country}" "${site}" "${location}" \
        "${direction}" "${sdepth}" "${replicate}" "${run}" \
        "${READS[${s}]:-}" "${MAPPED[${s}]:-}" \
        "${DEPTH[${s}]:-}" "${BREADTH[${s}]:-}" \
        "${SCO_N[${s}]:-}" "${SCO_PCT[${s}]:-}" \
        "$([[ -n "${INGROUP[${s}]:-}" ]] && echo TRUE || echo FALSE)" \
        "${MLG[${s}]:-}" "${MLG_N[${s}]:-}" "${STATUS[${s}]:-}"
done
} > "${OUT}"

### clean up
rm -rf "${TMP}"

### report
echo "Rows written: $(( $(wc -l < "${OUT}") - 1 ))"
echo ""
echo "  project:"
awk -F'\t' 'NR>1 {print "    " $2}' "${OUT}" | sort | uniq -c | sort -rn
echo "  country:"
awk -F'\t' 'NR>1 {print "    " ($3 == "" ? "(none)" : $3)}' "${OUT}" | sort | uniq -c | sort -rn
echo "  run_id:"
awk -F'\t' 'NR>1 {print "    " ($9 == "" ? "(none)" : $9)}' "${OUT}" | sort | uniq -c | sort -rn
echo "  in ingroup set:"
awk -F'\t' 'NR>1 {print "    " $16}' "${OUT}" | sort | uniq -c | sort -rn
echo "  clone status:"
awk -F'\t' 'NR>1 {print "    " ($19 == "" ? "(not in clone detection)" : $19)}' "${OUT}" | sort | uniq -c | sort -rn
echo ""

echo "  mapping and coverage, ten lowest depth:"
awk -F'\t' 'NR>1 {printf "    %-14s %10s reads  %6s%% mapped  %7sx depth  %6s%% breadth\n", $1, $10, $11, $12, $13}' \
    "${OUT}" | sort -k6 -g | head -10
echo ""

### empty fields in columns that should be complete
awk -F'\t' 'NR>1 && $9 == "" {print "  no run_id: " $1 " (" $2 ")"}' "${OUT}"
for C in 10:n_reads 11:primary_mapped_pct 12:mean_seq_depth 14:n_scos_present; do
    printf "  empty %-20s %s sample(s)\n" "${C#*:}:" \
        "$(awk -F'\t' -v c="${C%%:*}" 'NR>1 && $c == ""' "${OUT}" | wc -l)"
done

### output files (copy sample_metadata.tsv to local/01_variant_calling/data/103_data_mating_type_collect as sample_metadata_edit.tsv, with the PCR mating type column (mat) added by hand, for local/01_variant_calling/103_mating_type.R)
echo ""
echo "=== Output file in ${RESULTS} ==="
echo "  sample_metadata.tsv"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Completed in ${RUNTIME} seconds."
