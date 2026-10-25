#!/bin/bash
#SBATCH --job-name=mat_idiomorphs   ### job name
#SBATCH --cpus-per-task=4           ### cores
#SBATCH --mem-per-cpu=8G            ### memory per core
#SBATCH --time=01:00:00             ### walltime
#SBATCH --output /dev/null          ### stdout (logged by tee)
#SBATCH --error /dev/null           ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/101_results_mat_idiomorphs"
LOGS="${ROOT}/01_variant_calling/logs/101_logs_mat_idiomorphs"

### external: genome assemblies, one mating type each (copy manually into ${ROOT}/01_variant_calling/data/references/assemblies)
ASM_DIR="${ROOT}/01_variant_calling/data/references/assemblies"
ASM_GLOB="Cenge*.fasta"

### external: MAT query proteins and MAT domain HMMs (copy manually into ${ROOT}/01_variant_calling/data/references/mating_type)
MAT_DIR="${ROOT}/01_variant_calling/data/references/mating_type"
MAT_PROTEINS="${MAT_DIR}/mat_query_proteins.fasta"
MAT_HMM="${MAT_DIR}/mat_domains.hmm"

### working folders and output files
WIN_DIR="${RESULTS}/windows"
TMP="${RESULTS}/tmp"

REF="${RESULTS}/mat_idiomorphs.fasta"
REPORT="${RESULTS}/mat_reference_report.txt"
SCAN="${RESULTS}/mat_assembly_scan.tsv"
CHOICE="${RESULTS}/mat_reference_choice.tsv"
ANNOT="${RESULTS}/mat_gene_annotation.tsv"
CORE="${RESULTS}/mat_core.bed"
MASK="${RESULTS}/mat_mask.bed"
CONSERVED="${RESULTS}/mat_conserved.bed"
PRIMERS="${RESULTS}/mat_primer_sites.tsv"
PRIMERS_BED="${RESULTS}/mat_primer_sites.bed"
DIVERGENCE="${RESULTS}/mat_assembly_divergence.tsv"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"

### primers and fragment sizes (sequences withheld, replace XXXXXX with the primer sequences upon publication)
MAT1_F="XXXXXX"
MAT1_R="XXXXXX"
MAT2_F="XXXXXX"
MAT2_R="XXXXXX"

LEN1_EXPECTED=406
LEN2_EXPECTED=353

### mismatches allowed per primer
PRIMER_MM=2

### tolerance on the amplicon length (bp)
LEN_TOL=30

### flank kept either side of the idiomorph (bp)
FLANK=300

### window cut either side of the amplicon (bp)
WINDOW=25000

### shortest blast hit counted as shared sequence (bp)
MIN_HSP=200

### min. identity of a hit in the other mating type (%)
MASK_PIDENT=80

### assemblies the reference is cut from (empty: chosen below)
REF1_ASM=""
REF2_ASM=""

### create results folders
mkdir -p "${RESULTS}" "${WIN_DIR}" "${TMP}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_mat_idiomorphs.log"
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

echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"

: > "${REPORT}"

## check if the query proteins and domain models exist
if [[ ! -s "${MAT_PROTEINS}" || ! -s "${MAT_HMM}" ]]; then
    echo "ERROR: query proteins or domain models missing, see ${MAT_DIR}/README.txt"
    exit 1
fi


### =====================================================================
### step 1: idiomorph and amplicon window per assembly
### =====================================================================

echo ""
echo "=== Step 1: amplicon windows in all assemblies ==="

printf "assembly\tmat1_amplicon_bp\tmat2_amplicon_bp\tidiomorph\n" > "${SCAN}"

for ASM in "${ASM_DIR}"/${ASM_GLOB}; do

    [[ -s "${ASM}" ]] || continue

    BASE=$(basename "${ASM}" .fasta)
    L1=NA
    L2=NA

    for PAIR in "1 ${MAT1_F} ${MAT1_R} ${LEN1_EXPECTED}" \
                "2 ${MAT2_F} ${MAT2_R} ${LEN2_EXPECTED}"; do

        set -- ${PAIR}
        IDX="$1"; FWD="$2"; REV="$3"; EXP="$4"

        BED="${TMP}/${BASE}_amplicon${IDX}.bed"

        ### in silico PCR
        seqkit amplicon -F "${FWD}" -R "${REV}" -m "${PRIMER_MM}" \
            --bed "${ASM}" 2>/dev/null | cut -f1-6 > "${BED}"

        N_HIT=$(wc -l < "${BED}")

        ## check if the primer pair gives a single amplicon
        if [[ "${N_HIT}" -gt 1 ]]; then
            echo "WARNING: MAT1-${IDX} matches ${N_HIT} places in ${BASE}, skipped" \
              | tee -a "${REPORT}"
            continue
        fi

        [[ "${N_HIT}" -eq 1 ]] || continue

        ## check if the amplicon length is within the tolerance
        read -r _ AMP_START AMP_END _ < "${BED}"
        AMP_LEN=$(( AMP_END - AMP_START ))

        if [[ "${AMP_LEN}" -lt $(( EXP - LEN_TOL )) || "${AMP_LEN}" -gt $(( EXP + LEN_TOL )) ]]; then
            echo "WARNING: MAT1-${IDX} amplicon is ${AMP_LEN} bp in ${BASE}," \
                 "expected ${EXP} bp, skipped" | tee -a "${REPORT}"
            continue
        fi

        ### window around the amplicon
        seqkit subseq --bed "${BED}" -u "${WINDOW}" -d "${WINDOW}" \
            "${ASM}" 2>/dev/null \
          | seqkit replace -p '.*' -r "${BASE}_mat${IDX}" \
          > "${WIN_DIR}/${BASE}_mat${IDX}.fasta"

        if [[ "${IDX}" -eq 1 ]]; then L1="${AMP_LEN}"; else L2="${AMP_LEN}"; fi
    done

    ### idiomorph carried
    if   [[ "${L1}" != NA && "${L2}" == NA ]]; then TYPE=MAT1-1
    elif [[ "${L2}" != NA && "${L1}" == NA ]]; then TYPE=MAT1-2
    elif [[ "${L1}" != NA && "${L2}" != NA ]]; then TYPE=both
    else                                           TYPE=none
    fi

    printf "%s\t%s\t%s\t%s\n" "${BASE}" "${L1}" "${L2}" "${TYPE}" >> "${SCAN}"
    echo "  ${BASE}  ${TYPE}"
done

### assemblies per idiomorph
N_ASM=$(tail -n +2 "${SCAN}" | wc -l)
N_M1=$(cut -f4 "${SCAN}" | grep -cx "MAT1-1")
N_M2=$(cut -f4 "${SCAN}" | grep -cx "MAT1-2")
N_BOTH=$(cut -f4 "${SCAN}" | grep -cx "both")
N_NONE=$(cut -f4 "${SCAN}" | grep -cx "none")

{
    echo ""
    echo "assemblies scanned: ${N_ASM}"
    echo "  MAT1-1: ${N_M1}"
    echo "  MAT1-2: ${N_M2}"
    echo "  both:   ${N_BOTH}"
    echo "  none:   ${N_NONE}"
} | tee -a "${REPORT}"

## check if both idiomorphs were found
if [[ "${N_M1}" -lt 1 || "${N_M2}" -lt 1 ]]; then
    echo "ERROR: MAT1-1 windows: ${N_M1}, MAT1-2 windows: ${N_M2}"
    echo "       raise PRIMER_MM, raise LEN_TOL or widen ASM_GLOB"
    exit 1
fi

## check if any assembly carries both primer pairs
if [[ "${N_BOTH}" -gt 0 ]]; then
    echo "WARNING: ${N_BOTH} assemblies carry both primer pairs, see ${SCAN}" \
      | tee -a "${REPORT}"
fi

### reference assembly per mating type
### window length and N count per candidate
printf "assembly\tidiomorph\twindow_bp\texpected_bp\tfull_flank\tn_bases\n" > "${CHOICE}"

while IFS=$'\t' read -r BASE L1 L2 T; do

    case "${T}" in
        MAT1-1) IDX=1; AMPL="${L1}" ;;
        MAT1-2) IDX=2; AMPL="${L2}" ;;
        *)      continue ;;
    esac

    W="${WIN_DIR}/${BASE}_mat${IDX}.fasta"
    [[ -s "${W}" ]] || continue

    read -r WLEN N_BASES < <(seqkit fx2tab -n -l -C N "${W}" | cut -f2,3)

    EXP=$(( 2 * WINDOW + AMPL ))
    FULL=no
    [[ "${WLEN}" -eq "${EXP}" ]] && FULL=yes

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
           "${BASE}" "${T}" "${WLEN}" "${EXP}" "${FULL}" "${N_BASES}" >> "${CHOICE}"

done < <(tail -n +2 "${SCAN}")

### full flank first, then fewest N, then longest window, then name
pick_ref () {
    grep -P "\t$1\t" "${CHOICE}" \
      | sort -t$'\t' -k5,5r -k6,6n -k3,3nr -k1,1 | head -1 | cut -f1
}

[[ -n "${REF1_ASM}" ]] || REF1_ASM=$(pick_ref MAT1-1)
[[ -n "${REF2_ASM}" ]] || REF2_ASM=$(pick_ref MAT1-2)

## check if each mating type has a reference window
if [[ -z "${REF1_ASM}" || -z "${REF2_ASM}" ]]; then
    echo "ERROR: no reference window, MAT1-1: '${REF1_ASM}', MAT1-2: '${REF2_ASM}'"
    exit 1
fi

## check if the reference windows carry the full flank
for R in "${REF1_ASM}" "${REF2_ASM}"; do
    grep -qP "^${R}\t.*\tyes\t" "${CHOICE}" \
      || echo "WARNING: the window of ${R} lacks the full flank" | tee -a "${REPORT}"
done

cp "${WIN_DIR}/${REF1_ASM}_mat1.fasta" "${RESULTS}/mat1_region.fasta"
cp "${WIN_DIR}/${REF2_ASM}_mat2.fasta" "${RESULTS}/mat2_region.fasta"

### stale indices of the region windows
rm -f "${RESULTS}"/mat[12]_region.fasta.fai

{
    echo ""
    echo "reference windows:"
    echo "  MAT1-1  ${REF1_ASM}"
    echo "  MAT1-2  ${REF2_ASM}"
    echo "candidates per mating type and full flank:"
    tail -n +2 "${CHOICE}" | cut -f2,5 | sort | uniq -c
} | tee -a "${REPORT}"


### =====================================================================
### step 2: idiomorphs
### =====================================================================

echo ""
echo "=== Step 2: idiomorphs ==="

### window lengths and forward primer sites
for IDX in 1 2; do
    seqkit fx2tab -n -l "${RESULTS}/mat${IDX}_region.fasta" > "${TMP}/mat${IDX}.genome"
done

seqkit locate -i -P -m "${PRIMER_MM}" -p "${MAT1_F}" --bed \
    "${RESULTS}/mat1_region.fasta" > "${TMP}/mat1_primer.bed"
seqkit locate -i -P -m "${PRIMER_MM}" -p "${MAT2_F}" --bed \
    "${RESULTS}/mat2_region.fasta" > "${TMP}/mat2_primer.bed"

conda deactivate
conda activate cg_bioinfo_v1_blast
echo "Running in environment: $CONDA_DEFAULT_ENV"

### sequence shared by the two windows
for PAIR in "1 2" "2 1"; do

    set -- ${PAIR}
    Q="$1"; S="$2"

    blastn -task blastn -query "${RESULTS}/mat${Q}_region.fasta" \
           -subject "${RESULTS}/mat${S}_region.fasta" \
           -dust no -evalue 1e-10 -outfmt "6 qseqid qstart qend length" 2>/dev/null \
      | awk -v OFS='\t' -v m="${MIN_HSP}" '$4 >= m {print $1, $2 - 1, $3}' \
      > "${RESULTS}/mat${Q}_shared.bed"
done

echo "shared blocks: MAT1-1 $(wc -l < "${RESULTS}/mat1_shared.bed")," \
     "MAT1-2 $(wc -l < "${RESULTS}/mat2_shared.bed")" | tee -a "${REPORT}"

## check if the two windows share sequence
if [[ ! -s "${RESULTS}/mat1_shared.bed" || ! -s "${RESULTS}/mat2_shared.bed" ]]; then
    echo "ERROR: the two windows share no sequence, so no boundary can be drawn"
    echo "       try another pair of assemblies in REF1_ASM and REF2_ASM"
    exit 1
fi

conda deactivate
conda activate cg_bioinfo_v1_bedtools
echo "Running in environment: $CONDA_DEFAULT_ENV"

rm -f "${REF}"

for PAIR in "1 MAT1_1" "2 MAT1_2"; do

    set -- ${PAIR}
    IDX="$1"; NAME="$2"

    G="${TMP}/mat${IDX}.genome"

    ### unshared stretch holding the forward primer
    bedtools sort -i "${RESULTS}/mat${IDX}_shared.bed" \
      | bedtools merge -i - \
      | bedtools complement -i - -g "${G}" \
      | bedtools intersect -a - -b "${TMP}/mat${IDX}_primer.bed" -u \
      > "${TMP}/mat${IDX}_idiomorph.bed"

    ## check if one unshared stretch holds the primer
    if [[ "$(wc -l < "${TMP}/mat${IDX}_idiomorph.bed")" -ne 1 ]]; then
        echo "ERROR: no single unshared block holds the ${NAME} primer"
        exit 1
    fi

    read -r _ B_START B_END < "${TMP}/mat${IDX}_idiomorph.bed"
    read -r _ WLEN < "${G}"

    echo "${NAME} ${B_START}-${B_END} in the window, $(( B_END - B_START )) bp unshared" \
      | tee -a "${REPORT}"

    ## check if the idiomorph reaches the window edge
    if [[ "${B_START}" -eq 0 || "${B_END}" -eq "${WLEN}" ]]; then
        echo "WARNING: ${NAME} reaches the window edge, raise WINDOW" \
          | tee -a "${REPORT}"
    fi

    ### idiomorph with flanks
    bedtools slop -i "${TMP}/mat${IDX}_idiomorph.bed" -g "${G}" -b "${FLANK}" \
      | sed "s/$/\t${NAME}/" \
      | bedtools getfasta -fi "${RESULTS}/mat${IDX}_region.fasta" -bed - -nameOnly \
      > "${TMP}/${NAME}.fasta"

    cat "${TMP}/${NAME}.fasta" >> "${REF}"
done


### =====================================================================
### step 3: index
### =====================================================================

echo ""
echo "=== Step 3: indexing ==="

conda deactivate
conda activate cg_bioinfo_v1_bwa
echo "Running in environment: $CONDA_DEFAULT_ENV"

bwa index "${REF}" 2>/dev/null
samtools faidx "${REF}"

## check if both idiomorphs are in the reference
if [[ "$(grep -c '^>' "${REF}")" -ne 2 ]]; then
    echo "ERROR: the reference holds $(grep -c '^>' "${REF}") sequences, expected 2"
    exit 1
fi

echo "" | tee -a "${REPORT}"
cut -f1,2 "${REF}.fai" | tee -a "${REPORT}"


### =====================================================================
### step 4: gene annotation
### =====================================================================

echo ""
echo "=== Step 4: gene annotation ==="

conda deactivate
conda activate cg_bioinfo_v1_blast
echo "Running in environment: $CONDA_DEFAULT_ENV"

### MAT query proteins against each window
printf "region\tquery\tpident\tlength\tevalue\tbitscore\tsstart\tsend\n" > "${ANNOT}"

for IDX in 1 2; do
    tblastn -query "${MAT_PROTEINS}" -subject "${RESULTS}/mat${IDX}_region.fasta" \
            -outfmt "6 qseqid pident length evalue bitscore sstart send" \
            -evalue 1e-05 -num_threads "${THREADS}" 2>/dev/null \
      | sort -t$'\t' -k5,5gr \
      | sed "s/^/region${IDX}\t/" >> "${ANNOT}"
done

### best hit per window
echo "" | tee -a "${REPORT}"
for IDX in 1 2; do
    grep -m1 -P "^region${IDX}\t" "${ANNOT}" | cut -f1-3,5
done | tee -a "${REPORT}"

conda deactivate
conda activate cg_bioinfo_v1_seqkit
echo "Running in environment: $CONDA_DEFAULT_ENV"

### six frame translation of each window
for IDX in 1 2; do
    seqkit translate -f 6 -F -x "${RESULTS}/mat${IDX}_region.fasta" \
      > "${TMP}/mat${IDX}_frames.faa" 2>/dev/null
done

conda deactivate
conda activate cg_bioinfo_v1_hmmer
echo "Running in environment: $CONDA_DEFAULT_ENV"

### MAT domains in each window
for IDX in 1 2; do
    hmmsearch --domtblout "${RESULTS}/mat${IDX}_domains.txt" -E 1e-05 \
              "${MAT_HMM}" "${TMP}/mat${IDX}_frames.faa" > /dev/null 2>&1
done

### domain hits with HMM and frame coordinates
{
    echo ""
    printf "frame\tdomain\ti_evalue\thmm_from\thmm_to\tali_from\tali_to\n"
    grep -hv '^#' "${RESULTS}"/mat[12]_domains.txt | tr -s ' ' '\t' | cut -f1,4,13,16-19
} | tee -a "${REPORT}"


### =====================================================================
### step 5: divergence and conserved sequence among assemblies
### =====================================================================

echo ""
echo "=== Step 5: divergence among assemblies ==="

conda deactivate
conda activate cg_bioinfo_v1_blast
echo "Running in environment: $CONDA_DEFAULT_ENV"

### reference idiomorph against the window of every assembly of its type
printf "assembly\tidiomorph\tpident\taligned_bp\tquery_cov\n" > "${DIVERGENCE}"

for PAIR in "MAT1_1 MAT1-1 1" "MAT1_2 MAT1-2 2"; do

    set -- ${PAIR}
    NAME="$1"; TYPE="$2"; IDX="$3"

    mkdir -p "${TMP}/hits_${NAME}"

    for BASE in $(grep -P "\t${TYPE}$" "${SCAN}" | cut -f1); do

        W="${WIN_DIR}/${BASE}_mat${IDX}.fasta"
        [[ -s "${W}" ]] || continue

        blastn -task blastn -query "${TMP}/${NAME}.fasta" -subject "${W}" \
               -dust no -evalue 1e-10 -outfmt "6 qseqid qstart qend pident length qcovs" \
               2>/dev/null > "${TMP}/${BASE}_hits.tsv"

        ### hits as BED
        awk -v OFS='\t' '{print $1, $2 - 1, $3}' "${TMP}/${BASE}_hits.tsv" \
          > "${TMP}/hits_${NAME}/${BASE}.bed"

        ### best hit identity and length, coverage over all hits
        HIT=$(head -1 "${TMP}/${BASE}_hits.tsv" | cut -f4-6)
        [[ -n "${HIT}" ]] || HIT=$'NA\tNA\tNA'

        printf "%s\t%s\t%s\n" "${BASE}" "${NAME}" "${HIT}" >> "${DIVERGENCE}"
    done
done

{
    echo ""
    echo "lowest identity per idiomorph:"
    for NAME in MAT1_1 MAT1_2; do
        grep -P "\t${NAME}\t" "${DIVERGENCE}" | sort -t$'\t' -k3,3g | head -1
    done
    echo "lowest coverage per idiomorph:"
    for NAME in MAT1_1 MAT1_2; do
        grep -P "\t${NAME}\t" "${DIVERGENCE}" | sort -t$'\t' -k5,5g | head -1
    done
} | tee -a "${REPORT}"

conda deactivate
conda activate cg_bioinfo_v1_bedtools
echo "Running in environment: $CONDA_DEFAULT_ENV"

### reference sequence hit in every assembly of its type
: > "${CONSERVED}"

for NAME in MAT1_1 MAT1_2; do

    N_ASM=$(ls "${TMP}/hits_${NAME}"/*.bed | wc -l)

    for F in "${TMP}/hits_${NAME}"/*.bed; do
        bedtools sort -i "${F}" | bedtools merge -i - > "${F}.merged"
    done

    bedtools multiinter -i "${TMP}/hits_${NAME}"/*.bed.merged \
      | awk -v n="${N_ASM}" '$4 == n' | cut -f1-3 \
      | bedtools merge -i - >> "${CONSERVED}"
done


### =====================================================================
### step 6: core and primer sites
### =====================================================================

echo ""
echo "=== Step 6: core and primer sites ==="

conda deactivate
conda activate cg_bioinfo_v1_seqkit
echo "Running in environment: $CONDA_DEFAULT_ENV"

### primer sites in the reference
printf ">MAT1_F\n%s\n>MAT1_R\n%s\n>MAT2_F\n%s\n>MAT2_R\n%s\n" \
       "${MAT1_F}" "${MAT1_R}" "${MAT2_F}" "${MAT2_R}" > "${TMP}/primers.fasta"

seqkit locate -i -m "${PRIMER_MM}" -f "${TMP}/primers.fasta" "${REF}" > "${PRIMERS}"
seqkit locate -i -m "${PRIMER_MM}" -f "${TMP}/primers.fasta" --bed "${REF}" > "${PRIMERS_BED}"

## check if every primer sits once in the reference
N_SITES=$(wc -l < "${PRIMERS_BED}")
N_PRIMERS=$(cut -f4 "${PRIMERS_BED}" | sort -u | wc -l)

if [[ "${N_SITES}" -ne 4 || "${N_PRIMERS}" -ne 4 ]]; then
    echo "WARNING: ${N_SITES} sites for ${N_PRIMERS} primers in the reference," \
         "expected 4, see ${PRIMERS}" | tee -a "${REPORT}"
fi

conda deactivate
conda activate cg_bioinfo_v1_blast
echo "Running in environment: $CONDA_DEFAULT_ENV"

### idiomorph sequence found in genomes of the other mating type
: > "${TMP}/mask.bed"

for PAIR in "MAT1_1 MAT1-2" "MAT1_2 MAT1-1"; do

    set -- ${PAIR}
    NAME="$1"; OTHER="$2"

    for BASE in $(grep -P "\t${OTHER}$" "${SCAN}" | cut -f1); do

        makeblastdb -in "${ASM_DIR}/${BASE}.fasta" -dbtype nucl \
                    -out "${TMP}/blastdb/${BASE}" > /dev/null

        blastn -task blastn -query "${TMP}/${NAME}.fasta" -db "${TMP}/blastdb/${BASE}" \
               -dust no -evalue 1e-10 -perc_identity "${MASK_PIDENT}" \
               -num_threads "${THREADS}" \
               -outfmt "6 qseqid qstart qend" 2>/dev/null \
          | awk -v OFS='\t' '{print $1, $2 - 1, $3}' >> "${TMP}/mask.bed"
    done
done

conda deactivate
conda activate cg_bioinfo_v1_bedtools
echo "Running in environment: $CONDA_DEFAULT_ENV"

### mask: merged hits in the other mating type
bedtools sort -i "${TMP}/mask.bed" | bedtools merge -i - > "${MASK}"

### core: idiomorph without flanks and mask, present in all assemblies of its type
cut -f1,2 "${REF}.fai" > "${TMP}/ref.genome"

sed 's/\t/\t0\t/' "${TMP}/ref.genome" \
  | bedtools slop -i - -g "${TMP}/ref.genome" -b "-${FLANK}" \
  | bedtools subtract -a - -b "${MASK}" \
  | bedtools intersect -a - -b "${CONSERVED}" \
  > "${CORE}"

{
    echo ""
    echo "conserved bp per idiomorph:"
    bedtools summary -i "${CONSERVED}" -g "${TMP}/ref.genome" | cut -f1,4
    echo "core bp per idiomorph:"
    bedtools summary -i "${CORE}" -g "${TMP}/ref.genome" | cut -f1,4
} | tee -a "${REPORT}"

## check if both idiomorphs keep a core
if [[ "$(cut -f1 "${CORE}" | sort -u | wc -l)" -ne 2 ]]; then
    echo "ERROR: the core holds $(cut -f1 "${CORE}" | sort -u | wc -l) idiomorphs, expected 2"
    exit 1
fi

## check if the primers sit inside the core
bedtools intersect -a "${PRIMERS_BED}" -b "${CORE}" -f 1.0 -v \
  | cut -f4 | sed 's/^/WARNING: outside the core: /' | tee -a "${REPORT}"

### primer bp inside the core
{
    echo "primer bp inside the core:"
    bedtools intersect -a "${PRIMERS_BED}" -b "${CORE}" -wao | cut -f4,10 | sort -u
} | tee -a "${REPORT}"

rm -rf "${TMP}"


### output files (copy mat_primer_sites.tsv to local/01_variant_calling/data/103_data_mating_type_collect for local/01_variant_calling/103_mating_type.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  mat_idiomorphs.fasta        the two idiomorphs with ${FLANK} bp flanks"
echo "  mat_idiomorphs.fasta.fai    the index"
echo "  mat_core.bed                idiomorphs without flanks and mask, present in all assemblies of their type"
echo "  mat_mask.bed                idiomorph sequence found in the other mating type"
echo "  mat_conserved.bed           reference sequence found in all assemblies of its type"
echo "  mat_primer_sites.tsv        primer sites in the reference"
echo "  mat_primer_sites.bed        primer sites in the reference, as BED"
echo "  mat_assembly_divergence.tsv reference idiomorphs against every assembly of their type"
echo "  mat_reference_report.txt    which assembly gave what, and the lengths"
echo "  mat_assembly_scan.tsv       the idiomorph carried by every assembly"
echo "  mat_reference_choice.tsv    every candidate window"
echo "  mat1_region.fasta           the MAT1-1 window the reference was cut from"
echo "  mat2_region.fasta           the MAT1-2 window the reference was cut from"
echo "  mat1_shared.bed             MAT1-1 window sequence shared with MAT1-2"
echo "  mat2_shared.bed             MAT1-2 window sequence shared with MAT1-1"
echo "  mat_gene_annotation.tsv     curated MAT proteins against both windows"
echo "  mat1_domains.txt            Pfam hits in the MAT1-1 window"
echo "  mat2_domains.txt            Pfam hits in the MAT1-2 window"
echo "  windows/                    the window from every assembly"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Idiomorph reference built in ${RUNTIME} seconds."
