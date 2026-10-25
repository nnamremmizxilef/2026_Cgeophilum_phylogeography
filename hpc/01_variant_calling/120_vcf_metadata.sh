#!/bin/bash
#SBATCH --job-name=vcf_metadata     ### job name
#SBATCH --cpus-per-task=4           ### cores
#SBATCH --mem-per-cpu=4G            ### memory per core
#SBATCH --time=08:00:00             ### walltime
#SBATCH --output /dev/null          ### stdout (logged by tee)
#SBATCH --error /dev/null           ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/120_results_vcf_metadata"
LOGS="${ROOT}/01_variant_calling/logs/120_logs_vcf_metadata"

### VCFs in pipeline order, from hpc/01_variant_calling/070_variant_calling_consensus.sh, 080_filtering.sh, 091_decolned.sh, 092_decloned_filters.sh and 110_pruning_sco_consensus.sh
INPUT_DIRS=(
    "${ROOT}/01_variant_calling/results/070_results_variant_calling_consensus"
    "${ROOT}/01_variant_calling/results/080_results_filtering"
    "${ROOT}/01_variant_calling/results/091_results_decloned"
    "${ROOT}/01_variant_calling/results/092_results_decloned_filters"
    "${ROOT}/01_variant_calling/results/110_results_pruning_sco_consensus"
)

### parameters
### thresholds of hpc/01_variant_calling/080_filtering.sh, for the filter column (not read from the VCFs)
QUAL=30
MQ=50
GQ=20
DP=3
SITE_MISS=0.10
SAMPLE_MISS=0.50
MAF_V1=0.01
MAF_V2=0.05

### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_vcf_metadata.log"
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

echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"
echo "bcftools version: $(bcftools --version | head -1)"
echo ""


### =====================================================================
### function: describe (what a file is, read off its name)
### =====================================================================

### every step appends a suffix (the position of _decloned gives the order of clone correction and site filters)
describe () {

    local NAME=$1
    local core="${NAME}"
    local pruned=no poly=no declone=no
    local parent set_label maf maf_expr maf_order site_set clone filters

    ### the pruning and the polymorphic subset are always the last steps
    [[ "${core}" == *_1snp_per_sco ]] && { pruned=yes; core="${core%_1snp_per_sco}"; }
    [[ "${core}" == *_polymorphic  ]] && { poly=yes;   core="${core%_polymorphic}"; }

    ### peel the clone step and the MAF level off in whatever order they sit
    maf=none
    while : ; do
        case "${core}" in
            *_decloned) declone=yes; core="${core%_decloned}" ;;
            *_maf01)    maf="${MAF_V1}"; core="${core%_maf01}" ;;
            *_maf05)    maf="${MAF_V2}"; core="${core%_maf05}" ;;
            *)          break ;;
        esac
    done

    ### polymorphic subset as the MAF > 0 level of the same filter
    [[ "${poly}" == yes && "${maf}" == none ]] && maf="> 0"

    ### where the site filters were evaluated (_decloned_maf05 after clone correction, _maf05_decloned before)
    case "${NAME}" in
        *_decloned_maf0[15]*|*_decloned_polymorphic*)
            maf_order="after clone correction" ;;
        *_maf0[15]_decloned*)
            maf_order="before clone correction" ;;
        *)  [[ "${maf}" == none ]] && maf_order="none" || maf_order="no clone correction" ;;
    esac

    ### the threshold as it was written
    if   [[ "${maf}" == "> 0" ]]; then maf_expr="MAF>0"
    elif [[ "${maf}" != none  ]]; then maf_expr="MAF>=${maf}"
    else                               maf_expr=""
    fi

    case "${core}" in
        *_ingroup|*_ingroup_*)           set_label="ingroup" ;;
        *_plusoutgroup|*_plusoutgroup_*) set_label="ingroup_plus_outgroup" ;;
        *)                               set_label="all_called" ;;
    esac

    ### the parent is this file with its last step undone
    case "${NAME}" in
        *_1snp_per_sco)   parent="${NAME%_1snp_per_sco}" ;;
        consensus_sitefiltered_base_decloned_polymorphic)
                          parent="consensus_sitefiltered_base_decloned" ;;
        *_decloned_polymorphic)
                          parent="${NAME%_decloned_polymorphic}_noMAF_decloned" ;;
        *_polymorphic)    parent="${NAME%_polymorphic}" ;;
        *_decloned_maf01) parent="${NAME%_decloned_maf01}_noMAF_decloned" ;;
        *_decloned_maf05) parent="${NAME%_decloned_maf05}_noMAF_decloned" ;;
        *_decloned)       parent="${NAME%_decloned}" ;;
        consensus_raw)    parent="-" ;;
        consensus_sitefiltered_base) parent="consensus_raw" ;;
        *_noMAF)          parent="consensus_sitefiltered_base" ;;
        *_maf01)          parent="${NAME%_maf01}_noMAF" ;;
        *_maf05)          parent="${NAME%_maf05}_noMAF" ;;
    esac

    ### which sites the file holds
    if   [[ "${pruned}" == yes ]];         then site_set="one random SNP per SCO"
    elif [[ "${poly}"   == yes ]];         then site_set="polymorphic only"
    elif [[ "${core}" == consensus_raw ]]; then site_set="unfiltered"
    else                                        site_set="all filtered sites"
    fi

    [[ "${declone}" == yes ]] && clone="clonemates removed" || clone="none"

    ### the filters, in the order they were applied
    if [[ "${core}" == consensus_raw ]]; then
        filters="mpileup -q 5, no indels, ploidy 1, variant sites only"
    else
        filters="biallelic SNP; QUAL>=${QUAL}; MQ>=${MQ}; genotypes GQ>=${GQ} and DP>=${DP} set to missing"
        [[ "${core}" != consensus_sitefiltered_base ]] && \
            filters="${filters}; sample missingness<=${SAMPLE_MISS}; site F_MISSING<=${SITE_MISS}"

        if [[ "${maf_order}" == "before clone correction" ]]; then
            filters="${filters}; ${maf_expr}; clonemates removed, INFO recalculated"
        elif [[ "${maf_order}" == "after clone correction" ]]; then
            ### site filters re-evaluated on the retained samples
            filters="${filters}; clonemates removed, INFO recalculated"
            if [[ "${core}" == consensus_sitefiltered_base ]]; then
                filters="${filters}; ${maf_expr} re-evaluated on the retained samples"
            else
                filters="${filters}; ${maf_expr} and F_MISSING<=${SITE_MISS} re-evaluated on the retained samples"
            fi
        else
            [[ "${declone}" == yes ]] && filters="${filters}; clonemates removed, INFO recalculated"
            [[ -n "${maf_expr}" ]]    && filters="${filters}; ${maf_expr}"
        fi
    fi
    [[ "${pruned}" == yes ]] && filters="${filters}; one random SNP per contig"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s' \
        "${parent}" "${set_label}" "${clone}" "${site_set}" \
        "${maf}" "${maf_order}" "${filters}"
}


### =====================================================================
### step 1: one row per VCF
### =====================================================================

echo "=== Step 1: Reading the VCFs ==="
echo ""

TABLE="${RESULTS}/tmp_rows.tsv"
: > "${TABLE}"

for INPUT_DIR in "${INPUT_DIRS[@]}"; do

  STAGE=$(basename "${INPUT_DIR}" | cut -d_ -f1)

  for VCF_FILE in "${INPUT_DIR}"/*.vcf.gz; do

      NAME=$(basename "${VCF_FILE}" .vcf.gz)
      SITES="${RESULTS}/tmp_sites.txt"

      ### one pass over the file (invariant: no alternate allele carried or, at a biallelic site, no reference allele)
      bcftools query -f '%CHROM\t%POS\t%INFO/AC\t%INFO/AN\n' "${VCF_FILE}" > "${SITES}"

      read -r N_VARIANTS N_CONTIGS N_MONO N_MULTI <<<"$(awk -F'\t' '
      {
          n++
          contig[$1] = 1
          nalt = split($3, a, ",")
          if (nalt > 1) multi++
          if ($3 == "." || $4 == ".") next
          sum = 0
          for (i = 1; i <= nalt; i++) sum += a[i]
          if (sum == 0 || (nalt == 1 && sum == $4 + 0)) mono++
      }
      END {print n + 0, length(contig), mono + 0, multi + 0}
      ' "${SITES}")"

      SITES_MD5=$(cut -f1,2 "${SITES}" | md5sum | cut -d" " -f1)
      SAMPLES_MD5=$(bcftools query -l "${VCF_FILE}" | md5sum | cut -d" " -f1)
      N_SAMPLES=$(bcftools query -l "${VCF_FILE}" | wc -l)

      ### fields: name, stage, seven of describe, five counts, two checksums (fields 15 and 16)
      printf '%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n' \
          "${NAME}" "${STAGE}" "$(describe "${NAME}")" \
          "${N_SAMPLES}" "${N_VARIANTS}" "${N_CONTIGS}" "${N_MONO}" "${N_MULTI}" \
          "${SITES_MD5}" "${SAMPLES_MD5}" >> "${TABLE}"

      printf "  %-62s %8d samples  %10d sites  %7d invariant\n" \
          "${NAME}" "${N_SAMPLES}" "${N_VARIANTS}" "${N_MONO}"

      rm -f "${SITES}"

  done
done

echo ""
echo "VCFs read: $(wc -l < "${TABLE}")"


### =====================================================================
### step 2: which files share a site set and which share a sample set
### =====================================================================

echo ""
echo "=== Step 2: Grouping ==="

### site and sample groups from checksums of CHROM/POS and of the sample names (numbered in order of appearance)
OUT="${RESULTS}/vcf_metadata.tsv"

{
printf "file\tstage\tparent\tsample_set\tclone_correction\tsite_set\tmaf_level\tmaf_order\tfilters"
printf "\tn_samples\tn_variants\tn_contigs\tn_monomorphic\tn_multiallelic"
printf "\tsites_group\tsamples_group\tsites_md5\tsamples_md5\n"

awk -F'\t' 'BEGIN {OFS = "\t"}
{
    if (!($15 in sg)) sg[$15] = ++nsg
    if (!($16 in mg)) mg[$16] = ++nmg
    print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, \
          sg[$15], mg[$16], $15, $16
}' "${TABLE}"
} > "${OUT}"

rm -f "${TABLE}"

echo ""
echo "Site sets, files sharing one hold exactly the same sites:"
awk -F'\t' 'NR>1 {n[$15]++; f[$15] = f[$15] "\n      " $1}
END {for (g in n) if (n[g] > 1) print "  group " g " (" n[g] " files)" f[g]}' "${OUT}"

echo ""
echo "Sample sets:"
awk -F'\t' 'NR>1 {n[$16]++; s[$16] = $10}
END {for (g in n) printf "  group %s: %s samples in %d file(s)\n", g, s[g], n[g]}' "${OUT}" \
| sort -k3 -n

echo ""
echo "MAF order:"
awk -F'\t' 'NR>1 {print "  " $8}' "${OUT}" | sort | uniq -c | sort -rn


### =====================================================================
### step 3: the pipeline by stage
### =====================================================================

echo ""
echo "=== The pipeline, by stage ==="
printf "%-62s %8s %10s %9s %10s\n" "file" "samples" "sites" "invariant" "parent"
awk -F'\t' 'NR>1 {printf "%-62s %8d %10d %9d   %s\n", $1, $10, $11, $13, $3}' "${OUT}"


### output files
echo ""
echo "=== Output file in ${RESULTS} ==="
echo "  vcf_metadata.tsv"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Completed in ${RUNTIME} seconds ($((RUNTIME / 60)) min)."
