#!/bin/bash
#SBATCH --job-name=snp_density      ### job name
#SBATCH --cpus-per-task=4           ### cores
#SBATCH --mem-per-cpu=4G            ### memory per core
#SBATCH --time=04:00:00             ### walltime
#SBATCH --output /dev/null          ### stdout (logged by tee)
#SBATCH --error /dev/null           ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/01_variant_calling/results/130_results_snp_density"
LOGS="${ROOT}/01_variant_calling/logs/130_logs_snp_density"

### from hpc/01_variant_calling/080_filtering.sh
FILTERED="${ROOT}/01_variant_calling/results/080_results_filtering"
### plink.genome, from hpc/01_variant_calling/090_clone_id.sh
CLONE="${ROOT}/01_variant_calling/results/090_results_clone_id"
### consensus reference from hpc/01_variant_calling/052_consensus_build.sh (index from 053_indexing.sh)
REFERENCE="${ROOT}/01_variant_calling/results/050_results_sco_consensus/consensus_sco_reference.fasta"

### ingroup set without MAF filter (the set clone detection ran on), from hpc/01_variant_calling/080_filtering.sh
VCF="${FILTERED}/consensus_ingroup_noMAF.vcf.gz"

### parameters
### clone threshold (clone_dist printed by local/01_variant_calling/090_clone_id.R)
THRESHOLD=7.117085267241543e-4

### extra pairs for the difference counts (step 2 adds the two pairs bracketing the threshold)
PAIRS=(
    "DE-REX-W1B:DE-REX-W1D"   ### the least divergent clonemate pair, 19 sites
    "FR-REX-N2B:DE-MAX-E1G"   ### two representatives from different countries
)

### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-4}"

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-$(date +%Y%m%d_%H%M%S)}_snp_density.log"
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
echo "Input: ${VCF}"
echo "Clone threshold: ${THRESHOLD}"
echo ""


### =====================================================================
### step 0: the distribution of pairwise distances
### =====================================================================

### split into within site, within country and between countries (gawk for asort)
echo "=== Step 0: Pairwise distances ==="

gawk -v thr="${THRESHOLD}" '
NR==1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
{
    d = $(c["IBS0"]) / ($(c["IBS0"]) + $(c["IBS1"]) + $(c["IBS2"]))
    n++; v[n] = d; tot += d
    if (d <= thr) below++
    split($(c["IID1"]), a, "-"); split($(c["IID2"]), b, "-")
    if (a[1] == b[1] && a[2] == b[2]) {ns++; ss += d}
    else if (a[1] == b[1])            {nc++; sc += d}
    else                              {nb++; sb += d}
}
END {
    asort(v)
    printf "  pairs                          %d\n", n
    printf "  mean                           %.4e\n", tot / n
    printf "  median                         %.4e\n", v[int(n * 0.50)]
    printf "  first quartile                 %.4e\n", v[int(n * 0.25)]
    printf "  third quartile                 %.4e\n", v[int(n * 0.75)]
    printf "  min                            %.4e\n", v[1]
    printf "  max                            %.4e\n", v[n]
    printf "  pairs at or below the threshold %d  (%.2f%%)\n", below + 0, below / n * 100
    printf "  threshold below the mean       %.0f fold\n", (tot / n) / thr
    printf "  threshold below the median     %.0f fold\n", v[int(n * 0.50)] / thr
    printf "\n  mean within site               %.4e  (%d pairs)\n", ss / ns, ns
    printf "  mean within country            %.4e  (%d pairs)\n", sc / nc, nc
    printf "  mean between countries         %.4e  (%d pairs)\n", sb / nb, nb
}' "${CLONE}/plink.genome"

echo ""


### =====================================================================
### step 1: SNPs per orthologue
### =====================================================================

echo "=== Step 1: SNPs per orthologue ==="

### SNPs per contig (VCF sorted, so uniq -c counts in one pass; orthologues without a variant added from the reference index)
PER_SCO="${RESULTS}/snps_per_sco.tsv"

bcftools query -f '%CHROM\n' "${VCF}" | uniq -c \
    | awk '{print $2 "\t" $1}' | sort -k1,1 > "${RESULTS}/tmp_counts.tsv"

sort -k1,1 "${REFERENCE}.fai" | cut -f1,2 > "${RESULTS}/tmp_lengths.tsv"

printf "sco\tlength\tn_snps\tsnps_per_kb\n" > "${PER_SCO}"
join -a1 -e0 -o '0,1.2,2.2' -t$'\t' \
    "${RESULTS}/tmp_lengths.tsv" "${RESULTS}/tmp_counts.tsv" \
| awk 'BEGIN {OFS = "\t"} {printf "%s\t%d\t%d\t%.2f\n", $1, $2, $3, $3 / ($2 / 1000)}' \
>> "${PER_SCO}"

rm -f "${RESULTS}/tmp_counts.tsv" "${RESULTS}/tmp_lengths.tsv"

### distribution of SNP density and its concentration (gawk for asort)
gawk -F'\t' -v thr="${THRESHOLD}" '
NR>1 {n++; len += $2; snp += $3; d[n] = $4 + 0; s[n] = $3 + 0; if ($3 == 0) zero++}
END {
    asort(d); asort(s)
    printf "  orthologues                    %d\n", n
    printf "  total length                   %d bp\n", len
    printf "  total SNPs                     %d\n", snp
    printf "  genome wide density            %.1f per kb\n", snp / (len / 1000)
    printf "  orthologues without a SNP      %d\n", zero + 0
    printf "  density  min                   %.1f per kb\n", d[1]
    printf "  density  first quartile        %.1f per kb\n", d[int(n * 0.25)]
    printf "  density  median                %.1f per kb\n", d[int(n * 0.50)]
    printf "  density  third quartile        %.1f per kb\n", d[int(n * 0.75)]
    printf "  density  max                   %.1f per kb\n", d[n]
    ### SNPs carried by the most variable genes
    top5 = 0; for (i = int(n * 0.95) + 1; i <= n; i++) top5 += s[i]
    top1 = 0; for (i = int(n * 0.99) + 1; i <= n; i++) top1 += s[i]
    printf "  SNPs in the top 5%% of genes    %d  (%.1f%% of all)\n", top5, top5 / snp * 100
    printf "  SNPs in the top 1%% of genes    %d  (%.1f%% of all)\n", top1, top1 / snp * 100
    ### threshold for one gene of mean length
    printf "\n  mean orthologue length         %.0f bp\n", len / n
    printf "  sites differing at the threshold %.2f (gene of mean length)\n", \
        thr * (len / n)
}' "${PER_SCO}"

echo ""
echo "Written: $(basename "${PER_SCO}")"
echo ""


### =====================================================================
### step 2: the two pairs that bracket the threshold
### =====================================================================

echo "=== Step 2: The pairs closest to the threshold ==="

### closest pair above the threshold and most distant pair below it
read -r BELOW ABOVE <<<"$(awk -v thr="${THRESHOLD}" '
NR==1 {for (i = 1; i <= NF; i++) c[$i] = i; next}
{
    d = $(c["IBS0"]) / ($(c["IBS0"]) + $(c["IBS1"]) + $(c["IBS2"]))
    if (d <= thr && d > lo) {lo = d; lop = $(c["IID1"]) ":" $(c["IID2"])}
    if (d >  thr && (hi == 0 || d < hi)) {hi = d; hip = $(c["IID1"]) ":" $(c["IID2"])}
}
END {print lop, hip}' "${CLONE}/plink.genome")"

echo "  most distant pair still called clonal: ${BELOW}"
echo "  closest pair called distinct:          ${ABOVE}"
echo ""

PAIRS+=("${BELOW}" "${ABOVE}")


### =====================================================================
### step 3: differences per orthologue, for each pair
### =====================================================================

echo "=== Step 3: Differences per orthologue ==="
echo ""

PAIR_TABLE="${RESULTS}/pair_differences_per_sco.tsv"
printf "pair\tsco\tn_diff\n" > "${PAIR_TABLE}"

printf "  %-30s %10s %10s %12s %8s %8s\n" \
    "pair" "compared" "differing" "distance" "SCOs" "max/SCO"

for PAIR in "${PAIRS[@]}"; do

    A="${PAIR%%:*}"
    B="${PAIR##*:}"

    ### one pass per pair (compared: both called; differing: both called and unequal)
    bcftools view --threads "${THREADS}" -s "${A},${B}" "${VCF}" -Ou \
    | bcftools query -f '%CHROM[\t%GT]\n' \
    | awk -v pair="${A}:${B}" -v out="${PAIR_TABLE}" 'BEGIN {OFS = "\t"}
    {
        if ($2 == "." || $3 == ".") next
        cmp++
        if ($2 != $3) {diff++; per[$1]++}
    }
    END {
        nsco = 0; mx = 0
        for (s in per) {nsco++; if (per[s] > mx) mx = per[s]; print pair, s, per[s] >> out}
        printf "  %-30s %10d %10d %12.3e %8d %8d\n", \
            pair, cmp, diff + 0, (cmp > 0 ? diff / cmp : 0), nsco, mx
    }'

done

echo ""
echo "Written: $(basename "${PAIR_TABLE}")"


### output files (copy snps_per_sco.tsv and pair_differences_per_sco.tsv to local/01_variant_calling/data/130_data_snp_density for local/01_variant_calling/130_snp_density.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  snps_per_sco.tsv               one row per orthologue"
echo "  pair_differences_per_sco.tsv   one row per orthologue and pair"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Completed in ${RUNTIME} seconds ($((RUNTIME / 60)) min)."
