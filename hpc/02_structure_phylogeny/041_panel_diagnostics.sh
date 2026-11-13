#!/bin/bash
#SBATCH --job-name=panel_diagnostics    ### job name
#SBATCH --cpus-per-task=8               ### cores
#SBATCH --mem-per-cpu=4G                ### memory per core
#SBATCH --time=04:00:00                 ### walltime
#SBATCH --output /dev/null              ### stdout (logged by tee)
#SBATCH --error /dev/null               ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/02_structure_phylogeny/results/041_results_panel_diagnostics"
LOGS="${ROOT}/02_structure_phylogeny/logs/041_logs_panel_diagnostics"
### from hpc/02_structure_phylogeny/020_coverage_mask.sh
BED_DIR="${ROOT}/02_structure_phylogeny/results/020_results_coverage_mask/per_sample_bed"
### from hpc/02_structure_phylogeny/040_merge_alignments.sh
FINAL_DIR="${ROOT}/02_structure_phylogeny/results/040_results_merge_alignments/per_og_final"
MASK_TABLE="${RESULTS}/mask_fraction.tsv"
GENE_TABLE="${RESULTS}/alignment_stats.tsv"
TIP_TABLE="${RESULTS}/tip_missingness.tsv"
EMPTY_SAMPLE="${RESULTS}/empty_genes_per_sample.tsv"
EMPTY_GENE="${RESULTS}/empty_tips_per_gene.tsv"

### parameters
### cores (also outside SLURM)
THREADS="${SLURM_CPUS_PER_TASK:-8}"
### length of the consensus reference in bp (denominator of the masked fraction; sum of the lengths in consensus_sco_reference.fasta.fai of hpc/01_variant_calling/053_indexing.sh)
REF_SIZE=9229711

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_panel_diagnostics.log"
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
### step 1: masked fraction per sample
### =====================================================================

echo ""
echo "=== Step 1: Coverage mask ==="

for F in "${BED_DIR}"/*_lowcov.bed; do
    awk -v s="$(basename "${F}" _lowcov.bed)" -v r="${REF_SIZE}" \
        '{n += $3 - $2} END {printf "%s\t%.4f\n", s, (n + 0) / r * 100}' "${F}"
done | sort -k2 -g > "${MASK_TABLE}.tmp"

{ echo -e "sample\tmasked"; cat "${MASK_TABLE}.tmp"; } > "${MASK_TABLE}"
rm -f "${MASK_TABLE}.tmp"

echo "Samples: $(($(wc -l < "${MASK_TABLE}") - 1))"
echo "Median:  $(awk 'NR > 1 {v[NR - 1] = $2} END {printf "%.2f%%", v[int(NR / 2)]}' "${MASK_TABLE}")"
echo "Maximum: $(tail -1 "${MASK_TABLE}" | awk '{printf "%.2f%% (%s)", $2, $1}')"


### =====================================================================
### step 2: alignment length per orthogroup
### =====================================================================

echo ""
echo "=== Step 2: Alignment lengths ==="

seqkit stats -T -j "${THREADS}" "${FINAL_DIR}"/*_final_alignment.fasta \
    > "${GENE_TABLE}"

N_GENES=$(($(wc -l < "${GENE_TABLE}") - 1))

echo "Orthogroups: ${N_GENES}"
echo "Columns:     $(awk 'NR > 1 {n += $6} END {print n}' "${GENE_TABLE}")"


### =====================================================================
### step 3: missing data per tip
### =====================================================================

echo ""
echo "=== Step 3: Tip missingness ==="

### N content per tip over the concatenation (length and N percent per sequence, summed over the genes)
seqkit fx2tab -n -i -l -B N "${FINAL_DIR}"/*_final_alignment.fasta \
| awk 'BEGIN {OFS = "\t"}
    {
        sub(/^OG[0-9]+_/, "", $1)
        len[$1] += $2
        n[$1]   += $2 * $3 / 100
    }
    END {
        print "sample", "columns", "n_count", "n_percent"
        for (s in len) print s, len[s], int(n[s] + 0.5), 100 * n[s] / len[s]
    }' > "${TIP_TABLE}"

echo "Tips:    $(($(wc -l < "${TIP_TABLE}") - 1))"
echo "Median:  $(awk 'NR > 1 {print $4}' "${TIP_TABLE}" | sort -g | awk '{v[NR] = $1} END {printf "%.2f%%", v[int((NR + 1) / 2)]}')"
echo "Maximum: $(awk 'NR > 1 {print $4, $1}' "${TIP_TABLE}" | sort -gr | head -1 | awk '{printf "%.2f%% (%s)", $1, $2}')"


### =====================================================================
### step 4: genes with no called bases per tip
### =====================================================================

echo ""
echo "=== Step 4: Empty gene sequences ==="

### A, C, G and T content per sequence (a tip without any has no data for that gene, and iqtree2 refuses such an alignment)
EMPTY_RAW="${RESULTS}/.empty_scan.tmp"

seqkit fx2tab -n -i -B ACGTacgt "${FINAL_DIR}"/*_final_alignment.fasta \
    > "${EMPTY_RAW}"

### per sample, over all genes
awk -v total="${N_GENES}" 'BEGIN {OFS = "\t"}
    {
        s = $1
        sub(/^OG[0-9]+_/, "", s)
        seen[s]
        if ($2 <= 0.001) empty[s]++
    }
    END {
        print "sample", "empty_genes", "total_genes", "percent"
        for (s in seen) print s, empty[s] + 0, total, 100 * (empty[s] + 0) / total
    }' "${EMPTY_RAW}" > "${EMPTY_SAMPLE}.tmp"

{ head -1 "${EMPTY_SAMPLE}.tmp"; tail -n +2 "${EMPTY_SAMPLE}.tmp" | sort -k2 -gr; } \
    > "${EMPTY_SAMPLE}"
rm -f "${EMPTY_SAMPLE}.tmp"

### per gene, the tips dropped before the gene tree
awk 'BEGIN {OFS = "\t"}
    {
        og = $1
        sub(/_.*$/, "", og)
        seen[og]
        if ($2 <= 0.001) empty[og]++
    }
    END {
        print "orthogroup", "empty_tips"
        for (g in seen) print g, empty[g] + 0
    }' "${EMPTY_RAW}" > "${EMPTY_GENE}.tmp"

{ head -1 "${EMPTY_GENE}.tmp"; tail -n +2 "${EMPTY_GENE}.tmp" | sort -k2 -gr; } \
    > "${EMPTY_GENE}"
rm -f "${EMPTY_GENE}.tmp" "${EMPTY_RAW}"

echo "Genes with at least one empty tip: $(awk 'NR > 1 && $2 > 0' "${EMPTY_GENE}" | wc -l) of ${N_GENES}"
echo "Tips with at least one empty gene:  $(awk 'NR > 1 && $2 > 0' "${EMPTY_SAMPLE}" | wc -l)"
echo "Worst tips:"
awk 'NR > 1 && NR <= 11 {printf "  %-16s %5d genes (%.2f%%)\n", $1, $2, $4}' "${EMPTY_SAMPLE}"

### output files (copy the five tables to local/02_structure_phylogeny/data/041_data_panel_diagnostics for local/02_structure_phylogeny/041_panel_diagnostics.R, and tip_missingness.tsv also to local/02_structure_phylogeny/data/050_data_species_tree and local/02_structure_phylogeny/data/052_data_species_tree_decloned for local/02_structure_phylogeny/050_species_tree.R and 052_species_tree_decloned.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  mask_fraction.tsv            masked fraction of the reference per sample"
echo "  alignment_stats.tsv          seqkit stats per orthogroup alignment"
echo "  tip_missingness.tsv          N content per tip over all genes"
echo "  empty_genes_per_sample.tsv   genes with no called bases, per tip"
echo "  empty_tips_per_gene.tsv      tips with no called bases, per orthogroup"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Panel diagnostics completed in ${RUNTIME} seconds ($((RUNTIME / 60)) min)."
