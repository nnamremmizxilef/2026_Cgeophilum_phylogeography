#!/bin/bash
#SBATCH --job-name=alignment_divergence    ### job name
#SBATCH --cpus-per-task=4                  ### cores
#SBATCH --mem-per-cpu=8G                   ### memory per core
#SBATCH --time=08:00:00                    ### walltime
#SBATCH --output /dev/null                 ### stdout (logged by tee)
#SBATCH --error /dev/null                  ### stderr (logged by tee)

### define paths
## check if ROOT is set
: "${ROOT:?set ROOT to the project folder}"
RESULTS="${ROOT}/02_structure_phylogeny/results/051_results_alignment_divergence"
LOGS="${ROOT}/02_structure_phylogeny/logs/051_logs_alignment_divergence"
### from hpc/02_structure_phylogeny/050_species_tree.sh
CONCAT="${ROOT}/02_structure_phylogeny/results/050_results_species_tree/cenococcum_orthologs_concatenated.fasta"
TIP_TABLE="${RESULTS}/tip_divergence.tsv"
COL_TABLE="${RESULTS}/callability_bins.tsv"

### parameters
### reference tip (every distance measured against it)
REF_TIP=consensus
### core columns (callable in at least this share of the population samples)
CORE_MIN=0.95
### tip groups by name
PAT_GENOME='^(Cenge|Pseflo)'
PAT_OUTGROUP='^M4-'

### create results folders
mkdir -p "${RESULTS}" "${LOGS}"

### create log file
LOGFILE="${LOGS}/${SLURM_JOB_ID:-local}_alignment_divergence.log"
exec > >(tee "${LOGFILE}") 2>&1

### print the script to the log
echo "########## Start of job script ##########"
cat "$0"
echo "########## End of job script ##########"

### conda (python with numpy for the column scan over 9.2 million sites)
eval "$(conda shell.bash hook)"
conda deactivate
conda activate cg_bioinfo_v1_python

START_TIME=$(date +%s)

### run info
echo "Date: $(date)"
echo "Running in environment: $CONDA_DEFAULT_ENV"
echo "python version: $(python3 --version)"
echo "numpy version: $(python3 -c 'import numpy; print(numpy.__version__)')"

## check if the concatenation is there
if [[ ! -s "${CONCAT}" ]]; then
    echo "ERROR: concatenation not found at ${CONCAT}"
    exit 1
fi

echo "Concatenation: ${CONCAT}"
echo "Sequences:     $(grep -c '^>' "${CONCAT}")"


### =====================================================================
### step 1: scan the concatenation
### =====================================================================

echo ""
echo "=== Step 1: Column and tip scan ==="

### two passes (per column: called population samples and those differing from the reference; per tip: distance to the reference over all comparable columns and over the core, and divergence at the columns kept and lost)
python3 - "${CONCAT}" "${REF_TIP}" "${CORE_MIN}" "${PAT_GENOME}" \
         "${PAT_OUTGROUP}" "${TIP_TABLE}" "${COL_TABLE}" <<'PYTHON'
import re
import sys
 
import numpy as np
 
concat, ref_tip, core_min, pat_genome, pat_outgroup, tip_out, col_out = sys.argv[1:]
core_min = float(core_min)
 
re_genome = re.compile(pat_genome)
re_outgroup = re.compile(pat_outgroup)
 
BASES = np.zeros(256, dtype=bool)
for b in b"ACGTacgt":
    BASES[b] = True
 
UPPER = np.arange(256, dtype=np.uint8)
for b in range(ord("a"), ord("z") + 1):
    UPPER[b] = b - 32
 
 
def group_of(name):
    if name == ref_tip:
        return "reference"
    if re_outgroup.search(name):
        return "outgroup"
    if re_genome.search(name):
        return "genome"
    return "sample"
 
 
def records(path):
    """yield (name, sequence as uint8 array), one record at a time"""
    name, chunks = None, []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    yield name, np.frombuffer(b"".join(chunks), dtype=np.uint8)
                name = line[1:].split()[0].decode()
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, np.frombuffer(b"".join(chunks), dtype=np.uint8)
 
 
# pass 1 ----------------------------------------------------------------------
 
ref = None
n_callable = None
n_alt = None
n_samples = 0
 
for name, seq in records(concat):
 
    seq = UPPER[seq]
 
    if n_callable is None:
        n_callable = np.zeros(seq.size, dtype=np.int32)
        n_alt = np.zeros(seq.size, dtype=np.int32)
 
    if name == ref_tip:
        ref = seq.copy()
 
    if group_of(name) != "sample":
        continue
 
    ok = BASES[seq]
    n_callable += ok
    n_samples += 1
 
if ref is None:
    sys.exit("reference tip %s not in the alignment" % ref_tip)
 
print("population samples: %d" % n_samples, flush=True)
print("columns:            %d" % ref.size, flush=True)
 
# the reference itself is not called everywhere, and a column without a
# reference base has nothing to be measured against
ref_ok = BASES[ref]
 
for name, seq in records(concat):
    if group_of(name) != "sample":
        continue
    seq = UPPER[seq]
    n_alt += BASES[seq] & ref_ok & (seq != ref)
 
callability = n_callable / max(n_samples, 1)
### the share of the samples that have a base and differ from the reference,
### which is the divergence a column carries
alt_frac = np.where(n_callable > 0, n_alt / np.maximum(n_callable, 1), np.nan)
 
core = (callability >= core_min) & ref_ok
print("core columns:       %d (%.2f%%)"
      % (core.sum(), 100 * core.sum() / ref.size), flush=True)
 
 
# pass 2 ----------------------------------------------------------------------
 
rows = []
 
for name, seq in records(concat):
 
    seq = UPPER[seq]
    ok = BASES[seq] & ref_ok
    diff = ok & (seq != ref)
 
    ok_core = ok & core
    diff_core = diff & core
 
    n_ok = int(ok.sum())
    n_ok_core = int(ok_core.sum())
 
    ### the divergence the rest of the panel sees at the columns this tip kept,
    ### and at the columns it lost. a tip masked at random has the same value
    ### twice, a tip masked where it is divergent has a higher value at the
    ### columns it lost
    lost = ref_ok & ~BASES[seq]
 
    rows.append((
        name,
        group_of(name),
        n_ok,
        100 * n_ok / ref.size,
        int(diff.sum()) / max(n_ok, 1),
        n_ok_core,
        int(diff_core.sum()) / max(n_ok_core, 1),
        float(np.nanmean(alt_frac[ok])) if n_ok else float("nan"),
        float(np.nanmean(alt_frac[lost])) if lost.any() else float("nan"),
    ))
 
with open(tip_out, "w") as fh:
    fh.write("tip\tgroup\tn_callable\tcallable_pct\tp_dist_all"
             "\tn_callable_core\tp_dist_core\talt_frac_kept\talt_frac_lost\n")
    for r in rows:
        fh.write("%s\t%s\t%d\t%.4f\t%.6g\t%d\t%.6g\t%.6g\t%.6g\n" % r)
 
 
# columns binned by callability ------------------------------------------------
 
### if masking were independent of divergence every bin would carry the same
### alternative allele fraction
edges = np.array([0, .5, .7, .8, .9, .95, .99, .999, 1.0001])
 
with open(col_out, "w") as fh:
    fh.write("callability_from\tcallability_to\tcolumns\tmean_alt_frac\n")
    for i in range(edges.size - 1):
        sel = (callability >= edges[i]) & (callability < edges[i + 1]) & ref_ok
        n = int(sel.sum())
        m = float(np.nanmean(alt_frac[sel])) if n else float("nan")
        fh.write("%.3f\t%.3f\t%d\t%.6g\n" % (edges[i], edges[i + 1], n, m))
PYTHON


### =====================================================================
### step 2: divergence by group
### =====================================================================

echo ""
echo "=== Step 2: Divergence by group ==="

awk -F'\t' 'NR > 1 {
        g = $2
        n[g]++
        cal[g] += $4
        all[g] += $5
        core[g] += $7
        kept[g] += $8
        lost[g] += $9
    }
    END {
        printf "%-10s %5s %10s %12s %12s %12s %12s\n",
               "group", "n", "callable%", "p_all", "p_core", "alt_kept", "alt_lost"
        for (g in n)
            printf "%-10s %5d %10.2f %12.6f %12.6f %12.6f %12.6f\n",
                   g, n[g], cal[g]/n[g], all[g]/n[g], core[g]/n[g],
                   kept[g]/n[g], lost[g]/n[g]
    }' "${TIP_TABLE}"


### =====================================================================
### step 3: divergence against column callability
### =====================================================================

echo ""
echo "=== Step 3: Divergence against column callability ==="

### alternative allele fraction per callability bin (higher in poorly callable bins if masking removes divergent sites)
awk -F'\t' 'NR == 1 {next}
    {printf "  %.3f to %.3f  %10d columns  alt frac %.6f\n", $1, $2, $3, $4}' \
    "${COL_TABLE}"

### output files (copy tip_divergence.tsv and callability_bins.tsv to local/02_structure_phylogeny/data/051_data_alignment_divergence for local/02_structure_phylogeny/051_allignment_divergence.R)
echo ""
echo "=== Output files in ${RESULTS} ==="
echo "  tip_divergence.tsv    per tip callability, distance to the reference and"
echo "                        the divergence at the columns it kept and lost"
echo "  callability_bins.tsv  divergence per column callability bin"

### runtime
END_TIME=$(date +%s)
RUNTIME=$((END_TIME - START_TIME))
echo ""
echo "Alignment divergence completed in ${RUNTIME} seconds ($((RUNTIME / 3600))h $(((RUNTIME % 3600) / 60))m)."
