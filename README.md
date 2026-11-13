# 2026_Cgeophilum_phylogeography

> **Work in progress.** Scripts are added phase by phase as the analyses are finalised and may still change.

Felix Zimmermann, Swiss Federal Institute for Forest, Snow and Landscape Research WSL

## Study

Population structure and phylogeography of the ectomycorrhizal fungus *Cenococcum geophilum* across Europe. Whole-genome short reads of several hundred isolates are mapped to a single-copy orthologue reference to call variants. The variants are then used to identify clones, build phylogenies, delimit cryptic lineages and describe the soil and climate niches of those lineages.

## Repository

```
hpc/<phase>/NNN_topic.sh     cluster scripts (SLURM)
local/<phase>/NNN_topic.R    R scripts
```

Phases (further phases to follow):

- `01_variant_calling`: read trimming and QC; single-copy orthologue reference (OrthoFinder); mapping, variant calling and filtering; clone identification; sample metadata and mating types.
- `02_structure_phylogeny`: call sets for phylogeny and clustering; population clustering (DAPC, sNMF, fineRADstructure); per-gene consensus sequences and alignments; concatenation (IQ-TREE), coalescent (ASTER) and site-based (CASTER-site) species trees; site and gene concordance factors; topology tests; monophyly of the clusters in each tree.
- `03_delimitation`: lineage labels at K = 4 to 6; unit trees (ASTER), Neighbour-Net, recombination (PHI, index of association), gene flow (Dsuite), mPTP and BPP (A00, A11, gdi); one table of evidence per unit.

Scripts run in numerical order within a phase; exceptions are named in the scripts.

- Cluster scripts: `export ROOT=/path/to/project` (the folder holding the phase folders), then `sbatch` or `bash` the script. Add account, QOS and partition at submission.
- R scripts: open `local/cenococcum_europe.Rproj`; paths are relative to `local/`.
- Each script names the files it needs and where they come from.

## License

[CC0 1.0 Universal](LICENSE)
