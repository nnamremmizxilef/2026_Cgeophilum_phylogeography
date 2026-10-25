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

Phases: `01_variant_calling` (further phases to follow).

- Cluster scripts: `export ROOT=/path/to/project` (the folder holding the phase folders), then `sbatch` or `bash` the script. Add account, QOS and partition at submission.
- R scripts: open `local/cenococcum_europe.Rproj`; paths are relative to `local/`.
- Each script names the files it needs and where they come from.

## License

[CC0 1.0 Universal](LICENSE)
