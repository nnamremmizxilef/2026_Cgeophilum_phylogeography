## ============================================================================
## 001_lineage_labels.R
## Lineage labels at K = 4, 5 and 6 in one table, named by nesting and frozen
## with a checksum
##
## Files needed:
##   lineage_assignment_K4.tsv        DAPC lineage of every ingroup sample at K = 4
##                                    (clonemates inherit their MLG)
##   lineage_assignment_K5.tsv        the same at K = 5
##   lineage_assignment_K6.tsv        the same at K = 6
##                                    (all three from local/02_structure_phylogeny/014_clustering.R,
##                                     results/014_results_clustering/summary)
##   sample_metadata_ALL_FIXED.tsv    per-sample table (in_ingroup, mlg, clone_status, project,
##                                    country, site), from local/01_variant_calling/103_mating_type.R
## ============================================================================

library(tidyverse)



# files ------------------------------------------------------------------------

clusters_root <- "02_structure_phylogeny/results/014_results_clustering/summary"
meta_file     <- "01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv"
out_root      <- "03_delimitation/results/001_results_lineage_labels"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### candidate levels, coarsest first (the coarsest names the lineages)
k_levels <- c(4, 5, 6)

### expected number of ingroup samples
n_ingroup_expected <- 585



# read -------------------------------------------------------------------------

### lineage assignment per level from 014
assign <- map(set_names(k_levels), function(k) {
  read_tsv(file.path(clusters_root, sprintf("lineage_assignment_K%d.tsv", k)),
           col_types = cols(.default = "c", clustered = "l"))
})

### per-sample table, ingroup only
meta <- read_tsv(meta_file, col_types = cols(.default = "c"))

## check if the metadata carries the columns used here
need <- c("sample", "in_ingroup", "mlg", "clone_status", "project", "country", "site")

if (!all(need %in% names(meta))) {
  stop(basename(meta_file), " lacks: ", paste(setdiff(need, names(meta)), collapse = ", "))
}

meta <- filter(meta, in_ingroup == "TRUE")



# unit names -------------------------------------------------------------------

### a cluster takes the name of the unit above it, plus a letter by size where
### that unit splits
name_level <- function(child, parent) {
  key <- child %>%
    inner_join(parent, by = "sample") %>%
    count(cluster, parent = name) %>%
    slice_max(n, n = 1, by = cluster, with_ties = FALSE) %>%
    select(cluster, parent) %>%
    left_join(count(child, cluster, name = "size"), by = "cluster") %>%
    arrange(parent, desc(size), cluster) %>%
    mutate(name = if (n() == 1) parent else paste0(parent, letters[row_number()]),
           .by = parent)
  child %>%
    left_join(select(key, cluster, name), by = "cluster") %>%
    select(sample, cluster, name)
}

### coarsest level keeps the DAPC numbering (named on the clustered samples)
k_first <- as.character(k_levels[1])

units <- list()
units[[k_first]] <- assign[[k_first]] %>%
  filter(clustered) %>%
  transmute(sample, cluster = lineage, name = paste0("L", lineage))

for (i in seq_along(k_levels)[-1]) {
  k  <- as.character(k_levels[i])
  kp <- as.character(k_levels[i - 1])
  units[[k]] <- name_level(assign[[k]] %>%
                             filter(clustered) %>%
                             select(sample, cluster = lineage),
                           select(units[[kp]], sample, name))
}

### DAPC cluster number and unit name per level
unit_key <- imap_dfr(units, function(u, k) {
  u %>%
    distinct(cluster, name) %>%
    mutate(k = as.integer(k), .before = 1)
}) %>%
  rename(dapc = cluster, unit = name)



# labels -----------------------------------------------------------------------

### every ingroup sample, named through its level's DAPC cluster
labels <- reduce(imap(assign, function(a, kk) {
  a %>%
    left_join(unit_key %>% filter(k == as.integer(kk)) %>% select(dapc, unit),
              by = c("lineage" = "dapc")) %>%
    select(sample, !!paste0("k", kk) := unit)
}), left_join, by = "sample")

labels <- meta %>%
  select(sample, project, country, site, mlg, clone_status) %>%
  left_join(select(assign[[k_first]], sample, clustered), by = "sample") %>%
  left_join(labels, by = "sample") %>%
  arrange(sample)

### size per unit and level
unit_sizes <- map_dfr(k_levels, function(k) {
  col <- paste0("k", k)
  labels %>%
    summarise(n_mlg          = sum(clustered %in% TRUE),
              n_samples      = n(),
              pct_magic_mlg  = round(100 * sum(clustered %in% TRUE & project == "magic") /
                                       sum(clustered %in% TRUE), 1),
              n_sites        = n_distinct(site),
              .by = all_of(col)) %>%
    rename(unit = all_of(col)) %>%
    mutate(k = k, .before = 1)
}) %>%
  left_join(unit_key, by = c("k", "unit")) %>%
  relocate(dapc, .after = unit) %>%
  arrange(k, unit)



# checks -----------------------------------------------------------------------

## check if the three tables hold the same samples
if (!all(map_lgl(assign, ~ setequal(.x$sample, assign[[1]]$sample)))) {
  warning("the lineage assignment tables hold different samples")
}

## check if 014 labelled the ingroup of the metadata
if (!setequal(assign[[k_first]]$sample, meta$sample)) {
  warning(length(setdiff(meta$sample, assign[[k_first]]$sample)),
          " metadata ingroup samples not in 014, ",
          length(setdiff(assign[[k_first]]$sample, meta$sample)),
          " 014 samples not in the metadata ingroup")
}

## check if the clustered samples are the MLG representatives of the metadata
n_status <- sum((labels$clustered %in% TRUE) !=
                  (labels$clone_status %in% c("Unique", "Representative")))

if (n_status > 0) {
  warning(n_status, " samples differ between clustered (014) and clone_status (metadata)")
}

## check if every sample has a lineage at every level
if (anyNA(select(labels, matches("^k\\d+$")))) {
  warning(sum(!complete.cases(select(labels, matches("^k\\d+$")))),
          " samples lack a label at some level")
}

## check if every cluster sits inside one unit of the level above
crossing <- map_int(seq_along(k_levels)[-1], function(i) {
  x <- labels[[paste0("k", k_levels[i])]]
  y <- labels[[paste0("k", k_levels[i - 1])]]
  sum(!str_starts(x, fixed(y)), na.rm = TRUE)
})

if (any(crossing > 0)) {
  warning("samples cross between levels: ", paste(crossing, collapse = ", "))
}

## check if the table holds the expected number of samples
if (nrow(labels) != n_ingroup_expected) {
  warning(nrow(labels), " samples labelled, expected ", n_ingroup_expected)
}

## check if no sample appears twice
if (anyDuplicated(labels$sample) > 0) {
  warning("samples labelled twice: ",
          paste(unique(labels$sample[duplicated(labels$sample)]), collapse = ", "))
}



# tables -----------------------------------------------------------------------

### frozen table and its checksum (md5sum -c format)
labels_file <- file.path(out_root, "lineage_labels.tsv")
write_tsv(labels, labels_file, na = "NA")
write_lines(sprintf("%s  %s", unname(tools::md5sum(labels_file)), basename(labels_file)),
            paste0(labels_file, ".md5"))

### unit names against the DAPC numbering, with sizes
write_tsv(unit_sizes, file.path(out_root, "unit_sizes.tsv"))

cat("\n=== Lineage labels ===\n")
cat(sprintf("samples:  %d\n", nrow(labels)))
cat(sprintf("MLGs:     %d\n", sum(labels$clustered %in% TRUE)))
cat(sprintf("md5:      %s\n", unname(tools::md5sum(labels_file))))

cat("\n=== Units per level ===\n")
print(unit_sizes, n = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
### copy lineage_labels.tsv and lineage_labels.tsv.md5 to ${ROOT}/03_delimitation/supplementary
### for hpc/02_structure_phylogeny/052_species_tree_decloned.sh and hpc/03_delimitation/010,
### 020, 021, 030, 040 and 050
cat("  lineage_labels.tsv\n")
cat("  lineage_labels.tsv.md5\n")
cat("  unit_sizes.tsv\n")
