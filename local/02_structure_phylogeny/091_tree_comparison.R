## ============================================================================
## 091_tree_comparison.R
## Distances, clade agreement and AU test of the four topology estimates
##
## Files needed:
##   candidate_topologies.treefile  concat, aster_collapsed, aster and caster on shared tips
##   topology_distances.rfdist      pairwise Robinson-Foulds distances
##   topology_test.iqtree           AU test report
##                                  (all three from hpc/02_structure_phylogeny/091_tree_comparison.sh;
##                                   copy them from ${ROOT}/02_structure_phylogeny/results/091_results_tree_comparison
##                                   to local/02_structure_phylogeny/data/091_data_tree_comparison)
##   sample_metadata.tsv            sample metadata, from hpc/01_variant_calling/100_metadata.sh
##                                  (copy ${ROOT}/01_variant_calling/results/100_results_metadata/sample_metadata.tsv
##                                   to local/02_structure_phylogeny/data/091_data_tree_comparison)
## ============================================================================

library(ape)
library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

data_root <- "02_structure_phylogeny/data/091_data_tree_comparison"
out_root  <- "02_structure_phylogeny/results/091_results_tree_comparison"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

f_trees <- file.path(data_root, "candidate_topologies.treefile")
f_rf    <- file.path(data_root, "topology_distances.rfdist")
f_test  <- file.path(data_root, "topology_test.iqtree")
f_meta  <- file.path(data_root, "sample_metadata.tsv")



# parameters -------------------------------------------------------------------

### topologies in the order of the tree file
tree_names <- c("concat", "aster_collapsed", "aster", "caster")

### outgroup (assemblies and mapped samples; the assemblies alone are not a clade)
outgroup <- c("PsefloM402_1", "PsefloM405_1", "PsefloM406_1",
              "M4-01", "M4-02", "M4-03", "M4-04", "M4-05", "M4-06")

### assemblies (tip name pattern)
genome_pattern <- "^(Cenge|Pseflo)"

### AU p value below which a topology is rejected
au_alpha <- 0.05

### smallest clade counted as large (shallow clades counted separately)
clade_min <- 5

### clades printed
n_show <- 25



# style ------------------------------------------------------------------------

### topology colours
col_tree <- c(concat = "#EE6677", aster_collapsed = "#4477AA", aster = "#228833",
              caster = "#AA3377")

### clade agreement (number of topologies containing the clade)
col_agree <- c("4 of 4" = "#228833", "3 of 4" = "#CCBB44", "2 of 4" = "#EE7733",
               "1 of 4" = "#EE6677")

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

trees <- read.tree(f_trees)

## check if all four topologies were read
if (length(trees) != length(tree_names)) {
  stop("expected ", length(tree_names), " topologies, found ", length(trees))
}

names(trees) <- tree_names

meta <- read_tsv(f_meta, show_col_types = FALSE)

### all four rooted on the outgroup (via an ingroup tip first, root() misplaces
### the root when the outgroup spans the stored root)
trees <- lapply(trees, function(tr) {
  og <- intersect(outgroup, tr$tip.label)
  if (length(og) > 0) {
    tr <- root(tr, outgroup = setdiff(tr$tip.label, og)[1],
               resolve.root = TRUE, edgelabel = TRUE)
    tr <- root(tr, outgroup = og, resolve.root = TRUE, edgelabel = TRUE)
  }
  ladderize(tr)
})

n_tips <- length(trees[[1]]$tip.label)
n_ref  <- length(trees)

message("topologies: ", length(trees), "   shared tips: ", n_tips)



# checks -----------------------------------------------------------------------

## check if the four topologies share one tip set
if (!all(vapply(trees, function(tr) setequal(tr$tip.label, trees[[1]]$tip.label),
                logical(1)))) {
  stop("the topologies do not share one tip set")
}

## check if the outgroup is monophyletic in each topology
for (n in names(trees)) {
  og <- intersect(outgroup, trees[[n]]$tip.label)
  if (length(og) > 1 && !is.monophyletic(trees[[n]], og)) {
    warning("the outgroup is not monophyletic in ", n)
  }
}



# distances --------------------------------------------------------------------

### Robinson-Foulds distances, raw and normalised (1: no clade in common)
pairs <- combn(names(trees), 2, simplify = FALSE)

max_rf <- 2 * (n_tips - 3)

distances <- map_dfr(pairs, function(p) {
  ### on the unrooted trees (the normalisation below is the unrooted maximum)
  rf <- as.numeric(dist.topo(unroot(trees[[p[1]]]), unroot(trees[[p[2]]]),
                             method = "PH85"))
  tibble(tree_a        = p[1],
         tree_b        = p[2],
         rf_distance   = rf,
         rf_max        = max_rf,
         rf_normalised = rf / max_rf)
})

write_tsv(distances, file.path(out_root, "topology_distances.tsv"))

cat("\n=== pairwise distances ===\n")
print(distances, n = Inf)

p_dist <- distances %>%
  mutate(pair = paste(tree_a, "vs", tree_b)) %>%
  ggplot(aes(x = reorder(pair, rf_normalised), y = rf_normalised)) +
  geom_col(fill = "grey40", width = 0.6) +
  coord_flip() +
  labs(x = NULL, y = "Normalised Robinson-Foulds distance")

save_fig("topology_distances.pdf", p_dist, 7.1, 4)



# clades -----------------------------------------------------------------------

### clades of each topology as sorted tip sets (independent of node numbering)
clades_of <- function(tr) {
  pp   <- prop.part(tr)
  labs <- attr(pp, "labels")
  tibble(clade = vapply(pp, function(i) paste(sort(labs[i]), collapse = "|"),
                        character(1)),
         size  = vapply(pp, length, integer(1)))
}

clades <- imap_dfr(trees, ~ mutate(clades_of(.x), tree = .y))

### one row per clade, with which topologies contain it
clade_table <- clades %>%
  distinct(clade, size, tree) %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = tree, values_from = present, values_fill = FALSE) %>%
  mutate(n_trees = rowSums(across(all_of(names(trees)))),
         support = factor(sprintf("%d of %d", n_trees, n_ref),
                          levels = sprintf("%d of %d", n_ref:1, n_ref)))

write_tsv(arrange(clade_table, desc(size)),
          file.path(out_root, "clade_presence.tsv"))

clade_summary <- clade_table %>%
  group_by(support) %>%
  summarise(clades       = n(),
            large_clades = sum(size >= clade_min),
            median_size  = median(size),
            max_size     = max(size),
            .groups      = "drop")

write_tsv(clade_summary, file.path(out_root, "clade_summary.tsv"))

cat("\n=== clades by how many topologies contain them ===\n")
print(clade_summary, n = Inf)

### clade agreement among the large clades
p_clades <- clade_table %>%
  filter(size >= clade_min) %>%
  ggplot(aes(x = size, fill = support)) +
  geom_histogram(bins = 40, position = "stack") +
  scale_x_log10() +
  scale_fill_manual(values = col_agree, name = "Topologies with the clade", drop = FALSE) +
  labs(x = "N tips in the clade", y = "N clades")

save_fig("clade_agreement.pdf", p_clades, 7.1, 5)


cat("\n=== the", n_show, "largest clades not found in all", n_ref, "===\n")
clade_table %>%
  filter(n_trees < n_ref, size >= clade_min) %>%
  arrange(desc(size)) %>%
  mutate(tips = str_trunc(clade, 60)) %>%
  select(size, all_of(names(trees)), tips) %>%
  print(n = n_show)



# topology test ----------------------------------------------------------------

### AU test (does the concatenated alignment reject a topology)
lines <- read_lines(f_test)

start <- which(str_detect(lines, "^Tree\\s+logL"))

## check if the report holds the USER TREES block
if (length(start) == 0) {
  warning("no USER TREES block found in ", f_test)
  au <- NULL
} else {

  block <- lines[(start[1] + 1):length(lines)]
  block <- block[str_detect(block, "^\\s*[0-9]+\\s")]

  au <- tibble(raw = block) %>%
    mutate(fields = str_split(str_trim(raw), "\\s+")) %>%
    mutate(tree   = tree_names[as.integer(map_chr(fields, 1))],
           logL   = as.numeric(map_chr(fields, 2)),
           deltaL = as.numeric(map_chr(fields, 3))) %>%
    select(tree, logL, deltaL)

  ### p values by pattern (column positions depend on the tests run)
  au <- au %>%
    mutate(p_AU = map_dbl(block, function(l) {
      m <- str_match_all(l, "([0-9.eE+-]+)\\s*[-+]")[[1]]
      suppressWarnings(as.numeric(tail(m[, 2], 1)))
    }))

  write_tsv(au, file.path(out_root, "topology_test.tsv"))

  cat("\n=== topology test ===\n")
  print(au, n = Inf)

  p_au <- ggplot(au, aes(x = reorder(tree, deltaL), y = deltaL, fill = tree)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = col_tree, guide = "none") +
    coord_flip() +
    labs(x = NULL, y = "Log likelihood difference from the best topology")

  save_fig("topology_test.pdf", p_au, 7.1, 4)
}



# tables -----------------------------------------------------------------------

comparison_summary <- bind_rows(
  tibble(quantity = "shared tips",             value = as.character(n_tips)),
  tibble(quantity = sprintf("clades in %d of %d", n_ref:1, n_ref),
         value    = as.character(map_int(n_ref:1, ~ sum(clade_table$n_trees == .x)))),
  tibble(quantity = "largest disputed clade",
         value = as.character(max(c(0, clade_table$size[clade_table$n_trees < n_ref])))),
  distances %>% transmute(quantity = paste0("normalised RF, ", tree_a, " vs ", tree_b),
                          value    = sprintf("%.4f", rf_normalised)),
  if (!is.null(au)) {
    au %>% transmute(quantity = paste0("AU p value, ", tree), value = sprintf("%.4f", p_AU))
  })

write_tsv(comparison_summary, file.path(out_root, "comparison_summary.tsv"))

cat("\n=== summary ===\n")
print(comparison_summary, n = Inf)



# figures ----------------------------------------------------------------------

save_fig("tree_comparison_combined.pdf",
         (p_dist + p_clades) +
           plot_annotation(tag_levels = "A"), 7.1, 3)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  topology_distances.pdf\n")
cat("  clade_agreement.pdf\n")
cat("  topology_test.pdf\n")
cat("  tree_comparison_combined.pdf\n")
cat("  topology_distances.tsv\n")
cat("  clade_presence.tsv\n")
cat("  clade_summary.tsv\n")
cat("  topology_test.tsv\n")
cat("  comparison_summary.tsv\n")
