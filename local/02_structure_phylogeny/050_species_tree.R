## ============================================================================
## 050_species_tree.R
## Species tree figures and the missing-data diagnostic (does the N content of a
## tip predict where that tip sits in the topology)
##
## Files needed:
##   species_tree.treefile          ML tree with UFBoot support,
##                                  from hpc/02_structure_phylogeny/050_species_tree.sh
##                                  (copy ${ROOT}/02_structure_phylogeny/results/050_results_species_tree/species_tree.treefile
##                                   to local/02_structure_phylogeny/data/050_data_species_tree)
##   tip_missingness.tsv            N content per tip,
##                                  from hpc/02_structure_phylogeny/041_panel_diagnostics.sh
##                                  (copy ${ROOT}/02_structure_phylogeny/results/041_results_panel_diagnostics/tip_missingness.tsv
##                                   to local/02_structure_phylogeny/data/050_data_species_tree)
##   sample_metadata_ALL_FIXED.tsv  metadata and clone status,
##                                  from local/01_variant_calling/103_mating_type.R
##   lineage_labels.tsv             lineage per sample at K = 4 to 6,
##                                  from local/03_delimitation/001_lineage_labels.R
## ============================================================================

library(ape)
library(ggtree)
library(ggnewscale)
library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

data_root <- "02_structure_phylogeny/data/050_data_species_tree"
out_root  <- "02_structure_phylogeny/results/050_results_species_tree"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

f_tree   <- file.path(data_root, "species_tree.treefile")
f_n      <- file.path(data_root, "tip_missingness.tsv")
f_meta   <- "01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv"

f_lineage <- "03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv"



# parameters -------------------------------------------------------------------

### levels carried forward from 014 (all tested in the diagnostic)
k_focus   <- c(4, 5, 6)

### level the tree is coloured by
k_lineage <- 6

### define outgroup
outgroup <- c("PsefloM402_1", "PsefloM405_1", "PsefloM406_1",
              "M4-01", "M4-02", "M4-03", "M4-04", "M4-05", "M4-06")

### assemblies (tips without a mapping history, so without missing data)
genome_pattern <- "^(Cenge|Pseflo)"

### reference sequence carried along as a tip
reference_tip <- "consensus"

### tips listed in the missingness table, seed
n_high <- 20
seed   <- 42

### support labels printed up to this value (100: every branch)
boot_label_cut <- 100
label_size     <- 0.9
label_hjust    <- 2.2

### tree figures (every tip drawn)
tip_size    <- 1.1
tree_width  <- 24
tree_height <- 60
label_space <- 0.06

### clonemates in one grey (no independent information)
clone_colour <- "grey70"



# style ------------------------------------------------------------------------

### low to high N
col_n <- c("#4575B4", "#FFFFBF", "#D73027")

### bootstrap support bands
col_boot_band <- c("< 70"     = "#D73027",
                   "70 to 94" = "#FDAE61",
                   ">= 95"    = "#1A9850")

### lineage and project colours
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

pal_project <- c(magic = "white", phytoakmeter = "grey25")
lab_project <- c(magic = "MAGIC", phytoakmeter = "PhytOakmeter")

### variables in fixed colours (outlined points, black labels)
fill_vars <- c("cluster", "project")

pal_discrete <- viridisLite::viridis

### qualitative palette for the categorical tip colours (viridis beyond twelve
### levels)
pal_cat_base <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                  "#A65628", "#F781BF", "#00868B", "#B8860B", "#666666",
                  "#7570B3", "#8C564B")

pal_cat <- function(n) {
  if (n <= length(pal_cat_base)) pal_cat_base[seq_len(n)] else pal_discrete(n)
}

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf, limitsize = FALSE)
}



# read -------------------------------------------------------------------------

tree_raw <- read.tree(f_tree)

### rooted on the outgroup, first on an ingroup tip (root() misplaces the root
### when the outgroup spans the stored root; edgelabel keeps the support on its
### branch)
tree <- root(tree_raw, outgroup = setdiff(tree_raw$tip.label, outgroup)[1],
             resolve.root = TRUE, edgelabel = TRUE)
tree <- root(tree, outgroup = outgroup, resolve.root = TRUE,
             edgelabel = TRUE)
tree <- ladderize(tree)

meta <- read_tsv(f_meta, show_col_types = FALSE)

### clustered samples at every focal K (as in 014)
clusters <- read_tsv(f_lineage, show_col_types = FALSE) %>%
  filter(clustered) %>%
  pivot_longer(c(k4, k5, k6), names_to = "K", values_to = "cluster") %>%
  transmute(tip = sample, K = as.integer(str_remove(K, "^k")), cluster) %>%
  filter(K %in% k_focus)

n_raw <- read_tsv(f_n, show_col_types = FALSE) %>%
  select(tip = sample, columns, n_frac = n_percent)

### clonemates (clone_status from the metadata)
## check if the metadata carries the clone status
if (!"clone_status" %in% names(meta)) stop("no clone_status column in ", f_meta)

clone_tips <- meta$sample[meta$clone_status %in% "Clonemate"]

## check if the clonemates are in the tree
if (!any(clone_tips %in% tree$tip.label)) stop("no clonemate of the metadata is a tip of the tree")

tip_data <- tibble(tip = tree$tip.label) %>%
  left_join(n_raw, by = "tip") %>%
  left_join(select(meta, tip = sample, project, country, site), by = "tip") %>%
  mutate(is_genome    = str_detect(tip, genome_pattern),
         is_reference = tip == reference_tip,
         is_sample    = !is_genome & !is_reference,
         is_clone     = tip %in% clone_tips) %>%
  left_join(clusters %>%
              filter(K == k_lineage) %>%
              select(tip, cluster), by = "tip")

### plotting column per tip variable (clonemates as their own grey level)
cat_vars <- c("cluster", "country", "site", "project")

for (v in cat_vars) {
  tip_data[[paste0(v, "_plot")]] <-
    if_else(tip_data$is_clone, "clonemate", as.character(tip_data[[v]]))
}

### palette per plotting column (clonemate level grey and last)
cat_scale <- function(v, name = NULL, aesthetics = "colour") {
  lv <- setdiff(sort(unique(na.omit(tip_data[[paste0(v, "_plot")]]))), "clonemate")
  ### lineages and projects in their fixed colours (outgroup samples grey),
  ### other variables from pal_cat
  pal <- if (v == "cluster") {
    pal_unit[lv]
  } else if (v == "project") {
    pal_project[intersect(names(pal_project), lv)]
  } else {
    setNames(pal_cat(length(lv)), lv)
  }
  lab <- if (v == "project") c(lab_project[names(pal)], clonemate = "clonemate") else waiver()
  vals <- c(pal, clonemate = clone_colour)
  scale_colour_manual(values = vals, name = name, na.value = "grey85",
                      breaks = c(names(pal), "clonemate"), labels = lab,
                      aesthetics = aesthetics)
}



# checks -----------------------------------------------------------------------

## check if the N table covers every tip
if (anyNA(tip_data$n_frac)) {
  stop(sum(is.na(tip_data$n_frac)), " tips without an N value, check that ",
       "tip_missingness.tsv was written from the alignments of this tree")
}

n_columns <- unique(tip_data$columns)

## check if the tips share one alignment length
if (length(n_columns) != 1) {
  stop("the tips do not share one alignment length")
}

## check if the metadata covers the sample tips
missing_meta <- sum(is.na(tip_data$project) & tip_data$is_sample)

if (missing_meta > 0) {
  warning(missing_meta, " sample tips without metadata")
}

## check if the outgroup is a clade (the root was placed on it)
if (!is.monophyletic(tree, outgroup)) {
  warning("the outgroup is not monophyletic in the rooted tree")
}

### mapped outgroup samples (reference bias)
m4 <- grep("^M4-", tree$tip.label, value = TRUE)

if (length(m4) > 0) {
  message("mapped outgroup tips: ", length(m4),
          "   monophyletic: ", is.monophyletic(tree, m4),
          "   with the assemblies: ", is.monophyletic(tree, c(m4, outgroup)))
}

## check if every lineage label has a colour
if (any(!is.na(tip_data$cluster) & !tip_data$cluster %in% names(pal_unit))) {
  warning("lineage labels without a colour in pal_unit: ",
          paste(setdiff(na.omit(tip_data$cluster), names(pal_unit)), collapse = ", "))
}

## check if every focal K was read
if (!setequal(unique(clusters$K), k_focus)) {
  warning("cluster files missing for K = ",
          paste(setdiff(k_focus, unique(clusters$K)), collapse = ", "))
}

message("tips: ", Ntip(tree), "   samples: ", sum(tip_data$is_sample),
        "   genomes: ", sum(tip_data$is_genome),
        "   clonemates: ", sum(tip_data$is_clone))



# support ----------------------------------------------------------------------

### UFBoot per internal node, in three bands
annot <- tibble(node   = Ntip(tree) + seq_len(Nnode(tree)),
                ufboot = suppressWarnings(as.numeric(tree$node.label))) %>%
  mutate(boot_band = case_when(is.na(ufboot) ~ NA_character_,
                               ufboot < 70   ~ "< 70",
                               ufboot < 95   ~ "70 to 94",
                               TRUE          ~ ">= 95"),
         boot_band = factor(boot_band, levels = names(col_boot_band)))

support_summary <- annot %>%
  filter(!is.na(ufboot)) %>%
  summarise(branches = n(),
            median   = median(ufboot),
            below_95 = sum(ufboot < 95),
            below_70 = sum(ufboot < 70))



# tree -------------------------------------------------------------------------

plot_data <- bind_rows(
  tibble(node = seq_len(Ntip(tree)), label = tree$tip.label),
  tibble(node = annot$node, label = NA_character_)) %>%
  left_join(select(annot, node, ufboot, boot_band), by = "node") %>%
  left_join(tip_data, by = c("label" = "tip"))

### band per branch (weakest band on the path from the root, carried downstream)
band_levels <- c(">= 95", "70 to 94", "< 70")

inherit_band <- function(tr, band_num) {
  out <- band_num
  ### root starts in the best band
  out[Ntip(tr) + 1] <- 1
  ### cladewise order (parents before children)
  ed <- reorder(tr, "cladewise")$edge
  for (i in seq_len(nrow(ed))) {
    out[ed[i, 2]] <- max(out[ed[i, 1]], out[ed[i, 2]], na.rm = TRUE)
  }
  out
}

band_num <- rep(NA_integer_, Ntip(tree) + Nnode(tree))
band_num[annot$node] <- match(as.character(annot$boot_band), band_levels)

plot_data <- plot_data %>%
  mutate(branch_band = factor(band_levels[inherit_band(tree, band_num)[node]],
                              levels = names(col_boot_band)))

base_tree <- function() {
  ggtree(tree, linewidth = 0.2) %<+% plot_data +
    geom_nodepoint(aes(subset = !is.na(boot_band) & boot_band != ">= 95",
                       fill = boot_band),
                   shape = 21, size = 1.4, stroke = 0.2, colour = "black") +
    scale_fill_manual(values = col_boot_band, name = "UFBoot",
                      na.translate = FALSE, drop = FALSE) +
    geom_treescale(width = 0.01, fontsize = 3, linesize = 0.3, offset = 2) +
    hexpand(label_space)
}


### one figure per tip variable, all on the same tree
tip_figures <- list(
  cluster = list(var = "cluster", name = sprintf("Lineage (K = %d)", k_lineage)),
  country = list(var = "country", name = "Country"),
  site    = list(var = "site",    name = "Site"),
  project = list(var = "project", name = "Project")
)

for (nm in names(tip_figures)) {

  v <- tip_figures[[nm]]$var
  if (all(is.na(tip_data[[v]]))) {
    message("no values for ", v, ", that figure is skipped")
    next
  }

  vp <- paste0(v, "_plot")

  ### fixed colours as outlined points (a second fill scale next to the UFBoot one)
  p <- if (v %in% fill_vars) {
    base_tree() +
      new_scale_fill() +
      geom_tippoint(aes(fill = .data[[vp]]), shape = 21, colour = "black",
                    stroke = 0.1, size = 1.1, na.rm = TRUE) +
      geom_tiplab(size = tip_size, align = TRUE, linesize = 0.1,
                  linetype = "dotted", colour = "black") +
      cat_scale(v, name = tip_figures[[nm]]$name, aesthetics = "fill") +
      guides(fill = guide_legend(override.aes = list(size = 3)))
  } else {
    base_tree() +
      new_scale_fill() +
      geom_tippoint(aes(fill = .data[[vp]]), shape = 21, colour = "black",
                    stroke = 0.2, size = 0.9, na.rm = TRUE) +
      geom_tiplab(size = tip_size, align = TRUE, linesize = 0.1,
                  linetype = "dotted", colour = "black") +
      cat_scale(v, name = tip_figures[[nm]]$name, aesthetics = "fill")
  }

  p <- p +
    theme(legend.text = element_text(size = 7))

  save_fig(sprintf("species_tree_%s.pdf", nm), p, tree_width, tree_height)
}


### support tree (value on every branch, branch in its inherited band)
p_tree_support <- ggtree(tree, linewidth = 0.3) %<+% plot_data +
  aes(colour = branch_band) +
  geom_text2(aes(subset = !isTip & !is.na(ufboot) & ufboot <= boot_label_cut,
                 label = sprintf("%.0f", ufboot)),
             size = label_size, hjust = label_hjust, vjust = -0.4,
             colour = "black") +
  scale_colour_manual(values = col_boot_band, name = "UFBoot",
                      drop = FALSE, na.value = "grey70") +
  geom_tiplab(size = tip_size, align = TRUE, linesize = 0.1,
              linetype = "dotted", colour = "black") +
  hexpand(label_space)

save_fig("species_tree_support.pdf", p_tree_support, tree_width, tree_height)


### missing data on the tree (continuous scale, clonemates not greyed)
p_tree_n <- base_tree() +
  new_scale_fill() +
  geom_tippoint(aes(fill = n_frac), shape = 21, colour = "black",
                stroke = 0.2, size = 0.9) +
  geom_tiplab(size = tip_size, align = TRUE, linesize = 0.1,
              linetype = "dotted", colour = "black") +
  scale_fill_gradientn(colours = col_n, name = "N (%)")

save_fig("species_tree_missingness.pdf", p_tree_n, tree_width, tree_height)


### unrooted equal-angle layout
tree_net <- tree

if (is.null(tree_net$edge.length)) {
  tree_net$edge.length <- rep(1, nrow(tree_net$edge))
} else {
  tree_net$edge.length[is.na(tree_net$edge.length)] <- 0
}

p_net <- ggtree(tree_net, layout = "equal_angle", linewidth = 0.2) %<+% plot_data +
  geom_tippoint(aes(fill = cluster_plot), shape = 21, colour = "black",
                stroke = 0.1, size = 1.3, na.rm = TRUE) +
  cat_scale("cluster", name = sprintf("Lineage (K = %d)", k_lineage),
            aesthetics = "fill") +
  guides(fill = guide_legend(override.aes = list(size = 3)))

save_fig("species_tree_network.pdf", p_net, 7.1, 7.1)



# missingness ------------------------------------------------------------------

### missing-data diagnostic: N along the tree, N against branch length, N
### between clusters
sample_tips <- filter(tip_data, is_sample)

### Moran's I of N along the tree (patristic weights)
patristic <- cophenetic.phylo(tree)
patristic <- patristic[sample_tips$tip, sample_tips$tip]

### zero weight for pairs at distance zero (clonemates)
w <- 1 / patristic
w[!is.finite(w)] <- 0
diag(w) <- 0

set.seed(seed)
moran <- Moran.I(sample_tips$n_frac, w)

### root-to-tip distance against N (long branches of mostly unknown tips)
root_dist <- node.depth.edgelength(tree)[seq_len(Ntip(tree))]
names(root_dist) <- tree$tip.label

sample_tips <- mutate(sample_tips, root_to_tip = root_dist[tip])

rt_cor <- cor.test(sample_tips$n_frac, sample_tips$root_to_tip,
                   method = "spearman", exact = FALSE)

p_rt <- ggplot(sample_tips, aes(x = n_frac, y = root_to_tip)) +
  geom_point(shape = 21, colour = "black", fill = "grey55", size = 1,
             alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, colour = "#D55E00", linewidth = 0.6) +
  labs(x = "N in the concatenation (%)", y = "Root to tip distance")

### N between clusters at every focal K
cluster_tips <- sample_tips %>%
  select(tip, n_frac) %>%
  inner_join(clusters, by = "tip")

kw <- cluster_tips %>%
  group_by(K) %>%
  summarise(chisq = kruskal.test(n_frac ~ cluster)$statistic,
            df    = kruskal.test(n_frac ~ cluster)$parameter,
            p     = kruskal.test(n_frac ~ cluster)$p.value,
            .groups = "drop")

cluster_n <- cluster_tips %>%
  group_by(K, cluster) %>%
  summarise(n_tips = n(), median_n = median(n_frac), max_n = max(n_frac),
            .groups = "drop") %>%
  arrange(K, desc(median_n))

p_cluster_n <- ggplot(cluster_tips, aes(x = reorder(cluster, n_frac, median),
                                        y = n_frac)) +
  geom_boxplot(outlier.size = 0.6, linewidth = 0.4) +
  facet_wrap(~ K, scales = "free_x", labeller = label_both) +
  labs(x = "Cluster", y = "N in the concatenation (%)")

p_support <- ggplot(filter(annot, !is.na(ufboot)), aes(x = ufboot)) +
  geom_histogram(binwidth = 2, fill = "grey40") +
  labs(x = "Ultrafast bootstrap support", y = "Internal branches")

save_fig("species_tree_summary.pdf",
         (p_support + p_rt) / p_cluster_n +
           plot_annotation(tag_levels = "A"), 7.1, 5.9)



# tables -----------------------------------------------------------------------

write.tree(tree, file.path(out_root, "species_tree_rooted.tre"))

write_tsv(arrange(tip_data, desc(n_frac)), file.path(out_root, "tip_missingness.tsv"))
write_tsv(cluster_n, file.path(out_root, "missingness_by_cluster.tsv"))
write_tsv(support_summary, file.path(out_root, "support_summary.tsv"))

diagnostics <- bind_rows(
  tibble(test      = "Moran's I of N along the tree",
         statistic = moran$observed,
         p_value   = moran$p.value),
  tibble(test      = "Spearman, N vs root to tip distance",
         statistic = as.numeric(rt_cor$estimate),
         p_value   = rt_cor$p.value),
  kw %>% transmute(test      = sprintf("Kruskal-Wallis, N between clusters at K = %d", K),
                   statistic = as.numeric(chisq),
                   p_value   = p))

write_tsv(diagnostics, file.path(out_root, "missingness_diagnostics.tsv"))

cat("\n=== branch support ===\n")
print(support_summary, n = Inf)

cat("\n=== missing data diagnostics ===\n")
print(diagnostics, n = Inf)

cat("\n=== the", n_high, "tips with the most N ===\n")
tip_data %>%
  filter(is_sample) %>%
  arrange(desc(n_frac)) %>%
  select(tip, project, country, cluster, n_frac) %>%
  print(n = n_high)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  species_tree_cluster.pdf      tips coloured by lineage, K =", k_lineage, "\n")
cat("  species_tree_country.pdf      tips coloured by country\n")
cat("  species_tree_site.pdf         tips coloured by sampling site\n")
cat("  species_tree_project.pdf      tips coloured by project\n")
cat("  species_tree_support.pdf      UFBoot printed and coloured per branch\n")
cat("  species_tree_missingness.pdf  tips coloured by N content\n")
cat("  species_tree_network.pdf      unrooted layout, tips coloured by lineage\n")
cat("  species_tree_summary.pdf      support, branch length and cluster diagnostics\n")
cat("  species_tree_rooted.tre       rooted on the outgroup\n")
cat("  tip_missingness.tsv\n")
cat("  missingness_by_cluster.tsv\n")
cat("  support_summary.tsv\n")
cat("  missingness_diagnostics.tsv\n")
