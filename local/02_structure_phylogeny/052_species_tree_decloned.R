## ============================================================================
## 052_species_tree_decloned.R
## The concatenation tree on one isolate per MLG with its ultrafast bootstrap
## support, the clonemates drawn on their representative, and the missing data
## diagnostic
##
## Files needed:
##   species_tree_decloned.treefile  ML tree with UFBoot support on the MLG representatives,
##                                   from hpc/02_structure_phylogeny/052_species_tree_decloned.sh
##                                   (copy species_tree_decloned.treefile of
##                                    ${ROOT}/02_structure_phylogeny/results/052_results_species_tree_decloned
##                                    to local/02_structure_phylogeny/data/052_data_species_tree_decloned)
##   tip_missingness.tsv             N content per tip, from hpc/02_structure_phylogeny/041_panel_diagnostics.sh
##                                   (copy tip_missingness.tsv of
##                                    ${ROOT}/02_structure_phylogeny/results/041_results_panel_diagnostics
##                                    to local/02_structure_phylogeny/data/052_data_species_tree_decloned)
##   sample_metadata_ALL_FIXED.tsv   metadata and clone status, from local/01_variant_calling/103_mating_type.R
##   lineage_labels.tsv              MLG, representative and lineage per sample at K = 4 to 6,
##                                   from local/03_delimitation/001_lineage_labels.R
## ============================================================================

library(ape)
library(ggtree)
library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

data_root <- "02_structure_phylogeny/data/052_data_species_tree_decloned"
out_root  <- "02_structure_phylogeny/results/052_results_species_tree_decloned"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

f_tree    <- file.path(data_root, "species_tree_decloned.treefile")
f_n       <- file.path(data_root, "tip_missingness.tsv")
f_meta    <- "01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv"
f_cluster <- "03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv"



# parameters -------------------------------------------------------------------

### define outgroup
outgroup <- c("PsefloM402_1", "PsefloM405_1", "PsefloM406_1",
              "M4-01", "M4-02", "M4-03", "M4-04", "M4-05", "M4-06")
genome_pattern <- "^(Cenge|Pseflo)"
reference_tip  <- "consensus"

### expected tips (412 MLG representatives, 6 mapped outgroup isolates, 26 assemblies, reference; 618 with the clonemates)
n_tips_expected      <- 445
n_tips_full_expected <- 618

### UFBoot bands (>= 95 supported)
support_high <- 95
support_low  <- 70

### lineage levels (branches among and within lineages, and N between lineages, per level)
k_levels  <- c("k4", "k5", "k6")
lab_level <- c(k4 = "K = 4", k5 = "K = 5", k6 = "K = 6")

### tips listed with the most N, largest clades listed
n_high <- 20
n_show <- 20

### jitter of the points over the boxplots
seed <- 42



# style ------------------------------------------------------------------------

### annotations
support_label_cut <- 100
label_size        <- 0.9
label_hjust       <- 2.2
tip_size          <- 1.1
tree_width        <- 24
tree_height       <- 60
label_space       <- 0.06
scale_width       <- 0.01
clone_colour      <- "grey70"

### UFBoot colours
col_boot_band <- c("< 70"     = "#D73027",
                   "70 to 94" = "#FDAE61",
                   ">= 95"    = "#1A9850")

### branch classes
col_class <- c("among lineages" = "grey15", "within a lineage" = "grey80")

### N content, low to high
col_n <- c("#4575B4", "#FFFFBF", "#D73027")

### lineage and project colours, shared with 03
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

pal_project <- c(magic = "white", phytoakmeter = "grey25")
lab_project <- c(magic = "MAGIC", phytoakmeter = "PhytOakmeter")

### variables in fixed colours (larger points)
fill_vars <- c("cluster", "project")

pal_discrete <- viridisLite::viridis

### categorical tip colours
pal_cat_base <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                  "#A65628", "#F781BF", "#00868B", "#B8860B", "#666666",
                  "#7570B3", "#8C564B")

pal_cat <- function(n) {
  if (n <= length(pal_cat_base)) pal_cat_base[seq_len(n)] else pal_discrete(n)
}

### base style
theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf, limitsize = FALSE)
}

### p values as in the publication tables
fmt_p <- function(p) if_else(p < 0.001, "< 0.001", sprintf("%.3f", p))



# read -------------------------------------------------------------------------

tree_raw <- read.tree(f_tree)

### root on the outgroup (via an ingroup tip first, root() misplaces the root
### when the outgroup spans the stored root)
tree <- root(tree_raw, outgroup = setdiff(tree_raw$tip.label, outgroup)[1],
             resolve.root = TRUE, edgelabel = TRUE)
tree <- root(tree, outgroup = outgroup, resolve.root = TRUE, edgelabel = TRUE)
tree <- ladderize(tree)

meta <- read_tsv(f_meta, show_col_types = FALSE)

### lineage labels at every level (clonemates carry the label of their MLG)
labels <- read_tsv(f_cluster, show_col_types = FALSE)

n_raw <- read_tsv(f_n, show_col_types = FALSE) %>%
  select(tip = sample, columns, n_pct = n_percent)

### clonemates (clone_status from the metadata)
## check if the metadata carries the clone status
if (!"clone_status" %in% names(meta)) stop("no clone_status column in ", f_meta)

clone_tips <- meta$sample[meta$clone_status %in% "Clonemate"]

### representative of each clonemate (same MLG, clustered in the labels)
grafts <- labels %>%
  filter(!clustered) %>%
  select(clonemate = sample, mlg) %>%
  left_join(labels %>%
              filter(clustered) %>%
              select(representative = sample, mlg),
            by = "mlg")

### every tip of the figures (the decloned tree and the clonemates)
tip_data <- tibble(tip = c(tree$tip.label, grafts$clonemate)) %>%
  left_join(n_raw, by = "tip") %>%
  left_join(select(meta, tip = sample, project, country, site), by = "tip") %>%
  left_join(labels %>%
              filter(clustered) %>%
              select(tip = sample, cluster = k6), by = "tip") %>%
  mutate(is_genome    = str_detect(tip, genome_pattern),
         is_reference = tip == reference_tip,
         is_sample    = !is_genome & !is_reference,
         is_clone     = tip %in% clone_tips)

### one plotting column per tip variable
cat_vars <- c("cluster", "country", "site", "project")

for (v in cat_vars) {
  tip_data[[paste0(v, "_plot")]] <-
    if_else(tip_data$is_clone, "clonemate", as.character(tip_data[[v]]))
}

### column colours, clonemate level set to grey
cat_scale <- function(v, name = NULL, aesthetics = "colour") {
  lv <- setdiff(sort(unique(na.omit(tip_data[[paste0(v, "_plot")]]))), "clonemate")
  ### lineages and projects in their fixed colours (outgroup samples grey), the
  ### other variables from the palette
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

message("tips: ", Ntip(tree), "   with the clonemates: ", nrow(tip_data),
        "   clonemates: ", nrow(grafts),
        "   MLGs with clonemates: ", n_distinct(grafts$representative))



# checks -----------------------------------------------------------------------

## check if the support values were parsed
if (all(is.na(suppressWarnings(as.numeric(tree$node.label))))) {
  stop("no numeric support labels on the tree")
}

## check if the tree holds the decloned panel
if (Ntip(tree) != n_tips_expected) {
  warning("tree has ", Ntip(tree), " tips, expected ", n_tips_expected)
}

## check if the tree is fully resolved
if (!is.binary(tree)) warning("the ML tree has polytomies")

## check if the outgroup is a clade
if (!all(outgroup %in% tree$tip.label)) {
  warning("outgroup tips missing from the tree")
} else if (!is.monophyletic(tree, outgroup)) {
  warning("the outgroup is not monophyletic in the rooted tree")
}

### placement of the mapped outgroup samples
m4 <- grep("^M4-", tree$tip.label, value = TRUE)

if (length(m4) > 0) {
  message("mapped outgroup tips: ", length(m4),
          "   monophyletic: ", is.monophyletic(tree, m4),
          "   with the assemblies: ",
          is.monophyletic(tree, c(m4, intersect(outgroup, tree$tip.label))))
}

## check if every clonemate has exactly one representative in the tree
if (anyNA(grafts$representative) || anyDuplicated(grafts$clonemate) > 0) {
  stop("clonemates without exactly one representative of their MLG")
}

## check if every representative is a tip of the tree
if (!all(grafts$representative %in% tree$tip.label)) {
  stop("representatives missing from the tree: ",
       paste(head(setdiff(grafts$representative, tree$tip.label), 10), collapse = ", "))
}

## check if the tree holds no clonemate
if (any(grafts$clonemate %in% tree$tip.label)) {
  stop(sum(grafts$clonemate %in% tree$tip.label), " clonemates are tips of the decloned tree")
}

## check if the clonemates of the labels are those of the metadata
if (!setequal(grafts$clonemate, intersect(clone_tips, labels$sample))) {
  warning("clonemates differ between the labels and the metadata: ",
          length(setdiff(grafts$clonemate, clone_tips)), " only in the labels, ",
          length(setdiff(intersect(clone_tips, labels$sample), grafts$clonemate)),
          " only in the metadata")
}

## check if the N table covers every tip
if (anyNA(tip_data$n_pct)) {
  stop(sum(is.na(tip_data$n_pct)), " tips without an N value in ", f_n)
}

## check if the tips share one alignment length
if (n_distinct(tip_data$columns) != 1) stop("the tips do not share one alignment length")

## check if every ingroup sample has a lineage at every level
n_unlabelled <- sum(tip_data$is_sample & !tip_data$tip %in% outgroup &
                      !tip_data$tip %in% filter(labels, if_all(all_of(k_levels), ~ !is.na(.x)))$sample)
if (n_unlabelled > 0) {
  warning(n_unlabelled, " ingroup sample tips without a lineage at every level")
}

## check if every lineage label has a colour
if (any(!is.na(tip_data$cluster) & !tip_data$cluster %in% names(pal_unit))) {
  warning("lineage labels without a colour in pal_unit: ",
          paste(setdiff(na.omit(tip_data$cluster), names(pal_unit)), collapse = ", "))
}

## check if the metadata covers the sample tips
missing_meta <- sum(is.na(tip_data$project) & tip_data$is_sample)

if (missing_meta > 0) {
  warning(missing_meta, " sample tips without metadata")
}



# support ----------------------------------------------------------------------

### UFBoot band of a support value
boot_band_of <- function(x) {
  factor(case_when(is.na(x)          ~ NA_character_,
                   x < support_low   ~ "< 70",
                   x < support_high  ~ "70 to 94",
                   TRUE              ~ ">= 95"),
         levels = names(col_boot_band))
}

### UFBoot and clade size per internal node of the decloned tree (root excluded)
parts <- prop.part(tree)

### lineages below each node, per level (assemblies and the reference carry none)
n_units <- function(kk) {
  lab <- setNames(labels[[kk]], labels$sample)
  map_int(parts, ~ n_distinct(na.omit(lab[attr(parts, "labels")[.x]])))
}

annot <- tibble(node    = Ntip(tree) + seq_len(Nnode(tree)),
                support = suppressWarnings(as.numeric(tree$node.label)),
                n_tips  = lengths(parts),
                in_og   = map_lgl(parts, ~ all(attr(parts, "labels")[.x] %in% outgroup)),
                !!!set_names(map(k_levels, n_units), paste0("units_", k_levels))) %>%
  filter(node != Ntip(tree) + 1) %>%
  mutate(boot_band = boot_band_of(support))

### ingroup branches (the outgroup and the root split excluded)
ingroup_support <- annot %>%
  filter(!in_og, n_tips < Ntip(tree) - length(outgroup), !is.na(support))

### branch class per level (among lineages: the clade holds samples of more than one lineage)
support_classes <- map_dfr(k_levels, function(kk) {
  ingroup_support %>%
    mutate(level = kk,
           class = if_else(.data[[paste0("units_", kk)]] > 1, "among lineages", "within a lineage"))
})

### support summary for all ingroup branches and per level among and within lineages
summ <- function(d) {
  d %>%
    summarise(branches = n(),
              median   = median(support),
              q25      = quantile(support, 0.25, names = FALSE),
              q75      = quantile(support, 0.75, names = FALSE),
              high     = sum(support >= support_high),
              below_70 = sum(support < support_low),
              lowest   = min(support),
              .by = c(level, class))
}

support_summary <- bind_rows(
  summ(mutate(ingroup_support, level = "all", class = "all ingroup branches")),
  summ(support_classes) %>% arrange(level, class))



# clonemates -------------------------------------------------------------------

### clonemates bound at zero length on their representative (figures only,
### the UFBoot values stay those of the decloned tree)
graft <- function(tr, rep_tip, clones) {
  star <- stree(length(clones) + 1, tip.label = c(rep_tip, clones))
  star$edge.length <- rep(0, Nedge(star))
  star$node.label  <- ""
  bind.tree(tr, star, where = which(tr$tip.label == rep_tip))
}

groups <- split(grafts$clonemate, grafts$representative)

tree_full <- reduce2(names(groups), groups, graft, .init = tree) %>%
  ladderize()

## check if the grafted tree holds every tip once
if (Ntip(tree_full) != n_tips_full_expected || anyDuplicated(tree_full$tip.label) > 0) {
  warning("grafted tree has ", Ntip(tree_full), " tips, expected ", n_tips_full_expected)
}

## check if grafting left the decloned topology unchanged
if (!all.equal.phylo(drop.tip(tree_full, grafts$clonemate), tree, use.edge.length = FALSE)) {
  stop("the grafted tree pruned of the clonemates differs from the decloned tree")
}



# missingness ------------------------------------------------------------------

### ingroup samples of the decloned tree (the outgroup sits across the root and
### its N reflects divergence from the reference)
sample_tips <- tip_data %>%
  filter(tip %in% tree$tip.label, is_sample, !tip %in% outgroup)

### phylogenetic autocorrelation of N (inverse patristic distance as weight)
patristic <- cophenetic.phylo(tree)[sample_tips$tip, sample_tips$tip]

### pairs at zero distance carry no weight (infinite inverse distance)
w <- 1 / patristic
w[!is.finite(w)] <- 0
diag(w) <- 0

moran <- Moran.I(sample_tips$n_pct, w)

### root to tip distance against N
root_dist <- setNames(node.depth.edgelength(tree)[seq_len(Ntip(tree))], tree$tip.label)

sample_tips <- mutate(sample_tips, root_to_tip = root_dist[tip])

rt_cor <- cor.test(sample_tips$n_pct, sample_tips$root_to_tip,
                   method = "spearman", exact = FALSE)

### N between lineages, per level (MLG representatives)
cluster_tips <- labels %>%
  filter(clustered) %>%
  select(tip = sample, all_of(k_levels)) %>%
  pivot_longer(all_of(k_levels), names_to = "level", values_to = "cluster") %>%
  filter(!is.na(cluster)) %>%
  inner_join(select(sample_tips, tip, n_pct), by = "tip")

kw <- cluster_tips %>%
  nest(d = -level) %>%
  mutate(test  = map(d, ~ kruskal.test(.x$n_pct, factor(.x$cluster))),
         n     = map_int(d, nrow),
         chisq = map_dbl(test, ~ unname(.x$statistic)),
         df    = map_dbl(test, ~ unname(.x$parameter)),
         p     = map_dbl(test, ~ .x$p.value)) %>%
  select(level, n, chisq, df, p) %>%
  arrange(level)

cluster_n <- cluster_tips %>%
  summarise(n_tips   = n(),
            median_n = median(n_pct),
            max_n    = max(n_pct),
            .by = c(level, cluster)) %>%
  arrange(level, desc(median_n))

diagnostics <- bind_rows(
  tibble(test      = "Moran's I, N along the tree",
         level     = "",
         n         = nrow(sample_tips),
         statistic = moran$observed,
         expected  = moran$expected,
         df        = NA_real_,
         p_value   = moran$p.value),
  tibble(test      = "Spearman, N against root to tip distance",
         level     = "",
         n         = nrow(sample_tips),
         statistic = unname(rt_cor$estimate),
         expected  = 0,
         df        = NA_real_,
         p_value   = rt_cor$p.value),
  kw %>% transmute(test      = "Kruskal-Wallis, N between lineages",
                   level     = lab_level[level],
                   n,
                   statistic = chisq,
                   expected  = NA_real_,
                   df,
                   p_value   = p))



# figures ----------------------------------------------------------------------

### UFBoot of the decloned tree on its nodes in the grafted tree (clonemate nodes carry none)
annot_full <- tibble(node    = Ntip(tree_full) + seq_len(Nnode(tree_full)),
                     support = suppressWarnings(as.numeric(tree_full$node.label))) %>%
  mutate(boot_band = boot_band_of(support))

plot_data <- bind_rows(
  tibble(node = seq_len(Ntip(tree_full)), label = tree_full$tip.label),
  tibble(node = annot_full$node, label = NA_character_)) %>%
  left_join(annot_full, by = "node") %>%
  left_join(tip_data, by = c("label" = "tip"))

### weakest band on the path from the root is carried downstream
band_levels <- c(">= 95", "70 to 94", "< 70")

inherit_band <- function(tr, band_num) {
  out <- band_num
  ### the root has no branch of its own
  out[Ntip(tr) + 1] <- 1
  ### cladewise order visits a parent before its children
  ed <- reorder(tr, "cladewise")$edge
  for (i in seq_len(nrow(ed))) {
    out[ed[i, 2]] <- max(out[ed[i, 1]], out[ed[i, 2]], na.rm = TRUE)
  }
  out
}

band_num <- rep(NA_integer_, Ntip(tree_full) + Nnode(tree_full))
band_num[annot_full$node] <- match(as.character(annot_full$boot_band), band_levels)

plot_data <- plot_data %>%
  mutate(branch_band = factor(band_levels[inherit_band(tree_full, band_num)[node]],
                              levels = names(col_boot_band)))

### the tip variable carries the figure, the support has its own tree
base_tree <- function() {
  ggtree(tree_full, linewidth = 0.2) %<+% plot_data +
    ### tip names aligned as a column
    geom_tiplab(size = tip_size, align = TRUE, linesize = 0.1,
                linetype = "dotted") +
    geom_treescale(width = scale_width, fontsize = 3, linesize = 0.3, offset = 2) +
    hexpand(label_space)
}

### one figure per tip variable
tip_figures <- list(
  cluster = list(var = "cluster", name = "Lineage (K = 6)"),
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

  p <- if (v %in% fill_vars) {
    base_tree() +
      geom_tippoint(aes(fill = .data[[vp]]), shape = 21, colour = "black",
                    stroke = 0.1, size = 1.1, na.rm = TRUE) +
      cat_scale(v, name = tip_figures[[nm]]$name, aesthetics = "fill") +
      guides(fill = guide_legend(override.aes = list(size = 3)))
  } else {
    base_tree() +
      geom_tippoint(aes(fill = .data[[vp]]), shape = 21, colour = "black",
                    stroke = 0.2, size = 0.9, na.rm = TRUE) +
      cat_scale(v, name = tip_figures[[nm]]$name, aesthetics = "fill")
  }

  p <- p + theme(legend.text = element_text(size = 7))

  save_fig(sprintf("species_tree_%s.pdf", nm), p, tree_width, tree_height)
}


### support tree, value printed and branch coloured by its inherited band
p_tree_support <- ggtree(tree_full, linewidth = 0.3) %<+% plot_data +
  aes(colour = branch_band) +
  ### support values
  geom_text2(aes(subset = !isTip & !is.na(support) & support <= support_label_cut,
                 label = sprintf("%.0f", support)),
             size = label_size, hjust = label_hjust, vjust = -0.4, colour = "black") +
  scale_colour_manual(values = col_boot_band, name = "UFBoot",
                      drop = FALSE, na.value = "grey70") +
  geom_tiplab(size = tip_size, align = TRUE, linesize = 0.1,
              linetype = "dotted", colour = "black") +
  geom_treescale(width = scale_width, fontsize = 3, linesize = 0.3, offset = 2) +
  hexpand(label_space)

save_fig("species_tree_support.pdf", p_tree_support, tree_width, tree_height)


### N content per tip (all tips coloured, the value is the point of the figure)
p_tree_n <- base_tree() +
  geom_tippoint(aes(fill = n_pct), shape = 21, colour = "black", stroke = 0.1,
                size = 1.1) +
  scale_fill_gradientn(colours = col_n, name = "N (%)")

save_fig("species_tree_missingness.pdf", p_tree_n, tree_width, tree_height)


### rooted fan (branch lengths dropped, the outgroup branch would take the whole radius)
p_fan <- ggtree(tree_full, layout = "fan", open.angle = 10, branch.length = "none",
                linewidth = 0.2) %<+% plot_data +
  geom_tippoint(aes(fill = cluster_plot), shape = 21, colour = "black",
                stroke = 0.1, size = 1.4, na.rm = TRUE) +
  ### labels follow the radius, the canvas carries 618 of them
  geom_tiplab(size = 1.6, offset = 0.6, colour = "black") +
  hexpand(0.2) +
  cat_scale("cluster", name = "Lineage (K = 6)", aesthetics = "fill") +
  guides(fill = guide_legend(override.aes = list(size = 3)))

save_fig("species_tree_fan.pdf", p_fan, 30, 30)


### unrooted equal angle layout of the ingroup (the lineages relative to each
### other, the outgroup branch would take the whole canvas)
tree_net <- drop.tip(tree_full, outgroup)

p_net <- ggtree(tree_net, layout = "equal_angle", linewidth = 0.2) %<+%
  select(tip_data, label = tip, cluster_plot) +
  geom_tippoint(aes(fill = cluster_plot), shape = 21, colour = "black",
                stroke = 0.1, size = 1.3, na.rm = TRUE) +
  cat_scale("cluster", name = "Lineage (K = 6)", aesthetics = "fill") +
  guides(fill = guide_legend(override.aes = list(size = 3)))

save_fig("species_tree_network.pdf", p_net, 7.1, 7.1)


### support distribution, and support against clade size among and within lineages
p_hist <- ggplot(ingroup_support, aes(x = support)) +
  geom_histogram(binwidth = 2, fill = "grey40") +
  geom_vline(xintercept = support_high, colour = "#D73027", linetype = "dashed",
             linewidth = 0.5) +
  labs(x = "Ultrafast bootstrap support (%)", y = "N internal branches")

### among-lineage branches drawn last (on top)
p_size <- support_classes %>%
  filter(level == "k6") %>%
  arrange(class == "among lineages") %>%
  ggplot(aes(x = n_tips, y = support, fill = class)) +
  geom_hline(yintercept = support_high, colour = "grey40", linetype = "dashed",
             linewidth = 0.4) +
  geom_point(shape = 21, colour = "black", stroke = 0.2, size = 1.4) +
  scale_fill_manual(values = col_class, name = lab_level[["k6"]]) +
  scale_x_log10() +
  labs(x = "N tips in the clade", y = "Ultrafast bootstrap support (%)")

save_fig("species_tree_summary.pdf",
         (p_hist + p_size) + plot_annotation(tag_levels = "A"), 7.1, 3.8)


### N against root to tip distance, and N between lineages per level
p_rt <- sample_tips %>%
  mutate(cluster = replace_na(cluster, "none")) %>%
  ggplot(aes(x = n_pct, y = root_to_tip, fill = cluster)) +
  geom_point(shape = 21, colour = "black", stroke = 0.15, size = 1.3) +
  annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.4, size = 3,
           label = sprintf("rho = %.2f, p %s", rt_cor$estimate,
                           if_else(rt_cor$p.value < 0.001, "< 0.001",
                                   sprintf("= %.3f", rt_cor$p.value)))) +
  scale_fill_manual(values = c(pal_unit, none = "grey85"), name = lab_level[["k6"]]) +
  labs(x = "N in the concatenation (%)", y = "Root to tip distance") +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 2.5)))

kw_lab <- kw %>%
  mutate(label = sprintf("p %s", if_else(p < 0.001, "< 0.001", sprintf("= %.3f", p))))

p_cluster_n <- ggplot(cluster_tips, aes(x = cluster, y = n_pct)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3, width = 0.6) +
  geom_point(aes(fill = cluster), shape = 21, colour = "black", stroke = 0.1,
             size = 0.9, position = position_jitter(width = 0.15, height = 0, seed = seed)) +
  geom_text(data = kw_lab, aes(x = -Inf, y = Inf, label = label), hjust = -0.1,
            vjust = 1.4, size = 3, inherit.aes = FALSE) +
  facet_wrap(~ level, scales = "free_x", labeller = as_labeller(lab_level)) +
  scale_fill_manual(values = pal_unit, guide = "none") +
  labs(x = "Lineage", y = "N in the concatenation (%)") +
  theme(axis.text.x = element_text(size = 7))

save_fig("species_tree_missingness_tests.pdf",
         (p_rt + p_cluster_n) +
           plot_layout(widths = c(1, 2.6), guides = "collect") +
           plot_annotation(tag_levels = "A") &
           theme(legend.position = "bottom"), 7.1, 3.8)



# tables -----------------------------------------------------------------------

publication_support <- support_summary %>%
  transmute(Level    = if_else(level == "all", "", lab_level[level]),
            Branches = class,
            N        = branches,
            `UFBoot, median [IQR]` = sprintf("%.1f [%.1f, %.1f]", median, q25, q75),
            !!sprintf("UFBoot >= %d", support_high) :=
              sprintf("%d (%.1f%%)", high, 100 * high / branches),
            Lowest   = sprintf("%.0f", lowest))

publication_missingness <- diagnostics %>%
  transmute(Test      = test,
            Level     = level,
            N         = n,
            Statistic = sprintf("%.3f", statistic),
            Expected  = if_else(is.na(expected), "", sprintf("%.3f", expected)),
            df        = if_else(is.na(df), "", sprintf("%.0f", df)),
            p         = fmt_p(p_value))

### rooted decloned tree with UFBoot, and the same with the clonemates grafted (for 100)
write.tree(tree, file.path(out_root, "species_tree_decloned_rooted.tre"))
write.tree(tree_full, file.path(out_root, "species_tree_grafted_rooted.tre"))

write_tsv(arrange(grafts, representative, clonemate), file.path(out_root, "clonemate_grafts.tsv"))
write_tsv(select(annot, node, n_tips, in_og, starts_with("units_"), support, boot_band),
          file.path(out_root, "species_tree_support.tsv"))
write_tsv(support_summary, file.path(out_root, "species_tree_support_summary.tsv"))
write_tsv(publication_support, file.path(out_root, "publication_species_tree_support.tsv"))
write_tsv(select(tip_data, -ends_with("_plot")) %>% arrange(desc(n_pct)),
          file.path(out_root, "tip_missingness.tsv"))
write_tsv(cluster_n, file.path(out_root, "missingness_by_cluster.tsv"))
write_tsv(diagnostics, file.path(out_root, "missingness_diagnostics.tsv"))
write_tsv(publication_missingness, file.path(out_root, "publication_missingness_diagnostics.tsv"))

cat("\n=== UFBoot, ingroup ===\n")
print(publication_support, n = Inf, width = Inf)

cat("\n=== missing data diagnostics ===\n")
print(publication_missingness, n = Inf, width = Inf)

cat("\n=== the", n_show, "largest ingroup clades ===\n")
ingroup_support %>%
  arrange(desc(n_tips)) %>%
  select(node, n_tips, support) %>%
  print(n = n_show)

cat("\n=== the", n_high, "ingroup samples with the most N ===\n")
sample_tips %>%
  arrange(desc(n_pct)) %>%
  select(tip, project, country, cluster, n_pct, root_to_tip) %>%
  print(n = n_high)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  species_tree_cluster.pdf                     tips coloured by lineage, K = 6\n")
cat("  species_tree_country.pdf                     tips coloured by country\n")
cat("  species_tree_site.pdf                        tips coloured by sampling site\n")
cat("  species_tree_project.pdf                     tips coloured by project\n")
cat("  species_tree_support.pdf                     UFBoot printed and coloured per branch\n")
cat("  species_tree_missingness.pdf                 tips coloured by N content\n")
cat("  species_tree_fan.pdf                         rooted fan, tips coloured by lineage\n")
cat("  species_tree_network.pdf                     unrooted layout, tips coloured by lineage\n")
cat("  species_tree_summary.pdf                     support distribution, and support against clade size among and within lineages\n")
cat("  species_tree_missingness_tests.pdf           N against root to tip distance, and N between lineages per level\n")
cat("  species_tree_decloned_rooted.tre             rooted on the outgroup, UFBoot, MLG representatives\n")
cat("  species_tree_grafted_rooted.tre              the same with the clonemates at zero length on their representative, for 100\n")
cat("  clonemate_grafts.tsv                         clonemate, MLG and representative\n")
cat("  species_tree_support.tsv                     UFBoot, clade size and lineages below per branch\n")
cat("  species_tree_support_summary.tsv             support for all ingroup branches and per level among and within lineages\n")
cat("  publication_species_tree_support.tsv         the same, formatted\n")
cat("  tip_missingness.tsv                          N content and metadata per tip\n")
cat("  missingness_by_cluster.tsv                   N per lineage and level\n")
cat("  missingness_diagnostics.tsv                  Moran's I, Spearman and Kruskal-Wallis\n")
cat("  publication_missingness_diagnostics.tsv      the same, formatted\n")
