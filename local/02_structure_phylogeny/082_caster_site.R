## ============================================================================
## 082_caster_site.R
## Species tree from CASTER-site and its local bootstrap support
##
## Files needed:
##   caster_site_species_tree.tre   site-based species tree with local bootstrap,
##                                  from hpc/02_structure_phylogeny/082_caster_site.sh
##                                  (copy ${ROOT}/02_structure_phylogeny/results/082_results_caster_site/caster_site_species_tree.tre
##                                   to local/02_structure_phylogeny/data/082_data_caster_site)
##   sample_metadata_ALL_FIXED.tsv  metadata and clone status, from local/01_variant_calling/103_mating_type.R
##   lineage_labels.tsv             lineage per sample at K = 4 to 6, from local/03_delimitation/001_lineage_labels.R
## ============================================================================

library(ape)
library(ggtree)
library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

data_root <- "02_structure_phylogeny/data/082_data_caster_site"
out_root  <- "02_structure_phylogeny/results/082_results_caster_site"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

f_tree    <- file.path(data_root, "caster_site_species_tree.tre")
f_meta    <- "01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv"
f_cluster <- "03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv"



# parameters -------------------------------------------------------------------

### outgroup
outgroup <- c("PsefloM402_1", "PsefloM405_1", "PsefloM406_1",
              "M4-01", "M4-02", "M4-03", "M4-04", "M4-05", "M4-06")

### assemblies (tip name pattern) and the consensus reference tip
genome_pattern <- "^(Cenge|Pseflo)"
reference_tip  <- "consensus"

### expected tips (591 samples, 26 assemblies, reference)
n_tips_expected <- 618

### local bootstrap taken as resolved
support_high <- 95

### lineage levels (branches among and within lineages, per level)
k_levels  <- c("k4", "k5", "k6")
lab_level <- c(k4 = "K = 4", k5 = "K = 5", k6 = "K = 6")

### largest clades printed
n_show <- 20



# style ------------------------------------------------------------------------

### annotations
support_label_cut <- 100
label_size        <- 0.9
label_hjust       <- 2.2
tip_size          <- 1.1
tree_width        <- 24
tree_height       <- 60
label_space       <- 0.06
clone_colour      <- "grey70"

### local bootstrap colours (bands as in 050)
col_boot_band <- c("< 70"     = "#D73027",
                   "70 to 94" = "#FDAE61",
                   ">= 95"    = "#1A9850")

### branch classes
col_class <- c("among lineages" = "grey15", "within a lineage" = "grey80")

### lineage and project colours (shared with 03)
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

pal_project <- c(magic = "white", phytoakmeter = "grey25")
lab_project <- c(magic = "MAGIC", phytoakmeter = "PhytOakmeter")

### variables in fixed colours (larger points)
fill_vars <- c("cluster", "project")

pal_discrete <- viridisLite::viridis

### categorical tip colours (viridis above 12 levels)
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



# read -------------------------------------------------------------------------

tree_raw <- read.tree(f_tree)

### root on the outgroup (via an ingroup tip first, root() misplaces the root
### when the outgroup spans the stored root)
tree <- root(tree_raw, outgroup = setdiff(tree_raw$tip.label, outgroup)[1],
             resolve.root = TRUE, edgelabel = TRUE)
tree <- root(tree, outgroup = outgroup, resolve.root = TRUE, edgelabel = TRUE)
tree <- ladderize(tree)

### terminal branch lengths (not estimated by CASTER-site) drawn as 0
n_no_length <- sum(is.na(tree$edge.length))
tree$edge.length[is.na(tree$edge.length)] <- 0

meta <- read_tsv(f_meta, show_col_types = FALSE)

### lineage labels at every level (clonemates carry the label of their MLG)
labels <- read_tsv(f_cluster, show_col_types = FALSE)

### lineages for colouring (MLG representatives, clonemates grey)
clusters <- labels %>%
  filter(clustered) %>%
  transmute(tip = sample, cluster = k6)

tip_data <- tibble(tip = tree$tip.label) %>%
  left_join(select(meta, tip = sample, project, country, site), by = "tip") %>%
  mutate(is_genome    = str_detect(tip, genome_pattern),
         is_reference = tip == reference_tip,
         is_sample    = !is_genome & !is_reference) %>%
  left_join(clusters, by = "tip")

### clonemates (clone_status from the metadata)
## check if the metadata carries the clone status
if (!"clone_status" %in% names(meta)) stop("no clone_status column in ", f_meta)

clone_tips <- meta$sample[meta$clone_status %in% "Clonemate"]

## check if the clonemates are in the tree
if (!any(clone_tips %in% tree$tip.label)) stop("no clonemate of the metadata is a tip of the tree")

tip_data <- mutate(tip_data, is_clone = tip %in% clone_tips)

### one plotting column per tip variable
cat_vars <- c("cluster", "country", "site", "project")

for (v in cat_vars) {
  tip_data[[paste0(v, "_plot")]] <-
    if_else(tip_data$is_clone, "clonemate", as.character(tip_data[[v]]))
}

### column colours, clonemate level set to grey
cat_scale <- function(v, name = NULL, aesthetics = "colour") {
  lv <- setdiff(sort(unique(na.omit(tip_data[[paste0(v, "_plot")]]))), "clonemate")
  ### fixed colours for lineages and projects (outgroup samples grey), palette
  ### for the other variables
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

message("tips: ", Ntip(tree), "   samples: ", sum(tip_data$is_sample),
        "   genomes: ", sum(tip_data$is_genome),
        "   clonemates: ", sum(tip_data$is_clone))



# annotations ------------------------------------------------------------------

### local bootstrap and clade size per internal node (root excluded)
desc <- prop.part(tree)

### lineages below each node, per level (assemblies and the reference carry none)
n_units <- function(kk) {
  lab <- setNames(labels[[kk]], labels$sample)
  map_int(desc, ~ n_distinct(na.omit(lab[attr(desc, "labels")[.x]])))
}

annot <- tibble(node    = Ntip(tree) + seq_len(Nnode(tree)),
                support = suppressWarnings(as.numeric(tree$node.label)),
                n_tips  = lengths(desc),
                in_og   = map_lgl(desc, ~ all(attr(desc, "labels")[.x] %in% outgroup)),
                !!!set_names(map(k_levels, n_units), paste0("units_", k_levels))) %>%
  filter(node != Ntip(tree) + 1) %>%
  mutate(boot_band   = case_when(is.na(support) ~ NA_character_,
                                 support < 70   ~ "< 70",
                                 support < 95   ~ "70 to 94",
                                 TRUE           ~ ">= 95"),
         boot_band   = factor(boot_band, levels = names(col_boot_band)))



# checks -----------------------------------------------------------------------

## check if the support values were parsed
if (all(is.na(annot$support))) stop("no numeric support labels on the tree")

## check if the tree holds the full panel
if (Ntip(tree) != n_tips_expected) {
  warning("tree has ", Ntip(tree), " tips, expected ", n_tips_expected)
}

## check if the tree is fully resolved
if (!is.binary(tree)) warning("CASTER-site tree has polytomies")

## check if only the terminal branches lack a length
if (n_no_length > Ntip(tree) + 1) {
  warning(n_no_length - Ntip(tree), " internal branches without a length (drawn as 0)")
}

## check if the outgroup is a clade
if (!all(outgroup %in% tree$tip.label)) {
  warning("outgroup tips missing from the tree")
} else if (!is.monophyletic(tree, outgroup)) {
  warning("the outgroup is not monophyletic in the rooted tree")
}

### position of the mapped outgroup samples
m4 <- grep("^M4-", tree$tip.label, value = TRUE)

if (length(m4) > 0) {
  message("mapped outgroup tips: ", length(m4),
          "   monophyletic: ", is.monophyletic(tree, m4),
          "   with the assemblies: ",
          is.monophyletic(tree, c(m4, intersect(outgroup, tree$tip.label))))
}

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



# tree -------------------------------------------------------------------------

### annotation per node of the drawn tree (nothing pruned)
plot_data <- bind_rows(
  tibble(node = seq_len(Ntip(tree)), label = tree$tip.label),
  tibble(node = annot$node, label = NA_character_)) %>%
  left_join(select(annot, node, support, boot_band), by = "node") %>%
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

band_num <- rep(NA_integer_, Ntip(tree) + Nnode(tree))
band_num[annot$node] <- match(as.character(annot$boot_band), band_levels)

plot_data <- plot_data %>%
  mutate(branch_band = factor(band_levels[inherit_band(tree, band_num)[node]],
                              levels = names(col_boot_band)))

### base tree for the tip-variable figures
base_tree <- function() {
  ggtree(tree, linewidth = 0.2) %<+% plot_data +
    ### tip names aligned as a column
    geom_tiplab(size = tip_size, align = TRUE, linesize = 0.1,
                linetype = "dotted") +
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

  save_fig(sprintf("caster_site_tree_%s.pdf", nm), p, tree_width, tree_height)
}


### support tree, value printed and branch coloured by its inherited band
p_tree_support <- ggtree(tree, linewidth = 0.3) %<+% plot_data +
  aes(colour = branch_band) +
  ### support values
  geom_text2(aes(subset = !isTip & !is.na(support) & support <= support_label_cut,
                 label = sprintf("%.0f", support)),
             size = label_size, hjust = label_hjust, vjust = -0.4, colour = "black") +
  scale_colour_manual(values = col_boot_band, name = "Local bootstrap",
                      drop = FALSE, na.value = "grey70") +
  geom_tiplab(size = tip_size, align = TRUE, linesize = 0.1,
              linetype = "dotted", colour = "black") +
  hexpand(label_space)

save_fig("caster_site_tree_support.pdf", p_tree_support, tree_width, tree_height)


### rooted fan (branch lengths dropped: the outgroup branch would take the
### whole radius)
p_fan <- ggtree(tree, layout = "fan", open.angle = 10, branch.length = "none",
                linewidth = 0.2) %<+% plot_data +
  geom_tippoint(aes(fill = cluster_plot), shape = 21, colour = "black",
                stroke = 0.1, size = 1.4, na.rm = TRUE) +
  ### tip labels along the radius
  geom_tiplab(size = 1.6, offset = 0.6, colour = "black") +
  hexpand(0.2) +
  cat_scale("cluster", name = "Lineage (K = 6)", aesthetics = "fill") +
  guides(fill = guide_legend(override.aes = list(size = 3)))

save_fig("caster_site_tree_fan.pdf", p_fan, 30, 30)



# summary figures --------------------------------------------------------------

### ingroup branches (the outgroup and the root split excluded)
ingroup_support <- annot %>%
  filter(!in_og, n_tips < Ntip(tree) - length(outgroup), !is.na(support))

### branch class per level (among lineages: samples of more than one lineage in the clade)
support_classes <- map_dfr(k_levels, function(kk) {
  ingroup_support %>%
    mutate(level = kk,
           class = if_else(.data[[paste0("units_", kk)]] > 1, "among lineages", "within a lineage"))
})

p_hist <- ggplot(ingroup_support, aes(x = support)) +
  geom_histogram(binwidth = 2, fill = "grey40") +
  geom_vline(xintercept = support_high, colour = "#D73027", linetype = "dashed",
             linewidth = 0.5) +
  labs(x = "Local bootstrap support (%)", y = "N internal branches")

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
  labs(x = "N tips in the clade", y = "Local bootstrap support (%)")

save_fig("caster_site_summary.pdf",
         (p_hist + p_size) + plot_annotation(tag_levels = "A"), 7.1, 3.8)



# tables -----------------------------------------------------------------------

### support summary for all ingroup branches and per level among and within lineages
summ <- function(d) {
  d %>%
    summarise(branches = n(),
              median   = median(support),
              q25      = quantile(support, 0.25, names = FALSE),
              q75      = quantile(support, 0.75, names = FALSE),
              high     = sum(support >= support_high),
              below_70 = sum(support < 70),
              lowest   = min(support),
              .by = c(level, class))
}

support_summary <- bind_rows(
  summ(mutate(ingroup_support, level = "all", class = "all ingroup branches")),
  summ(support_classes) %>% arrange(level, class))

publication_support <- support_summary %>%
  transmute(Level    = if_else(level == "all", "", lab_level[level]),
            Branches = class,
            N        = branches,
            `Local bootstrap, median [IQR]` = sprintf("%.1f [%.1f, %.1f]", median, q25, q75),
            !!sprintf("Local bootstrap >= %d", support_high) :=
              sprintf("%d (%.1f%%)", high, 100 * high / branches),
            Lowest   = sprintf("%.1f", lowest))

### topology and support only (CASTER-site leaves terminal branch lengths undefined)
tree_out <- tree
tree_out$edge.length <- NULL
write.tree(tree_out, file.path(out_root, "caster_site_tree_rooted.tre"))

write_tsv(select(annot, node, n_tips, in_og, starts_with("units_"), support, boot_band),
          file.path(out_root, "caster_site_support.tsv"))
write_tsv(support_summary, file.path(out_root, "caster_site_support_summary.tsv"))
write_tsv(publication_support, file.path(out_root, "publication_caster_site_support.tsv"))

cat("\n=== local bootstrap, ingroup ===\n")
print(publication_support, n = Inf, width = Inf)

cat("\n=== the", n_show, "largest ingroup clades ===\n")
ingroup_support %>%
  arrange(desc(n_tips)) %>%
  select(node, n_tips, support) %>%
  print(n = n_show)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  caster_site_tree_cluster.pdf           tips coloured by lineage, K = 6\n")
cat("  caster_site_tree_country.pdf           tips coloured by country\n")
cat("  caster_site_tree_site.pdf              tips coloured by sampling site\n")
cat("  caster_site_tree_project.pdf           tips coloured by project\n")
cat("  caster_site_tree_support.pdf           local bootstrap printed and coloured per branch\n")
cat("  caster_site_tree_fan.pdf               rooted fan, tips coloured by lineage\n")
cat("  caster_site_summary.pdf                support distribution, and support against clade size among and within lineages\n")
cat("  caster_site_tree_rooted.tre            rooted on the outgroup, support, no branch lengths, for 100_cluster_monophylly.R\n")
cat("  caster_site_support.tsv                local bootstrap, clade size and lineages below per branch\n")
cat("  caster_site_support_summary.tsv        support for all ingroup branches and per level among and within lineages\n")
cat("  publication_caster_site_support.tsv    the same, formatted\n")
