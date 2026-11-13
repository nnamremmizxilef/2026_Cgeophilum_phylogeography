## ============================================================================
## 081_aster.R
## The coalescent tree from ASTER, its branch annotations and the polytomy test
##
## Files needed:
##   aster_species_tree.tre              coalescent tree with annotations
##   aster_polytomy_test.tre             per branch p value against a polytomy
##                                       both from hpc/02_structure_phylogeny/081_aster.sh
##                                       (copy both files of
##                                        ${ROOT}/02_structure_phylogeny/results/081_results_aster
##                                        to local/02_structure_phylogeny/data/081_data_aster)
##   sample_metadata_ALL_FIXED.tsv       metadata and clone status, from local/01_variant_calling/103_mating_type.R
##   lineage_labels.tsv                  lineage per sample, from local/03_delimitation/001_lineage_labels.R
## ============================================================================

library(ape)
library(ggtree)
library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

data_root <- "02_structure_phylogeny/data/081_data_aster"
out_root  <- "02_structure_phylogeny/results/081_results_aster"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

f_tree    <- file.path(data_root, "aster_species_tree.tre")
f_poly    <- file.path(data_root, "aster_polytomy_test.tre")
f_meta    <- "01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv"
f_cluster <- "03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv"



# parameters -------------------------------------------------------------------

### define outgroup
outgroup <- c("PsefloM402_1", "PsefloM405_1", "PsefloM406_1",
              "M4-01", "M4-02", "M4-03", "M4-04", "M4-05", "M4-06")
genome_pattern <- "^(Cenge|Pseflo)"
reference_tip  <- "consensus"

### local posterior probability at or above this is taken as resolved
pp_high <- 0.95

### polytomy test threshold (branches not separable from a polytomy)
poly_alpha <- 0.05

### asymmetry test threshold for incomplete lineage sorting
asym_flag <- 0.05

### q1 of a branch without signal (one of three quartet resolutions)
q_null <- 1 / 3

### most asymmetric branches listed
n_show <- 20



# style ------------------------------------------------------------------------

### annotations
pp_label_cut <- 1
label_size   <- 0.9
label_hjust  <- 2.2
tip_size     <- 1.1
tree_width   <- 24
tree_height  <- 60
label_space  <- 0.06
clone_colour <- "grey70"

### posterior colors
col_pp_band <- c("< 0.7"      = "#D73027",
                 "0.7 to 0.9" = "#FDAE61",
                 "> 0.9"      = "#1A9850")

### polytomy test colors
col_poly_band <- c("not separable" = "#D73027",
                   "weak"          = "#FDAE61",
                   "resolved"      = "#1A9850")

col_flag <- c(`FALSE` = "grey60", `TRUE` = "#D73027")

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

### categorical tip colors
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

tree <- read.tree(f_tree)

### root on the outgroup (via an ingroup tip first, root() misplaces the root
### when the outgroup spans the stored root)
tree <- root(tree, outgroup = setdiff(tree$tip.label, outgroup)[1],
             resolve.root = TRUE, edgelabel = TRUE)
tree <- root(tree, outgroup = outgroup, resolve.root = TRUE, edgelabel = TRUE)
tree <- ladderize(tree)
meta <- read_tsv(f_meta, show_col_types = FALSE)

### lineages (K = 6)
clusters <- read_tsv(f_cluster, show_col_types = FALSE) %>%
  filter(clustered) %>%
  transmute(tip = sample, cluster = k6)

tip_data <- tibble(tip = tree$tip.label) %>%
  left_join(select(meta, tip = sample, project, country, site), by = "tip") %>%
  mutate(is_genome    = str_detect(tip, genome_pattern),
         is_reference = tip == reference_tip,
         is_sample    = !is_genome & !is_reference)

### lineage per tip
tip_data <- left_join(tip_data, clusters, by = "tip")

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

### column colors, clonemate level set to grey
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

message("tips: ", Ntip(tree), "   samples: ", sum(tip_data$is_sample),
        "   genomes: ", sum(tip_data$is_genome),
        "   clonemates: ", sum(tip_data$is_clone))



# annotations ------------------------------------------------------------------

### read ASTER blocks by name
field <- function(x, key) {
  m <- str_match(x, paste0("(?<![A-Za-z])", key,
                           "=(-?[0-9.]+(?:[eE][+-]?[0-9]+)?)"))[, 2]
  suppressWarnings(as.numeric(m))
}

annot <- tibble(node  = Ntip(tree) + seq_len(Nnode(tree)),
                label = tree$node.label) %>%
  mutate(q1  = field(label, "q1"),
         q2  = field(label, "q2"),
         q3  = field(label, "q3"),
         pp1 = field(label, "pp1"),
         f1  = field(label, "f1"),
         ### coalescent branch lengths (population-scaled time units)
         cu  = field(label, "CULength"),
         su  = field(label, "SULength"),
         asym = abs(q2 - q3),
         ### styling
         pp_band = case_when(pp1 <  0.70 ~ "< 0.7",
                             pp1 <= 0.90 ~ "0.7 to 0.9",
                             TRUE        ~ "> 0.9"),
         pp_band = factor(pp_band, levels = names(col_pp_band)))



# polytomy ---------------------------------------------------------------------

### match p values by unrooted split, keyed by the smaller side
split_key <- function(tr) {
  pp   <- prop.part(tr)
  labs <- attr(pp, "labels")
  vapply(pp, function(i) {
    a <- sort(labs[i])
    b <- sort(labs[-i])
    ### a single tip side carries no topology
    if (length(a) < 2 || length(b) < 2) return(NA_character_)
    ka <- paste(a, collapse = "|")
    kb <- paste(b, collapse = "|")
    if (length(a) < length(b)) ka else if (length(b) < length(a)) kb else min(ka, kb)
  }, character(1))
}

poly_tree <- read.tree(f_poly)

## check if the two trees describe the same set of tips
if (!setequal(poly_tree$tip.label, tree$tip.label)) {
  warning("the polytomy tree and the ASTER tree carry different tips, ",
          length(intersect(poly_tree$tip.label, tree$tip.label)), " shared")
}

poly_raw <- tibble(key     = split_key(poly_tree),
                   p_value = suppressWarnings(as.numeric(poly_tree$node.label))) %>%
  filter(!is.na(key)) %>%
  group_by(key) %>%
  summarise(p_value = p_value[which(!is.na(p_value))[1]], .groups = "drop")

annot <- annot %>%
  mutate(key = split_key(tree)) %>%
  left_join(poly_raw, by = "key") %>%
  select(-key)

annot <- annot %>%
  mutate(poly_band = case_when(is.na(p_value)        ~ NA_character_,
                               p_value >= poly_alpha ~ "not separable",
                               p_value >= 0.01       ~ "weak",
                               TRUE                  ~ "resolved"),
         poly_band = factor(poly_band, levels = names(col_poly_band)))

message("polytomy p values matched: ", sum(!is.na(annot$p_value)),
        " of ", nrow(annot), " internal branches")



# checks -----------------------------------------------------------------------

## check if the annotations were parsed
if (all(is.na(annot$pp1))) {
  stop("no pp1 values parsed, check that ASTER was run with -u 2")
}

## check if the coalescent branch lengths are there (CASTLES-2)
if (all(is.na(annot$cu))) {
  warning("no CULength values parsed, the coalescent branch length statistics ",
          "will be empty")
}

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

### annotations on the nodes of the drawn tree (nothing pruned)
plot_data <- bind_rows(
  tibble(node = seq_len(Ntip(tree)), label = tree$tip.label),
  tibble(node = annot$node, label = NA_character_)) %>%
  left_join(select(annot, node, pp1, pp_band, q1, asym, cu, p_value, poly_band),
            by = "node") %>%
  left_join(tip_data, by = c("label" = "tip"))

### weakest band on the path from the root is carried downstream
band_levels <- c("> 0.9", "0.7 to 0.9", "< 0.7")

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
band_num[annot$node] <- match(as.character(annot$pp_band), band_levels)

plot_data <- plot_data %>%
  mutate(branch_band = factor(band_levels[inherit_band(tree, band_num)[node]],
                              levels = names(col_pp_band)))

### the same for the polytomy test (weakest result on the path from the root)
poly_levels <- c("resolved", "weak", "not separable")

poly_num <- rep(NA_integer_, Ntip(tree) + Nnode(tree))
poly_num[annot$node] <- match(as.character(annot$poly_band), poly_levels)

plot_data <- plot_data %>%
  mutate(branch_poly = factor(poly_levels[inherit_band(tree, poly_num)[node]],
                              levels = names(col_poly_band)))

### the tip variable carries the figure, the polytomy test has its own tree
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

  save_fig(sprintf("aster_tree_%s.pdf", nm), p, tree_width, tree_height)
}


### support tree, value printed and branch coloured by its inherited band
p_tree_support <- ggtree(tree, linewidth = 0.3) %<+% plot_data +
  aes(colour = branch_band) +
  ### posterior values
  geom_text2(aes(subset = !isTip & !is.na(pp1) & pp1 <= pp_label_cut,
                 label = sprintf("%.2f", pp1)),
             size = label_size, hjust = label_hjust, vjust = -0.4, colour = "black") +
  scale_colour_manual(values = col_pp_band, name = "Local posterior",
                      drop = FALSE, na.value = "grey70") +
  geom_tiplab(size = tip_size, align = TRUE, linesize = 0.1,
              linetype = "dotted", colour = "black") +
  hexpand(label_space)

save_fig("aster_tree_support.pdf", p_tree_support, tree_width, tree_height)


### polytomy tree, p value printed and branch coloured by its inherited band
p_tree_poly <- ggtree(tree, linewidth = 0.3) %<+% plot_data +
  aes(colour = branch_poly) +
  ### p values
  geom_text2(aes(subset = !isTip & !is.na(p_value),
                 label = sprintf("%.2f", p_value)),
             size = label_size, hjust = label_hjust, vjust = -0.4, colour = "black") +
  scale_colour_manual(values = col_poly_band, name = "Polytomy test",
                      drop = FALSE, na.value = "grey70") +
  geom_tiplab(size = tip_size, align = TRUE, linesize = 0.1,
              linetype = "dotted", colour = "black") +
  hexpand(label_space)

save_fig("aster_tree_polytomy.pdf", p_tree_poly, tree_width, tree_height)


### rooted fan (branch lengths dropped, the outgroup branch would take the whole radius)
p_fan <- ggtree(tree, layout = "fan", open.angle = 10, branch.length = "none",
                linewidth = 0.2) %<+% plot_data +
  geom_tippoint(aes(fill = cluster_plot), shape = 21, colour = "black",
                stroke = 0.1, size = 1.4, na.rm = TRUE) +
  ### labels follow the radius, the canvas carries 618 of them
  geom_tiplab(size = 1.6, offset = 0.6, colour = "black") +
  hexpand(0.2) +
  cat_scale("cluster", name = "Lineage (K = 6)", aesthetics = "fill") +
  guides(fill = guide_legend(override.aes = list(size = 3)))

save_fig("aster_tree_fan.pdf", p_fan, 30, 30)



# summary figures --------------------------------------------------------------

pp_summary <- annot %>%
  filter(!is.na(pp1)) %>%
  summarise(branches = n(),
            median   = median(pp1),
            below_95 = sum(pp1 < pp_high),
            below_70 = sum(pp1 < 0.70))

p_pp <- ggplot(filter(annot, !is.na(pp1)), aes(x = pp1)) +
  geom_histogram(binwidth = 0.02, fill = "grey40") +
  geom_vline(xintercept = pp_high, colour = "#D73027", linetype = "dashed",
             linewidth = 0.5) +
  labs(x = "Local posterior probability", y = "N internal branches")

### quartet frequency against local posterior
p_q1 <- ggplot(filter(annot, !is.na(q1)), aes(x = q1, y = pp1)) +
  geom_hline(yintercept = pp_high, colour = "grey50", linetype = "dotted") +
  geom_vline(xintercept = q_null, colour = "#D73027", linetype = "dashed") +
  geom_point(shape = 21, colour = "black", fill = "#4575B4", size = 1, alpha = 0.5) +
  labs(x = "q1", y = "Local posterior probability")

### asymmetry of the alternative resolutions (not produced by incomplete lineage sorting alone)
p_asym <- ggplot(filter(annot, !is.na(q2)), aes(x = q2, y = q3)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dashed") +
  geom_point(aes(fill = asym > asym_flag), shape = 21, colour = "black", size = 1,
             alpha = 0.6) +
  scale_fill_manual(values = col_flag, name = NULL,
                    labels = c(`FALSE` = "symmetric",
                               `TRUE` = sprintf("|q2 - q3| > %.2f", asym_flag))) +
  labs(x = "q2", y = "q3")

p_poly <- ggplot(filter(annot, !is.na(p_value)), aes(x = p_value)) +
  geom_histogram(binwidth = 0.02, fill = "grey40") +
  geom_vline(xintercept = poly_alpha, colour = "#D73027",
             linetype = "dashed", linewidth = 0.5) +
  labs(x = "p value against a polytomy", y = "N internal branches")

save_fig("aster_summary.pdf",
         (p_pp + p_q1) / (p_asym + p_poly) +
           plot_annotation(tag_levels = "A"), 7.1, 5.9)



# tables -----------------------------------------------------------------------

cu_ok <- filter(annot, !is.na(cu), cu > 0)

poly_summary <- annot %>%
  filter(!is.na(p_value)) %>%
  summarise(branches     = n(),
            rejected     = sum(p_value < poly_alpha),
            not_rejected = sum(p_value >= poly_alpha))

write.tree(tree, file.path(out_root, "aster_tree_rooted.tre"))

write_tsv(select(annot, node, pp1, q1, q2, q3, asym, f1, cu, su, p_value),
          file.path(out_root, "aster_branch_annotations.tsv"))

write_tsv(pp_summary, file.path(out_root, "aster_pp_summary.tsv"))

aster_statistics <- bind_rows(
  tibble(quantity = "internal branches",            value = pp_summary$branches),
  tibble(quantity = "median local posterior",       value = round(pp_summary$median, 4)),
  tibble(quantity = "branches below 0.95",          value = pp_summary$below_95),
  tibble(quantity = "median q1",                    value = round(median(annot$q1, na.rm = TRUE), 4)),
  tibble(quantity = "branches with |q2 - q3| > flag", value = sum(annot$asym > asym_flag, na.rm = TRUE)),
  if (nrow(cu_ok) > 0) {
    tibble(quantity = "median coalescent branch length", value = round(median(cu_ok$cu), 4))
  },
  if (nrow(cu_ok) > 0) {
    tibble(quantity = "branches below 0.1 coalescent units", value = sum(cu_ok$cu < 0.1))
  },
  tibble(quantity = "branches with a polytomy p value", value = poly_summary$branches),
  tibble(quantity = "branches not separable from a polytomy", value = poly_summary$not_rejected))

write_tsv(aster_statistics, file.path(out_root, "aster_statistics.tsv"))

cat("\n=== branch support ===\n")
print(pp_summary, n = Inf)

cat("\n=== statistics ===\n")
print(aster_statistics, n = Inf)

cat("\n=== polytomy test ===\n")
print(poly_summary, n = Inf)

cat("\n=== the", n_show, "most asymmetric branches ===\n")
annot %>%
  filter(!is.na(asym)) %>%
  arrange(desc(asym)) %>%
  select(node, pp1, q1, q2, q3, asym, cu, p_value) %>%
  print(n = n_show)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  aster_tree_cluster.pdf     tips coloured by lineage, K = 6\n")
cat("  aster_tree_country.pdf     tips coloured by country\n")
cat("  aster_tree_site.pdf        tips coloured by sampling site\n")
cat("  aster_tree_project.pdf     tips coloured by project\n")
cat("  aster_tree_support.pdf     local posterior printed and coloured per branch\n")
cat("  aster_tree_polytomy.pdf    polytomy p value printed and coloured per branch\n")
cat("  aster_tree_fan.pdf         rooted fan, tips coloured by lineage\n")
cat("  aster_summary.pdf          support, quartets and the polytomy test\n")
cat("  aster_tree_rooted.tre      rooted on the outgroup, for 100\n")
cat("  aster_branch_annotations.tsv\n")
cat("  aster_pp_summary.tsv\n")
cat("  aster_statistics.tsv\n")
