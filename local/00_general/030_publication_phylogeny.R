## ============================================================================
## 030_publication_phylogeny.R
## Publication figure of the phylogeny and genetic clustering: concatenation tree
## against coalescent tree on their shared tips, with the support of unit stems
## and of branches between units (A), DAPC BIC and sNMF cross-entropy against K
## (B), PCA of the clustering SNP set (C) and sNMF ancestry at K = 4, 5 and 6 (D)
##
## Files needed:
##   species_tree_grafted_rooted.tre  rooted concatenation tree of the MLG representatives with UFBoot,
##                                    clonemates grafted at zero length (every tip drawn),
##                                    from local/02_structure_phylogeny/052_species_tree_decloned.R
##   aster_tree_rooted.tre            rooted ASTER tree with local posterior,
##                                    from local/02_structure_phylogeny/081_aster.R
##   aster_collapsed_tree_rooted.tre  rooted ASTER tree of the collapsed gene trees (option),
##                                    from local/02_structure_phylogeny/080_aster_collapsed.R
##   caster_site_tree_rooted.tre      rooted CASTER-site tree (option),
##                                    from local/02_structure_phylogeny/082_caster_site.R
##   site_concordance.cf.stat, site_concordance.cf.branch
##                                    sCF per branch of the concatenation tree,
##                                    from hpc/02_structure_phylogeny/060_site_concordance.sh
##                                    (copied to local/02_structure_phylogeny/data/060_data_site_concordance
##                                     for local/02_structure_phylogeny/100_cluster_monophylly.R)
##   gene_concordance_<tree>.cf.stat, gene_concordance_<tree>.cf.branch
##                                    gCF per branch of each tree,
##                                    from hpc/02_structure_phylogeny/090_gene_concordance.sh
##                                    (copied to local/02_structure_phylogeny/data/090_data_gene_concordance
##                                     for local/02_structure_phylogeny/100_cluster_monophylly.R)
##   k_diagnostics.tsv                BIC and mean cross-entropy per K of the primary draw
##   replicate_curves.tsv             the same for every replicate draw
##   clusters_K<k>.tsv                sNMF max q per sample at K = k, in the sample order of the sNMF input
##   input.geno                       haploid genotypes of the sNMF input (loci x samples, 9 missing)
##   input.snmf/K<k>/run<r>/*.Q       sNMF ancestry per K and run
##                                    (all five from local/02_structure_phylogeny/014_clustering.R,
##                                     results/014_results_clustering)
##   lineage_labels.tsv               unit at K = 4, 5 and 6 and project per sample,
##                                    from local/03_delimitation/001_lineage_labels.R
## ============================================================================

library(ape)
library(ggtree)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(ggrepel)



# files ------------------------------------------------------------------------

res_02   <- "02_structure_phylogeny/results"
res_03   <- "03_delimitation/results"
dir_014  <- file.path(res_02, "014_results_clustering")
out_root <- "00_general/results/030_results_publication_phylogeny"

f_trees <- c(concat          = file.path(res_02, "052_results_species_tree_decloned", "species_tree_grafted_rooted.tre"),
             aster           = file.path(res_02, "081_results_aster", "aster_tree_rooted.tre"),
             aster_collapsed = file.path(res_02, "080_results_aster_collapsed", "aster_collapsed_tree_rooted.tre"),
             caster          = file.path(res_02, "082_results_caster_site", "caster_site_tree_rooted.tre"))
f_scf    <- file.path("02_structure_phylogeny/data/060_data_site_concordance", "site_concordance")
f_gcf    <- setNames(file.path("02_structure_phylogeny/data/090_data_gene_concordance",
                               paste0("gene_concordance_", names(f_trees))), names(f_trees))
f_kdiag  <- file.path(dir_014, "summary", "k_diagnostics.tsv")
f_kreps  <- file.path(dir_014, "summary", "replicate_curves.tsv")
f_geno   <- file.path(dir_014, "snmf", "input.geno")
dir_q    <- file.path(dir_014, "snmf", "input.snmf")
f_labels <- file.path(res_03, "001_results_lineage_labels", "lineage_labels.tsv")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### trees left and right (names of f_trees)
tree_left  <- "concat"
tree_right <- "aster"

### level of the colours and the ordering; levels with stem support and ancestry
k_main   <- "k6"
k_levels <- c("k4", "k5", "k6")

### K range of the selection curves and the levels shaded there
k_show  <- 2:12
k_shade <- c(4, 6)

### gap between the two trees (share of the tree depth)
tree_gap <- 0.55

### tree labels below the trees and stem support format per tree (support, sCF, gCF)
lab_tree <- c(concat          = "Concatenation (IQ-TREE): UFBoot / sCF / gCF",
              aster           = "Coalescent (ASTER): LPP / gCF",
              aster_collapsed = "Coalescent, collapsed gene trees (ASTER): LPP / gCF",
              caster          = "Coalescent, sites (CASTER): local bootstrap / gCF")
fmt_support <- list(concat          = function(s, scf, gcf) sprintf("%s/%s/%s", f0(s), f0(scf), f0(gcf)),
                    aster           = function(s, scf, gcf) sprintf("%s/%s", f2(s), f0(gcf)),
                    aster_collapsed = function(s, scf, gcf) sprintf("%s/%s", f2(s), f0(gcf)),
                    caster          = function(s, scf, gcf) sprintf("%s/%s", f0(s), f0(gcf)))

### largest mean difference in max q between a run and 014 accepted as the same run
q_tol <- 1e-3



# style ------------------------------------------------------------------------

### unit colours
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

### project labels and symbols; edges above tips of more than one unit
lab_project   <- c(magic = "MAGIC", phytoakmeter = "PhytOakmeter")
shape_project <- c(magic = 21, phytoakmeter = 24)
col_mixed     <- "grey60"

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}

### level labels; support values (-: not applicable)
lab_k <- function(k) sub("k", "K = ", k)
f0    <- function(x) if_else(is.na(x), "-", sprintf("%.0f", x))
f2    <- function(x) if_else(is.na(x), "-", sprintf("%.2f", x))



# read -------------------------------------------------------------------------

## check if every input exists
f_need <- c(f_trees[c(tree_left, tree_right)], f_kdiag, f_geno, f_labels)
miss   <- f_need[!file.exists(f_need)]
if (length(miss) > 0) stop("missing inputs: ", paste(miss, collapse = ", "))

labels <- read_tsv(f_labels, col_types = cols(.default = "c"))
kdiag  <- read_tsv(f_kdiag, show_col_types = FALSE)
kreps  <- read_tsv(f_kreps, show_col_types = FALSE)

### sNMF max q per level in the sample order of the sNMF input
qmax <- map(set_names(k_levels), function(kk) {
  read_tsv(file.path(dir_014, "summary", sprintf("clusters_K%s.tsv", sub("k", "", kk))), show_col_types = FALSE) %>%
    select(sample, snmf_max_q)
})
samples_snmf <- qmax[[1]]$sample

unit_of <- setNames(labels[[k_main]], labels$sample)

### concordance factors per branch (branch IDs are the node labels of the .cf.branch tree)
read_cf <- function(prefix) {
  list(stat = read_tsv(paste0(prefix, ".cf.stat"), comment = "#", col_types = cols(.default = "c")),
       tree = read.tree(paste0(prefix, ".cf.branch")))
}
cf_s <- read_cf(f_scf)
cf_g <- map(f_gcf[c(tree_left, tree_right)], read_cf)



# checks -----------------------------------------------------------------------

## check if the support formats cover both trees
if (!all(c(tree_left, tree_right) %in% names(fmt_support))) stop("no support format for ", tree_right)

## check if the clustering tables hold the same samples in the same order
if (!all(map_lgl(qmax, ~ identical(.x$sample, samples_snmf)))) stop("clusters_K tables differ in samples or order")

## check if every clustered sample has a unit
if (anyNA(unit_of[samples_snmf])) warning(sum(is.na(unit_of[samples_snmf])), " clustered samples without a unit")



# trees (A) --------------------------------------------------------------------

### both trees on their shared labelled tips; the left one ladderized, the right one rotated to its tip order
t1_all <- read.tree(f_trees[[tree_left]])
t2_all <- read.tree(f_trees[[tree_right]])
keep <- Reduce(intersect, list(t1_all$tip.label, t2_all$tip.label, labels$sample))
t1 <- drop.tip(t1_all, setdiff(t1_all$tip.label, keep))
t2 <- drop.tip(t2_all, setdiff(t2_all$tip.label, keep))

### tips drawn (and tips of each tree)
message("tips drawn: ", length(keep), " (left ", Ntip(t1_all), ", right ", Ntip(t2_all), ")")

### bipartition key over a tip set (the smaller side), to match branches between a tree and its .cf.branch tree
split_key <- function(x, tips) {
  a  <- sort(intersect(x, tips))
  b  <- sort(setdiff(tips, a))
  ka <- paste(a, collapse = "|")
  kb <- paste(b, collapse = "|")
  if (length(a) < length(b) || (length(a) == length(b) && ka < kb)) ka else kb
}

### support per clade: node label of the tree (local posterior pp1 for ASTER, bootstrap otherwise), sCF and gCF
### (no sCF for the coalescent trees)
support_of <- function(tr, tree_name) {
  pp     <- prop.part(tr)
  labs   <- attr(pp, "labels")
  clades <- map(pp, ~ labs[.x])
  nl     <- if (is.null(tr$node.label)) rep(NA_character_, tr$Nnode) else tr$node.label
  pp1    <- str_match(nl, "pp1=(-?[0-9.]+(?:[eE][+-]?[0-9]+)?)")[, 2]
  cf_val <- function(cf, col) {
    if (is.null(cf) || !col %in% names(cf$stat)) return(rep(NA_real_, length(clades)))
    common <- intersect(labs, cf$tree$tip.label)
    cpp    <- prop.part(cf$tree)
    ck     <- map_chr(cpp, ~ split_key(attr(cpp, "labels")[.x], common))
    v      <- setNames(as.numeric(cf$stat[[col]][match(cf$tree$node.label, cf$stat$ID)]), ck)
    unname(v[map_chr(clades, ~ split_key(.x, common))])
  }
  tibble(key     = map_chr(clades, ~ paste(sort(.x), collapse = "|")),
         support = suppressWarnings(as.numeric(coalesce(pp1, str_extract(nl, "^[0-9.]+")))),
         sCF     = cf_val(if (tree_name == "concat") cf_s else NULL, "sCF"),
         gCF     = cf_val(cf_g[[tree_name]], "gCF"))
}
sup1 <- support_of(t1, tree_left)
sup2 <- support_of(t2, tree_right)

t1 <- ladderize(t1)

d1  <- ggtree(t1, branch.length = "none", ladderize = FALSE)$data
ord <- d1 %>% filter(isTip) %>% arrange(y) %>% pull(label)
t2  <- rotateConstr(t2, ord)
d2  <- ggtree(t2, branch.length = "none", ladderize = FALSE)$data

### edge colour: the unit when all tips below the edge belong to one unit
node_unit <- function(tr) {
  u  <- unit_of[tr$tip.label]
  pp <- prop.part(tr)
  tibble(node = seq_len(Ntip(tr) + tr$Nnode),
         grp  = c(coalesce(u, "mixed"),
                  map_chr(pp, function(i) { v <- unique(u[i]); if (length(v) == 1 && !is.na(v)) v else "mixed" })))
}
d1 <- left_join(d1, node_unit(t1), by = "node")
d2 <- left_join(d2, node_unit(t2), by = "node")

### right tree mirrored next to the left one
x1  <- max(d1$x)
gap <- tree_gap * x1
d2  <- mutate(d2, x = x1 + gap + (max(x) - x))
xs1 <- x1 + 1.2
xs2 <- x1 + gap - 1.2

tips <- d1 %>%
  filter(isTip) %>%
  select(label, y1 = y) %>%
  inner_join(d2 %>% filter(isTip) %>% select(label, y2 = y), by = "label") %>%
  mutate(unit = coalesce(unit_of[label], "mixed"))

### labelled branches: the stem of every unit (K = 4 to 6) that is a clade, and every branch between units of the
### main level (root excluded)
stems_of <- function(tr, d, tree_name, sup) {
  pp    <- prop.part(tr)
  labs  <- attr(pp, "labels")
  n     <- Ntip(tr)
  n_lin <- map_int(pp, ~ n_distinct(unit_of[labs[.x]]))
  among <- setdiff(n + which(n_lin >= 2), n + 1)
  stems <- unlist(map(k_levels, function(kk) {
    map_int(unique(labels[[kk]]), function(uu) {
      m <- intersect(labels$sample[labels[[kk]] == uu], tr$tip.label)
      if (length(m) < 2) return(NA_integer_)
      nd <- getMRCA(tr, m)
      if (length(pp[[nd - n]]) == length(m)) nd else NA_integer_
    })
  }))
  tibble(node = sort(unique(c(among, na.omit(stems))))) %>%
    mutate(key = map_chr(node, ~ paste(sort(labs[pp[[.x - n]]]), collapse = "|"))) %>%
    left_join(sup, by = "key") %>%
    mutate(label = fmt_support[[tree_name]](support, sCF, gCF)) %>%
    left_join(select(d, node, x, y), by = "node")
}
st1 <- stems_of(t1, d1, tree_left, sup1)
st2 <- stems_of(t2, d2, tree_right, sup2)

pal_edge <- c(pal_unit, mixed = col_mixed)
y_lab    <- -0.03 * max(d1$y)

p_trees <- ggplot() +
  geom_tree(data = d1, aes(x, y, colour = grp), linewidth = 0.15) +
  geom_tree(data = d2, aes(x, y, colour = grp), linewidth = 0.15) +
  geom_segment(data = tips, aes(x = xs1 + 0.8, xend = xs2 - 0.8, y = y1, yend = y2, colour = unit),
               linewidth = 0.12, alpha = 0.7) +
  geom_tile(data = tips, aes(x = xs1, y = y1, fill = unit), width = 1.4, height = 1) +
  geom_tile(data = tips, aes(x = xs2, y = y2, fill = unit), width = 1.4, height = 1) +
  geom_label_repel(data = st1, aes(x, y, label = label), size = 5.5 / .pt, label.size = 0.1,
                   label.padding = unit(0.05, "lines"), box.padding = 0.06, direction = "y", nudge_x = -1.5,
                   hjust = 1, min.segment.length = 0, segment.size = 0.2, segment.colour = "grey40", seed = 1) +
  geom_label_repel(data = st2, aes(x, y, label = label), size = 5.5 / .pt, label.size = 0.1,
                   label.padding = unit(0.05, "lines"), box.padding = 0.06, direction = "y", nudge_x = 1.5,
                   hjust = 0, min.segment.length = 0, segment.size = 0.2, segment.colour = "grey40", seed = 1) +
  annotate("text", x = c(x1 / 2, x1 + gap + x1 / 2), y = y_lab, label = lab_tree[c(tree_left, tree_right)],
           size = 6.5 / .pt) +
  scale_colour_manual(values = pal_edge, guide = "none") +
  scale_fill_manual(values = pal_edge, breaks = intersect(names(pal_unit), tips$unit), name = NULL) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(colour = "black", linewidth = 0.2))) +
  coord_cartesian(clip = "off") +
  theme_tree() +
  theme(legend.position = "bottom", legend.key.size = unit(3, "mm"), legend.text = element_text(size = 7))



# K selection (B) --------------------------------------------------------------

lab_stat <- c(BIC = "DAPC, BIC", mean_ce = "sNMF, cross-entropy")

curve_long <- function(d) {
  d %>%
    filter(K %in% k_show) %>%
    pivot_longer(c(BIC, mean_ce), names_to = "stat", values_to = "value") %>%
    mutate(stat = factor(lab_stat[stat], levels = lab_stat))
}
prim <- curve_long(select(kdiag, K, BIC, mean_ce))
reps <- curve_long(filter(kreps, draw != "primary"))

p_k <- ggplot(prim, aes(K, value)) +
  annotate("rect", xmin = k_shade[1] - 0.4, xmax = k_shade[2] + 0.4, ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_line(data = reps, aes(group = draw), colour = "grey65", linewidth = 0.2) +
  geom_line(linewidth = 0.5) +
  geom_point(shape = 21, size = 1.2, stroke = 0.3, fill = "white", colour = "black") +
  facet_wrap(~ stat, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_x_continuous(breaks = k_show) +
  labs(x = "Number of clusters (K)", y = NULL) +
  theme(strip.placement = "outside", strip.background = element_blank())



# PCA (C) ----------------------------------------------------------------------

### genotypes as samples x loci, missing at the locus mean, centred
geno <- read_lines(f_geno)
G <- do.call(rbind, strsplit(geno, "", fixed = TRUE))
storage.mode(G) <- "integer"
G[G == 9L] <- NA_integer_

## check if the genotype matrix holds the clustered samples
if (ncol(G) != length(samples_snmf)) stop("input.geno has ", ncol(G), " samples, clusters_K ", length(samples_snmf))

X <- t(G)
X <- apply(X, 2, function(v) replace(v, is.na(v), mean(v, na.rm = TRUE)))
X <- X[, apply(X, 2, var) > 0, drop = FALSE]
pc <- prcomp(X, center = TRUE, scale. = FALSE, rank. = 2)
pct <- 100 * pc$sdev^2 / sum(pc$sdev^2)

pca <- tibble(sample = samples_snmf, PC1 = pc$x[, 1], PC2 = pc$x[, 2]) %>%
  left_join(select(labels, sample, project, unit = all_of(k_main)), by = "sample") %>%
  mutate(unit = factor(unit, levels = intersect(names(pal_unit), unit)))

p_pca <- ggplot(pca, aes(PC1, PC2, fill = unit, shape = project)) +
  geom_point(size = 1.4, stroke = 0.2, colour = "black") +
  scale_fill_manual(values = pal_unit, guide = "none") +
  scale_shape_manual(values = shape_project, labels = lab_project, name = NULL) +
  labs(x = sprintf("PC1 (%.1f%%)", pct[1]), y = sprintf("PC2 (%.1f%%)", pct[2])) +
  theme(legend.position = "inside", legend.position.inside = c(0.99, 0.99), legend.justification = c(1, 1),
        legend.background = element_blank(), legend.key.size = unit(3, "mm"))



# ancestry (D) -----------------------------------------------------------------

### Q matrix of the run matching the max q of 014, components named after the unit with the highest mean share
q_level <- function(kk) {
  k  <- as.integer(sub("k", "", kk))
  fq <- list.files(file.path(dir_q, paste0("K", k)), pattern = "\\.Q$", recursive = TRUE, full.names = TRUE)
  if (length(fq) == 0) stop("no sNMF Q files for K = ", k, " in ", dir_q)
  qs   <- map(fq, ~ as.matrix(read.table(.x)))
  dev  <- map_dbl(qs, ~ mean(abs(apply(.x, 1, max) - qmax[[kk]]$snmf_max_q)))
  best <- which.min(dev)

  ## check if the run matches the max q written by 014
  if (dev[best] > q_tol) warning("K = ", k, ": no sNMF run matches clusters_K", k, " (mean difference ",
                                 signif(dev[best], 2), ")")

  q    <- qs[[best]]
  u    <- labels[[kk]][match(samples_snmf, labels$sample)]
  comp <- map_chr(seq_len(ncol(q)), function(j) names(which.max(tapply(q[, j], u, mean))))

  ## check if every component maps to a different unit
  if (anyDuplicated(comp)) warning("K = ", k, ": sNMF components map to the same unit: ", paste(comp, collapse = ", "))

  colnames(q) <- comp
  as_tibble(q) %>%
    mutate(sample = samples_snmf, level = lab_k(kk)) %>%
    pivot_longer(-c(sample, level), names_to = "comp", values_to = "share")
}
q_long <- map_dfr(k_levels, q_level)

### samples by unit of the main level, then by their share of that unit
q_main  <- filter(q_long, level == lab_k(k_main))
order_s <- tibble(sample = samples_snmf, unit = unit_of[samples_snmf]) %>%
  left_join(q_main, by = c("sample", "unit" = "comp")) %>%
  mutate(unit = factor(unit, levels = intersect(names(pal_unit), unit))) %>%
  arrange(unit, desc(share)) %>%
  mutate(i = row_number()) %>%
  select(sample, i, unit)
brk <- summarise(order_s, lo = min(i), hi = max(i), .by = unit) %>% mutate(mid = (lo + hi) / 2)

q_long <- q_long %>%
  inner_join(select(order_s, sample, i), by = "sample") %>%
  mutate(level = factor(level, levels = lab_k(k_levels)), comp = factor(comp, levels = names(pal_unit)))

p_q <- ggplot(q_long, aes(i, share, fill = comp)) +
  geom_col(width = 1, linewidth = 0) +
  geom_vline(xintercept = brk$hi[-nrow(brk)] + 0.5, colour = "white", linewidth = 0.5) +
  facet_grid(level ~ ., switch = "y") +
  scale_fill_manual(values = pal_unit, guide = "none") +
  scale_x_continuous(breaks = brk$mid, labels = brk$unit, expand = c(0, 0), position = "top") +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  labs(x = NULL, y = "Ancestry (sNMF)") +
  theme(strip.placement = "outside", strip.background = element_blank(), axis.ticks.x = element_blank(),
        axis.line.x = element_blank(), panel.spacing = unit(1, "mm"))



# figures ----------------------------------------------------------------------

p_fig <- p_trees / (p_k | p_pca) / p_q +
  plot_layout(heights = c(4.2, 2.1, 1.9)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold"))

save_fig(sprintf("publication_phylogeny_%s_%s.pdf", tree_left, tree_right), p_fig, 7.1, 9.6)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat(sprintf("  publication_phylogeny_%s_%s.pdf\n", tree_left, tree_right))
