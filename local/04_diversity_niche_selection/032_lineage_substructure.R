## ============================================================================
## 032_lineage_substructure.R
## Substructure within L2a and L3b at K = 6 (the units with depth-structured relatedness in
## 031): principal coordinates of the sequence distances between their MLGs, clusters by
## partitioning around medoids (number by the average silhouette width), NJ trees, and the
## association of clusters and principal coordinates with sampling depth inside the
## PhytOakmeter sites, tested by permuting depth labels within site and direction
##
## Files needed:
##   classes.tsv                    representatives (one per MLG), site, units
##   pair_differences.tsv.gz        differences per pair of representatives over all loci
##   pair_sites.tsv.gz              compared sites per pair of representatives over all loci
##                                  (all three from hpc/04_diversity_niche_selection/010_sequence_distances.sh;
##                                   copy ${ROOT}/04_diversity_niche_selection/results/010_results_sequence_distances
##                                   to local/04_diversity_niche_selection/data/010_data_sequence_distances)
##   lineage_labels.tsv             unit at K = 4, 5 and 6 per ingroup sample,
##                                  from local/03_delimitation/001_lineage_labels.R
##   sample_metadata_ALL_FIXED.tsv  project, site, direction, sampling_depth and mlg per sample,
##                                  from local/01_variant_calling/103_mating_type.R
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)
library(cluster)
library(ape)



# files ------------------------------------------------------------------------

data_root   <- "04_diversity_niche_selection/data/010_data_sequence_distances"
labels_file <- "03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv"
meta_file   <- "01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv"
out_root    <- "04_diversity_niche_selection/results/032_results_lineage_substructure"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### level and units examined, project sampled at several depths and directions
k_focal     <- "k6"
units_focal <- c("L2a", "L3b")
depth_proj  <- "phytoakmeter"

### numbers of clusters compared
k_range <- 2:12

### tolerance below the largest average silhouette width (the smallest number of clusters within it is chosen)
sil_tol <- 0.02

### smallest average silhouette width for substantial structure
sil_min <- 0.25

### principal coordinates tested against depth
n_axes  <- 3

### permutations, seed
n_perm <- 999
seed   <- 42

### depth intervals
lab_depth <- c(`1` = "5-10 cm", `2` = "25-30 cm", `3` = "45-50 cm")



# style ------------------------------------------------------------------------

### depth class colours of the MLGs (where they were found)
pal_depth <- c(`5-10 cm` = "#C6DBEF", `25-30 cm` = "#6BAED6", `45-50 cm` = "#08519C",
               `Several depths` = "grey55", `MAGIC (0-10 cm)` = "white")

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}

### p values
fmt_p <- function(p) if_else(p < 0.001, "< 0.001", sprintf("%.3f", p))



# read -------------------------------------------------------------------------

classes <- read_tsv(file.path(data_root, "classes.tsv"), col_types = cols(.default = "c", mixed = "l"))
lab     <- read_tsv(labels_file, col_types = cols(.default = "c"))
meta    <- read_tsv(meta_file, col_types = cols(.default = "c"))

### distance matrices (representatives in rows and columns)
read_mat <- function(f) {
  m <- read.delim(gzfile(f), row.names = 1, check.names = FALSE)
  as.matrix(m)
}
diffs     <- read_mat(file.path(data_root, "pair_differences.tsv.gz"))
cmp_sites <- read_mat(file.path(data_root, "pair_sites.tsv.gz"))

### p-distance between MLGs (via their representatives)
rep_mlg <- meta$mlg[match(rownames(diffs), meta$sample)]
p_dist  <- diffs / cmp_sites
diag(p_dist) <- 0
dimnames(p_dist) <- list(rep_mlg, rep_mlg)

### isolates of the examined units with project, site, direction, depth and MLG
iso <- lab %>%
  select(sample, unit = all_of(k_focal)) %>%
  filter(unit %in% units_focal) %>%
  inner_join(select(meta, sample, project, site, direction, depth = sampling_depth, mlg), by = "sample")



# checks -----------------------------------------------------------------------

## check if the matrices hold the representatives of classes.tsv in the same order
if (!identical(rownames(diffs), classes$sample) || !identical(colnames(diffs), classes$sample)) {
  stop("distance matrices and classes.tsv differ in samples or order")
}

## check if every representative has one MLG
if (anyNA(rep_mlg) || anyDuplicated(rep_mlg)) {
  stop("representatives without MLG or MLGs with several representatives")
}

## check if every MLG of the examined units has a representative
no_rep <- setdiff(unique(iso$mlg), rep_mlg)

if (length(no_rep) > 0) {
  warning(length(no_rep), " MLGs without representative are left out: ", paste(no_rep, collapse = ", "))
  iso <- filter(iso, mlg %in% rep_mlg)
}

## check if the isolates of the depth project have direction and depth
if (anyNA(filter(iso, project == depth_proj)[c("direction", "depth")])) {
  stop("isolates of ", depth_proj, " without direction or depth")
}



# analysis per unit ------------------------------------------------------------

### depth class per MLG (the interval it was found in, several intervals, or MAGIC)
depth_of_mlg <- function(d) {
  d %>%
    summarise(proj = if (all(project != depth_proj)) "magic" else depth_proj,
              depths = list(unique(depth[project == depth_proj])), .by = mlg) %>%
    mutate(depth_class = case_when(proj == "magic" ~ "MAGIC (0-10 cm)",
                                   lengths(depths) > 1 ~ "Several depths",
                                   .default = lab_depth[map_chr(depths, ~ .x[1])])) %>%
    select(mlg, depth_class)
}

### depth labels permuted among the cores of a site and direction (as in 031)
permute_depth <- function(cores) {
  strat <- paste(cores$site, cores$direction)
  lab   <- cores$depth
  for (s in unique(strat)) {
    j <- which(strat == s)
    lab[j] <- lab[j][sample.int(length(j))]
  }
  setNames(lab, cores$core)
}

### Cramer's V between cluster and depth over the MLG occurrences per core
cramer_v <- function(cl, dep) {
  tab <- table(cl, dep)
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  if (min(dim(tab)) < 2) return(NA_real_)
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  unname(sqrt(chi / (sum(tab) * (min(dim(tab)) - 1))))
}

### slope of a principal coordinate (in standard deviations) on depth (intervals 1 to 3) within sites, over the
### MLG occurrences per core
slope_depth <- function(score, dep, site) {
  X <- cbind(model.matrix(~ factor(site) - 1), depth = as.numeric(dep))
  unname(lm.fit(X, as.numeric(scale(score)))$coefficients["depth"])
}

substructure <- function(uu) {
  d    <- filter(iso, unit == uu)
  mlgs <- sort(unique(d$mlg))
  D    <- p_dist[mlgs, mlgs]

  ### principal coordinates (share of the positive eigenvalues)
  pc   <- cmdscale(as.dist(D), k = n_axes, eig = TRUE)
  prop <- pc$eig[seq_len(n_axes)] / sum(pc$eig[pc$eig > 0])

  ### clusters for every number in k_range (the smallest number within sil_tol of the largest width chosen)
  sel  <- tibble(k = k_range, silhouette = map_dbl(k_range, ~ pam(as.dist(D), .x, diss = TRUE)$silinfo$avg.width))
  k_bt <- min(sel$k[sel$silhouette >= max(sel$silhouette) - sil_tol])
  cl   <- setNames(pam(as.dist(D), k_bt, diss = TRUE)$clustering, mlgs)

  mlg_tbl <- tibble(mlg = mlgs, cluster = factor(cl[mlgs])) %>%
    bind_cols(setNames(as.data.frame(pc$points), paste0("PCo", seq_len(n_axes)))) %>%
    left_join(depth_of_mlg(d), by = "mlg")

  ### MLG occurrences per core at the depth project, with cluster and principal coordinates
  occ <- d %>%
    filter(project == depth_proj) %>%
    distinct(site, direction, depth, mlg) %>%
    mutate(core = paste(site, direction, depth, sep = "|")) %>%
    left_join(mlg_tbl, by = "mlg")
  cores <- distinct(occ, core, site, direction, depth)

  ### observed association with depth and the permutation null (depth labels of the cores within site and direction)
  stats <- function(dep) {
    c(cramer_v = cramer_v(occ$cluster, dep),
      setNames(map_dbl(paste0("PCo", seq_len(n_axes)), ~ slope_depth(occ[[.x]], dep, occ$site)),
               paste0("slope_PCo", seq_len(n_axes))))
  }
  obs  <- stats(occ$depth)
  null <- replicate(n_perm, {
    pd <- permute_depth(cores)
    stats(pd[occ$core])
  })
  assoc <- tibble(unit = uu, test = names(obs), observed = obs,
                  null_mean = rowMeans(null, na.rm = TRUE),
                  null_lo = apply(null, 1, quantile, 0.025, na.rm = TRUE),
                  null_hi = apply(null, 1, quantile, 0.975, na.rm = TRUE),
                  p = map_dbl(seq_along(obs), function(i) {
                    x <- null[i, ]
                    if (names(obs)[i] == "cramer_v") (1 + sum(x >= obs[i], na.rm = TRUE)) / (1 + sum(!is.na(x)))
                    else (1 + sum(abs(x - mean(x, na.rm = TRUE)) >= abs(obs[i] - mean(x, na.rm = TRUE)),
                                  na.rm = TRUE)) / (1 + sum(!is.na(x)))
                  }))

  ### clusters: MLGs, MAGIC MLGs, occurrences per depth and sites; sites holding several clusters
  within_d <- tibble(cluster = factor(sort(unique(cl))),
                     mean_within = map_dbl(sort(unique(cl)), function(cc) {
                       m <- names(cl)[cl == cc]
                       if (length(m) < 2) NA_real_ else mean(D[m, m][upper.tri(D[m, m])])
                     }))
  clusters <- mlg_tbl %>%
    summarise(mlgs = n(), mlgs_magic = sum(depth_class == "MAGIC (0-10 cm)"), .by = cluster) %>%
    left_join(within_d, by = "cluster") %>%
    left_join(occ %>% count(cluster, depth) %>% mutate(depth = lab_depth[depth]) %>%
                pivot_wider(names_from = depth, values_from = n, values_fill = 0), by = "cluster") %>%
    left_join(distinct(d, mlg, site) %>% left_join(select(mlg_tbl, mlg, cluster), by = "mlg") %>%
                summarise(sites = n_distinct(site), .by = cluster), by = "cluster") %>%
    arrange(cluster)
  co_sites <- distinct(d, mlg, site) %>%
    mutate(cluster = cl[mlg]) %>%
    summarise(n_clusters = n_distinct(cluster), .by = site)

  list(unit = uu, D = D, mlg_tbl = mlg_tbl, prop = prop, sel = mutate(sel, unit = uu, .before = 1), k_best = k_bt,
       occ = occ, assoc = assoc, clusters = mutate(clusters, unit = uu, .before = 1),
       co_sites = tibble(unit = uu, sites = nrow(co_sites), sites_several = sum(co_sites$n_clusters > 1)))
}

set.seed(seed)
res <- map(units_focal, substructure)

selection <- map_dfr(res, "sel")
assoc     <- map_dfr(res, "assoc")
clusters  <- map_dfr(res, "clusters")
co_sites  <- map_dfr(res, "co_sites")



# figures ----------------------------------------------------------------------

### per unit: principal coordinates of the MLGs by depth class with clusters as convex hulls, average silhouette
### width per number of clusters (chosen number filled), depth of the occurrences per cluster (PhytOakmeter)
plot_unit <- function(r) {
  m <- r$mlg_tbl %>% mutate(depth_class = factor(depth_class, levels = names(pal_depth)))
  hulls <- m %>%
    group_by(cluster) %>%
    filter(n() >= 3) %>%
    slice(chull(PCo1, PCo2)) %>%
    ungroup()
  centres <- summarise(m, PCo1 = mean(PCo1), PCo2 = mean(PCo2), .by = cluster)
  p_a <- ggplot(m, aes(PCo1, PCo2)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_polygon(data = hulls, aes(group = cluster), fill = NA, colour = "grey60", linewidth = 0.3) +
    geom_point(aes(fill = depth_class), shape = 21, size = 2, stroke = 0.3, colour = "black") +
    geom_text(data = centres, aes(label = cluster), size = 7 / .pt, fontface = "bold", colour = "black") +
    scale_fill_manual(values = pal_depth, drop = FALSE, name = NULL) +
    guides(fill = guide_legend(nrow = 2)) +
    scale_x_continuous(labels = scales::label_number(drop0trailing = TRUE)) +
    scale_y_continuous(labels = scales::label_number(drop0trailing = TRUE)) +
    labs(x = sprintf("PCo1 (%.1f%%)", 100 * r$prop[1]), y = sprintf("PCo2 (%.1f%%)", 100 * r$prop[2]))
  p_b <- ggplot(r$sel, aes(k, silhouette)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = max(r$sel$silhouette) - sil_tol, ymax = max(r$sel$silhouette),
             fill = "grey90") +
    geom_hline(yintercept = sil_min, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
    geom_line(linewidth = 0.4, colour = "grey30") +
    geom_point(aes(fill = k == r$k_best), shape = 21, size = 2, stroke = 0.3, colour = "black") +
    scale_fill_manual(values = c(`TRUE` = "grey15", `FALSE` = "white"), guide = "none") +
    scale_x_continuous(breaks = k_range) +
    labs(x = "Clusters", y = "Average silhouette width")
  p_c <- r$occ %>%
    count(cluster, depth) %>%
    mutate(depth = factor(lab_depth[depth], levels = lab_depth), share = n / sum(n), .by = cluster) %>%
    ggplot(aes(cluster, share, fill = depth)) +
    geom_col(width = 0.7, colour = "grey20", linewidth = 0.2) +
    geom_text(data = ~ summarise(.x, n = sum(n), .by = cluster), aes(cluster, 1.05, label = n), size = 6 / .pt,
              inherit.aes = FALSE) +
    scale_fill_manual(values = pal_depth[lab_depth], guide = "none") +
    scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Cluster", y = "Share of occurrences")
  (p_a | p_b | p_c) + plot_layout(widths = c(1.6, 1, 1), guides = "collect") + plot_annotation(tag_levels = "A") &
    theme(legend.position = "bottom", legend.box = "horizontal")
}

walk(res, ~ save_fig(sprintf("substructure_%s.pdf", .x$unit), plot_unit(.x), 7.1, 3.2))


### per unit: neighbour-joining tree of the MLGs (unrooted), tips by depth class, cluster numbers at the tips
walk(res, function(r) {
  tr <- nj(as.dist(r$D))
  m  <- r$mlg_tbl[match(tr$tip.label, r$mlg_tbl$mlg), ]
  cairo_pdf(file.path(out_root, sprintf("tree_%s.pdf", r$unit)), width = 5, height = 5)
  par(mar = c(1, 1, 1, 1))
  plot(tr, type = "unrooted", show.tip.label = FALSE, edge.width = 0.6, edge.color = "grey40")
  tiplabels(pch = 21, bg = pal_depth[m$depth_class], col = "black", cex = 0.9)
  tiplabels(as.character(m$cluster), frame = "none", cex = 0.45, adj = c(-0.8, 0.5))
  add.scale.bar(cex = 0.6)
  legend("bottomright", legend = names(pal_depth), pt.bg = pal_depth, pch = 21, cex = 0.6, bty = "n")
  dev.off()
})



# tables -----------------------------------------------------------------------

### test labels
lab_test <- c(cramer_v = "Cluster and depth (Cramer's V)",
              setNames(sprintf("PCo%d on depth (SD per interval)", seq_len(n_axes)),
                       paste0("slope_PCo", seq_len(n_axes))))

### average silhouette width per unit and number of clusters
pub_selection <- selection %>%
  left_join(tibble(unit = map_chr(res, "unit"), k_best = map_int(res, "k_best")), by = "unit") %>%
  transmute(Unit = unit, Clusters = k, `Average silhouette width` = sprintf("%.3f", silhouette),
            Chosen = if_else(k == k_best,
                             if_else(silhouette >= sil_min, "yes", "yes (below 0.25: no substantial structure)"), ""))

### association with depth per unit and test
pub_assoc <- assoc %>%
  transmute(Unit = unit, Test = lab_test[test], Observed = sprintf("%.3f", observed),
            `Null mean [95% range]` = sprintf("%.3f [%.3f, %.3f]", null_mean, null_lo, null_hi), p = fmt_p(p))

### clusters per unit
pub_clusters <- clusters %>%
  left_join(co_sites, by = "unit", suffix = c("", "_unit")) %>%
  transmute(Unit = unit, Cluster = cluster, MLGs = mlgs, `MAGIC MLGs` = mlgs_magic, Sites = sites,
            `Mean p-distance within` = sprintf("%.4f", mean_within),
            !!!syms(intersect(unname(lab_depth), names(clusters))),
            `Sites with several clusters` = sprintf("%d of %d", sites_several, sites_unit))

write_tsv(map_dfr(res, ~ mutate(.x$mlg_tbl, unit = .x$unit, .before = 1)),
          file.path(out_root, "substructure_mlgs.tsv"))
write_tsv(selection, file.path(out_root, "substructure_selection.tsv"))
write_tsv(assoc, file.path(out_root, "substructure_depth_tests.tsv"))
write_tsv(pub_selection, file.path(out_root, "publication_substructure_selection.tsv"))
write_tsv(pub_assoc, file.path(out_root, "publication_substructure_depth_tests.tsv"))
write_tsv(pub_clusters, file.path(out_root, "publication_substructure_clusters.tsv"))



# summary ----------------------------------------------------------------------

cat("\n=== Number of clusters ===\n")
print(pub_selection, n = Inf, width = Inf)

cat("\n=== Clusters ===\n")
print(pub_clusters, n = Inf, width = Inf)

cat("\n=== Association with depth ===\n")
print(pub_assoc, n = Inf, width = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  substructure_<unit>.pdf\n")
cat("  tree_<unit>.pdf\n")
cat("  substructure_mlgs.tsv\n")
cat("  substructure_selection.tsv\n")
cat("  substructure_depth_tests.tsv\n")
cat("  publication_substructure_selection.tsv\n")
cat("  publication_substructure_depth_tests.tsv\n")
cat("  publication_substructure_clusters.tsv\n")
