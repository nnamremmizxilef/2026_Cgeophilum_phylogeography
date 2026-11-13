## ============================================================================
## 031_within_site_relatedness.R
## Relatedness of distinct genets within the units at K = 4, 5 and 6 inside the PhytOakmeter
## sites: sequence distance between distinct MLGs of a unit within a site against their depth
## and direction, without and with data quality as covariates, tested by permuting the labels
## of the soil cores; per unit and pooled, all isolates and without mixed-MAT isolates
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
##   sample_metadata_ALL_FIXED.tsv  project, site, direction, sampling_depth, mlg, mat_genome,
##                                  mean_seq_depth and pct_scos_present per sample,
##                                  from local/01_variant_calling/103_mating_type.R
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

data_root   <- "04_diversity_niche_selection/data/010_data_sequence_distances"
labels_file <- "03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv"
meta_file   <- "01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv"
out_root    <- "04_diversity_niche_selection/results/031_results_within_site_relatedness"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### levels, level of the figure, sample sets, project sampled at several depths and directions
k_levels <- c("k4", "k5", "k6")
k_main   <- "k6"
sets     <- c("all", "no_mixed")
project  <- "phytoakmeter"

### mating type calls of possible mixed cultures (sample set without them)
mat_exclude <- c("MAT12", "ambiguous")

### smallest number of pairs of distinct MLGs of a unit for its own test
n_min_pairs <- 30

### permutations, bootstrap replicates over sites (figure intervals), seed
n_perm <- 999
n_boot <- 999
seed   <- 42

### depth intervals and depth combinations (against the outermost intervals 5-10 and 45-50 cm)
lab_depth <- c(`1` = "5-10 cm", `2` = "25-30 cm", `3` = "45-50 cm")
lab_class <- c(same_1 = "Both 5-10 cm", same_2 = "Both 25-30 cm", same_3 = "Both 45-50 cm",
               apart_1 = "Adjacent intervals")

### separation classes in depth and around the tree
lab_vsep  <- c(`0` = "Same depth", `1` = "20 cm apart", `2` = "40 cm apart")
lab_hsep  <- c(`2` = "Same pit", `1` = "Adjacent pits (7 m)", `0` = "Opposite pits (10 m)")

### sampling design: depth interval midpoints (cm), pits 5 m from the tree in four directions (7.07 or 10 m apart)
depth_mid <- c(`1` = 7.5, `2` = 27.5, `3` = 47.5)
dir_pos   <- c(N = 0, E = 1, S = 2, W = 3)
pit_dist  <- outer(dir_pos, dir_pos, function(a, b) c(0, 5 * sqrt(2), 10, 5 * sqrt(2))[abs(a - b) + 1])

### separations at which the trends are expressed (cm, m)
v_max     <- 40
h_max     <- 10

### effects: separation trends, same depth and direction (without and with data quality), depth combinations
lab_term  <- c(v_trend = "Depth separation (effect at 40 cm)", h_trend = "Horizontal separation (effect at 10 m)",
               depth = "Same depth", direction = "Same direction",
               depth_adj = "Same depth", direction_adj = "Same direction", lab_class)

### effects tested with permuted depth labels, with permuted direction labels; trends (p from the upper tail)
dep_terms <- c("v_trend", "depth", "depth_adj", names(lab_class))
dir_terms <- c("h_trend", "direction", "direction_adj")
pos_terms <- c("v_trend", "h_trend")



# style ------------------------------------------------------------------------

### unit colours (also the units tested and their order; all: units pooled)
pal_unit  <- c(L1  = "#E69F00",
               L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
               L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
               L4  = "#CC79A7", all = "grey60")

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}

### p values and level labels
fmt_p <- function(p) if_else(p < 0.001, "< 0.001", sprintf("%.3f", p))
lab_k <- function(k) sub("k", "K = ", k)



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

### p-distance and compared sites between MLGs (via their representatives), sequencing depth per representative
rep_sample <- rownames(diffs)
rep_mlg    <- meta$mlg[match(rep_sample, meta$sample)]
p_dist     <- diffs / cmp_sites
diag(p_dist) <- 0
dimnames(p_dist) <- dimnames(cmp_sites) <- list(rep_mlg, rep_mlg)
rep_seq <- setNames(as.numeric(meta$mean_seq_depth[match(rep_sample, meta$sample)]), rep_mlg)

### isolates of the project with site, direction, depth, MLG, units and data quality
iso <- lab %>%
  select(sample, all_of(k_levels)) %>%
  inner_join(meta %>% filter(project == !!project) %>%
               select(sample, site, direction, depth = sampling_depth, mlg, mat_genome, mean_seq_depth,
                      pct_scos_present),
             by = "sample") %>%
  mutate(mixed = mat_genome %in% mat_exclude, across(c(mean_seq_depth, pct_scos_present), as.numeric))



# checks -----------------------------------------------------------------------

## check if the matrices hold the representatives of classes.tsv in the same order
if (!identical(rep_sample, classes$sample) || !identical(colnames(diffs), classes$sample)) {
  stop("distance matrices and classes.tsv differ in samples or order")
}

## check if every pair has compared sites
if (any(cmp_sites[upper.tri(cmp_sites)] == 0)) {
  stop("pairs of representatives without compared sites")
}

## check if every representative has one MLG, a sequencing depth, and every MLG one representative
if (anyNA(rep_mlg) || anyDuplicated(rep_mlg) || anyNA(rep_seq)) {
  stop("representatives without MLG or sequencing depth, or MLGs with several representatives")
}

## check if every isolate has a direction, a depth and an MLG
if (anyNA(iso$direction) || anyNA(iso$depth) || anyNA(iso$mlg)) {
  stop("isolates without direction, depth or MLG")
}

## check if every MLG has one unit per level
if (any(map_lgl(k_levels, ~ any(count(distinct(iso, mlg, .data[[.x]]), mlg)$n > 1)))) {
  stop("MLGs with several units at a level")
}

## check if every MLG has a representative in the distance matrices
no_rep <- setdiff(unique(iso$mlg), rep_mlg)

if (length(no_rep) > 0) {
  warning(length(no_rep), " MLGs without representative are left out: ", paste(no_rep, collapse = ", "))
  iso <- filter(iso, mlg %in% rep_mlg)
}



# data -------------------------------------------------------------------------

### MLG occurrences per soil core (site, direction, depth) for a sample set
occurrences <- function(ss) {
  d <- if (ss == "all") iso else filter(iso, !mixed)
  distinct(d, site, direction, depth, mlg, !!!syms(k_levels)) %>%
    mutate(core = paste(site, direction, depth, sep = "|"))
}

### pairs of distinct MLGs within a site (clonemates excluded): p-distance, compared sites, lower sequencing depth
pairs_of <- function(occ) {
  occ <- mutate(occ, i = row_number())
  pr  <- inner_join(occ, occ, by = "site", suffix = c("_1", "_2"), relationship = "many-to-many") %>%
    filter(i_1 < i_2, mlg_1 != mlg_2)
  idx <- cbind(match(pr$mlg_1, rownames(p_dist)), match(pr$mlg_2, rownames(p_dist)))
  mutate(pr, dist = p_dist[idx], compared = cmp_sites[idx], seq_min = pmin(rep_seq[pr$mlg_1], rep_seq[pr$mlg_2]))
}

### depth combination of a pair (same interval, adjacent intervals 20 cm apart, outermost intervals 40 cm apart)
depth_class <- function(a, b) {
  case_when(a == b ~ paste0("same_", a),
            abs(as.integer(a) - as.integer(b)) == 1 ~ "apart_1",
            .default = "apart_2")
}

### core labels permuted within strata (depth within site and direction, direction within site and depth)
permute_cores <- function(cores, dim) {
  strat <- if (dim == "depth") paste(cores$site, cores$direction) else paste(cores$site, cores$depth)
  lab   <- cores[[dim]]
  for (s in unique(strat)) {
    j <- which(strat == s)
    lab[j] <- lab[j][sample.int(length(j))]
  }
  setNames(lab, cores$core)
}

### test groups per level (every unit with enough pairs, all units pooled) with strata and standardised data quality
groups_of <- function(pr, kk) {
  pr <- filter(pr, .data[[paste0(kk, "_1")]] == .data[[paste0(kk, "_2")]])
  pr$unit <- pr[[paste0(kk, "_1")]]
  units <- intersect(names(pal_unit), names(which(table(pr$unit) >= n_min_pairs)))
  c(setNames(map(units, ~ which(pr$unit == .x)), units), list(all = seq_len(nrow(pr)))) %>%
    map(function(i) {
      g <- pr[i, ]
      list(pr = g, S = model.matrix(~ factor(paste(g$site, g$unit)) - 1),
           Q = scale(cbind(compared = g$compared, seq_min = g$seq_min)))
    })
}

### columns with variation (constant test and covariate columns of a group dropped)
varying <- function(M) M[, apply(M, 2, function(x) all(is.finite(x)) && var(x) > 0), drop = FALSE]

### effects on the p-distance of distinct MLGs (% of the mean distance; site, and unit when pooled, as strata)
fit_terms <- function(g, dep, dirn) {
  pr <- g$pr
  d1 <- dep[pr$core_1]
  d2 <- dep[pr$core_2]
  s_dep <- as.numeric(d1 == d2)
  s_dir <- as.numeric(dirn[pr$core_1] == dirn[pr$core_2])
  cls   <- depth_class(d1, d2)
  C     <- sapply(names(lab_class), function(cl) as.numeric(cls == cl))
  coef_of <- function(X, terms) {
    b <- lm.fit(cbind(g$S, varying(X)), pr$dist)$coefficients
    setNames(100 * unname(b[terms]) / mean(pr$dist), terms)
  }
  v_sep <- abs(depth_mid[d1] - depth_mid[d2])
  h_sep <- pit_dist[cbind(dirn[pr$core_1], dirn[pr$core_2])]
  tr    <- coef_of(cbind(v_trend = v_sep, h_trend = h_sep, g$Q), c("v_trend", "h_trend")) * c(v_max, h_max)
  c(tr,
    coef_of(cbind(depth = s_dep, direction = s_dir), c("depth", "direction")),
    setNames(coef_of(cbind(depth = s_dep, direction = s_dir, g$Q), c("depth", "direction")),
             c("depth_adj", "direction_adj")),
    coef_of(cbind(C, direction = s_dir, g$Q), names(lab_class)))
}



# tests ------------------------------------------------------------------------

### observed effects and permutation nulls per sample set, level, group and effect (one-sided p)
run_set <- function(ss) {
  occ   <- occurrences(ss)
  pr    <- pairs_of(occ)
  cores <- distinct(occ, core, site, direction, depth)
  obs_d <- setNames(cores$depth, cores$core)
  obs_r <- setNames(cores$direction, cores$core)

  grp <- map_dfr(k_levels, function(kk) {
    g <- groups_of(pr, kk)
    tibble(k = kk, unit = names(g), g = unname(g))
  })
  obs <- map(grp$g, ~ fit_terms(.x, obs_d, obs_r))

  null_dep <- replicate(n_perm, {
    pd <- permute_cores(cores, "depth")
    unlist(map(grp$g, ~ fit_terms(.x, pd, obs_r)[dep_terms]))
  })
  null_dir <- replicate(n_perm, {
    pr_dir <- permute_cores(cores, "direction")
    unlist(map(grp$g, ~ fit_terms(.x, obs_d, pr_dir)[dir_terms]))
  })

  grp %>%
    mutate(i = row_number(), pairs = map_int(g, ~ nrow(.x$pr)), sites = map_int(g, ~ n_distinct(.x$pr$site)),
           mean_dist = map_dbl(g, ~ mean(.x$pr$dist))) %>%
    select(-g) %>%
    expand_grid(term = names(lab_term)) %>%
    mutate(observed = map2_dbl(i, term, ~ obs[[.x]][[.y]]),
           null = map2(i, term, function(i, t) {
             if (t %in% dep_terms) null_dep[(i - 1) * length(dep_terms) + match(t, dep_terms), ]
             else null_dir[(i - 1) * length(dir_terms) + match(t, dir_terms), ]
           }),
           null_mean = map_dbl(null, mean, na.rm = TRUE),
           null_lo   = map_dbl(null, ~ quantile(.x, 0.025, na.rm = TRUE)),
           null_hi   = map_dbl(null, ~ quantile(.x, 0.975, na.rm = TRUE)),
           p = pmap_dbl(list(null, observed, term), function(x, o, t) {
             if (is.na(o)) return(NA_real_)
             hits <- if (t %in% pos_terms) x >= o else x <= o
             (1 + sum(hits, na.rm = TRUE)) / (1 + sum(!is.na(x)))
           }),
           set = ss, .before = 1) %>%
    select(-null, -i)
}

set.seed(seed)
res <- map_dfr(sets, run_set)



# data quality by depth --------------------------------------------------------

### data quality per isolate by sampling depth, compared sites per pair by depth combination (main level,
### all units pooled, all isolates)
qual_iso <- iso %>%
  select(depth, `Mean sequencing depth` = mean_seq_depth, `Single-copy orthologues present (%)` = pct_scos_present) %>%
  pivot_longer(-depth, names_to = "measure")

pairs_main <- pairs_of(occurrences("all")) %>%
  filter(.data[[paste0(k_main, "_1")]] == .data[[paste0(k_main, "_2")]]) %>%
  mutate(class = depth_class(depth_1, depth_2))

quality <- bind_rows(
  qual_iso %>%
    summarise(n = n(), median = median(value), q1 = quantile(value, 0.25), q3 = quantile(value, 0.75),
              .by = c(measure, depth)) %>%
    rename(group = depth) %>%
    mutate(group = lab_depth[group]),
  pairs_main %>%
    summarise(n = n(), median = median(compared), q1 = quantile(compared, 0.25), q3 = quantile(compared, 0.75),
              .by = class) %>%
    mutate(measure = "Compared sites per pair", group = coalesce(lab_class[class], "Outermost intervals")) %>%
    select(-class))



# figures ----------------------------------------------------------------------

### distance between distinct MLGs of a unit at a tree by separation in depth and around the tree (main level,
### all isolates), as % of the mean of the site (and unit when pooled)
sep_pairs <- function(pr) {
  pr %>%
    mutate(rel   = 100 * dist / ave(dist, paste(site, unit)),
           v_sep = factor(lab_vsep[as.character(abs(as.integer(depth_1) - as.integer(depth_2)))], levels = lab_vsep),
           h_sep = factor(lab_hsep[as.character(abs(dir_pos[direction_1] - dir_pos[direction_2]) %% 2 +
                                                  2 * (direction_1 == direction_2))], levels = lab_hsep))
}

### mean with a 95% interval from bootstrapping sites
boot_ci <- function(d, cat) {
  lv <- levels(d[[cat]])
  by_site <- split(d, d$site)
  bs <- replicate(n_boot, {
    b <- bind_rows(by_site[sample(names(by_site), replace = TRUE)])
    tapply(b$rel, factor(b[[cat]], levels = lv), mean)
  })
  tibble(level = factor(lv, levels = lv), mean = as.numeric(tapply(d$rel, factor(d[[cat]], levels = lv), mean)),
         lo = apply(bs, 1, quantile, 0.025, na.rm = TRUE), hi = apply(bs, 1, quantile, 0.975, na.rm = TRUE))
}

### one panel per unit tested at the main level and one for all units pooled
pr_main <- pairs_of(occurrences("all")) %>%
  filter(.data[[paste0(k_main, "_1")]] == .data[[paste0(k_main, "_2")]]) %>%
  mutate(unit = .data[[paste0(k_main, "_1")]])
panels  <- res %>% filter(set == "all", k == k_main) %>% distinct(unit) %>% pull(unit)
fig_pairs <- bind_rows(map(setdiff(panels, "all"), ~ filter(pr_main, unit == .x) %>% mutate(panel = .x)),
                       mutate(pr_main, panel = "all")) %>%
  sep_pairs() %>%
  mutate(panel = factor(panel, levels = intersect(names(pal_unit), panels)))

set.seed(seed)
fig_sep <- fig_pairs %>%
  group_by(panel) %>%
  group_modify(~ bind_rows(mutate(boot_ci(.x, "v_sep"), axis = "depth"),
                           mutate(boot_ci(.x, "h_sep"), axis = "direction"))) %>%
  ungroup()

### p of the trend with the separation (permutation test) in the panel corner
fig_p <- res %>%
  filter(set == "all", k == k_main, term %in% pos_terms) %>%
  transmute(panel = factor(unit, levels = levels(fig_sep$panel)),
            axis  = if_else(term == "v_trend", "depth", "direction"),
            label = paste0("p ", if_else(p < 0.001, "< 0.001", sprintf("= %.3f", p))))

plot_sep <- function(ax, x_lab) {
  d <- filter(fig_sep, axis == ax)
  ggplot(d, aes(level, mean, group = 1)) +
    geom_hline(yintercept = 100, linetype = "dashed", linewidth = 0.3, colour = "grey60") +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.4, colour = "grey40") +
    geom_line(linewidth = 0.4, colour = "grey40") +
    geom_point(aes(fill = panel), shape = 21, size = 2.4, stroke = 0.3, colour = "black") +
    geom_text(data = filter(fig_p, axis == ax), aes(x = -Inf, y = Inf, label = label), hjust = -0.1, vjust = 1.4,
              size = 6.5 / .pt, inherit.aes = FALSE) +
    facet_wrap(~ panel, nrow = 1, labeller = as_labeller(function(x) if_else(x == "all", "All", x))) +
    scale_fill_manual(values = pal_unit, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    labs(x = x_lab, y = "Genetic distance (% of site mean)") +
    theme(strip.background = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
}

p_rel <- (plot_sep("depth", "Depth separation") / plot_sep("direction", "Horizontal separation")) +
  plot_annotation(tag_levels = "A")

save_fig("within_site_relatedness.pdf", p_rel, 7.1, 5.6)


### data quality per isolate by sampling depth, compared sites per pair by depth combination
p_qi <- qual_iso %>%
  mutate(depth = factor(lab_depth[depth], levels = lab_depth)) %>%
  ggplot(aes(depth, value)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.3, fill = "grey90") +
  geom_jitter(width = 0.15, shape = 21, colour = "black", fill = "grey55", size = 0.5, stroke = 0.2,
              alpha = 0.5) +
  facet_wrap(~ measure, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = NULL) +
  theme(strip.background = element_blank())

p_qp <- pairs_main %>%
  mutate(class = factor(coalesce(lab_class[class], "Outermost intervals"),
                        levels = c(lab_class, "Outermost intervals"))) %>%
  ggplot(aes(class, compared)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.3, fill = "grey90") +
  geom_jitter(width = 0.15, shape = 21, colour = "black", fill = "grey55", size = 0.4, stroke = 0.2,
              alpha = 0.3) +
  labs(x = NULL, y = "Compared sites per pair") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_quality <- (p_qi | p_qp) + plot_layout(widths = c(2, 1.3)) + plot_annotation(tag_levels = "A")

save_fig("data_quality_by_depth.pdf", p_quality, 7.1, 3)



# tables -----------------------------------------------------------------------

### one row per level, group and effect (all isolates), the p without mixed-MAT isolates alongside
pub_res <- res %>%
  filter(set == "all") %>%
  left_join(res %>% filter(set == "no_mixed") %>% select(k, unit, term, p_nm = p), by = c("k", "unit", "term")) %>%
  mutate(model = case_when(term %in% pos_terms ~ "Trend with separation, with data quality",
                           term %in% c("depth", "direction") ~ "Without data quality",
                           term %in% c("depth_adj", "direction_adj") ~ "With data quality",
                           .default = "Depth combination (against the outermost intervals), with data quality")) %>%
  arrange(k, match(unit, names(pal_unit)), match(term, names(lab_term))) %>%
  transmute(Level = lab_k(k), Unit = if_else(unit == "all", "All", unit), Model = model, Effect = lab_term[term],
            Pairs = pairs, Sites = sites, `Mean p-distance` = sprintf("%.4f", mean_dist),
            `Effect (% of mean)`    = sprintf("%.1f", observed),
            `Null mean [95% range]` = sprintf("%.1f [%.1f, %.1f]", null_mean, null_lo, null_hi),
            p                       = if_else(is.na(p), "-", fmt_p(p)),
            `Without mixed MAT: p`  = if_else(is.na(p_nm), "-", fmt_p(p_nm)))

### median and interquartile range per measure and group
pub_quality <- quality %>%
  transmute(Measure = measure, Group = group, N = n, `Median [IQR]` = sprintf("%.1f [%.1f, %.1f]", median, q1, q3))

write_tsv(res, file.path(out_root, "within_site_relatedness.tsv"))
write_tsv(fig_sep, file.path(out_root, "within_site_relatedness_separation.tsv"))
write_tsv(quality, file.path(out_root, "data_quality_by_depth.tsv"))
write_tsv(pub_res, file.path(out_root, "publication_within_site_relatedness.tsv"))
write_tsv(pub_quality, file.path(out_root, "publication_data_quality_by_depth.tsv"))



# summary ----------------------------------------------------------------------

cat("\n=== Data quality by depth ===\n")
print(pub_quality, n = Inf, width = Inf)

cat("\n=== Relatedness of distinct genets within sites ===\n")
print(pub_res, n = Inf, width = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  within_site_relatedness.pdf\n")
cat("  data_quality_by_depth.pdf\n")
cat("  within_site_relatedness.tsv\n")
cat("  within_site_relatedness_separation.tsv\n")
cat("  data_quality_by_depth.tsv\n")
cat("  publication_within_site_relatedness.tsv\n")
cat("  publication_data_quality_by_depth.tsv\n")
