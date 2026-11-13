## ============================================================================
## 011_nucleotide_diversity.R
## Nucleotide diversity within units (pi), absolute and net divergence between
## units (d_xy, d_a) and Hudson's F_ST at K = 4, 5 and 6, from all sites of the
## single-copy orthologues, with 95% confidence intervals from bootstrapping
## loci; all representatives and without mixed-MAT samples
##
## Files needed:
##   classes.tsv              representatives, units, mixed-MAT flag, class
##   class_pair_sums.tsv.gz   per locus and class pair: differences and compared sites
##   loci.tsv                 per locus: length, called share, variable sites
##                            (all three from hpc/04_diversity_niche_selection/010_sequence_distances.sh;
##                             copy them from ${ROOT}/04_diversity_niche_selection/results/010_results_sequence_distances
##                             to local/04_diversity_niche_selection/data/010_data_sequence_distances)
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

data_root <- "04_diversity_niche_selection/data/010_data_sequence_distances"
out_root  <- "04_diversity_niche_selection/results/011_results_nucleotide_diversity"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### levels
k_levels <- c("k4", "k5", "k6")

### bootstrap replicates over loci and seed
n_boot <- 1000
seed   <- 42



# style ------------------------------------------------------------------------

### unit colours
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

### F_ST colour scale
pal_fst <- c("white", "#B2182B")

### sample set labels and shapes
lab_set   <- c(all = "All representatives", no_mixed = "Without mixed-MAT samples")
shape_set <- c(all = 21, no_mixed = 24)

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}

### estimate [95% CI]
fmt_ci <- function(m, lo, hi, digits = 4) sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
                                                  m, lo, hi)



# read -------------------------------------------------------------------------

classes <- read_tsv(file.path(data_root, "classes.tsv"), col_types = cols(.default = "c", mixed = "l"))
sums    <- read_tsv(file.path(data_root, "class_pair_sums.tsv.gz"), col_types = "cccdd")
loci    <- read_tsv(file.path(data_root, "loci.tsv"), col_types = "cddd")

cls <- distinct(classes, class, across(all_of(k_levels)), mixed)

### classes per sample set
sets <- list(all = cls$class, no_mixed = cls$class[!cls$mixed])

### representatives per level, unit and sample set
n_rep <- expand_grid(k = k_levels, set = names(sets)) %>%
  pmap_dfr(function(k, set) classes %>% filter(class %in% sets[[set]]) %>% count(unit = .data[[k]]) %>%
             mutate(k = k, set = set, .before = 1))



# per locus and unit pair ------------------------------------------------------

### class pairs summed into unit pairs (pairs of classes in one unit count as within the unit)
agg <- expand_grid(k = k_levels, set = names(sets)) %>%
  mutate(data = map2(k, set, function(kk, ss) {
    key <- setNames(cls[[kk]], cls$class)
    sums %>%
      filter(class_a %in% sets[[ss]], class_b %in% sets[[ss]]) %>%
      mutate(unit_a = pmin(key[class_a], key[class_b]), unit_b = pmax(key[class_a], key[class_b])) %>%
      summarise(differences = sum(differences), sites = sum(sites), .by = c(locus, unit_a, unit_b))
  })) %>%
  unnest(data)



# estimates and bootstrap ------------------------------------------------------

### pi, d_xy, d_a and Hudson's F_ST from summed differences over summed compared sites
derive <- function(est) {
  pi_u <- est %>% filter(unit_a == unit_b) %>% select(unit = unit_a, value)
  est %>%
    filter(unit_a != unit_b) %>%
    left_join(rename(pi_u, pi_a = value), by = c(unit_a = "unit")) %>%
    left_join(rename(pi_u, pi_b = value), by = c(unit_b = "unit")) %>%
    transmute(unit_a, unit_b, dxy = value, da = value - (pi_a + pi_b) / 2,
              fst = 1 - ((pi_a + pi_b) / 2) / value) %>%
    list(pi = pi_u, pairs = .)
}

set.seed(seed)

res <- agg %>%
  nest(.by = c(k, set)) %>%
  mutate(out = map(data, function(d) {
    pairs <- distinct(d, unit_a, unit_b)
    D <- d %>% select(locus, unit_a, unit_b, differences) %>%
      pivot_wider(names_from = c(unit_a, unit_b), values_from = differences, values_fill = 0)
    S <- d %>% select(locus, unit_a, unit_b, sites) %>%
      pivot_wider(names_from = c(unit_a, unit_b), values_from = sites, values_fill = 0)
    Dm <- as.matrix(D[, -1])
    Sm <- as.matrix(S[, colnames(D)[-1]])
    key <- tibble(col = colnames(Dm)) %>%
      left_join(mutate(pairs, col = paste(unit_a, unit_b, sep = "_")), by = "col")

    ### one estimate from given locus weights (1 = every locus once)
    est <- function(w) derive(mutate(key, value = as.numeric(crossprod(w, Dm) / crossprod(w, Sm))))

    point <- est(rep(1, nrow(Dm)))
    boot  <- map(seq_len(n_boot), ~ est(tabulate(sample(nrow(Dm), replace = TRUE), nrow(Dm))))

    pi_ci <- map_dfr(boot, "pi", .id = "b") %>%
      summarise(pi_lo = quantile(value, 0.025), pi_hi = quantile(value, 0.975), .by = unit)
    pair_ci <- map_dfr(boot, "pairs", .id = "b") %>%
      summarise(across(c(dxy, da, fst), list(lo = ~ quantile(.x, 0.025), hi = ~ quantile(.x, 0.975))),
                .by = c(unit_a, unit_b))

    list(pi    = point$pi %>% rename(pi = value) %>% left_join(pi_ci, by = "unit"),
         pairs = point$pairs %>% left_join(pair_ci, by = c("unit_a", "unit_b")))
  }))

pi_tab <- res %>% select(k, set, out) %>% mutate(out = map(out, "pi")) %>% unnest(out) %>%
  left_join(n_rep, by = c("k", "set", "unit"))
pair_tab <- res %>% select(k, set, out) %>% mutate(out = map(out, "pairs")) %>% unnest(out)



# checks -----------------------------------------------------------------------

## check if every locus was summed for every level and sample set
n_loc <- agg %>% summarise(n = n_distinct(locus), .by = c(k, set))
if (any(n_loc$n != nrow(loci))) {
  warning("a level or sample set lacks loci")
}

## check if every unit and pair of units has an estimate at every level and in every sample set
if (any(is.na(pi_tab$pi)) || any(is.na(pair_tab$dxy))) {
  warning("units or pairs without an estimate")
}

## check if no unit is more diverse than its divergence from any other (negative d_a)
n_neg <- sum(pair_tab$da < 0)
if (n_neg > 0) {
  message(n_neg, " pairs with negative net divergence (d_a < 0)")
}



# figures ----------------------------------------------------------------------

### pi per unit, one panel per level
plot_pi <- function(kk, y_title) {
  d <- pi_tab %>% filter(k == kk) %>% mutate(set = factor(set, levels = names(lab_set)))

  ggplot(d, aes(x = unit, y = pi, fill = unit, shape = set, group = set)) +
    geom_errorbar(aes(ymin = pi_lo, ymax = pi_hi), width = 0, linewidth = 0.3,
                  colour = "grey30", position = position_dodge(width = 0.6)) +
    geom_point(size = 2, stroke = 0.3, colour = "black", position = position_dodge(width = 0.6)) +
    scale_fill_manual(values = pal_unit, guide = "none") +
    scale_shape_manual(values = shape_set, labels = lab_set, name = NULL) +
    scale_y_continuous(limits = c(0, max(pi_tab$pi_hi)), expand = expansion(mult = c(0, 0.05))) +
    guides(shape = guide_legend(override.aes = list(fill = "grey70"))) +
    labs(x = sprintf("Units at K = %s", sub("k", "", kk)),
         y = if (y_title) expression(pi ~ "(per site)") else NULL)
}

widths <- table(distinct(pi_tab, k, unit)$k)[k_levels]

p_pi <- wrap_plots(map2(k_levels, seq_along(k_levels) == 1, plot_pi), nrow = 1, widths = widths) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

save_fig("nucleotide_diversity.pdf", p_pi, 7.1, 3)

### Hudson's F_ST between units (all representatives), one panel per level
plot_fst <- function(kk) {
  d <- pair_tab %>% filter(k == kk, set == "all")
  u <- sort(unique(c(d$unit_a, d$unit_b)))
  d <- d %>% mutate(unit_a = factor(unit_a, levels = u), unit_b = factor(unit_b, levels = rev(u)))

  ggplot(d, aes(x = unit_a, y = unit_b, fill = fst)) +
    geom_tile(colour = "grey80", linewidth = 0.2) +
    geom_text(aes(label = sprintf("%.2f", fst)), size = 6 / .pt) +
    scale_fill_gradient(low = pal_fst[1], high = pal_fst[2], limits = c(0, max(pair_tab$fst)),
                        name = expression(F[ST])) +
    coord_equal() +
    labs(x = sprintf("Units at K = %s", sub("k", "", kk)), y = NULL) +
    theme(axis.line = element_blank(), axis.ticks = element_blank())
}

p_fst <- wrap_plots(map(k_levels, plot_fst), nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right")

save_fig("fst.pdf", p_fst, 7.1, 2.8)



# tables -----------------------------------------------------------------------

### pi per level and unit
pub_pi <- pi_tab %>%
  filter(set == "all") %>%
  left_join(pi_tab %>%
              filter(set == "no_mixed") %>%
              select(k, unit, n_nm = n, pi_nm = pi, lo_nm = pi_lo, hi_nm = pi_hi),
            by = c("k", "unit")) %>%
  arrange(k, unit) %>%
  transmute(Level = sub("k", "K = ", k), Unit = unit, n = n,
            `pi [95% CI]` = fmt_ci(pi, pi_lo, pi_hi),
            `Without mixed MAT: n` = n_nm,
            `Without mixed MAT: pi [95% CI]` = fmt_ci(pi_nm, lo_nm, hi_nm))

### divergence per level and pair of units
pub_pairs <- pair_tab %>%
  filter(set == "all") %>%
  left_join(pair_tab %>%
              filter(set == "no_mixed") %>%
              select(k, unit_a, unit_b, fst_nm = fst, lo_nm = fst_lo, hi_nm = fst_hi),
            by = c("k", "unit_a", "unit_b")) %>%
  arrange(k, unit_a, unit_b) %>%
  transmute(Level = sub("k", "K = ", k), Pair = paste(unit_a, unit_b, sep = " - "),
            `d_xy [95% CI]` = fmt_ci(dxy, dxy_lo, dxy_hi),
            `d_a [95% CI]`  = fmt_ci(da, da_lo, da_hi),
            `F_ST [95% CI]` = fmt_ci(fst, fst_lo, fst_hi, 3),
            `Without mixed MAT: F_ST [95% CI]` = fmt_ci(fst_nm, lo_nm, hi_nm, 3))

write_tsv(pi_tab, file.path(out_root, "nucleotide_diversity_units.tsv"))
write_tsv(pair_tab, file.path(out_root, "nucleotide_diversity_pairs.tsv"))
write_tsv(pub_pi, file.path(out_root, "publication_nucleotide_diversity.tsv"))
write_tsv(pub_pairs, file.path(out_root, "publication_divergence.tsv"))



# summary ----------------------------------------------------------------------

cat("\n=== pi per level and unit ===\n")
print(pub_pi, n = Inf, width = Inf)

cat("\n=== Divergence per level and pair ===\n")
print(pub_pairs, n = Inf, width = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  nucleotide_diversity.pdf\n")
cat("  fst.pdf\n")
cat("  nucleotide_diversity_units.tsv\n")
cat("  nucleotide_diversity_pairs.tsv\n")
cat("  publication_nucleotide_diversity.tsv\n")
cat("  publication_divergence.tsv\n")
