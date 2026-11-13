## ============================================================================
## 013_tajimas_d.R
## Tajima's D per unit at K = 4, 5 and 6 from the single-copy orthologues, every
## site projected to the same sample size; per-gene D summarised over genes,
## with Watterson's theta and pi per site; all representatives and without
## mixed-MAT samples
##
## Files needed:
##   groups.tsv                level, sample set, unit, size and group
##   locus_group_sums.tsv.gz   per locus and group: sites, pi, expected segregating sites
##                             (both from hpc/04_diversity_niche_selection/013_tajimas_d.sh;
##                              copy them from ${ROOT}/04_diversity_niche_selection/results/013_results_tajimas_d
##                              to local/04_diversity_niche_selection/data/013_data_tajimas_d)
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

data_root <- "04_diversity_niche_selection/data/013_data_tajimas_d"
out_root  <- "04_diversity_niche_selection/results/013_results_tajimas_d"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### levels
k_levels <- c("k4", "k5", "k6")

### projected sample size (N_PROJ of hpc/04_diversity_niche_selection/013_tajimas_d.sh)
n_proj <- 16

### bootstrap replicates over genes and seed
n_boot <- 1000
seed   <- 42



# style ------------------------------------------------------------------------

### unit colours
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

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
fmt_ci <- function(m, lo, hi, digits = 2) sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
                                                  m, lo, hi)



# read -------------------------------------------------------------------------

groups <- read_tsv(file.path(data_root, "groups.tsv"), col_types = "cccdcc")
sums   <- read_tsv(file.path(data_root, "locus_group_sums.tsv.gz"), col_types = "ccdddd")



# Tajima's D -------------------------------------------------------------------

### constants for the projected sample size
i  <- seq_len(n_proj - 1)
a1 <- sum(1 / i)
a2 <- sum(1 / i^2)
b1 <- (n_proj + 1) / (3 * (n_proj - 1))
b2 <- 2 * (n_proj^2 + n_proj + 3) / (9 * n_proj * (n_proj - 1))
c1 <- b1 - 1 / a1
c2 <- b2 - (n_proj + 2) / (a1 * n_proj) + a2 / a1^2
e1 <- c1 / a1
e2 <- c2 / (a1^2 + a2)

tajima <- function(pi, S) if_else(S > 0, (pi - S / a1) / sqrt(e1 * S + e2 * S * (S - 1)), NA_real_)

### per gene and group
per_locus <- sums %>%
  mutate(D = tajima(pi_sum, seg_proj))



# summaries per group ----------------------------------------------------------

### per level, sample set and unit: genes with a D, mean D with a bootstrap CI over genes, median D, share of genes with D < 0, Watterson's theta and pi per site over all genes
boot_ci <- function(d) quantile(replicate(n_boot, mean(sample(d, replace = TRUE))), c(0.025, 0.975))

set.seed(seed)

per_group <- per_locus %>%
  summarise(genes  = sum(!is.na(D)),
            mean_D = mean(D, na.rm = TRUE),
            ci     = list(boot_ci(D[!is.na(D)])),
            median_D = median(D, na.rm = TRUE),
            neg      = mean(D[!is.na(D)] < 0),
            theta_w  = sum(seg_proj) / a1 / sum(sites),
            pi       = sum(pi_sum) / sum(sites),
            sites    = sum(sites),
            .by = group) %>%
  mutate(lo = map_dbl(ci, 1), hi = map_dbl(ci, 2), .after = mean_D) %>%
  select(-ci)

res <- groups %>%
  select(k, set, unit, n, group) %>%
  left_join(per_group, by = "group")



# checks -----------------------------------------------------------------------

## check if every group has genes with a D
if (any(res$genes == 0)) {
  warning("units without any gene with segregating sites")
}

## check if every gene was summed for every group
n_loc <- count(sums, group)
if (n_distinct(n_loc$n) != 1) {
  warning("groups differ in the number of genes")
}



# figures ----------------------------------------------------------------------

### per-gene D (violin, all representatives) with the mean and 95% CI per sample set, one panel per level
plot_level <- function(kk, y_title) {
  g <- filter(res, k == kk)
  v <- per_locus %>%
    inner_join(filter(g, set == "all") %>% select(group, unit), by = "group") %>%
    filter(!is.na(D))
  d <- g %>% mutate(set = factor(set, levels = names(lab_set)))

  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
    geom_violin(data = v, aes(x = unit, y = D), fill = "grey92", colour = "grey70", linewidth = 0.2,
                scale = "width", width = 0.8) +
    geom_errorbar(data = d, aes(x = unit, ymin = lo, ymax = hi, group = set), width = 0, linewidth = 0.3,
                  colour = "grey20", position = position_dodge(width = 0.4)) +
    geom_point(data = d, aes(x = unit, y = mean_D, fill = unit, shape = set, group = set), size = 2,
               stroke = 0.3, colour = "black", position = position_dodge(width = 0.4)) +
    scale_fill_manual(values = pal_unit, guide = "none") +
    scale_shape_manual(values = shape_set, labels = lab_set, name = NULL) +
    guides(shape = guide_legend(override.aes = list(fill = "grey70"))) +
    labs(x = sprintf("Units at K = %s", sub("k", "", kk)), y = if (y_title) "Tajima's D (per gene)" else NULL)
}

widths <- table(distinct(res, k, unit)$k)[k_levels]

p_d <- wrap_plots(map2(k_levels, seq_along(k_levels) == 1, plot_level), nrow = 1, widths = widths) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

save_fig("tajimas_d.pdf", p_d, 7.1, 3)



# tables -----------------------------------------------------------------------

### Tajima's D per level and unit (all representatives and without mixed-MAT samples)
pub <- res %>%
  filter(set == "all") %>%
  left_join(res %>%
              filter(set == "no_mixed") %>%
              select(k, unit, n_nm = n, mean_nm = mean_D, lo_nm = lo, hi_nm = hi),
            by = c("k", "unit")) %>%
  arrange(k, unit) %>%
  transmute(Level = sub("k", "K = ", k), Unit = unit, n = n, Genes = genes,
            `theta_W (per site)` = sprintf("%.4f", theta_w), `pi (per site)` = sprintf("%.4f", pi),
            `Tajima's D, mean [95% CI]` = fmt_ci(mean_D, lo, hi),
            `Median D` = sprintf("%.2f", median_D), `Genes with D < 0 (%)` = sprintf("%.0f", 100 * neg),
            `Without mixed MAT: n, mean D [95% CI]` = sprintf("%d, %s", n_nm, fmt_ci(mean_nm, lo_nm, hi_nm)))

write_tsv(per_locus, file.path(out_root, "tajimas_d_genes.tsv.gz"))
write_tsv(res, file.path(out_root, "tajimas_d_units.tsv"))
write_tsv(pub, file.path(out_root, "publication_tajimas_d.tsv"))



# summary ----------------------------------------------------------------------

cat("\n=== Tajima's D per level and unit ===\n")
print(pub, n = Inf, width = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  tajimas_d.pdf\n")
cat("  tajimas_d_genes.tsv.gz\n")
cat("  tajimas_d_units.tsv\n")
cat("  publication_tajimas_d.tsv\n")
