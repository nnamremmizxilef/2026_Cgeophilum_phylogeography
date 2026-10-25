## ============================================================================
## 130_snp_density.R
## SNP density over the single-copy orthologues, and differences per orthologue
## for pairs on either side of the clone threshold
##
## Files needed:
##   snps_per_sco.tsv               SNPs per orthologue
##   pair_differences_per_sco.tsv   differences per orthologue and pair
##                                  (both from hpc/01_variant_calling/130_snp_distribution.sh;
##                                   copy them from ${ROOT}/01_variant_calling/results/130_results_snp_density
##                                   to local/01_variant_calling/data/130_data_snp_density)
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

data_130 <- "01_variant_calling/data/130_data_snp_density"
out_root <- "01_variant_calling/results/130_results_snp_density"

f_sco  <- file.path(data_130, "snps_per_sco.tsv")
f_pair <- file.path(data_130, "pair_differences_per_sco.tsv")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### clone threshold (clone_dist of local/01_variant_calling/090_clone_id.R)
clone_dist <- 7.117085267241543e-4



# style ------------------------------------------------------------------------

col_cut  <- "grey20"
col_fill <- "grey60"
col_band <- "grey88"

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

### one row per orthologue (orthologues without a SNP included); pair table:
### orthologues at which a pair differs
sco  <- read_tsv(f_sco,  show_col_types = FALSE)
pair <- read_tsv(f_pair, show_col_types = FALSE)



# SNPs over the orthologues ----------------------------------------------------

density_all <- sum(sco$n_snps) / (sum(sco$length) / 1000)
len_mean    <- mean(sco$length)
n_zero      <- sum(sco$n_snps == 0)

dens_med <- median(sco$snps_per_kb)
dens_q25 <- unname(quantile(sco$snps_per_kb, 0.25))
dens_q75 <- unname(quantile(sco$snps_per_kb, 0.75))

message("orthologues: ", nrow(sco),
        "   length: ", sum(sco$length), " bp",
        "   SNPs: ", sum(sco$n_snps))
message("density: ", round(density_all, 1), " per kb",
        "   without a SNP: ", n_zero,
        "   range: ", round(min(sco$snps_per_kb), 1), " - ",
        round(max(sco$snps_per_kb), 1))

summary_table <- tibble(
  n_sco          = nrow(sco),
  total_length   = sum(sco$length),
  mean_length    = len_mean,
  total_snps     = sum(sco$n_snps),
  density_per_kb = density_all,
  n_without_snp  = n_zero,
  density_q25    = dens_q25,
  density_median = dens_med,
  density_q75    = dens_q75,
  density_max    = max(sco$snps_per_kb)
)

write_tsv(summary_table, file.path(out_root, "summary_statistics.tsv"))

cat("\n=== SNPs over the orthologue set ===\n")
print(summary_table, width = Inf)

cat(sprintf("\nsites differing at clone_dist (%.3e), mean-length orthologue  %.2f\n",
            clone_dist, clone_dist * len_mean))



# differences per pair ---------------------------------------------------------

### total differences, orthologues differing and maximum per orthologue
pair_sum <- pair %>%
  group_by(pair) %>%
  summarise(total_diff = sum(n_diff),
            n_sco      = n(),
            max_diff   = max(n_diff),
            .groups = "drop") %>%
  mutate(f_sco   = n_sco / nrow(sco),
         per_sco = total_diff / n_sco) %>%
  arrange(total_diff)

write_tsv(pair_sum, file.path(out_root, "pair_summary.tsv"))

cat("\n=== differences per pair ===\n")
print(pair_sum, n = Inf, width = Inf)



# figures ----------------------------------------------------------------------

### reference line, median and interquartile band (mapped for the legend)
ref_line <- tibble(slope = density_all, intercept = 0, lab = "Mean SNP density")
ref_med  <- tibble(xintercept = dens_med, lab = "Median")
ref_band <- tibble(xmin = dens_q25, xmax = dens_q75, lab = "50% of orthologues")

### interquartile wedge (the same half of the orthologues as the band in B)
ref_wedge <- tibble(length_kb = c(0, max(sco$length) / 1000)) %>%
  mutate(lo = dens_q25 * length_kb, hi = dens_q75 * length_kb,
         lab = ref_band$lab)

### A: SNP count against orthologue length
p_scale <- ggplot(sco, aes(x = length / 1000, y = n_snps)) +
  geom_ribbon(data = ref_wedge, inherit.aes = FALSE,
              aes(x = length_kb, ymin = lo, ymax = hi, fill = lab)) +
  geom_point(shape = 21, size = 0.5, stroke = 0.2, alpha = 0.25,
             fill = col_fill, colour = "black") +
  geom_abline(data = ref_line, inherit.aes = FALSE, key_glyph = "vline",
              aes(slope = slope, intercept = intercept, linetype = lab),
              colour = col_cut, linewidth = 0.5) +
  scale_fill_manual(values = setNames(col_band, ref_wedge$lab[1]), name = NULL) +
  scale_linetype_manual(values = setNames("dashed", ref_line$lab), name = NULL) +
  guides(fill = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
  labs(x = "Orthologue length (kb)",
       y = "N segregating sites")

### B: SNP density
p_density <- ggplot(sco, aes(x = snps_per_kb)) +
  geom_rect(data = ref_band, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = lab),
            colour = NA) +
  geom_histogram(binwidth = 5, fill = col_fill, colour = NA) +
  geom_vline(data = ref_med, inherit.aes = FALSE,
             aes(xintercept = xintercept, linetype = lab),
             colour = col_cut, linewidth = 0.5) +
  scale_fill_manual(values = setNames(col_band, ref_band$lab), name = NULL) +
  scale_linetype_manual(values = setNames("dashed", ref_med$lab), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = "SNPs per kb",
       y = "Orthologues")

save_fig("snp_distribution.pdf",
         p_scale / p_density + plot_annotation(tag_levels = "A"), 7.1, 8.5)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  summary_statistics.tsv\n")
cat("  pair_summary.tsv\n")
cat("  snp_distribution.pdf\n")
