## ============================================================================
## 021_phi_test.R
## Share of genes with PHI evidence for recombination per unit at K = 4, 5 and
## 6, for all representatives, without mixed-MAT samples and in random subsets
## of equal size
##
## Files needed:
##   phi_per_gene.tsv    PHI p values per sample set, unit and gene
##   unit_levels.tsv     units of each level
##                       (both from hpc/03_delimitation/021_phi_test.sh;
##                        copy them from ${ROOT}/03_delimitation/results/021_results_phi_test
##                        to local/03_delimitation/data/021_data_phi_test)
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

data_root <- "03_delimitation/data/021_data_phi_test"
out_root  <- "03_delimitation/results/021_results_phi_test"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### significance per gene (also the share of genes expected without recombination)
alpha <- 0.05



# style ------------------------------------------------------------------------

### unit colours (a K4 or K5 unit takes the colour of its larger K6 part)
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

### sample set labels and shapes
lab_set   <- c(all = "All representatives", no_mixed = "Without mixed-MAT samples",
               equal_n = "Random subsets of equal size")
shape_set <- c(all = 21, no_mixed = 24, equal_n = 22)

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

phi <- read_tsv(file.path(data_root, "phi_per_gene.tsv"), na = c("NA", "--"),
                col_types = cols(set = "c", unit = "c", gene = "c", .default = "d"))

lev <- read_tsv(file.path(data_root, "unit_levels.tsv"), col_types = "cc")



# summary per sample set and unit ----------------------------------------------

### genes tested, genes with p < alpha, BH count and binomial test against alpha
per_set <- phi %>%
  summarise(n_seq              = first(n_seq),
            n_seq_distinct     = n_distinct(n_seq),
            genes              = n(),
            tested             = sum(!is.na(p_perm)),
            median_informative = median(n_informative[!is.na(p_perm)]),
            n_sig              = sum(p_perm < alpha, na.rm = TRUE),
            n_fdr              = sum(p.adjust(p_perm[!is.na(p_perm)], "BH") < alpha),
            .by = c(set, unit)) %>%
  mutate(group   = str_remove(set, "_[0-9]+$"),
         pct_sig = 100 * n_sig / tested,
         p_binom = map2_dbl(n_sig, tested,
                            ~ binom.test(.x, .y, alpha, alternative = "greater")$p.value))

### per sample-set group (mean and range over the equal-size draws)
per_group <- per_set %>%
  summarise(n_seq   = first(n_seq),
            draws   = n(),
            pct     = mean(pct_sig),
            pct_min = min(pct_sig),
            pct_max = max(pct_sig),
            .by = c(group, unit))

n_sub <- per_group$n_seq[per_group$group == "equal_n"][1]

lab_set["equal_n"] <- sprintf("Random subsets of %d representatives", n_sub)



# checks -----------------------------------------------------------------------

## check if every sample set and unit holds every gene
if (n_distinct(per_set$genes) != 1) {
  warning("sample sets or units differ in the number of genes")
}

## check if every gene subset holds the same samples within a sample set and unit
if (any(per_set$n_seq_distinct != 1)) {
  warning(sum(per_set$n_seq_distinct != 1),
          " sample sets and units with genes of differing sample size")
}

## check if every unit was tested in every sample set
if (nrow(per_set) != n_distinct(per_set$set) * n_distinct(lev$unit)) {
  warning("not every unit was tested in every sample set")
}

## check if the equal-size draws have one size across units
if (n_distinct(per_set$n_seq[per_set$group == "equal_n"]) != 1) {
  warning("the equal-size draws differ in size between units")
}

## check if dropping mixed-MAT samples never added samples
n_cmp <- per_group %>%
  select(group, unit, n_seq) %>%
  pivot_wider(names_from = group, values_from = n_seq)

if (any(n_cmp$no_mixed > n_cmp$all)) {
  warning("a unit holds more samples without the mixed-MAT samples than with them")
}



# figures ----------------------------------------------------------------------

### one panel per level
plot_level <- function(kk, y_title) {
  d <- per_group %>%
    inner_join(filter(lev, k == kk), by = "unit") %>%
    mutate(group = factor(group, levels = names(lab_set)))

  ggplot(d, aes(x = unit, y = pct, fill = unit, shape = group, group = group)) +
    geom_hline(yintercept = 100 * alpha, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
    geom_errorbar(aes(ymin = pct_min, ymax = pct_max), width = 0, linewidth = 0.3,
                  colour = "grey30", position = position_dodge(width = 0.6)) +
    geom_point(size = 2, stroke = 0.3, colour = "black", position = position_dodge(width = 0.6)) +
    scale_fill_manual(values = pal_unit, guide = "none") +
    scale_shape_manual(values = shape_set, labels = lab_set, name = NULL, drop = FALSE) +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
    guides(shape = guide_legend(override.aes = list(fill = "grey70"))) +
    labs(x = sprintf("Units at K = %s", sub("k", "", kk)),
         y = if (y_title) "Genes with PHI p < 0.05 (%)" else NULL)
}

levels_k <- sort(unique(lev$k))
p_phi <- wrap_plots(map2(levels_k, seq_along(levels_k) == 1, plot_level), nrow = 1,
                    widths = table(lev$k)[levels_k]) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

save_fig("phi_recombination.pdf", p_phi, 7.1, 3)



# tables -----------------------------------------------------------------------

### one row per level and unit
pub <- lev %>%
  left_join(filter(per_set, set == "all"), by = "unit") %>%
  left_join(per_group %>%
              filter(group == "no_mixed") %>%
              select(unit, n_no_mixed = n_seq, pct_no_mixed = pct),
            by = "unit") %>%
  left_join(per_group %>%
              filter(group == "equal_n") %>%
              select(unit, pct_equal = pct, pct_equal_min = pct_min, pct_equal_max = pct_max),
            by = "unit") %>%
  transmute(Level                                   = sub("k", "K = ", k),
            Unit                                    = unit,
            n                                       = n_seq,
            `Genes tested`                          = tested,
            `Informative sites (median per gene)`   = median_informative,
            `Genes p < 0.05`                        = sprintf("%d (%.1f%%)", n_sig, pct_sig),
            `Genes q < 0.05 (BH)`                   = n_fdr,
            `Binomial p`                            = format.pval(p_binom, digits = 2, eps = 1e-16),
            `Without mixed MAT: n`                  = n_no_mixed,
            `Without mixed MAT: %`                  = sprintf("%.1f", pct_no_mixed),
            !!sprintf("Equal size (n = %d): %% mean [range]", n_sub) :=
              sprintf("%.1f [%.1f-%.1f]", pct_equal, pct_equal_min, pct_equal_max))

write_tsv(per_set %>% select(-n_seq_distinct), file.path(out_root, "phi_summary_sets.tsv"))
write_tsv(pub, file.path(out_root, "publication_phi.tsv"))

cat("\n=== Genes with PHI p < 0.05 (%) per level and unit ===\n")
per_group %>%
  select(group, unit, pct) %>%
  pivot_wider(names_from = group, values_from = pct) %>%
  inner_join(lev, by = "unit") %>%
  relocate(k) %>%
  arrange(k, unit) %>%
  mutate(across(where(is.double), ~ round(.x, 1))) %>%
  print(n = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  phi_recombination.pdf\n")
cat("  phi_summary_sets.tsv\n")
cat("  publication_phi.tsv\n")
