## ============================================================================
## 051_alignment_divergence.R
## Divergence from the reference and the coverage mask: are the masked columns
## the divergent ones (per column, per tip and per lineage)
##
## Files needed:
##   tip_divergence.tsv     callability and distance to the reference per tip
##   callability_bins.tsv   alternative allele fraction per callability bin
##                          (both from hpc/02_structure_phylogeny/051_alignment_divergence.sh;
##                           copy them from ${ROOT}/02_structure_phylogeny/results/051_results_alignment_divergence
##                           to local/02_structure_phylogeny/data/051_data_alignment_divergence)
##   lineage_labels.tsv     lineage per sample at K = 4 to 6,
##                          from local/03_delimitation/001_lineage_labels.R
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

res_root  <- "02_structure_phylogeny/results"
data_root <- "02_structure_phylogeny/data/051_data_alignment_divergence"
out_root  <- file.path(res_root, "051_results_alignment_divergence")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

f_tips <- file.path(data_root, "tip_divergence.tsv")
f_bins <- file.path(data_root, "callability_bins.tsv")

### lineage labels
f_lineage <- "03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv"



# parameters -------------------------------------------------------------------

k_focus <- c(4, 5, 6)

### level the tips are coloured by (as in the trees)
k_lineage <- 6

### core threshold of hpc 051 (callable in at least this share of the samples)
core_min <- 0.95

### columns callable in fewer than this share of the samples
low_max <- 0.5

### samples listed in the console
n_show <- 20



# style ------------------------------------------------------------------------

col_group <- c("sample"          = "#666666",
               "genome"          = "#E41A1C",
               "outgroup genome" = "#377EB8",
               "outgroup mapped" = "#4DAF4A",
               "reference"       = "#B8860B")

### lineage colours
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

### groups of hpc 051, assemblies split into Cenococcum and the outgroup genus
tips <- read_tsv(f_tips, show_col_types = FALSE) %>%
  mutate(group = case_when(group == "genome" & str_detect(tip, "^Pseflo") ~ "outgroup genome",
                           group == "outgroup"                            ~ "outgroup mapped",
                           TRUE                                           ~ group),
         alt_shift = alt_frac_lost - alt_frac_kept)

bins <- read_tsv(f_bins, show_col_types = FALSE) %>%
  mutate(bin = sprintf("%.3f-%.3f", callability_from, pmin(callability_to, 1)),
         bin = fct_inorder(bin))

### one row per sample and K (clonemates carry their representative's label)
lineages <- read_tsv(f_lineage, show_col_types = FALSE) %>%
  pivot_longer(c(k4, k5, k6), names_to = "K", values_to = "lineage") %>%
  transmute(tip = sample, clustered, K = as.integer(str_remove(K, "^k")), lineage) %>%
  filter(K %in% k_focus)

### population samples (clonemates included, each is its own sequence)
samples <- filter(tips, group == "sample") %>%
  left_join(filter(lineages, K == k_lineage) %>% select(tip, lineage), by = "tip")

### clustered samples (clone representatives) at every K
samples_k <- inner_join(select(samples, -lineage),
                        filter(lineages, clustered) %>% select(-clustered), by = "tip")



# checks -----------------------------------------------------------------------

## check if every group is there
missing_groups <- setdiff(names(col_group), unique(tips$group))

if (length(missing_groups) > 0) {
  warning("groups missing from tip_divergence.tsv: ", paste(missing_groups, collapse = ", "))
}

## check if the labelled samples are in the alignment
missing_labelled <- setdiff(unique(lineages$tip), samples$tip)

if (length(missing_labelled) > 0) {
  warning(length(missing_labelled), " labelled samples not in tip_divergence.tsv")
}

## check if every sample has a lineage
if (any(is.na(samples$lineage))) {
  warning(sum(is.na(samples$lineage)), " samples without a lineage label")
}

## check if every lineage label has a colour
if (any(!lineages$lineage %in% names(pal_unit))) {
  warning("lineage labels without a colour in pal_unit: ",
          paste(setdiff(lineages$lineage, names(pal_unit)), collapse = ", "))
}

## check if the bins are the eight of 051
if (nrow(bins) != 8) {
  warning("expected 8 callability bins, found ", nrow(bins))
}

message("tips: ", nrow(tips), "   samples: ", nrow(samples),
        "   clustered: ", n_distinct(samples_k$tip))



# columns ----------------------------------------------------------------------

### alternative allele fraction in poorly callable columns against the core
### (flat if the mask is independent of divergence)
col_summary <- bins %>%
  mutate(class = case_when(callability_from >= core_min ~ "core",
                           callability_to   <= low_max  ~ "low",
                           TRUE                         ~ "intermediate")) %>%
  group_by(class) %>%
  summarise(alt_frac = weighted.mean(mean_alt_frac, columns, na.rm = TRUE),
            columns  = sum(columns),
            .groups  = "drop")

alt_core <- col_summary$alt_frac[col_summary$class == "core"]

p_bins <- ggplot(bins, aes(x = bin, y = mean_alt_frac, group = 1)) +
  geom_line(colour = "grey60") +
  geom_point(aes(size = columns), shape = 21, colour = "black", fill = "grey30") +
  geom_hline(yintercept = alt_core, colour = "#D73027", linetype = "dashed",
             linewidth = 0.4) +
  scale_size_continuous(trans = "log10", name = "Columns") +
  labs(x = "Column callability (share of samples with a base)",
       y = "Alternative allele fraction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# tips -------------------------------------------------------------------------

### divergence at the columns a sample lost against the columns it kept (above
### the 1:1 line the mask took divergent columns)
shift <- filter(samples, !is.na(alt_frac_lost), !is.na(alt_frac_kept))

shift_test <- wilcox.test(shift$alt_frac_lost, shift$alt_frac_kept,
                          paired = TRUE, exact = FALSE)

p_shift <- ggplot(shift, aes(x = alt_frac_kept, y = alt_frac_lost, fill = lineage)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dashed",
              linewidth = 0.4) +
  geom_point(shape = 21, colour = "black", size = 1, alpha = 0.6) +
  scale_fill_manual(values = pal_unit, na.value = "grey75",
                    name = sprintf("Lineage (K = %d)", k_lineage)) +
  labs(x = "Alternative allele fraction, columns kept",
       y = "Alternative allele fraction, columns lost")

### distance to the reference over the core against callability (negative:
### poorly covered samples look closer to the reference)
dist_cor <- cor.test(samples$callable_pct, samples$p_dist_core,
                     method = "spearman", exact = FALSE)

p_dist <- ggplot(samples, aes(x = callable_pct, y = p_dist_core, fill = lineage)) +
  geom_point(shape = 21, colour = "black", size = 1, alpha = 0.6) +
  scale_fill_manual(values = pal_unit, na.value = "grey75",
                    name = sprintf("Lineage (K = %d)", k_lineage)) +
  labs(x = "Callable columns (%)", y = "p distance to the reference, core columns")



# groups and lineages ----------------------------------------------------------

by_group <- tips %>%
  group_by(group) %>%
  summarise(n                  = n(),
            callable_median    = median(callable_pct),
            p_dist_all_median  = median(p_dist_all),
            p_dist_core_median = median(p_dist_core),
            alt_shift_median   = median(alt_shift, na.rm = TRUE),
            .groups            = "drop")

by_lineage <- samples_k %>%
  group_by(K, lineage) %>%
  summarise(n                  = n(),
            callable_median    = median(callable_pct),
            callable_min       = min(callable_pct),
            p_dist_core_median = median(p_dist_core),
            alt_shift_median   = median(alt_shift, na.rm = TRUE),
            .groups            = "drop")

### callability and distance between lineages, per K
kw <- samples_k %>%
  group_by(K) %>%
  summarise(p_callable = kruskal.test(callable_pct ~ lineage)$p.value,
            p_dist     = kruskal.test(p_dist_core ~ lineage)$p.value,
            .groups    = "drop")

p_lineage <- samples_k %>%
  filter(K == k_lineage) %>%
  pivot_longer(c(callable_pct, p_dist_core), names_to = "measure", values_to = "value") %>%
  mutate(measure = recode(measure,
                          callable_pct = "Callable columns (%)",
                          p_dist_core  = "p distance to the reference, core")) %>%
  ggplot(aes(x = lineage, y = value, fill = lineage)) +
  geom_boxplot(outlier.size = 0.6, linewidth = 0.3) +
  facet_wrap(~ measure, scales = "free_y") +
  scale_fill_manual(values = pal_unit, guide = "none") +
  labs(x = sprintf("Lineage (K = %d)", k_lineage), y = NULL)



# tables -----------------------------------------------------------------------

write_tsv(col_summary, file.path(out_root, "divergence_by_callability.tsv"))
write_tsv(by_group,    file.path(out_root, "divergence_by_group.tsv"))
write_tsv(by_lineage,  file.path(out_root, "divergence_by_lineage.tsv"))
write_tsv(kw,          file.path(out_root, "divergence_lineage_tests.tsv"))

tests <- tibble(
  test      = c("paired Wilcoxon, alternative fraction lost vs kept (samples)",
                "Spearman, callable columns vs p distance over the core (samples)"),
  n         = c(nrow(shift), nrow(samples)),
  statistic = c(as.numeric(shift_test$statistic), as.numeric(dist_cor$estimate)),
  p_value   = c(shift_test$p.value, dist_cor$p.value))

write_tsv(tests, file.path(out_root, "divergence_tests.tsv"))

cat("\n=== columns by callability ===\n")
print(col_summary, n = Inf)

cat("\n=== tips by group ===\n")
print(by_group, n = Inf)

cat("\n=== samples by lineage ===\n")
print(by_lineage, n = Inf)
print(kw, n = Inf)

cat("\n=== tests ===\n")
print(tests, n = Inf)

cat("\n=== the", n_show, "samples whose lost columns were the most divergent ===\n")
shift %>%
  arrange(desc(alt_shift)) %>%
  select(tip, lineage, callable_pct, p_dist_core, alt_frac_kept, alt_frac_lost, alt_shift) %>%
  print(n = n_show)



# figures ----------------------------------------------------------------------

save_fig("divergence_by_callability.pdf", p_bins, 7.1, 5)
save_fig("divergence_lost_vs_kept.pdf", p_shift, 7.1, 6)
save_fig("divergence_vs_callability.pdf", p_dist, 7.1, 6)
save_fig("divergence_by_lineage.pdf", p_lineage, 7.1, 5)

save_fig("divergence_combined.pdf",
         (p_bins + p_shift) / (p_dist + p_lineage) +
           plot_annotation(tag_levels = "A"), 7.1, 6)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  divergence_by_callability.pdf   alternative allele fraction per callability bin\n")
cat("  divergence_lost_vs_kept.pdf     per sample, columns lost against kept\n")
cat("  divergence_vs_callability.pdf   distance to the reference against callability\n")
cat("  divergence_by_lineage.pdf       callability and distance per lineage\n")
cat("  divergence_combined.pdf\n")
cat("  divergence_by_callability.tsv\n")
cat("  divergence_by_group.tsv\n")
cat("  divergence_by_lineage.tsv\n")
cat("  divergence_lineage_tests.tsv\n")
cat("  divergence_tests.tsv\n")
