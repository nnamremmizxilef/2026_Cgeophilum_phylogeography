## ============================================================================
## 080_filtering.R
## Hard-filtering diagnostics: site metrics, cut-off sensitivity, sample
## missingness, orthologue recovery and collapsed paralogues
##
## Files needed:
##   consensus_raw_all_variants_metrics.txt
##                                  per-variant metrics of the raw call set, from
##                                  hpc/01_variant_calling/070_variant_calling_consensus.sh
##                                  (copy it from
##                                   ${ROOT}/01_variant_calling/results/070_results_variant_calling_consensus
##                                   to local/01_variant_calling/data/070_data_variant_calling)
##   consensus_sco_reference.fasta.fai
##                                  index of the consensus reference, from
##                                  hpc/01_variant_calling/053_indexing.sh
##                                  (copy it from
##                                   ${ROOT}/01_variant_calling/results/050_results_sco_consensus
##                                   to local/01_variant_calling/data/070_data_variant_calling)
##   consensus_ingroup_noMAF_metrics.txt.gz
##                                  per-variant metrics of the filtered ingroup set,
##                                  from hpc/01_variant_calling/080_filtering.sh
##   variant_counts.tsv             variants and samples per VCF, from
##                                  hpc/01_variant_calling/080_filtering.sh
##   consensus_sample_missingness.tsv
##                                  missing genotypes per sample, from
##                                  hpc/01_variant_calling/080_filtering.sh
##   sco_presence.tsv               orthologues present per sample, from
##                                  hpc/01_variant_calling/080_filtering.sh
##                                  (copy these four files from
##                                   ${ROOT}/01_variant_calling/results/080_results_filtering
##                                   to local/01_variant_calling/data/080_data_filtering)
##
## What each choice rests on:
##   site, genotype and sample      Danecek et al. 2021, GigaScience 10:giab008
##   filters (bcftools)
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

data_070 <- "01_variant_calling/data/070_data_variant_calling"
data_080 <- "01_variant_calling/data/080_data_filtering"
out_root <- "01_variant_calling/results/080_results_filtering"

f_raw      <- file.path(data_070, "consensus_raw_all_variants_metrics.txt")
f_filt     <- file.path(data_080, "consensus_ingroup_noMAF_metrics.txt.gz")
f_counts   <- file.path(data_080, "variant_counts.tsv")
f_smiss    <- file.path(data_080, "consensus_sample_missingness.tsv")
f_sco      <- file.path(data_080, "sco_presence.tsv")
f_fai      <- file.path(data_070, "consensus_sco_reference.fasta.fai")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### thresholds of hpc/01_variant_calling/080_filtering.sh
qual        <- 30     ### QUAL, site level
mq          <- 50     ### INFO/MQ, site level
site_miss   <- 0.10   ### INFO/F_MISSING, site level
maf01       <- 0.01   ### lenient MAF level
maf05       <- 0.05   ### stringent MAF level
sample_miss <- 0.50   ### per-sample missingness
sco_breadth <- 50     ### breadth at which an orthologue counts as covered

### orthologue share for a second sample criterion (shown, not applied)
sco_min_pct <- 50

### outgroup samples
outgroup_pattern <- "^M4-"



# style ------------------------------------------------------------------------

col_stage  <- c(Raw = "grey65", Filtered = "#228833")
col_status <- c(Kept = "grey55", Removed = "#CC3311", Outgroup = "#EE7733")
col_cut    <- "#D55E00"

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

num_cols <- c("QUAL", "DP", "MQ", "MAF", "F_MISSING", "AF", "SCR", "SGB")

read_metrics <- function(path, stage) {
  read_tsv(path, show_col_types = FALSE, na = c(".", "NA", "")) %>%
    mutate(across(any_of(num_cols), as.numeric),
           across(any_of(c("NS", "AC", "AN")), as.integer),
           stage = stage)
}

raw_metrics  <- read_metrics(f_raw,  "Raw")
filt_metrics <- read_metrics(f_filt, "Filtered")

metrics <- bind_rows(raw_metrics, filt_metrics) %>%
  mutate(stage = factor(stage, levels = c("Raw", "Filtered")))

message("raw: ", nrow(raw_metrics), "   filtered: ", nrow(filt_metrics))

counts <- read_tsv(f_counts, show_col_types = FALSE) %>%
  mutate(maf_level  = factor(maf_level, levels = c("noMAF", "maf01", "maf05")),
         sample_set = factor(sample_set, levels = c("ingroup", "plusoutgroup")))

sco <- read_tsv(f_sco, show_col_types = FALSE)

### status per sample (outgroup exempt from the missingness filter)
smiss <- read_tsv(f_smiss, show_col_types = FALSE) %>%
  mutate(is_outgroup = str_detect(sample, outgroup_pattern),
         removed = !is_outgroup & f_missing > sample_miss,
         status  = case_when(is_outgroup ~ "Outgroup",
                             removed     ~ "Removed",
                             TRUE        ~ "Kept")) %>%
  left_join(sco, by = "sample") %>%
  arrange(f_missing) %>%
  mutate(rank = row_number())

message("samples: ", nrow(smiss),
        "   outgroup: ", sum(smiss$is_outgroup),
        "   removed: ", sum(smiss$removed))



# site metrics and variant counts ----------------------------------------------

summary_table <- metrics %>%
  group_by(stage) %>%
  summarise(
    n_variants       = n(),
    median_qual      = median(QUAL, na.rm = TRUE),
    median_mq        = median(MQ, na.rm = TRUE),
    median_site_miss = median(F_MISSING, na.rm = TRUE),
    median_maf       = median(MAF, na.rm = TRUE),
    pct_pass_qual    = mean(QUAL >= qual, na.rm = TRUE) * 100,
    pct_pass_mq      = mean(MQ >= mq, na.rm = TRUE) * 100,
    pct_pass_miss    = mean(F_MISSING <= site_miss, na.rm = TRUE) * 100,
    pct_pass_maf01   = mean(MAF >= maf01, na.rm = TRUE) * 100,
    pct_pass_maf05   = mean(MAF >= maf05, na.rm = TRUE) * 100,
    .groups = "drop"
  )

write_tsv(summary_table, file.path(out_root, "summary_statistics.tsv"))
print(summary_table, width = Inf)

print(counts, n = Inf, width = Inf)



# site cut-off sensitivity -----------------------------------------------------

raw_long <- raw_metrics %>%
  transmute(QUAL, MQ, `Site missingness` = F_MISSING, MAF) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
  filter(!is.na(value))

cut_spec <- tribble(
  ~metric,             ~keep,    ~from, ~to,  ~by,
  "QUAL",              "above",  0,     1000, 5,
  "MQ",                "above",  0,     60,   1,
  "Site missingness",  "below",  0,     1,    0.01,
  "MAF",               "above",  0,     0.5,  0.005
)

retention <- cut_spec %>%
  pmap_dfr(function(metric, keep, from, to, by) {
    v <- raw_long$value[raw_long$metric == metric]
    x <- seq(from, to, by = by)
    frac <- if (keep == "above") map_dbl(x, ~ mean(v >= .x)) else
      map_dbl(x, ~ mean(v <= .x))
    tibble(metric = metric, x = x, frac = frac)
  }) %>%
  mutate(metric = factor(metric, levels = cut_spec$metric))

write_tsv(retention, file.path(out_root, "site_cutoff_sensitivity.tsv"))



# sample missingness -----------------------------------------------------------

write_tsv(smiss, file.path(out_root, "per_sample_missingness.tsv"))

smiss %>%
  filter(removed) %>%
  select(sample, f_missing, pct_scos_present) %>%
  arrange(desc(f_missing)) %>%
  write_tsv(file.path(out_root, "samples_removed_missingness.tsv"))

### ingroup samples removed per cut-off
smiss_sens <- smiss %>%
  filter(!is_outgroup) %>%
  select(sample, f_missing) %>%
  expand_grid(x = seq(0, 1, by = 0.005)) %>%
  group_by(x) %>%
  summarise(n_removed = sum(f_missing > x), .groups = "drop")

write_tsv(smiss_sens, file.path(out_root, "sample_cutoff_sensitivity.tsv"))



# orthologue recovery ----------------------------------------------------------

sco_dat <- filter(smiss, !is.na(pct_scos_present))

sco_dat %>%
  select(sample, n_scos_present, n_scos_total, pct_scos_present,
         f_missing, status) %>%
  arrange(pct_scos_present) %>%
  write_tsv(file.path(out_root, "sco_presence_per_sample.tsv"))

cat("\n=== twenty samples with the fewest orthologues recovered ===\n")
sco_dat %>%
  arrange(pct_scos_present) %>%
  select(sample, n_scos_present, n_scos_total, pct_scos_present,
         f_missing, status) %>%
  head(20) %>%
  print(n = Inf)

### missingness criterion against orthologue criterion, ingroup only
compare <- sco_dat %>%
  filter(!is_outgroup) %>%
  transmute(
    sample, f_missing, pct_scos_present, n_scos_present, n_scos_total,
    fail_miss = f_missing > sample_miss,
    fail_sco  = pct_scos_present < sco_min_pct,
    category  = case_when(fail_miss &  fail_sco ~ "both",
                          fail_miss & !fail_sco ~ "missingness only",
                          !fail_miss & fail_sco ~ "orthologues only",
                          TRUE                  ~ "neither")
  )

write_tsv(compare, file.path(out_root, "sample_criteria_flags.tsv"))



# collapsed paralogues ---------------------------------------------------------

### variant density against depth per contig (collapsed paralogues carry
### excess depth and excess variants)
contig_len <- read_tsv(f_fai, col_names = c("CHROM", "length", "offset",
                                            "linebases", "linewidth"),
                       show_col_types = FALSE) %>%
  select(CHROM, length)

dp_median <- median(filt_metrics$DP, na.rm = TRUE)

per_contig <- filt_metrics %>%
  group_by(CHROM) %>%
  summarise(n_variants = n(),
            median_dp  = median(DP, na.rm = TRUE),
            median_maf = median(MAF, na.rm = TRUE),
            .groups = "drop") %>%
  right_join(contig_len, by = "CHROM") %>%
  replace_na(list(n_variants = 0)) %>%
  mutate(per_kb      = n_variants / (length / 1000),
         depth_ratio = median_dp / dp_median)

write_tsv(per_contig, file.path(out_root, "per_contig_density_depth.tsv"))

message("contigs: ", nrow(per_contig),
        "   median depth: ", dp_median,
        "   median density: ", round(median(per_contig$per_kb, na.rm = TRUE), 1),
        " per kb")

### share of all variants on contigs with excess depth
walk(c(1.5, 2, 3), function(k) {
  d <- filter(per_contig, depth_ratio >= k)
  message("depth >= ", k, "x median: ", nrow(d), " contigs, ",
          sum(d$n_variants), " variants (",
          round(sum(d$n_variants) / sum(per_contig$n_variants) * 100, 1), "%)")
})

cat("\n=== twenty contigs with the highest density ===\n")
per_contig %>%
  arrange(desc(per_kb)) %>%
  select(CHROM, length, n_variants, per_kb, depth_ratio, median_maf) %>%
  head(20) %>%
  print(n = Inf)



# figures ----------------------------------------------------------------------

### site metrics before and after filtering, with cut-offs
metric_levels <- c("QUAL (log10)", "MQ", "Site missingness", "MAF")

metrics_long <- metrics %>%
  transmute(stage,
            `QUAL (log10)`     = log10(QUAL),
            `MQ`               = MQ,
            `Site missingness` = F_MISSING,
            `MAF`              = MAF) %>%
  pivot_longer(-stage, names_to = "metric", values_to = "value") %>%
  filter(!is.na(value), is.finite(value)) %>%
  mutate(metric = factor(metric, levels = metric_levels))

### cut-offs (two MAF levels)
cuts_site <- tibble(
  metric = factor(c("QUAL (log10)", "MQ", "Site missingness", "MAF", "MAF"),
                  levels = metric_levels),
  cutoff = c(log10(qual), mq, site_miss, maf01, maf05)
)

p_dist <- ggplot(metrics_long, aes(x = value, y = after_stat(density),
                                   fill = stage)) +
  geom_histogram(bins = 50, alpha = 0.6, position = "identity", colour = NA) +
  geom_vline(data = cuts_site, aes(xintercept = cutoff),
             colour = col_cut, linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~ metric, scales = "free", nrow = 2) +
  scale_fill_manual(values = col_stage, name = NULL) +
  labs(x = NULL, y = "Density")

save_fig("site_metric_distributions.pdf", p_dist, 7.1, 4.5)


### site cut-off sensitivity
cuts_ret <- tibble(
  metric = factor(c("QUAL", "MQ", "Site missingness", "MAF", "MAF"),
                  levels = levels(retention$metric)),
  cutoff = c(qual, mq, site_miss, maf01, maf05)
)

p_ret <- ggplot(retention, aes(x = x, y = frac)) +
  geom_line(linewidth = 0.7, colour = "#228833") +
  geom_vline(data = cuts_ret, aes(xintercept = cutoff),
             colour = col_cut, linetype = "dashed", linewidth = 0.5) +
  facet_wrap(~ metric, scales = "free_x", nrow = 2) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Cutoff value", y = "Raw variants retained")

save_fig("site_cutoff_sensitivity.pdf", p_ret, 7.1, 4.5)


### variants per set
p_counts <- ggplot(counts, aes(x = maf_level, y = n_variants / 1e6,
                               fill = sample_set)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = format(n_variants, big.mark = ",")),
            position = position_dodge(width = 0.8),
            vjust = -0.4, size = 2.6) +
  geom_text(aes(label = paste0("n=", n_samples), y = 0),
            position = position_dodge(width = 0.8),
            vjust = -0.6, size = 2.4, colour = "white") +
  scale_fill_manual(values = c(ingroup = "#228833", plusoutgroup = "#4477AA"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "Variants (millions)")

save_fig("variant_counts_per_set.pdf", p_counts, 7.1, 6)


### per-sample missingness
p_smiss <- ggplot(smiss, aes(x = rank, y = f_missing)) +
  geom_point(aes(fill = status), shape = 21, colour = "black", size = 1) +
  geom_hline(yintercept = sample_miss, colour = col_cut,
             linetype = "dashed", linewidth = 0.5) +
  geom_text(data = filter(smiss, removed | is_outgroup),
            aes(label = sample), hjust = 1.15, size = 2.2,
            check_overlap = TRUE) +
  scale_fill_manual(values = col_status, name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Samples, ordered by missingness", y = "Missing genotypes")

save_fig("per_sample_missingness.pdf", p_smiss, 7.1, 4.3)


### sample cut-off sensitivity
p_smiss_sens <- ggplot(smiss_sens, aes(x = x, y = n_removed)) +
  geom_step(linewidth = 0.7, colour = "#228833") +
  geom_vline(xintercept = sample_miss, colour = col_cut,
             linetype = "dashed", linewidth = 0.5) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Sample missingness cutoff", y = "Ingroup samples removed")

save_fig("sample_cutoff_sensitivity.pdf", p_smiss_sens, 7.1, 6)


### orthologues recovered against missingness
### (labels for the outgroup and for samples failing either criterion)
label_df <- sco_dat %>%
  filter(is_outgroup | f_missing > sample_miss | pct_scos_present < sco_min_pct)

p_scatter <- ggplot(sco_dat, aes(x = f_missing, y = pct_scos_present)) +
  geom_point(aes(fill = status), shape = 21, colour = "black", size = 1.4,
             alpha = 0.8) +
  geom_vline(xintercept = sample_miss, colour = col_cut,
             linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = sco_min_pct, colour = col_cut,
             linetype = "dashed", linewidth = 0.5) +
  geom_text(data = label_df, aes(label = sample),
            hjust = -0.15, size = 2.2, check_overlap = TRUE) +
  scale_fill_manual(values = col_status, name = NULL) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "Missing genotypes", y = "Orthologues recovered (%)")

p_sco_rank <- sco_dat %>%
  arrange(pct_scos_present) %>%
  mutate(sco_rank = row_number()) %>%
  ggplot(aes(x = sco_rank, y = pct_scos_present)) +
  geom_point(aes(fill = status), shape = 21, colour = "black", size = 1) +
  geom_hline(yintercept = sco_min_pct, colour = col_cut,
             linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(values = col_status, guide = "none") +
  labs(x = "Samples, ordered by orthologue recovery",
       y = "Orthologues recovered (%)")

save_fig("sco_presence_vs_missingness.pdf",
         p_scatter / p_sco_rank + plot_layout(heights = c(2, 1)) +
           plot_annotation(tag_levels = "A"),
         7.1, 7.8)


### variant density against depth per contig
p_paralog <- ggplot(per_contig, aes(x = depth_ratio, y = per_kb)) +
  geom_point(shape = 21, colour = "black", fill = "grey55", size = 0.5,
             stroke = 0.2, alpha = 0.3) +
  geom_vline(xintercept = c(1.5, 2), colour = col_cut,
             linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = median(per_contig$per_kb, na.rm = TRUE),
             colour = "#228833", linetype = "dashed", linewidth = 0.5) +
  scale_x_continuous(trans = "log2",
                     breaks = c(0.25, 0.5, 1, 2, 4, 8)) +
  labs(x = "Median site depth / genome-wide median",
       y = "Variants per kb")

p_paralog_hist <- ggplot(per_contig, aes(x = per_kb)) +
  geom_histogram(bins = 60, fill = "grey55", colour = NA) +
  geom_vline(xintercept = median(per_contig$per_kb, na.rm = TRUE),
             colour = "#228833", linetype = "dashed", linewidth = 0.5) +
  labs(x = "Variants per kb", y = "Contigs")

save_fig("per_contig_density_depth.pdf",
         p_paralog / p_paralog_hist + plot_layout(heights = c(2, 1)) +
           plot_annotation(tag_levels = "A"),
         7.1, 7.1)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  summary_statistics.tsv\n")
cat("  site_cutoff_sensitivity.tsv\n")
cat("  per_sample_missingness.tsv\n")
cat("  samples_removed_missingness.tsv\n")
cat("  sample_cutoff_sensitivity.tsv\n")
cat("  sco_presence_per_sample.tsv\n")
cat("  sample_criteria_flags.tsv\n")
cat("  per_contig_density_depth.tsv\n")
cat("  site_metric_distributions.pdf\n")
cat("  site_cutoff_sensitivity.pdf\n")
cat("  variant_counts_per_set.pdf\n")
cat("  per_sample_missingness.pdf\n")
cat("  sample_cutoff_sensitivity.pdf\n")
cat("  sco_presence_vs_missingness.pdf\n")
cat("  per_contig_density_depth.pdf\n")
