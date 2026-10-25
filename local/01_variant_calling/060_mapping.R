## ============================================================================
## 060_mapping.R
## Mapping QC per sample against the single-copy orthologue consensus
##
## Files needed:
##   flagstat/*_flagstat.txt        samtools flagstat per sample, from
##                                  hpc/01_variant_calling/060_mapping.sh
##   coverage/*_coverage.txt        samtools coverage per sample, from
##                                  hpc/01_variant_calling/060_mapping.sh
##                                  (copy flagstat and coverage of
##                                   ${ROOT}/01_variant_calling/results/060_results_mapping/consensus
##                                   to local/01_variant_calling/data/060_data_mapping)
##
## What each choice rests on:
##   flagstat, coverage             Danecek et al. 2021, GigaScience 10:giab008
## ============================================================================

library(tidyverse)
library(ggrepel)
library(ggpubr)



# files ------------------------------------------------------------------------

data_root <- "01_variant_calling/data/060_data_mapping"
out_root  <- "01_variant_calling/results/060_results_mapping"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### outgroup samples
outgroup_pattern <- "^M4-"

### breadth at which an orthologue counts as present (%)
sco_present_pct <- 50

### plotting cut-off for near-empty libraries (not a filter)
min_primary_reads <- 1e6



# style ------------------------------------------------------------------------

group_colors <- c(Ingroup = "grey40", Outgroup = "#EE7733")

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

### one row per sample from a samtools flagstat file
parse_flagstat <- function(f) {

  lines  <- readLines(f)
  sample <- sub("_flagstat\\.txt$", "", basename(f))

  extract_count <- function(pattern) {
    line <- grep(pattern, lines, value = TRUE)[1]
    if (is.na(line)) return(NA_real_)
    as.numeric(sub(" \\+.*", "", line))
  }

  primary         <- extract_count("primary$")
  primary_mapped  <- extract_count("primary mapped")
  primary_dup     <- extract_count("primary duplicates$")
  paired          <- extract_count("paired in sequencing")
  properly_paired <- extract_count("properly paired")

  ### NA instead of NaN for an empty total
  safe_pct <- function(x, total) {
    if (is.na(total) || total == 0) return(NA_real_)
    x / total * 100
  }

  tibble(
    sample          = sample,
    primary_reads   = primary,
    primary_mapped  = primary_mapped,
    duplicates      = primary_dup,
    ### sequenced reads placed on the reference
    mapping_rate    = safe_pct(primary_mapped, primary),
    ### pairs in the expected orientation and distance
    proper_pair_rate = safe_pct(properly_paired, paired),
    ### PCR and optical duplicates among the mapped reads
    duplicate_rate  = safe_pct(primary_dup, primary_mapped)
  )
}


flagstat_files <- list.files(file.path(data_root, "flagstat"),
                             pattern = "_flagstat\\.txt$", full.names = TRUE)

## check if there are flagstat files
if (length(flagstat_files) == 0) stop("no flagstat files in ", data_root)

flagstat_df <- map_dfr(flagstat_files, parse_flagstat)

message("flagstat files read: ", nrow(flagstat_df))


### coverage columns rname, startpos, endpos, numreads, covbases, coverage,
### meandepth, meanbaseq, meanmapq (one row per orthologue)
coverage_files <- list.files(file.path(data_root, "coverage"),
                             pattern = "_coverage\\.txt$", full.names = TRUE)

## check if there are coverage files
if (length(coverage_files) == 0) stop("no coverage files in ", data_root)

coverage_df <- map_dfr(coverage_files, function(f) {

  d <- read_tsv(f, show_col_types = FALSE) %>%
    rename_with(~ sub("^#", "", .x)) %>%
    mutate(contig_len = endpos - startpos + 1)

  tibble(
    sample         = sub("_coverage\\.txt$", "", basename(f)),
    n_scos_total   = nrow(d),
    ### orthologues carrying reads
    n_scos_present = sum(d$coverage >= sco_present_pct),
    reference_size = sum(d$contig_len),
    ### depth over the whole reference (weighted by contig length)
    mean_depth     = sum(d$meandepth * d$contig_len) / sum(d$contig_len),
    ### reference positions covered by at least one read
    breadth_1x     = sum(d$covbases) / sum(d$contig_len) * 100,
    ### mean mapping quality of the placed reads
    mean_mapq      = sum(d$meanmapq * d$numreads) / sum(d$numreads)
  )
})

message("coverage files read: ", nrow(coverage_df))



# per-sample summary -----------------------------------------------------------

sample_summary <- flagstat_df %>%
  full_join(coverage_df, by = "sample") %>%
  mutate(
    is_outgroup = str_detect(sample, outgroup_pattern),
    group       = if_else(is_outgroup, "Outgroup", "Ingroup"),
    ### flag for the log-scale plots (not an exclusion)
    low_reads   = is.na(primary_reads) | primary_reads < min_primary_reads,
    pct_scos_present = n_scos_present / n_scos_total * 100
  ) %>%
  arrange(mapping_rate)

write_tsv(sample_summary, file.path(out_root, "sample_summary.tsv"))

message("\nsamples in the summary: ", nrow(sample_summary),
        "  (", sum(sample_summary$is_outgroup), " outgroup)")



# checks -----------------------------------------------------------------------

## check if every sample has a flagstat and a coverage file
incomplete <- sample_summary %>%
  filter(is.na(mapping_rate) | is.na(mean_depth)) %>%
  select(sample, mapping_rate, mean_depth)

if (nrow(incomplete) > 0) {
  warning("samples with a missing flagstat or coverage file: ",
          paste(incomplete$sample, collapse = ", "))
  print(incomplete, n = Inf)
}

## check if every coverage file lists every orthologue of the reference
if (n_distinct(sample_summary$n_scos_total) > 1) {
  warning("the coverage files do not all list the same number of contigs")
  print(count(sample_summary, n_scos_total))
}



# lowest samples ---------------------------------------------------------------

cat("\n=== ten samples with the lowest mapping rate ===\n")
sample_summary %>%
  select(sample, group, mapping_rate, mean_depth, breadth_1x,
         pct_scos_present) %>%
  head(10) %>%
  print(n = Inf)

cat("\n=== ten samples with the fewest orthologues covered ===\n")
sample_summary %>%
  arrange(pct_scos_present) %>%
  select(sample, group, pct_scos_present, mapping_rate, mean_depth) %>%
  head(10) %>%
  print(n = Inf)



# figures ----------------------------------------------------------------------

### every metric, samples ranked
ranked <- sample_summary %>%
  select(sample, group,
         `Primary mapping rate (%)` = mapping_rate,
         `Mean depth (x)`           = mean_depth,
         `Breadth >=1x (%)`         = breadth_1x,
         `Orthologues covered (%)`  = pct_scos_present,
         `Duplicate rate (%)`       = duplicate_rate,
         `Mean MAPQ`                = mean_mapq) %>%
  pivot_longer(-c(sample, group), names_to = "metric", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(metric = factor(metric,
                         levels = c("Primary mapping rate (%)",
                                    "Mean depth (x)",
                                    "Breadth >=1x (%)",
                                    "Orthologues covered (%)",
                                    "Duplicate rate (%)",
                                    "Mean MAPQ"))) %>%
  group_by(metric) %>%
  arrange(value, .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup()

p_ranked <- ggplot(ranked, aes(x = rank, y = value, fill = group)) +
  geom_point(shape = 21, colour = "black", size = 0.7, stroke = 0.2) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = group_colors, name = NULL) +
  labs(x = "Samples, ranked", y = NULL) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"))

save_fig("qc_ranked_metrics.pdf", p_ranked, 7.1, 4.5)


### mapping rate against depth (near-empty libraries left out)
scatter_df <- filter(sample_summary, !low_reads, !is.na(mean_depth))

p_scatter <- ggplot(scatter_df, aes(x = mapping_rate, y = mean_depth,
                                    fill = group)) +
  geom_point(shape = 21, colour = "black", size = 1.6, alpha = 0.8) +
  geom_text_repel(data = slice_min(scatter_df, mapping_rate, n = 15),
                  aes(label = sample), size = 2.2, max.overlaps = 20,
                  show.legend = FALSE) +
  scale_fill_manual(values = group_colors, name = NULL) +
  scale_y_log10() +
  labs(x = "Primary mapping rate (%)", y = "Mean depth (x, log scale)") +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

save_fig("qc_mapping_vs_depth.pdf", p_scatter, 7.1, 6)


### orthologues covered against the mapping rate
sco_df <- filter(sample_summary, !is.na(pct_scos_present), !is.na(mapping_rate))

p_sco <- ggplot(sco_df, aes(x = mapping_rate, y = pct_scos_present,
                            fill = group)) +
  geom_point(shape = 21, colour = "black", size = 1.6, alpha = 0.8) +
  geom_text_repel(data = slice_min(sco_df, pct_scos_present, n = 15),
                  aes(label = sample), size = 2.2, max.overlaps = 20,
                  show.legend = FALSE) +
  scale_fill_manual(values = group_colors, name = NULL) +
  labs(x = "Primary mapping rate (%)", y = "Orthologues covered (%)") +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

save_fig("qc_scos_vs_mapping_rate.pdf", p_sco, 7.1, 6)



# tables -----------------------------------------------------------------------

pub_table <- sample_summary %>%
  filter(!is_outgroup) %>%
  summarise(
    n = n(),
    n_depth = sum(!is.na(mean_depth)),
    `Reference size (Mb)`  = sprintf("%.2f", first(na.omit(reference_size)) / 1e6),
    `Orthologues`          = first(na.omit(n_scos_total)),
    `Primary mapped (%)`   = sprintf("%.2f ± %.2f", mean(mapping_rate, na.rm = TRUE),
                                     sd(mapping_rate, na.rm = TRUE)),
    ### duplicates among the reads mapped to the consensus
    `Duplicates (%)`       = sprintf("%.2f ± %.2f", mean(duplicate_rate, na.rm = TRUE),
                                     sd(duplicate_rate, na.rm = TRUE)),
    `Mean depth (x)`       = sprintf("%.2f ± %.2f", mean(mean_depth, na.rm = TRUE),
                                     sd(mean_depth, na.rm = TRUE)),
    ### median and IQR
    `Median depth (x) [IQR]` = sprintf("%.2f [%.2f-%.2f]",
                                       median(mean_depth, na.rm = TRUE),
                                       quantile(mean_depth, 0.25, na.rm = TRUE),
                                       quantile(mean_depth, 0.75, na.rm = TRUE)),
    `Breadth >=1x (%) [IQR]` = sprintf("%.2f [%.2f-%.2f]",
                                       median(breadth_1x, na.rm = TRUE),
                                       quantile(breadth_1x, 0.25, na.rm = TRUE),
                                       quantile(breadth_1x, 0.75, na.rm = TRUE)),
    `Orthologues covered (%) [IQR]` = sprintf("%.2f [%.2f-%.2f]",
                                              median(pct_scos_present, na.rm = TRUE),
                                              quantile(pct_scos_present, 0.25, na.rm = TRUE),
                                              quantile(pct_scos_present, 0.75, na.rm = TRUE)),
    `Mean MAPQ`            = sprintf("%.2f ± %.2f", mean(mean_mapq, na.rm = TRUE),
                                     sd(mean_mapq, na.rm = TRUE))
  )

write_tsv(pub_table, file.path(out_root, "publication_mapping_table.tsv"))

cat("\n=== mapping statistics, ingroup ===\n")
print(pub_table, n = Inf, width = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  sample_summary.tsv\n")
cat("  qc_ranked_metrics.pdf\n")
cat("  qc_mapping_vs_depth.pdf\n")
cat("  qc_scos_vs_mapping_rate.pdf\n")
cat("  publication_mapping_table.tsv\n")
