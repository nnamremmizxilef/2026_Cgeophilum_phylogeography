## ============================================================================
## 030_multiqc.R
## Read QC per library from the MultiQC reports of FastQC and fastp
##
## Files needed:
##   pre*/multiqc_data.json         FastQC on the raw reads, from
##                                  hpc/01_variant_calling/030_multiqc.sh
##   post*/multiqc_data.json        FastQC on the trimmed reads, from
##                                  hpc/01_variant_calling/030_multiqc.sh
##   fastp*/multiqc_data.json       fastp reports, from
##                                  hpc/01_variant_calling/030_multiqc.sh
##                                  (copy the three multiqc_report_*_data folders of
##                                   ${ROOT}/01_variant_calling/results/030_results_multiqc
##                                   to local/01_variant_calling/data/030_data_multiqc,
##                                   renamed to start with pre (raw), post (trimmed)
##                                   and fastp)
##
## What each choice rests on:
##   read trimming, pair counts     Chen et al. 2018, Bioinformatics 34:i884-i890
##   FastQC module verdicts         Andrews 2010, Babraham Bioinformatics (software)
##   MultiQC reports                Ewels et al. 2016, Bioinformatics 32:3047-3048
## ============================================================================

library(tidyverse)
library(jsonlite)
library(ggpubr)



# files ------------------------------------------------------------------------

data_root <- "01_variant_calling/data/030_data_multiqc"
out_root  <- "01_variant_calling/results/030_results_multiqc"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

### report folders by their first word
all_dirs   <- list.dirs(data_root, recursive = FALSE)
stage_dirs <- map_chr(c(raw = "^pre", trimmed = "^post", fastp = "^fastp"),
                      ~ all_dirs[grep(.x, basename(all_dirs))[1]])



# parameters -------------------------------------------------------------------

### FastQC modules whose verdict is kept
modules <- c("adapter_content", "per_base_sequence_quality",
             "per_sequence_quality_scores", "overrepresented_sequences")



# style ------------------------------------------------------------------------

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

### nested report entry as text values (fastp JSON is nested, FastQC read
### length is a number or a range)
flatten <- function(x, prefix = "") {

  out <- list()

  for (nm in names(x)) {
    value <- x[[nm]]
    key   <- if (prefix == "") nm else paste(prefix, nm, sep = "_")
    if (is.list(value)) {
      out <- c(out, flatten(value, key))
    } else if (is.atomic(value) && length(value) == 1) {
      out[[key]] <- as.character(value)
    }
  }
  out
}

### one row per entry of one module (read files for FastQC, libraries for fastp)
read_module <- function(dir, module) {

  raw <- fromJSON(file.path(dir, "multiqc_data.json"),
                  simplifyVector = FALSE)$report_saved_raw_data
  raw <- raw[[grep(module, names(raw))[1]]]

  imap_dfr(raw, ~ bind_cols(tibble(entry = .y), as_tibble(flatten(.x))))
}

### library id from the entry name (stage tokens anywhere, mate suffix only at
### the end)
library_id <- function(x) {
  x %>%
    str_remove("\\.(fastq|fq)(\\.gz)?$") %>%
    str_remove_all(regex("(_|\\.)(fastp|trimmed|trim|val|paired|filtered|clean)[0-9]*",
                         ignore_case = TRUE)) %>%
    str_remove("(_|\\.)R?[12]$")
}

### fastp: pair counts, retention and Q30
fastp <- read_module(stage_dirs[["fastp"]], "fastp") %>%
  transmute(
    library      = library_id(entry),
    ### pairs (fastp counts reads)
    pairs_raw    = as.numeric(summary_before_filtering_total_reads) / 2,
    pairs_trim   = as.numeric(summary_after_filtering_total_reads) / 2,
    q30_bases    = as.numeric(summary_after_filtering_q30_bases),
    bases_trim   = as.numeric(summary_after_filtering_total_bases),
    ### R1 only
    read_length  = as.numeric(summary_after_filtering_read1_mean_length),
    pct_retained = pairs_trim / pairs_raw * 100,
    ### post-trimming Q30 (over the bases fastp kept)
    pct_q30      = q30_bases / bases_trim * 100,
    ### adapter-trimmed reads over all raw reads (not fastp's pct_adapter)
    adapter_rate = as.numeric(adapter_cutting_adapter_trimmed_reads) /
      (pairs_raw * 2) * 100
  )

### FastQC: module verdicts and read counts, before and after trimming
flags <- map_dfr(c("raw", "trimmed"), function(stage) {

  df <- read_module(stage_dirs[[stage]], "fastqc") %>%
    ### one spelling across MultiQC versions
    rename_with(~ .x %>%
                  str_replace_all(" ", "_") %>%
                  str_to_lower() %>%
                  str_remove("^statuses_"))

  count_col <- str_subset(names(df), "total_sequences$")[1]

  df %>%
    transmute(stage   = stage,
              file    = entry,
              library = library_id(entry),
              mate    = if_else(str_detect(entry, "(_|\\.)R?2(_|\\.|$)"),
                                "R2", "R1"),
              n_reads = as.numeric(.data[[count_col]]),
              across(all_of(modules)))
})



# per-library summary ----------------------------------------------------------

status_counts <- flags %>%
  pivot_longer(all_of(modules), names_to = "module", values_to = "status") %>%
  count(stage, module, status) %>%
  pivot_wider(names_from = status, values_from = n, values_fill = 0) %>%
  arrange(module, stage)

### worst verdict over both mates
worst <- function(x) {
  ranks <- c(pass = 1, warn = 2, fail = 3)
  names(ranks)[max(ranks[x])]
}

summary_df <- fastp %>%
  left_join(
    flags %>%
      filter(stage == "trimmed") %>%
      group_by(library) %>%
      summarise(across(all_of(modules), worst), .groups = "drop"),
    by = "library"
  ) %>%
  arrange(pairs_trim)



# checks -----------------------------------------------------------------------

## check if every library has both mates at both stages
gaps <- flags %>%
  count(library, stage) %>%
  complete(library, stage, fill = list(n = 0)) %>%
  filter(n != 2)

if (nrow(gaps) > 0) {
  warning(nrow(gaps), " library/stage combinations do not have two read files, ",
          "first: ", paste(head(sprintf("%s (%s, n = %d)", gaps$library,
                                        gaps$stage, gaps$n), 5),
                           collapse = ", "))
}

## check if every library has FastQC verdicts
if (anyNA(summary_df$adapter_content)) {
  warning(sum(is.na(summary_df$adapter_content)), " libraries without FastQC ",
          "verdicts, check library_id() against the entry names")
}

## check if the FastQC count of the trimmed R1 file equals the fastp pair count
count_check <- flags %>%
  filter(stage == "trimmed", mate == "R1") %>%
  inner_join(select(summary_df, library, pairs_trim), by = "library") %>%
  filter(n_reads != pairs_trim)

if (nrow(count_check) > 0) {
  warning(nrow(count_check), " libraries where the FastQC read count and the ",
          "fastp pair count disagree, first: ",
          paste(head(count_check$library, 5), collapse = ", "))
}



# figures ----------------------------------------------------------------------

### retention against library size
p_retention <- ggplot(summary_df, aes(pairs_raw / 1e6, pct_retained)) +
  geom_point(shape = 21, colour = "black", fill = "grey30", size = 1.8,
             alpha = 0.7) +
  scale_x_log10() +
  labs(x = "Raw read pairs (millions, log scale)", y = "Pairs retained (%)")

### distributions of the reported metrics
p_dist <- summary_df %>%
  transmute(`Trimmed read pairs (millions)` = pairs_trim / 1e6,
            `Pairs retained (%)`            = pct_retained,
            `Bases at Q30 or above (%)`     = pct_q30) %>%
  pivot_longer(everything(), names_to = "metric") %>%
  ggplot(aes(value)) +
  geom_histogram(bins = 40, fill = "grey40") +
  facet_wrap(~ metric, nrow = 1, scales = "free") +
  labs(x = NULL, y = "Libraries")

save_fig("qc_retention_vs_size.pdf", p_retention, 7.1, 5)
save_fig("qc_metric_distributions.pdf", p_dist, 7.1, 2.6)



# tables -----------------------------------------------------------------------

pub_table <- summary_df %>%
  summarise(
    n                        = n(),
    ### totals
    `Total raw pairs (G)`    = sprintf("%.2f", sum(pairs_raw) / 1e9),
    `Total trimmed pairs (G)` = sprintf("%.2f", sum(pairs_trim) / 1e9),
    `Total trimmed bases (Tb)` = sprintf("%.2f", sum(bases_trim) / 1e12),
    `Raw read pairs (M)`     = sprintf("%.1f ± %.1f", mean(pairs_raw / 1e6),
                                       sd(pairs_raw / 1e6)),
    `Trimmed read pairs (M)` = sprintf("%.1f ± %.1f", mean(pairs_trim / 1e6),
                                       sd(pairs_trim / 1e6)),
    ### range over the whole cohort
    `Trimmed pairs range (M)` = sprintf("%.2f-%.1f", min(pairs_trim) / 1e6,
                                        max(pairs_trim) / 1e6),
    ### median and IQR (library sizes are right-skewed)
    `Median trimmed pairs (M) [IQR]` = sprintf("%.1f [%.1f-%.1f]",
                                               median(pairs_trim / 1e6),
                                               quantile(pairs_trim / 1e6, 0.25),
                                               quantile(pairs_trim / 1e6, 0.75)),
    `Pairs retained (%)`     = sprintf("%.1f ± %.1f", mean(pct_retained),
                                       sd(pct_retained)),
    ### pooled over all bases
    `Bases >=Q30 (%)`        = sprintf("%.1f", sum(q30_bases) /
                                         sum(bases_trim) * 100),
    ### mean over libraries
    `Bases >=Q30 per library (%)` = sprintf("%.1f ± %.1f", mean(pct_q30),
                                            sd(pct_q30)),
    ### median and IQR (a few libraries carry most of the adapter)
    `Adapter in reads (%) [IQR]` = sprintf("%.1f [%.1f-%.1f]",
                                           median(adapter_rate),
                                           quantile(adapter_rate, 0.25),
                                           quantile(adapter_rate, 0.75)),
    `Read length R1 (bp)`    = sprintf("%.0f ± %.0f", mean(read_length),
                                       sd(read_length))
  )

write_tsv(summary_df,    file.path(out_root, "read_qc_summary.tsv"))
write_tsv(status_counts, file.path(out_root, "fastqc_status_counts.tsv"))
write_tsv(pub_table,     file.path(out_root, "publication_read_table.tsv"))

### smallest libraries
cat(sprintf("\nsmallest library   %s, %.0f read pairs\nnext smallest      %.1f M read pairs\n",
            slice_min(summary_df, pairs_trim, n = 1, with_ties = FALSE)$library,
            min(summary_df$pairs_trim), sort(summary_df$pairs_trim)[2] / 1e6))



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  read_qc_summary.tsv\n")
cat("  fastqc_status_counts.tsv\n")
cat("  publication_read_table.tsv\n")
cat("  qc_retention_vs_size.pdf\n")
cat("  qc_metric_distributions.pdf\n")
