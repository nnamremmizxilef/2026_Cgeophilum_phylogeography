## ============================================================================
## 070_variant_calling.R
## Summary of the raw SCO consensus call set
##
## Files needed:
##   consensus_raw_all_variants_metrics.txt
##                                  per-variant INFO metrics of the raw call set, from
##                                  hpc/01_variant_calling/070_variant_calling_consensus.sh
##                                  (copy it from
##                                   ${ROOT}/01_variant_calling/results/070_results_variant_calling_consensus
##                                   to local/01_variant_calling/data/070_data_variant_calling)
##
## What each choice rests on:
##   variant calling (bcftools)     Danecek et al. 2021, GigaScience 10:giab008
## ============================================================================

library(tidyverse)
library(hexbin)
library(ggpubr)



# files ------------------------------------------------------------------------

data_root <- "01_variant_calling/data/070_data_variant_calling"
out_root  <- "01_variant_calling/results/070_results_variant_calling"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

metrics_file <- file.path(data_root,
                          "consensus_raw_all_variants_metrics.txt")



# style ------------------------------------------------------------------------

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

metrics <- read_tsv(metrics_file, show_col_types = FALSE,
                    na = c(".", "NA", ""))

message("loaded ", nrow(metrics), " variants (complete table)")



# call set summary -------------------------------------------------------------

summary_table <- tibble(
  n_variants     = nrow(metrics),
  n_contigs      = n_distinct(metrics$CHROM),
  max_an         = max(metrics$AN, na.rm = TRUE),
  median_qual    = median(metrics$QUAL, na.rm = TRUE),
  median_mq      = median(metrics$MQ, na.rm = TRUE),
  median_info_dp = median(metrics$DP, na.rm = TRUE),
  median_missing = median(metrics$F_MISSING, na.rm = TRUE),
  median_maf     = median(metrics$MAF, na.rm = TRUE)
)

write_tsv(summary_table, file.path(out_root, "summary_statistics.tsv"))

cat("\n=== raw call set ===\n")
print(summary_table, width = Inf)



# figures ----------------------------------------------------------------------

### site missingness against minor allele frequency
p_bivar <- ggplot(metrics, aes(x = MAF, y = F_MISSING)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(name = "Variants", trans = "log10") +
  labs(caption = paste0("all ", format(nrow(metrics), big.mark = " "),
                        " raw variants"),
       x = "MAF", y = "Fraction of missing genotypes") +
  theme(panel.grid.minor = element_blank(),
        plot.caption = element_text(size = 8, colour = "grey40"))

save_fig("missingness_vs_maf.pdf", p_bivar, 7.1, 6)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  summary_statistics.tsv\n")
cat("  missingness_vs_maf.pdf\n")
