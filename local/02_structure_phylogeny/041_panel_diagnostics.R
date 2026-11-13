## ============================================================================
## 041_panel_diagnostics.R
## Alignment panel of the tree: coverage mask, gene lengths, missing data per tip
## and genes without data
##
## Files needed:
##   mask_fraction.tsv                masked fraction per sample
##   alignment_stats.tsv              seqkit stats per orthogroup
##   tip_missingness.tsv              N content per tip
##   empty_genes_per_sample.tsv       genes without called bases per tip
##   empty_tips_per_gene.tsv          tips without called bases per gene
##                                    (all five from hpc/02_structure_phylogeny/041_panel_diagnostics.sh;
##                                     copy them from ${ROOT}/02_structure_phylogeny/results/041_results_panel_diagnostics
##                                     to local/02_structure_phylogeny/data/041_data_panel_diagnostics)
##   sample_metadata_ALL_FIXED.tsv    metadata, from local/01_variant_calling/103_mating_type.R
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

data_root <- "02_structure_phylogeny/data/041_data_panel_diagnostics"
out_root  <- "02_structure_phylogeny/results/041_results_panel_diagnostics"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

f_mask    <- file.path(data_root, "mask_fraction.tsv")
f_genes   <- file.path(data_root, "alignment_stats.tsv")
f_tips    <- file.path(data_root, "tip_missingness.tsv")
f_empty_s <- file.path(data_root, "empty_genes_per_sample.tsv")
f_empty_g <- file.path(data_root, "empty_tips_per_gene.tsv")
f_meta    <- "01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv"



# parameters -------------------------------------------------------------------

### assemblies and reference sequence (tips without a mapping history)
genome_pattern <- "^(Cenge|Pseflo)"
reference_tip  <- "consensus"

### mapped outgroup samples (other genus, gene loss expected, kept apart from
### the ingroup)
outgroup_pattern <- "^M4-"

### tips listed in the console
n_high <- 20



# style ------------------------------------------------------------------------

col_in_tree <- c(`TRUE` = "black", `FALSE` = "grey65")

col_group <- c(ingroup = "#4575B4", outgroup = "#D73027")

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

meta <- read_tsv(f_meta, show_col_types = FALSE)

tips <- read_tsv(f_tips, show_col_types = FALSE) %>%
  mutate(is_genome    = str_detect(sample, genome_pattern),
         is_reference = sample == reference_tip,
         is_outgroup  = str_detect(sample, outgroup_pattern),
         is_sample    = !is_genome & !is_reference)

genes <- read_tsv(f_genes, show_col_types = FALSE) %>%
  select(length = min_len, n_seqs = num_seqs)

### mean depth from the metadata
mask <- read_tsv(f_mask, show_col_types = FALSE) %>%
  left_join(select(meta, sample, mean_depth = mean_seq_depth), by = "sample") %>%
  mutate(in_tree = sample %in% tips$sample)

empty_sample <- read_tsv(f_empty_s, show_col_types = FALSE)
empty_gene   <- read_tsv(f_empty_g, show_col_types = FALSE)

n_columns <- sum(genes$length)

message("tips: ", nrow(tips), "   genes: ", nrow(genes),
        "   columns: ", n_columns)



# checks -----------------------------------------------------------------------

## check if the metadata covers every masked sample
if (anyNA(mask$mean_depth)) {
  warning(sum(is.na(mask$mean_depth)), " samples without a depth value")
}

## check if every alignment carries the same number of sequences
if (any(genes$n_seqs != genes$n_seqs[1])) {
  stop("the alignments do not all carry the same number of sequences")
}

## check if the two empty tables agree on how much was lost
if (sum(empty_sample$empty_genes) != sum(empty_gene$empty_tips)) {
  stop("the per sample and per gene tables of empty sequences do not agree")
}



# mask -------------------------------------------------------------------------

### masked fraction against depth (mask concentrated in thin samples)
mask_cor <- cor.test(mask$mean_depth, mask$masked, method = "spearman",
                     exact = FALSE)

p_mask <- ggplot(mask, aes(x = mean_depth, y = masked, fill = in_tree)) +
  geom_point(shape = 21, colour = "black", size = 1.1, alpha = 0.7) +
  scale_x_log10() +
  scale_fill_manual(values = col_in_tree, name = NULL,
                    labels = c(`TRUE` = "in the tree",
                               `FALSE` = "removed by filtering")) +
  labs(x = "Mean sequencing depth (x)", y = "Masked reference positions (%)")

save_fig("mask_vs_depth.pdf", p_mask, 7.1, 5)



# genes ------------------------------------------------------------------------

p_genes <- ggplot(genes, aes(x = length)) +
  geom_histogram(bins = 50, fill = "grey40") +
  scale_x_log10() +
  labs(x = "Alignment length (bp)", y = "N orthogroups")

save_fig("gene_lengths.pdf", p_genes, 7.1, 5)



# tips -------------------------------------------------------------------------

sample_tips <- filter(tips, is_sample)

p_tips <- ggplot(sample_tips, aes(x = n_percent)) +
  geom_histogram(binwidth = 1, fill = "grey40") +
  labs(x = "Ambiguous sites in the concatenation (%)", y = "N sample tips")

save_fig("tip_missingness.pdf", p_tips, 7.1, 5)



# empty genes ------------------------------------------------------------------

### diffuse missing data against whole genes lost, per tip
empty <- empty_sample %>%
  left_join(select(tips, sample, n_percent, is_genome, is_reference,
                   is_outgroup, is_sample), by = "sample") %>%
  left_join(select(mask, sample, masked, mean_depth), by = "sample") %>%
  mutate(group = if_else(is_outgroup, "outgroup", "ingroup"))

### mapped tips only (the genomes define the orthogroups and lose nothing)
empty_reads <- filter(empty, !is_genome, !is_reference)

empty_cor <- cor.test(empty_reads$empty_genes, empty_reads$n_percent,
                      method = "spearman", exact = FALSE)

p_empty_vs_n <- ggplot(empty_reads, aes(x = n_percent, y = empty_genes,
                                        fill = group)) +
  geom_point(shape = 21, colour = "black", size = 1.2, alpha = 0.7) +
  scale_x_log10() +
  scale_fill_manual(values = col_group, name = NULL) +
  labs(x = "Ambiguous sites in the concatenation (%)",
       y = "N orthogroups without called bases")

save_fig("empty_genes_vs_missingness.pdf", p_empty_vs_n, 7.1, 5)


p_empty_gene <- ggplot(filter(empty_gene, empty_tips > 0), aes(x = empty_tips)) +
  geom_histogram(bins = 40, fill = "grey40") +
  scale_x_log10() +
  labs(x = "N tips without called bases", y = "N orthogroups")

save_fig("empty_tips_per_gene.pdf", p_empty_gene, 7.1, 5)



# tables -----------------------------------------------------------------------

panel_summary <- tibble(
  quantity = c("samples with a mask", "samples in the tree", "orthogroups",
               "tips", "columns", "masked %, median", "masked %, maximum",
               "tip N %, median", "tip N %, maximum",
               "genes with an empty tip", "tips with an empty gene",
               "empty genes per tip, median", "empty genes per tip, maximum"),
  value = c(nrow(mask), sum(mask$in_tree), nrow(genes), nrow(tips), n_columns,
            round(median(mask$masked), 2), round(max(mask$masked), 2),
            round(median(sample_tips$n_percent), 2),
            round(max(sample_tips$n_percent), 2),
            sum(empty_gene$empty_tips > 0),
            sum(empty_sample$empty_genes > 0),
            median(empty_reads$empty_genes),
            max(empty_reads$empty_genes)))

write_tsv(panel_summary, file.path(out_root, "panel_summary.tsv"))

### test statistics
panel_statistics <- tribble(
  ~test,                                          ~statistic,                    ~p_value,
  "Spearman, masked fraction vs mean depth",      as.numeric(mask_cor$estimate), mask_cor$p.value,
  "Spearman, empty genes vs tip N",               as.numeric(empty_cor$estimate), empty_cor$p.value)

write_tsv(panel_statistics, file.path(out_root, "panel_statistics.tsv"))

cat("\n=== panel ===\n")
print(panel_summary, n = Inf)

cat("\n=== statistics ===\n")
print(panel_statistics, n = Inf)

cat("\n=== the", n_high, "tips with the most N ===\n")
sample_tips %>%
  left_join(select(meta, sample, project, country), by = "sample") %>%
  arrange(desc(n_percent)) %>%
  select(sample, project, country, n_percent) %>%
  print(n = n_high)

cat("\n=== the", n_high, "tips losing the most whole genes ===\n")
empty_reads %>%
  left_join(select(meta, sample, project, country), by = "sample") %>%
  arrange(desc(empty_genes)) %>%
  select(sample, project, country, empty_genes, n_percent, masked) %>%
  print(n = n_high)



# figures ----------------------------------------------------------------------

save_fig("panel_diagnostics.pdf",
         (p_mask + p_genes) / (p_tips + p_empty_vs_n) +
           plot_annotation(tag_levels = "A"), 7.1, 5.9)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  mask_vs_depth.pdf\n")
cat("  gene_lengths.pdf\n")
cat("  tip_missingness.pdf\n")
cat("  empty_genes_vs_missingness.pdf\n")
cat("  empty_tips_per_gene.pdf\n")
cat("  panel_diagnostics.pdf\n")
cat("  panel_summary.tsv\n")
cat("  panel_statistics.tsv\n")
