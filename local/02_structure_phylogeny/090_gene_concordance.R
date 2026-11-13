## ============================================================================
## 090_gene_concordance.R
## Gene concordance factors on the four reference topologies
##
## Files needed:
##   gene_concordance_concat.cf.stat            gCF per branch on the 052 tree
##   gene_concordance_aster_collapsed.cf.stat   gCF per branch on the 080 tree
##   gene_concordance_aster.cf.stat             gCF per branch on the 081 tree
##   gene_concordance_caster.cf.stat            gCF per branch on the 082 tree
##   gene_models.tsv                            rate model per gene, counted
##                                  (all five from hpc/02_structure_phylogeny/090_gene_concordance.sh;
##                                   copy them from ${ROOT}/02_structure_phylogeny/results/090_results_gene_concordance
##                                   to local/02_structure_phylogeny/data/090_data_gene_concordance)
##   site_concordance.cf.stat                   sCF per branch of the 052 tree,
##                                  from hpc/02_structure_phylogeny/060_site_concordance.sh
##                                  (copy ${ROOT}/02_structure_phylogeny/results/060_results_site_concordance/site_concordance.cf.stat
##                                   to local/02_structure_phylogeny/data/090_data_gene_concordance)
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

data_root <- "02_structure_phylogeny/data/090_data_gene_concordance"
out_root  <- "02_structure_phylogeny/results/090_results_gene_concordance"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

### one entry per reference topology
f_cf <- c(concat          = file.path(data_root, "gene_concordance_concat.cf.stat"),
          aster_collapsed = file.path(data_root, "gene_concordance_aster_collapsed.cf.stat"),
          aster           = file.path(data_root, "gene_concordance_aster.cf.stat"),
          caster          = file.path(data_root, "gene_concordance_caster.cf.stat"))

f_models <- file.path(data_root, "gene_models.tsv")
f_scf    <- file.path(data_root, "site_concordance.cf.stat")



# parameters -------------------------------------------------------------------

### gCF without signal (three quartet resolutions)
gcf_null <- 100 / 3

### high branch support
support_high <- 95

### few decisive genes
gn_low <- 100

### branches printed per table
n_show <- 20



# style ------------------------------------------------------------------------

### reference colours (as in 091)
col_ref <- c(concat = "#EE6677", aster_collapsed = "#4477AA", aster = "#228833",
             caster = "#AA3377")

col_null <- "#D73027"

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

read_cf <- function(path, reference) {
  read_tsv(path, comment = "#", show_col_types = FALSE) %>%
    mutate(reference = reference,
           gCF     = as.numeric(gCF),
           gN      = as.numeric(gN),
           gDF1    = as.numeric(gDF1),
           gDF2    = as.numeric(gDF2),
           Length  = as.numeric(Length),
           support = suppressWarnings(as.numeric(Label)))
}

cf <- imap(f_cf, ~ read_cf(.x, .y)) %>% bind_rows()

models <- read_tsv(f_models, show_col_types = FALSE)



# checks -----------------------------------------------------------------------

## check if the expected columns are present (cf.stat layout differs between IQ-TREE versions)
required <- c("ID", "gCF", "gN", "Label", "Length")

if (!all(required %in% names(cf))) {
  stop("columns missing from the cf.stat file: ",
       paste(setdiff(required, names(cf)), collapse = ", "))
}

## check if the branches have enough decisive genes
low_gn <- cf %>%
  group_by(reference) %>%
  summarise(below = sum(gN < gn_low, na.rm = TRUE), .groups = "drop")

if (any(low_gn$below > 0)) {
  warning("branches with fewer than ", gn_low, " decisive genes: ",
          paste(sprintf("%s %d", low_gn$reference, low_gn$below), collapse = ", "))
}



# distribution -----------------------------------------------------------------

cf_ok <- filter(cf, !is.na(gCF))

gcf_summary <- cf_ok %>%
  group_by(reference) %>%
  summarise(branches   = n(),
            min        = min(gCF),
            q1         = quantile(gCF, 0.25),
            median     = median(gCF),
            mean       = mean(gCF),
            q3         = quantile(gCF, 0.75),
            max        = max(gCF),
            below_null = sum(gCF < gcf_null),
            gN_median  = median(gN),
            .groups    = "drop")

p_dist <- ggplot(cf_ok, aes(x = gCF, fill = reference)) +
  geom_histogram(binwidth = 1, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = gcf_null, colour = col_null, linetype = "dashed",
             linewidth = 0.6) +
  scale_fill_manual(values = col_ref, name = NULL) +
  labs(x = "Gene concordance factor (%)", y = "Internal branches")

save_fig("gcf_distribution.pdf", p_dist, 7.1, 5)



# decisive genes ---------------------------------------------------------------

### gCF against decisive genes (few: no information; many: conflict)
p_gn <- ggplot(cf_ok, aes(x = gN, y = gCF, fill = reference)) +
  geom_hline(yintercept = gcf_null, colour = col_null, linetype = "dashed",
             linewidth = 0.5) +
  geom_vline(xintercept = gn_low, colour = "grey50", linetype = "dotted",
             linewidth = 0.5) +
  geom_point(shape = 21, colour = "black", size = 1, alpha = 0.5) +
  scale_x_log10() +
  scale_fill_manual(values = col_ref, name = NULL) +
  labs(x = "Decisive gene trees per branch", y = "Gene concordance factor (%)")

save_fig("gcf_vs_decisive_genes.pdf", p_gn, 7.1, 5)



# support ----------------------------------------------------------------------

### gCF against branch support
cf_sup <- filter(cf_ok, !is.na(support))

p_support <- ggplot(cf_sup, aes(x = support, y = gCF, fill = reference)) +
  geom_hline(yintercept = gcf_null, colour = col_null, linetype = "dashed",
             linewidth = 0.5) +
  geom_point(shape = 21, colour = "black", size = 1, alpha = 0.5) +
  scale_fill_manual(values = col_ref, name = NULL) +
  labs(x = "Branch support (UFBoot or local bootstrap, %)", y = "Gene concordance factor (%)")

save_fig("gcf_vs_support.pdf", p_support, 7.1, 5)



# site concordance -------------------------------------------------------------

### gCF against sCF on the concatenation tree
scf <- read_tsv(f_scf, comment = "#", show_col_types = FALSE) %>%
  select(ID, sCF) %>%
  mutate(sCF = as.numeric(sCF))

plane <- cf_ok %>%
  filter(reference == "concat") %>%
  inner_join(scf, by = "ID") %>%
  filter(!is.na(sCF))

if (nrow(plane) > 0) {

  plane_cor <- cor.test(plane$gCF, plane$sCF, method = "spearman",
                        exact = FALSE)

  p_plane <- ggplot(plane, aes(x = gCF, y = sCF, fill = support)) +
    geom_hline(yintercept = gcf_null, colour = col_null, linetype = "dashed",
               linewidth = 0.4) +
    geom_vline(xintercept = gcf_null, colour = col_null, linetype = "dashed",
               linewidth = 0.4) +
    geom_point(shape = 21, colour = "black", size = 1.1, alpha = 0.7) +
    scale_fill_viridis_c(name = "Support") +
    labs(x = "Gene concordance factor (%)",
         y = "Site concordance factor (%)")

  save_fig("gcf_vs_scf.pdf", p_plane, 7.1, 5)

  cat(sprintf("\n=== gCF against sCF, concatenation tree: %d branches, Spearman rho = %.3f ===\n",
              nrow(plane), plane_cor$estimate))

  write_tsv(plane, file.path(out_root, "gcf_scf_per_branch.tsv"))
}



# models -----------------------------------------------------------------------

### rate model per gene (substitution matrix fixed to GTR)
p_models <- ggplot(models, aes(x = reorder(model, genes), y = genes)) +
  geom_col(fill = "grey40") +
  coord_flip() +
  labs(x = NULL, y = "Orthogroups")

save_fig("gene_models.pdf", p_models, 7.1, 4)



# tables -----------------------------------------------------------------------

write_tsv(gcf_summary, file.path(out_root, "gcf_summary.tsv"))
write_tsv(arrange(cf_ok, reference, gCF),
          file.path(out_root, "gcf_per_branch.tsv"))

cat("\n=== gene concordance ===\n")
print(gcf_summary, n = Inf)

cat("\n=== branches with high support and no gene signal ===\n")
cf_sup %>%
  filter(support >= support_high, gCF < gcf_null) %>%
  arrange(gCF) %>%
  select(reference, ID, support, gCF, gN, Length) %>%
  print(n = n_show)

### decisive genes that disagree (conflict, not absence of information)
cat("\n=== conflict, not absence of information ===\n")
cf_ok %>%
  filter(gN >= gn_low, gCF < gcf_null) %>%
  arrange(gCF) %>%
  select(reference, ID, gCF, gDF1, gDF2, gN) %>%
  print(n = n_show)



# figures ----------------------------------------------------------------------

save_fig("gcf_combined.pdf",
         (p_dist + p_gn) / (p_support + plot_spacer()) +
           plot_annotation(tag_levels = "A"), 7.1, 5.9)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  gcf_distribution.pdf\n")
cat("  gcf_vs_decisive_genes.pdf\n")
cat("  gcf_vs_support.pdf\n")
cat("  gcf_vs_scf.pdf\n")
cat("  gene_models.pdf\n")
cat("  gcf_combined.pdf\n")
cat("  gcf_summary.tsv\n")
cat("  gcf_per_branch.tsv\n")
cat("  gcf_scf_per_branch.tsv\n")
