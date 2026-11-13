## ============================================================================
## 060_site_concordance.R
## Site concordance factors (sCF) per branch of the decloned species tree
##
## Files needed:
##   site_concordance.cf.stat   sCF per branch, from hpc/02_structure_phylogeny/060_site_concordance.sh
##   site_concordance.cf.tree   the tree carrying sCF, from hpc/02_structure_phylogeny/060_site_concordance.sh
##                              (copy both files of
##                               ${ROOT}/02_structure_phylogeny/results/060_results_site_concordance
##                               to local/02_structure_phylogeny/data/060_data_site_concordance)
## ============================================================================

library(ape)
library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

data_root <- "02_structure_phylogeny/data/060_data_site_concordance"
out_root  <- "02_structure_phylogeny/results/060_results_site_concordance"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

f_stat    <- file.path(data_root, "site_concordance.cf.stat")
f_tree    <- file.path(data_root, "site_concordance.cf.tree")



# parameters -------------------------------------------------------------------

### sCF of a branch without signal (one of three quartet resolutions)
scf_null <- 100 / 3

### support threshold
support_high <- 95

### seed
seed <- 42



# style ------------------------------------------------------------------------

col_point <- "#4575B4"
col_null  <- "#D73027"

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

### sCF per branch (comment header skipped)
cf <- read_tsv(f_stat, comment = "#", show_col_types = FALSE) %>%
  rename_with(~ str_replace(.x, "^sCF_N$", "sCF_n")) %>%
  mutate(across(where(is.character), ~ na_if(.x, "NA")))

### branch support of 052 (branch label)
cf <- cf %>%
  mutate(support = suppressWarnings(as.numeric(Label)),
         sCF     = as.numeric(sCF),
         sN      = as.numeric(sN),
         Length  = as.numeric(Length))

tree <- read.tree(f_tree)

message("branches with sCF: ", sum(!is.na(cf$sCF)),
        "   of ", nrow(cf))



# checks -----------------------------------------------------------------------

## check if the expected columns are present (layout differs by IQ-TREE version)
required <- c("ID", "sCF", "sN", "Label", "Length")

if (!all(required %in% names(cf))) {
  stop("columns missing from the cf.stat file: ",
       paste(setdiff(required, names(cf)), collapse = ", "))
}

## check if any branch sits below the quartet floor
below_null <- sum(cf$sCF < scf_null, na.rm = TRUE)

if (below_null > 0) {
  warning(below_null, " branches have sCF below the ",
          sprintf("%.1f", scf_null), "% floor")
}



# distribution -----------------------------------------------------------------

cf_ok <- filter(cf, !is.na(sCF))

scf_summary <- cf_ok %>%
  summarise(branches   = n(),
            min        = min(sCF),
            q1         = quantile(sCF, 0.25),
            median     = median(sCF),
            mean       = mean(sCF),
            q3         = quantile(sCF, 0.75),
            max        = max(sCF),
            below_null = sum(sCF < scf_null))

p_dist <- ggplot(cf_ok, aes(x = sCF)) +
  geom_histogram(binwidth = 1, fill = "grey40") +
  geom_vline(xintercept = scf_null, colour = col_null, linetype = "dashed",
             linewidth = 0.6) +
  labs(x = "Site concordance factor (%)", y = "Internal branches")

save_fig("scf_distribution.pdf", p_dist, 7.1, 5)



# support ----------------------------------------------------------------------

### sCF against branch support (bootstrap saturates on long alignments, sCF does not)
cf_sup <- filter(cf_ok, !is.na(support))

p_support <- ggplot(cf_sup, aes(x = support, y = sCF)) +
  geom_hline(yintercept = scf_null, colour = col_null, linetype = "dashed",
             linewidth = 0.5) +
  geom_vline(xintercept = support_high, colour = "grey50", linetype = "dotted",
             linewidth = 0.5) +
  geom_point(shape = 21, colour = "black", fill = col_point, size = 1, alpha = 0.5) +
  labs(x = "Branch support", y = "Site concordance factor (%)")

save_fig("scf_vs_support.pdf", p_support, 7.1, 5)



# branch length ----------------------------------------------------------------

### sCF against branch length (short branches carry few informative sites)
len_cor <- cor.test(cf_ok$Length, cf_ok$sCF, method = "spearman", exact = FALSE)

p_length <- ggplot(cf_ok, aes(x = Length, y = sCF)) +
  geom_hline(yintercept = scf_null, colour = col_null, linetype = "dashed",
             linewidth = 0.5) +
  geom_point(shape = 21, colour = "black", fill = col_point, size = 1, alpha = 0.5) +
  scale_x_log10() +
  labs(x = "Branch length", y = "Site concordance factor (%)")

save_fig("scf_vs_branchlength.pdf", p_length, 7.1, 5)



# informative sites ------------------------------------------------------------

### sCF against the number of decisive sites (sN)
p_sites <- ggplot(cf_ok, aes(x = sN, y = sCF)) +
  geom_hline(yintercept = scf_null, colour = col_null, linetype = "dashed",
             linewidth = 0.5) +
  geom_point(shape = 21, colour = "black", fill = col_point, size = 1, alpha = 0.5) +
  scale_x_log10() +
  labs(x = "Decisive sites per branch", y = "Site concordance factor (%)")

save_fig("scf_vs_sites.pdf", p_sites, 7.1, 5)



# tables -----------------------------------------------------------------------

write_tsv(scf_summary, file.path(out_root, "scf_summary.tsv"))
write_tsv(arrange(cf_ok, sCF), file.path(out_root, "scf_per_branch.tsv"))

cat("\n=== site concordance ===\n")
print(scf_summary, n = Inf)

cat("\n=== sCF against support, branch length and decisive sites ===\n")
cat(sprintf("%-42s %d\n", "branches with support", nrow(cf_sup)))
cat(sprintf("%-42s %d\n", sprintf("support >= %d and sCF below the floor", support_high),
            sum(cf_sup$support >= support_high & cf_sup$sCF < scf_null)))
cat(sprintf("%-42s %.3f\n", "Spearman rho, sCF against branch length", len_cor$estimate))
cat(sprintf("%-42s %.3g\n", "p", len_cor$p.value))
cat(sprintf("%-42s %s\n", "median decisive sites per branch",
            format(median(cf_ok$sN), big.mark = " ")))

cat("\n=== branches with high support and no site signal ===\n")
cf_sup %>%
  filter(support >= support_high, sCF < scf_null) %>%
  arrange(sCF) %>%
  select(ID, support, sCF, sN, Length) %>%
  print(n = 20)



# figures ----------------------------------------------------------------------

save_fig("scf_combined.pdf",
         (p_dist + p_support) / (p_length + p_sites) +
           plot_annotation(tag_levels = "A"), 7.1, 5.9)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  scf_distribution.pdf\n")
cat("  scf_vs_support.pdf\n")
cat("  scf_vs_branchlength.pdf\n")
cat("  scf_vs_sites.pdf\n")
cat("  scf_combined.pdf\n")
cat("  scf_summary.tsv\n")
cat("  scf_per_branch.tsv\n")
