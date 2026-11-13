## ============================================================================
## 015_clustering_publication.R
## Publication figures of the clustering: choice of K, sNMF ancestry at the
## levels carried forward and a PCA of the primary draw, by lineage
##
## Files needed:
##   replicate_curves.tsv       BIC and cross-entropy, primary and replicate draws
##   replicate_ari.tsv          replicate draws against the primary
##   assignment_quality.tsv     assignment per sample and K, row order of the sNMF input
##   input.snmfProject          the sNMF runs
##                              (all four from local/02_structure_phylogeny/014_clustering.R)
##   consensus_phylogeny_ingroup_decloned_maf01_1snp_per_sco.vcf.gz
##                              primary draw, from hpc/02_structure_phylogeny/011_filtering_clustering.sh
##                              (in local/02_structure_phylogeny/data/014_data_clustering, copied for 014)
##   lineage_labels.tsv         lineage names at K = 4 to 6,
##                              from local/03_delimitation/001_lineage_labels.R
## ============================================================================

library(vcfR)
library(LEA)
library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

res_014  <- "02_structure_phylogeny/results/014_results_clustering"
data_014 <- "02_structure_phylogeny/data/014_data_clustering"
out_root <- "02_structure_phylogeny/results/015_results_clustering_publication"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

f_curves  <- file.path(res_014, "summary", "replicate_curves.tsv")
f_rep_ari <- file.path(res_014, "summary", "replicate_ari.tsv")
f_quality <- file.path(res_014, "summary", "assignment_quality.tsv")
f_snmf    <- file.path(res_014, "snmf", "input.snmfProject")
f_vcf     <- file.path(data_014, "consensus_phylogeny_ingroup_decloned_maf01_1snp_per_sco.vcf.gz")
f_lineage <- "03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv"



# parameters -------------------------------------------------------------------

### levels carried forward
k_show <- 4:6

### level the samples are ordered and coloured by
k_lineage <- 6

### K range of the curves
k_curve <- 2:12

### PCs drawn
n_pc <- 3



# style ------------------------------------------------------------------------

### lineage colours
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

### projects by shape (colour is reserved for the lineages)
shape_project <- c(magic = 21, phytoakmeter = 24)
lab_project   <- c(magic = "MAGIC", phytoakmeter = "PhytOakmeter")

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

labels <- read_tsv(f_lineage, show_col_types = FALSE)

### one row per clustered sample and K
lineage_long <- labels %>%
  filter(clustered) %>%
  pivot_longer(c(k4, k5, k6), names_to = "K", values_to = "lineage") %>%
  transmute(sample, project, K = as.integer(str_remove(K, "^k")), lineage)

quality <- read_tsv(f_quality, show_col_types = FALSE)

### row order of the sNMF input (014 writes the samples in that order)
sample_ids <- quality %>% filter(K == min(K)) %>% pull(sample)

### curves of the primary and the replicate draws
curves <- read_tsv(f_curves, show_col_types = FALSE)

rep_ari <- read_tsv(f_rep_ari, show_col_types = FALSE)

### ancestry of the run with the lowest cross-entropy, as in 014
snmf_project <- load.snmfProject(f_snmf)

q_long <- map_dfr(k_show, function(k) {
  best <- which.min(cross.entropy(snmf_project, K = k))
  qm   <- Q(snmf_project, K = k, run = best)
  colnames(qm) <- paste0("Anc", seq_len(k))
  as_tibble(qm) %>%
    mutate(sample = sample_ids, K = k, .before = 1) %>%
    pivot_longer(starts_with("Anc"), names_to = "component", values_to = "q")
})

### genotypes of the primary draw (haploid, 0 or 1)
vcf <- read.vcfR(f_vcf, verbose = FALSE)
gt  <- extract.gt(vcf, element = "GT")
geno <- t(matrix(suppressWarnings(as.numeric(gt)), nrow = nrow(gt),
                 dimnames = dimnames(gt)))



# checks -----------------------------------------------------------------------

## check if the sNMF rows are the samples of 014 (hard call of the best run
## against the snmf column)
hard <- q_long %>%
  group_by(K, sample) %>%
  slice_max(q, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(K, sample, snmf_q = as.integer(str_remove(component, "^Anc"))) %>%
  inner_join(select(quality, K, sample, snmf), by = c("K", "sample"))

if (any(hard$snmf_q != hard$snmf)) {
  stop("the sNMF runs do not match assignment_quality.tsv, check that the ",
       "project is the one of the last 014 run")
}

## check if every clustered sample has a lineage
no_label <- setdiff(sample_ids, lineage_long$sample)

if (length(no_label) > 0) {
  warning(length(no_label), " clustered samples without a lineage label")
}

## check if every lineage has a colour
if (any(!lineage_long$lineage %in% names(pal_unit))) {
  warning("lineages without a colour in pal_unit: ",
          paste(setdiff(lineage_long$lineage, names(pal_unit)), collapse = ", "))
}

## check if the VCF carries the clustered samples
if (!setequal(rownames(geno), sample_ids)) {
  warning("the VCF and the sNMF input differ in their samples")
}



# choice of K ------------------------------------------------------------------

curve_long <- curves %>%
  filter(K %in% k_curve) %>%
  pivot_longer(c(BIC, mean_ce), names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric, BIC = "DAPC, BIC", mean_ce = "sNMF, cross-entropy"))

p_k <- ggplot(curve_long, aes(x = K, y = value, group = draw)) +
  annotate("rect", xmin = min(k_show) - 0.4, xmax = max(k_show) + 0.4,
           ymin = -Inf, ymax = Inf, fill = "grey90") +
  geom_line(data = function(d) filter(d, draw != "primary"),
            colour = "grey70", linewidth = 0.3) +
  geom_line(data = function(d) filter(d, draw == "primary"),
            colour = "black", linewidth = 0.6) +
  geom_point(data = function(d) filter(d, draw == "primary"),
             shape = 21, colour = "black", fill = "black", size = 1) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1, strip.position = "left") +
  scale_x_continuous(breaks = k_curve) +
  labs(x = "Number of clusters (K)", y = NULL) +
  theme(strip.placement = "outside", strip.background = element_blank())



# ancestry ---------------------------------------------------------------------

### every sNMF component named by the lineage whose samples carry the most of it
comp_lineage <- q_long %>%
  inner_join(select(lineage_long, sample, K, lineage), by = c("sample", "K")) %>%
  group_by(K, component, lineage) %>%
  summarise(q = mean(q), .groups = "drop_last") %>%
  slice_max(q, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(K, component, comp_lineage = lineage)

## check if every component names a different lineage (one to one at each K)
if (any(count(comp_lineage, K, comp_lineage)$n > 1)) {
  warning("two sNMF components name the same lineage at one K, the colours ",
          "of the barplot are ambiguous there")
}

### samples ordered by lineage, then by the share of their own component
sample_order <- q_long %>%
  filter(K == k_lineage) %>%
  inner_join(comp_lineage, by = c("K", "component")) %>%
  inner_join(filter(lineage_long, K == k_lineage) %>% select(sample, lineage),
             by = "sample") %>%
  filter(comp_lineage == lineage) %>%
  arrange(lineage, desc(q)) %>%
  pull(sample)

bar_dat <- q_long %>%
  inner_join(comp_lineage, by = c("K", "component")) %>%
  inner_join(filter(lineage_long, K == k_lineage) %>% select(sample, lineage),
             by = "sample") %>%
  mutate(sample = factor(sample, levels = sample_order),
         K      = factor(paste("K =", K), levels = paste("K =", k_show)))

p_bar <- ggplot(bar_dat, aes(x = sample, y = q, fill = comp_lineage)) +
  geom_col(width = 1) +
  facet_grid(K ~ lineage, scales = "free_x", space = "free_x", switch = "y") +
  scale_fill_manual(values = pal_unit, guide = "none") +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0), position = "right") +
  labs(x = NULL, y = "Ancestry proportion") +
  theme(axis.text.x      = element_blank(),
        axis.ticks.x     = element_blank(),
        axis.line        = element_blank(),
        strip.background = element_blank(),
        strip.placement  = "outside",
        panel.spacing.x  = unit(0.6, "mm"),
        panel.spacing.y  = unit(1, "mm"))



# PCA --------------------------------------------------------------------------

### mean imputation (as the DAPC input of 014), monomorphic columns dropped
geno_imp <- apply(geno, 2, function(v) {
  v[is.na(v)] <- mean(v, na.rm = TRUE)
  v
})
geno_imp <- geno_imp[, apply(geno_imp, 2, var) > 0]

pca    <- prcomp(geno_imp, center = TRUE, scale. = FALSE)
var_pc <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pc_dat <- as_tibble(pca$x[, seq_len(n_pc)]) %>%
  mutate(sample = rownames(geno_imp)) %>%
  left_join(filter(lineage_long, K == k_lineage) %>% select(sample, project, lineage),
            by = "sample")

pc_plot <- function(a, b) {
  ggplot(pc_dat, aes(x = .data[[a]], y = .data[[b]], fill = lineage,
                     shape = project)) +
    geom_point(size = 1.6, stroke = 0.2, colour = "black", alpha = 0.9) +
    scale_fill_manual(values = pal_unit, name = sprintf("Lineage (K = %d)", k_lineage)) +
    scale_shape_manual(values = shape_project, labels = lab_project, name = NULL) +
    guides(fill  = guide_legend(override.aes = list(shape = 21, size = 2.2), order = 1),
           shape = guide_legend(override.aes = list(size = 2, fill = "grey85"), order = 2)) +
    theme(legend.position = "right", legend.box = "vertical") +
    labs(x = sprintf("%s (%.1f%%)", a, var_pc[as.integer(str_remove(a, "PC"))]),
         y = sprintf("%s (%.1f%%)", b, var_pc[as.integer(str_remove(b, "PC"))]))
}

p_pc12 <- pc_plot("PC1", "PC2")
p_pc13 <- pc_plot("PC1", "PC3")



# replicate draws --------------------------------------------------------------

p_ari <- rep_ari %>%
  pivot_longer(-c(draw, K), names_to = "comparison", values_to = "ARI") %>%
  mutate(comparison = factor(comparison,
                             levels = c("dapc_vs_primary", "snmf_vs_primary", "dapc_vs_snmf"),
                             labels = c("DAPC vs primary", "sNMF vs primary", "DAPC vs sNMF")),
         K = factor(paste("K =", K), levels = paste("K =", sort(unique(K))))) %>%
  ggplot(aes(x = comparison, y = ARI)) +
  geom_jitter(width = 0.12, height = 0, shape = 21, colour = "black",
              fill = "grey55", size = 0.9, stroke = 0.2, alpha = 0.7) +
  stat_summary(fun = median, fun.min = median, fun.max = median,
               geom = "crossbar", width = 0.4, linewidth = 0.3, colour = "#D73027") +
  facet_wrap(~ K, nrow = 1) +
  labs(x = NULL, y = "Adjusted Rand Index") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))



# tables -----------------------------------------------------------------------

### lineage sizes per K (all ingroup samples and clone corrected)
lineage_summary <- labels %>%
  pivot_longer(c(k4, k5, k6), names_to = "K", values_to = "lineage") %>%
  mutate(K = as.integer(str_remove(K, "^k"))) %>%
  group_by(K, lineage) %>%
  summarise(samples      = n(),
            genotypes    = sum(clustered),
            magic        = sum(project == "magic"),
            phytoakmeter = sum(project == "phytoakmeter"),
            sites        = n_distinct(site),
            countries    = n_distinct(country),
            .groups      = "drop")

write_tsv(lineage_summary, file.path(out_root, "lineage_summary.tsv"))
write_tsv(comp_lineage, file.path(out_root, "snmf_component_lineage.tsv"))
write_tsv(select(pc_dat, sample, lineage, project, starts_with("PC")),
          file.path(out_root, "pca_scores.tsv"))
write_tsv(tibble(pc = seq_along(var_pc), variance_pct = var_pc)[seq_len(10), ],
          file.path(out_root, "pca_variance.tsv"))

cat("\n=== lineages per K ===\n")
print(lineage_summary, n = Inf)

cat("\n=== sNMF components and the lineages they name ===\n")
print(comp_lineage, n = Inf)

cat("\n=== PCA ===\n")
cat(sprintf("SNPs           %d\n", ncol(geno_imp)))
cat(sprintf("samples        %d\n", nrow(geno_imp)))
cat(sprintf("PC1, PC2, PC3  %.1f%%, %.1f%%, %.1f%%\n", var_pc[1], var_pc[2], var_pc[3]))



# figures ----------------------------------------------------------------------

save_fig("Fig_clustering_k_choice.pdf", p_k, 3.4, 3.5)
save_fig("Fig_clustering_ancestry.pdf", p_bar, 7.1, 2.8)
save_fig("Fig_clustering_pca.pdf",
         (p_pc12 + p_pc13) + plot_layout(guides = "collect") +
           plot_annotation(tag_levels = "A"), 7.1, 3.1)

save_fig("FigS_clustering_replicates.pdf", p_ari, 7.1, 2.8)

### main figure: K choice and PCA above, ancestry below
p_main <- (p_k | p_pc12) / p_bar +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(tag_levels = "A")

save_fig("Fig_clustering_combined.pdf", p_main, 7.1, 6.7)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  Fig_clustering_k_choice.pdf     BIC and cross-entropy, primary and replicate draws\n")
cat("  Fig_clustering_ancestry.pdf     sNMF ancestry at K =", paste(k_show, collapse = ", "), "\n")
cat("  Fig_clustering_pca.pdf          PC1 to PC3, coloured by lineage\n")
cat("  Fig_clustering_combined.pdf     the three panels together\n")
cat("  FigS_clustering_replicates.pdf  replicate draws against the primary\n")
cat("  lineage_summary.tsv             samples, genotypes, projects, sites per lineage\n")
cat("  snmf_component_lineage.tsv      the lineage each sNMF component stands for\n")
cat("  pca_scores.tsv\n")
cat("  pca_variance.tsv\n")
