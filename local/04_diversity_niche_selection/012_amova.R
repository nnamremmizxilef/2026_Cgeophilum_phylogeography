## ============================================================================
## 012_amova.R
## Hierarchical AMOVA (among units, among sites within units, within sites) at
## K = 4, 5 and 6 on the pairwise sequence distances of the single-copy
## orthologues; all representatives and without mixed-MAT samples
##
## Files needed:
##   classes.tsv               representatives, site, units, mixed-MAT flag
##   pair_differences.tsv.gz   differences per pair of representatives over all loci
##   pair_sites.tsv.gz         compared sites per pair of representatives over all loci
##                             (all three from hpc/04_diversity_niche_selection/010_sequence_distances.sh;
##                              copy them from ${ROOT}/04_diversity_niche_selection/results/010_results_sequence_distances
##                              to local/04_diversity_niche_selection/data/010_data_sequence_distances)
## ============================================================================

library(tidyverse)
library(ggpubr)
library(ade4)



# files ------------------------------------------------------------------------

data_root <- "04_diversity_niche_selection/data/010_data_sequence_distances"
out_root  <- "04_diversity_niche_selection/results/012_results_amova"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### levels
k_levels <- c("k4", "k5", "k6")

### permutations and seed
n_perm <- 999
seed   <- 42



# style ------------------------------------------------------------------------

### variance components (dark to light from the highest level)
lab_source <- c(`among units` = "Among units", `among sites within units` = "Among sites within units",
                `within sites` = "Within sites")
pal_source <- c(`among units` = "grey25", `among sites within units` = "grey60", `within sites` = "grey90")

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

classes <- read_tsv(file.path(data_root, "classes.tsv"), col_types = cols(.default = "c", mixed = "l"))

read_mat <- function(f) {
  m <- read.delim(gzfile(f), row.names = 1, check.names = FALSE)
  as.matrix(m)
}
diffs <- read_mat(file.path(data_root, "pair_differences.tsv.gz"))
sites <- read_mat(file.path(data_root, "pair_sites.tsv.gz"))

### p-distance per pair of representatives (squared Euclidean distance of the AMOVA)
p_dist <- diffs / sites
diag(p_dist) <- 0

### representatives per sample set
sets <- list(all = classes$sample, no_mixed = classes$sample[!classes$mixed])



# checks -----------------------------------------------------------------------

## check if the matrices hold the representatives of classes.tsv in the same order
if (!identical(rownames(p_dist), classes$sample) || !identical(colnames(p_dist), classes$sample)) {
  stop("distance matrices and classes.tsv differ in samples or order")
}

## check if every pair has compared sites
if (any(sites[upper.tri(sites)] == 0)) {
  stop("pairs of representatives without compared sites")
}



# AMOVA ------------------------------------------------------------------------

### units as groups, sites within units as populations, representatives as haplotypes
run_amova <- function(kk, ss) {
  keep <- classes$sample %in% sets[[ss]]
  cl   <- classes[keep, ]
  pop  <- paste(cl[[kk]], cl$site, sep = "|")

  samples    <- as.data.frame.matrix(table(cl$sample, pop)[cl$sample, ])
  structures <- data.frame(unit = factor(str_remove(colnames(samples), "\\|.*$")))

  d <- as.dist(sqrt(p_dist[keep, keep]))
  euclid <- is.euclid(d)
  if (!euclid) d <- quasieuclid(d)

  am <- amova(samples, d, structures)
  rt <- randtest(am, nrepet = n_perm)

  comp <- am$componentsofcovariance
  tibble(k = kk, set = ss, n = nrow(cl), units = n_distinct(cl[[kk]]), sites = n_distinct(pop),
         source = c("among units", "among sites within units", "within sites"),
         df     = am$results$Df[1:3],
         sigma2 = comp$Sigma[1:3],
         pct    = comp$`%`[1:3],
         phi    = rev(am$statphi$Phi),
         phi_name = c("Phi_CT", "Phi_SC", "Phi_ST"),
         p      = rt$pvalue[c(3, 2, 1)],
         euclidean = euclid)
}

set.seed(seed)
res <- expand_grid(k = k_levels, set = names(sets)) %>%
  mutate(out = map2(k, set, run_amova)) %>%
  select(out) %>%
  unnest(out)



# figures ----------------------------------------------------------------------

### share of the molecular variance per level (all representatives)
p_amova <- res %>%
  filter(set == "all") %>%
  mutate(Level  = factor(sub("k", "K = ", k), levels = rev(sub("k", "K = ", k_levels))),
         source = factor(source, levels = rev(names(lab_source)))) %>%
  ggplot(aes(x = pct, y = Level, fill = source)) +
  geom_col(width = 0.6, colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = if_else(pct >= 3, sprintf("%.0f%%", pct), ""), colour = source == "among units"),
            position = position_stack(vjust = 0.5), size = 7 / .pt) +
  scale_fill_manual(values = pal_source, labels = lab_source, breaks = names(lab_source), name = NULL) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
  labs(x = "Molecular variance (%)", y = NULL)

save_fig("amova.pdf", p_amova, 3.4, 2.2)



# tables -----------------------------------------------------------------------

### variance components per level (all representatives and without mixed-MAT samples)
pub <- res %>%
  filter(set == "all") %>%
  left_join(res %>% filter(set == "no_mixed") %>% select(k, source, pct_nm = pct, phi_nm = phi),
            by = c("k", "source")) %>%
  transmute(Level = sub("k", "K = ", k), Source = source, df = df,
            `Variance` = signif(sigma2, 3), `% of total` = sprintf("%.1f", pct),
            Statistic = phi_name, Phi = sprintf("%.3f", phi),
            p = if_else(p < 0.001, "< 0.001", sprintf("%.3f", p)),
            `Without mixed MAT: % (Phi)` = sprintf("%.1f (%.3f)", pct_nm, phi_nm))

write_tsv(res, file.path(out_root, "amova.tsv"))
write_tsv(pub, file.path(out_root, "publication_amova.tsv"))



# summary ----------------------------------------------------------------------

cat("\n=== AMOVA ===\n")
res %>% mutate(across(c(sigma2, pct, phi), ~ signif(.x, 3))) %>%
  select(k, set, source, df, sigma2, pct, phi_name, phi, p) %>%
  print(n = Inf, width = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  amova.pdf\n")
cat("  amova.tsv\n")
cat("  publication_amova.tsv\n")
