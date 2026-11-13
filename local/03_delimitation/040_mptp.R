## ============================================================================
## 040_mptp.R
## mPTP delimitation of the MLG representatives against the units at K = 4, 5 and 6:
## species number, support per unit and sensitivity to the minimum branch length
##
## Files needed:
##   mptp_ml_minbrF.txt           ML delimitation per minimum branch length factor
##   mptp_ml_sensitivity.tsv      species number per minimum branch length
##   mptp_mcmc.<seed>.stats       posterior of the species number per MCMC run
##   mptp_mcmc.42.combined.tree   probability that a node lies in the speciation process,
##                                runs combined
##                                (all four from hpc/03_delimitation/040_mptp.sh;
##                                 copy them from ${ROOT}/03_delimitation/results/040_results_mptp
##                                 to local/03_delimitation/data/040_data_mptp)
##   lineage_labels.tsv           units per sample at K = 4, 5 and 6,
##                                from local/03_delimitation/001_lineage_labels.R
## ============================================================================

library(ape)
library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

data_root   <- "03_delimitation/data/040_data_mptp"
out_root    <- "03_delimitation/results/040_results_mptp"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

labels_file <- "03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv"



# parameters -------------------------------------------------------------------

### levels
k_levels <- c("k4", "k5", "k6")

### seed of the MCMC call (names the combined files)
seed <- 42



# style ------------------------------------------------------------------------

### unit colours
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



# functions --------------------------------------------------------------------

### species blocks of an mPTP delimitation file
read_mptp <- function(path) {
  l <- read_lines(path)
  i <- str_which(l, "^Species [0-9]+:")
  map_dfr(seq_along(i), function(j) {
    e    <- if (j < length(i)) i[j + 1] - 1 else length(l)
    tips <- str_trim(l[(i[j] + 1):e])
    tibble(sample = tips[tips != ""], species = j)
  })
}

### posterior of the species number and its 95% credible interval, one MCMC run
read_stats <- function(path) {
  l <- read_lines(path)
  d <- str_subset(l, "^[0-9]+,")
  m <- str_match(str_subset(l, "^CCI"), "\\((\\d+),(\\d+)\\)")
  tibble(run     = basename(path),
         species = as.numeric(str_extract(d, "^[0-9]+")),
         pct     = as.numeric(str_extract(d, "(?<=,)[0-9.eE+-]+")),
         cci_lo  = as.numeric(m[1, 2]),
         cci_hi  = as.numeric(m[1, 3]))
}

### adjusted Rand index
ari <- function(x, y) {
  tab <- table(x, y)
  s   <- sum(choose(tab, 2))
  a   <- sum(choose(rowSums(tab), 2))
  b   <- sum(choose(colSums(tab), 2))
  e   <- a * b / choose(sum(tab), 2)
  (s - e) / ((a + b) / 2 - e)
}



# read -------------------------------------------------------------------------

labels <- read_tsv(labels_file, col_types = cols(.default = "c"))
reps   <- filter(labels, clustered == "TRUE")

sens <- read_tsv(file.path(data_root, "mptp_ml_sensitivity.tsv"), col_types = "dddd")

ml <- map_dfr(sens$factor, ~ read_mptp(file.path(data_root, sprintf("mptp_ml_minbr%g.txt", .x))) %>%
                mutate(factor = .x, .before = 1))

f_stats <- list.files(data_root, pattern = "^mptp_mcmc\\.[0-9]+\\.stats$", full.names = TRUE)
f_stats <- f_stats[!str_detect(basename(f_stats), sprintf("\\.%d\\.stats$", seed))]
post    <- map_dfr(f_stats, read_stats)
runs    <- distinct(post, run, cci_lo, cci_hi)

support_tree <- read.tree(file.path(data_root, sprintf("mptp_mcmc.%d.combined.tree", seed)))

### main minimum branch length (first factor)
factor_main <- first(sens$factor)
ml_main     <- filter(ml, factor == factor_main)



# checks -----------------------------------------------------------------------

## check if every representative was delimited
miss <- setdiff(reps$sample, ml_main$sample)

if (length(miss) > 0) warning(length(miss), " representatives missing from the ML delimitation")

## check if the delimitation holds only representatives (outgroup cropped)
extra <- setdiff(ml_main$sample, reps$sample)

if (length(extra) > 0) warning(length(extra), " delimited tips are not representatives: ",
                               paste(head(extra), collapse = ", "))

## check if there are MCMC runs to compare
if (length(f_stats) < 2) warning("fewer than two MCMC runs found")



# agreement with the units -----------------------------------------------------

### species number and adjusted Rand index per minimum branch length and level
levels_tbl <- expand_grid(factor = sens$factor, k = k_levels) %>%
  mutate(res = map2(factor, k, function(f, kk) {
    d <- filter(ml, factor == f) %>% inner_join(select(reps, sample, unit = all_of(kk)), by = "sample")
    tibble(units = n_distinct(d$unit), species = n_distinct(d$species), ari = ari(d$unit, d$species))
  })) %>%
  unnest(res) %>%
  left_join(select(sens, factor, minbr), by = "factor")

### mPTP species per unit at the main minimum branch length, species shared with other units
units_tbl <- map_dfr(k_levels, function(kk) {
  d <- ml_main %>% inner_join(select(reps, sample, unit = all_of(kk)), by = "sample")
  shared <- d %>% distinct(species, unit) %>% count(species, name = "units_in_species")
  d %>%
    left_join(shared, by = "species") %>%
    summarise(n = n(), shared = n_distinct(species[units_in_species > 1]), species = n_distinct(species),
              merged = any(units_in_species > 1), .by = unit) %>%
    relocate(shared, .after = species) %>%
    mutate(k = kk, .before = 1)
})

### MCMC support per unit (one species: 1 - support at its MRCA; separate from the sister: support at its parent)
node_sup <- function(tr, node) suppressWarnings(as.numeric(tr$node.label[node - Ntip(tr)]))

mcmc_tbl <- map_dfr(k_levels, function(kk) {
  key <- setNames(reps[[kk]], reps$sample)
  map_dfr(sort(unique(key)), function(u) {
    m    <- intersect(names(key)[key == u], support_tree$tip.label)
    mrca <- getMRCA(support_tree, m)
    mono <- setequal(intersect(extract.clade(support_tree, mrca)$tip.label, names(key)), m)
    par  <- support_tree$edge[support_tree$edge[, 2] == mrca, 1]
    tibble(k = kk, unit = u, monophyletic = mono,
           p_one_species = if (mono) 1 - node_sup(support_tree, mrca) else NA_real_,
           p_separate    = if (mono && length(par) == 1) node_sup(support_tree, par) else NA_real_)
  })
})

units_tbl <- left_join(units_tbl, mcmc_tbl, by = c("k", "unit"))

### posterior of the species number, mean and range over runs
post_sum <- post %>%
  select(run, species, pct) %>%
  complete(run, species, fill = list(pct = 0)) %>%
  summarise(lo = min(pct), hi = max(pct), pct = mean(pct), .by = species) %>%
  filter(hi > 0) %>%
  arrange(species)



# figures ----------------------------------------------------------------------

n_ml    <- n_distinct(ml_main$species)
n_units <- levels_tbl %>% distinct(k, units)

p_post <- ggplot(post_sum, aes(x = species, y = pct)) +
  geom_col(fill = "grey60", width = 0.8) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.3) +
  geom_vline(xintercept = n_ml, linetype = "dashed", linewidth = 0.4) +
  labs(x = "Number of species (mPTP)", y = "Posterior (%)")

p_units <- units_tbl %>%
  mutate(level = sub("k", "K = ", k)) %>%
  ggplot(aes(x = unit, y = species, fill = unit)) +
  geom_col(colour = "grey15", linewidth = 0.2, width = 0.7) +
  facet_grid(. ~ level, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal_unit, guide = "none") +
  labs(x = NULL, y = "mPTP species within the unit")

save_fig("mptp.pdf",
         (p_post | p_units) + plot_layout(widths = c(1, 1.6)) + plot_annotation(tag_levels = "A"),
         7.1, 2.8)



# tables -----------------------------------------------------------------------

pub_units <- units_tbl %>%
  transmute(Level = sub("k", "K = ", k), Unit = unit, Representatives = n, `mPTP species` = species,
            `mPTP species shared with other units` = shared,
            `P(one species)` = sprintf("%.2f", p_one_species),
            `P(separate from sister)` = sprintf("%.2f", p_separate))

pub_levels <- levels_tbl %>%
  transmute(`Minimum branch length` = sprintf("%.2e", minbr), Level = sub("k", "K = ", k), Units = units,
            `mPTP species` = species, ARI = sprintf("%.3f", ari))

write_tsv(ml_main, file.path(out_root, "mptp_species.tsv"))
write_tsv(units_tbl, file.path(out_root, "mptp_units.tsv"))
write_tsv(levels_tbl, file.path(out_root, "mptp_levels.tsv"))
write_tsv(post_sum, file.path(out_root, "mptp_species_posterior.tsv"))
write_tsv(runs, file.path(out_root, "mptp_mcmc_runs.tsv"))
write_tsv(pub_units, file.path(out_root, "publication_mptp_units.tsv"))
write_tsv(pub_levels, file.path(out_root, "publication_mptp_levels.tsv"))



# summary ----------------------------------------------------------------------

cat("\n=== Species number and agreement with the units, per minimum branch length ===\n")
print(pub_levels, n = Inf, width = Inf)

cat("\n=== Units at the main minimum branch length ===\n")
print(pub_units, n = Inf, width = Inf)

cat("\n=== 95% credible interval of the species number per MCMC run ===\n")
print(runs, n = Inf)

cat("\n=== Posterior of the species number (mean over", length(f_stats), "runs) ===\n")
print(post_sum, n = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  mptp.pdf\n")
cat("  mptp_species.tsv\n")
cat("  mptp_units.tsv\n")
cat("  mptp_levels.tsv\n")
cat("  mptp_species_posterior.tsv\n")
cat("  mptp_mcmc_runs.tsv\n")
cat("  publication_mptp_units.tsv\n")
cat("  publication_mptp_levels.tsv\n")
