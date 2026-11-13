## ============================================================================
## 060_delimitation_summary.R
## Lines of evidence per unit at K = 4, 5 and 6 in one table, and the units used from
## here on (popgen_units.tsv)
##
## Files needed:
##   lineage_labels.tsv             units per sample at K = 4, 5 and 6
##   unit_sizes.tsv                 DAPC cluster number of each unit
##                                  (both from local/03_delimitation/001_lineage_labels.R)
##   sample_metadata_ALL_FIXED.tsv  mating type per sample,
##                                  from local/01_variant_calling/103_mating_type.R
##   aster_units.tsv                sister unit, quartet support and terminal length in coalescent units,
##                                  from local/03_delimitation/010_aster_lineages.R
##   unit_distances.tsv             mean 1 - IBS within and between units,
##                                  from local/03_delimitation/020_neighbournet.R
##   phi_summary_sets.tsv           share of genes with PHI evidence for recombination,
##                                  from local/03_delimitation/021_phi_test.R
##   index_association_sets.tsv     rbarD per unit and sample set,
##                                  from local/03_delimitation/022_index_association.R
##   dsuite_trios.tsv               D and f4-ratio per trio
##   dsuite_fbranch.tsv             f-branch per branch and P3 unit
##                                  (both from local/03_delimitation/030_dsuite.R)
##   cluster_monophyly.tsv          monophyly and stem support of the 014 clusters in the four trees,
##                                  from local/02_structure_phylogeny/100_cluster_monophylly.R
##   mptp_units.tsv                 mPTP species per unit and MCMC support,
##                                  from local/03_delimitation/040_mptp.R
##   bpp_gdi.tsv                    gdi per unit and locus set
##   publication_bpp_a11.tsv        A11 posterior per delimited species
##                                  (both from local/03_delimitation/051_bpp.R)
##   depth_profiles.tsv             depth profile test per unit, from the soil niche script 040
##                                  of local/04_diversity_niche_selection (results/040_results_soil_niche;
##                                  a later phase, run it before this script)
## ============================================================================

library(tidyverse)



# files ------------------------------------------------------------------------

res_03   <- "03_delimitation/results"
out_root <- file.path(res_03, "060_results_delimitation_summary")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

labels_file <- file.path(res_03, "001_results_lineage_labels", "lineage_labels.tsv")
key_file    <- file.path(res_03, "001_results_lineage_labels", "unit_sizes.tsv")
meta_file   <- "01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv"

### one entry per line of evidence
f_in <- c(aster  = file.path(res_03, "010_results_aster_lineages", "aster_units.tsv"),
          nnet   = file.path(res_03, "020_results_neighbournet", "unit_distances.tsv"),
          phi    = file.path(res_03, "021_results_phi_test", "phi_summary_sets.tsv"),
          ia     = file.path(res_03, "022_results_index_association", "index_association_sets.tsv"),
          trios  = file.path(res_03, "030_results_dsuite", "dsuite_trios.tsv"),
          fb     = file.path(res_03, "030_results_dsuite", "dsuite_fbranch.tsv"),
          trees  = "02_structure_phylogeny/results/100_results_cluster_monophyly/cluster_monophyly.tsv",
          mptp   = file.path(res_03, "040_results_mptp", "mptp_units.tsv"),
          gdi    = file.path(res_03, "051_results_bpp", "bpp_gdi.tsv"),
          a11    = file.path(res_03, "051_results_bpp", "publication_bpp_a11.tsv"),
          depth  = "04_diversity_niche_selection/results/040_results_soil_niche/depth_profiles.tsv")



# parameters -------------------------------------------------------------------

### levels
k_levels <- c("k4", "k5", "k6")

### decision: level of the units used from here on, and units merged at that level (empty: none)
final_level <- "k6"
merge_units <- list()

### significance level (Dsuite q values)
alpha <- 0.05

### gdi thresholds (one lineage below, separate species above)
gdi_low  <- 0.2
gdi_high <- 0.7



# read -------------------------------------------------------------------------

labels <- read_tsv(labels_file, col_types = cols(.default = "c"))
meta   <- read_tsv(meta_file, col_types = cols(.default = "c"))
reps   <- filter(labels, clustered == "TRUE") %>%
  left_join(select(meta, sample, mat_genome), by = "sample")

### lines of evidence
inp <- map(f_in, ~ read_tsv(.x, show_col_types = FALSE, guess_max = 10000))

### levels written as k4, k5, k6 everywhere
norm_k <- function(d) if ("k" %in% names(d)) mutate(d, k = paste0("k", parse_number(as.character(k)))) else d
inp    <- map(inp, norm_k)



# checks -----------------------------------------------------------------------

## check if every representative has a mating type call
if (anyNA(reps$mat_genome)) warning(sum(is.na(reps$mat_genome)), " representatives without a mating type call")

## check if the final level is one of the levels
if (!final_level %in% k_levels) stop("final_level must be one of ", paste(k_levels, collapse = ", "))

## check if merged units exist at the final level
bad <- setdiff(unlist(merge_units), unique(labels[[final_level]]))

if (length(bad) > 0) stop("units to merge not found at ", final_level, ": ", paste(bad, collapse = ", "))



# units ------------------------------------------------------------------------

### size, spread, project share and mating types of the representatives
units <- map_dfr(k_levels, function(kk) {
  reps %>%
    mutate(unit = .data[[kk]]) %>%
    summarise(n_mlg     = n(),
              n_sites   = n_distinct(site),
              n_country = n_distinct(country),
              pct_magic = 100 * mean(project == "magic"),
              mat1      = sum(mat_genome == "MAT1", na.rm = TRUE),
              mat2      = sum(mat_genome == "MAT2", na.rm = TRUE),
              .by = unit) %>%
    left_join(count(labels, unit = .data[[kk]], name = "n_samples"), by = "unit") %>%
    mutate(k = kk, both_mat = mat1 > 0 & mat2 > 0, .before = 1)
})



# lines of evidence ------------------------------------------------------------

ev <- list(units)

### 010: sister unit, quartet support of the clade with the sister, terminal length (coalescent units)
ev$aster <- select(inp$aster, k, unit, sister, q1_sister, p_sister, cu_terminal)

### 020: mean distance within the unit and to the closest other unit
d <- inp$nnet
within  <- d %>% filter(within) %>% transmute(k, unit = unit_1, d_within = mean_d)
between <- d %>%
  filter(!within) %>%
  (\(x) bind_rows(transmute(x, k, unit = unit_1, other = unit_2, mean_d),
                  transmute(x, k, unit = unit_2, other = unit_1, mean_d)))() %>%
  slice_min(mean_d, n = 1, by = c(k, unit), with_ties = FALSE) %>%
  rename(nearest = other, d_nearest = mean_d)
ev$nnet <- left_join(within, between, by = c("k", "unit")) %>% mutate(d_ratio = d_nearest / d_within)

### 021: share of genes with recombination signal (all, equal size, without mixed-MAT samples)
ev$phi <- inp$phi %>%
  summarise(pct = mean(pct_sig), p = if (n() == 1) p_binom else NA_real_, .by = c(group, unit)) %>%
  pivot_wider(names_from = group, values_from = c(pct, p), names_glue = "phi_{.value}_{group}") %>%
  select(unit, any_of(c("phi_pct_all", "phi_p_all", "phi_pct_equal_n", "phi_pct_no_mixed")))

### 022: standardised index of association (all, equal size)
ev$ia <- inp$ia %>%
  summarise(rbarD = mean(rbarD), p = if (n() == 1) p_rD else NA_real_, .by = c(group, unit)) %>%
  pivot_wider(names_from = group, values_from = c(rbarD, p), names_glue = "ia_{.value}_{group}") %>%
  select(unit, any_of(c("ia_rbarD_all", "ia_p_all", "ia_rbarD_equal_n")))

### 030: significant trios involving the unit and the largest f-branch into or out of it (all samples)
t <- filter(inp$trios, set == "all")
ev$trios <- map_dfr(k_levels, function(kk) {
  tk <- filter(t, k == kk)
  tibble(unit = unique(c(tk$P1, tk$P2, tk$P3))) %>%
    mutate(k = kk,
           n_trios     = map_int(unit, ~ sum(tk$P1 == .x | tk$P2 == .x | tk$P3 == .x)),
           n_trios_sig = map_int(unit, ~ sum((tk$P1 == .x | tk$P2 == .x | tk$P3 == .x) & tk$q < alpha)),
           f4_max_sig  = map_dbl(unit, ~ max(c(0, tk$f4[(tk$P2 == .x | tk$P3 == .x) & tk$q < alpha]))))
})
f <- filter(inp$fb, set == "all", !is.na(fb))
ev$fb <- full_join(f %>% rename(unit = branch_descendants) %>% summarise(fb_into = max(fb), .by = c(k, unit)),
                   f %>% rename(unit = P3) %>% summarise(fb_from = max(fb), .by = c(k, unit)), by = c("k", "unit"))

### 100 (02): monophyly in the four trees and stem support, DAPC cluster numbers translated to units
key <- read_tsv(key_file, col_types = cols(.default = "c")) %>%
  transmute(k = paste0("k", k), dapc, unit, n_mlg = as.integer(n_mlg))
tr <- inp$trees %>%
  transmute(tree, k = paste0("k", K), dapc = as.character(cluster), n_tips, monophyletic,
            across(c(support, sCF, gCF), as.numeric)) %>%
  inner_join(key, by = c("k", "dapc"))

## check if the clusters of 100 are the units of 001 (same number of representatives)
if (any(tr$n_tips != tr$n_mlg)) warning("cluster sizes in 100 differ from the unit sizes in 001")

### stem support of monophyletic units (concatenation UFBoot, sCF, gCF; CASTER local bootstrap; ASTER local posterior)
ev$trees <- tr %>%
  summarise(mono = sprintf("%d/%d", sum(monophyletic), n()), .by = c(k, unit)) %>%
  left_join(tr %>% filter(tree == "concat") %>%
              transmute(k, unit, ufboot = if_else(monophyletic, support, NA), scf = if_else(monophyletic, sCF, NA),
                        gcf = if_else(monophyletic, gCF, NA)), by = c("k", "unit")) %>%
  left_join(tr %>% filter(tree == "caster") %>% transmute(k, unit, lbs = if_else(monophyletic, support, NA)),
            by = c("k", "unit")) %>%
  left_join(tr %>% filter(tree == "aster") %>% transmute(k, unit, lpp = if_else(monophyletic, support, NA)),
            by = c("k", "unit"))

### 040: mPTP species within the unit and shared with others, support as one species and as separate from the sister
ev$mptp <- transmute(inp$mptp, k, unit, mptp_species = species,
                     mptp_shared = shared,
                     mptp_merged = merged, mptp_p_one = p_one_species, mptp_p_separate = p_separate)

### 051: gdi per locus set and the call (thresholds met in all sets)
ev$gdi <- inp$gdi %>%
  summarise(gdi_min = min(gdi_mean), gdi_max = max(gdi_mean),
            gdi_call = case_when(all(gdi_mean > gdi_high) ~ "separate",
                                 all(gdi_mean < gdi_low)  ~ "one lineage",
                                 TRUE                     ~ "ambiguous"),
            .by = c(k, unit))

### 051: lowest A11 posterior of the delimited unit over runs (finest level)
ev$a11 <- inp$a11 %>%
  select(species = `P per delimited species`) %>%
  separate_longer_delim(species, ", ") %>%
  separate_wider_delim(species, " ", names = c("unit", "pp")) %>%
  summarise(a11_pp = min(as.numeric(pp)), .by = unit) %>%
  mutate(k = last(k_levels))

### 040 (04): depth profile test (PhytOakmeter, all samples)
ev$depth <- inp$depth %>%
  filter(set == "all") %>%
  summarise(depth_q = min(q), .by = c(k, unit))

### one row per level and unit
join_ev <- function(x, y) left_join(x, y, by = intersect(c("k", "unit"), names(y)))
evidence <- reduce(ev, join_ev) %>%
  arrange(k, unit)



# final units ------------------------------------------------------------------

### units at the final level, merged where decided
merge_key <- unlist(unname(imap(merge_units, ~ setNames(rep(.y, length(.x)), .x))))
popgen <- labels %>%
  transmute(sample, mlg, clustered, unit = .data[[final_level]])
if (length(merge_key) > 0) {
  popgen <- mutate(popgen, unit = if_else(unit %in% names(merge_key), unname(merge_key[unit]), unit))
}



# tables -----------------------------------------------------------------------

### fixed decimals, "-" for NA (numeric first: formatC does not take an all-NA logical column)
fmt <- function(x, d) {
  x <- as.numeric(x)
  if_else(is.na(x), "-", formatC(x, format = "f", digits = d))
}

pub <- evidence %>%
  transmute(Level                            = sub("k", "K = ", k),
            Unit                             = unit,
            `MLGs (sites)`                   = sprintf("%d (%d)", n_mlg, n_sites),
            `MAGIC (%)`                      = fmt(pct_magic, 0),
            `MAT1:MAT2`                      = sprintf("%d:%d", mat1, mat2),
            Sister                           = sister,
            `Monophyletic (trees)`           = mono,
            `UFBoot/sCF/gCF`                 = sprintf("%s/%s/%s", fmt(ufboot, 0), fmt(scf, 0), fmt(gcf, 0)),
            `CASTER LBS`                     = fmt(lbs, 0),
            `ASTER LPP`                      = fmt(lpp, 2),
            `Recombining genes (%)`          = fmt(phi_pct_all, 1),
            rbarD                            = fmt(ia_rbarD_all, 3),
            `Max f-branch`                   = fmt(pmax(fb_into, fb_from, na.rm = TRUE), 3),
            `mPTP species (within / shared)` = sprintf("%d / %s", mptp_species,
                                                       if_else(is.na(mptp_shared), "-", as.character(mptp_shared))),
            gdi                              = sprintf("%s-%s (%s)", fmt(gdi_min, 2), fmt(gdi_max, 2), gdi_call),
            `A11 PP`                         = fmt(a11_pp, 2)) %>%
  select(where(~ !all(is.na(.x))))

write_tsv(evidence, file.path(out_root, "delimitation_evidence.tsv"))
write_tsv(pub, file.path(out_root, "publication_delimitation_evidence.tsv"))
write_tsv(popgen, file.path(out_root, "popgen_units.tsv"))



# summary ----------------------------------------------------------------------

cat("\n=== Evidence per level and unit ===\n")
print(pub, n = Inf, width = Inf)

cat("\n=== Final units (", final_level, if (length(merge_units) > 0) ", merged" else "", ") ===\n", sep = "")
popgen %>%
  summarise(samples = n(), mlgs = sum(clustered == "TRUE"), .by = unit) %>%
  arrange(unit) %>%
  print(n = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  delimitation_evidence.tsv\n")
cat("  publication_delimitation_evidence.tsv\n")
cat("  popgen_units.tsv\n")
