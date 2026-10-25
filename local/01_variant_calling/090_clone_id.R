## ============================================================================
## 090_clone_id.R
## Clone detection: multilocus genotypes from the plink distance matrix
##
## Files needed:
##   plink.genome                   pairwise IBS counts
##   plink.imiss                    per-sample missingness
##                                  (both from hpc/01_variant_calling/090_clone_id.sh;
##                                   copy them from ${ROOT}/01_variant_calling/results/090_results_clone_id
##                                   to local/01_variant_calling/data/090_data_clone_id)
##
## What each choice rests on:
##   plink IBS counts               Chang et al. 2015, GigaScience 4:7
##   poppr                          Kamvar et al. 2014, PeerJ 2:e281
##   filter_stats, cutoff_predictor Kamvar et al. 2015, Front Genet 6:208
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)
library(poppr)



# files ------------------------------------------------------------------------

data_090 <- "01_variant_calling/data/090_data_clone_id"
out_root <- "01_variant_calling/results/090_results_clone_id"

f_genome <- file.path(data_090, "plink.genome")
f_imiss  <- file.path(data_090, "plink.imiss")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### clone threshold (NA: farthest neighbour prediction of cutoff_predictor)
clone_dist <- NA

### mlg.filter algorithms (farthest neighbour: default and strictest)
algos <- c(farthest = "complete", nearest = "single", average = "average")

### tolerance of the genotype-count plateau around clone_dist
tol <- 0.01



# style ------------------------------------------------------------------------

### Okabe-Ito colours, neutral reference lines
col_algo   <- c(farthest = "#D55E00", nearest = "#0072B2", average = "#009E73")
col_status <- c(Unique = "grey70", Representative = "#0072B2",
                Clonemate = "#D55E00")
col_cut    <- "grey20"
col_fill   <- "grey60"

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

### pdf device (fonts embedded)
save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}

lab_log <- scales::label_log()



# read -------------------------------------------------------------------------

### pairwise counts (IBS0 differing, IBS2 matching, IBS1 zero for haploid
### calls; denominator per pair)
pairs <- read_table(f_genome, show_col_types = FALSE) %>%
  transmute(sample_1 = IID1, sample_2 = IID2,
            n_diff = IBS0, n_het = IBS1, n_same = IBS2,
            n_compared = IBS0 + IBS1 + IBS2,
            dist = n_diff / n_compared)

### per-sample missingness, country and site (1-xxx samples: Swiss, no site code)
stats <- read_table(f_imiss, show_col_types = FALSE) %>%
  transmute(sample = IID, n_missing = N_MISS, f_missing = F_MISS) %>%
  mutate(country = if_else(str_detect(sample, "^1-"), "CH",
                           str_extract(sample, "^[A-Z]{2}")),
         site    = replace_na(str_extract(sample, "^[A-Z]{2}-[A-Z]{3}"), "other"))

site_of <- setNames(stats$site, stats$sample)

message("pairs: ", nrow(pairs), "   samples: ", nrow(stats),
        "   IBS1 > 0 in ", sum(pairs$n_het > 0), " pairs")
message("distance: ", signif(min(pairs$dist), 3), " - ",
        signif(max(pairs$dist), 3),
        "   sites compared: ", min(pairs$n_compared), " - ",
        max(pairs$n_compared))



# distance matrix and trees ----------------------------------------------------

### distance matrix from the exact counts (pairwise complete)
ids <- sort(unique(c(pairs$sample_1, pairs$sample_2)))

dmat <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
dmat[cbind(pairs$sample_1, pairs$sample_2)] <- pairs$dist
dmat[cbind(pairs$sample_2, pairs$sample_1)] <- pairs$dist

d     <- as.dist(dmat)
trees <- map(algos, ~ hclust(d, method = .x))

### merge heights (the THRESHOLDS of poppr filter_stats)
heights <- map(trees, ~ sort(.x$height))



# threshold --------------------------------------------------------------------

predicted <- tibble(algorithm = names(algos),
                    predicted = map_dbl(heights, cutoff_predictor))

if (is.na(clone_dist)) clone_dist <- predicted$predicted[1]

write_tsv(predicted, file.path(out_root, "predicted_thresholds.tsv"))

### clone_dist is mirrored in hpc/01_variant_calling/130_snp_distribution.sh
### (THRESHOLD) and local/01_variant_calling/130_snp_density.R
cat("\n=== predicted thresholds (poppr cutoff_predictor) ===\n")
print(predicted, width = Inf)
cat(sprintf("\nclone_dist %.3e\n", clone_dist))
cat(sprintf("pairs below it: %d of %d\n",
            sum(pairs$dist <= clone_dist), nrow(pairs)))



# genotypes against threshold --------------------------------------------------

### genotype count after each merge (poppr filter_stats curve)
collapse <- imap_dfr(heights, ~ tibble(algorithm = .y,
                                       threshold = c(0, .x),
                                       n_mlg = length(ids) - seq(0, length(.x))))

write_tsv(collapse, file.path(out_root, "filter_stats.tsv"))



# sampling-design bound --------------------------------------------------------

### genotypes spanning two countries (upper bound on the threshold)
country_of <- setNames(stats$country, stats$sample)

### grid over the observed distance range
x_range <- range(pairs$dist)
x_grid  <- 10^seq(log10(x_range[1]), log10(x_range[2]), length.out = 400)

cross_country <- map_dfr(x_grid, function(h) {
  tibble(sample = trees$farthest$labels,
         grp    = cutree(trees$farthest, h = h)) %>%
    mutate(country = country_of[sample]) %>%
    group_by(grp) %>%
    summarise(n_countries = n_distinct(country), .groups = "drop") %>%
    summarise(threshold = h, n_cross = sum(n_countries > 1))
})

write_tsv(cross_country, file.path(out_root, "cross_country_bound.tsv"))

cross_bound <- cross_country %>% filter(n_cross > 0) %>% slice(1) %>%
  pull(threshold)

n_final <- n_distinct(cutree(trees$farthest, h = clone_dist))

### last grid point without a cross-country genotype
cross_below <- cross_country %>% filter(threshold < cross_bound) %>%
  slice_tail(n = 1) %>% pull(threshold)

### plateau (genotype count within tol of the count at clone_dist)
in_band <- collapse %>%
  filter(algorithm == "farthest", threshold > 0,
         abs(n_mlg - n_final) <= tol * n_final)

band_lo <- min(in_band$threshold)
band_hi <- max(in_band$threshold)

cat(sprintf("\ncross-country bound                %.3e to %.3e\n",
            cross_below, cross_bound))
cat(sprintf("  relative to clone_dist           %.2f to %.2f\n",
            cross_below / clone_dist, cross_bound / clone_dist))
cat(sprintf("plateau (count within %.0f%% of %d)  %.3e to %.3e (%.0f-fold)\n",
            tol * 100, n_final, band_lo, band_hi, band_hi / band_lo))
cat(sprintf("  genotype count on the plateau    %d to %d\n",
            min(in_band$n_mlg), max(in_band$n_mlg)))



# multilocus genotypes ---------------------------------------------------------

### genotypes at clone_dist (cutree(h = x) equals mlg.filter at
### x + .Machine$double.eps^0.5)
assign <- tibble(sample   = trees$farthest$labels,
                 farthest = cutree(trees$farthest, h = clone_dist),
                 nearest  = cutree(trees$nearest,  h = clone_dist),
                 average  = cutree(trees$average,  h = clone_dist))

n_joint <- n_distinct(paste(assign$farthest, assign$nearest))
agree   <- n_joint == n_distinct(assign$farthest) &&
  n_joint == n_distinct(assign$nearest)

cat("\n=== multilocus genotypes at", signif(clone_dist, 4), "===\n")
cat("farthest neighbour:", n_distinct(assign$farthest), "\n")
cat("nearest neighbour: ", n_distinct(assign$nearest), "\n")
cat("average neighbour: ", n_distinct(assign$average), "\n")
cat("farthest and nearest agree:", agree, "\n")

mlg <- assign %>%
  select(sample, mlg = farthest) %>%
  left_join(stats, by = "sample") %>%
  group_by(mlg) %>%
  mutate(n_members = n()) %>%
  ungroup()



# clone groups -----------------------------------------------------------------

### worst pair per genotype (most differing sites, fewest compared sites)
within <- pairs %>%
  left_join(select(mlg, sample, mlg_1 = mlg), by = c("sample_1" = "sample")) %>%
  left_join(select(mlg, sample, mlg_2 = mlg), by = c("sample_2" = "sample")) %>%
  filter(mlg_1 == mlg_2) %>%
  group_by(mlg = mlg_1) %>%
  summarise(n_pairs = n(), max_n_diff = max(n_diff),
            min_n_compared = min(n_compared), .groups = "drop")

groups <- mlg %>%
  filter(n_members > 1) %>%
  group_by(mlg, n_members) %>%
  summarise(members     = paste(sort(sample), collapse = ","),
            n_sites     = n_distinct(site),
            n_countries = n_distinct(country),
            countries   = paste(sort(unique(country)), collapse = ", "),
            .groups = "drop") %>%
  left_join(within, by = "mlg") %>%
  arrange(desc(n_members))

write_tsv(groups, file.path(out_root, "clone_groups.tsv"))

cat("\n=== clone groups ===\n")
groups %>%
  select(mlg, n_members, n_sites, n_countries, max_n_diff, min_n_compared,
         countries) %>%
  print(n = Inf)



# representatives --------------------------------------------------------------

### representative: the member with the lowest missingness
mlg <- mlg %>%
  group_by(mlg) %>%
  arrange(f_missing, .by_group = TRUE) %>%
  mutate(is_representative = row_number() == 1,
         status = case_when(n_members == 1    ~ "Unique",
                            is_representative ~ "Representative",
                            TRUE              ~ "Clonemate"),
         status = factor(status, levels = names(col_status))) %>%
  ungroup()

write_tsv(mlg, file.path(out_root, "per_sample_mlg.tsv"))

write_lines(filter(mlg, status == "Clonemate")$sample,
            file.path(out_root, "samples_excluded_clones.txt"))
write_lines(filter(mlg, status != "Clonemate")$sample,
            file.path(out_root, "samples_retained_clones.txt"))

cat("\n=== samples ===\n")
cat("total:                ", nrow(mlg), "\n")
cat("multilocus genotypes: ", n_distinct(mlg$mlg), "\n")
cat("clone groups:         ", nrow(groups), "\n")
cat("in a clone group:     ", sum(mlg$n_members > 1), "\n")
cat("clonemates to remove: ", sum(mlg$status == "Clonemate"), "\n")
cat("retained:             ", sum(mlg$status != "Clonemate"), "\n")



# sites compared per clone group -----------------------------------------------

cat("\n=== sites compared, worst pair per clone group ===\n")
groups %>%
  summarise(min      = min(min_n_compared),
            median   = median(min_n_compared),
            max      = max(min_n_compared),
            below_2M = sum(min_n_compared < 2e6),
            worst_n_diff = max(max_n_diff)) %>%
  print(width = Inf)

cat("\n=== five groups on the fewest compared sites ===\n")
groups %>%
  arrange(min_n_compared) %>%
  select(mlg, n_members, min_n_compared, max_n_diff, countries) %>%
  head(5) %>%
  print(n = Inf)



# figures ----------------------------------------------------------------------

### threshold diagnostics (shared log x scale, threshold line, plateau band)
band <- annotate("rect", xmin = band_lo, xmax = band_hi,
                 ymin = -Inf, ymax = Inf, fill = "grey88", colour = NA)

cut_line <- geom_vline(xintercept = clone_dist, colour = col_cut,
                       linetype = "dashed", linewidth = 0.4)

x_scale <- scale_x_continuous(trans = "log10", labels = lab_log,
                              limits = x_range,
                              expand = expansion(mult = c(0.01, 0.01)))

### A: distance distribution (log y)
n_zero <- sum(pairs$dist == 0)
## check if any pair is at distance 0
if (n_zero > 0) message(n_zero, " pairs at distance 0, not shown in figure 1A")

p_dist <- ggplot(filter(pairs, dist > 0), aes(x = dist)) +
  band + geom_histogram(bins = 90, fill = col_fill, colour = NA) + cut_line +
  x_scale +
  scale_y_continuous(trans = "log10", labels = scales::label_comma(),
                     expand = expansion(mult = c(0.08, 0.12))) +
  labs(x = NULL, y = "Sample pairs")

### B: genotypes against threshold (threshold 0 dropped on the log axis)
p_filter <- ggplot(filter(collapse, threshold > 0),
                   aes(x = threshold, y = n_mlg, colour = algorithm)) +
  band + geom_step(linewidth = 0.6) + cut_line +
  x_scale +
  scale_y_continuous(limits = c(0, length(ids)),
                     expand = expansion(mult = c(0, 0.06))) +
  scale_colour_manual(values = col_algo, name = "Clustering algorithm") +
  labs(x = NULL, y = "Multilocus genotypes")

### C: genotypes spanning two countries
p_cross <- ggplot(cross_country, aes(x = threshold, y = n_cross)) +
  band + geom_step(linewidth = 0.6, colour = col_cut) + cut_line +
  x_scale +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = "Pairwise distance / threshold",
       y = "Genotypes spanning\ntwo countries")

save_fig("threshold_diagnostics.pdf",
         p_dist / p_filter / p_cross +
           plot_layout(heights = c(1.3, 1.5, 1), axes = "collect") +
           plot_annotation(tag_levels = "A"), 7.1, 8.5)

### missingness check (A: distance against compared sites, B: worst pair per
### genotype)
p_comp <- ggplot(pairs, aes(x = n_compared / 1e6, y = dist)) +
  geom_hex(bins = 60) +
  geom_hline(yintercept = clone_dist, colour = col_cut,
             linetype = "dashed", linewidth = 0.4) +
  scale_y_continuous(trans = "log10", labels = lab_log) +
  scale_fill_viridis_c(trans = "log10", name = "Sample pairs",
                       labels = scales::label_comma()) +
  labs(x = "Sites compared (millions)", y = "Pairwise distance")

p_group <- ggplot(groups, aes(x = min_n_compared / 1e6, y = max_n_diff)) +
  geom_point(aes(size = n_members), shape = 21, alpha = 0.7,
             fill = col_algo[["farthest"]], colour = "black") +
  scale_size_continuous(range = c(1.5, 5), breaks = c(2, 4, 6, 8),
                        name = "Members") +
  labs(x = "Fewest sites compared within the genotype (millions)",
       y = "Most differing sites\nwithin the genotype")

save_fig("missingness_check.pdf",
         p_comp / p_group + plot_annotation(tag_levels = "A"), 7.1, 8)

### distance heatmap (ordered by the farthest neighbour tree)
ord <- trees$farthest$labels[trees$farthest$order]

heat <- as_tibble(dmat, rownames = "sample_1") %>%
  pivot_longer(-sample_1, names_to = "sample_2", values_to = "dist") %>%
  mutate(sample_1 = factor(sample_1, levels = ord),
         sample_2 = factor(sample_2, levels = ord))

p_heat <- ggplot(heat, aes(x = sample_1, y = sample_2, fill = dist)) +
  geom_raster() +
  scale_fill_viridis_c(name = "Pairwise distance", option = "magma",
                       direction = -1) +
  coord_fixed() +
  labs(x = paste0("Samples, ordered by the clustering (n = ", length(ids), ")"),
       y = NULL) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.key.width = unit(1.2, "cm"))

save_fig("distance_heatmap.pdf", p_heat, 7.1, 7)

### per-sample missingness
p_miss <- mlg %>%
  arrange(f_missing) %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(x = rank, y = f_missing, fill = status)) +
  geom_point(shape = 21, colour = "black", size = 1.2, alpha = 0.9) +
  scale_fill_manual(values = col_status, name = NULL) +
  scale_y_continuous(trans = "log10", labels = scales::label_percent()) +
  labs(x = "Samples, ordered by missingness", y = "Missing genotypes")

save_fig("sample_qc.pdf", p_miss, 7.1, 4.5)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  predicted_thresholds.tsv\n")
cat("  filter_stats.tsv\n")
cat("  cross_country_bound.tsv\n")
cat("  clone_groups.tsv\n")
### copy per_sample_mlg.tsv to ${ROOT}/01_variant_calling/supplementary/091_supplementary_decloned
### for hpc/01_variant_calling/091_decolned.sh and 100_metadata.sh
cat("  per_sample_mlg.tsv\n")
cat("  samples_excluded_clones.txt\n")
cat("  samples_retained_clones.txt\n")
cat("  threshold_diagnostics.pdf\n")
cat("  missingness_check.pdf\n")
cat("  distance_heatmap.pdf\n")
cat("  sample_qc.pdf\n")
