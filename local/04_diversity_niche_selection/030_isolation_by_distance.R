## ============================================================================
## 030_isolation_by_distance.R
## Isolation by distance within every unit at K = 4, 5 and 6: Mantel test of the pairwise
## sequence distance of the single-copy orthologues against the log geographic distance
## between representatives; all representatives and without mixed-MAT samples
##
## Files needed:
##   classes.tsv                    representatives, site, units, mixed-MAT flag
##   pair_differences.tsv.gz        differences per pair of representatives over all loci
##   pair_sites.tsv.gz              compared sites per pair of representatives over all loci
##                                  (all three from hpc/04_diversity_niche_selection/010_sequence_distances.sh;
##                                   copy ${ROOT}/04_diversity_niche_selection/results/010_results_sequence_distances
##                                   to local/04_diversity_niche_selection/data/010_data_sequence_distances)
##   sample_metadata_ALL_FIXED.tsv  coordinates per sample,
##                                  from local/01_variant_calling/103_mating_type.R
## ============================================================================

library(tidyverse)
library(ggpubr)
library(vegan)



# files ------------------------------------------------------------------------

data_root <- "04_diversity_niche_selection/data/010_data_sequence_distances"
meta_file <- "01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv"
out_root  <- "04_diversity_niche_selection/results/030_results_isolation_by_distance"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### levels
k_levels <- c("k4", "k5", "k6")

### coordinate columns of the metadata (decimal degrees, WGS84)
lat_col <- "latitude"
lon_col <- "longitude"

### smallest number of sites in a unit for a Mantel test
site_min <- 3

### permutations and seed
n_perm <- 999
seed   <- 42



# style ------------------------------------------------------------------------

### unit colours (also the order of the panels)
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

### p values
fmt_p <- function(p) if_else(p < 0.001, "< 0.001", sprintf("%.3f", p))



# read -------------------------------------------------------------------------

classes <- read_tsv(file.path(data_root, "classes.tsv"), col_types = cols(.default = "c", mixed = "l"))
meta    <- read_tsv(meta_file, col_types = cols(.default = "c"))

## check if the metadata hold the coordinate columns
if (!all(c(lat_col, lon_col) %in% names(meta))) {
  stop("coordinate columns not in the metadata; columns are: ", paste(names(meta), collapse = ", "))
}

### distance matrices (representatives in rows and columns)
read_mat <- function(f) {
  m <- read.delim(gzfile(f), row.names = 1, check.names = FALSE)
  as.matrix(m)
}
diffs <- read_mat(file.path(data_root, "pair_differences.tsv.gz"))
sites <- read_mat(file.path(data_root, "pair_sites.tsv.gz"))

### p-distance per pair of representatives
p_dist <- diffs / sites
diag(p_dist) <- 0

### coordinates per representative
reps <- classes %>%
  left_join(select(meta, sample, lat = all_of(lat_col), lon = all_of(lon_col)), by = "sample") %>%
  mutate(across(c(lat, lon), as.numeric))



# checks -----------------------------------------------------------------------

## check if the matrices hold the representatives of classes.tsv in the same order
if (!identical(rownames(p_dist), classes$sample) || !identical(colnames(p_dist), classes$sample)) {
  stop("distance matrices and classes.tsv differ in samples or order")
}

## check if every representative has coordinates in range
bad <- reps %>% filter(is.na(lat) | is.na(lon) | abs(lat) > 90 | abs(lon) > 180)

if (nrow(bad) > 0) {
  warning(nrow(bad), " representatives without valid coordinates are left out: ",
          paste(bad$sample, collapse = ", "))
}

## check if the representatives of a site share coordinates (sites with spread are reported)
spread <- reps %>%
  filter(!is.na(lat), !is.na(lon)) %>%
  summarise(lat_range = diff(range(lat)), lon_range = diff(range(lon)), .by = site) %>%
  filter(lat_range > 0.01 | lon_range > 0.01)

if (nrow(spread) > 0) {
  message(nrow(spread),
          " sites with representatives more than 0.01 degrees apart (coordinates are used per representative)")
}

## check if a unit name holds the same representatives at every level
same <- map_dfr(k_levels, function(kk) tibble(k = kk, unit = classes[[kk]], sample = classes$sample)) %>%
  summarise(members = paste(sort(sample), collapse = ","), .by = c(k, unit)) %>%
  summarise(n = n_distinct(members), .by = unit)

if (any(same$n > 1)) {
  stop("unit names with different members between levels: ", paste(same$unit[same$n > 1], collapse = ", "))
}

## check if every pair has compared sites
if (any(sites[upper.tri(sites)] == 0)) {
  stop("pairs of representatives without compared sites")
}



# geographic distance ----------------------------------------------------------

### great-circle distance in km between all representatives (haversine)
haversine <- function(lat, lon) {
  phi <- lat * pi / 180
  lam <- lon * pi / 180
  h <- sin(outer(phi, phi, "-") / 2)^2 + outer(cos(phi), cos(phi)) * sin(outer(lam, lam, "-") / 2)^2
  2 * 6371.0088 * asin(pmin(1, sqrt(h)))
}

ok <- !is.na(reps$lat) & !is.na(reps$lon) & abs(reps$lat) <= 90 & abs(reps$lon) <= 180
geo <- matrix(NA_real_, nrow(reps), nrow(reps), dimnames = dimnames(p_dist))
geo[ok, ok] <- haversine(reps$lat[ok], reps$lon[ok])

### sample sets (representatives with coordinates)
sets <- list(all = reps$sample[ok], no_mixed = reps$sample[ok & !reps$mixed])



# Mantel tests -----------------------------------------------------------------

### per level, sample set and unit: p-distance against ln(1 + km), Pearson, permutations of representatives
run_ibd <- function(kk, ss, uu) {
  s  <- reps$sample[reps$sample %in% sets[[ss]] & reps[[kk]] == uu]
  g  <- geo[s, s]
  d  <- p_dist[s, s]
  ns <- n_distinct(reps$site[reps$sample %in% s])

  out <- tibble(k = kk, set = ss, unit = uu, n = length(s), sites = ns, pairs = choose(length(s), 2),
                max_km = max(g), r = NA_real_, p = NA_real_, slope = NA_real_, intercept = NA_real_)
  if (ns < site_min) return(out)

  x  <- log1p(g[upper.tri(g)])
  y  <- d[upper.tri(d)]
  mt <- mantel(as.dist(d), as.dist(log1p(g)), method = "pearson", permutations = n_perm)
  fit <- coef(lm(y ~ x))

  mutate(out, r = mt$statistic, p = mt$signif, slope = fit[2], intercept = fit[1])
}

### units per level; a unit present at several levels (same members, checked above) is tested once
units    <- map_dfr(k_levels, ~ tibble(k = .x, unit = sort(unique(reps[[.x]]))))
level_of <- distinct(units, unit, .keep_all = TRUE)
units_fig <- names(pal_unit)[names(pal_unit) %in% level_of$unit]

set.seed(seed)
per_unit <- expand_grid(set = names(sets), level_of) %>%
  pmap_dfr(function(set, k, unit) run_ibd(k, set, unit)) %>%
  select(-k)

res <- units %>%
  left_join(per_unit, by = "unit", relationship = "many-to-many") %>%
  relocate(set, .after = k)

### pairs of representatives within every unit (all representatives), for the figure
pairs <- pmap_dfr(level_of, function(k, unit) {
  s <- reps$sample[reps$sample %in% sets$all & reps[[k]] == unit]
  i <- which(upper.tri(diag(length(s))), arr.ind = TRUE)
  tibble(unit = unit, sample_a = s[i[, 1]], sample_b = s[i[, 2]],
         km = geo[s, s][i], p_distance = p_dist[s, s][i])
})



# figures ----------------------------------------------------------------------

### p-distance against geographic distance per unit (all representatives), least-squares fit on ln(1 + km),
### Mantel r and p
lab_fig <- per_unit %>%
  filter(set == "all") %>%
  mutate(label = if_else(is.na(r), sprintf("%d sites", sites),
                         sprintf("r = %.2f, p %s", r, if_else(p < 0.001, "< 0.001", sprintf("= %.3f", p)))))

fit_units <- lab_fig$unit[!is.na(lab_fig$r)]

p_ibd <- pairs %>%
  mutate(unit = factor(unit, levels = units_fig)) %>%
  ggplot(aes(x = km, y = p_distance)) +
  geom_point(size = 0.4, alpha = 0.3, colour = "grey55", stroke = 0) +
  geom_smooth(data = ~ filter(.x, unit %in% fit_units), method = "lm", formula = y ~ x, se = FALSE,
              colour = "grey15", linewidth = 1.1) +
  geom_smooth(data = ~ filter(.x, unit %in% fit_units), aes(colour = unit), method = "lm", formula = y ~ x,
              se = FALSE, linewidth = 0.6) +
  geom_text(data = mutate(lab_fig, unit = factor(unit, levels = units_fig)), aes(label = label),
            x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, size = 6 / .pt, inherit.aes = FALSE) +
  facet_wrap(~ unit, nrow = 2, scales = "free") +
  scale_x_continuous(trans = "log1p", breaks = c(0, 10, 100, 1000), limits = c(0, max(pairs$km))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_colour_manual(values = pal_unit, guide = "none") +
  labs(x = "Geographic distance (km)", y = "Sequence distance (p-distance)") +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))

save_fig("isolation_by_distance.pdf", p_ibd, 7.1, 3.6)



# tables -----------------------------------------------------------------------

### one row per level and unit (all representatives), the test without mixed-MAT samples alongside
pub <- res %>%
  filter(set == "all") %>%
  left_join(res %>% filter(set == "no_mixed") %>% select(k, unit, n_nm = n, r_nm = r, p_nm = p),
            by = c("k", "unit")) %>%
  arrange(k, unit) %>%
  transmute(Level = sub("k", "K = ", k), Unit = unit, n = n, Sites = sites,
            `Maximum distance (km)`       = sprintf("%.0f", max_km),
            `Mantel r`                    = if_else(is.na(r), "–", sprintf("%.3f", r)),
            p                             = if_else(is.na(p), "–", fmt_p(p)),
            `Slope (per ln km)`           = if_else(is.na(slope), "–", sprintf("%.2e", slope)),
            `Without mixed MAT: n, r (p)` = if_else(is.na(r_nm), sprintf("%d, –", n_nm),
                                                    sprintf("%d, %.3f (%s)", n_nm, r_nm, fmt_p(p_nm))))

write_tsv(res, file.path(out_root, "ibd.tsv"))
write_tsv(pairs, file.path(out_root, "ibd_pairs.tsv.gz"))
write_tsv(pub, file.path(out_root, "publication_ibd.tsv"))



# summary ----------------------------------------------------------------------

cat("\n=== Isolation by distance per level and unit ===\n")
print(pub, n = Inf, width = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  isolation_by_distance.pdf\n")
cat("  ibd.tsv\n")
cat("  ibd_pairs.tsv.gz\n")
cat("  publication_ibd.tsv\n")
