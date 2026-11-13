## ============================================================================
## 020_genotypic_diversity.R
## Genotypic diversity per unit at K = 4, 5 and 6, per unit within each sampling
## project and per sampling site within units: isolates, multilocus genotypes
## (MLG), clonal fraction, rarefied number of MLGs (eMLG), Shannon H, Stoddart &
## Taylor G, corrected Simpson index and evenness E.5; all ingroup isolates
## (clonemates included, since clonality is what is measured); sites shared by
## pairs of units
##
## Files needed:
##   lineage_labels.tsv             unit at K = 4, 5 and 6 per ingroup sample,
##                                  from local/03_delimitation/001_lineage_labels.R
##   sample_metadata_ALL_FIXED.tsv  per-sample table (mlg, site, project, mat_genome, coordinates),
##                                  from local/01_variant_calling/103_mating_type.R
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)
library(vegan)



# files ------------------------------------------------------------------------

labels_file <- "03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv"
meta_file   <- "01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv"
out_root    <- "04_diversity_niche_selection/results/020_results_genotypic_diversity"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### levels
k_levels <- c("k4", "k5", "k6")

### mating type calls of possible mixed cultures (sample set without them)
mat_exclude <- c("MAT12", "ambiguous")

### projects (sampling designs differ, so clonality is compared within a project)
lab_project <- c(magic = "MAGIC", phytoakmeter = "PhytOakmeter")

### smallest group rarefied to within a project, and smallest site with a Simpson index
n_rarefy_min <- 10
site_min     <- 5



# style ------------------------------------------------------------------------

### unit colours
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

### project fills
pal_project <- c(magic = "white", phytoakmeter = "grey25")

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

lab  <- read_tsv(labels_file, col_types = cols(.default = "c"))
meta <- read_tsv(meta_file, col_types = cols(.default = "c"))

### labelled samples without sampling coordinates (kept in the phylogeny and delimitation, not here)
excluded <- lab %>%
  select(sample, clustered) %>%
  inner_join(select(meta, sample, project, site, latitude, longitude), by = "sample") %>%
  filter(is.na(as.numeric(latitude)) | is.na(as.numeric(longitude))) %>%
  transmute(sample, project, site, clustered, reason = "no sampling coordinates")

### isolates with unit, MLG, site, project and mixed-MAT flag
iso <- lab %>%
  filter(!sample %in% excluded$sample) %>%
  select(sample, all_of(k_levels)) %>%
  inner_join(select(meta, sample, mlg, site, project, mat_genome), by = "sample") %>%
  mutate(mixed = mat_genome %in% mat_exclude)

### isolates per sample set
sets <- list(all = iso, no_mixed = filter(iso, !mixed))



# checks -----------------------------------------------------------------------

## check if every labelled sample has metadata
if (nrow(iso) != nrow(lab) - nrow(excluded)) {
  warning(nrow(lab) - nrow(excluded) - nrow(iso), " labelled samples without metadata")
}

## check if every isolate has an MLG
if (anyNA(iso$mlg) || any(iso$mlg == "")) {
  stop("isolates without an MLG")
}

## check if every site belongs to one project
n_proj <- iso %>% distinct(site, project) %>% count(site) %>% filter(n > 1)
if (nrow(n_proj) > 0) {
  warning("sites in more than one project: ", paste(n_proj$site, collapse = ", "))
}

## check if every isolate has a known project
if (!all(iso$project %in% names(lab_project))) {
  warning("isolates of other projects: ", paste(unique(setdiff(iso$project, names(lab_project))), collapse = ", "))
}

## check if every MLG sits in one unit at every level
n_split <- map_int(k_levels, function(kk) iso %>% distinct(mlg, .data[[kk]]) %>% count(mlg) %>% filter(n > 1) %>% nrow())
if (any(n_split > 0)) {
  warning("MLGs split between units: ", paste(sprintf("%s %d", k_levels, n_split), collapse = ", "))
}



# diversity --------------------------------------------------------------------

### indices from MLG counts of one group (Simpson corrected by N / (N - 1))
diversity_row <- function(counts) {
  N <- sum(counts)
  p <- counts / N
  H <- -sum(p * log(p))
  G <- 1 / sum(p^2)
  tibble(N = N, MLG = length(counts), clonal_fraction = 1 - length(counts) / N,
         H = H, G = G,
         lambda = if (N > 1) N / (N - 1) * (1 - sum(p^2)) else NA_real_,
         E5 = if (length(counts) > 1) (G - 1) / (exp(H) - 1) else NA_real_)
}

### per unit of a level, with eMLG rarefied to the smallest unit holding at least n_min isolates
unit_diversity <- function(d, kk, n_min = 1) {
  tab  <- table(d[[kk]], d$mlg)
  n_g  <- rowSums(tab)
  n_r  <- min(n_g[n_g >= n_min])
  rare <- rarefy(tab[n_g >= n_min, , drop = FALSE], sample = n_r, se = TRUE)
  map_dfr(rownames(tab), function(u) diversity_row(tab[u, tab[u, ] > 0]) %>% mutate(unit = u, .before = 1)) %>%
    mutate(eMLG    = rare[1, match(unit, colnames(rare))],
           eMLG_se = rare[2, match(unit, colnames(rare))],
           eMLG_n  = n_r)
}

### per unit, all isolates and without mixed-MAT isolates
per_unit <- expand_grid(k = k_levels, set = names(sets)) %>%
  mutate(data = map2(k, set, function(kk, ss) unit_diversity(sets[[ss]], kk))) %>%
  unnest(data)

### per unit within each project (eMLG for units with at least n_rarefy_min isolates in the project)
per_project <- expand_grid(k = k_levels, project = names(lab_project)) %>%
  mutate(data = map2(k, project, function(kk, pp) unit_diversity(filter(iso, project == pp), kk, n_rarefy_min))) %>%
  unnest(data)

### per site within unit (all isolates)
per_site <- map_dfr(k_levels, function(kk) {
  iso %>%
    group_by(unit = .data[[kk]], project, site) %>%
    group_modify(~ diversity_row(as.numeric(table(.x$mlg)))) %>%
    ungroup() %>%
    mutate(k = kk, .before = 1)
})



# co-occurrence ----------------------------------------------------------------

### sampling sites holding both units of a pair (all isolates), per level
cooc <- map_dfr(k_levels, function(kk) {
  s <- distinct(iso, unit = .data[[kk]], site)
  combn(sort(unique(s$unit)), 2, simplify = FALSE) %>%
    map_dfr(function(p) {
      a <- s$site[s$unit == p[1]]
      b <- s$site[s$unit == p[2]]
      tibble(k = kk, unit_a = p[1], unit_b = p[2], sites_a = length(a), sites_b = length(b),
             shared = length(intersect(a, b)), shared_sites = paste(sort(intersect(a, b)), collapse = ","))
    })
})



# figures ----------------------------------------------------------------------

### rarefaction curves: expected MLGs against isolates sampled, per unit (all isolates)
curves <- map_dfr(k_levels, function(kk) {
  tab <- table(iso[[kk]], iso$mlg)
  map_dfr(rownames(tab), function(u) {
    x <- tab[u, tab[u, ] > 0]
    tibble(k = kk, unit = u, n = seq_len(sum(x)), eMLG = as.numeric(rarefy(x, sample = seq_len(sum(x)))))
  })
})

plot_curves <- function(kk, y_title) {
  d   <- filter(curves, k == kk)
  end <- d %>% slice_max(n, by = unit)
  n_r <- per_unit %>% filter(k == kk, set == "all") %>% pull(eMLG_n) %>% first()

  ggplot(d, aes(x = n, y = eMLG, group = unit)) +
    geom_vline(xintercept = n_r, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
    geom_line(colour = "grey20", linewidth = 0.9) +
    geom_line(aes(colour = unit), linewidth = 0.5) +
    geom_text(data = end, aes(label = unit), hjust = -0.15, size = 7 / .pt) +
    scale_colour_manual(values = pal_unit, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = sprintf("Isolates at K = %s", sub("k", "", kk)), y = if (y_title) "Expected MLGs" else NULL)
}

### clonal fraction per unit and project (groups of at least n_rarefy_min isolates)
plot_clonal <- function(kk, y_title) {
  d <- per_project %>% filter(k == kk, N >= n_rarefy_min)

  ggplot(d, aes(x = unit, y = clonal_fraction, fill = project)) +
    geom_point(shape = 21, size = 2.2, stroke = 0.3, colour = "black",
               position = position_dodge(width = 0.5)) +
    scale_fill_manual(values = pal_project, labels = lab_project, name = NULL) +
    scale_y_continuous(limits = c(0, max(per_project$clonal_fraction[per_project$N >= n_rarefy_min])),
                       expand = expansion(mult = c(0, 0.05))) +
    labs(x = sprintf("Units at K = %s", sub("k", "", kk)), y = if (y_title) "Clonal fraction" else NULL)
}

first_col <- seq_along(k_levels) == 1
widths    <- table(distinct(per_unit, k, unit)$k)[k_levels]

p_geno <- (wrap_plots(map2(k_levels, first_col, plot_curves), nrow = 1) /
             wrap_plots(map2(k_levels, first_col, plot_clonal), nrow = 1, widths = widths)) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

save_fig("genotypic_diversity.pdf", p_geno, 7.1, 5)



# tables -----------------------------------------------------------------------

### fixed decimals (two by default), "-" for NA
fmt <- function(x, d = 2) if_else(is.na(x), "-", sprintf(paste0("%.", d, "f"), x))

### per level and unit (all isolates and without mixed-MAT isolates)
pub_unit <- per_unit %>%
  filter(set == "all") %>%
  left_join(per_unit %>%
              filter(set == "no_mixed") %>%
              select(k, unit, N_nm = N, MLG_nm = MLG, lambda_nm = lambda),
            by = c("k", "unit")) %>%
  arrange(k, unit) %>%
  transmute(Level = sub("k", "K = ", k), Unit = unit, N = N, MLG = MLG,
            `Clonal fraction` = fmt(clonal_fraction),
            eMLG = sprintf("%.1f (%.1f)", eMLG, eMLG_se),
            H = fmt(H), G = fmt(G, 1), `lambda (corrected)` = fmt(lambda, 3), E.5 = fmt(E5),
            `Without mixed MAT: N, MLG, lambda` = sprintf("%d, %d, %s", N_nm, MLG_nm, fmt(lambda_nm, 3)))

### per level, unit and project
pub_project <- per_project %>%
  arrange(k, unit, project) %>%
  transmute(Level = sub("k", "K = ", k), Unit = unit, Project = lab_project[project], N = N, MLG = MLG,
            `Clonal fraction` = fmt(clonal_fraction),
            eMLG = if_else(is.na(eMLG), "-", sprintf("%.1f (%.1f)", eMLG, eMLG_se)),
            `lambda (corrected)` = fmt(lambda, 3), E.5 = fmt(E5))

### sites shared by pairs of units
pub_cooc <- cooc %>%
  arrange(k, unit_a, unit_b) %>%
  transmute(Level = sub("k", "K = ", k), Pair = paste(unit_a, unit_b, sep = " - "),
            `Sites (first)` = sites_a, `Sites (second)` = sites_b, `Shared sites` = shared,
            `Shared, % of the smaller` = sprintf("%.0f", 100 * shared / pmin(sites_a, sites_b)))

### per site within the K = 6 units
pub_site <- per_site %>%
  filter(k == "k6") %>%
  arrange(unit, project, site) %>%
  transmute(Unit = unit, Project = lab_project[project], Site = site, N = N, MLG = MLG,
            `Clonal fraction` = fmt(clonal_fraction),
            `lambda (corrected)` = fmt(if_else(N >= site_min, lambda, NA_real_), 3))

write_tsv(excluded, file.path(out_root, "samples_excluded_no_phylogeny.txt"))
write_tsv(per_unit, file.path(out_root, "genotypic_diversity_units.tsv"))
write_tsv(per_project, file.path(out_root, "genotypic_diversity_projects.tsv"))
write_tsv(per_site, file.path(out_root, "genotypic_diversity_sites.tsv"))
write_tsv(pub_unit, file.path(out_root, "publication_genotypic_diversity.tsv"))
write_tsv(pub_project, file.path(out_root, "publication_genotypic_diversity_projects.tsv"))
write_tsv(cooc, file.path(out_root, "site_cooccurrence.tsv"))
write_tsv(pub_cooc, file.path(out_root, "publication_site_cooccurrence.tsv"))
write_tsv(pub_site, file.path(out_root, "publication_genotypic_diversity_sites.tsv"))



# summary ----------------------------------------------------------------------

cat("\n=== Genotypic diversity per level and unit ===\n")
print(pub_unit, n = Inf, width = Inf)

cat("\n=== Genotypic diversity per level, unit and project ===\n")
print(pub_project, n = Inf, width = Inf)

cat("\n=== Sites shared by pairs of units ===\n")
print(pub_cooc, n = Inf, width = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  genotypic_diversity.pdf\n")
cat("  samples_excluded_no_phylogeny.txt\n")
cat("  genotypic_diversity_units.tsv\n")
cat("  genotypic_diversity_projects.tsv\n")
cat("  genotypic_diversity_sites.tsv\n")
cat("  publication_genotypic_diversity.tsv\n")
cat("  publication_genotypic_diversity_projects.tsv\n")
cat("  site_cooccurrence.tsv\n")
cat("  publication_site_cooccurrence.tsv\n")
cat("  publication_genotypic_diversity_sites.tsv\n")
