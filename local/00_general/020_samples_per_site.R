## ============================================================================
## 020_samples_per_site.R
## Isolates and multilocus genotypes per sampling site within each unit at
## K = 4, 5 and 6, and the number of sites per unit that reach a minimum sample
## size for unit-specific population genetics
##
## Files needed:
##   lineage_labels.tsv  unit at K = 4, 5 and 6, project, site and MLG per ingroup sample,
##                       from local/03_delimitation/001_lineage_labels.R
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

labels_file <- "03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv"
out_root    <- "00_general/results/020_results_samples_per_site"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### levels
k_levels <- c("k4", "k5", "k6")

### smallest minimum sample size per site (first grid column, dashed line in the bars)
n_min_from <- 6

### minimum sample sizes in the publication table
n_show <- c(6, 8, 10, 15, 20)

### expected ingroup samples
n_ingroup_expected <- 585



# style ------------------------------------------------------------------------

### unit colours
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

### count and level labels
lab_count <- c(isolates = "Isolates", mlg = "Multilocus genotypes")
lab_axis  <- c(isolates = "isolates", mlg = "multilocus genotypes")
lab_short <- c(isolates = "Isolates", mlg = "MLGs")
lab_level <- c(k4 = "K = 4", k5 = "K = 5", k6 = "K = 6")

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

labels <- read_tsv(labels_file, col_types = cols(.default = "c"))



# checks -----------------------------------------------------------------------

## check if the labels hold the columns needed
need <- c("sample", "project", "site", "mlg", k_levels)
if (!all(need %in% names(labels))) {
  stop(basename(labels_file), " lacks: ", paste(setdiff(need, names(labels)), collapse = ", "))
}

## check if the labels hold the full ingroup
if (nrow(labels) != n_ingroup_expected) {
  warning(nrow(labels), " samples in the labels, expected ", n_ingroup_expected)
}

## check if every sample has a site
if (anyNA(labels$site)) {
  warning(sum(is.na(labels$site)), " samples without a site (left out)")
}

## check if every sample has an MLG
if (anyNA(labels$mlg)) {
  warning(sum(is.na(labels$mlg)), " samples without an MLG")
}

## check if every site belongs to one project
site_proj <- labels %>% distinct(site, project) %>% count(site) %>% filter(n > 1)
if (nrow(site_proj) > 0) {
  warning("sites in more than one project: ", paste(site_proj$site, collapse = ", "))
}



# samples per site -------------------------------------------------------------

### isolates and MLGs per unit and site (an MLG counts once per site)
per_site <- map_dfr(k_levels, function(kk) {
  labels %>%
    filter(!is.na(site)) %>%
    rename(unit = all_of(kk)) %>%
    summarise(isolates = n(), mlg = n_distinct(mlg), .by = c(unit, project, site)) %>%
    mutate(k = kk, .before = 1)
}) %>%
  arrange(k, unit, desc(isolates), site)



# sites per minimum sample size ------------------------------------------------

### sites per unit with at least n_min isolates or MLGs
thresholds <- seq(n_min_from, max(n_show, per_site$isolates))

sites_by_n <- per_site %>%
  pivot_longer(c(isolates, mlg), names_to = "count", values_to = "n") %>%
  cross_join(tibble(n_min = thresholds)) %>%
  summarise(sites = sum(n >= n_min), .by = c(k, unit, count, n_min))



# figures ----------------------------------------------------------------------

### sites per unit at each minimum sample size (one row per unit, one block per level)
plot_grid <- function(cnt) {
  d <- sites_by_n %>%
    filter(count == cnt) %>%
    mutate(k     = factor(k, names(lab_level), lab_level),
           unit  = fct_rev(factor(unit, names(pal_unit))),
           label = if_else(sites > 0, as.character(sites), ""),
           dark  = sites > max(sites) / 2)
  marks <- distinct(d, k, unit)
  ggplot(d, aes(x = n_min, y = unit)) +
    geom_tile(aes(fill = sites), colour = "grey70", linewidth = 0.2) +
    geom_text(aes(label = label, colour = dark), size = 2.3, show.legend = FALSE) +
    geom_point(data = marks, aes(x = n_min_from - 1, y = unit), shape = 22, size = 2.8,
               colour = "black", fill = pal_unit[as.character(marks$unit)]) +
    facet_grid(k ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient(low = "white", high = "grey20", limits = c(0, NA), name = "Sites") +
    scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white")) +
    scale_x_continuous(breaks = thresholds, limits = c(n_min_from - 1.6, max(thresholds) + 0.5),
                       expand = c(0, 0)) +
    labs(x = sprintf("Minimum %s per site", lab_axis[cnt]), y = NULL) +
    theme(axis.line = element_blank(), axis.ticks.y = element_blank(),
          legend.key.height = unit(0.3, "cm"))
}

### samples per site, sites ranked by size within each unit (one row per level)
plot_bars <- function(cnt) {
  d <- per_site %>%
    mutate(n = .data[[cnt]]) %>%
    arrange(k, unit, desc(n), site) %>%
    mutate(rank = row_number(), .by = c(k, unit)) %>%
    mutate(k = factor(k, names(lab_level), lab_level), unit = factor(unit, names(pal_unit)))
  y_max <- max(d$n)
  lv    <- levels(droplevels(d$k))
  rows  <- map(seq_along(lv), function(i) {
    ggplot(filter(d, k == lv[i]), aes(x = rank, y = n, fill = unit)) +
      geom_hline(yintercept = n_min_from, colour = "grey30", linetype = "dashed", linewidth = 0.3) +
      geom_col(colour = "black", linewidth = 0.15, width = 0.85) +
      facet_grid(k ~ unit, scales = "free_x", space = "free_x") +
      scale_fill_manual(values = pal_unit, guide = "none") +
      scale_x_continuous(breaks = NULL, expand = expansion(add = 0.6)) +
      scale_y_continuous(limits = c(0, y_max), expand = expansion(mult = c(0, 0.05)),
                         breaks = function(l) seq(0, floor(l[2]), by = max(1, round(l[2] / 5)))) +
      labs(x = if (i == length(lv)) "Sites, ranked by sample size within unit" else NULL,
           y = if (i == ceiling(length(lv) / 2)) sprintf("%s per site", lab_short[cnt]) else NULL) +
      theme(panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3))
  })
  wrap_plots(rows, ncol = 1) +
    plot_annotation(tag_levels = "A")
}

for (cnt in names(lab_count)) {
  save_fig(sprintf("sites_by_min_samples_%s.pdf", cnt), plot_grid(cnt), 7.1, 4.6)
  save_fig(sprintf("samples_per_site_%s.pdf", cnt), plot_bars(cnt), 7.1, 5.5)
}



# tables -----------------------------------------------------------------------

### sites per unit at the shown minimum sample sizes (isolates / MLGs)
pub_sites <- sites_by_n %>%
  filter(n_min %in% n_show) %>%
  pivot_wider(names_from = count, values_from = sites) %>%
  mutate(cell = sprintf("%d / %d", isolates, mlg)) %>%
  select(k, unit, n_min, cell) %>%
  pivot_wider(names_from = n_min, values_from = cell, names_prefix = ">= ") %>%
  left_join(per_site %>% summarise(sites = n(), .by = c(k, unit)), by = c("k", "unit")) %>%
  transmute(Level = lab_level[k], Unit = unit, `Sites (any)` = sites, across(starts_with(">= "))) %>%
  arrange(Level, Unit)

write_tsv(per_site, file.path(out_root, "samples_per_site.tsv"))
write_tsv(sites_by_n, file.path(out_root, "sites_by_min_samples.tsv"))
write_tsv(pub_sites, file.path(out_root, "publication_sites_by_min_samples.tsv"))

cat("\n=== sites per unit with at least n samples (isolates / MLGs) ===\n")
print(pub_sites, n = Inf, width = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  sites_by_min_samples_isolates.pdf\n")
cat("  samples_per_site_isolates.pdf\n")
cat("  sites_by_min_samples_mlg.pdf\n")
cat("  samples_per_site_mlg.pdf\n")
cat("  samples_per_site.tsv\n")
cat("  sites_by_min_samples.tsv\n")
cat("  publication_sites_by_min_samples.tsv\n")
