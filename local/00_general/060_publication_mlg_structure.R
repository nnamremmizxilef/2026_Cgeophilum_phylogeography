## ============================================================================
## 060_publication_mlg_structure.R
## Publication figure of clonality and fine-scale structure: molecular variance
## at K = 4, 5 and 6 (A), clonal fraction per unit and project (B), MAT1-1 share
## per unit (C), PhytOakmeter sampling design (D), clonal sample pairs by pit and
## depth separation (E), genets by soil pits and depth intervals occupied (F),
## distance between distinct MLGs of a unit within a site against depth and
## horizontal separation (G, H)
##
## Files needed:
##   amova.tsv                               variance components per level,
##                                           from local/04_diversity_niche_selection/012_amova.R
##   genotypic_diversity_projects.tsv        clonal fraction per unit and project,
##                                           from local/04_diversity_niche_selection/020_genotypic_diversity.R
##   mating_type_by_lineage.tsv              MAT1-1 share and 95% CI per lineage (one isolate per MLG),
##                                           from local/01_variant_calling/103_mating_type.R
##   clonal_pairs_by_separation.tsv          clonal sample pairs per pit and depth separation
##   genet_extent.tsv                        soil pits and depth intervals per genet
##                                           (both from local/01_variant_calling/091_genet_extent.R)
##   within_site_relatedness_separation.tsv  distance between distinct MLGs per separation class with
##                                           bootstrap interval
##   within_site_relatedness.tsv             permutation p of the trends with separation
##                                           (both from local/04_diversity_niche_selection/031_within_site_realtedness.R)
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

res_04   <- "04_diversity_niche_selection/results"
res_01   <- "01_variant_calling/results"
out_root <- "00_general/results/060_results_publication_mlg_structure"

f_in <- c(amova  = file.path(res_04, "012_results_amova", "amova.tsv"),
          geno   = file.path(res_04, "020_results_genotypic_diversity", "genotypic_diversity_projects.tsv"),
          mat    = file.path(res_01, "103_results_mating_type_collect", "mating_type_by_lineage.tsv"),
          pairs  = file.path(res_01, "091_results_genet_extent", "clonal_pairs_by_separation.tsv"),
          extent = file.path(res_01, "091_results_genet_extent", "genet_extent.tsv"),
          sep    = file.path(res_04, "031_results_within_site_relatedness", "within_site_relatedness_separation.tsv"),
          rel    = file.path(res_04, "031_results_within_site_relatedness", "within_site_relatedness.tsv"))

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### level of panels B, C, G and H; levels of panel A
k_main   <- "k6"
k_levels <- c("k4", "k5", "k6")

### smallest number of isolates of a unit in a project for its clonal fraction
### (as in local/04_diversity_niche_selection/020_genotypic_diversity.R)
n_min_project <- 10

### sample set of the mating type and relatedness panels
set_main <- "all"

### sampling design: pits 5 m from the tree in the four directions, depth intervals 5-10, 25-30 and 45-50 cm
### at their midpoints (m; as in local/01_variant_calling/091_genet_extent.R)
radius     <- 5
directions <- c("N", "E", "S", "W")
depth_dz   <- c(`1` = 0.075, `2` = 0.275, `3` = 0.475)

### projection of the design (x against y, flatness, depth and pit exaggeration)
stretch    <- sqrt(3)
squash     <- 0.45
z_exag     <- 28
blow_up    <- 3

### separation classes of the relatedness panels
### (as in local/04_diversity_niche_selection/031_within_site_realtedness.R)
lab_depth <- c(`1` = "5-10 cm", `2` = "25-30 cm", `3` = "45-50 cm")
lab_vsep  <- c(`0` = "Same depth", `1` = "20 cm apart", `2` = "40 cm apart")
lab_hsep  <- c(`2` = "Same pit", `1` = "Adjacent pits (7 m)", `0` = "Opposite pits (10 m)")

### variance components (as in local/04_diversity_niche_selection/012_amova.R)
lab_source <- c(`among units` = "Among lineages", `among sites within units` = "Among sites",
                `within sites` = "Within sites")
pal_source <- c(`among units` = "grey25", `among sites within units` = "grey60", `within sites` = "grey90")



# style ------------------------------------------------------------------------

### unit colours
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

### project symbols and their offset on the unit axis; pooled units; empty and full cells of the clonality panels
lab_project    <- c(magic = "MAGIC", phytoakmeter = "PhytOakmeter")
shape_project  <- c(magic = 21, phytoakmeter = 24)
offset_project <- c(magic = -0.15, phytoakmeter = 0.15)
col_all        <- "grey60"
col_none       <- "grey92"
col_full       <- "grey20"

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}

lab_k <- function(k) sub("k", "K = ", k)
fmt_p <- function(p) if_else(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))



# read -------------------------------------------------------------------------

## check if every input exists
miss <- f_in[!file.exists(f_in)]
if (length(miss) > 0) stop("missing inputs: ", paste(miss, collapse = ", "))

inp <- map(f_in, ~ read_tsv(.x, show_col_types = FALSE, guess_max = 10000))

### levels written as k4, k5, k6 everywhere
inp$mat <- mutate(inp$mat, k = tolower(k))



# checks -----------------------------------------------------------------------

## check if the per-level tables hold the main level
for (x in c("amova", "geno", "mat", "rel")) {
  if (!k_main %in% inp[[x]]$k) stop(x, " has no results at ", k_main)
}

## check if panel A has every level
if (!all(k_levels %in% inp$amova$k)) warning("AMOVA without ", paste(setdiff(k_levels, inp$amova$k), collapse = ", "))



# AMOVA (A) --------------------------------------------------------------------

p_amova <- inp$amova %>%
  filter(set == set_main) %>%
  mutate(level  = factor(lab_k(k), levels = rev(lab_k(k_levels))),
         source = factor(source, levels = rev(names(lab_source)))) %>%
  ggplot(aes(pct, level, fill = source)) +
  geom_col(width = 0.65, colour = "grey20", linewidth = 0.2) +
  geom_text(aes(label = if_else(pct >= 5, sprintf("%.0f", pct), ""), colour = source == "among units"),
            position = position_stack(vjust = 0.5), size = 6.5 / .pt, show.legend = FALSE) +
  scale_fill_manual(values = pal_source, labels = lab_source, breaks = names(lab_source), name = NULL) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_x_continuous(breaks = seq(0, 100, 25), expand = expansion(mult = c(0, 0.05))) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(x = "Molecular variance (%)", y = NULL) +
  theme(legend.key.size = unit(3, "mm"))



# clonal fraction (B) ----------------------------------------------------------

### units with enough isolates in a project, and all units of the project pooled (an MLG belongs to one unit)
geno_main <- filter(inp$geno, k == k_main)
units_cf  <- intersect(names(pal_unit), geno_main$unit)

clonal <- bind_rows(
  filter(geno_main, N >= n_min_project) %>% select(project, unit, N, MLG, clonal_fraction),
  summarise(geno_main, N = sum(N), MLG = sum(MLG), .by = project) %>%
    mutate(unit = "all", clonal_fraction = 1 - MLG / N)) %>%
  mutate(x = match(unit, c(units_cf, "all")) + offset_project[project])

p_clonal <- ggplot(clonal, aes(x, clonal_fraction, fill = unit, shape = project)) +
  geom_vline(xintercept = length(units_cf) + 0.5, colour = "grey80", linewidth = 0.3) +
  geom_point(size = 2, stroke = 0.3, colour = "black") +
  scale_fill_manual(values = c(pal_unit, all = col_all), guide = "none") +
  scale_shape_manual(values = shape_project, labels = lab_project, name = NULL) +
  scale_x_continuous(breaks = seq_len(length(units_cf) + 1), labels = c(units_cf, "All"),
                     expand = expansion(add = 0.5)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  guides(shape = guide_legend(override.aes = list(fill = "grey70"))) +
  labs(x = NULL, y = "Clonal fraction") +
  theme(legend.key.size = unit(3, "mm"))



# mating types (C) -------------------------------------------------------------

mat <- inp$mat %>%
  filter(set == set_main, k == k_main) %>%
  mutate(lineage = factor(lineage, levels = intersect(names(pal_unit), lineage)),
         label   = paste0(n, if_else(skewed, "*", "")))

p_mat <- ggplot(mat, aes(lineage, f_mat1)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_linerange(aes(ymin = lo, ymax = hi), colour = "grey30", linewidth = 0.4) +
  geom_point(aes(fill = lineage), shape = 21, size = 2.2, stroke = 0.3, colour = "black") +
  geom_text(aes(y = 1.07, label = label), size = 6.5 / .pt) +
  scale_fill_manual(values = pal_unit, guide = "none") +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "MAT1-1 share (95% CI)")



# sampling design (D) ----------------------------------------------------------

### pits at the corners of a square seen at an angle, depth planes stacked below; depth and pits drawn larger
### than they are (z_exag, blow_up)
design <- tibble(direction = directions, px = c(0, radius, 0, -radius), py = c(radius, 0, -radius, 0))
depths <- tibble(depth = names(depth_dz), dz = unname(depth_dz))
pit_r  <- 0.5 * blow_up / sqrt(2)

sx_of <- function(px, py) px * stretch
sy_of <- function(px, py, dz) py * squash - dz * z_exag

### horizontal square per row of d with half diagonal r
square <- function(d, r) {
  d %>%
    mutate(id = row_number()) %>%
    expand_grid(tibble(ox = c(0, 1, 0, -1), oy = c(1, 0, -1, 0))) %>%
    mutate(sx = sx_of(px + ox * r, py + oy * r), sy = sy_of(px + ox * r, py + oy * r, dz))
}

dz0     <- min(depths$dz)
planes  <- square(tibble(px = 0, py = 0, dz = rev(depths$dz)), radius)
blocks  <- square(arrange(expand_grid(design, depths), sy_of(px, py, dz)), pit_r)
columns <- mutate(design, sx = sx_of(px, py), sy0 = sy_of(px, py, dz0), sy1 = sy_of(px, py, max(depths$dz)))
east    <- design %>% filter(direction == "E") %>% mutate(sx = sx_of(px, py), sy = sy_of(px, py, dz0))
dir_lab <- mutate(design, sx = sx_of(px, py), sy = sy_of(px, py, max(depths$dz)))

p_design <- ggplot(mapping = aes(x = sx, y = sy)) +
  geom_segment(data = columns, aes(xend = sx, y = sy0, yend = sy1), colour = "grey70", linewidth = 0.3) +
  geom_polygon(data = planes, aes(group = id), alpha = 0.55, fill = "grey93", colour = "grey75", linewidth = 0.3) +
  geom_polygon(data = blocks, aes(group = id), fill = col_none, colour = "grey45", linewidth = 0.3) +
  geom_segment(data = east, aes(x = 0, xend = sx, yend = sy), colour = "grey45", linewidth = 0.35) +
  geom_point(data = bind_rows(east, mutate(east, sx = 0)), shape = 3, size = 1.6, stroke = 0.5, colour = "grey30") +
  geom_text(data = east, aes(x = sx / 2), label = "5 m", vjust = -0.7, size = 6.5 / .pt) +
  geom_text(data = filter(blocks, direction == "E", ox == 1), aes(label = lab_depth[depth]),
            hjust = -0.25, size = 6.5 / .pt) +
  geom_text(data = dir_lab, aes(label = direction), size = 6.5 / .pt) +
  coord_fixed(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(),
        plot.margin = margin(4, 30, 4, 4))



# clonal pairs by separation (E) -----------------------------------------------

pairs <- inp$pairs %>%
  mutate(d_h = factor(sprintf("%.0f", d_h), levels = sprintf("%.0f", sort(unique(d_h)))),
         d_v = factor(sprintf("%.0f", 100 * d_v), levels = sprintf("%.0f", 100 * sort(unique(d_v)))))
p_mid <- max(pairs$p_clonal) / 2

p_pairs <- ggplot(pairs, aes(d_h, d_v, fill = p_clonal)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%d/%d", n_clonal, n_pairs), colour = p_clonal > p_mid), size = 6 / .pt,
            show.legend = FALSE) +
  scale_fill_gradient(low = col_none, high = col_full, limits = c(0, NA), name = "Clonal share") +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white")) +
  guides(fill = guide_colourbar(barwidth = unit(18, "mm"), barheight = unit(1.5, "mm"), title.vjust = 1)) +
  labs(x = "Distance between pits (m)", y = "Distance between depths (cm)") +
  theme(axis.line = element_blank(), axis.ticks = element_blank())



# genet extent (F) -------------------------------------------------------------

grid <- inp$extent %>%
  count(n_pits, n_depths) %>%
  complete(n_pits = seq_along(directions), n_depths = seq_along(depth_dz), fill = list(n = 0))
g_mid <- max(grid$n) / 2

p_grid <- ggplot(grid, aes(factor(n_pits), factor(n_depths), fill = n)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = n, colour = n > g_mid), size = 7 / .pt, show.legend = FALSE) +
  scale_fill_gradient(low = col_none, high = col_full, name = "Genets") +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white")) +
  guides(fill = guide_colourbar(barwidth = unit(18, "mm"), barheight = unit(1.5, "mm"), title.vjust = 1)) +
  labs(x = "Soil pits occupied", y = "Depth intervals occupied") +
  theme(axis.line = element_blank(), axis.ticks = element_blank())



# relatedness of distinct MLGs (G, H) ------------------------------------------

panel_lv <- intersect(c(names(pal_unit), "all"), unique(inp$sep$panel))

sep <- inp$sep %>%
  mutate(panel = factor(panel, levels = panel_lv))

### one distance axis for both separations
sep_lim <- range(c(sep$lo, sep$hi, sep$mean), na.rm = TRUE)

sep_p <- inp$rel %>%
  filter(set == set_main, k == k_main, term %in% c("v_trend", "h_trend"), unit %in% panel_lv) %>%
  transmute(panel = factor(unit, levels = panel_lv), axis = if_else(term == "v_trend", "depth", "direction"),
            label = fmt_p(p))

## check if every relatedness panel has its p value
if (nrow(sep_p) != 2 * length(panel_lv)) warning("relatedness panels without p value at ", k_main)

plot_sep <- function(ax, lv, x_lab) {
  d <- filter(sep, axis == ax) %>% mutate(level = factor(level, levels = lv))
  ggplot(d, aes(level, mean, group = 1)) +
    geom_hline(yintercept = 100, linetype = "dashed", linewidth = 0.3, colour = "grey60") +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.4, colour = "grey40") +
    geom_line(linewidth = 0.4, colour = "grey40") +
    geom_point(aes(fill = panel), shape = 21, size = 2, stroke = 0.3, colour = "black") +
    geom_text(data = filter(sep_p, axis == ax), aes(x = -Inf, y = Inf, label = label), hjust = -0.08, vjust = 1.4,
              size = 6 / .pt, inherit.aes = FALSE) +
    facet_wrap(~ panel, nrow = 1, labeller = as_labeller(function(x) if_else(x == "all", "All", x))) +
    scale_fill_manual(values = c(pal_unit, all = col_all), guide = "none") +
    scale_y_continuous(limits = sep_lim, expand = expansion(mult = c(0.05, 0.2))) +
    labs(x = x_lab, y = "Distance (% of site mean)") +
    theme(strip.background = element_blank(), axis.text.x = element_text(angle = 35, hjust = 1))
}

p_vsep <- plot_sep("depth", lab_vsep, NULL)
p_hsep <- plot_sep("direction", lab_hsep, NULL)



# figures ----------------------------------------------------------------------

row_1 <- (p_amova | p_clonal | p_mat) + plot_layout(widths = c(1.2, 1, 1))
row_2 <- (p_design | p_pairs | p_grid) + plot_layout(widths = c(1.15, 1, 0.9))

p_fig <- row_1 / row_2 / p_vsep / p_hsep +
  plot_layout(heights = c(1, 1.15, 1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold"))

save_fig(sprintf("publication_mlg_structure_%s.pdf", k_main), p_fig, 7.1, 9.4)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat(sprintf("  publication_mlg_structure_%s.pdf\n", k_main))
