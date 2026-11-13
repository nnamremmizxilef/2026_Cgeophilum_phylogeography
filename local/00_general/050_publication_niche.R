## ============================================================================
## 050_publication_niche.R
## Publication niche figure: sampling sites with the lineage composition (MLGs)
## per site, MAGIC sites in an inset at the same pie scale (A); soil partial RDA
## within the PhytOakmeter and within the MAGIC sites (B, C); climate partial RDA
## across all sites (D); variance partitioning into soil, depth, climate and site
## context at K = 4, 5 and 6 (E)
##
## Files needed:
##   site_climate.tsv               site coordinates, project and host,
##                                  from local/04_diversity_niche_selection/041_climate_data.R
##   soil_niche_scores.tsv          soil triplot scores and model fit per project,
##                                  from local/04_diversity_niche_selection/040_soil_niche.R
##   climate_niche_scores.tsv       climate triplot scores and model fit
##   soil_climate_varpart.tsv       soil, depth, climate and site context fractions
##                                  (both from local/04_diversity_niche_selection/042_climate_niche.R)
##   lineage_labels.tsv             unit at K = 4, 5 and 6 per ingroup sample,
##                                  from local/03_delimitation/001_lineage_labels.R
##   sample_metadata_ALL_FIXED.tsv  site and MLG per sample,
##                                  from local/01_variant_calling/103_mating_type.R
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)
library(ggrepel)
library(ggforce)
library(sf)



# files ------------------------------------------------------------------------

res_root    <- "04_diversity_niche_selection/results"
labels_file <- "03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv"
meta_file   <- "01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv"
out_root    <- "00_general/results/050_results_publication_niche"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### level of panels A to D; levels of panel E
k_main   <- "k6"
k_levels <- c("k4", "k5", "k6")

### map projection
crs_map <- "EPSG:3035"

### map window: margin around all sites (m), widened to this aspect (width / height of the map row)
map_margin <- 250000
map_aspect <- 1.7

### MAGIC inset: extent (degrees), height and distance to the top and left border (shares of the map height)
ext_ch       <- c(xmin = 6.2, xmax = 10.3, ymin = 45.6, ymax = 48.2)
inset_h      <- 0.52
inset_margin <- 0.03

### pie radius per square root of the MLG count on the main map (m; the inset at the same paper scale)
r_eu <- 22000

### gap between neighbouring pies (share of the radii)
pie_gap <- 0.15

### MLG counts in the size legend
n_leg <- c(5, 10, 20)



# style ------------------------------------------------------------------------

### unit colours
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

### symbols of the PhytOakmeter depths, MAGIC hosts and site groups (project and host); level fills
lab_layer   <- c(phytoakmeter_1 = "5-10 cm", phytoakmeter_2 = "25-30 cm", phytoakmeter_3 = "45-50 cm")
shape_layer <- c(phytoakmeter_1 = 16, phytoakmeter_2 = 1, phytoakmeter_3 = 4)
lab_host    <- c(Quercus_robur = "Oak", Fagus_sylvatica = "Beech")
shape_host  <- c(Quercus_robur = 17, Fagus_sylvatica = 2)
lab_group   <- c(phytoakmeter_Quercus_robur = "PhytOakmeter, oak", magic_Quercus_robur = "MAGIC, oak",
                 magic_Fagus_sylvatica = "MAGIC, beech")
shape_group <- c(phytoakmeter_Quercus_robur = 16, magic_Quercus_robur = 17, magic_Fagus_sylvatica = 2)
pal_level   <- c(k4 = "grey80", k5 = "grey55", k6 = "grey25")

### variance fractions
lab_frac <- c(soil_only = "Soil only", depth_only = "Depth only", climate_only = "Climate only",
              context_only = "Site context only", shared = "Shared")

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}

lab_k <- function(k) sub("k", "K = ", k)



# read -------------------------------------------------------------------------

site_env <- read_tsv(file.path(res_root, "041_results_climate_data/site_climate.tsv"),
                     col_types = cols(site = "c", project = "c", host = "c", group = "c", .default = "d"))
soil_sc  <- read_tsv(file.path(res_root, "040_results_soil_niche/soil_niche_scores.tsv"),
                     col_types = cols(depth = "c", .default = col_guess()))
clim_sc  <- read_tsv(file.path(res_root, "042_results_climate_niche/climate_niche_scores.tsv"), show_col_types = FALSE)
vp4      <- read_tsv(file.path(res_root, "042_results_climate_niche/soil_climate_varpart.tsv"), show_col_types = FALSE)

lab  <- read_tsv(labels_file, col_types = cols(.default = "c"))
meta <- read_tsv(meta_file, col_types = cols(.default = "c"))

### MLGs per unit and site at the main level
comp <- lab %>%
  select(sample, unit = all_of(k_main)) %>%
  inner_join(select(meta, sample, site, mlg), by = "sample") %>%
  inner_join(select(site_env, site, project, lat, lon), by = "site") %>%
  distinct(project, site, lat, lon, mlg, unit) %>%
  count(project, site, lat, lon, unit)



# checks -----------------------------------------------------------------------

## check if every panel has data at the main level
for (x in list(soil = soil_sc, clim = clim_sc, vp = vp4)) {
  if (!k_main %in% x$k) stop("no results at ", k_main)
}



# map (A) ----------------------------------------------------------------------

### box of an extent in the map projection
box <- function(e) st_as_sfc(st_bbox(e, crs = 4326)) %>% st_segmentize(10000) %>% st_transform(crs_map)
lim <- function(e) st_bbox(box(e))

world <- st_as_sf(maps::map("world", plot = FALSE, fill = TRUE)) %>%
  st_make_valid() %>%
  st_geometry() %>%
  st_crop(xmin = -30, xmax = 50, ymin = 30, ymax = 75) %>%
  st_transform(crs_map)

### map window: all sites with a margin, widened to map_aspect
xy_sites <- st_coordinates(st_transform(st_as_sf(site_env, coords = c("lon", "lat"), crs = 4326), crs_map))
xr <- range(xy_sites[, 1]) + c(-1, 1) * map_margin
yr <- range(xy_sites[, 2]) + c(-1, 1) * map_margin
if (diff(xr) / diff(yr) < map_aspect) {
  xr <- mean(xr) + c(-1, 1) * diff(yr) * map_aspect / 2
} else {
  yr <- mean(yr) + c(-1, 1) * diff(xr) / map_aspect / 2
}

### inset window and its box on the map (shares of the panel), shaped like the inset (same distance to the top
### and left border on paper); pie radius of the inset for the paper scale of the main map
l_ch  <- lim(ext_ch)
ch_w  <- as.numeric(l_ch["xmax"] - l_ch["xmin"])
ch_h  <- as.numeric(l_ch["ymax"] - l_ch["ymin"])
asp   <- diff(xr) / diff(yr)
inset_box <- c(left = inset_margin / asp, bottom = 1 - inset_margin - inset_h,
               right = (inset_margin + inset_h * ch_w / ch_h) / asp, top = 1 - inset_margin)
r_ch  <- r_eu * ch_h / (inset_h * diff(yr))

### sites in the map projection with pie radius r_scale * sqrt(MLGs)
site_xy <- function(d, r_scale) {
  s  <- summarise(d, n = sum(n), .by = c(site, lat, lon))
  xy <- st_coordinates(st_transform(st_as_sf(s, coords = c("lon", "lat"), crs = 4326), crs_map))
  mutate(s, x0 = xy[, 1], y0 = xy[, 2], r = r_scale * sqrt(n))
}

### pie centres pushed apart until no two pies overlap, kept inside the map window (sites at the same point
### split by direction)
repel <- function(s, xlim, ylim, n_iter = 2000) {
  x <- s$x0
  y <- s$y0
  if (nrow(s) > 1) {
    for (it in seq_len(n_iter)) {
      moved <- FALSE
      for (i in seq_len(nrow(s) - 1)) for (j in (i + 1):nrow(s)) {
        dx <- x[j] - x[i]
        dy <- y[j] - y[i]
        d  <- sqrt(dx^2 + dy^2)
        need <- (s$r[i] + s$r[j]) * (1 + pie_gap)
        if (d < need) {
          if (d < 1) {
            dx <- cos(j)
            dy <- sin(j)
            d  <- 1
          }
          push <- (need - d) / 2 / d
          x[i] <- x[i] - push * dx
          y[i] <- y[i] - push * dy
          x[j] <- x[j] + push * dx
          y[j] <- y[j] + push * dy
          moved <- TRUE
        }
      }
      x <- pmin(pmax(x, xlim[1] + s$r), xlim[2] - s$r)
      y <- pmin(pmax(y, ylim[1] + s$r), ylim[2] - s$r)
      if (!moved) break
    }
  }
  mutate(s, x = x, y = y, moved = sqrt((x - x0)^2 + (y - y0)^2) > 0.3 * r)
}

### pie layers: leader line and point at the sampling site for moved pies, wedges per unit
pie_layers <- function(d, s) {
  w <- d %>%
    left_join(select(s, site, x, y, r), by = "site") %>%
    mutate(unit = factor(unit, levels = names(pal_unit))) %>%
    arrange(site, unit) %>%
    mutate(end = 2 * pi * cumsum(n) / sum(n), start = lag(end, default = 0), .by = site)
  list(
    geom_segment(data = filter(s, moved), aes(x = x0, y = y0, xend = x, yend = y), linewidth = 0.25,
                 colour = "grey30"),
    geom_point(data = filter(s, moved), aes(x = x0, y = y0), shape = 21, size = 0.5, stroke = 0.2,
               colour = "black", fill = "grey15"),
    geom_arc_bar(data = w, aes(x0 = x, y0 = y, r0 = 0, r = r, start = start, end = end, fill = unit),
                 colour = "grey20", linewidth = 0.3)
  )
}

s_eu <- repel(site_xy(filter(comp, project == "phytoakmeter"), r_eu), xr, yr)
s_ch <- repel(site_xy(filter(comp, project == "magic"), r_ch), as.numeric(l_ch[c("xmin", "xmax")]),
              as.numeric(l_ch[c("ymin", "ymax")]))

### size legend: circles side by side in the lower left corner of the main map, MLG counts below
gap_leg <- 0.01 * diff(xr)
leg <- tibble(n = n_leg, r = r_eu * sqrt(n_leg)) %>%
  mutate(x = xr[1] + 0.03 * diff(xr) + cumsum(2 * r + gap_leg) - r - gap_leg,
         y = yr[1] + 0.07 * diff(yr) + r, y_lab = yr[1] + 0.035 * diff(yr))

p_map <- ggplot() +
  geom_sf(data = world, fill = "grey93", colour = "grey70", linewidth = 0.15) +
  geom_sf(data = box(ext_ch), fill = NA, colour = "grey20", linewidth = 0.3, linetype = "dashed") +
  pie_layers(filter(comp, project == "phytoakmeter"), s_eu) +
  geom_circle(data = leg, aes(x0 = x, y0 = y, r = r), fill = "white", colour = "grey20", linewidth = 0.3) +
  geom_text(data = leg, aes(x = x, y = y_lab, label = n), size = 6.5 / .pt) +
  annotate("text", x = mean(range(leg$x)), y = yr[1] + 0.1 * diff(yr) + 2 * max(leg$r), label = "MLGs",
           size = 6.5 / .pt) +
  scale_fill_manual(values = pal_unit, name = NULL, drop = TRUE) +
  coord_sf(crs = crs_map, xlim = xr, ylim = yr, expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(legend.position = "inside", legend.position.inside = c(0.99, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", colour = "grey60", linewidth = 0.2),
        panel.grid = element_blank(), legend.key.size = unit(3, "mm"))

p_ch <- ggplot() +
  geom_sf(data = world, fill = "grey93", colour = "grey70", linewidth = 0.15) +
  pie_layers(filter(comp, project == "magic"), s_ch) +
  scale_fill_manual(values = pal_unit, guide = "none") +
  coord_sf(crs = crs_map, xlim = l_ch[c("xmin", "xmax")], ylim = l_ch[c("ymin", "ymax")], expand = FALSE) +
  theme_void() +
  theme(panel.border = element_rect(fill = NA, colour = "grey20", linewidth = 0.3),
        panel.background = element_rect(fill = "white", colour = NA))

p_a <- p_map + inset_element(p_ch, left = inset_box["left"], bottom = inset_box["bottom"],
                             right = inset_box["right"], top = inset_box["top"], ignore_tag = TRUE)



# triplots (B to D) ------------------------------------------------------------

### sites or cells as grey symbols, units as filled points, variables as arrows
triplot <- function(d, shape_col, shapes, shape_labs) {
  st <- filter(d, type %in% c("site", "cell"))
  sp <- filter(d, type == "unit")
  bp <- filter(d, type == "variable")
  ggplot() +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_point(data = st, aes(A1, A2, shape = .data[[shape_col]]), colour = "grey55", size = 1.1) +
    geom_segment(data = bp, aes(0, 0, xend = A1, yend = A2), arrow = arrow(length = unit(1.5, "mm")),
                 linewidth = 0.35, colour = "grey15") +
    geom_text_repel(data = bp, aes(A1, A2, label = label), size = 6.5 / .pt, colour = "black", seed = 1) +
    geom_point(data = sp, aes(A1, A2, fill = id), shape = 21, size = 2.6, stroke = 0.4, colour = "black") +
    geom_text_repel(data = sp, aes(A1, A2, label = id), size = 7 / .pt, fontface = "bold", seed = 2,
                    box.padding = 0.4) +
    scale_fill_manual(values = pal_unit, guide = "none") +
    scale_shape_manual(values = shapes, labels = shape_labs, breaks = names(shape_labs), name = NULL) +
    annotate("text", x = -Inf, y = Inf, label = d$fit[1], hjust = -0.05, vjust = 1.4, size = 6.5 / .pt) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    labs(x = d$axis1[1], y = d$axis2[1]) +
    theme(legend.key.width = unit(2, "mm"))
}

p_b <- triplot(filter(soil_sc, k == k_main, model == "phytoakmeter"), "layer", shape_layer, lab_layer)
p_c <- triplot(filter(soil_sc, k == k_main, model == "magic"), "host", shape_host, lab_host)
p_d <- triplot(filter(clim_sc, k == k_main), "group", shape_group, lab_group)



# variance partitioning (E) ----------------------------------------------------

p_e <- vp4 %>%
  filter(set == "all", part != "residual") %>%
  mutate(part = factor(part, levels = rev(names(lab_frac))), k = factor(k, levels = k_levels)) %>%
  ggplot(aes(x = 100 * adj_r2, y = part, fill = k)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, colour = "grey20", linewidth = 0.2) +
  scale_y_discrete(labels = lab_frac) +
  scale_fill_manual(values = pal_level, labels = lab_k(k_levels), breaks = k_levels, name = NULL) +
  labs(x = "Explained (adjusted R², %)", y = NULL) +
  theme(legend.position = "inside", legend.position.inside = c(1, 0.5), legend.justification = c(1, 0.5),
        legend.direction = "vertical", legend.key.size = unit(3, "mm"))



# figures ----------------------------------------------------------------------

p_niche <- wrap_elements(full = p_a) / (p_b | p_c) / (p_d | p_e) +
  plot_layout(heights = c(5.6, 2.4, 2.4)) +
  plot_annotation(tag_levels = "A")

save_fig(sprintf("publication_niche_%s.pdf", k_main), p_niche, 7.1, 9.4)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat(sprintf("  publication_niche_%s.pdf\n", k_main))
