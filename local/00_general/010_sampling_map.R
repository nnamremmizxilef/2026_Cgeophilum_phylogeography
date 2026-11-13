## ============================================================================
## 010_sampling_map.R
## Sampling maps: European background, EEA environmental zones, the Quercus
## robur distribution range, the sampling sites and the lineages per site
##
## Files needed:
##   env-zones*.tif                 environmental zones raster (EEA)
##   Quercus_robur_plg*.shp         Quercus robur distribution polygons
##                                  (both external: copy manually into
##                                   local/00_general/data/010_data_sampling_map)
##   sample_metadata_ALL_FIXED.tsv  project, site and coordinates per sample,
##                                  from local/01_variant_calling/103_mating_type.R
##   lineage_labels.tsv             units at K = 4, 5 and 6 and MLG per sample,
##                                  from local/03_delimitation/001_lineage_labels.R
## ============================================================================

library(tidyverse)
library(sf)
library(terra)
library(ggpattern)
library(ggspatial)
library(rnaturalearth)
library(ggrepel)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

data_root <- "00_general/data/010_data_sampling_map"
out_root  <- "00_general/results/010_results_sampling_map"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

### input files matched by pattern (date stamps in the file names)
input <- function(pattern) list.files(data_root, pattern, full.names = TRUE)[1]

f_meta    <- "01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv"
f_lineage <- "03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv"
f_zones   <- input("env-zones.*\\.tif$")



# parameters -------------------------------------------------------------------

### map projection (ETRS89 / LAEA Europe: equal area, the projection of the zone raster)
crs_map <- 3035

### map windows (degrees, converted once): Europe and the Swiss window of the MAGIC sites
view_eu <- st_bbox(c(xmin = -13, ymin = 33, xmax = 47, ymax = 71), crs = 4326) %>%
  st_as_sfc() %>% st_transform(crs_map) %>% st_bbox()

view_ch <- st_bbox(c(xmin = 6.3, ymin = 45.7, xmax = 10.4, ymax = 48.1), crs = 4326) %>%
  st_as_sfc() %>% st_transform(crs_map) %>% st_bbox()

### lineage level of the pies (column of lineage_labels.tsv)
k_map <- "k6"

### pie radius in map metres, smallest and largest site (area by genotypes)
pie_eu <- c(min = 25e3, max = 90e3)
pie_ch <- c(min = 2.5e3, max = 8e3)



# style ------------------------------------------------------------------------

### project colours of the site markers (dark enough against the zone colours)
col_project <- c(magic = "#D7191C", phytoakmeter = "#2C7BB6")

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

### PDF and PNG of a map
save_map <- function(name, plot, width = 7.1, height = 5.1) {
  save_fig(paste0(name, ".pdf"), plot, width, height)
  ggsave(file.path(out_root, paste0(name, ".png")), plot,
         width = width, height = height, dpi = 400)
}



# read -------------------------------------------------------------------------

meta <- read_tsv(f_meta, show_col_types = FALSE)
lin  <- read_tsv(f_lineage, show_col_types = FALSE)



# background -------------------------------------------------------------------

land_all <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(crs_map)

land_for <- function(view) suppressWarnings(st_crop(land_all, view))



# environmental zones ----------------------------------------------------------

### zone names and colours (raster codes 1 to 15 in this order)
zone_colours <- c(
  "Alpine North"            = "#584060",
  "Boreal"                  = "#6C609C",
  "Nemoral"                 = "#7D80C3",
  "Atlantic North"          = "#8CA5D4",
  "Alpine South"            = "#8BD0F0",
  "Continental"             = "#96D0C9",
  "Atlantic Central"        = "#A5D39E",
  "Pannonian"               = "#B5D06A",
  "Lusitanian"              = "#D2DC77",
  "Anatolian"               = "#EFDF81",
  "Mediterranean Mountains" = "#F7CA7E",
  "Mediterranean North"     = "#F9A366",
  "Mediterranean South"     = "#F88057",
  "Macaronesia"             = "#D66E4D",
  "Arctic"                  = "#B9BFC4"
)

### zone raster of a window as a data frame: 100 m cells cropped and aggregated to about 2000 cells
### across (most frequent class per block), zones as a factor in gradient order
zone_df_for <- function(view) {
  win <- rast(f_zones) %>%
    crop(ext(view[c("xmin", "xmax", "ymin", "ymax")]))

  agg <- aggregate(win, fact = max(1, floor(max(dim(win)) / 2000)),
                   fun = "modal", na.rm = TRUE)

  out <- agg %>%
    as.data.frame(xy = TRUE, na.rm = TRUE) %>%
    rename(code = 3) %>%
    mutate(zone = factor(names(zone_colours)[code],
                         levels = names(zone_colours)))

  ## check if every zone code has a label
  extra <- setdiff(unique(out$code), seq_along(zone_colours))
  if (length(extra) > 0) {
    warning("zone codes without a label: ", paste(extra, collapse = ", "))
  }
  out
}

zone_eu <- zone_df_for(view_eu)
zone_ch <- zone_df_for(view_ch)



# distribution range -----------------------------------------------------------

### missing .shx index rebuilt by GDAL on the first read
Sys.setenv(SHAPE_RESTORE_SHX = "YES")

### range polygons (WGS84 where the CRS is missing), dissolved into one
range_sf <- st_read(input("^Quercus_robur_plg.*\\.shp$"), quiet = TRUE)
if (is.na(st_crs(range_sf))) st_crs(range_sf) <- 4326

range_sf <- range_sf %>%
  st_make_valid() %>% st_transform(crs_map) %>% st_union() %>% st_as_sf()

range_for <- function(view) suppressWarnings(st_crop(range_sf, view))



# sites ------------------------------------------------------------------------

### samples without coordinates (not on the map)
no_coord <- meta %>%
  filter(project != "outgroup", is.na(latitude) | is.na(longitude))

if (nrow(no_coord) > 0) {
  message(nrow(no_coord), " samples without coordinates, not on the map: ",
          paste(sort(unique(no_coord$sample)), collapse = ", "))
}

### sites in the call set and samples per site (from the metadata)
site_n <- meta %>%
  filter(project != "outgroup", !is.na(latitude), !is.na(longitude)) %>%
  group_by(project, site, site_name, country) %>%
  summarise(n = n(), lat = first(latitude), lon = first(longitude),
            .groups = "drop")

## check if every site has one coordinate
spread <- meta %>%
  filter(!is.na(latitude)) %>%
  group_by(site) %>%
  summarise(n_coord = n_distinct(paste(latitude, longitude)), .groups = "drop") %>%
  filter(n_coord > 1)

if (nrow(spread) > 0) {
  warning("sites with more than one coordinate: ",
          paste(spread$site, collapse = ", "))
}

sites_sf <- st_as_sf(site_n, coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs_map)



# lineages per site ------------------------------------------------------------

### one per genotype and site (clonemates are copies of one genotype)
lin_site <- lin %>%
  select(sample, mlg, lineage = all_of(k_map)) %>%
  inner_join(select(meta, sample, site, latitude, longitude), by = "sample") %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  distinct(site, mlg, .keep_all = TRUE) %>%
  group_by(site) %>%
  mutate(lat = first(latitude), lon = first(longitude)) %>%
  ungroup() %>%
  count(site, lat, lon, lineage, name = "n")

### labelled samples without coordinates (not in the pies)
no_place <- anti_join(lin, filter(meta, !is.na(latitude), !is.na(longitude)),
                      by = "sample")

if (nrow(no_place) > 0) {
  message(nrow(no_place), " labelled samples without coordinates, not in the pies")
}

## check if every lineage has a colour
if (any(!lin_site$lineage %in% names(pal_unit))) {
  warning("lineages without a colour in pal_unit: ",
          paste(setdiff(lin_site$lineage, names(pal_unit)), collapse = ", "))
}

### one point per site, in the map projection
lin_sites <- lin_site %>%
  group_by(site, lat, lon) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs_map)

lin_sites <- lin_sites %>%
  mutate(x     = st_coordinates(lin_sites)[, 1],
         y     = st_coordinates(lin_sites)[, 2],
         in_ch = lengths(st_within(lin_sites, st_as_sfc(view_ch))) > 0)

### pie wedges as polygons in the map projection (clockwise from north, scaled with the map)
pie_wedges <- function(x, y, r, n, lineage, n_arc = 72) {
  f   <- n / sum(n)
  ang <- 2 * pi * c(0, cumsum(f))
  rings <- map(seq_along(n), function(i) {
    th  <- seq(ang[i], ang[i + 1], length.out = max(3, ceiling(n_arc * f[i])))
    arc <- cbind(x + r * sin(th), y + r * cos(th))
    if (f[i] == 1) rbind(arc, arc[1, ]) else rbind(c(x, y), arc, c(x, y))
  })
  st_sf(lineage  = lineage,
        geometry = st_sfc(map(rings, ~ st_polygon(list(.x))), crs = crs_map))
}

### pies of a set of sites, radius by the square root of the genotypes
pies_for <- function(sites, r) {
  n_max <- max(lin_sites$n)
  map(seq_len(nrow(sites)), function(i) {
    d <- filter(lin_site, site == sites$site[i])
    pie_wedges(sites$x[i], sites$y[i],
               r[["min"]] + (r[["max"]] - r[["min"]]) * sqrt(sites$n[i] / n_max),
               d$n, d$lineage)
  }) %>%
    do.call(rbind, .)
}



# figures ----------------------------------------------------------------------

### map builder (by_project: project markers and legend; zones: zone raster; pies: pies in place of markers)
draw_map <- function(pts, view, zone_df = NULL, by_project = FALSE,
                     zones = TRUE, label_size = 2.1, point_size = 2.6,
                     scalebar = TRUE, pies = NULL) {

  p <- ggplot() +
    geom_sf(data = land_for(view), fill = "grey92", colour = NA)

  if (zones) {
    p <- p +
      geom_raster(data = zone_df, aes(x, y, fill = zone)) +
      scale_fill_manual(values = zone_colours, na.value = "transparent",
                        na.translate = FALSE, name = "Environmental zone")
  }

  ### borders and the range hatching (no fill, the zones stay visible)
  p <- p +
    geom_sf(data = land_for(view), fill = NA, colour = "grey45", linewidth = 0.15) +
    geom_sf_pattern(data = range_for(view), aes(pattern = "range"),
                    fill = NA, colour = "grey10", linewidth = 0.35,
                    pattern_fill = "grey15", pattern_colour = NA,
                    pattern_angle = 45, pattern_density = 0.06,
                    pattern_spacing = 0.016, pattern_key_scale_factor = 0.5) +
    scale_pattern_manual(values = c(range = "stripe"), name = NULL,
                         labels = expression(italic("Quercus robur")~"range"))

  ### site markers: white fill, dark ring (sample count in the label); pies in their place where given
  if (!is.null(pies)) {
    p <- p +
      geom_sf(data = pies, aes(fill = lineage), colour = "white", linewidth = 0.1) +
      scale_fill_manual(values = pal_unit, name = "Lineage")
  } else if (by_project) {
    p <- p +
      geom_sf(data = pts, aes(colour = project), shape = 21, fill = "white",
              stroke = 0.9, size = point_size) +
      scale_colour_manual(values = col_project, name = "Project",
                          labels = c(magic = "MAGIC",
                                     phytoakmeter = "PhytOakmeter"))
  } else {
    p <- p +
      geom_sf(data = pts, aes(shape = "Sampling site"), fill = "white",
              colour = "black", stroke = 0.5, size = point_size) +
      scale_shape_manual(values = c("Sampling site" = 21), name = NULL)
  }

  ### site labels (fixed seed for a stable layout, white halo over dark zones)
  p <- p +
    geom_text_repel(data = pts,
                    aes(geometry = geometry,
                        label = sprintf("%s (%d)", site, n)),
                    stat = "sf_coordinates", size = label_size, colour = "black",
                    bg.color = "white", bg.r = 0.12, segment.size = 0.2,
                    segment.colour = "grey30", min.segment.length = 0,
                    box.padding = 0.35, point.padding = 0.2,
                    max.overlaps = Inf, seed = 42)

  ### scale bar (orientation from the labelled graticule)
  if (scalebar) {
    p <- p + annotation_scale(location = "bl", width_hint = 0.22, style = "bar",
                              bar_cols = c("grey15", "white"),
                              height = unit(0.15, "cm"),
                              text_cex = 0.6, line_width = 0.4)
  }

  p +
    coord_sf(xlim = view[c("xmin", "xmax")], ylim = view[c("ymin", "ymax")],
             crs = st_crs(crs_map), datum = st_crs(4326), expand = FALSE) +
    guides(shape   = guide_legend(order = 1, override.aes = list(size = 2.8)),
           colour  = guide_legend(order = 1, override.aes = list(size = 2.8)),
           pattern = guide_legend(order = 2),
           fill    = guide_legend(order = 3, override.aes = list(colour = NA))) +
    labs(x = NULL, y = NULL) +
    theme_pubr(base_size = 10, border = TRUE, legend = "right") +
    theme(
      ### sea
      panel.background = element_rect(fill = "#EAF1F6", colour = NA),
      panel.grid.major = element_line(colour = "grey80", linewidth = 0.15),
      panel.border     = element_rect(colour = "grey20", fill = NA,
                                      linewidth = 0.4),
      axis.text        = element_text(size = 6.5, colour = "grey30"),
      ### legend keys without background (hatching visible)
      legend.key       = element_blank(),
      legend.title     = element_text(size = 8, face = "bold"),
      legend.text      = element_text(size = 7),
      legend.key.size  = unit(0.35, "cm"),
      plot.margin      = margin(3, 3, 3, 3)
    )
}


### one map per project, without the zone raster
p_phyt <- draw_map(filter(sites_sf, project == "phytoakmeter"), view_eu,
                   zones = FALSE)

p_magic <- draw_map(filter(sites_sf, project == "magic"), view_ch,
                    zones = FALSE, label_size = 2.4)

save_map("Fig_sampling_map_phytoakmeter", p_phyt)
save_map("Fig_sampling_map_magic", p_magic)


### both projects with the zones, the Swiss window as a box and as an inset
box_ch <- st_as_sfc(view_ch)

p_both <- draw_map(sites_sf, view_eu, zone_eu, by_project = TRUE) +
  geom_sf(data = box_ch, fill = NA, colour = "grey10", linewidth = 0.4,
          linetype = "dashed") +
  coord_sf(xlim = view_eu[c("xmin", "xmax")], ylim = view_eu[c("ymin", "ymax")],
           crs = st_crs(crs_map), datum = st_crs(4326), expand = FALSE)

p_inset <- draw_map(filter(sites_sf, project == "magic"), view_ch, zone_ch,
                    label_size = 1.7, point_size = 1.6, scalebar = FALSE) +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = margin(0, 0, 0, 0))

p_combined <- p_both +
  inset_element(p_inset, left = 0.02, bottom = 0.52, right = 0.40, top = 0.98,
                align_to = "panel")

save_map("Fig_sampling_map_combined", p_combined, width = 7.1, height = 5.9)


### lineage pies of the genotypes per site, without the zones; the Swiss sites in the inset only
p_lin <- draw_map(filter(lin_sites, !in_ch), view_eu, zones = FALSE,
                  pies = pies_for(filter(lin_sites, !in_ch), pie_eu)) +
  geom_sf(data = box_ch, fill = NA, colour = "grey10", linewidth = 0.4,
          linetype = "dashed") +
  coord_sf(xlim = view_eu[c("xmin", "xmax")], ylim = view_eu[c("ymin", "ymax")],
           crs = st_crs(crs_map), datum = st_crs(4326), expand = FALSE)

p_lin_ch <- draw_map(filter(lin_sites, in_ch), view_ch, zones = FALSE,
                     pies = pies_for(filter(lin_sites, in_ch), pie_ch),
                     label_size = 2.4)

p_lin_inset <- draw_map(filter(lin_sites, in_ch), view_ch, zones = FALSE,
                        pies = pies_for(filter(lin_sites, in_ch), pie_ch),
                        label_size = 1.7, scalebar = FALSE) +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        plot.margin = margin(0, 0, 0, 0))

p_lin_combined <- p_lin +
  inset_element(p_lin_inset, left = 0.02, bottom = 0.52, right = 0.40, top = 0.98,
                align_to = "panel")

save_map("Fig_lineage_map_combined", p_lin_combined, width = 7.1, height = 5.9)
save_map("Fig_lineage_map_ch", p_lin_ch)



# tables -----------------------------------------------------------------------

write_tsv(st_drop_geometry(sites_sf), file.path(out_root, "sites_sampled.tsv"))
write_tsv(lin_site, file.path(out_root, "site_lineages.tsv"))

cat("\n=== genotypes per site and lineage (", k_map, ") ===\n", sep = "")
lin_site %>%
  select(site, lineage, n) %>%
  pivot_wider(names_from = lineage, values_from = n, values_fill = 0) %>%
  print(n = Inf)

cat("\n=== samples per site ===\n")
site_n %>%
  arrange(project, desc(n)) %>%
  print(n = Inf)

cat("\nsites                ", nrow(site_n), "\n")
cat("samples              ", sum(site_n$n), "\n")
cat("sites in MAGIC       ", sum(site_n$project == "magic"), "\n")
cat("sites in PhytOakmeter", sum(site_n$project == "phytoakmeter"), "\n")



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  Fig_sampling_map_phytoakmeter.pdf\n")
cat("  Fig_sampling_map_phytoakmeter.png\n")
cat("  Fig_sampling_map_magic.pdf\n")
cat("  Fig_sampling_map_magic.png\n")
cat("  Fig_sampling_map_combined.pdf\n")
cat("  Fig_sampling_map_combined.png\n")
cat("  Fig_lineage_map_combined.pdf\n")
cat("  Fig_lineage_map_combined.png\n")
cat("  Fig_lineage_map_ch.pdf\n")
cat("  Fig_lineage_map_ch.png\n")
cat("  sites_sampled.tsv\n")
cat("  site_lineages.tsv\n")
