## ============================================================================
## 041_climate_data.R
## Climate at all sampling sites (CHELSA V2.1, 1981-2010), with project and host per site;
## climate candidates screened for correlation in order of priority
##
## Files needed:
##   sample_metadata_ALL_FIXED.tsv        site, project, host and coordinates per sample,
##                                        from local/01_variant_calling/103_mating_type.R
##   CHELSA_<var>_1981-2010_V.2.1.tif     CHELSA V2.1 bio1, bio4, bio12 and cmi_mean,
##                                        external: downloaded by this script (once) into
##                                        local/04_diversity_niche_selection/supplementary/041_supplementary_rasters/chelsa
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)
library(terra)



# files ------------------------------------------------------------------------

meta_file  <- "01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv"
raster_dir <- "04_diversity_niche_selection/supplementary/041_supplementary_rasters"
out_root   <- "04_diversity_niche_selection/results/041_results_climate_data"

dir.create(file.path(raster_dir, "chelsa"), recursive = TRUE, showWarnings = FALSE)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### CHELSA V2.1 climatologies 1981-2010 (bioclim and BIOCLIM+ files share one folder)
chelsa_base <- "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio"

### climate candidates in order of priority, scale and offset of the stored integers (CHELSA V2.1 technical
### specification; applied only if the file does not carry them)
clim_vars <- tribble(
  ~var,       ~label,                      ~unit,         ~scale, ~offset,
  "bio1",     "Mean annual temperature",   "°C",          0.1,    -273.15,
  "bio12",    "Annual precipitation",      "mm",          0.1,    0,
  "bio4",     "Temperature seasonality",   "°C x 100",    0.1,    0,
  "cmi_mean", "Climate moisture index",    "mm month⁻¹",  0.1,    0)

### correlation threshold, buffer for sites on a missing pixel (m), download timeout (s)
r_max     <- 0.7
buffer_m  <- 2000
timeout_s <- 3600

### site groups (project and host)
lab_group   <- c(phytoakmeter_Quercus_robur = "PhytOakmeter, oak", magic_Quercus_robur = "MAGIC, oak",
                 magic_Fagus_sylvatica = "MAGIC, beech")



# style ------------------------------------------------------------------------

### site group shapes
shape_group <- c(phytoakmeter_Quercus_robur = 16, magic_Quercus_robur = 17, magic_Fagus_sylvatica = 2)

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

meta <- read_tsv(meta_file, col_types = cols(.default = "c"))

### sampling sites of the ingroup isolates with coordinates, project and host (the first species of dominant_host)
sites <- meta %>%
  filter(in_ingroup == "TRUE", !is.na(latitude), !is.na(longitude)) %>%
  mutate(across(c(latitude, longitude), as.numeric), host = str_trim(str_remove(dominant_host, ";.*"))) %>%
  summarise(lat = mean(latitude), lon = mean(longitude),
            spread = max(diff(range(latitude)), diff(range(longitude))), isolates = n(),
            hosts = n_distinct(host), host = first(host), .by = c(project, site)) %>%
  mutate(group = paste(project, host, sep = "_")) %>%
  arrange(project, site)



# checks -----------------------------------------------------------------------

## check if every site has one coordinate
if (any(sites$spread > 0.01)) {
  warning("sites with isolates more than 0.01 degrees apart (mean coordinate used): ",
          paste(sites$site[sites$spread > 0.01], collapse = ", "))
}

## check if every site belongs to one project
if (anyDuplicated(sites$site)) {
  stop("sites in more than one project: ", paste(sites$site[duplicated(sites$site)], collapse = ", "))
}

## check if every site has one host
if (any(sites$hosts > 1) || anyNA(sites$host)) {
  stop("sites without exactly one host: ", paste(sites$site[sites$hosts != 1 | is.na(sites$host)], collapse = ", "))
}

## check if every site falls into a defined site group
if (!all(sites$group %in% names(lab_group))) {
  stop("site groups without a label: ", paste(setdiff(sites$group, names(lab_group)), collapse = ", "))
}



# downloads --------------------------------------------------------------------

### CHELSA files (downloaded once into raster_dir)
options(timeout = timeout_s)
chelsa_file <- function(v) file.path(raster_dir, "chelsa", sprintf("CHELSA_%s_1981-2010_V.2.1.tif", v))
for (v in clim_vars$var) {
  if (!file.exists(chelsa_file(v))) {
    url <- sprintf("%s/CHELSA_%s_1981-2010_V.2.1.tif", chelsa_base, v)
    ok  <- tryCatch(download.file(url, chelsa_file(v), mode = "wb", quiet = TRUE), error = function(e) 1)
    if (ok != 0) {
      unlink(chelsa_file(v))
      stop("CHELSA download failed for ", v, " (", url, "); update chelsa_base if the files have moved")
    }
  }
}



# extraction -------------------------------------------------------------------

### value at the site, or the mean within buffer_m where the site falls on a missing pixel
extract_sites <- function(r, d) {
  p <- vect(d, geom = c("lon", "lat"), crs = "EPSG:4326")
  x <- terra::extract(r, p, ID = FALSE)[[1]]
  miss <- is.na(x)
  if (any(miss)) {
    x[miss] <- terra::extract(r, buffer(p[miss], buffer_m), fun = mean, na.rm = TRUE, ID = FALSE)[[1]]
  }
  x
}

### climate at every site (scale and offset applied when the file stores raw integers, temperature in 0.1 K)
clim_raw <- map(setNames(clim_vars$var, clim_vars$var), ~ extract_sites(rast(chelsa_file(.x)), sites))
raw <- median(clim_raw$bio1, na.rm = TRUE) > 150
site_env <- sites %>%
  select(project, host, group, site, lat, lon, isolates) %>%
  bind_cols(as_tibble(map2(clim_raw, seq_along(clim_raw), function(x, i) {
    if (raw) x * clim_vars$scale[i] + clim_vars$offset[i] else x
  })))



# checks on the extracted values -----------------------------------------------

## check if every site has climate values
if (anyNA(site_env[clim_vars$var])) {
  stop("sites without climate values: ",
       paste(site_env$site[!complete.cases(site_env[clim_vars$var])], collapse = ", "))
}

## check if the temperatures are plausible after scaling
if (any(site_env$bio1 < -10 | site_env$bio1 > 30)) {
  stop("mean annual temperature outside -10 to 30 °C; check the CHELSA scale and offset")
}



# climate variables ------------------------------------------------------------

### correlations across all sites; candidates kept in order of priority unless |r| > r_max with a kept one
clim_cor <- cor(as.matrix(site_env[clim_vars$var]))
clim_keep <- character(0)
for (v in clim_vars$var) {
  if (length(clim_keep) == 0 || all(abs(clim_cor[v, clim_keep]) <= r_max)) clim_keep <- c(clim_keep, v)
}

clim_sel <- clim_vars %>%
  mutate(priority = row_number(), kept = var %in% clim_keep, scaled_here = raw,
         max_r_kept = map_dbl(var, function(v) {
           others <- setdiff(clim_keep, v)
           if (length(others) == 0) NA_real_ else max(abs(clim_cor[v, others]))
         }))



# figures ----------------------------------------------------------------------

### climate variable selection: correlations across sites in order of priority, kept variables outlined
cor_long <- as_tibble(clim_cor, rownames = "v1") %>%
  pivot_longer(-v1, names_to = "v2", values_to = "r") %>%
  mutate(kept = v1 == v2 & v1 %in% clim_keep,
         v1 = factor(v1, levels = clim_vars$var), v2 = factor(v2, levels = rev(clim_vars$var)))

p_clim <- ggplot(cor_long, aes(v1, v2, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_tile(data = filter(cor_long, kept), fill = NA, colour = "black", linewidth = 0.5, width = 1, height = 1) +
  geom_text(aes(label = sprintf("%.2f", r)), size = 6 / .pt) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1), name = "r") +
  coord_equal() +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "right", legend.key.height = unit(4, "mm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

### climate of the sites by site group (first two kept variables)
p_sites <- site_env %>%
  ggplot(aes(x = .data[[clim_keep[1]]], y = .data[[clim_keep[min(2, length(clim_keep))]]], shape = group)) +
  geom_point(size = 1.6, colour = "black") +
  scale_shape_manual(values = shape_group, labels = lab_group, breaks = names(lab_group), name = NULL) +
  labs(x = sprintf("%s (%s)", clim_vars$label[clim_vars$var == clim_keep[1]],
                   clim_vars$unit[clim_vars$var == clim_keep[1]]),
       y = sprintf("%s (%s)", clim_vars$label[clim_vars$var == clim_keep[min(2, length(clim_keep))]],
                   clim_vars$unit[clim_vars$var == clim_keep[min(2, length(clim_keep))]]))

p_climate <- (p_clim | p_sites) + plot_annotation(tag_levels = "A")
save_fig("climate_variable_selection.pdf", p_climate, 7.1, 3)



# tables -----------------------------------------------------------------------

### climate candidates in order of priority
pub_clim <- clim_sel %>%
  transmute(Priority = priority, Variable = sprintf("%s (%s)", var, label), Unit = unit,
            Kept = if_else(kept, "yes", "no"),
            `Largest |r| with a kept variable` = if_else(is.na(max_r_kept), "-", sprintf("%.2f", max_r_kept)))

### climate per site
pub_sites <- site_env %>%
  transmute(`Site group` = lab_group[group], Site = site, Latitude = round(lat, 4), Longitude = round(lon, 4),
            Isolates = isolates, !!!setNames(map(clim_vars$var, ~ round(site_env[[.x]], 1)), clim_vars$var))

write_tsv(site_env, file.path(out_root, "site_climate.tsv"))
write_tsv(clim_sel, file.path(out_root, "climate_variables.tsv"))
write_tsv(as_tibble(clim_cor, rownames = "var"), file.path(out_root, "climate_correlations.tsv"))
write_tsv(pub_clim, file.path(out_root, "publication_climate_variables.tsv"))
write_tsv(pub_sites, file.path(out_root, "publication_site_climate.tsv"))



# summary ----------------------------------------------------------------------

cat("\nCHELSA values  ",
    if (raw) "scaled from the stored integers" else "read with the scale and offset of the files", "\n")

cat("\n=== Climate candidates ===\n")
print(pub_clim, width = Inf)
print(round(clim_cor, 2))

cat("\n=== Site climate ===\n")
print(pub_sites, n = Inf, width = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  climate_variable_selection.pdf\n")
cat("  site_climate.tsv\n")
cat("  climate_variables.tsv\n")
cat("  climate_correlations.tsv\n")
cat("  publication_climate_variables.tsv\n")
cat("  publication_site_climate.tsv\n")
cat("\nCHELSA rasters in", file.path(raster_dir, "chelsa"), "\n")
cat("  CHELSA_<var>_1981-2010_V.2.1.tif\n")
