## ============================================================================
## 042_climate_niche.R
## Climate, project and host of the units at K = 4, 5 and 6 across all sampling
## sites (partial RDA, variance partitioning, host within MAGIC, composition per
## site group), and soil, depth, climate and site context at the sites with in
## situ soil data
##
## Files needed:
##   site_climate.tsv               project, host and climate per site
##   climate_variables.tsv          climate candidates and the kept ones
##                                  (both from local/04_diversity_niche_selection/041_climate_data.R)
##   soildata_ALL_FIXED.txt         external: in situ soil data per PhytOakmeter site and depth and per
##                                  MAGIC site (copy manually into local/04_diversity_niche_selection/supplementary)
##   lineage_labels.tsv             unit at K = 4, 5 and 6 per ingroup sample,
##                                  from local/03_delimitation/001_lineage_labels.R
##   sample_metadata_ALL_FIXED.tsv  per-sample table (site, sampling_depth, mlg, mat_genome, in_ingroup,
##                                  latitude), from local/01_variant_calling/103_mating_type.R
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)
library(ggrepel)
library(vegan)



# files ------------------------------------------------------------------------

env_root    <- "04_diversity_niche_selection/results/041_results_climate_data"
soil_file   <- "04_diversity_niche_selection/supplementary/soildata_ALL_FIXED.txt"
labels_file <- "03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv"
meta_file   <- "01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv"
out_root    <- "04_diversity_niche_selection/results/042_results_climate_niche"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### levels and sample sets
k_levels <- c("k4", "k5", "k6")
sets     <- c("all", "no_mixed")

### depth code of the single MAGIC sample (0-10 cm, empty in the metadata)
depth_single <- c(magic = "1")

### mating type calls of possible mixed cultures (sample set without them)
mat_exclude <- c("MAT12", "ambiguous")

### kept in situ soil variables of 040 (column, name, ln-transformed)
soil_main <- tribble(
  ~column,     ~var,   ~log,
  "pH_H20",    "pH",   FALSE,
  "C",         "C",    TRUE,
  "C_N",       "CN",   FALSE,
  "Prop_Clay", "Clay", FALSE,
  "P2O5",      "P2O5", TRUE)

### short labels of the climate variables (arrows)
lab_clim <- c(bio1 = "Temperature", bio12 = "Precipitation", bio4 = "Seasonality", cmi_mean = "Moisture")

### permutations, seed
n_perm <- 999
seed   <- 42



# style ------------------------------------------------------------------------

### unit colours
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

### site groups (project and host): labels and shapes; level fills (variance partitioning)
lab_group   <- c(phytoakmeter_Quercus_robur = "PhytOakmeter, oak", magic_Quercus_robur = "MAGIC, oak",
                 magic_Fagus_sylvatica = "MAGIC, beech")
shape_group <- c(phytoakmeter_Quercus_robur = 16, magic_Quercus_robur = 17, magic_Fagus_sylvatica = 2)
pal_level   <- c(k4 = "grey80", k5 = "grey55", k6 = "grey25")

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}

### p values in tables
fmt_p <- function(p) if_else(p < 0.001, "< 0.001", sprintf("%.3f", p))

### model fit for the triplot corner (adjusted R2 and p of the overall test), level labels
lab_fit <- function(r2, p) sprintf("Adj. R² = %.1f%%, %s", 100 * r2,
                                   if_else(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))
lab_k <- function(k) sub("k", "K = ", k)



# read -------------------------------------------------------------------------

site_env  <- read_tsv(file.path(env_root, "site_climate.tsv"),
                      col_types = cols(site = "c", project = "c", host = "c", group = "c", .default = "d"))
clim_vars <- read_tsv(file.path(env_root, "climate_variables.tsv"), show_col_types = FALSE)
clim_kept <- clim_vars$var[clim_vars$kept]

### in situ soil per site and depth (site codes with "-" as in the metadata), kept variables (ln where flagged)
soil <- read_tsv(soil_file, col_types = cols(Site = "c", .default = "d")) %>%
  mutate(site = str_replace_all(Site, "_", "-"), depth = as.character(Depth))

soil_x <- soil %>%
  transmute(site, depth,
            !!!setNames(map2(soil_main$column, soil_main$log,
                             function(cl, lg) if (lg) log(soil[[cl]]) else soil[[cl]]),
                        soil_main$var))

lab  <- read_tsv(labels_file, col_types = cols(.default = "c"))
meta <- read_tsv(meta_file, col_types = cols(.default = "c"))

### ingroup isolates with site, project, host, depth and MLG
iso <- lab %>%
  select(sample, all_of(k_levels)) %>%
  inner_join(meta %>% select(sample, site, depth = sampling_depth, mlg, mat_genome), by = "sample") %>%
  inner_join(select(site_env, site, project, host), by = "site") %>%
  mutate(depth = coalesce(depth, depth_single[project]), mixed = mat_genome %in% mat_exclude)



# checks -----------------------------------------------------------------------

## check if climate variables were kept
if (length(clim_kept) == 0) {
  stop("no kept climate variables in climate_variables.tsv")
}

## check if every isolate site has climate
if (!all(unique(meta$site[meta$in_ingroup == "TRUE" & !is.na(meta$latitude)]) %in% site_env$site)) {
  warning("sites with isolates but without climate are left out")
}

## check if both hosts occur within a project (host and project are separable)
if (!any(count(distinct(site_env, project, host), project)$n > 1)) {
  stop("no project with both hosts; check dominant_host in the metadata")
}

## check if the kept soil variables are complete and every soil site has a project
if (anyNA(soil_x[soil_main$var]) || !all(soil_x$site %in% site_env$site)) {
  stop("missing soil values or soil sites without project in site_climate.tsv")
}

## check if every isolate at a site with soil data has a depth
if (anyNA(iso$depth[iso$site %in% soil_x$site])) {
  stop("isolates at sites with soil data without a sampling depth")
}



# data -------------------------------------------------------------------------

### MLGs per unit per site (all sites) or per site and depth (sites with soil data), for a level and sample set
comp <- function(kk, ss, by) {
  d <- if (ss == "all") iso else filter(iso, !mixed)
  tibble(site = d$site, depth = if (by == "site") "all" else d$depth, mlg = d$mlg, unit = d[[kk]]) %>%
    distinct() %>%
    count(site, depth, unit) %>%
    pivot_wider(names_from = unit, values_from = n, values_fill = 0) %>%
    left_join(site_env, by = "site") %>%
    arrange(site, depth)
}

### dummy columns of a factor (none when it has one level)
dummies <- function(x) {
  if (n_distinct(x) > 1) model.matrix(~ x)[, -1, drop = FALSE] else matrix(0, length(x), 0)
}

### standardisation within project and one column per variable and project (zero outside the project), as in 040
z_within <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))
by_project <- function(X, g) {
  out <- do.call(cbind, map(sort(unique(g)), function(p) {
    x <- X
    x[g != p, ] <- 0
    colnames(x) <- paste(colnames(X), p, sep = "@")
    x
  }))
  out[, order(match(sub("@.*", "", colnames(out)), colnames(X))), drop = FALSE]
}

### test of a partial RDA (response R, explanatory X, conditions Z)
test_rda <- function(R, X, Z, perm) {
  m <- rda(R, X, Z)
  a <- anova(m, permutations = perm)
  tibble(adj_r2 = RsquareAdj(m)$adj.r.squared, F = a$F[1], p = a$`Pr(>F)`[1])
}

### unique fractions [a] to [d] of a four-table vegan variance partitioning, the rest shared
four_parts <- function(fr, keys) {
  part <- case_when(str_detect(rownames(fr), "^\\[a\\]") ~ keys[1],
                    str_detect(rownames(fr), "^\\[b\\]") ~ keys[2],
                    str_detect(rownames(fr), "^\\[c\\]") ~ keys[3],
                    str_detect(rownames(fr), "^\\[d\\]") ~ keys[4],
                    str_detect(rownames(fr), "^\\[p\\]") ~ "residual",
                    .default = "shared")
  tibble(fraction = rownames(fr), part = part, adj_r2 = fr$Adj.R.square)
}



# climate, geography, project and host across all sites ------------------------

### climate given geography, project and host: overall and marginal tests, VIF; variance partitioning into
### climate, geography, project and host with tests of the unique fractions; host within the MAGIC sites,
### given geography
fit_sites <- function(kk, ss) {
  dat   <- comp(kk, ss, "site")
  units <- intersect(names(pal_unit), names(dat))
  R <- decostand(as.matrix(dat[units]), "hellinger")
  C <- scale(as.matrix(dat[clim_kept]))
  G <- scale(as.matrix(dat[c("lat", "lon")]))
  P <- dummies(dat$project)
  H <- dummies(dat$host)
  free <- how(nperm = n_perm)

  m    <- rda(R, C, cbind(G, P, H))
  a    <- test_rda(R, C, cbind(G, P, H), free)
  marg <- map_dfr(clim_kept, ~ test_rda(R, C[, .x, drop = FALSE],
                                        cbind(G, P, H, C[, setdiff(clim_kept, .x), drop = FALSE]), free))
  vif  <- vif.cca(m)[clim_kept]

  vp     <- varpart(R, C, G, P, H)
  t_geo  <- test_rda(R, G, cbind(C, P, H), free)
  t_proj <- test_rda(R, P, cbind(C, G, H), free)
  t_host <- test_rda(R, H, cbind(C, G, P), free)

  mg   <- dat$project == "magic"
  h_mg <- test_rda(R[mg, , drop = FALSE], dummies(dat$host[mg]), scale(as.matrix(dat[mg, c("lat", "lon")])), free)

  list(test = tibble(k = kk, set = ss, sites = nrow(dat), mlg = sum(dat[units]), adj_r2 = a$adj_r2, F = a$F, p = a$p,
                     vif_max = max(vif)),
       margin = tibble(k = kk, set = ss, var = clim_kept, F = marg$F, p = marg$p, vif = vif),
       varpart = four_parts(vp$part$indfract, c("climate_only", "geography_only", "project_only", "host_only")) %>%
         mutate(k = kk, set = ss, .before = 1) %>%
         mutate(p = case_when(part == "climate_only" ~ a$p, part == "geography_only" ~ t_geo$p,
                              part == "project_only" ~ t_proj$p, part == "host_only" ~ t_host$p,
                              .default = NA_real_)),
       host_magic = tibble(k = kk, set = ss, sites = sum(mg), sites_oak = sum(mg & dat$host == "Quercus_robur"),
                           adj_r2 = h_mg$adj_r2, F = h_mg$F, p = h_mg$p),
       model = m, group = dat$group)
}

set.seed(seed)
fits <- expand_grid(k = k_levels, set = sets) %>% mutate(fit = map2(k, set, fit_sites))

tests      <- map_dfr(fits$fit, "test")
margins    <- map_dfr(fits$fit, "margin")
vparts     <- map_dfr(fits$fit, "varpart")
host_magic <- map_dfr(fits$fit, "host_magic")

vparts_sum <- vparts %>%
  summarise(adj_r2 = sum(adj_r2), p = first(na.omit(p), default = NA_real_), .by = c(k, set, part))

### share of every unit among the MLGs of a site group (MLGs counted once per site; all isolates)
groups <- map_dfr(k_levels, function(kk) {
  d <- distinct(iso, group = paste(project, host, sep = "_"), site, mlg, unit = .data[[kk]])
  d %>%
    count(group, unit) %>%
    complete(group, unit, fill = list(n = 0)) %>%
    mutate(share = n / sum(n), mlgs = sum(n), .by = group) %>%
    left_join(summarise(d, sites = n_distinct(site), .by = group), by = "group") %>%
    mutate(k = kk, .before = 1)
})



# soil, depth, climate and site context at the sites with soil data ------------

### variance partitioning of the composition per site and depth: soil standardised within project with
### responses per project (as in 040), site context = geography, project and host; unique soil fraction tested
### by free permutation (climate varies between sites only and is tested at site level above)
fit_four <- function(kk, ss) {
  dat   <- comp(kk, ss, "site_depth") %>% inner_join(soil_x, by = c("site", "depth"))
  units <- intersect(names(pal_unit), names(dat))
  R <- decostand(as.matrix(dat[units]), "hellinger")
  S <- by_project(sapply(soil_main$var, function(v) z_within(dat[[v]], dat$project)), dat$project)
  D <- dummies(dat$depth)
  C <- scale(as.matrix(dat[clim_kept]))
  X <- cbind(scale(as.matrix(dat[c("lat", "lon")])), dummies(dat$project), dummies(dat$host))
  vp <- varpart(R, S, D, C, X)
  t_soil <- test_rda(R, S, cbind(D, C, X), how(nperm = n_perm))
  four_parts(vp$part$indfract, c("soil_only", "depth_only", "climate_only", "context_only")) %>%
    mutate(k = kk, set = ss, cells = nrow(dat), sites = n_distinct(dat$site), .before = 1) %>%
    mutate(p = if_else(part == "soil_only", t_soil$p, NA_real_))
}

set.seed(seed)
four <- map2_dfr(rep(k_levels, each = length(sets)), rep(sets, length(k_levels)), fit_four)

four_sum <- four %>%
  summarise(adj_r2 = sum(adj_r2), p = first(na.omit(p), default = NA_real_), .by = c(k, set, cells, sites, part))



# ordination scores ------------------------------------------------------------

### triplot scores per level (all isolates): sites, units, climate variables
xy <- function(x, id) {
  x <- as.matrix(x)
  if (ncol(x) == 1) x <- cbind(x, 0)
  colnames(x) <- c("A1", "A2")
  as_tibble(x, rownames = id)
}

scores_tbl <- map_dfr(k_levels, function(kk) {
  f   <- fits$fit[[which(fits$k == kk & fits$set == "all")]]
  m   <- f$model
  eig <- c(m$CCA$eig, m$CA$eig)
  ax  <- sprintf("%s (%.1f%%)", names(eig)[1:2], 100 * eig[1:2] / m$tot.chi)
  st  <- xy(scores(m, display = "sites", choices = 1:2, scaling = 2), "id") %>% mutate(type = "site", group = f$group)
  sp  <- xy(scores(m, display = "species", choices = 1:2, scaling = 2), "id") %>% mutate(type = "unit")
  bp  <- xy(scores(m, display = "bp", choices = 1:2, scaling = 2), "id") %>% mutate(type = "variable")
  bp  <- mutate(bp, across(c(A1, A2), ~ .x * 0.9 * max(abs(c(st$A1, st$A2))) / max(abs(c(bp$A1, bp$A2)))))
  bind_rows(st, sp, bp) %>% mutate(k = kk, axis1 = ax[1], axis2 = ax[2],
                                   fit = lab_fit(f$test$adj_r2, f$test$p), .before = 1)
}) %>%
  mutate(label = if_else(type == "variable", coalesce(lab_clim[id], id), id))



# figures ----------------------------------------------------------------------

### triplot per level (all isolates): sites by site group, climate arrows, unit scores
plot_triplot <- function(kk) {
  d  <- filter(scores_tbl, k == kk)
  st <- filter(d, type == "site")
  sp <- filter(d, type == "unit")
  bp <- filter(d, type == "variable")
  ggplot() +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_point(data = st, aes(A1, A2, shape = group), colour = "grey55", size = 1.2) +
    geom_segment(data = bp, aes(0, 0, xend = A1, yend = A2), arrow = arrow(length = unit(1.5, "mm")),
                 linewidth = 0.35, colour = "grey15") +
    geom_text_repel(data = bp, aes(A1, A2, label = label), size = 6.5 / .pt, colour = "black", seed = 1) +
    geom_point(data = sp, aes(A1, A2, fill = id), shape = 21, size = 2.6, stroke = 0.4, colour = "black") +
    geom_text_repel(data = sp, aes(A1, A2, label = id), size = 7 / .pt, fontface = "bold", seed = 2,
                    box.padding = 0.4) +
    scale_fill_manual(values = pal_unit, guide = "none") +
    scale_shape_manual(values = shape_group, labels = lab_group, breaks = names(lab_group), name = NULL) +
    annotate("text", x = -Inf, y = Inf, label = d$fit[1], hjust = -0.05, vjust = 1.4, size = 6.5 / .pt) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    labs(x = d$axis1[1], y = d$axis2[1])
}

### variance partitioning at the sites with soil data per level (all isolates)
lab_frac <- c(soil_only = "Soil only", depth_only = "Depth only", climate_only = "Climate only",
              context_only = "Site context only", shared = "Shared")
p_vp <- four_sum %>%
  filter(set == "all", part != "residual") %>%
  mutate(part = factor(part, levels = rev(names(lab_frac))), k = factor(k, levels = k_levels)) %>%
  ggplot(aes(x = 100 * adj_r2, y = part, fill = k)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, colour = "grey20", linewidth = 0.2) +
  scale_y_discrete(labels = lab_frac) +
  scale_fill_manual(values = pal_level, labels = lab_k(k_levels), breaks = k_levels, name = NULL) +
  labs(x = "Lineage composition explained (adjusted R², %)", y = NULL)

p_clim <- (wrap_plots(map(k_levels, plot_triplot), nrow = 1) + plot_layout(guides = "collect")) / p_vp +
  plot_layout(heights = c(1.3, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

save_fig("climate_niche.pdf", p_clim, 7.1, 5.6)

### composition per site group (supplement): share of the units among the MLGs, one panel per level
p_groups <- groups %>%
  mutate(group = factor(group, levels = names(lab_group)), unit = factor(unit, levels = rev(names(pal_unit))),
         k = lab_k(k)) %>%
  ggplot(aes(x = group, y = share, fill = unit)) +
  geom_col(width = 0.7, colour = "grey20", linewidth = 0.2) +
  geom_text(data = ~ distinct(.x, k, group, mlgs, sites),
            aes(x = group, y = 1.04, label = sprintf("%d (%d)", mlgs, sites)),
            size = 6 / .pt, inherit.aes = FALSE) +
  facet_wrap(~ k, nrow = 1) +
  scale_x_discrete(labels = lab_group) +
  scale_fill_manual(values = pal_unit, breaks = names(pal_unit), name = NULL) +
  scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Share of MLGs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.background = element_blank(),
        legend.position = "right")

save_fig("site_group_composition.pdf", p_groups, 7.1, 3)



# tables -----------------------------------------------------------------------

pub_test <- tests %>%
  filter(set == "all") %>%
  left_join(tests %>% filter(set == "no_mixed") %>% select(k, adj_r2_nm = adj_r2, p_nm = p), by = "k") %>%
  transmute(Level = lab_k(k), Sites = sites, MLGs = mlg,
            `Climate variables` = paste(clim_kept, collapse = ", "),
            `Adjusted R² (%)` = sprintf("%.1f", 100 * adj_r2), F = sprintf("%.2f", F), p = fmt_p(p),
            `Largest VIF` = sprintf("%.2f", vif_max),
            `Without mixed MAT: adjusted R² (%), p` = sprintf("%.1f, %s", 100 * adj_r2_nm, fmt_p(p_nm)))

pub_margin <- margins %>%
  filter(set == "all") %>%
  left_join(margins %>% filter(set == "no_mixed") %>% select(k, var, p_nm = p), by = c("k", "var")) %>%
  transmute(Level = lab_k(k), Variable = var, F = sprintf("%.2f", F), p = fmt_p(p), VIF = sprintf("%.2f", vif),
            `Without mixed MAT: p` = fmt_p(p_nm))

lab_parts <- c(climate_only = "Climate only", geography_only = "Geography only", project_only = "Project only",
               host_only = "Host only", shared = "Shared", residual = "Residual")
pub_varpart <- vparts_sum %>%
  filter(set == "all") %>%
  left_join(vparts_sum %>% filter(set == "no_mixed") %>% select(k, part, adj_r2_nm = adj_r2, p_nm = p),
            by = c("k", "part")) %>%
  mutate(part = factor(part, levels = names(lab_parts))) %>%
  arrange(k, part) %>%
  transmute(Level = lab_k(k), Fraction = lab_parts[as.character(part)],
            `Adjusted R² (%)` = sprintf("%.1f", 100 * adj_r2), p = if_else(is.na(p), "", fmt_p(p)),
            `Without mixed MAT: adjusted R² (%), p` = sprintf("%.1f%s", 100 * adj_r2_nm,
                                                             if_else(is.na(p_nm), "", paste0(", ", fmt_p(p_nm)))))

pub_host <- host_magic %>%
  filter(set == "all") %>%
  left_join(host_magic %>% filter(set == "no_mixed") %>% select(k, adj_r2_nm = adj_r2, p_nm = p), by = "k") %>%
  transmute(Level = lab_k(k), `MAGIC sites (oak)` = sprintf("%d (%d)", sites, sites_oak),
            `Adjusted R² (%)` = sprintf("%.1f", 100 * adj_r2), F = sprintf("%.2f", F), p = fmt_p(p),
            `Without mixed MAT: adjusted R² (%), p` = sprintf("%.1f, %s", 100 * adj_r2_nm, fmt_p(p_nm)))

pub_groups <- groups %>%
  mutate(cell = sprintf("%.2f", share)) %>%
  select(k, group, sites, mlgs, unit, cell) %>%
  pivot_wider(names_from = unit, values_from = cell, values_fill = "") %>%
  mutate(group = factor(group, levels = names(lab_group))) %>%
  arrange(k, group) %>%
  transmute(Level = lab_k(k), `Site group` = lab_group[as.character(group)], Sites = sites, MLGs = mlgs,
            !!!syms(intersect(names(pal_unit), unique(groups$unit))))

pub_four <- four_sum %>%
  filter(set == "all") %>%
  left_join(four_sum %>% filter(set == "no_mixed") %>% select(k, part, adj_r2_nm = adj_r2), by = c("k", "part")) %>%
  transmute(Level = lab_k(k), `Site and depth combinations` = cells, Sites = sites,
            Fraction = coalesce(lab_frac[part], "Residual"),
            `Adjusted R² (%)` = sprintf("%.1f", 100 * adj_r2), p = if_else(is.na(p), "", fmt_p(p)),
            `Without mixed MAT: adjusted R² (%)` = sprintf("%.1f", 100 * adj_r2_nm))

write_tsv(tests, file.path(out_root, "climate_niche_tests.tsv"))
write_tsv(margins, file.path(out_root, "climate_niche_margins.tsv"))
write_tsv(vparts, file.path(out_root, "climate_niche_varpart_fractions.tsv"))
write_tsv(vparts_sum, file.path(out_root, "climate_niche_varpart.tsv"))
write_tsv(host_magic, file.path(out_root, "host_within_magic.tsv"))
write_tsv(groups, file.path(out_root, "site_group_composition.tsv"))
write_tsv(four, file.path(out_root, "soil_climate_varpart_fractions.tsv"))
write_tsv(four_sum, file.path(out_root, "soil_climate_varpart.tsv"))
write_tsv(scores_tbl, file.path(out_root, "climate_niche_scores.tsv"))
write_tsv(pub_test, file.path(out_root, "publication_climate_niche.tsv"))
write_tsv(pub_margin, file.path(out_root, "publication_climate_niche_margins.tsv"))
write_tsv(pub_varpart, file.path(out_root, "publication_climate_niche_varpart.tsv"))
write_tsv(pub_host, file.path(out_root, "publication_host_within_magic.tsv"))
write_tsv(pub_groups, file.path(out_root, "publication_site_group_composition.tsv"))
write_tsv(pub_four, file.path(out_root, "publication_soil_climate_varpart.tsv"))



# summary ----------------------------------------------------------------------

cat("\n=== Climate across all sites ===\n")
print(pub_test, width = Inf)

cat("\n=== Marginal tests ===\n")
print(pub_margin, n = Inf, width = Inf)

cat("\n=== Climate, geography, project and host ===\n")
print(pub_varpart, n = Inf, width = Inf)

cat("\n=== Host within the MAGIC sites ===\n")
print(pub_host, width = Inf)

cat("\n=== Composition per site group ===\n")
print(pub_groups, n = Inf, width = Inf)

cat("\n=== Soil, depth, climate and site context at the sites with soil data ===\n")
print(pub_four, n = Inf, width = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  climate_niche.pdf\n")
cat("  site_group_composition.pdf\n")
cat("  climate_niche_tests.tsv\n")
cat("  climate_niche_margins.tsv\n")
cat("  climate_niche_varpart_fractions.tsv\n")
cat("  climate_niche_varpart.tsv\n")
cat("  host_within_magic.tsv\n")
cat("  site_group_composition.tsv\n")
cat("  soil_climate_varpart_fractions.tsv\n")
cat("  soil_climate_varpart.tsv\n")
cat("  climate_niche_scores.tsv\n")
cat("  publication_climate_niche.tsv\n")
cat("  publication_climate_niche_margins.tsv\n")
cat("  publication_climate_niche_varpart.tsv\n")
cat("  publication_host_within_magic.tsv\n")
cat("  publication_site_group_composition.tsv\n")
cat("  publication_soil_climate_varpart.tsv\n")
