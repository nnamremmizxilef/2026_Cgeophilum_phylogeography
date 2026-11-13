## ============================================================================
## 040_soil_niche.R
## Soil niche of the units at K = 4, 5 and 6 at the sites with in situ soil data:
## partial RDA across sites and per project with variance partitioning, depth and
## soil within the PhytOakmeter sites, and depth profiles per unit
##
## Files needed:
##   soildata_ALL_FIXED.txt         external: in situ soil data per PhytOakmeter site and depth and per
##                                  MAGIC site (copy manually into local/04_diversity_niche_selection/supplementary)
##   lineage_labels.tsv             unit at K = 4, 5 and 6 per ingroup sample,
##                                  from local/03_delimitation/001_lineage_labels.R
##   sample_metadata_ALL_FIXED.tsv  per-sample table (project, site, dominant_host, sampling_depth, mlg,
##                                  mat_genome), from local/01_variant_calling/103_mating_type.R
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)
library(ggrepel)
library(vegan)
library(lme4)



# files ------------------------------------------------------------------------

soil_file   <- "04_diversity_niche_selection/supplementary/soildata_ALL_FIXED.txt"
labels_file <- "03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv"
meta_file   <- "01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv"
out_root    <- "04_diversity_niche_selection/results/040_results_soil_niche"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### levels and sample sets
k_levels <- c("k4", "k5", "k6")
sets     <- c("all", "no_mixed")

### projects (all sites), projects sampled at three depths (within-site analyses and depth profiles), depth
### code of the single MAGIC sample (0-10 cm, empty in the metadata)
projects       <- c("magic", "phytoakmeter")
depth_projects <- "phytoakmeter"
depth_single   <- c(magic = "1")

### mating type calls of possible mixed cultures (sample set without them)
mat_exclude <- c("MAT12", "ambiguous")

### candidate soil variables: column, model name, label, ln (right-skewed concentrations), group of correlated
### variables, kept (one per group), persistent (water content and mineral N: snapshots of the sampling day)
soil_vars <- tribble(
  ~column,         ~var,      ~label,    ~log,  ~group,                   ~main, ~persistent,
  "pH_H20",        "pH",      "pH",      FALSE, "Acidity",                TRUE,  TRUE,
  "CaO",           "CaO",     "CaO",     TRUE,  "Acidity",                FALSE, TRUE,
  "C",             "C",       "C",       TRUE,  "Organic matter",         TRUE,  TRUE,
  "N_total",       "N",       "N",       TRUE,  "Organic matter",         FALSE, TRUE,
  "C_N",           "CN",      "C:N",     FALSE, "Organic matter quality", TRUE,  TRUE,
  "Prop_Clay",     "Clay",    "Clay",    FALSE, "Texture and exchange",   TRUE,  TRUE,
  "Cations",       "Cations", "Cations", TRUE,  "Texture and exchange",   FALSE, TRUE,
  "sand",          "Sand",    "Sand",    FALSE, "Texture and exchange",   FALSE, TRUE,
  "K2O",           "K2O",     "K2O",     TRUE,  "Texture and exchange",   FALSE, TRUE,
  "MgO",           "MgO",     "MgO",     TRUE,  "Texture and exchange",   FALSE, TRUE,
  "Na2O",          "Na2O",    "Na2O",    TRUE,  "Texture and exchange",   FALSE, TRUE,
  "P2O5",          "P2O5",    "P2O5",    TRUE,  "Phosphorus",             TRUE,  TRUE,
  "Water_Content", "Water",   "Water",   FALSE, "Sampling day",           FALSE, FALSE,
  "NO3",           "NO3",     "NO3",     TRUE,  "Sampling day",           FALSE, FALSE,
  "NH4",           "NH4",     "NH4",     TRUE,  "Sampling day",           FALSE, FALSE)

### soil table columns not used: pH_KCl (second pH), MOT (C x 1.724), Lime_total and Lime_active (zero
### everywhere), silt (clay + silt + sand = 100), K2IAB (humus mineralisation coefficient computed from clay), R

### kept variables and the sensitivity set without P2O5 (phosphorus extracted differently in the projects)
main_vars <- soil_vars$var[soil_vars$main]
nop_vars  <- setdiff(main_vars, "P2O5")

### depth codes and labels (PhytOakmeter)
lab_depth <- c(`1` = "5-10 cm", `2` = "25-30 cm", `3` = "45-50 cm")

### smallest number of MLGs of a unit for its depth profile; draws for the interval of the share
n_min_mlg <- 15
n_draws   <- 1000

### permutations, significance, correlation threshold of the groups, histogram bins of the diagnostics, seed
n_perm  <- 999
alpha   <- 0.05
r_max   <- 0.7
n_bins  <- 15
seed    <- 42



# style ------------------------------------------------------------------------

### unit colours
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

### project fills and labels, level fills (variance partitioning)
pal_project <- c(magic = "white", phytoakmeter = "grey25")
lab_project <- c(magic = "MAGIC", phytoakmeter = "PhytOakmeter")
pal_level   <- c(k4 = "grey80", k5 = "grey55", k6 = "grey25")

### sampled layers (project and depth code) and hosts: labels and shapes
lab_layer   <- c(phytoakmeter_1 = "5-10 cm", phytoakmeter_2 = "25-30 cm", phytoakmeter_3 = "45-50 cm",
                 magic_1 = "0-10 cm (MAGIC)")
shape_layer <- c(phytoakmeter_1 = 16, phytoakmeter_2 = 1, phytoakmeter_3 = 4, magic_1 = 17)
lab_host    <- c(Quercus_robur = "Oak", Fagus_sylvatica = "Beech")
shape_host  <- c(Quercus_robur = 17, Fagus_sylvatica = 2)

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}

### p values in tables
fmt_p   <- function(p) if_else(p < 0.001, "< 0.001", sprintf("%.3f", p))

### model fit for the triplot corner (adjusted R2 and p of the overall test), level and variable labels
lab_fit <- function(r2, p) sprintf("Adj. R² = %.1f%%, %s", 100 * r2,
                                   if_else(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))
lab_k   <- function(k) sub("k", "K = ", k)
lab_var <- function(v) paste(soil_vars$label[match(v, soil_vars$var)], collapse = ", ")

### moment skewness
skew <- function(x) {
  x <- x[!is.na(x)]
  mean((x - mean(x))^3) / mean((x - mean(x))^2)^1.5
}



# read -------------------------------------------------------------------------

lab  <- read_tsv(labels_file, col_types = cols(.default = "c"))
meta <- read_tsv(meta_file, col_types = cols(.default = "c"))

### project and host of every site (the sampled tree species: first species listed in dominant_host)
site_info <- meta %>%
  filter(project %in% projects, !is.na(site)) %>%
  mutate(host = str_trim(str_remove(dominant_host, ";.*"))) %>%
  distinct(site, project, host)

### soil per site and depth (site codes with "-" as in the metadata), with project and host of the site
soil <- read_tsv(soil_file, col_types = cols(Site = "c", .default = "d")) %>%
  mutate(site  = str_replace_all(Site, "_", "-"),
         depth = as.character(Depth),
         sand  = Prop_Sand_fine + Prop_Sand_coarse) %>%
  left_join(site_info, by = "site")

### isolates with project, site, depth and MLG (PhytOakmeter: four directions pooled per site and depth, as the
### one soil sample per depth; MAGIC: one sample per site; clonemates counted once per site and depth)
iso <- lab %>%
  select(sample, all_of(k_levels)) %>%
  inner_join(meta %>% filter(project %in% projects) %>%
               select(sample, project, site, depth = sampling_depth, mlg, mat_genome), by = "sample") %>%
  mutate(depth = coalesce(depth, depth_single[project]), mixed = mat_genome %in% mat_exclude)



# checks -----------------------------------------------------------------------

## check if every candidate column is in the soil table
need <- c(soil_vars$column, "Lat", "Lon", "Depth")

if (!all(need %in% names(soil))) {
  stop("soil columns missing: ", paste(setdiff(need, names(soil)), collapse = ", "))
}

## check if every site has one project and one host
dup <- site_info$site[duplicated(site_info$site)]

if (length(dup) > 0) {
  stop("sites with more than one project or host: ", paste(unique(dup), collapse = ", "))
}

## check if every soil site has a project and a host in the metadata
if (anyNA(soil$project) || anyNA(soil$host)) {
  stop("soil sites without project or host in the metadata: ",
       paste(soil$site[is.na(soil$project) | is.na(soil$host)], collapse = ", "))
}

## check if every host has a label
if (!all(soil$host %in% names(lab_host))) {
  stop("hosts without a label: ", paste(setdiff(soil$host, names(lab_host)), collapse = ", "))
}

## check if every site and depth is in the soil table once
if (anyDuplicated(soil[c("site", "depth")])) {
  stop("site and depth combinations repeated in the soil table")
}

## check if the kept variables are complete at all sites
if (anyNA(soil[soil_vars$column[soil_vars$main]])) {
  stop("missing values in the kept soil variables")
}

## check if all candidates are complete at the sites sampled at three depths
if (anyNA(filter(soil, project %in% depth_projects)[soil_vars$column])) {
  stop("missing soil values at the ", paste(depth_projects, collapse = ", "), " sites")
}

### candidates not measured at every site
not_all <- soil_vars$var[map_lgl(soil_vars$column, ~ anyNA(soil[[.x]]))]

if (length(not_all) > 0) {
  message("candidates not measured at every site (left out of the all-site analyses): ", lab_var(not_all))
}

## check if the ln-transformed variables are positive
neg <- soil_vars$column[soil_vars$log & map_lgl(soil_vars$column, ~ any(soil[[.x]] <= 0, na.rm = TRUE))]

if (length(neg) > 0) {
  stop("zero or negative values in ln-transformed variables: ", paste(neg, collapse = ", "))
}

## check if the texture fractions add up to 100 (where measured)
tex <- soil$Prop_Clay + soil$Prop_Silt_fine + soil$Prop_Silt_coarse + soil$sand

if (any(abs(tex - 100) > 1, na.rm = TRUE)) {
  warning("texture fractions not adding up to 100 at ", sum(abs(tex - 100) > 1, na.rm = TRUE),
          " site and depth combinations")
}

## check if every isolate has a depth and an MLG
if (anyNA(iso$depth) || anyNA(iso$mlg)) {
  stop("isolates without a sampling depth or MLG")
}

## check if every isolate site and depth has soil data
no_soil <- anti_join(distinct(iso, site, depth), soil, by = c("site", "depth"))

if (nrow(no_soil) > 0) {
  warning(nrow(no_soil), " site and depth combinations with isolates but without soil data are left out: ",
          paste(no_soil$site, no_soil$depth, collapse = ", "))
}

### site and depth combinations with soil data but without isolates
no_iso <- anti_join(soil, distinct(iso, site, depth), by = c("site", "depth"))

if (nrow(no_iso) > 0) {
  message(nrow(no_iso), " site and depth combinations without isolates: ",
          paste(no_iso$site, no_iso$depth, collapse = ", "))
}

### sites per project and host
print(count(distinct(soil, project, host, site), project, host))



# data -------------------------------------------------------------------------

### soil variables in the models (ln where flagged) and site coordinates
soil_x <- soil %>%
  transmute(project, host, site, depth, lat = Lat, lon = Lon,
            !!!setNames(map2(soil_vars$column, soil_vars$log,
                             function(cl, lg) if (lg) log(soil[[cl]]) else soil[[cl]]),
                        soil_vars$var))

### candidates measured at every site; the persistent ones for the forward selection (sensitivity)
vars_all_sites <- soil_vars$var[map_lgl(soil_vars$column, ~ !anyNA(soil[[.x]]))]
fwd_vars       <- intersect(soil_vars$var[soil_vars$persistent], vars_all_sites)

### standardisation within project (removes differences in level and scale between the soil analyses of the
### projects, e.g. pH in water or CaCl2, phosphorus extraction)
z_within <- function(x, g) ave(x, g, FUN = function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE))

### one column per variable and project, zero outside the project (soil responses free to differ between
### projects; columns named variable@project)
by_project <- function(X, g) {
  out <- do.call(cbind, map(sort(unique(g)), function(p) {
    x <- X
    x[g != p, ] <- 0
    colnames(x) <- paste(colnames(X), p, sep = "@")
    x
  }))
  out[, order(match(sub("@.*", "", colnames(out)), colnames(X))), drop = FALSE]
}

### columns of variables (all projects) in a matrix of by_project
cols_of <- function(X, v) colnames(X)[sub("@.*", "", colnames(X)) %in% v]

### dummy columns of a factor (none when it has one level)
dummies <- function(x) {
  if (n_distinct(x) > 1) model.matrix(~ x)[, -1, drop = FALSE] else matrix(0, length(x), 0)
}

### site coordinates and great-circle distances between the sites (km)
site_xy <- distinct(soil_x, project, host, site, lat, lon) %>% arrange(site)
gc_dist <- function(lat, lon) {
  r <- pi / 180
  a <- outer(lat, lat, function(x, y) sin((y - x) * r / 2)^2) +
    outer(cos(lat * r), cos(lat * r)) * outer(lon, lon, function(x, y) sin((y - x) * r / 2)^2)
  6371 * 2 * asin(pmin(sqrt(a), 1))
}

### distance-based Moran's eigenvector maps of the sites (as adespatial::dbmem: distances truncated at the
### longest edge of the minimum spanning tree, weights 1 - (d / 4t)^2, eigenvectors of the doubly centred
### weights with positive eigenvalues)
dbmem <- function(d) {
  t <- max(spantree(as.dist(d))$dist)
  w <- ifelse(d <= t, 1 - (d / (4 * t))^2, 0)
  diag(w) <- 0
  h <- diag(nrow(w)) - 1 / nrow(w)
  e <- eigen(h %*% w %*% h, symmetric = TRUE)
  m <- e$vectors[, e$values > 1e-8 * max(e$values), drop = FALSE]
  dimnames(m) <- list(rownames(d), paste0("MEM", seq_len(ncol(m))))
  m
}
site_d <- gc_dist(site_xy$lat, site_xy$lon)
dimnames(site_d) <- list(site_xy$site, site_xy$site)
mem_sites <- dbmem(site_d)

### MLGs per unit per site and depth, with soil, for a level, sample set and projects
comp <- function(kk, ss, prj = projects) {
  d <- if (ss == "all") iso else filter(iso, !mixed)
  d <- filter(d, project %in% prj)
  tibble(site = d$site, depth = d$depth, mlg = d$mlg, unit = d[[kk]]) %>%
    distinct() %>%
    count(site, depth, unit) %>%
    pivot_wider(names_from = unit, values_from = n, values_fill = 0) %>%
    inner_join(soil_x, by = c("site", "depth")) %>%
    arrange(project, site, depth)
}

### response (Hellinger) and explanatory matrices: soil standardised within project (Xw) and per project (Xp),
### geography (standardised coordinates), project, host, dbMEMs, depth and site (dummies); matrix interface of
### rda (models independent of formula environments)
model_data <- function(kk, ss, prj = projects) {
  dat   <- comp(kk, ss, prj)
  units <- intersect(names(pal_unit), names(dat))
  Xw    <- sapply(soil_vars$var, function(v) z_within(dat[[v]], dat$project))
  depth <- factor(dat$depth, levels = names(lab_depth))
  list(dat = dat, units = units, depth = depth, site = factor(dat$site),
       layer = paste(dat$project, dat$depth, sep = "_"),
       R  = decostand(as.matrix(dat[units]), "hellinger"),
       Xw = Xw,
       Xp = by_project(Xw, dat$project),
       G  = scale(as.matrix(dat[c("lat", "lon")])),
       P  = dummies(dat$project),
       H  = dummies(dat$host),
       M  = mem_sites[dat$site, , drop = FALSE],
       D  = model.matrix(~ depth)[, -1, drop = FALSE],
       S  = dummies(dat$site))
}

### test of a partial RDA (response R, explanatory X, conditions Z)
test_rda <- function(R, X, Z, perm) {
  m <- rda(R, X, Z)
  a <- anova(m, permutations = perm)
  tibble(adj_r2 = RsquareAdj(m)$adj.r.squared, F = a$F[1], p = a$`Pr(>F)`[1])
}

### forward selection over blocks of columns (a soil variable in every project, or one dbMEM) after a
### significant global test, stopping at alpha or at the adjusted R2 of the global model
forward <- function(R, X, blocks, Z, perm) {
  cols <- function(b) unlist(blocks[b], use.names = FALSE)
  glob <- test_rda(R, X[, cols(names(blocks)), drop = FALSE], Z, perm)
  out  <- list(global = glob, sel = character(0))
  if (glob$p >= alpha) return(out)
  r2_cur <- -Inf
  repeat {
    cand <- setdiff(names(blocks), out$sel)
    if (length(cand) == 0) break
    r2 <- map_dbl(cand, ~ RsquareAdj(rda(R, X[, cols(c(out$sel, .x)), drop = FALSE], Z))$adj.r.squared)
    best <- cand[which.max(r2)]
    if (max(r2) <= r2_cur || max(r2) > glob$adj_r2) break
    t <- test_rda(R, X[, cols(best), drop = FALSE], cbind(Z, X[, cols(out$sel), drop = FALSE]), perm)
    if (t$p >= alpha) break
    out$sel <- c(out$sel, best)
    r2_cur  <- max(r2)
  }
  out
}



# transformations and variable groups ------------------------------------------

### skewness of every candidate before and after the ln-transformation (all site and depth combinations)
transforms <- soil_vars %>%
  transmute(var, label, log, skew_raw = map_dbl(column, ~ skew(soil[[.x]])),
            skew_ln = if_else(log, map_dbl(column, ~ skew(log(soil[[.x]]))), NA_real_))

### correlations among the candidates within projects (standardised within project, pairs measured at both)
### and share of the variation of each variable between sites
soil_w   <- sapply(soil_vars$var, function(v) z_within(soil_x[[v]], soil_x$project))
soil_cor <- cor(soil_w, use = "pairwise.complete.obs")
kept_cor <- soil_cor[main_vars, main_vars]

between <- map_dfr(soil_vars$var, function(v) {
  a <- anova(lm(soil_x[[v]] ~ factor(soil_x$site)))
  tibble(var = v, between_sites = a$`Sum Sq`[1] / sum(a$`Sum Sq`))
})



# checks on the kept variables -------------------------------------------------

## check if the kept variables are below the correlation threshold
r_kept <- max(abs(kept_cor[upper.tri(kept_cor)]))

if (r_kept > r_max) {
  warning("kept soil variables correlated above ", r_max, " (largest |r| = ", round(r_kept, 2), ")")
}

## check if every group of correlated variables holds one kept variable
grp_kept <- soil_vars %>% filter(persistent) %>% summarise(n = sum(main), .by = group)

if (any(grp_kept$n != 1)) {
  warning("groups without exactly one kept variable: ", paste(grp_kept$group[grp_kept$n != 1], collapse = ", "))
}



# across sites -----------------------------------------------------------------

### dbMEMs per level and sample set: forward selection on the composition per site (MLGs summed over the
### depths), conditioned on latitude, longitude, project and host (spatial structure beyond the linear trend)
select_mem <- function(kk, ss) {
  md <- model_data(kk, ss)
  cs <- rowsum(as.matrix(md$dat[md$units]), md$dat$site)
  sx <- site_xy[match(rownames(cs), site_xy$site), ]
  M  <- mem_sites[rownames(cs), , drop = FALSE]
  f  <- forward(decostand(cs, "hellinger"), M, setNames(as.list(colnames(M)), colnames(M)),
                cbind(scale(as.matrix(sx[c("lat", "lon")])), dummies(sx$project), dummies(sx$host)),
                how(nperm = n_perm))
  tibble(k = kk, set = ss, mem_adj_r2 = f$global$adj_r2, mem_p = f$global$p, mem_sel = list(f$sel))
}

set.seed(seed)
mems <- map2_dfr(rep(k_levels, each = length(sets)), rep(sets, length(k_levels)), select_mem)
mem_of <- function(kk, ss) mems$mem_sel[[which(mems$k == kk & mems$set == ss)]]

### kept soil variables with responses per project (Xs) given the site context (Z: geography, project, host):
### overall (free and within sites) and marginal tests, VIF; common response (Xc) and the project-specific
### part beyond it; variance partitioning into soil, depth and site context; without P2O5; with the dbMEMs
fit_across <- function(kk, ss) {
  md   <- model_data(kk, ss)
  Xs   <- md$Xp[, cols_of(md$Xp, main_vars), drop = FALSE]
  Xc   <- md$Xw[, main_vars, drop = FALSE]
  Z    <- cbind(md$G, md$P, md$H)
  free <- how(nperm = n_perm)

  m      <- rda(md$R, Xs, Z)
  a_free <- test_rda(md$R, Xs, Z, free)
  a_site <- test_rda(md$R, Xs, Z, how(nperm = n_perm, blocks = md$site))
  a_com  <- test_rda(md$R, Xc, Z, free)
  ### project-specific part: columns of all but one project given the common response (same space as all
  ### project columns, without aliased constraints)
  Xd     <- Xs[, sub(".*@", "", colnames(Xs)) != last(sort(unique(md$dat$project))), drop = FALSE]
  a_spec <- test_rda(md$R, Xd, cbind(Z, Xc), free)
  a_nop  <- test_rda(md$R, md$Xp[, cols_of(md$Xp, nop_vars), drop = FALSE], Z, free)
  mem    <- mem_of(kk, ss)
  a_mem  <- test_rda(md$R, Xs, cbind(Z, md$M[, mem, drop = FALSE]), free)
  a_marg <- map_dfr(main_vars, ~ test_rda(md$R, Xs[, cols_of(Xs, .x), drop = FALSE],
                                          cbind(Z, Xs[, cols_of(Xs, setdiff(main_vars, .x)), drop = FALSE]), free))
  vif    <- vif.cca(m)

  vp     <- varpart(md$R, Xs, md$D, Z)
  t_soil <- test_rda(md$R, Xs, cbind(md$D, Z), free)
  t_dep  <- test_rda(md$R, md$D, cbind(Xs, Z), free)
  t_ctx  <- test_rda(md$R, Z, cbind(Xs, md$D), free)

  list(
    test = tibble(k = kk, set = ss, cells = nrow(md$dat), sites = n_distinct(md$dat$site),
                  sites_magic = n_distinct(md$dat$site[md$dat$project == "magic"]), mlg = sum(md$dat[md$units]),
                  adj_r2 = a_free$adj_r2, F = a_free$F, p = a_free$p, p_within_sites = a_site$p,
                  adj_r2_common = a_com$adj_r2, p_common = a_com$p,
                  adj_r2_specific = a_spec$adj_r2, p_specific = a_spec$p,
                  adj_r2_nop = a_nop$adj_r2, p_nop = a_nop$p, n_mem = length(mem),
                  mems = paste(mem, collapse = ", "), adj_r2_mem = a_mem$adj_r2, p_mem = a_mem$p,
                  vif_max = max(vif[colnames(Xs)])),
    margin = tibble(k = kk, set = ss, var = main_vars, F = a_marg$F, p = a_marg$p,
                    vif = map_dbl(main_vars, ~ max(vif[cols_of(Xs, .x)]))),
    varpart = tibble(k = kk, set = ss,
                     fraction = c("soil_only", "depth_only", "context_only", "soil_depth", "depth_context",
                                  "soil_context", "all_three", "residual"),
                     adj_r2 = vp$part$indfract$Adj.R.square,
                     p = c(t_soil$p, t_dep$p, t_ctx$p, rep(NA, 5))))
}

set.seed(seed)
fits <- expand_grid(k = k_levels, set = sets) %>% mutate(fit = map2(k, set, fit_across))

tests   <- map_dfr(fits$fit, "test")
margins <- map_dfr(fits$fit, "margin")
vparts  <- map_dfr(fits$fit, "varpart")

### per project: kept soil variables given geography and host (where it varies), free permutations and within
### sites (where sites hold several depths); the models give the triplots
fit_project <- function(kk, ss, pp) {
  md <- model_data(kk, ss, pp)
  X  <- md$Xw[, main_vars, drop = FALSE]
  Z  <- cbind(md$G, md$H)
  a  <- test_rda(md$R, X, Z, how(nperm = n_perm))
  p_site <- if (any(table(md$site) > 1)) test_rda(md$R, X, Z, how(nperm = n_perm, blocks = md$site))$p else NA_real_
  list(test = tibble(k = kk, set = ss, project = pp, cells = nrow(md$dat), sites = n_distinct(md$dat$site),
                     mlg = sum(md$dat[md$units]), adj_r2 = a$adj_r2, F = a$F, p = a$p, p_within_sites = p_site),
       model = rda(md$R, X, Z), dat = md$dat, layer = md$layer)
}

set.seed(seed)
fits_proj <- expand_grid(k = k_levels, set = sets, project = projects) %>%
  mutate(fit = pmap(list(k, set, project), fit_project))
tests_proj <- map_dfr(fits_proj$fit, "test")

### sensitivity: forward selection over the persistent candidates measured at all sites, responses per project
### (all isolates)
set.seed(seed)
selection <- tibble(k = k_levels) %>%
  mutate(res = map(k, function(kk) {
    md <- model_data(kk, "all")
    X  <- md$Xp[, cols_of(md$Xp, fwd_vars), drop = FALSE]
    bl <- setNames(map(fwd_vars, ~ cols_of(X, .x)), fwd_vars)
    Z  <- cbind(md$G, md$P, md$H)
    f  <- forward(md$R, X, bl, Z, how(nperm = n_perm))
    f$vif <- if (length(f$sel) > 0) vif.cca(rda(md$R, X[, cols_of(X, f$sel), drop = FALSE], Z)) else numeric(0)
    f
  }))



# within sites -----------------------------------------------------------------

### depth and the kept soil variables per level and sample set at the sites sampled at three depths,
### conditioned on site, permutations within sites
fit_within <- function(kk, ss) {
  md   <- model_data(kk, ss, depth_projects)
  perm <- how(nperm = n_perm, blocks = md$site)
  t_depth <- test_rda(md$R, md$D, md$S, perm)
  t_soil  <- test_rda(md$R, md$Xw[, main_vars, drop = FALSE], md$S, perm)
  tibble(k = kk, set = ss, sites = n_distinct(md$dat$site),
         depth_adj_r2 = t_depth$adj_r2, depth_F = t_depth$F, depth_p = t_depth$p,
         soil_adj_r2 = t_soil$adj_r2, soil_F = t_soil$F, soil_p = t_soil$p)
}

set.seed(seed)
within_res <- map2_dfr(rep(k_levels, each = length(sets)), rep(sets, length(k_levels)), fit_within)



# depth profiles per unit ------------------------------------------------------

### share of a unit among the MLGs of a site and depth at the sites sampled at three depths: binomial mixed
### model (depth fixed, site random), likelihood-ratio test against the model without depth
profile_unit <- function(kk, ss, uu) {
  d <- comp(kk, ss, depth_projects)
  units <- intersect(names(pal_unit), names(d))
  d <- d %>% mutate(y = .data[[uu]], n = rowSums(across(all_of(units))),
                    depth = factor(depth, levels = names(lab_depth)))
  ctrl <- glmerControl(optimizer = "bobyqa")
  m1 <- glmer(cbind(y, n - y) ~ depth + (1 | site), family = binomial, data = d, control = ctrl)
  m0 <- glmer(cbind(y, n - y) ~ 1 + (1 | site), family = binomial, data = d, control = ctrl)
  lrt <- anova(m0, m1)

  ### share per depth averaged over sites (random effect integrated out), 95% interval from draws of the
  ### fixed effects
  X     <- model.matrix(~ depth, data = tibble(depth = factor(names(lab_depth), levels = names(lab_depth))))
  z     <- qnorm(ppoints(200)) * sqrt(VarCorr(m1)$site[1])
  draws <- MASS::mvrnorm(n_draws, fixef(m1), as.matrix(vcov(m1)))
  marg  <- function(b) rowMeans(plogis(outer(drop(X %*% b), z, "+")))
  sim   <- apply(draws, 1, marg)
  tibble(k = kk, set = ss, unit = uu, mlg = sum(d$y), chisq = lrt$Chisq[2], df = lrt$Df[2],
         p = lrt$`Pr(>Chisq)`[2], depth = names(lab_depth), share = marg(fixef(m1)),
         lo = apply(sim, 1, quantile, 0.025), hi = apply(sim, 1, quantile, 0.975),
         observed = as.numeric(tapply(d$y, d$depth, sum) / tapply(d$n, d$depth, sum)))
}

### units with at least n_min_mlg MLGs per level and sample set
set.seed(seed)
profiles <- expand_grid(k = k_levels, set = sets) %>%
  mutate(units = map2(k, set, function(kk, ss) {
    d <- comp(kk, ss, depth_projects)
    u <- intersect(names(pal_unit), names(d))
    u[colSums(d[u]) >= n_min_mlg]
  })) %>%
  unnest(units) %>%
  pmap_dfr(function(k, set, units) profile_unit(k, set, units))

### q values (BH within level and sample set)
profiles <- profiles %>%
  left_join(profiles %>% distinct(k, set, unit, p) %>% mutate(q = p.adjust(p, "BH"), .by = c(k, set)) %>%
              select(-p),
            by = c("k", "set", "unit"))



# ordination scores ------------------------------------------------------------

### triplot scores per level and project (all isolates): site and depth combinations, units, soil variables
xy <- function(x, id) {
  x <- as.matrix(x)
  if (ncol(x) == 1) x <- cbind(x, 0)
  colnames(x) <- c("A1", "A2")
  as_tibble(x, rownames = id)
}

scores_tbl <- fits_proj %>%
  filter(set == "all") %>%
  pmap_dfr(function(k, set, project, fit) {
    m   <- fit$model
    eig <- c(m$CCA$eig, m$CA$eig)
    ax  <- sprintf("%s (%.1f%%)", names(eig)[1:2], 100 * eig[1:2] / m$tot.chi)
    st  <- xy(scores(m, display = "sites", choices = 1:2, scaling = 2), "id") %>%
      mutate(type = "cell", depth = fit$dat$depth, host = fit$dat$host, layer = fit$layer)
    sp  <- xy(scores(m, display = "species", choices = 1:2, scaling = 2), "id") %>% mutate(type = "unit")
    bp  <- xy(scores(m, display = "bp", choices = 1:2, scaling = 2), "id") %>% mutate(type = "variable")
    bp  <- mutate(bp, across(c(A1, A2), ~ .x * 0.9 * max(abs(c(st$A1, st$A2))) / max(abs(c(bp$A1, bp$A2)))))
    bind_rows(st, sp, bp) %>% mutate(k = k, model = project, axis1 = ax[1], axis2 = ax[2],
                                     fit = lab_fit(fit$test$adj_r2, fit$test$p), .before = 1)
  }) %>%
  mutate(label = if_else(type == "variable", soil_vars$label[match(id, soil_vars$var)], id))



# figures ----------------------------------------------------------------------

### transformation diagnostics (supplement): every candidate on the raw scale and, where transformed, on the
### ln scale, with skewness; projects stacked
diag_long <- soil_vars %>%
  select(var, label, column, log) %>%
  mutate(data = map2(column, log, function(cl, lg) {
    raw <- tibble(project = soil$project, scale = "raw", value = soil[[cl]])
    if (lg) bind_rows(raw, tibble(project = soil$project, scale = "ln", value = log(soil[[cl]]))) else raw
  })) %>%
  unnest(data) %>%
  filter(!is.na(value)) %>%
  left_join(select(transforms, var, skew_raw, skew_ln), by = "var") %>%
  mutate(order = match(var, soil_vars$var[order(!soil_vars$log)]),
         panel = sprintf("%s, %s\nskewness %.2f", label, if_else(log, scale, "untransformed"),
                         if_else(scale == "ln", skew_ln, skew_raw))) %>%
  mutate(panel = fct_reorder(panel, order + (scale == "ln") / 2))

p_diag <- ggplot(diag_long, aes(value, fill = project)) +
  geom_histogram(bins = n_bins, colour = "grey20", linewidth = 0.2) +
  facet_wrap(~ panel, scales = "free", ncol = 6) +
  scale_fill_manual(values = pal_project, labels = lab_project, name = NULL) +
  labs(x = NULL, y = "Site and depth combinations") +
  theme(strip.background = element_blank(), strip.text = element_text(size = 6), axis.text = element_text(size = 5))

save_fig("soil_transformations.pdf", p_diag, 7.1, 6.5)

### triplot of a project and level (all isolates): cells by sampled layer (PhytOakmeter) or host (MAGIC), soil
### arrows, unit scores
plot_triplot <- function(kk, pp) {
  d  <- filter(scores_tbl, k == kk, model == pp)
  st <- filter(d, type == "cell")
  sp <- filter(d, type == "unit")
  bp <- filter(d, type == "variable")
  shp <- if (pp == "magic") list(col = "host", values = shape_host, labels = lab_host) else
    list(col = "layer", values = shape_layer, labels = lab_layer[names(lab_layer) != "magic_1"])
  ggplot() +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_point(data = st, aes(A1, A2, shape = .data[[shp$col]]), colour = "grey55", size = 1.1) +
    geom_segment(data = bp, aes(0, 0, xend = A1, yend = A2), arrow = arrow(length = unit(1.5, "mm")),
                 linewidth = 0.35, colour = "grey15") +
    geom_text_repel(data = bp, aes(A1, A2, label = label), size = 6.5 / .pt, colour = "black", seed = 1) +
    geom_point(data = sp, aes(A1, A2, fill = id), shape = 21, size = 2.6, stroke = 0.4, colour = "black") +
    geom_text_repel(data = sp, aes(A1, A2, label = id), size = 7 / .pt, fontface = "bold", seed = 2,
                    box.padding = 0.4) +
    scale_fill_manual(values = pal_unit, guide = "none") +
    scale_shape_manual(values = shp$values, labels = shp$labels, breaks = names(shp$labels), name = NULL) +
    annotate("text", x = -Inf, y = Inf, label = d$fit[1], hjust = -0.05, vjust = 1.4, size = 6.5 / .pt) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
    labs(x = d$axis1[1], y = d$axis2[1]) +
    theme(legend.key.width = unit(2, "mm"))
}

### variance partitioning across sites per level (all isolates): unique fractions and the shared remainder
lab_frac <- c(soil_only = "Soil only", depth_only = "Depth only", context_only = "Site context only",
              shared = "Shared")
vp_plot <- vparts %>%
  filter(set == "all", fraction != "residual") %>%
  mutate(fraction = if_else(fraction %in% names(lab_frac), fraction, "shared")) %>%
  summarise(adj_r2 = sum(adj_r2), .by = c(k, fraction)) %>%
  mutate(fraction = factor(fraction, levels = rev(names(lab_frac))), k = factor(k, levels = k_levels))

p_vp <- ggplot(vp_plot, aes(x = 100 * adj_r2, y = fraction, fill = k)) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, colour = "grey20", linewidth = 0.2) +
  scale_y_discrete(labels = lab_frac) +
  scale_fill_manual(values = pal_level, labels = lab_k(k_levels), breaks = k_levels, name = NULL) +
  labs(x = "Lineage composition explained (adjusted R², %)", y = NULL)

p_soil <- (wrap_plots(map(k_levels, ~ plot_triplot(.x, "phytoakmeter")), nrow = 1) + plot_layout(guides = "collect")) /
  (wrap_plots(map(k_levels, ~ plot_triplot(.x, "magic")), nrow = 1) + plot_layout(guides = "collect")) /
  p_vp +
  plot_layout(heights = c(1.3, 1.3, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

save_fig("soil_niche.pdf", p_soil, 7.1, 8.4)

### depth profiles per level (all isolates): share per depth averaged over sites with 95% interval; solid lines
### for units with q < alpha
plot_profile <- function(kk, y_title) {
  d <- profiles %>%
    filter(k == kk, set == "all") %>%
    mutate(unit = factor(unit, levels = intersect(names(pal_unit), unit)),
           x = match(depth, names(lab_depth)) + 0.24 * (as.integer(unit) - (nlevels(unit) + 1) / 2) / nlevels(unit),
           sig = q < alpha)
  ggplot(d, aes(x = x, y = share, group = unit)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.3, colour = "grey40") +
    geom_line(aes(linetype = sig), colour = "grey20", linewidth = 0.4) +
    geom_point(aes(fill = unit), shape = 21, size = 2, stroke = 0.3, colour = "black") +
    geom_text_repel(data = filter(d, depth == last(names(lab_depth))), aes(label = unit), nudge_x = 0.3,
                    direction = "y", hjust = 0, size = 6.5 / .pt, segment.size = 0.2, segment.colour = "grey60",
                    seed = 1) +
    scale_fill_manual(values = pal_unit, guide = "none") +
    scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dashed"), guide = "none") +
    scale_x_continuous(breaks = seq_along(lab_depth), labels = lab_depth, limits = c(0.7, 3.8)) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(x = sprintf("Depth (K = %s)", sub("k", "", kk)), y = if (y_title) "Share of MLGs" else NULL)
}

p_prof <- wrap_plots(map2(k_levels, seq_along(k_levels) == 1, plot_profile), nrow = 1) +
  plot_annotation(tag_levels = "A")

save_fig("depth_profiles.pdf", p_prof, 7.1, 2.6)

### variable selection (supplement): correlations within projects among the candidates ordered by group, kept
### variables outlined on the diagonal; share of the variation between sites; variables in the model per level
ord <- soil_vars$var

cor_long <- as_tibble(soil_cor, rownames = "v1") %>%
  pivot_longer(-v1, names_to = "v2", values_to = "r") %>%
  mutate(kept = v1 == v2 & v1 %in% main_vars,
         v1 = factor(v1, levels = ord), v2 = factor(v2, levels = rev(ord)))

p_cor <- ggplot(cor_long, aes(v1, v2, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_tile(data = filter(cor_long, kept), fill = NA, colour = "black", linewidth = 0.5, width = 1, height = 1) +
  geom_text(aes(label = if_else(abs(r) > r_max & v1 != v2, sprintf("%.2f", r), "")), size = 5 / .pt) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1), name = "r") +
  scale_x_discrete(labels = setNames(soil_vars$label, soil_vars$var)) +
  scale_y_discrete(labels = setNames(soil_vars$label, soil_vars$var)) +
  coord_equal() +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "right", legend.key.height = unit(4, "mm"))

p_between <- between %>%
  left_join(select(soil_vars, var, label, main, group), by = "var") %>%
  mutate(var = factor(var, levels = rev(ord))) %>%
  ggplot(aes(x = 100 * between_sites, y = var, fill = main)) +
  geom_col(width = 0.7, colour = "grey20", linewidth = 0.2) +
  scale_fill_manual(values = c(`TRUE` = "grey30", `FALSE` = "grey85"), guide = "none") +
  scale_y_discrete(labels = setNames(soil_vars$label, soil_vars$var)) +
  scale_x_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Between sites (%)", y = NULL)

lab_appr <- c(kept = "Kept", forward = "Forward selection")
in_model <- bind_rows(
  expand_grid(k = k_levels, var = main_vars) %>% mutate(approach = "kept"),
  selection %>% mutate(var = map(res, "sel")) %>% select(k, var) %>% unnest(var) %>% mutate(approach = "forward"))

p_sel <- expand_grid(k = k_levels, approach = names(lab_appr), var = ord) %>%
  left_join(mutate(in_model, used = TRUE), by = c("k", "approach", "var")) %>%
  mutate(used = replace_na(used, FALSE), approach = factor(approach, levels = names(lab_appr)),
         var = factor(var, levels = rev(ord)), k = lab_k(k)) %>%
  ggplot(aes(approach, var, fill = used)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  facet_grid(. ~ k) +
  scale_fill_manual(values = c(`TRUE` = "grey25", `FALSE` = "grey93"), guide = "none") +
  scale_x_discrete(labels = lab_appr) +
  scale_y_discrete(labels = setNames(soil_vars$label, soil_vars$var)) +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.background = element_blank())

p_selection <- (p_cor | p_between | p_sel) +
  plot_layout(widths = c(1.3, 0.5, 0.8)) +
  plot_annotation(tag_levels = "A")

save_fig("soil_variable_selection.pdf", p_selection, 7.1, 3.4)



# tables -----------------------------------------------------------------------

pub_test <- tests %>%
  filter(set == "all") %>%
  left_join(tests %>% filter(set == "no_mixed") %>% select(k, adj_r2_nm = adj_r2, p_nm = p), by = "k") %>%
  transmute(Level = lab_k(k), `Site and depth combinations` = cells,
            `Sites (MAGIC)` = sprintf("%d (%d)", sites, sites_magic),
            MLGs = mlg, `Adjusted R² (%)` = sprintf("%.1f", 100 * adj_r2), F = sprintf("%.2f", F), p = fmt_p(p),
            `p (permutation within sites)` = fmt_p(p_within_sites), `Largest VIF` = sprintf("%.2f", vif_max),
            `Common response: adjusted R² (%), p` = sprintf("%.1f, %s", 100 * adj_r2_common, fmt_p(p_common)),
            `Project-specific part: adjusted R² (%), p` = sprintf("%.1f, %s", 100 * adj_r2_specific,
                                                                 fmt_p(p_specific)),
            `Without P2O5: adjusted R² (%), p` = sprintf("%.1f, %s", 100 * adj_r2_nop, fmt_p(p_nop)),
            `With dbMEMs: MEMs, adjusted R² (%), p` = sprintf("%d, %.1f, %s", n_mem, 100 * adj_r2_mem, fmt_p(p_mem)),
            `Without mixed MAT: adjusted R² (%), p` = sprintf("%.1f, %s", 100 * adj_r2_nm, fmt_p(p_nm)))

pub_project <- tests_proj %>%
  filter(set == "all") %>%
  left_join(tests_proj %>% filter(set == "no_mixed") %>% select(k, project, adj_r2_nm = adj_r2, p_nm = p),
            by = c("k", "project")) %>%
  transmute(Level = lab_k(k), Project = lab_project[project], `Site and depth combinations` = cells, Sites = sites,
            MLGs = mlg, `Adjusted R² (%)` = sprintf("%.1f", 100 * adj_r2), F = sprintf("%.2f", F), p = fmt_p(p),
            `p (permutation within sites)` = if_else(is.na(p_within_sites), "-", fmt_p(p_within_sites)),
            `Without mixed MAT: adjusted R² (%), p` = sprintf("%.1f, %s", 100 * adj_r2_nm, fmt_p(p_nm)))

pub_mem <- mems %>%
  transmute(Level = lab_k(k), Set = set, `dbMEMs (positive)` = ncol(mem_sites),
            `Global test given geography, project and host: adjusted R² (%), p` =
              sprintf("%.1f, %s", 100 * mem_adj_r2, fmt_p(mem_p)),
            `Selected` = map_chr(mem_sel, ~ if (length(.x)) paste(.x, collapse = ", ") else "-"))

pub_margin <- margins %>%
  filter(set == "all") %>%
  left_join(margins %>% filter(set == "no_mixed") %>% select(k, var, p_nm = p), by = c("k", "var")) %>%
  transmute(Level = lab_k(k), Variable = soil_vars$label[match(var, soil_vars$var)], F = sprintf("%.2f", F),
            p = fmt_p(p), `Largest VIF` = sprintf("%.2f", vif), `Without mixed MAT: p` = fmt_p(p_nm))

pub_varpart <- vparts %>%
  filter(set == "all") %>%
  transmute(Level = lab_k(k),
            Fraction = recode(fraction, soil_only = "Soil only", depth_only = "Depth only",
                              context_only = "Site context only", soil_depth = "Soil and depth",
                              depth_context = "Depth and site context", soil_context = "Soil and site context",
                              all_three = "All three", residual = "Residual"),
            `Adjusted R² (%)` = sprintf("%.1f", 100 * adj_r2), p = if_else(is.na(p), "", fmt_p(p)))

pub_within <- within_res %>%
  filter(set == "all") %>%
  left_join(within_res %>% filter(set == "no_mixed") %>% select(k, depth_p_nm = depth_p, soil_p_nm = soil_p),
            by = "k") %>%
  transmute(Level = lab_k(k), Sites = sites,
            `Depth: adjusted R² (%), p` = sprintf("%.1f, %s", 100 * depth_adj_r2, fmt_p(depth_p)),
            `Soil: adjusted R² (%), p` = sprintf("%.1f, %s", 100 * soil_adj_r2, fmt_p(soil_p)),
            `Without mixed MAT: p depth, p soil` = sprintf("%s, %s", fmt_p(depth_p_nm), fmt_p(soil_p_nm)))

pub_profiles <- profiles %>%
  filter(set == "all") %>%
  select(k, unit, mlg, chisq, p, q, depth, share, lo, hi) %>%
  mutate(cell = sprintf("%.2f [%.2f, %.2f]", share, lo, hi), depth = lab_depth[depth]) %>%
  select(-share, -lo, -hi) %>%
  pivot_wider(names_from = depth, values_from = cell) %>%
  left_join(profiles %>% filter(set == "no_mixed") %>% distinct(k, unit, q_nm = q), by = c("k", "unit")) %>%
  transmute(Level = lab_k(k), Unit = unit, MLGs = mlg, !!!syms(unname(lab_depth)),
            `LRT chi² (df 2)` = sprintf("%.2f", chisq), p = fmt_p(p), q = fmt_p(q),
            `Without mixed MAT: q` = fmt_p(q_nm))

pub_sens <- selection %>%
  transmute(Level = lab_k(k), `Candidate variables` = lab_var(fwd_vars),
            `Global adjusted R² (%)` = map_chr(res, ~ sprintf("%.1f", 100 * .x$global$adj_r2)),
            `Global p` = map_chr(res, ~ fmt_p(.x$global$p)),
            `Selected (in order)` = map_chr(res, ~ if (length(.x$sel)) lab_var(.x$sel) else "-"),
            `Largest VIF` = map_chr(res, ~ if (length(.x$vif)) sprintf("%.2f", max(.x$vif)) else "-"))

pub_vars <- soil_vars %>%
  left_join(between, by = "var") %>%
  left_join(select(transforms, var, skew_raw, skew_ln), by = "var") %>%
  transmute(Variable = label, `Soil table column` = column, Transformation = if_else(log, "ln", "none"),
            `Skewness raw` = sprintf("%.2f", skew_raw),
            `Skewness ln` = if_else(is.na(skew_ln), "-", sprintf("%.2f", skew_ln)),
            Group = group, Kept = if_else(main, "yes", "no"),
            `Measured at` = if_else(var %in% vars_all_sites, "all sites", "PhytOakmeter sites"),
            `Variation between sites (%)` = sprintf("%.0f", 100 * between_sites))

write_tsv(map_dfr(k_levels, ~ comp(.x, "all") %>% mutate(k = .x, .before = 1)),
          file.path(out_root, "soil_niche_cells.tsv"))
write_tsv(transforms, file.path(out_root, "soil_transformations.tsv"))
write_tsv(as_tibble(soil_cor, rownames = "var"), file.path(out_root, "soil_correlations.tsv"))
write_tsv(between, file.path(out_root, "soil_between_sites.tsv"))
write_tsv(tests, file.path(out_root, "soil_niche_tests.tsv"))
write_tsv(tests_proj, file.path(out_root, "soil_niche_project_tests.tsv"))
write_tsv(margins, file.path(out_root, "soil_niche_margins.tsv"))
write_tsv(vparts, file.path(out_root, "soil_niche_varpart.tsv"))
write_tsv(within_res, file.path(out_root, "soil_niche_within_sites.tsv"))
write_tsv(profiles, file.path(out_root, "depth_profiles.tsv"))
write_tsv(scores_tbl, file.path(out_root, "soil_niche_scores.tsv"))
write_tsv(pub_test, file.path(out_root, "publication_soil_niche.tsv"))
write_tsv(pub_project, file.path(out_root, "publication_soil_niche_projects.tsv"))
write_tsv(pub_margin, file.path(out_root, "publication_soil_niche_margins.tsv"))
write_tsv(pub_varpart, file.path(out_root, "publication_soil_niche_varpart.tsv"))
write_tsv(pub_within, file.path(out_root, "publication_soil_niche_within_sites.tsv"))
write_tsv(pub_profiles, file.path(out_root, "publication_depth_profiles.tsv"))
write_tsv(pub_sens, file.path(out_root, "publication_soil_niche_forward_selection.tsv"))
write_tsv(pub_vars, file.path(out_root, "publication_soil_variables.tsv"))
write_tsv(pub_mem, file.path(out_root, "publication_soil_niche_dbmem.tsv"))
write_tsv(as_tibble(mem_sites, rownames = "site"), file.path(out_root, "soil_niche_dbmem_sites.tsv"))



# summary ----------------------------------------------------------------------

cat("\n=== Soil variables: transformations, groups and variation between sites ===\n")
print(pub_vars, n = Inf, width = Inf)
cat("largest |r| among kept variables:", round(r_kept, 2), "\n")

cat("\n=== Across sites: kept soil variables, responses per project ===\n")
print(pub_test, width = Inf)

cat("\n=== Across sites: marginal tests ===\n")
print(pub_margin, n = Inf, width = Inf)

cat("\n=== Across sites: variance partitioning (all isolates) ===\n")
print(pub_varpart, n = Inf, width = Inf)

cat("\n=== Per project ===\n")
print(pub_project, n = Inf, width = Inf)

cat("\n=== Within the PhytOakmeter sites: depth and soil ===\n")
print(pub_within, width = Inf)

cat("\n=== Depth profiles per unit (PhytOakmeter sites) ===\n")
print(pub_profiles, n = Inf, width = Inf)

cat("\n=== Sensitivity: forward selection ===\n")
print(pub_sens, n = Inf, width = Inf)

cat("\n=== Sensitivity: dbMEM selection ===\n")
print(pub_mem, n = Inf, width = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  soil_transformations.pdf\n")
cat("  soil_niche.pdf\n")
cat("  depth_profiles.pdf\n")
cat("  soil_variable_selection.pdf\n")
cat("  soil_niche_cells.tsv\n")
cat("  soil_transformations.tsv\n")
cat("  soil_correlations.tsv\n")
cat("  soil_between_sites.tsv\n")
cat("  soil_niche_tests.tsv\n")
cat("  soil_niche_project_tests.tsv\n")
cat("  soil_niche_margins.tsv\n")
cat("  soil_niche_varpart.tsv\n")
cat("  soil_niche_within_sites.tsv\n")
### depth_profiles.tsv is read by local/03_delimitation/060_delimitation_summary.R (same machine, no copy)
cat("  depth_profiles.tsv\n")
cat("  soil_niche_scores.tsv\n")
cat("  publication_soil_niche.tsv\n")
cat("  publication_soil_niche_projects.tsv\n")
cat("  publication_soil_niche_margins.tsv\n")
cat("  publication_soil_niche_varpart.tsv\n")
cat("  publication_soil_niche_within_sites.tsv\n")
cat("  publication_depth_profiles.tsv\n")
cat("  publication_soil_niche_forward_selection.tsv\n")
cat("  publication_soil_variables.tsv\n")
cat("  publication_soil_niche_dbmem.tsv\n")
cat("  soil_niche_dbmem_sites.tsv\n")
