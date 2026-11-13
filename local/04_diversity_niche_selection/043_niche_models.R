## ============================================================================
## 043_niche_models.R
## Soil niche per unit at K = 4, 5 and 6 within each project: probability that a
## genotype found at a site and depth belongs to the unit, Bayesian model
## averaging over all subsets of the in situ soil variables (binomial GLMMs)
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
library(lme4)



# files ------------------------------------------------------------------------

soil_file   <- "04_diversity_niche_selection/supplementary/soildata_ALL_FIXED.txt"
labels_file <- "03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv"
meta_file   <- "01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv"
out_root    <- "04_diversity_niche_selection/results/043_results_niche_models"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### levels and sample sets
k_levels <- c("k4", "k5", "k6")
sets     <- c("all", "no_mixed")

### projects; depth code of the single MAGIC sample (0-10 cm, empty in the metadata)
projects     <- c("magic", "phytoakmeter")
depth_single <- c(magic = "1")

### mating type calls of possible mixed cultures (sample set without them)
mat_exclude <- c("MAT12", "ambiguous")

### in situ soil variables: the kept set of the soil pRDA in 040 (column, name, ln-transformed)
soil_is <- tribble(
  ~column,     ~var,   ~log,
  "pH_H20",    "pH",   FALSE,
  "C",         "C",    TRUE,
  "C_N",       "CN",   FALSE,
  "Prop_Clay", "Clay", FALSE,
  "P2O5",      "P2O5", TRUE)

### smallest number of MLGs of a unit in a project for its models (about six per predictor with host)
n_min_model <- 40

### histogram bins of the diagnostics, seed
n_bins <- 15
seed   <- 42



# style ------------------------------------------------------------------------

### unit colours
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

### project labels and fills, predictor labels, colours of negative and positive effects
lab_project <- c(magic = "MAGIC", phytoakmeter = "PhytOakmeter")
pal_project <- c(magic = "white", phytoakmeter = "grey25")
lab_pred    <- c(pH = "pH", C = "ln C", CN = "C:N", Clay = "Clay", P2O5 = "ln P2O5")
pal_sign    <- c(negative = "#762A83", positive = "#1B7837")

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}

### p values in tables, level labels
fmt_p <- function(p) if_else(p < 0.001, "< 0.001", sprintf("%.3f", p))
lab_k <- function(k) sub("k", "K = ", k)

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

### in situ soil per site and depth (site codes with "-" as in the metadata)
soil <- read_tsv(soil_file, col_types = cols(Site = "c", .default = "d")) %>%
  mutate(site = str_replace_all(Site, "_", "-"), depth = as.character(Depth))

in_situ <- soil %>%
  transmute(site, depth,
            !!!setNames(map2(soil_is$column, soil_is$log,
                             function(cl, lg) if (lg) log(soil[[cl]]) else soil[[cl]]),
                        soil_is$var))

### one row per MLG and site and depth (PhytOakmeter directions pooled), with project, host and soil
genotypes <- lab %>%
  select(sample, all_of(k_levels)) %>%
  inner_join(meta %>% filter(project %in% projects) %>%
               select(sample, project, site, depth = sampling_depth, mlg, mat_genome), by = "sample") %>%
  mutate(depth = coalesce(depth, depth_single[project]), mixed = mat_genome %in% mat_exclude) %>%
  distinct(site, depth, mlg, .keep_all = TRUE) %>%
  left_join(select(site_info, site, host), by = "site") %>%
  inner_join(in_situ, by = c("site", "depth"))



# checks -----------------------------------------------------------------------

## check if every in situ variable is in the soil table
if (!all(soil_is$column %in% names(soil))) {
  stop("soil columns missing: ", paste(setdiff(soil_is$column, names(soil)), collapse = ", "))
}

## check if the in situ variables are complete
if (anyNA(soil[soil_is$column])) {
  stop("missing values in the kept soil variables")
}

## check if the ln-transformed variables are positive
if (any(map_lgl(soil_is$column[soil_is$log], ~ any(soil[[.x]] <= 0)))) {
  stop("zero or negative values in ln-transformed soil variables")
}

## check if every site has one project and one host
if (anyDuplicated(site_info$site) || anyNA(genotypes$host)) {
  stop("sites without exactly one project and host in the metadata")
}

## check if the smallest modelled unit has at least five genotypes per predictor (soil and host)
if (n_min_model / (nrow(soil_is) + 1) < 5) {
  warning("fewer than five genotypes of the unit per predictor")
}



# models -----------------------------------------------------------------------

### binomial GLMM with random intercepts for site and depth within site (where sites hold several depths;
### genotypes of one site and depth share their predictor values), host where it varies; fits with warnings flagged
fit_glmm <- function(d, x) {
  re   <- if (any(count(distinct(d, site, depth), site)$n > 1)) "(1 | site/depth)" else "(1 | site)"
  host <- if (n_distinct(d$host) > 1) "host" else character(0)
  warn <- FALSE
  m <- withCallingHandlers(
    glmer(reformulate(c(host, if (length(x) > 0) x else "1", re), response = "y"), family = binomial,
          data = d, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))),
    warning = function(w) { warn <<- TRUE; invokeRestart("muffleWarning") },
    message = function(m) invokeRestart("muffleMessage"))
  list(model = m, warning = warn)
}

### all subsets of the soil variables (the last one holds all of them)
bma_vars <- soil_is$var
subsets  <- unlist(map(0:length(bma_vars), ~ combn(bma_vars, .x, simplify = FALSE)), recursive = FALSE)

### Bayesian model averaging for one unit: BIC weights as posterior model probabilities, posterior inclusion
### probability (PIP) = summed weight of the models holding a variable, model-averaged effect over all models
### (zero where absent); effects, intervals and p values of the full model
bma_unit <- function(d) {
  fits <- map(subsets, ~ fit_glmm(d, .x))
  bic  <- map_dbl(fits, ~ BIC(.x$model))
  w    <- exp(-0.5 * (bic - min(bic)))
  w    <- w / sum(w)
  b    <- map(fits, ~ fixef(.x$model))
  full <- fits[[length(fits)]]$model
  s    <- summary(full)$coefficients
  vc   <- as.data.frame(VarCorr(full))
  tibble(var      = bma_vars,
         pip      = map_dbl(bma_vars, function(v) sum(w[map_lgl(subsets, ~ v %in% .x)])),
         est_avg  = map_dbl(bma_vars, function(v) sum(w * map_dbl(b, ~ if (v %in% names(.x)) .x[[v]] else 0))),
         est      = s[bma_vars, 1],
         se       = s[bma_vars, 2],
         p        = s[bma_vars, 4],
         best     = paste(subsets[[which.max(w)]], collapse = " + "),
         w_best   = max(w),
         sd_site  = vc$sdcor[vc$grp == "site"],
         sd_depth = if (any(vc$grp == "depth:site")) vc$sdcor[vc$grp == "depth:site"] else NA_real_,
         singular = isSingular(full),
         warnings = sum(map_lgl(fits, "warning")))
}

### genotypes of a project and sample set, soil standardised within the project
model_rows <- function(pp, ss) {
  d <- filter(genotypes, project == pp)
  if (ss == "no_mixed") d <- filter(d, !mixed)
  mutate(d, across(all_of(bma_vars), ~ as.numeric(scale(.x))))
}

### units with enough MLGs per project, level and sample set
targets <- expand_grid(project = projects, k = k_levels, set = sets) %>%
  mutate(unit = pmap(list(project, k, set), function(pp, kk, ss) {
    n <- table(model_rows(pp, ss)[[kk]])
    intersect(names(pal_unit), names(n)[n >= n_min_model])
  })) %>%
  unnest(unit)

set.seed(seed)
bma <- targets %>%
  mutate(res = pmap(list(project, k, set, unit), function(pp, kk, ss, uu) {
    d <- model_rows(pp, ss)
    bma_unit(mutate(d, y = as.integer(d[[kk]] == uu)))
  })) %>%
  unnest(res)



# figures ----------------------------------------------------------------------

### transformation diagnostics (supplement): ln-transformed predictors on the raw and ln scale at the genotype
### rows of the models, with skewness per project
diag_long <- soil_is %>%
  filter(log) %>%
  mutate(data = map(column, function(cl) {
    x <- soil[[cl]][match(paste(genotypes$site, genotypes$depth), paste(soil$site, soil$depth))]
    bind_rows(tibble(project = genotypes$project, scale = "raw", value = x),
              tibble(project = genotypes$project, scale = "ln", value = log(x)))
  })) %>%
  unnest(data) %>%
  mutate(scale = factor(scale, levels = c("raw", "ln")))

diag_skew <- diag_long %>%
  summarise(skew = skew(value), .by = c(var, scale, project)) %>%
  summarise(text = paste(sprintf("%s %.2f", lab_project[project], skew), collapse = ", "), .by = c(var, scale))

p_diag <- diag_long %>%
  left_join(diag_skew, by = c("var", "scale")) %>%
  mutate(panel = sprintf("%s, %s\nskewness %s", var, scale, text),
         panel = fct_reorder(panel, match(var, soil_is$var) + (scale == "ln") / 2)) %>%
  ggplot(aes(value, fill = project)) +
  geom_histogram(bins = n_bins, colour = "grey20", linewidth = 0.2) +
  facet_wrap(~ panel, scales = "free", ncol = 2) +
  scale_fill_manual(values = pal_project, labels = lab_project, name = NULL) +
  labs(x = NULL, y = "Genotype rows") +
  theme(strip.background = element_blank(), strip.text = element_text(size = 6))

save_fig("niche_models_transformations.pdf", p_diag, 7.1, 4)

### soil niche per unit and project (all isolates; supplement): full-model effect with 95% interval and
### significance, point size the PIP, colour the sign; one row per level
p_bma <- bma %>%
  filter(set == "all") %>%
  mutate(predictor = factor(var, levels = rev(soil_is$var), labels = rev(lab_pred[soil_is$var])),
         sign = if_else(est > 0, "positive", "negative"),
         stars = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", .default = ""),
         panel = paste0(lab_k(k), ": ", unit, " (", lab_project[project], ")")) %>%
  mutate(panel = factor(panel, levels = unique(panel[order(k, project, match(unit, names(pal_unit)))]))) %>%
  ggplot(aes(x = est, y = predictor)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
  geom_errorbar(aes(xmin = est - 1.96 * se, xmax = est + 1.96 * se, colour = sign), width = 0, linewidth = 0.4,
                orientation = "y") +
  geom_point(aes(size = pip, fill = sign), shape = 21, stroke = 0.3, colour = "black") +
  geom_text(aes(label = stars), vjust = -0.3, size = 9 / .pt, show.legend = FALSE) +
  facet_wrap(~ panel, ncol = 5) +
  scale_colour_manual(values = pal_sign, labels = c(negative = "Negative", positive = "Positive"), name = NULL) +
  scale_fill_manual(values = pal_sign, labels = c(negative = "Negative", positive = "Positive"), name = NULL) +
  scale_size_area(limits = c(0, 1), breaks = c(0.25, 0.5, 0.75, 1), max_size = 3.5, name = "PIP") +
  guides(size = guide_legend(override.aes = list(fill = "grey70"))) +
  labs(x = "Standardised effect (log odds per SD within project)", y = NULL) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 6.5), legend.position = "right")

save_fig("niche_models_bma.pdf", p_bma, 7.1, 7)



# tables -----------------------------------------------------------------------

pub_bma <- bma %>%
  filter(set == "all") %>%
  left_join(bma %>% filter(set == "no_mixed") %>% select(project, k, unit, var, pip_nm = pip),
            by = c("project", "k", "unit", "var")) %>%
  mutate(n = pmap_int(list(project, k, unit),
                      function(pp, kk, uu) sum(genotypes$project == pp & genotypes[[kk]] == uu))) %>%
  transmute(Level = lab_k(k), Project = lab_project[project], Unit = unit, MLGs = n, Variable = lab_pred[var],
            PIP = sprintf("%.2f", pip), `Model-averaged effect` = sprintf("%.2f", est_avg),
            `Full model: effect [95% CI]` = sprintf("%.2f [%.2f, %.2f]", est, est - 1.96 * se, est + 1.96 * se),
            p = fmt_p(p),
            `Best model (weight)` = sprintf("%s (%.2f)", if_else(best == "", "no soil variable", best), w_best),
            `Without mixed MAT: PIP` = if_else(is.na(pip_nm), "-", sprintf("%.2f", pip_nm)))

pub_transform <- diag_long %>%
  summarise(skew = skew(value), .by = c(var, scale, project)) %>%
  pivot_wider(names_from = scale, values_from = skew, names_prefix = "skew_") %>%
  transmute(Variable = lab_pred[var], Project = lab_project[project],
            `Skewness raw` = sprintf("%.2f", skew_raw), `Skewness ln` = sprintf("%.2f", skew_ln))

write_tsv(genotypes, file.path(out_root, "niche_models_data.tsv"))
write_tsv(bma, file.path(out_root, "niche_models_bma.tsv"))
write_tsv(pub_bma, file.path(out_root, "publication_niche_models_bma.tsv"))
write_tsv(pub_transform, file.path(out_root, "publication_niche_models_transformations.tsv"))



# summary ----------------------------------------------------------------------

cat("\n=== Transformations (genotype rows) ===\n")
print(pub_transform, width = Inf)

cat("\n=== Soil niche per unit and project: Bayesian model averaging ===\n")
print(pub_bma, n = Inf, width = Inf)

fits_bma <- distinct(bma, project, k, set, unit, singular, warnings)
cat("\nFull models singular:", sum(fits_bma$singular), "of", nrow(fits_bma), "; GLMM fits with warnings:",
    sum(fits_bma$warnings), "of", nrow(fits_bma) * length(subsets), "\n")



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  niche_models_transformations.pdf\n")
cat("  niche_models_bma.pdf\n")
cat("  niche_models_data.tsv\n")
cat("  niche_models_bma.tsv\n")
cat("  publication_niche_models_bma.tsv\n")
cat("  publication_niche_models_transformations.tsv\n")
