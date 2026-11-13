## ============================================================================
## 051_bpp.R
## theta, tau and gdi per unit at K = 4, 5 and 6 (BPP A00) and posterior probabilities of
## delimitations and species trees (BPP A11), three disjoint locus sets and two chains each
##
## Files needed:
##   runs.tsv                   one row per BPP run
##   individuals.tsv            selected individuals, site and unit at every level
##   priors.tsv                 inverse-gamma priors for theta and root tau
##   loci.tsv                   set, locus, sites, variable sites, most missing per locus
##   tree_kX.nwk                unit tree per level (A00 topology)
##                              (all five from hpc/03_delimitation/050_bpp_input.sh;
##                               copy them from ${ROOT}/03_delimitation/results/050_results_bpp_input
##                               to local/03_delimitation/data/050_data_bpp_input)
##   <run>/<run>.ctl            control file
##   <run>/<run>.mcmc.txt       A00 MCMC sample
##   <run>/<run>.summary.txt    A11 posterior summary
##   <run>/extended_<samples>   marker of a finished extension
##                              (all four from hpc/03_delimitation/051_bpp.sh;
##                               copy every run folder from ${ROOT}/03_delimitation/results/051_results_bpp
##                               to local/03_delimitation/data/051_data_bpp)
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)
library(coda)



# files ------------------------------------------------------------------------

input_root <- "03_delimitation/data/050_data_bpp_input"
data_root  <- "03_delimitation/data/051_data_bpp"
out_root   <- "03_delimitation/results/051_results_bpp"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### gdi thresholds (below: one lineage, above: separate species, between: ambiguous)
gdi_low  <- 0.2
gdi_high <- 0.7

### convergence: smallest effective sample size of the pooled chains, largest PSRF between chains
ess_min  <- 200
psrf_max <- 1.1



# style ------------------------------------------------------------------------

### unit colours
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

### locus sets
shape_set <- c(`1` = 21, `2` = 24, `3` = 22)

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}

### mean [95% HPD]
fmt_hpd <- function(m, lo, hi, digits = 3) sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
                                                   m, lo, hi)
fmt_sci <- function(m, lo, hi) sprintf("%.2e [%.2e, %.2e]", m, lo, hi)



# read -------------------------------------------------------------------------

runs   <- read_tsv(file.path(input_root, "runs.tsv"), col_types = "ccddcd")
ind    <- read_tsv(file.path(input_root, "individuals.tsv"), col_types = cols(.default = "c"))
priors <- read_tsv(file.path(input_root, "priors.tsv"), col_types = "cddd")
loci   <- read_tsv(file.path(input_root, "loci.tsv"), col_types = "dcddd")

### MCMC settings from the first control file
ctl <- read_lines(file.path(data_root, runs$run[1], paste0(runs$run[1], ".ctl")))
ctl_value <- function(key) str_trim(str_remove(str_subset(ctl, paste0("^\\s*", key, "\\s*=")), ".*=\\s*"))
mcmc_set  <- as.numeric(c(burnin = ctl_value("burnin"), sampfreq = ctl_value("sampfreq"),
                          nsample = ctl_value("nsample")))
names(mcmc_set) <- c("burnin", "sampfreq", "nsample")

### A00 samples, one row per sample, run and parameter
a00_runs <- filter(runs, analysis == "A00")

a00 <- a00_runs %>%
  mutate(data = map(run, ~ read_tsv(file.path(data_root, .x, paste0(.x, ".mcmc.txt")),
                                    col_types = cols(.default = "d"), progress = FALSE))) %>%
  mutate(data = map(data, ~ .x %>%
                      select(Gen, starts_with("theta:"), starts_with("tau:")) %>%
                      pivot_longer(-Gen, names_to = "column", values_to = "value"))) %>%
  unnest(data) %>%
  separate_wider_delim(column, ":", names = c("parameter", "node", "members"), cols_remove = FALSE)

### A11 summaries
read_a11 <- function(path) {
  s <- read_lines(path)
  sec <- function(pattern) {
    i <- str_which(s, pattern)[1] + 1
    e <- which(str_trim(s[i:length(s)]) == "")
    j <- if (length(e) > 0) i + e[1] - 2 else length(s)
    str_squish(s[i:j])
  }
  best   <- sec("^\\(A\\) List of best models")[1]
  nsp    <- sec("^\\(D\\) Posterior probability for # of species")
  delim  <- sec("^\\(C\\) .*delimited species")
  list(best_tree  = str_extract(best, "\\(.*;$") %>% str_remove("^\\([^()]*\\)\\s*"),
       best_pp    = as.numeric(str_split_i(best, " ", 2)),
       n_species  = tibble(n = as.integer(str_extract(nsp, "(?<=P\\[)[0-9]+")),
                           pp = as.numeric(str_extract(nsp, "(?<=\\] = )[0-9.]+"))),
       species_pp = tibble(species = str_split_i(delim, " ", 3), pp = as.numeric(str_split_i(delim, " ", 2))))
}

a11_runs <- filter(runs, analysis == "A11") %>%
  mutate(res = map(run, ~ read_a11(file.path(data_root, .x, paste0(.x, ".summary.txt")))))



# checks -----------------------------------------------------------------------

## check if every A00 run has its full sample (nsample of the control file plus the extended_<samples> markers)
n_expected <- a00_runs %>%
  transmute(run, expected = mcmc_set[["nsample"]] +
              map_dbl(run, ~ sum(as.numeric(str_remove(list.files(file.path(data_root, .x),
                                                                  pattern = "^extended_[0-9]+$"),
                                                       "^extended_")))))
n_sample <- a00 %>% distinct(run, Gen) %>% count(run) %>% right_join(n_expected, by = "run")

if (any(is.na(n_sample$n) | n_sample$n != n_sample$expected)) {
  warning("A00 runs missing, incomplete or still extending: ",
          paste(n_sample$run[is.na(n_sample$n) | n_sample$n != n_sample$expected], collapse = ", "))
}

## check if every A00 run has a theta for every unit of its level
n_units <- ind %>% pivot_longer(starts_with("k"), names_to = "k", values_to = "unit") %>%
  distinct(k, unit) %>% count(k, name = "units")
n_theta <- a00 %>% filter(parameter == "theta", !str_detect(members, ",")) %>%
  distinct(run, k, members) %>% count(run, k) %>% left_join(n_units, by = "k")

if (any(n_theta$n != n_theta$units)) {
  warning("A00 runs without a theta for every unit")
}

## check if every A11 run was summarised
if (any(map_lgl(a11_runs$res, ~ is.na(.x$best_pp)))) {
  warning("A11 runs without a delimitation summary")
}



# convergence ------------------------------------------------------------------

### ESS per run and parameter
ess <- a00 %>%
  summarise(ess = as.numeric(effectiveSize(value)), .by = c(k, set, chain, run, parameter, members))

### PSRF between the two chains per level, locus set and parameter
psrf <- a00 %>%
  summarise(psrf = {
    ch <- split(value, chain)
    if (length(ch) < 2) NA_real_ else
      gelman.diag(mcmc.list(lapply(ch, mcmc)), autoburnin = FALSE, multivariate = FALSE)$psrf[1, 1]
  }, .by = c(k, set, parameter, members))

### ESS of the pooled chains (the sample used for all estimates; sum over chains)
ess_pooled <- ess %>%
  summarise(ess_pooled = sum(ess), .by = c(k, set, parameter, members))

### convergence per level and locus set
conv <- ess %>%
  summarise(ess_chain_min = min(ess), .by = c(k, set)) %>%
  left_join(ess_pooled %>% summarise(ess_min = min(ess_pooled), .by = c(k, set)), by = c("k", "set")) %>%
  left_join(psrf %>% summarise(psrf_max = max(psrf, na.rm = TRUE), .by = c(k, set)), by = c("k", "set"))

## check if the runs reach the convergence thresholds
if (any(conv$ess_min < ess_min) || any(conv$psrf_max > psrf_max)) {
  warning("convergence below the thresholds (ESS < ", ess_min, " or PSRF > ", psrf_max,
          "); extend the runs concerned with bpp --resume <last .chk> --extend <samples>")
}



# posterior parameters ---------------------------------------------------------

### posterior mean and 95% HPD per level, locus set and parameter (chains pooled)
hpd <- function(x) as.numeric(HPDinterval(mcmc(x), prob = 0.95))

params <- a00 %>%
  summarise(mean = mean(value), lo = hpd(value)[1], hi = hpd(value)[2],
            .by = c(k, set, parameter, node, members)) %>%
  mutate(size = str_count(members, ",") + 1)



# gdi --------------------------------------------------------------------------

### every unit with its parent node (the smallest tau clade holding it) and its sister clade
parents <- params %>%
  distinct(k, parameter, members, size) %>%
  (\(p) p %>%
     filter(parameter == "theta", size == 1) %>%
     select(k, unit = members) %>%
     left_join(filter(p, parameter == "tau") %>% select(k, parent = members, parent_size = size),
               by = "k", relationship = "many-to-many") %>%
     filter(map2_lgl(unit, parent, ~ .x %in% str_split_1(.y, ","))) %>%
     slice_min(parent_size, by = c(k, unit)) %>%
     mutate(sister = map2_chr(unit, parent, ~ paste(sort(setdiff(str_split_1(.y, ","), .x)), collapse = "+"))))()

### gdi = 1 - exp(-2 tau_parent / theta_unit), per MCMC sample (chains pooled)
gdi_samples <- a00 %>%
  filter(parameter == "theta", !str_detect(members, ",")) %>%
  select(k, set, chain, Gen, unit = members, theta = value) %>%
  left_join(select(parents, k, unit, parent), by = c("k", "unit")) %>%
  left_join(a00 %>% filter(parameter == "tau") %>% select(k, set, chain, Gen, parent = members, tau = value),
            by = c("k", "set", "chain", "Gen", "parent")) %>%
  mutate(gdi = 1 - exp(-2 * tau / theta))

gdi <- gdi_samples %>%
  summarise(gdi_mean = mean(gdi), gdi_lo = hpd(gdi)[1], gdi_hi = hpd(gdi)[2],
            theta_mean = mean(theta), theta_lo = hpd(theta)[1], theta_hi = hpd(theta)[2],
            tau_mean = mean(tau), tau_lo = hpd(tau)[1], tau_hi = hpd(tau)[2],
            .by = c(k, set, unit, parent)) %>%
  left_join(select(parents, k, unit, sister), by = c("k", "unit")) %>%
  arrange(k, unit, set)

### call per unit from all locus sets (thresholds on the posterior means)
gdi_call <- gdi %>%
  summarise(call = case_when(all(gdi_mean > gdi_high) ~ "separate",
                             all(gdi_mean < gdi_low)  ~ "one lineage",
                             TRUE                     ~ "ambiguous"),
            .by = c(k, unit))



# A11 --------------------------------------------------------------------------

### loci per A11 run (from its control file)
a11_loci <- function(run) {
  l <- read_lines(file.path(data_root, run, paste0(run, ".ctl")))
  as.integer(str_extract(str_subset(l, "^\\s*nloci\\s*="), "[0-9]+$"))
}

a11 <- a11_runs %>%
  transmute(set, chain,
            n_loci    = map_int(run, a11_loci),
            best_tree = map_chr(res, "best_tree"),
            best_pp   = map_dbl(res, "best_pp"),
            n_species = map_int(res, ~ .x$n_species$n[which.max(.x$n_species$pp)]),
            pp_n      = map_dbl(res, ~ max(.x$n_species$pp)),
            species   = map_chr(res, ~ paste(sprintf("%s %.2f", .x$species_pp$species, .x$species_pp$pp),
                                             collapse = ", ")))



# figures ----------------------------------------------------------------------

### one panel per level: gdi mean and 95% HPD per unit and locus set
plot_level <- function(kk, y_title) {
  d <- gdi %>% filter(k == kk) %>% mutate(set = factor(set))

  ggplot(d, aes(x = unit, y = gdi_mean, fill = unit, shape = set, group = set)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = gdi_low, ymax = gdi_high, fill = "grey93") +
    geom_errorbar(aes(ymin = gdi_lo, ymax = gdi_hi), width = 0, linewidth = 0.3,
                  colour = "grey30", position = position_dodge(width = 0.6)) +
    geom_point(size = 2, stroke = 0.3, colour = "black", position = position_dodge(width = 0.6)) +
    scale_fill_manual(values = pal_unit, guide = "none") +
    scale_shape_manual(values = shape_set, name = "Locus set") +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    guides(shape = guide_legend(override.aes = list(fill = "grey70"))) +
    labs(x = sprintf("Units at K = %s", sub("k", "", kk)), y = if (y_title) "gdi" else NULL)
}

levels_k <- sort(unique(gdi$k))
p_gdi <- wrap_plots(map2(levels_k, seq_along(levels_k) == 1, plot_level), nrow = 1,
                    widths = table(distinct(gdi, k, unit)$k)[levels_k]) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

save_fig("bpp_gdi.pdf", p_gdi, 7.1, 3)


### the parameter with the largest PSRF per level and locus set, both chains
worst <- psrf %>%
  slice_max(psrf, n = 1, with_ties = FALSE, by = c(k, set)) %>%
  mutate(label = sprintf("%s, set %s: %s %s (PSRF %.2f)", k, set, parameter, members, psrf))

p_trace <- a00 %>%
  inner_join(select(worst, k, set, parameter, members, label), by = c("k", "set", "parameter", "members")) %>%
  mutate(chain = factor(chain)) %>%
  ggplot(aes(x = Gen, y = value, colour = chain)) +
  geom_line(linewidth = 0.15, alpha = 0.7) +
  facet_wrap(~ label, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c(`1` = "grey20", `2` = "#E69F00"), name = "Chain") +
  labs(x = "MCMC step", y = "Parameter value") +
  theme(strip.text = element_text(size = 6))

save_fig("bpp_trace_worst_psrf.pdf", p_trace, 7.1, 6)


### unit trees dated by tau: node at the mean over locus sets, 95% HPD per locus set
tau_nodes <- params %>%
  filter(parameter == "tau") %>%
  mutate(key = map_chr(str_split(members, ","), ~ paste(sort(.x), collapse = ",")))

tau_max <- max(tau_nodes$hi)

plot_tree <- function(kk) {
  tr    <- ape::reorder.phylo(ape::ladderize(ape::read.tree(file.path(input_root, sprintf("tree_%s.nwk", kk)))),
                              "cladewise")
  n_tip <- ape::Ntip(tr)
  inner <- n_tip + seq_len(tr$Nnode)
  par   <- tr$edge[, 1]
  chi   <- tr$edge[, 2]

  ### node keys as in the BPP labels (sorted tips)
  node_key <- c(tr$tip.label, map_chr(inner, ~ paste(sort(ape::extract.clade(tr, .x)$tip.label), collapse = ",")))
  tk       <- filter(tau_nodes, k == kk)

  ## check if every inner node of the unit tree has a tau
  if (!all(node_key[inner] %in% tk$key)) warning("nodes of tree_", kk, " without a tau")

  node_age <- c(rep(0, n_tip), map_dbl(node_key[inner], ~ mean(tk$mean[tk$key == .x])))

  ### tips in plotting order, inner nodes at the mean of their children
  node_y <- rep(NA_real_, n_tip + tr$Nnode)
  node_y[chi[chi <= n_tip]] <- seq_len(n_tip)
  while (anyNA(node_y[inner])) {
    for (nd in inner[is.na(node_y[inner])]) {
      if (!anyNA(node_y[chi[par == nd]])) node_y[nd] <- mean(node_y[chi[par == nd]])
    }
  }

  seg_h <- tibble(x = node_age[par], xend = node_age[chi], y = node_y[chi])
  seg_v <- tibble(x    = node_age[inner],
                  ymin = map_dbl(inner, ~ min(node_y[chi[par == .x]])),
                  ymax = map_dbl(inner, ~ max(node_y[chi[par == .x]])))
  hpd_n <- tk %>%
    mutate(y = node_y[match(key, node_key)] + (set - 2) * 0.15, set = factor(set))
  tips  <- tibble(unit = tr$tip.label, y = node_y[seq_len(n_tip)])

  ggplot() +
    geom_segment(data = seg_h, aes(x = x, xend = xend, y = y, yend = y), linewidth = 0.4, colour = "grey30") +
    geom_segment(data = seg_v, aes(x = x, xend = x, y = ymin, yend = ymax), linewidth = 0.4, colour = "grey30") +
    geom_linerange(data = hpd_n, aes(xmin = lo, xmax = hi, y = y), linewidth = 0.3, colour = "grey45") +
    geom_point(data = hpd_n, aes(x = mean, y = y, shape = set), size = 1.4, stroke = 0.3,
               fill = "grey70", colour = "black") +
    geom_point(data = tips, aes(x = 0, y = y, fill = unit), shape = 21, size = 2.2, stroke = 0.3,
               colour = "black") +
    geom_text(data = tips, aes(x = 0, y = y, label = unit), hjust = 0, nudge_x = 0.035 * tau_max, size = 2.6) +
    scale_x_reverse(limits = c(tau_max * 1.02, -0.14 * tau_max), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(add = 0.5)) +
    scale_fill_manual(values = pal_unit, guide = "none") +
    scale_shape_manual(values = shape_set, name = "Locus set", drop = FALSE) +
    labs(x = expression(tau ~ "(expected substitutions per site)"), y = sprintf("K = %s", sub("k", "", kk))) +
    theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
}

n_units  <- table(distinct(gdi, k, unit)$k)[levels_k]
p_dated <- wrap_plots(map(levels_k, plot_tree), ncol = 1, heights = n_units) +
  plot_layout(guides = "collect", axes = "collect_x") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

save_fig("bpp_dated_trees.pdf", p_dated, 3.4, 6.5)



# tables -----------------------------------------------------------------------

### gdi per level and unit, one column per locus set
pub_gdi <- gdi %>%
  mutate(cell = fmt_hpd(gdi_mean, gdi_lo, gdi_hi, 2)) %>%
  select(k, unit, sister, set, cell) %>%
  pivot_wider(names_from = set, values_from = cell, names_prefix = "gdi, locus set ") %>%
  left_join(gdi_call, by = c("k", "unit")) %>%
  transmute(Level = sub("k", "K = ", k), Unit = unit, Sister = sister,
            across(starts_with("gdi, locus set")), Call = call)

### theta and tau per level and node, one column per locus set
pub_params <- params %>%
  mutate(cell = fmt_sci(mean, lo, hi), Node = map_chr(str_split(members, ","), ~ paste(sort(.x), collapse = "+"))) %>%
  select(k, parameter, Node, size, set, cell) %>%
  pivot_wider(names_from = set, values_from = cell, names_prefix = "Locus set ") %>%
  arrange(k, desc(parameter), size, Node) %>%
  transmute(Level = sub("k", "K = ", k), Parameter = parameter, Node, across(starts_with("Locus set")))

### A11 per locus set and chain
pub_a11 <- a11 %>%
  transmute(`Locus set` = set, Chain = chain, Loci = n_loci, `Species (posterior mode)` = n_species,
            `P(number of species)` = sprintf("%.3f", pp_n), `Best species tree` = best_tree,
            `P(best model)` = sprintf("%.3f", best_pp), `P per delimited species` = species)

write_tsv(gdi, file.path(out_root, "bpp_gdi.tsv"))
write_tsv(params, file.path(out_root, "bpp_parameters.tsv"))
write_tsv(ess %>% left_join(ess_pooled, by = c("k", "set", "parameter", "members")) %>%
            left_join(psrf, by = c("k", "set", "parameter", "members")),
          file.path(out_root, "bpp_convergence.tsv"))
write_tsv(pub_gdi, file.path(out_root, "publication_bpp_gdi.tsv"))
write_tsv(pub_params, file.path(out_root, "publication_bpp_parameters.tsv"))
write_tsv(pub_a11, file.path(out_root, "publication_bpp_a11.tsv"))



# summary ----------------------------------------------------------------------

cat("\n=== Convergence per level and locus set ===\n")
print(conv, n = Inf)

cat("\n=== Parameters above the PSRF threshold (in_gdi: tip theta or parent tau of a unit) ===\n")
psrf %>%
  filter(psrf > psrf_max) %>%
  mutate(in_gdi = (parameter == "theta" & !str_detect(members, ",")) |
           (parameter == "tau" & paste(k, members) %in% paste(parents$k, parents$parent))) %>%
  left_join(ess_pooled, by = c("k", "set", "parameter", "members")) %>%
  arrange(desc(psrf)) %>%
  print(n = Inf, width = Inf)

cat("\n=== gdi per level and unit (mean [95% HPD], one column per locus set) ===\n")
print(pub_gdi, n = Inf, width = Inf)

cat("\n=== A11 per locus set and chain ===\n")
print(select(a11, -species), n = Inf, width = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  bpp_gdi.pdf\n")
cat("  bpp_trace_worst_psrf.pdf\n")
cat("  bpp_dated_trees.pdf\n")
cat("  bpp_gdi.tsv\n")
cat("  bpp_parameters.tsv\n")
cat("  bpp_convergence.tsv\n")
cat("  publication_bpp_gdi.tsv\n")
cat("  publication_bpp_parameters.tsv\n")
cat("  publication_bpp_a11.tsv\n")
