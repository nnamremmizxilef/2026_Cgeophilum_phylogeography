## ============================================================================
## 040_publication_delimitation.R
## Publication figure of lineage delimitation: unit tree of the main level dated
## by BPP tau (A), gdi per unit at K = 4, 5 and 6 (B), BPP theta against
## nucleotide diversity (C), f-branch between units (D), recombination within
## units (E) and every line of evidence per unit (F)
##
## Files needed:
##   bpp_parameters.tsv              theta and tau per level, locus set and node
##   bpp_gdi.tsv                     gdi and tip theta per unit, level and locus set
##                                   (both from local/03_delimitation/051_bpp.R)
##   nucleotide_diversity_units.tsv  pi per unit with bootstrap interval,
##                                   from local/04_diversity_niche_selection/011_nucleotide_diversity.R
##   dsuite_fbranch.tsv              f-branch per level and sample set,
##                                   from local/03_delimitation/030_dsuite.R
##   phi_summary_sets.tsv            genes tested and with PHI evidence per unit and sample set,
##                                   from local/03_delimitation/021_phi_test.R
##   index_association_sets.tsv      rbarD and its permutation p per unit and sample set,
##                                   from local/03_delimitation/022_index_association.R
##   delimitation_evidence.tsv       lines of evidence per level and unit,
##                                   from local/03_delimitation/060_delimitation_summary.R
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)
library(ggrepel)



# files ------------------------------------------------------------------------

res_03   <- "03_delimitation/results"
res_04   <- "04_diversity_niche_selection/results"
out_root <- "00_general/results/040_results_publication_delimitation"

f_in <- c(params = file.path(res_03, "051_results_bpp", "bpp_parameters.tsv"),
          gdi    = file.path(res_03, "051_results_bpp", "bpp_gdi.tsv"),
          pi     = file.path(res_04, "011_results_nucleotide_diversity", "nucleotide_diversity_units.tsv"),
          fb     = file.path(res_03, "030_results_dsuite", "dsuite_fbranch.tsv"),
          phi    = file.path(res_03, "021_results_phi_test", "phi_summary_sets.tsv"),
          ia     = file.path(res_03, "022_results_index_association", "index_association_sets.tsv"))

f_ev <- file.path(res_03, "060_results_delimitation_summary", "delimitation_evidence.tsv")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### level of panels A and C to E; levels of panel B
k_main   <- "k6"
k_levels <- c("k4", "k5", "k6")

### gdi thresholds (below: one lineage, above: separate species)
gdi_low  <- 0.2
gdi_high <- 0.7

### sample set of panels C to E (all MLG representatives) and the equal-size draws of 021 and 022
set_main  <- "all"
set_equal <- "equal_n"

### significance level of the PHI test in 021 (share of genes expected without recombination) and of rbarD
alpha_phi <- 0.05
alpha_ia  <- 0.05

### smallest f-branch printed in its cell
fb_label_min <- 0.01

### stem support and A11 posterior counted as support of a unit
ufboot_min <- 95
lpp_min    <- 0.95
a11_min    <- 0.95



# style ------------------------------------------------------------------------

### unit colours
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

### locus set symbols, f-branch ramp, support fills of the evidence cells
shape_set <- c(`1` = 21, `2` = 24, `3` = 22)
pal_fb    <- c("white", "grey20")
pal_call  <- c(Supports = "grey45", Ambiguous = "grey80", Against = "grey97")

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}

lab_k <- function(k) sub("k", "K = ", k)



# read -------------------------------------------------------------------------

## check if every input exists
miss <- f_in[!file.exists(f_in)]
if (length(miss) > 0) stop("missing inputs: ", paste(miss, collapse = ", "))

inp <- map(f_in, ~ read_tsv(.x, show_col_types = FALSE, guess_max = 10000))

ev <- read_tsv(f_ev, show_col_types = FALSE, guess_max = 10000)

### units of the main level in palette order
units_main <- intersect(names(pal_unit), unique(filter(inp$gdi, k == k_main)$unit))



# checks -----------------------------------------------------------------------

## check if the BPP tables hold every level
if (!all(k_levels %in% inp$gdi$k) || !k_main %in% inp$params$k) {
  stop("BPP results missing for ", paste(setdiff(k_levels, inp$gdi$k), collapse = ", "))
}

## check if every unit of the main level has pi and is in the recombination tables
for (x in c("pi", "phi", "ia")) {
  d <- if (x == "pi") filter(inp[[x]], k == k_main, set == set_main) else filter(inp[[x]], group == set_main)
  lack <- setdiff(units_main, d$unit)
  if (length(lack) > 0) warning(x, " without ", paste(lack, collapse = ", "))
}

## check if the f-branch table holds the main level
if (!k_main %in% inp$fb$k) stop("f-branch has no results at ", k_main)



# dated unit tree (A) ----------------------------------------------------------

### clades carrying a tau at the main level (BPP node labels, tips sorted)
tau_main <- inp$params %>%
  filter(k == k_main, parameter == "tau") %>%
  mutate(key = map_chr(str_split(members, ","), ~ paste(sort(.x), collapse = ",")))
clades <- unique(tau_main$key)

### newick of a clade from the clades nested in it (the largest ones as children)
newick_of <- function(key) {
  m <- str_split_1(key, ",")
  if (length(m) == 1) return(m)
  sub_k <- clades[map_lgl(clades, ~ { s <- str_split_1(.x, ","); length(s) < length(m) && all(s %in% m) })]
  top_k <- sub_k[!map_lgl(sub_k, function(a) {
    sa <- str_split_1(a, ",")
    any(map_lgl(sub_k, function(b) { sb <- str_split_1(b, ","); length(sb) > length(sa) && all(sa %in% sb) }))
  })]
  kids <- c(top_k, setdiff(m, unlist(str_split(top_k, ","))))
  paste0("(", paste(map_chr(kids, newick_of), collapse = ","), ")")
}

root_key <- clades[which.max(str_count(clades, ",") + 1)]
tr <- ape::read.tree(text = paste0(newick_of(root_key), ";")) %>%
  ape::ladderize() %>%
  ape::reorder.phylo("cladewise")

## check if the tree holds every unit of the main level
if (!setequal(tr$tip.label, units_main)) warning("unit tree and gdi table differ in units")

n_tip <- ape::Ntip(tr)
inner <- n_tip + seq_len(tr$Nnode)
par   <- tr$edge[, 1]
chi   <- tr$edge[, 2]

### node ages at the mean over locus sets; tips in plotting order, inner nodes at the mean of their children
node_key <- c(tr$tip.label, map_chr(inner, ~ paste(sort(ape::extract.clade(tr, .x)$tip.label), collapse = ",")))
node_age <- c(rep(0, n_tip), map_dbl(node_key[inner], ~ mean(tau_main$mean[tau_main$key == .x])))
node_y   <- rep(NA_real_, n_tip + tr$Nnode)
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
hpd_n <- tau_main %>%
  mutate(y = node_y[match(key, node_key)] + (set - 2) * 0.13, set = factor(set))
tips  <- tibble(unit = tr$tip.label, y = node_y[seq_len(n_tip)])
tau_max <- max(hpd_n$hi)

p_tree <- ggplot() +
  geom_segment(data = seg_h, aes(x = x, xend = xend, y = y, yend = y), linewidth = 0.45, colour = "grey25") +
  geom_segment(data = seg_v, aes(x = x, xend = x, y = ymin, yend = ymax), linewidth = 0.45, colour = "grey25") +
  geom_linerange(data = hpd_n, aes(xmin = lo, xmax = hi, y = y), linewidth = 0.3, colour = "grey45") +
  geom_point(data = hpd_n, aes(x = mean, y = y, shape = set), size = 1.2, stroke = 0.3, fill = "grey70",
             colour = "black") +
  geom_point(data = tips, aes(x = 0, y = y, fill = unit), shape = 21, size = 2.4, stroke = 0.3, colour = "black") +
  geom_text(data = tips, aes(x = 0, y = y, label = unit), hjust = 0, nudge_x = 0.04 * tau_max, size = 7 / .pt) +
  scale_x_reverse(limits = c(tau_max * 1.02, -0.16 * tau_max), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  scale_fill_manual(values = pal_unit, guide = "none") +
  scale_shape_manual(values = shape_set, name = "Locus set") +
  labs(x = expression(tau ~ "(substitutions per site)"), y = NULL) +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.key.size = unit(3, "mm"))



# gdi (B) ----------------------------------------------------------------------

gdi <- inp$gdi %>%
  mutate(level = factor(lab_k(k), levels = lab_k(k_levels)),
         unit  = factor(unit, levels = intersect(names(pal_unit), unique(unit))),
         set   = factor(set))

p_gdi <- ggplot(gdi, aes(unit, gdi_mean, group = set)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = gdi_low, ymax = gdi_high, fill = "grey93") +
  geom_hline(yintercept = c(gdi_low, gdi_high), linetype = "dashed", linewidth = 0.25, colour = "grey40") +
  geom_linerange(aes(ymin = gdi_lo, ymax = gdi_hi), linewidth = 0.3, position = position_dodge(0.55)) +
  geom_point(aes(fill = unit, shape = set), size = 1.6, stroke = 0.3, colour = "black",
             position = position_dodge(0.55)) +
  facet_grid(~ level, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal_unit, guide = "none") +
  scale_shape_manual(values = shape_set, guide = "none") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = NULL, y = "gdi") +
  theme(strip.background = element_blank(), panel.spacing = unit(1.5, "mm"))



# theta against pi (C) ---------------------------------------------------------

### theta of the tips: mean over locus sets, range of the 95% HPDs
tp <- inp$gdi %>%
  filter(k == k_main) %>%
  summarise(theta = mean(theta_mean), theta_lo = min(theta_lo), theta_hi = max(theta_hi), .by = unit) %>%
  inner_join(inp$pi %>% filter(k == k_main, set == set_main) %>% select(unit, pi, pi_lo, pi_hi), by = "unit") %>%
  mutate(unit = factor(unit, levels = units_main))

p_theta <- ggplot(tp, aes(pi, theta)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.3, colour = "grey40") +
  geom_linerange(aes(ymin = theta_lo, ymax = theta_hi), linewidth = 0.3, colour = "grey30") +
  geom_linerange(aes(xmin = pi_lo, xmax = pi_hi), linewidth = 0.3, colour = "grey30") +
  geom_point(aes(fill = unit), shape = 21, size = 2.2, stroke = 0.3, colour = "black") +
  geom_text_repel(aes(label = unit), size = 7 / .pt, box.padding = 0.35, min.segment.length = 0.3, seed = 1) +
  scale_fill_manual(values = pal_unit, guide = "none") +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = expression("Nucleotide diversity, " * pi), y = expression(theta ~ "(BPP)"))



# f-branch (D) -----------------------------------------------------------------

### branches with at least one testable donor
fb_main <- inp$fb %>%
  filter(k == k_main, set == set_main) %>%
  filter(any(!is.na(fb)), .by = branch) %>%
  mutate(branch_lab = str_replace_all(branch_descendants, ",", "+"),
         branch_lab = fct_reorder(branch_lab, order, .desc = TRUE),
         P3 = factor(P3, levels = intersect(names(pal_unit), unique(P3))))

### cell text: value, dash where the branch and donor cannot be tested; white on dark cells
fb_main <- fb_main %>%
  mutate(txt = case_when(is.na(fb) ~ "-", fb >= fb_label_min ~ sprintf("%.2f", fb), .default = ""),
         dark = !is.na(fb) & fb > max(fb, na.rm = TRUE) / 2)

p_fb <- ggplot(fb_main, aes(P3, branch_lab, fill = fb)) +
  geom_tile(colour = "grey80", linewidth = 0.3) +
  geom_text(aes(label = txt, colour = dark), size = 5.5 / .pt, show.legend = FALSE) +
  scale_fill_gradient(low = pal_fb[1], high = pal_fb[2], limits = c(0, NA), na.value = "white",
                      name = expression(italic(f)[b])) +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white")) +
  labs(x = "Donor (P3)", y = "Recipient branch") +
  guides(fill = guide_colourbar(barwidth = unit(18, "mm"), barheight = unit(1.5, "mm"), title.vjust = 1)) +
  theme(legend.position = "bottom", axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))



# recombination (E) ------------------------------------------------------------

lab_metric <- c(phi = "PHI genes (%)", ia = "rbarD")

### all MLGs: PHI share with exact binomial interval, rbarD with its permutation p
rec_all <- bind_rows(
  inp$phi %>%
    filter(group == set_main) %>%
    mutate(ci = map2(n_sig, tested, ~ 100 * binom.test(.x, .y)$conf.int)) %>%
    transmute(unit, metric = "phi", value = pct_sig, lo = map_dbl(ci, 1), hi = map_dbl(ci, 2), sig = NA),
  inp$ia %>%
    filter(group == set_main) %>%
    transmute(unit, metric = "ia", value = rbarD, lo = NA_real_, hi = NA_real_, sig = p_rD < alpha_ia))

### equal-size draws: mean and range
rec_equal <- bind_rows(
  inp$phi %>% filter(group == set_equal) %>% transmute(unit, metric = "phi", value = pct_sig),
  inp$ia %>% filter(group == set_equal) %>% transmute(unit, metric = "ia", value = rbarD)) %>%
  summarise(mean = mean(value), lo = min(value), hi = max(value), .by = c(unit, metric))

rec_prep <- function(d) {
  d %>%
    filter(unit %in% units_main) %>%
    mutate(unit = factor(unit, levels = units_main), metric = factor(lab_metric[metric], levels = lab_metric))
}
rec_all   <- rec_prep(rec_all)
rec_equal <- rec_prep(rec_equal)
rec_ref   <- rec_prep(tibble(unit = units_main[1], metric = "phi", y = 100 * alpha_phi))

### significance marks above the bar and the equal-size range
rec_star <- rec_all %>%
  filter(sig %in% TRUE) %>%
  left_join(select(rec_equal, unit, metric, eq_hi = hi), by = c("unit", "metric")) %>%
  mutate(y = pmax(value, eq_hi, na.rm = TRUE))

p_rec <- ggplot(rec_all, aes(unit, value)) +
  geom_hline(data = rec_ref, aes(yintercept = y), linetype = "dashed", linewidth = 0.25, colour = "grey40") +
  geom_col(aes(fill = unit), width = 0.5, colour = "black", linewidth = 0.2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.25, na.rm = TRUE) +
  geom_text(data = rec_star, aes(y = y, label = "*"), vjust = -0.1, size = 9 / .pt) +
  geom_linerange(data = rec_equal, aes(y = mean, ymin = lo, ymax = hi), position = position_nudge(x = 0.36),
                 linewidth = 0.3, colour = "grey30") +
  geom_point(data = rec_equal, aes(y = mean), position = position_nudge(x = 0.36), shape = 23, size = 1.4,
             stroke = 0.3, fill = "white", colour = "black") +
  facet_wrap(~ metric, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_fill_manual(values = pal_unit, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = NULL) +
  theme(strip.placement = "outside", strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))



# lines of evidence (F) --------------------------------------------------------

### one cell per unit of the main level and line of evidence; fill: whether it supports the unit (empty: not
### scored); -: not applicable; mPTP: species within the unit / species shared with other units, scored by the
### shared ones
lab_crit <- c(clusters = "Genomic\nclusters", mono = "Monophyly\n(trees)", stem = "Stem support\n(UFBoot/LPP)",
              mptp = "mPTP species\n(within / shared)", gdi = "gdi\n(locus sets)", a11 = "A11 P\n(distinct species)",
              flow = "Max\nf-branch", mat = "MAT1:MAT2")

### levels at which a unit is found
lv <- ev %>%
  summarise(levels = paste0("K", paste(sort(parse_number(k)), collapse = ",")), .by = unit)

cells <- ev %>%
  filter(k == k_main, unit %in% units_main) %>%
  left_join(lv, by = "unit") %>%
  mutate(n_mono = parse_number(str_extract(mono, "^[0-9]+")), n_tree = parse_number(str_extract(mono, "[0-9]+$")))

cells <- bind_rows(
  transmute(cells, unit, crit = "clusters", txt = levels, call = NA_character_),
  transmute(cells, unit, crit = "mono", txt = coalesce(mono, "-"),
            call = case_when(n_mono == n_tree ~ "Supports", n_mono == 0 ~ "Against", !is.na(n_mono) ~ "Ambiguous")),
  transmute(cells, unit, crit = "stem",
            txt  = if_else(n_mono %in% 0, "-",
                           sprintf("%s/%s", if_else(is.na(ufboot), "-", sprintf("%.0f", ufboot)),
                                   if_else(is.na(lpp), "-", sprintf("%.2f", lpp)))),
            call = case_when(ufboot >= ufboot_min & lpp >= lpp_min ~ "Supports",
                             !is.na(ufboot) & !is.na(lpp) ~ "Ambiguous")),
  transmute(cells, unit, crit = "mptp",
            txt  = if_else(is.na(mptp_species) | is.na(mptp_shared), "-",
                           sprintf("%d / %d", mptp_species, mptp_shared)),
            call = case_when(mptp_shared == 0 ~ "Supports", mptp_shared > 0 ~ "Against")),
  transmute(cells, unit, crit = "gdi", txt = if_else(is.na(gdi_min), "-", sprintf("%.2f-%.2f", gdi_min, gdi_max)),
            call = recode(gdi_call, separate = "Supports", ambiguous = "Ambiguous", `one lineage` = "Against")),
  transmute(cells, unit, crit = "a11", txt = if_else(is.na(a11_pp), "-", sprintf("%.2f", a11_pp)),
            call = case_when(a11_pp >= a11_min ~ "Supports", !is.na(a11_pp) ~ "Ambiguous")),
  transmute(cells, unit, crit = "flow", txt = {
    f <- pmax(fb_into, fb_from, na.rm = TRUE)
    if_else(is.na(f), "-", sprintf("%.2f", f))
  }, call = NA_character_),
  transmute(cells, unit, crit = "mat", txt = if_else(is.na(mat1), "-", sprintf("%d:%d", mat1, mat2)),
            call = NA_character_)) %>%
  mutate(crit = factor(crit, levels = names(lab_crit), labels = lab_crit),
         unit = factor(unit, levels = rev(units_main)),
         call = factor(call, levels = names(pal_call)))

p_ev <- ggplot(cells, aes(crit, unit)) +
  geom_tile(aes(fill = call), colour = "grey70", linewidth = 0.3) +
  geom_text(aes(label = txt, colour = call %in% "Supports"), size = 6.5 / .pt, show.legend = FALSE) +
  geom_tile(data = distinct(cells, unit), aes(x = 0.35, y = unit, fill = NULL), width = 0.22,
            fill = pal_unit[levels(cells$unit)][as.integer(distinct(cells, unit)$unit)], colour = "black",
            linewidth = 0.2, inherit.aes = FALSE) +
  scale_fill_manual(values = pal_call, na.value = "white", na.translate = FALSE, drop = FALSE, name = NULL) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_x_discrete(position = "top", expand = expansion(add = c(0.75, 0.05))) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "bottom", legend.key.size = unit(3, "mm"), axis.line = element_blank(),
        axis.ticks = element_blank(), axis.text.x = element_text(size = 6.5, lineheight = 0.9))



# figures ----------------------------------------------------------------------

row_1 <- (p_tree | p_gdi) + plot_layout(widths = c(0.85, 1.25))
row_2 <- (p_theta | p_fb | p_rec) + plot_layout(widths = c(0.85, 1.15, 0.8))

p_fig <- row_1 / row_2 / p_ev + plot_layout(heights = c(1, 1.1, 0.7))
p_fig <- p_fig +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold"))

save_fig(sprintf("publication_delimitation_%s.pdf", k_main), p_fig, 7.1, 8.4)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat(sprintf("  publication_delimitation_%s.pdf\n", k_main))
