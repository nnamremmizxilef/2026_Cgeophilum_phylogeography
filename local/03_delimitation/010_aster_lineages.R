## ============================================================================
## 010_aster_lineages.R
## Unit trees at K = 4, 5 and 6 from ASTER with individuals mapped to units:
## quartet support and local posterior per branch, polytomy test, sister units
## and terminal branch lengths in coalescent units
##
## Files needed:
##   aster_cluster_tree_kX.tre           unit tree with q1-q3, f1-f3 and local PP per branch
##   aster_cluster_tree_kX_polytomy.tre  the same topology scored by ASTRAL-III: polytomy p per
##                                       branch, lengths in CU
##   individuals.tsv                     tip and its unit at every level
##                                       (all from hpc/03_delimitation/010_aster_lineages.sh, kX = k4, k5, k6;
##                                        copy them from ${ROOT}/03_delimitation/results/010_results_aster_lineages
##                                        to local/03_delimitation/data/010_data_aster_lineages)
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)
library(ape)



# files ------------------------------------------------------------------------

data_root <- "03_delimitation/data/010_data_aster_lineages"
out_root  <- "03_delimitation/results/010_results_aster_lineages"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### levels (as in hpc/03_delimitation/010_aster_lineages.sh)
k_levels <- c("k4", "k5", "k6")

### outgroup unit name in the mapping
outgroup <- "outgroup"

### polytomy test significance
alpha <- 0.05



# style ------------------------------------------------------------------------

### unit colours
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

### outgroup fill and label
fill_outgroup <- "grey60"
lab_outgroup  <- "Outgroup"

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# functions --------------------------------------------------------------------

### units below every internal node (nodes in tree order)
clade_units <- function(tr) map(prop.part(tr), ~ sort(tr$tip.label[.x]))

### a branch named by the units on its ingroup side (NA for terminal branches)
branch_key <- function(below, all_units) {
  side <- if (outgroup %in% below) setdiff(all_units, below) else below
  if (length(side) < 2 || length(all_units) - length(side) < 2) return(NA_character_)
  paste(sort(side), collapse = "+")
}

### length of the edge above every node (NA for the root)
edge_above <- function(tr) tr$edge.length[match(seq_len(Ntip(tr) + tr$Nnode), tr$edge[, 2])]

### ASTER tree: annotations moved out of the quoted labels, one row per internal branch
read_aster <- function(file) {
  txt <- str_trim(read_file(file))
  m   <- gregexpr("'\\[[^]]*\\]'", txt)
  ann <- regmatches(txt, m)[[1]]
  regmatches(txt, m) <- list(paste0("n", seq_along(ann)))
  tr  <- read.tree(text = txt)

  vals <- tibble(node_label = paste0("n", seq_along(ann)),
                 a          = str_remove_all(ann, "^'\\[|\\]'$")) %>%
    separate_longer_delim(a, ";") %>%
    separate_wider_delim(a, "=", names = c("name", "value")) %>%
    mutate(value = as.numeric(value)) %>%
    pivot_wider(names_from = name, values_from = value)

  units <- sort(tr$tip.label)
  branches <- tibble(node_label = tr$node.label,
                     branch     = map_chr(clade_units(tr), branch_key, all_units = units)) %>%
    filter(!is.na(branch)) %>%
    left_join(vals, by = "node_label") %>%
    filter(!is.na(q1)) %>%
    distinct(branch, .keep_all = TRUE) %>%
    select(branch, q1, q2, q3, f1, f2, f3, localPP, pp1, pp2, pp3)

  list(tree = tr, branches = branches)
}

### ASTRAL-III tree: polytomy p and length (CU) per internal branch, terminal
### length per unit
read_astral <- function(file) {
  tr    <- read.tree(file)
  n     <- Ntip(tr)
  units <- sort(tr$tip.label)
  len   <- edge_above(tr)
  root  <- n + 1

  branches <- tibble(branch  = map_chr(clade_units(tr), branch_key, all_units = units),
                     p_poly  = as.numeric(na_if(tr$node.label, "")),
                     cu      = len[root:(n + tr$Nnode)]) %>%
    filter(!is.na(branch)) %>%
    summarise(p_poly = first(na.omit(p_poly)), cu = sum(cu, na.rm = TRUE), .by = branch)

  ### a tip next to a two-way root takes the length of the other root edge as well
  root_kids <- tr$edge[tr$edge[, 1] == root, 2]
  terminal  <- len[seq_len(n)]
  if (length(root_kids) == 2) {
    tip_kid <- root_kids[root_kids <= n]
    if (length(tip_kid) == 1) terminal[tip_kid] <- terminal[tip_kid] +
      replace_na(len[setdiff(root_kids, tip_kid)], 0)
  }

  list(tree = tr, branches = branches,
       terminal = tibble(unit = tr$tip.label, cu_terminal = terminal))
}

### rectangular cladogram (tips aligned, nodes one step left of their deepest
### descendant)
layout_tree <- function(tr) {
  n  <- Ntip(tr)
  ed <- tr$edge
  nn <- n + tr$Nnode
  y  <- rep(NA_real_, nn)
  h  <- rep(0, nn)
  lo <- rep(Inf, nn)
  hi <- rep(-Inf, nn)

  tips_order <- ed[ed[, 2] <= n, 2]
  y[tips_order] <- seq_along(tips_order)

  for (e in rev(seq_len(nrow(ed)))) {
    p <- ed[e, 1]
    c <- ed[e, 2]
    if (c > n) y[c] <- (lo[c] + hi[c]) / 2
    h[p]  <- max(h[p], h[c] + 1)
    lo[p] <- min(lo[p], y[c])
    hi[p] <- max(hi[p], y[c])
  }
  y[n + 1] <- (lo[n + 1] + hi[n + 1]) / 2

  tibble(node = seq_len(nn), x = max(h) - h, y = y, lo = lo, hi = hi,
         unit = c(tr$tip.label, rep(NA_character_, tr$Nnode)))
}



# read -------------------------------------------------------------------------

ind <- read_tsv(file.path(data_root, "individuals.tsv"), col_types = cols(.default = "c"))

aster  <- map(set_names(k_levels),
              ~ read_aster(file.path(data_root, sprintf("aster_cluster_tree_%s.tre", .x))))
astral <- map(set_names(k_levels),
              ~ read_astral(file.path(data_root, sprintf("aster_cluster_tree_%s_polytomy.tre", .x))))



# branches ---------------------------------------------------------------------

### one row per level and internal branch (support from ASTER, polytomy p and
### length from ASTRAL-III; rounding noise below zero set to zero)
branches <- map_dfr(k_levels, ~ left_join(aster[[.x]]$branches, astral[[.x]]$branches,
                                          by = "branch") %>%
                      mutate(k = .x, .before = 1)) %>%
  mutate(across(c(q1:q3, f1:f3, pp1:pp3), ~ pmax(.x, 0)),
         EN       = f1 + f2 + f3,
         polytomy = p_poly >= alpha)



# units ------------------------------------------------------------------------

### rooted unit trees (topology only)
rooted <- map(aster, function(a) {
  tr <- root(a$tree, outgroup = outgroup, resolve.root = TRUE)
  tr$node.label  <- NULL
  tr$edge.length <- NULL
  tr
})

### same tip order at every level (finest level ladderized, a coarser unit at
### the place of its parts)
k_fine     <- last(k_levels)
tr_fine    <- ladderize(rooted[[k_fine]])
order_fine <- tr_fine$tip.label[tr_fine$edge[tr_fine$edge[, 2] <= Ntip(tr_fine), 2]]

rooted <- map(set_names(k_levels), function(kk) {
  coarse <- deframe(distinct(ind, .data[[k_fine]], .data[[kk]]))
  reorder(rotateConstr(rooted[[kk]], unique(coarse[order_fine])), "cladewise")
})

### sister of every unit and support of the clade with its sister
sisters <- map_dfr(k_levels, function(kk) {
  tr    <- rooted[[kk]]
  n     <- Ntip(tr)
  below <- clade_units(tr)
  map_dfr(setdiff(tr$tip.label, outgroup), function(u) {
    parent <- tr$edge[tr$edge[, 2] == match(u, tr$tip.label), 1]
    tibble(k      = kk,
           unit   = u,
           sister = paste(setdiff(below[[parent - n]], u), collapse = "+"),
           branch = branch_key(below[[parent - n]], sort(tr$tip.label)))
  })
})

n_unit <- ind %>%
  pivot_longer(all_of(k_levels), names_to = "k", values_to = "unit") %>%
  count(k, unit)

units <- sisters %>%
  left_join(n_unit, by = c("k", "unit")) %>%
  left_join(map_dfr(k_levels, ~ mutate(astral[[.x]]$terminal, k = .x)), by = c("k", "unit")) %>%
  left_join(select(branches, k, branch, q1_sister = q1, p_sister = p_poly),
            by = c("k", "branch")) %>%
  select(k, unit, n, sister, q1_sister, p_sister, cu_terminal) %>%
  arrange(k, unit)



# checks -----------------------------------------------------------------------

## check if ASTER and ASTRAL-III hold the same branches at every level
for (kk in k_levels) {
  if (!setequal(aster[[kk]]$branches$branch, astral[[kk]]$branches$branch)) {
    warning("ASTER and ASTRAL-III branches differ at ", kk)
  }
}

## check if every tree holds the units of its level and the outgroup
for (kk in k_levels) {
  if (!setequal(aster[[kk]]$tree$tip.label, unique(ind[[kk]]))) {
    warning("tree tips and units of individuals.tsv differ at ", kk)
  }
}

## check if every tree is fully resolved
n_res <- branches %>%
  count(k) %>%
  left_join(count(distinct(n_unit, k, unit), k, name = "units"), by = "k")

if (any(n_res$n != n_res$units - 3)) {
  warning("a unit tree is not fully resolved")
}

## check if the quartet supports of every branch sum to one
if (any(abs(branches$q1 + branches$q2 + branches$q3 - 1) > 1e-3)) {
  warning("quartet supports that do not sum to one")
}

## check if every branch has a polytomy p and every unit a terminal length
if (anyNA(branches$p_poly) || anyNA(units$cu_terminal)) {
  warning("branches without polytomy p or units without terminal length")
}



# figures ----------------------------------------------------------------------

### one rooted cladogram per level (q1 on internal branches, dashed where not
### distinguishable from a polytomy)
plot_tree <- function(kk) {
  tr   <- rooted[[kk]]
  n    <- Ntip(tr)
  lay  <- layout_tree(tr)
  keys <- c(rep(NA_character_, n),
            map_chr(clade_units(tr), branch_key, all_units = sort(tr$tip.label)))

  edges <- tibble(parent = tr$edge[, 1], child = tr$edge[, 2], branch = keys[tr$edge[, 2]]) %>%
    left_join(select(lay, parent = node, x0 = x), by = "parent") %>%
    left_join(select(lay, child = node, x1 = x, y1 = y), by = "child") %>%
    left_join(filter(branches, k == kk) %>% select(branch, q1, polytomy), by = "branch") %>%
    mutate(polytomy = replace_na(polytomy, FALSE))

  tips <- lay %>%
    filter(node <= n) %>%
    mutate(label = if_else(unit == outgroup, lab_outgroup, unit))

  root_xy <- filter(lay, node == n + 1)

  ggplot() +
    geom_segment(data = filter(lay, node > n), aes(x = x, xend = x, y = lo, yend = hi),
                 linewidth = 0.4) +
    geom_segment(data = edges, aes(x = x0, xend = x1, y = y1, yend = y1, linetype = polytomy),
                 linewidth = 0.4) +
    geom_segment(data = root_xy, aes(x = x - 0.3, xend = x, y = y, yend = y), linewidth = 0.4) +
    geom_text(data = filter(edges, !is.na(q1)),
              aes(x = (x0 + x1) / 2, y = y1, label = sprintf("%.2f", q1)),
              vjust = -0.5, size = 7 / .pt) +
    geom_point(data = tips, aes(x = x, y = y, fill = unit), shape = 21, colour = "black",
               size = 2.2, stroke = 0.3) +
    geom_text(data = tips, aes(x = x + 0.2, y = y, label = label), hjust = 0, size = 8 / .pt) +
    scale_fill_manual(values = c(pal_unit, setNames(fill_outgroup, outgroup)), guide = "none") +
    scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.35))) +
    scale_y_continuous(expand = expansion(add = 0.6)) +
    labs(x = sprintf("Units at K = %s", sub("k", "", kk))) +
    theme_void(base_size = 10) +
    theme(axis.title.x = element_text(margin = margin(t = 4)))
}

### panels tagged (tag style set here, theme_void replaces the global theme)
p_tree <- wrap_plots(map(k_levels, plot_tree), nrow = 1) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold"))

save_fig("aster_unit_trees.pdf", p_tree, 7.1, 2.4)



# tables -----------------------------------------------------------------------

fmt_p <- function(p) if_else(p < 0.001, "< 0.001", sprintf("%.3f", p))

### one row per level and internal branch
pub_branches <- branches %>%
  transmute(Level                         = sub("k", "K = ", k),
            Branch                        = branch,
            `Quartet support q1 (q2, q3)` = sprintf("%.2f (%.2f, %.2f)", q1, q2, q3),
            `Local PP`                    = sprintf("%.2f", localPP),
            `Length (CU)`                 = sprintf("%.2f", cu),
            `Polytomy p`                  = fmt_p(p_poly))

### one row per level and unit
pub_units <- units %>%
  transmute(Level                    = sub("k", "K = ", k),
            Unit                     = unit,
            n                        = n,
            Sister                   = sister,
            `q1 of unit and sister`  = if_else(is.na(q1_sister), "", sprintf("%.2f", q1_sister)),
            `Polytomy p`             = if_else(is.na(p_sister), "", fmt_p(p_sister)),
            `Terminal branch (CU)`   = sprintf("%.2f", cu_terminal))

write_tsv(branches, file.path(out_root, "aster_branches.tsv"))
write_tsv(units, file.path(out_root, "aster_units.tsv"))
write_tsv(pub_branches, file.path(out_root, "publication_aster_branches.tsv"))
write_tsv(pub_units, file.path(out_root, "publication_aster_units.tsv"))

cat("\n=== Unit trees (rooted) ===\n")
walk(k_levels, ~ cat(.x, write.tree(rooted[[.x]]), "\n"))

cat("\n=== Internal branches ===\n")
branches %>%
  mutate(across(c(q1, q2, q3, localPP, cu), ~ round(.x, 3)), p_poly = signif(p_poly, 3)) %>%
  select(k, branch, q1, q2, q3, localPP, cu, p_poly) %>%
  print(n = Inf)

cat("\n=== Units ===\n")
units %>% mutate(across(where(is.double), ~ signif(.x, 3))) %>% print(n = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  aster_unit_trees.pdf\n")
cat("  aster_branches.tsv\n")
cat("  aster_units.tsv\n")
cat("  publication_aster_branches.tsv\n")
cat("  publication_aster_units.tsv\n")
