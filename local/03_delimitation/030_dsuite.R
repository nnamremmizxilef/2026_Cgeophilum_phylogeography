## ============================================================================
## 030_dsuite.R
## D statistics, f4-ratios and f-branch between units at K = 4, 5 and 6 (all
## representatives and without mixed-MAT samples) and the most frequent site pattern per trio
##
## Files needed:
##   dtrios/kX_<set>_tree.txt     D and f4-ratio per trio arranged by the unit tree
##   dtrios/kX_<set>_fbranch.txt  f-branch matrix
##                                (both from hpc/03_delimitation/030_dsuite.sh;
##                                 copy dtrios/k*_<set>_tree.txt and dtrios/k*_<set>_fbranch.txt
##                                 from ${ROOT}/03_delimitation/results/030_results_dsuite
##                                 to local/03_delimitation/data/030_data_dsuite/dtrios)
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

data_root <- "03_delimitation/data/030_data_dsuite"
out_root  <- "03_delimitation/results/030_results_dsuite"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### levels and sample sets (as in hpc/03_delimitation/030_dsuite.sh)
k_levels <- c("k4", "k5", "k6")
sets     <- c("all", "no_mixed")

### significance (BH-adjusted within level and sample set)
alpha <- 0.05



# style ------------------------------------------------------------------------

### f-branch colour scale
pal_fb <- c("white", "#B2182B")

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

### D and f4-ratio per trio (arranged by the unit tree)
trios <- expand_grid(k = k_levels, set = sets) %>%
  mutate(data = map2(k, set, ~ read_tsv(file.path(data_root, "dtrios", sprintf("%s_%s_tree.txt", .x, .y)),
                                        col_types = cols(P1 = "c", P2 = "c", P3 = "c", .default = "d")))) %>%
  unnest(data) %>%
  rename(D = Dstatistic, Z = `Z-score`, p = `p-value`, f4 = `f4-ratio`) %>%
  mutate(q = p.adjust(p, "BH"), .by = c(k, set)) %>%
  mutate(sig = q < alpha)

### most frequent site pattern per trio (BBAA: P1 and P2 share the derived allele, as in the unit tree)
trios <- trios %>%
  mutate(bbaa_share = BBAA / (BBAA + ABBA + BABA),
         top        = c("BBAA", "ABBA", "BABA")[max.col(cbind(BBAA, ABBA, BABA), ties.method = "first")],
         tree_top   = top == "BBAA")

### f-branch matrices (branch, its descendants, one column per P3 unit)
fbranch <- expand_grid(k = k_levels, set = sets) %>%
  mutate(data = map2(k, set, ~ read_tsv(file.path(data_root, "dtrios", sprintf("%s_%s_fbranch.txt", .x, .y)),
                                        na = c("nan", "NaN", "NA"), col_types = cols(.default = "c")))) %>%
  mutate(data = map(data, ~ .x %>%
                      mutate(order = row_number()) %>%
                      pivot_longer(-c(branch, branch_descendants, order), names_to = "P3", values_to = "fb") %>%
                      mutate(fb = as.numeric(fb)))) %>%
  unnest(data) %>%
  filter(P3 != "Outgroup")



# checks -----------------------------------------------------------------------

## check if every level has all trios of its units in both sample sets
n_exp <- trios %>%
  summarise(n = n(), units = n_distinct(c(P1, P2, P3)), .by = c(k, set)) %>%
  mutate(expected = units * (units - 1) * (units - 2) / 6)

if (any(n_exp$n != n_exp$expected)) {
  warning("a level or sample set lacks trios")
}

## check if both sample sets hold the same trios
n_match <- trios %>% count(k, P1, P2, P3) %>% filter(n != length(sets))

if (nrow(n_match) > 0) {
  warning(nrow(n_match), " trios arranged differently or missing between the sample sets")
}

## check if every trio has a D statistic
if (anyNA(trios$D)) {
  warning(sum(is.na(trios$D)), " trios without a D statistic")
}

## check if the unit tree arrangement is the most frequent site pattern in every trio
n_conflict <- trios %>% filter(set == "all", !tree_top)

if (nrow(n_conflict) > 0) {
  message(nrow(n_conflict),
          " trios where a pattern other than BBAA is the most frequent (the unit tree does not fit these SNPs)")
}



# figures ----------------------------------------------------------------------

### largest f-branch over the levels (one colour scale)
fb_max <- max(fbranch$fb[fbranch$set == "all"], na.rm = TRUE)

### one f-branch panel per level (rows: branches from the tips to the root; empty: not testable)
plot_fbranch <- function(kk) {
  d <- fbranch %>%
    filter(k == kk, set == "all") %>%
    mutate(branch_lab = str_replace_all(branch_descendants, ",", "+"),
           branch_lab = fct_reorder(branch_lab, order, .desc = TRUE),
           P3         = factor(P3, levels = unique(P3)))

  ggplot(d, aes(x = P3, y = branch_lab, fill = fb)) +
    geom_tile(colour = "grey80", linewidth = 0.2) +
    scale_fill_gradient(low = pal_fb[1], high = pal_fb[2], limits = c(0, fb_max),
                        na.value = "grey95", name = expression(f[b])) +
    coord_equal() +
    labs(x = sprintf("P3 (K = %s)", sub("k", "", kk)), y = "Branch") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

p_fb <- wrap_plots(map(k_levels, plot_fbranch), nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

save_fig("fbranch.pdf", p_fb, 7.1, 3.2)



# tables -----------------------------------------------------------------------

### one row per level and trio (all representatives, with the run without mixed-MAT samples)
pub <- trios %>%
  filter(set == "all") %>%
  left_join(trios %>% filter(set == "no_mixed") %>% select(k, P1, P2, P3, D_nm = D, q_nm = q),
            by = c("k", "P1", "P2", "P3")) %>%
  transmute(Level                       = sub("k", "K = ", k),
            Trio                        = sprintf("(%s, %s), %s", P1, P2, P3),
            D                           = sprintf("%.3f", D),
            Z                           = sprintf("%.2f", Z),
            q                           = format.pval(q, digits = 2, eps = 1e-10),
            `f4-ratio`                  = sprintf("%.3f", f4),
            `BBAA, ABBA, BABA`          = sprintf("%.0f, %.0f, %.0f", BBAA, ABBA, BABA),
            `Most frequent pattern`     = top,
            `Without mixed MAT: D (q)`  = sprintf("%.3f (%s)", D_nm, format.pval(q_nm, digits = 2, eps = 1e-10)))

write_tsv(trios, file.path(out_root, "dsuite_trios.tsv"))
write_tsv(fbranch, file.path(out_root, "dsuite_fbranch.tsv"))
write_tsv(pub, file.path(out_root, "publication_dsuite.tsv"))



# summary ----------------------------------------------------------------------

cat("\n=== D per trio, all representatives (q < 0.05 marked) ===\n")
trios %>%
  filter(set == "all") %>%
  mutate(across(c(D, Z, f4, bbaa_share), ~ round(.x, 3)), q = signif(q, 2), sig = if_else(sig, "*", "")) %>%
  select(k, P1, P2, P3, D, Z, q, f4, sig, BBAA, ABBA, BABA, bbaa_share, top) %>%
  print(n = Inf, width = Inf)

cat("\n=== Trios with q < 0.05 per level and sample set ===\n")
trios %>% summarise(significant = sum(sig), trios = n(), .by = c(k, set)) %>% print()



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  fbranch.pdf\n")
cat("  dsuite_trios.tsv\n")
cat("  dsuite_fbranch.tsv\n")
cat("  publication_dsuite.tsv\n")
