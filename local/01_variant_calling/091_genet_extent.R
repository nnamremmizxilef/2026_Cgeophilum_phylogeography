## ============================================================================
## 091_genet_extent.R
## Extent of the clonal genotypes across the soil pits and depth intervals
##
## Files needed:
##   per_sample_mlg.tsv             multilocus genotype per sample,
##                                  from local/01_variant_calling/090_clone_id.R
##   plink.genome                   pairwise IBS counts, from hpc/01_variant_calling/090_clone_id.sh
##                                  (copy ${ROOT}/01_variant_calling/results/090_results_clone_id/plink.genome
##                                   to local/01_variant_calling/data/090_data_clone_id)
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

supp     <- "01_variant_calling/results/090_results_clone_id"
data_090 <- "01_variant_calling/data/090_data_clone_id"
out_root <- "01_variant_calling/results/091_results_genet_extent"

f_mlg    <- file.path(supp, "per_sample_mlg.tsv")
f_genome <- file.path(data_090, "plink.genome")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### sampling design (four pits 5 m from a centre, one per cardinal direction;
### depth intervals at 5, 25 and 45 cm)
radius     <- 5
directions <- c("N", "E", "S", "W")

design <- tibble(direction = directions,
                 px = c(0, radius, 0, -radius),
                 py = c(radius, 0, -radius, 0))

depths <- tibble(depth = c("1", "2", "3"),
                 dz    = c(0.05, 0.25, 0.45))

### projection
stretch <- sqrt(3)   ### x against y (view angle)
squash  <- 0.45      ### flatness of the view
z_exag  <- 28        ### depth exaggeration
blow_up <- 3         ### pit enlargement



# style ------------------------------------------------------------------------

col_none <- "grey88"

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# projection -------------------------------------------------------------------

pit_r   <- 0.5 * blow_up / sqrt(2)   ### half diagonal of a drawn pit

sx_of <- function(px, py)     px * stretch
sy_of <- function(px, py, dz) py * squash - dz * z_exag

### horizontal square per row of d, given its half diagonal
square <- function(d, r) {
  d %>%
    mutate(id = row_number()) %>%
    expand_grid(tibble(ox = c(0, 1, 0, -1), oy = c(1, 0, -1, 0))) %>%
    mutate(sx = sx_of(px + ox * r, py + oy * r),
           sy = sy_of(px + ox * r, py + oy * r, dz))
}

### one plot square per depth, deepest first
dz0 <- min(depths$dz)

planes <- square(tibble(px = 0, py = 0, dz = rev(depths$dz)), radius)

### pit blocks, drawn from the bottom of the page upwards
blocks <- square(arrange(expand_grid(design, depths), sy_of(px, py, dz)), pit_r)



# read -------------------------------------------------------------------------

### samples placed in the design (COUNTRY-SITE-<direction><depth><replicate>;
### 1-xxx samples carry no position)
mlg <- read_tsv(f_mlg, show_col_types = FALSE) %>%
  filter(str_detect(sample, "^[A-Z]{2}-[A-Z]{3}-[NESW]\\d[A-Z]$")) %>%
  mutate(direction = str_sub(sample, 8, 8),
         depth     = str_sub(sample, 9, 9)) %>%
  left_join(design, by = "direction") %>%
  left_join(depths, by = "depth")

message("samples placed: ", nrow(mlg),
        "   sites: ", n_distinct(mlg$site),
        "   genotypes with more than one sample: ",
        n_distinct(mlg$mlg[mlg$n_members > 1]))



# clonal pairs by separation ---------------------------------------------------

### same-site sample pairs by pit separation (0, 7.07, 10 m) and depth
### separation (0, 0.20, 0.40 m)
pair <- mlg %>%
  select(site, sample, mlg, px, py, dz) %>%
  inner_join(., ., by = "site", relationship = "many-to-many",
             suffix = c("_1", "_2")) %>%
  filter(sample_1 < sample_2) %>%
  transmute(site, sample_1, sample_2,
            d_h    = round(sqrt((px_2 - px_1)^2 + (py_2 - py_1)^2), 2),
            d_v    = round(abs(dz_2 - dz_1), 2),
            clonal = mlg_1 == mlg_2)

sep <- pair %>%
  group_by(d_h, d_v) %>%
  summarise(n_pairs  = n(),
            n_clonal = sum(clonal),
            p_clonal = n_clonal / n_pairs,
            .groups = "drop")

write_tsv(sep, file.path(out_root, "clonal_pairs_by_separation.tsv"))

cat("\n=== pairs of samples from one site ===\n")
print(sep, n = Inf, width = Inf)



# decay of clonality along the two axes ----------------------------------------

### logistic fit of clonality on both separations (anisotropy: ratio of the
### depth to the horizontal coefficient)
aniso_of <- function(d) {
  b <- coef(glm(clonal ~ d_h + d_v, binomial, data = d))
  b[["d_v"]] / b[["d_h"]]
}

fit     <- glm(clonal ~ d_h + d_v, binomial, data = pair)
by_site <- split(pair, pair$site)

### 95% interval from resampling whole sites
set.seed(1)
boot <- replicate(499, aniso_of(bind_rows(
  sample(by_site, length(by_site), replace = TRUE))))

cat("\n=== decay of clonality along the two axes ===\n")
print(round(summary(fit)$coefficients, 4))
cat(sprintf("\nanisotropy (depth / horizontal coefficient)  %.0f\n",
            aniso_of(pair)))
cat(sprintf("  site bootstrap, 95%%                         %.0f to %.0f\n",
            quantile(boot, 0.025, na.rm = TRUE),
            quantile(boot, 0.975, na.rm = TRUE)))



# missing data and the decay with depth ----------------------------------------

### missing genotypes by depth interval
cat("\n=== missing genotypes by depth interval ===\n")
mlg %>%
  group_by(depth, dz) %>%
  summarise(n_samples = n(),
            f_missing = mean(f_missing),
            .groups = "drop") %>%
  print(width = Inf)

### sites compared by depth separation
cat("\n=== sites compared by depth separation ===\n")
read_table(f_genome, show_col_types = FALSE) %>%
  transmute(sample_1 = IID1, sample_2 = IID2,
            n_compared = IBS0 + IBS1 + IBS2) %>%
  inner_join(select(mlg, sample_1 = sample, site_1 = site, dz_1 = dz),
             by = "sample_1") %>%
  inner_join(select(mlg, sample_2 = sample, site_2 = site, dz_2 = dz),
             by = "sample_2") %>%
  filter(site_1 == site_2) %>%
  group_by(d_v = round(abs(dz_2 - dz_1), 2)) %>%
  summarise(n_pairs = n(),
            sites_min    = min(n_compared),
            sites_mean   = mean(n_compared),
            .groups = "drop") %>%
  print(width = Inf)

### fit on the samples at or below the median missingness
cut_miss   <- median(mlg$f_missing)
keep       <- mlg$sample[mlg$f_missing <= cut_miss]
pair_clean <- filter(pair, sample_1 %in% keep, sample_2 %in% keep)

cat(sprintf("\nsamples at or below the median missingness (%.4f)  %d of %d\n",
            cut_miss, length(keep), nrow(mlg)))
print(round(summary(glm(clonal ~ d_h + d_v, binomial,
                        data = pair_clean))$coefficients, 4))
cat(sprintf("\nanisotropy, these samples  %.0f (all samples %.0f)\n",
            aniso_of(pair_clean), aniso_of(pair)))



# genotype extent --------------------------------------------------------------

### horizontal span (largest distance between occupied pits) and vertical span
### (shallowest to deepest interval) per genotype
span <- mlg %>%
  filter(n_members > 1) %>%
  group_by(site, mlg, n_members) %>%
  summarise(pits     = paste(sort(unique(direction)), collapse = ""),
            n_pits   = n_distinct(direction),
            n_depths = n_distinct(depth),
            span_v   = max(dz) - min(dz),
            .groups = "drop") %>%
  mutate(span_h = case_when(
    n_pits == 1 ~ 0,
    n_pits == 2 & pits %in% c("NS", "EW") ~ 2 * radius,
    n_pits == 2 ~ radius * sqrt(2),
    TRUE ~ 2 * radius))

write_tsv(span, file.path(out_root, "genet_extent.tsv"))

cat("\n=== genotypes found more than once ===\n")
cat(sprintf("genotypes                 %d\n", nrow(span)))
cat(sprintf("one pit, one depth        %d\n",
            sum(span$n_pits == 1 & span$n_depths == 1)))
cat(sprintf("one pit, several depths   %d\n",
            sum(span$n_pits == 1 & span$n_depths > 1)))
cat(sprintf("several pits, one depth   %d\n",
            sum(span$n_pits > 1 & span$n_depths == 1)))
cat(sprintf("several pits and depths   %d\n",
            sum(span$n_pits > 1 & span$n_depths > 1)))
cat(sprintf("widest genotype           %.1f m across, %.2f m deep\n",
            max(span$span_h), max(span$span_v)))

cat("\n=== the genotypes reaching beyond one pit ===\n")
span %>%
  filter(n_pits > 1) %>%
  arrange(desc(n_pits), desc(n_members)) %>%
  select(site, mlg, n_members, pits, n_depths, span_h, span_v) %>%
  print(n = Inf)



# figures ----------------------------------------------------------------------

### sampling design (depth planes behind the pit blocks, 5 m arm to the east
### pit, directions in the deepest blocks)
sites   <- unique(mlg$site)
columns <- design %>%
  mutate(sx = sx_of(px, py), sy0 = sy_of(px, py, dz0),
         sy1 = sy_of(px, py, max(depths$dz)))

centres <- design %>%
  mutate(sx = sx_of(px, py), sy = sy_of(px, py, dz0))

east <- filter(centres, direction == "E")

dir_lab <- design %>%
  mutate(sx = sx_of(px, py), sy = sy_of(px, py, max(depths$dz)))

p_design <- ggplot(mapping = aes(x = sx, y = sy)) +
  geom_segment(data = columns, aes(xend = sx, y = sy0, yend = sy1),
               colour = "grey70", linewidth = 0.3) +
  geom_polygon(data = planes, aes(group = id), alpha = 0.55,
               fill = "grey93", colour = "grey75", linewidth = 0.3) +
  geom_polygon(data = blocks, aes(group = id),
               fill = col_none, colour = "grey45", linewidth = 0.3) +
  geom_segment(data = east, aes(x = 0, xend = sx, yend = sy),
               colour = "grey45", linewidth = 0.35) +
  geom_point(data = bind_rows(east, mutate(east, sx = 0)),
             shape = 3, size = 2.2, stroke = 0.6, colour = "grey30") +
  geom_text(data = east, aes(x = sx / 2), label = "5 m",
            vjust = -0.7, size = 2.8, colour = "black") +
  geom_text(data = filter(blocks, direction == "E", ox == 1),
            aes(label = sprintf("%.0f cm", dz * 100)),
            hjust = -0.3, size = 2.8, colour = "black") +
  geom_text(data = dir_lab, aes(label = direction), size = 2.8,
            colour = "black") +
  coord_fixed(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(), plot.margin = margin(6, 26, 6, 6))

save_fig("sampling_design.pdf", p_design, 3.4, 3.8)

### genotypes in the design (one point per genotype, pit and depth, sized by
### samples; paths join the positions of a genotype; golden-angle hues per site)
gcol_of <- span %>%
  arrange(site, desc(n_members), mlg) %>%
  transmute(site, mlg, gcol = as.character(mlg))

col_genotype <- setNames(
  hcl(h = (seq_len(nrow(gcol_of)) * 137.508) %% 360, c = 70, l = 58),
  gcol_of$gcol)

pos <- mlg %>%
  count(site, direction, depth, mlg, px, py, dz, name = "n_samples") %>%
  left_join(gcol_of, by = c("site", "mlg")) %>%
  mutate(gcol = replace_na(gcol, "single")) %>%
  arrange(site, direction, depth, gcol == "single", mlg) %>%
  group_by(site, direction, depth) %>%
  mutate(sx = sx_of(px, py) + 0.55 * (row_number() - (n() + 1) / 2),
         sy = sy_of(px, py, dz)) %>%
  ungroup()

### genotypes in one position only are not joined
paths <- pos %>%
  filter(gcol != "single") %>%
  group_by(site, mlg) %>%
  filter(n() > 1) %>%
  ungroup()

### size key from one sample to the largest position
n_max <- max(pos$n_samples)
n_brk <- unique(round(seq(1, n_max, length.out = min(4, n_max))))

per_site <- function(d) {
  expand_grid(site = sites, d) %>% mutate(id = paste(site, id))
}

p_map <- ggplot(mapping = aes(x = sx, y = sy)) +
  geom_polygon(data = per_site(planes), aes(group = id),
               fill = "grey95", colour = "grey80", linewidth = 0.25) +
  geom_polygon(data = per_site(blocks), aes(group = id),
               fill = "grey92", colour = "grey70", linewidth = 0.2) +
  geom_line(data = paths, aes(group = mlg, colour = gcol), linewidth = 0.4) +
  geom_point(data = pos, aes(fill = gcol, size = n_samples),
             shape = 21, colour = "black", stroke = 0.25) +
  scale_colour_manual(values = col_genotype, guide = "none") +
  scale_fill_manual(values = c(col_genotype, single = "grey80"),
                    breaks = "single", labels = "genotype found once",
                    name = NULL) +
  scale_size_area(max_size = 4, breaks = n_brk, limits = c(0, n_max),
                  name = "Samples") +
  guides(fill = guide_legend(order = 1, override.aes = list(size = 3)),
         size = guide_legend(order = 2,
                             override.aes = list(fill = "grey65",
                                                 colour = "black"))) +
  facet_wrap(~ site, ncol = 5) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "grey90", fill = NA),
        strip.background = element_blank(),
        strip.text = element_text(size = 8))

save_fig("genet_maps.pdf", p_map, 7.1, 8.2)

### A: share of clonal pairs per pit and depth separation (counts per cell)
p_mid <- max(sep$p_clonal) / 2

p_sep <- ggplot(sep, aes(x = factor(d_h), y = factor(d_v), fill = p_clonal)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%s / %s", scales::comma(n_clonal),
                                scales::comma(n_pairs)),
                colour = p_clonal > p_mid),
            size = 2.6, show.legend = FALSE) +
  scale_fill_gradient(low = col_none, high = "#08519C",
                      name = "Proportion clonal pairs") +
  scale_colour_manual(values = c(`FALSE` = "grey20", `TRUE` = "white")) +
  guides(fill = guide_colourbar(barwidth = 8, barheight = 0.4)) +
  labs(x = "Distance between pits (m)", y = "Distance between depths (m)")

### B: genotypes by pits and depth intervals occupied (empty cells kept)
grid <- span %>%
  count(n_pits, n_depths) %>%
  complete(n_pits = seq_len(nrow(design)), n_depths = seq_len(nrow(depths)),
           fill = list(n = 0))

g_mid <- max(grid$n) / 2

p_grid <- ggplot(grid, aes(x = factor(n_pits), y = factor(n_depths))) +
  geom_tile(aes(fill = n), colour = "white", linewidth = 0.8) +
  geom_text(aes(label = n, colour = n > g_mid),
            size = 3.2, show.legend = FALSE) +
  scale_fill_gradient(low = col_none, high = "#08519C", name = "Genotypes") +
  scale_colour_manual(values = c(`FALSE` = "grey20", `TRUE` = "white")) +
  guides(fill = guide_colourbar(barwidth = 8, barheight = 0.4)) +
  labs(x = "Soil pits occupied", y = "Depth intervals occupied")

save_fig("genet_axes.pdf",
         p_sep / p_grid + plot_annotation(tag_levels = "A"), 7.1, 8)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  clonal_pairs_by_separation.tsv\n")
cat("  genet_extent.tsv\n")
cat("  sampling_design.pdf\n")
cat("  genet_maps.pdf\n")
cat("  genet_axes.pdf\n")
