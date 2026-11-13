## ============================================================================
## 020_neighbournet.R
## Neighbour-Net of the MLG representatives, coloured by lineage at K = 4, 5
## and 6 and by project, and distances within and between units
##
## Files needed:
##   neighbournet.nex        splits and equal-angle layout
##   run_log.json            fast-nnt settings and least-squares fit
##   network_ibs.mdist       1 - IBS distance matrix of the representatives
##   network_ibs.mdist.id    sample order of the matrix
##                           (all four from hpc/03_delimitation/020_neighbournet.sh;
##                            copy them from ${ROOT}/03_delimitation/results/020_results_neighbournet
##                            to local/03_delimitation/data/020_data_neighbournet)
##   lineage_labels.tsv      lineage at K = 4, 5 and 6 per sample,
##                           from local/03_delimitation/001_lineage_labels.R
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

net_root    <- "03_delimitation/data/020_data_neighbournet"
labels_file <- "03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv"
out_root    <- "03_delimitation/results/020_results_neighbournet"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### levels shown
k_levels <- c(4, 5, 6)



# style ------------------------------------------------------------------------

### unit colours (a K4 or K5 unit takes the colour of its largest K6 part)
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

### project fills (kept apart from the unit colours)
pal_project <- c(magic = "white", phytoakmeter = "grey25")
lab_project <- c(magic = "MAGIC", phytoakmeter = "PhytOakmeter")

### scale bar (1 - IBS)
scale_len <- 0.01

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

### network block of the nexus (taxa at their vertex, vertex coordinates, edges)
nex <- read_lines(file.path(net_root, "neighbournet.nex"))
net <- nex[which(str_to_upper(nex) == "BEGIN NETWORK;"):length(nex)]

i_tr  <- which(net == "TRANSLATE")
i_v   <- which(net == "VERTICES")
i_e   <- which(net == "EDGES")
i_end <- which(str_starts(net, "END;"))[1]

translate <- str_match(net[(i_tr + 1):(i_v - 1)], "^\\s*(\\d+)\\s+'([^']+)'") %>%
  as_tibble(.name_repair = ~ c("line", "vertex", "sample")) %>%
  filter(!is.na(vertex)) %>%
  transmute(vertex = as.integer(vertex), sample)

vertices <- str_match(net[(i_v + 1):(i_e - 1)], "^\\s*(\\d+)\\s+(\\S+)\\s+(\\S+)") %>%
  as_tibble(.name_repair = ~ c("line", "vertex", "x", "y")) %>%
  filter(!is.na(vertex)) %>%
  transmute(vertex = as.integer(vertex), x = as.numeric(x), y = as.numeric(y))

edges <- str_match(net[(i_e + 1):(i_end - 1)], "^\\s*\\d+\\s+(\\d+)\\s+(\\d+)\\s+s=") %>%
  as_tibble(.name_repair = ~ c("line", "from", "to")) %>%
  filter(!is.na(from)) %>%
  transmute(from = as.integer(from), to = as.integer(to)) %>%
  left_join(vertices, by = c("from" = "vertex")) %>%
  left_join(rename(vertices, xend = x, yend = y), by = c("to" = "vertex"))

### least-squares fit of the splits to the distances
fit <- as.numeric(str_match(read_file(file.path(net_root, "run_log.json")),
                            "\"fit_percent\":\\s*([0-9.]+)")[, 2])

### distance matrix with its sample order
ids <- read_table(file.path(net_root, "network_ibs.mdist.id"),
                  col_names = c("fid", "iid"), col_types = "cc")

m <- as.matrix(read.table(file.path(net_root, "network_ibs.mdist")))
dimnames(m) <- list(ids$iid, ids$iid)

### lineage labels
labels <- read_tsv(labels_file, col_types = cols(.default = "c"))

tips <- translate %>%
  left_join(vertices, by = "vertex") %>%
  left_join(labels, by = "sample")



# distances within and between units -------------------------------------------

dist_long <- as_tibble(m, rownames = "a") %>%
  pivot_longer(-a, names_to = "b", values_to = "d") %>%
  filter(a < b)

unit_dist <- map_dfr(k_levels, function(k) {
  key <- select(labels, sample, unit = all_of(paste0("k", k)))
  dist_long %>%
    left_join(rename(key, unit_a = unit), by = c("a" = "sample")) %>%
    left_join(rename(key, unit_b = unit), by = c("b" = "sample")) %>%
    mutate(unit_1 = pmin(unit_a, unit_b), unit_2 = pmax(unit_a, unit_b)) %>%
    summarise(n_pairs  = n(),
              mean_d   = mean(d),
              median_d = median(d),
              .by = c(unit_1, unit_2)) %>%
    mutate(k = k, within = unit_1 == unit_2, .before = 1) %>%
    arrange(unit_1, unit_2)
})



# checks -----------------------------------------------------------------------

## check if every sample of the matrix sits at a vertex of the network
if (!setequal(translate$sample, ids$iid)) {
  warning(length(setdiff(ids$iid, translate$sample)), " samples of the matrix are not in the network")
}

## check if every taxon vertex is a tip of the network
deg <- table(c(edges$from, edges$to))

if (any(deg[as.character(translate$vertex)] != 1)) {
  warning(sum(deg[as.character(translate$vertex)] != 1), " samples sit on internal vertices")
}

## check if the matrix is complete and symmetric
if (anyNA(m) || !isSymmetric(unname(m))) {
  warning("distance matrix has missing values or is not symmetric")
}

## check if any two representatives are identical over the SNPs
if (any(dist_long$d == 0)) {
  warning(sum(dist_long$d == 0), " pairs of representatives have a distance of 0")
}

## check if every sample has a label
if (anyNA(tips$k4)) {
  warning(sum(is.na(tips$k4)), " samples of the network have no label")
}

## check if every sample is an MLG representative
if (!all(tips$clustered == "TRUE", na.rm = TRUE)) {
  warning("the network holds samples that are not MLG representatives")
}



# figures ----------------------------------------------------------------------

### scale bar at the lower left of the layout
bar <- tibble(x = min(vertices$x), y = min(vertices$y))

### one panel per level, one by project
plot_net <- function(fill_col, pal, labs = waiver()) {
  ggplot() +
    geom_segment(data = edges, aes(x = x, y = y, xend = xend, yend = yend),
                 colour = "grey55", linewidth = 0.1) +
    geom_point(data = tips, aes(x = x, y = y, fill = .data[[fill_col]]),
               shape = 21, size = 1.2, stroke = 0.15, colour = "black") +
    geom_segment(data = bar, aes(x = x, y = y, xend = x + scale_len, yend = y),
                 linewidth = 0.3) +
    annotate("text", x = bar$x + scale_len / 2, y = bar$y, label = scale_len,
             vjust = -0.6, size = 2.5) +
    scale_fill_manual(values = pal, breaks = names(pal), labels = labs, name = NULL) +
    guides(fill = guide_legend(nrow = 1, override.aes = list(size = 2))) +
    coord_equal() +
    theme_void(base_size = 10) +
    theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
}

panels <- c(map(paste0("k", k_levels), ~ plot_net(.x, pal_unit)),
            list(plot_net("project", pal_project, lab_project)))

### panels tagged (tag style set here, theme_void replaces the global theme)
p_net <- wrap_plots(panels, ncol = 2) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold"))

save_fig("neighbournet.pdf", p_net, 7.1, 6.4)



# tables -----------------------------------------------------------------------

write_tsv(unit_dist, file.path(out_root, "unit_distances.tsv"))

cat("\n=== Neighbour-Net ===\n")
cat(sprintf("samples:  %d\n", nrow(translate)))
cat(sprintf("vertices: %d   edges: %d\n", nrow(vertices), nrow(edges)))
cat(sprintf("fit:      %.2f%%\n", fit))
cat(sprintf("1 - IBS between representatives: %.4f to %.4f\n",
            min(dist_long$d), max(dist_long$d)))

cat("\n=== Mean 1 - IBS within and between units, K = 6 ===\n")
unit_dist %>%
  filter(k == 6) %>%
  select(unit_1, unit_2, mean_d) %>%
  pivot_wider(names_from = unit_2, values_from = mean_d) %>%
  mutate(across(where(is.double), ~ round(.x, 4))) %>%
  print(n = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  neighbournet.pdf\n")
cat("  unit_distances.tsv\n")
