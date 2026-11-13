## ============================================================================
## 070_gene_trees.R
## The gene tree panel: sites, informative sites, rate model, tips and resolution
## per orthogroup tree
##
## Files needed:
##   per_og_tree/<OG>.treefile               ML gene tree
##   per_og_tree/<OG>.iqtree                 IQ-TREE report per gene
##   per_og_empty_tips/<OG>_empty_tips.txt   tips dropped as empty
##                                           all from hpc/02_structure_phylogeny/070_gene_trees.sh
##                                           (copy per_og_tree and per_og_empty_tips of
##                                            ${ROOT}/02_structure_phylogeny/results/070_results_gene_trees
##                                            to local/02_structure_phylogeny/data/070_data_gene_trees)
## ============================================================================

library(ape)
library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

data_root <- "02_structure_phylogeny/data/070_data_gene_trees"
out_root  <- "02_structure_phylogeny/results/070_results_gene_trees"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

tree_dir  <- file.path(data_root, "per_og_tree")
empty_dir <- file.path(data_root, "per_og_empty_tips")



# parameters -------------------------------------------------------------------

n_og    <- 5616   ### orthogroups in the panel
n_panel <- 618    ### tips of the panel (samples, genomes, reference)

### internal branch length counted as zero (IQ-TREE writes its minimum, 1e-6,
### for a branch the data cannot resolve)
bl_zero <- 1e-5

### IQ-TREE minimum branch length (-blmin default)
bl_floor <- 1e-6

### share of unresolved internal branches above which a gene tree is listed
unresolved_high <- 0.5

### least resolved gene trees listed
n_show <- 20



# style ------------------------------------------------------------------------

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

f_tree <- list.files(tree_dir, pattern = "\\.treefile$", full.names = TRUE)

### first match of a pattern in the report
iq_field <- function(lines, pattern) {
  m <- str_match(lines, pattern)[, 2]
  m <- m[!is.na(m)]
  if (length(m) == 0) NA_character_ else m[1]
}

read_report <- function(path) {
  if (!file.exists(path)) {
    return(tibble(n_seq = NA_integer_, n_sites = NA_integer_,
                  n_informative = NA_integer_, n_patterns = NA_integer_,
                  model = NA_character_))
  }
  l <- read_lines(path)
  tibble(n_seq         = as.integer(iq_field(l, "^Input data: ([0-9]+) sequences")),
         n_sites       = as.integer(iq_field(l, "^Input data: .* with ([0-9]+) nucleotide sites")),
         n_informative = as.integer(iq_field(l, "^Number of parsimony informative sites: ([0-9]+)")),
         n_patterns    = as.integer(iq_field(l, "^Number of distinct site patterns: ([0-9]+)")),
         model         = coalesce(iq_field(l, "^Best-fit model according to BIC: (\\S+)"),
                                  iq_field(l, "^Model of substitution: (\\S+)")))
}

read_tree_stats <- function(path) {
  tr  <- read.tree(path)
  int <- tr$edge[, 2] > Ntip(tr)
  bl  <- tr$edge.length
  tibble(n_tips          = Ntip(tr),
         tree_length     = sum(bl),
         internal_length = sum(bl[int]),
         n_internal      = sum(int),
         n_zero_internal = sum(bl[int] < bl_zero),
         bl_internal     = list(bl[int]))
}

n_empty_of <- function(og) {
  f <- file.path(empty_dir, paste0(og, "_empty_tips.txt"))
  if (!file.exists(f)) return(NA_integer_)
  x <- str_trim(read_lines(f))
  length(unique(x[x != ""]))
}

panel <- tibble(og = str_remove(basename(f_tree), "\\.treefile$"), f_tree = f_tree) %>%
  mutate(tree    = map(f_tree, read_tree_stats),
         report  = map(str_replace(f_tree, "\\.treefile$", ".iqtree"), read_report),
         n_empty = map_int(og, n_empty_of)) %>%
  unnest(c(tree, report)) %>%
  mutate(unresolved     = n_zero_internal / n_internal,
         internal_share = internal_length / tree_length,
         informative_pct = 100 * n_informative / n_sites) %>%
  select(-f_tree)

### internal branch lengths of all gene trees
bl_internal <- unlist(panel$bl_internal)
panel <- select(panel, -bl_internal)

message("gene trees: ", nrow(panel), " of ", n_og)



# checks -----------------------------------------------------------------------

## check if every orthogroup has a tree
if (nrow(panel) != n_og) {
  warning(n_og - nrow(panel), " orthogroups without a gene tree")
}

## check if every tree has its report
if (any(is.na(panel$model))) {
  warning(sum(is.na(panel$model)), " gene trees without a model in the report")
}

## check if the site counts were read from the reports
if (any(is.na(panel$n_informative) & !is.na(panel$model))) {
  warning(sum(is.na(panel$n_informative) & !is.na(panel$model)),
          " reports without site counts, check the patterns in read_report()")
}

## check if the tips plus the empty tips give the panel
bad_tips <- filter(panel, !is.na(n_empty), n_tips + n_empty != n_panel)

if (nrow(bad_tips) > 0) {
  warning(nrow(bad_tips), " gene trees whose tips plus empty tips are not ", n_panel)
}

## check if the report and the tree agree on the tips
if (any(panel$n_seq != panel$n_tips, na.rm = TRUE)) {
  warning(sum(panel$n_seq != panel$n_tips, na.rm = TRUE),
          " gene trees whose report lists a different number of sequences")
}

## check if no internal branch is below the IQ-TREE floor (else blmin differed)
if (any(bl_internal < 0.999 * bl_floor)) {
  warning(sum(bl_internal < 0.999 * bl_floor), " internal branches below ",
          bl_floor, ", check bl_floor")
}



# summary ----------------------------------------------------------------------

q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))

panel_summary <- panel %>%
  summarise(gene_trees              = n(),
            sites_median            = median(n_sites, na.rm = TRUE),
            sites_q1                = q(n_sites, 0.25),
            sites_q3                = q(n_sites, 0.75),
            informative_median      = median(n_informative, na.rm = TRUE),
            informative_q1          = q(n_informative, 0.25),
            informative_q3          = q(n_informative, 0.75),
            tips_median             = median(n_tips),
            trees_with_dropped_tips = sum(n_empty > 0, na.rm = TRUE),
            dropped_tips_max        = max(n_empty, na.rm = TRUE),
            unresolved_median       = median(unresolved),
            unresolved_q1           = q(unresolved, 0.25),
            unresolved_q3           = q(unresolved, 0.75),
            trees_mostly_unresolved = sum(unresolved > unresolved_high),
            tree_length_median      = median(tree_length))

models <- panel %>%
  count(model, name = "genes") %>%
  arrange(desc(genes))

### resolution against informative sites
res_cor <- cor.test(panel$n_informative, panel$unresolved,
                    method = "spearman", exact = FALSE)

### one substitution on a branch, in the longest and in the median gene
bl_one_max    <- 1 / max(panel$n_sites, na.rm = TRUE)
bl_one_median <- 1 / median(panel$n_sites, na.rm = TRUE)

### half an order of magnitude either side of the threshold
bl_window <- bl_zero * 10^c(-0.5, 0.5)

### internal branches collapsed and kept, and those in the window
bl_classes <- tibble(class = c("collapsed", "kept", "window_below", "window_above"),
                     from  = c(0, bl_zero, bl_window[1], bl_zero),
                     to    = c(bl_zero, Inf, bl_zero, bl_window[2])) %>%
  mutate(branches = map2_int(from, to, ~ sum(bl_internal >= .x & bl_internal < .y)),
         pct      = 100 * branches / length(bl_internal))

## check if the collapse threshold lies between the floor and one substitution in the longest gene
if (bl_zero <= bl_floor || bl_zero >= bl_one_max) {
  warning("bl_zero is not between the floor and one substitution in the longest gene")
}



# figures ----------------------------------------------------------------------

p_inf <- ggplot(panel, aes(x = n_informative)) +
  geom_histogram(bins = 60, fill = "grey40") +
  scale_x_log10() +
  labs(x = "Parsimony informative sites per gene", y = "Gene trees")

p_res <- ggplot(panel, aes(x = n_informative, y = 100 * unresolved)) +
  geom_point(shape = 21, colour = "black", fill = "#4575B4", stroke = 0.2,
             size = 0.6, alpha = 0.3) +
  geom_hline(yintercept = 100 * unresolved_high, colour = "#D73027",
             linetype = "dashed", linewidth = 0.4) +
  scale_x_log10() +
  labs(x = "Parsimony informative sites per gene",
       y = "Unresolved internal branches (%)")

p_empty <- ggplot(filter(panel, !is.na(n_empty)), aes(x = n_empty)) +
  geom_histogram(binwidth = 1, fill = "grey40") +
  labs(x = "Tips dropped as empty", y = "Gene trees")

p_model <- ggplot(models, aes(x = reorder(model, genes), y = genes)) +
  geom_col(fill = "grey40") +
  coord_flip() +
  labs(x = NULL, y = "Gene trees")

### internal branch lengths (square root counts)
p_bl <- ggplot(tibble(length = bl_internal), aes(x = length)) +
  annotate("rect", xmin = bl_window[1], xmax = bl_window[2], ymin = 0, ymax = Inf,
           fill = "#D73027", alpha = 0.1) +
  geom_histogram(bins = 120, fill = "grey40") +
  geom_vline(xintercept = bl_zero, colour = "#D73027", linetype = "dashed",
             linewidth = 0.4) +
  geom_vline(xintercept = c(bl_one_max, bl_one_median), colour = "grey20",
             linetype = "dotted", linewidth = 0.4) +
  annotate("text", x = bl_zero, y = Inf, label = "collapse threshold",
           angle = 90, hjust = 1.1, vjust = -0.5, size = 2.8, colour = "black") +
  annotate("text", x = c(bl_one_max, bl_one_median), y = Inf,
           label = c("one substitution, longest gene", "one substitution, median gene"),
           angle = 90, hjust = 1.1, vjust = -0.5, size = 2.8, colour = "black") +
  scale_x_log10(breaks = 10^(-6:0), labels = scales::label_math(10^.x, format = log10)) +
  scale_y_sqrt(labels = scales::label_comma()) +
  labs(x = "Internal branch length (substitutions per site)",
       y = "Internal branches (square root scale)")

save_fig("gene_tree_branch_lengths.pdf", p_bl, 7.1, 5)

save_fig("gene_tree_informative_sites.pdf", p_inf, 7.1, 5)
save_fig("gene_tree_resolution.pdf", p_res, 7.1, 5)
save_fig("gene_tree_dropped_tips.pdf", p_empty, 7.1, 5)
save_fig("gene_tree_models.pdf", p_model, 7.1, 4)

save_fig("gene_trees_combined.pdf",
         (p_inf + p_res) / (p_empty + p_model) +
           plot_annotation(tag_levels = "A"), 7.1, 5.9)



# tables -----------------------------------------------------------------------

write_tsv(panel, file.path(out_root, "gene_tree_panel.tsv"))
write_tsv(panel_summary, file.path(out_root, "gene_tree_summary.tsv"))
write_tsv(models, file.path(out_root, "gene_tree_models.tsv"))
write_tsv(bl_classes, file.path(out_root, "gene_tree_branch_length_classes.tsv"))

cat("\n=== gene tree panel ===\n")
print(pivot_longer(panel_summary, everything(), names_to = "quantity",
                   values_to = "value"), n = Inf)

cat("\n=== rate model per gene ===\n")
print(models, n = Inf)

cat("\n=== resolution against informative sites ===\n")
cat(sprintf("Spearman rho = %.3f, p = %.3g\n", res_cor$estimate, res_cor$p.value))

cat("\n=== internal branch lengths ===\n")
cat(sprintf("one substitution: longest gene %.2g, median gene %.2g\n",
            bl_one_max, bl_one_median))
print(bl_classes)

cat("\n=== the", n_show, "least resolved gene trees ===\n")
panel %>%
  arrange(desc(unresolved)) %>%
  select(og, n_tips, n_sites, n_informative, unresolved, model) %>%
  print(n = n_show)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  gene_tree_informative_sites.pdf   informative sites per gene\n")
cat("  gene_tree_resolution.pdf          unresolved internal branches against information\n")
cat("  gene_tree_dropped_tips.pdf        tips dropped as empty per gene\n")
cat("  gene_tree_models.pdf              rate model per gene\n")
cat("  gene_tree_branch_lengths.pdf      internal branch lengths with the collapse threshold\n")
cat("  gene_trees_combined.pdf\n")
cat("  gene_tree_panel.tsv               one row per gene tree\n")
cat("  gene_tree_summary.tsv\n")
cat("  gene_tree_models.tsv\n")
cat("  gene_tree_branch_length_classes.tsv   internal branches collapsed, kept, and within half a decade of the threshold\n")
