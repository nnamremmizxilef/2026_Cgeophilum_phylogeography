## ============================================================================
## 022_index_association.R
## Standardised index of association (rbarD) per unit at K = 4, 5 and 6 on one
## SNP per single-copy orthologue, for the sample sets of 021 (all, no_mixed,
## equal_n)
##
## Files needed:
##   consensus_phylogeny_ingroup_decloned_maf01_1snp_per_sco.vcf.gz
##                              one SNP per SCO, MLG representatives,
##                              from hpc/02_structure_phylogeny/011_filtering_clustering.sh
##                              (copied to local/02_structure_phylogeny/data/014_data_clustering
##                               for local/02_structure_phylogeny/014_clustering.R)
##   units/<set>/<unit>.txt     samples per sample set and unit
##   unit_levels.tsv            units of each level
##                              (both from hpc/03_delimitation/021_phi_test.sh;
##                               copy units/ and unit_levels.tsv from
##                               ${ROOT}/03_delimitation/results/021_results_phi_test
##                               to local/03_delimitation/data/022_data_index_association)
## ============================================================================

library(tidyverse)
library(ggpubr)
library(patchwork)
library(vcfR)
library(adegenet)
library(poppr)
library(parallel)



# files ------------------------------------------------------------------------

vcf_file   <- "02_structure_phylogeny/data/014_data_clustering/consensus_phylogeny_ingroup_decloned_maf01_1snp_per_sco.vcf.gz"
units_root <- "03_delimitation/data/022_data_index_association"
out_root   <- "03_delimitation/results/022_results_index_association"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



# parameters -------------------------------------------------------------------

### permutations of alleles among individuals (ia method 1)
n_perm <- 999

### significance
alpha <- 0.05

### cores and seed (parallel streams, reproducible for a given number of cores)
n_cores <- 4
seed    <- 42



# style ------------------------------------------------------------------------

### unit colours
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

### sample set labels and shapes
lab_set   <- c(all = "All representatives", no_mixed = "Without mixed-MAT samples",
               equal_n = "Random subsets of equal size")
shape_set <- c(all = 21, no_mixed = 24, equal_n = 22)

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf)
}



# read -------------------------------------------------------------------------

### haploid genotypes (as in 014)
vcf <- read.vcfR(vcf_file, verbose = FALSE)

gi <- vcfR2genind(vcf, return.alleles = TRUE)
ploidy(gi) <- 1L

gt <- extract.gt(vcf, element = "GT")

### sample sets and units of each level from 021
lev <- read_tsv(file.path(units_root, "unit_levels.tsv"), col_types = "cc")

sets <- tibble(path = list.files(file.path(units_root, "units"), pattern = "\\.txt$",
                                 recursive = TRUE, full.names = TRUE)) %>%
  mutate(set     = basename(dirname(path)),
         unit    = str_remove(basename(path), "\\.txt$"),
         group   = str_remove(set, "_[0-9]+$"),
         samples = map(path, read_lines)) %>%
  select(-path)

message("samples: ", nInd(gi), "   loci: ", nLoc(gi), "   sample sets x units: ", nrow(sets))



# index of association ---------------------------------------------------------

### rbarD of one unit in one sample set (loci with the minor allele in fewer
### than two individuals dropped)
run_ia <- function(s) {
  g <- gi[indNames(gi) %in% s, , drop = TRUE]
  g <- informloci(g, cutoff = 2 / nInd(g), quiet = TRUE)
  r <- ia(g, sample = n_perm, method = 1, quiet = TRUE, plot = FALSE)
  tibble(n = nInd(g), loci = nLoc(g), Ia = r[["Ia"]], p_Ia = r[["p.Ia"]],
         rbarD = r[["rbarD"]], p_rD = r[["p.rD"]])
}

RNGkind("L'Ecuyer-CMRG")
set.seed(seed)

t0  <- Sys.time()
res <- mclapply(sets$samples, run_ia, mc.cores = n_cores, mc.set.seed = TRUE)
message("index of association: ", round(difftime(Sys.time(), t0, units = "mins"), 1), " min")

## check if every analysis ran
failed <- map_lgl(res, inherits, "try-error")

if (any(failed)) stop("index of association failed for ", sum(failed),
                      " sample sets and units: ", res[[which(failed)[1]]])

per_set <- bind_cols(select(sets, set, group, unit), bind_rows(res)) %>%
  arrange(unit, set)

### per sample-set group (mean and range over the equal-size draws)
per_group <- per_set %>%
  summarise(n          = first(n),
            draws      = n(),
            loci       = mean(loci),
            rbarD_mean = mean(rbarD),
            rbarD_min  = min(rbarD),
            rbarD_max  = max(rbarD),
            n_sig      = sum(p_rD < alpha),
            .by = c(group, unit))

n_sub   <- per_group$n[per_group$group == "equal_n"][1]
n_draws <- per_group$draws[per_group$group == "equal_n"][1]

lab_set["equal_n"] <- sprintf("Random subsets of %d representatives", n_sub)



# checks -----------------------------------------------------------------------

## check if the genotypes are haploid calls
if (!all(gt %in% c("0", "1", NA))) {
  warning("genotypes other than 0, 1 and missing: ",
          paste(setdiff(unique(gt), c("0", "1", NA)), collapse = " "))
}

## check if every sample of the sample sets is in the VCF
n_absent <- length(setdiff(unlist(sets$samples), indNames(gi)))

if (n_absent > 0) {
  warning(n_absent, " samples of the sample sets are not in the VCF")
}

## check if every sample set and unit holds the samples of its list
if (any(per_set$n != lengths(sets$samples[match(paste(per_set$set, per_set$unit),
                                                paste(sets$set, sets$unit))]))) {
  warning("a sample set and unit lost samples in the subset")
}

## check if every unit of every level was analysed in every sample set
if (!all(lev$unit %in% per_set$unit) ||
    nrow(per_set) != n_distinct(per_set$set) * n_distinct(per_set$unit)) {
  warning("not every unit was analysed in every sample set")
}

## check if every analysis kept enough loci and returned a value
if (any(per_set$loci < 10 | is.na(per_set$rbarD))) {
  warning(sum(per_set$loci < 10 | is.na(per_set$rbarD)),
          " sample sets and units with fewer than 10 loci or no rbarD")
}



# figures ----------------------------------------------------------------------

### one panel per level
plot_level <- function(kk, y_title) {
  d <- per_group %>%
    inner_join(filter(lev, k == kk), by = "unit") %>%
    mutate(group = factor(group, levels = names(lab_set)))

  ggplot(d, aes(x = unit, y = rbarD_mean, fill = unit, shape = group, group = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey40") +
    geom_errorbar(aes(ymin = rbarD_min, ymax = rbarD_max), width = 0, linewidth = 0.3,
                  colour = "grey30", position = position_dodge(width = 0.6)) +
    geom_point(size = 2, stroke = 0.3, colour = "black", position = position_dodge(width = 0.6)) +
    scale_fill_manual(values = pal_unit, guide = "none") +
    scale_shape_manual(values = shape_set, labels = lab_set, name = NULL, drop = FALSE) +
    scale_y_continuous(limits = c(min(0, per_set$rbarD), max(per_set$rbarD))) +
    guides(shape = guide_legend(override.aes = list(fill = "grey70"))) +
    labs(x = sprintf("Units at K = %s", sub("k", "", kk)),
         y = if (y_title) expression(bar(r)[d]) else NULL)
}

levels_k <- sort(unique(lev$k))
p_ia <- wrap_plots(map2(levels_k, seq_along(levels_k) == 1, plot_level), nrow = 1,
                   widths = table(lev$k)[levels_k]) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

save_fig("index_association.pdf", p_ia, 7.1, 3)



# tables -----------------------------------------------------------------------

fmt_p <- function(p) sprintf("%.3f", p)

### one row per level and unit
pub <- lev %>%
  left_join(filter(per_set, set == "all"), by = "unit") %>%
  left_join(per_set %>%
              filter(set == "no_mixed") %>%
              select(unit, n_no_mixed = n, rbarD_no_mixed = rbarD, p_no_mixed = p_rD),
            by = "unit") %>%
  left_join(per_group %>%
              filter(group == "equal_n") %>%
              select(unit, rbarD_equal = rbarD_mean, rbarD_equal_min = rbarD_min,
                     rbarD_equal_max = rbarD_max, sig_equal = n_sig),
            by = "unit") %>%
  transmute(Level                        = sub("k", "K = ", k),
            Unit                         = unit,
            n                            = n,
            Loci                         = loci,
            rbarD                        = sprintf("%.4f", rbarD),
            p                            = fmt_p(p_rD),
            `Without mixed MAT: n`       = n_no_mixed,
            `Without mixed MAT: rbarD`   = sprintf("%.4f (p %s)", rbarD_no_mixed, fmt_p(p_no_mixed)),
            `Equal size: rbarD`          = sprintf("%.4f [%.4f to %.4f]", rbarD_equal,
                                                   rbarD_equal_min, rbarD_equal_max),
            `Equal size: draws p < 0.05` = sprintf("%d of %d", sig_equal, n_draws))

write_tsv(per_set, file.path(out_root, "index_association_sets.tsv"))
write_tsv(pub, file.path(out_root, "publication_index_association.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out_root, "sessionInfo.txt"))

cat("\n=== rbarD per level and unit ===\n")
per_group %>%
  select(group, unit, rbarD_mean) %>%
  pivot_wider(names_from = group, values_from = rbarD_mean) %>%
  inner_join(lev, by = "unit") %>%
  relocate(k) %>%
  arrange(k, unit) %>%
  mutate(across(where(is.double), ~ round(.x, 4))) %>%
  print(n = Inf)

cat("\n=== Permutation p, all representatives ===\n")
per_set %>% filter(set == "all") %>% select(unit, n, loci, rbarD, p_rD) %>% print(n = Inf)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  index_association.pdf\n")
cat("  index_association_sets.tsv\n")
cat("  publication_index_association.tsv\n")
cat("  sessionInfo.txt\n")
