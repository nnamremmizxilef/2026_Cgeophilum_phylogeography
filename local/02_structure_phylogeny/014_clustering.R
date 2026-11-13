## ============================================================================
## 014_clustering.R
## DAPC and sNMF clustering of the thinned call set, checked against coancestry
## clustering (fineRADstructure)
##
## Files needed:
##   consensus_phylogeny_ingroup_decloned_maf01_1snp_per_sco.vcf.gz   primary draw
##   replicates/..._1snp_per_sco_rep01 ... rep20.vcf.gz               replicate draws
##                     (both from hpc/02_structure_phylogeny/011_filtering_clustering.sh;
##                      copy the primary VCF and replicates/ from
##                      ${ROOT}/02_structure_phylogeny/results/011_results_filtering_clustering
##                      to local/02_structure_phylogeny/data/014_data_clustering)
##   fineradstructure_populations.tsv                                 MAP populations
##                     (from hpc/02_structure_phylogeny/014_clustering.sh; copy
##                      ${ROOT}/02_structure_phylogeny/results/014_results_clustering/fineradstructure_populations.tsv
##                      to local/02_structure_phylogeny/data/014_data_clustering)
##   sample_metadata_ALL_FIXED.tsv                                    site, project, clones
##                     (from local/01_variant_calling/103_mating_type.R)
## ============================================================================

library(vcfR)
library(adegenet)
library(LEA)
library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

data_root <- "02_structure_phylogeny/data/014_data_clustering"
out_root  <- "02_structure_phylogeny/results/014_results_clustering"

dir_dapc <- file.path(out_root, "dapc")
dir_snmf <- file.path(out_root, "snmf")
dir_sum  <- file.path(out_root, "summary")

walk(c(dir_dapc, dir_snmf, dir_sum),
     dir.create, recursive = TRUE, showWarnings = FALSE)

f_vcf   <- file.path(data_root, "consensus_phylogeny_ingroup_decloned_maf01_1snp_per_sco.vcf.gz")
f_meta  <- "01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv"
rep_dir <- file.path(data_root, "replicates")
f_frs   <- file.path(data_root, "fineradstructure_populations.tsv")



# parameters -------------------------------------------------------------------

### K range of the primary run
k_range <- 2:20

### ancestry barplot range
k_plot <- 2:10

### levels for the full test (from where BIC and cross-entropy level off to where
### the methods disagree; empty: only the diagnostics across k_range)
k_focus <- 4:7

### levels carried forward as lineages (the first also as lineage_assignment.tsv)
k_lineage <- 4:7

### replicate draws (curves at every K, partitions at k_focus)
run_replicates <- TRUE

### sNMF repetitions per K, k-means restarts inside find.clusters, seed
n_reps  <- 10
n_start <- 100
seed    <- 42

### ancestry reference (hybrid threshold for STRUCTURE), counted per K
q_ref <- 0.7

### clusters smaller than this counted as residual bins rather than groups
n_min_cluster <- 5

### per-sample missingness filter of the sensitivity run
miss_max <- 0.20



# style ------------------------------------------------------------------------

pal_cluster <- viridisLite::viridis

pal_project <- c(magic = "white", phytoakmeter = "grey25")
lab_project <- c(magic = "MAGIC", phytoakmeter = "PhytOakmeter")

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(path, plot, width, height) {
  ggsave(path, plot, width = width, height = height, device = cairo_pdf)
}

k_lab <- labeller(K = function(x) paste("K =", x))



# read -------------------------------------------------------------------------

meta <- read_tsv(f_meta, show_col_types = FALSE)

vcf <- read.vcfR(f_vcf, verbose = FALSE)

gi <- vcfR2genind(vcf, return.alleles = TRUE)
ploidy(gi) <- 1L

sample_ids <- indNames(gi)
n_samples  <- length(sample_ids)

### keep all principal components for the clustering step
n_pca_all <- n_samples - 1L

### upper end of the alpha-score scan
n_pca_max <- min(100L, n_pca_all)

### grouping by metadata
sample_info <- tibble(sample = sample_ids) %>%
  left_join(select(meta, sample, project, country, site), by = "sample") %>%
  arrange(project, country, site, sample) %>%
  mutate(order = row_number())

gt_matrix <- extract.gt(vcf, element = "GT")

### fineRADstructure MAP population per sample
frs     <- read_tsv(f_frs, show_col_types = FALSE)
frs_pop <- frs$population[match(sample_ids, frs$sample)]

message("samples: ", n_samples, "   loci: ", nLoc(gi))

writeLines(capture.output(sessionInfo()), file.path(dir_sum, "sessionInfo.txt"))



# checks -----------------------------------------------------------------------

## check if the metadata covers every sample
if (anyNA(sample_info$project)) {
  warning(sum(is.na(sample_info$project)), " samples without metadata, check ",
          "that sample_metadata_ALL_FIXED.tsv was copied from the same run as the VCF")
}

## check if the metadata carries clone columns
if (!all(c("in_ingroup", "mlg", "clone_status") %in% names(meta))) {
  stop("sample_metadata_ALL_FIXED.tsv lacks in_ingroup, mlg or clone_status")
}

## check if metadata carries run and depth for the missingness block
if (!all(c("run_id", "mean_seq_depth", "primary_mapped_pct") %in% names(meta))) {
  stop("sample_metadata_ALL_FIXED.tsv lacks run_id, mean_seq_depth or primary_mapped_pct")
}

## check if any sample has no site
if (anyNA(sample_info$site)) {
  warning(sum(is.na(sample_info$site)), " samples without a site, the per-site ",
          "tables will hold an NA column and the project comparison is not ",
          "interpretable")
}

## check if the genotypes are haploid
unexpected <- setdiff(na.omit(unique(as.vector(gt_matrix))), c("0", "1", "."))

if (length(unexpected) > 0) {
  stop("unexpected GT codes for haploid data: ", paste(unexpected, collapse = ", "))
}

## check if the top of the K range is identifiable
if (n_samples / max(k_range) < n_min_cluster) {
  warning("mean cluster size is below ", n_min_cluster, " at K = ", max(k_range),
          ", treat the top of the range as exploratory")
}

## check if K levels are in k_range
if (!all(k_focus %in% k_range) || !all(k_plot %in% k_range)) {
  stop("k_focus and k_plot have to lie inside k_range")
}

## check if the lineage levels are in k_range
if (!all(k_lineage %in% k_range)) {
  stop("every k_lineage has to lie inside k_range")
}

## check if fineRADstructure covers every sample
if (anyNA(frs_pop)) {
  stop(sum(is.na(frs_pop)), " samples without a fineRADstructure population")
}



# dapc -------------------------------------------------------------------------

### k-means on the principal components (BIC over the whole range, then the
### partition per K)
set.seed(seed)
dapc_scan <- find.clusters(gi, max.n.clust = max(k_range), n.pca = n_pca_all,
                           n.clust = NULL, stat = "BIC", n.start = n_start,
                           choose.n.clust = FALSE)

bic <- tibble(K = seq_along(dapc_scan$Kstat), BIC = as.numeric(dapc_scan$Kstat))

grp_dapc <- set_names(map(k_range, function(k) {
  set.seed(seed)
  as.integer(find.clusters(gi, n.pca = n_pca_all, n.clust = k, n.start = n_start,
                           choose.n.clust = FALSE)$grp)
}), as.character(k_range))



# snmf -------------------------------------------------------------------------

### .geno input: loci in rows, samples in columns, 0, 1 and 9 (missing) at
### ploidy 1
geno <- matrix(9L, nrow = nrow(gt_matrix), ncol = ncol(gt_matrix),
               dimnames = dimnames(gt_matrix))
geno[which(gt_matrix == "0")] <- 0L
geno[which(gt_matrix == "1")] <- 1L

stopifnot(identical(colnames(geno), sample_ids))

message("missing genotypes: ", sprintf("%.2f%%", 100 * mean(geno == 9L)))

f_geno <- file.path(dir_snmf, "input.geno")
write.table(geno, f_geno, col.names = FALSE, row.names = FALSE, sep = "",
            quote = FALSE)

snmf_project <- snmf(f_geno, K = k_range, entropy = TRUE, repetitions = n_reps,
                     project = "new", ploidy = 1, seed = seed, CPU = 4)

ce <- tibble(K = rep(k_range, each = n_reps),
             ce = as.numeric(unlist(map(k_range,
                                        ~ cross.entropy(snmf_project, K = .x)))))

ce_mean <- ce %>% group_by(K) %>% summarise(mean_ce = mean(ce), .groups = "drop")

### run with the lowest cross-entropy per K
snmf_best <- set_names(map_int(k_range,
                               ~ which.min(cross.entropy(snmf_project, K = .x))),
                       as.character(k_range))

q_best <- set_names(map(k_range, function(k) {
  qm <- Q(snmf_project, K = k, run = snmf_best[[as.character(k)]])
  colnames(qm) <- paste0("Anc", seq_len(k))
  rownames(qm) <- sample_ids
  qm
}), as.character(k_range))



# per-sample reproducibility ---------------------------------------------------

### hard assignment of every repetition at one K
hard_runs <- function(k) {
  map(seq_len(n_reps), ~ apply(Q(snmf_project, K = k, run = .x), 1, which.max))
}

### adjusted Rand index between repetitions (a K less reproducible than its
### neighbours forces a split)
run_ari <- function(hard) {
  combn(n_reps, 2, function(i) mclust::adjustedRandIndex(hard[[i[1]]], hard[[i[2]]]))
}

### consistency per sample across repetitions
sample_consistency <- function(hard, best) {
  ref      <- hard[[best]]
  ref_sets <- split(seq_along(ref), ref)
  acc      <- rep(0, length(ref))
  for (h in hard) {
    h_sets <- split(seq_along(h), h)
    a <- ref_sets[as.character(ref)]
    b <- h_sets[as.character(h)]
    acc <- acc + map2_dbl(a, b, ~ length(intersect(.x, .y)) / length(union(.x, .y)))
  }
  acc / length(hard)
}

### samples leaving their cluster when K grows by one
n_crossing <- function(prev, cur) {
  tab <- table(cur, prev)
  sum(rowSums(tab) - apply(tab, 1, max))
}

hard_all <- set_names(map(k_range, hard_runs), as.character(k_range))

grp_snmf <- set_names(map(as.character(k_range),
                          ~ hard_all[[.x]][[snmf_best[[.x]]]]),
                      as.character(k_range))



# per-sample table across the whole range --------------------------------------

quality <- map_dfr(k_range, function(k) {
  kk <- as.character(k)
  tibble(K                = k,
         sample           = sample_ids,
         dapc             = grp_dapc[[kk]],
         snmf             = as.integer(grp_snmf[[kk]]),
         snmf_max_q       = apply(q_best[[kk]], 1, max),
         snmf_consistency = sample_consistency(hard_all[[kk]], snmf_best[[kk]]))
}) %>%
  left_join(select(sample_info, sample, project, country, site), by = "sample")

write_tsv(quality, file.path(dir_sum, "assignment_quality.tsv"))



# K diagnostics ----------------------------------------------------------------

### diagnostics per K (ARI compares partitions regardless of cluster labels)
k_diag <- tibble(K = k_range) %>%
  mutate(
    rep_ARI_mean  = map_dbl(as.character(K), ~ mean(run_ari(hard_all[[.x]]))),
    rep_ARI_min   = map_dbl(as.character(K), ~ min(run_ari(hard_all[[.x]]))),
    ARI_dapc_snmf = map_dbl(as.character(K),
                            ~ mclust::adjustedRandIndex(grp_dapc[[.x]],
                                                        grp_snmf[[.x]])),
    min_cluster_n = map_int(as.character(K), ~ min(tabulate(grp_dapc[[.x]]))),
    n_small       = map_int(as.character(K),
                            ~ sum(tabulate(grp_dapc[[.x]]) < n_min_cluster)),
    crossing      = c(NA_real_,
                      map2_dbl(as.character(K[-1]), as.character(K[-length(K)]),
                               ~ n_crossing(grp_dapc[[.y]], grp_dapc[[.x]]))),
    ### fineRADstructure against DAPC and sNMF
    ### frs_crossing: samples of fineRADstructure populations split by the level
    ARI_dapc_frs  = map_dbl(as.character(K),
                            ~ mclust::adjustedRandIndex(grp_dapc[[.x]], frs_pop)),
    ARI_snmf_frs  = map_dbl(as.character(K),
                            ~ mclust::adjustedRandIndex(grp_snmf[[.x]], frs_pop)),
    frs_crossing  = map_dbl(as.character(K),
                            ~ n_crossing(grp_dapc[[.x]], frs_pop))) %>%
  left_join(bic, by = "K") %>%
  left_join(ce_mean, by = "K") %>%
  mutate(BIC_change = BIC - lag(BIC),
         ce_change  = mean_ce - lag(mean_ce)) %>%
  left_join(quality %>%
              group_by(K) %>%
              summarise(median_max_q       = median(snmf_max_q),
                        n_below_q          = sum(snmf_max_q < q_ref),
                        median_consistency = median(snmf_consistency),
                        .groups = "drop"),
            by = "K") %>%
  select(K, BIC, BIC_change, mean_ce, ce_change, rep_ARI_mean, rep_ARI_min,
         ARI_dapc_snmf, median_max_q, n_below_q, median_consistency,
         min_cluster_n, n_small, crossing, ARI_dapc_frs, ARI_snmf_frs, frs_crossing)

write_tsv(k_diag, file.path(dir_sum, "k_diagnostics.tsv"))

cat("\n=== diagnostics across K ===\n")
k_diag %>%
  mutate(across(where(is.double), ~ round(.x, 3))) %>%
  print(n = Inf)

cat("\n=== fineRADstructure ===\n")
cat("MAP populations:", n_distinct(frs_pop), "  sizes:",
    paste(sort(table(frs_pop), decreasing = TRUE), collapse = " "), "\n")

### one panel per diagnostic on the same K axis
diag_long <- k_diag %>%
  select(K, BIC, `cross-entropy` = mean_ce,
         `sNMF run-to-run ARI (min)` = rep_ARI_min,
         `DAPC vs sNMF ARI` = ARI_dapc_snmf,
         `median max ancestry` = median_max_q,
         `median consistency` = median_consistency,
         `smallest cluster (n)` = min_cluster_n,
         `samples crossing from K-1` = crossing,
         `DAPC vs fineRADstructure ARI` = ARI_dapc_frs,
         `fineRADstructure samples split` = frs_crossing) %>%
  pivot_longer(-K, names_to = "metric", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(metric = factor(metric, levels = c(
    "BIC", "cross-entropy",
    "sNMF run-to-run ARI (min)", "DAPC vs sNMF ARI",
    "median max ancestry", "median consistency",
    "smallest cluster (n)", "samples crossing from K-1",
    "DAPC vs fineRADstructure ARI", "fineRADstructure samples split")))

p_diag <- ggplot(diag_long, aes(x = K, y = value)) +
  geom_vline(xintercept = k_range, colour = "grey85", linetype = "dashed",
             linewidth = 0.25) +
  geom_line(linewidth = 0.6, na.rm = TRUE) +
  geom_point(shape = 21, colour = "black", fill = "black", size = 1.4,
             na.rm = TRUE) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2, axes = "all_x") +
  scale_x_continuous(breaks = k_range) +
  labs(x = "Number of clusters (K)", y = NULL)

save_fig(file.path(dir_sum, "k_diagnostics.pdf"), p_diag, 7.1, 9)



# ancestry barplots across K ---------------------------------------------------

### stacked ancestry barplots across k_plot
anc_levels <- paste0("Anc", seq_len(max(k_plot)))
anc_cols   <- set_names(pal_cluster(max(k_plot)), anc_levels)

qdf <- map_dfr(k_plot, function(k) {
  as_tibble(q_best[[as.character(k)]]) %>%
    mutate(sample = sample_ids, K = k) %>%
    pivot_longer(starts_with("Anc"), names_to = "ancestry",
                 values_to = "proportion")
}) %>%
  left_join(sample_info, by = "sample") %>%
  mutate(ancestry = factor(ancestry, levels = anc_levels))

p_bar <- ggplot(qdf, aes(x = reorder(sample, order), y = proportion,
                         fill = ancestry)) +
  geom_col(width = 1) +
  scale_fill_manual(values = anc_cols, name = NULL) +
  facet_grid(K ~ project + country, scales = "free_x", space = "free_x",
             labeller = labeller(K = function(x) paste("K =", x))) +
  labs(x = NULL, y = "Ancestry proportion") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text.x = element_text(size = 7),
        strip.text.y = element_text(size = 8),
        panel.spacing.x = unit(0.4, "mm"),
        panel.spacing.y = unit(1, "mm"))

save_fig(file.path(dir_snmf, "snmf_barplots.pdf"), p_bar,
         7.1, 0.65 * length(k_plot) + 1)



# replicate draws --------------------------------------------------------------

### same SNP pool with different random draws
f_rep <- sort(list.files(rep_dir, pattern = "_1snp_per_sco_rep[0-9]+\\.vcf\\.gz$",
                         full.names = TRUE))

### one draw: the two curves over k_range and the partitions at k_focus
cluster_draw <- function(f) {

  draw <- str_extract(basename(f), "rep[0-9]+")
  message("replicate ", draw)

  v  <- read.vcfR(f, verbose = FALSE)
  g  <- vcfR2genind(v, return.alleles = TRUE)
  ploidy(g) <- 1L
  gt <- extract.gt(v, element = "GT")

  ## check if the draw holds the same samples in the same order as the primary
  stopifnot(identical(indNames(g), sample_ids))

  ### dapc
  set.seed(seed)
  scan <- find.clusters(g, max.n.clust = max(k_range), n.pca = n_pca_all,
                        n.clust = NULL, stat = "BIC", n.start = n_start,
                        choose.n.clust = FALSE)

  dapc_k <- set_names(map(k_focus, function(k) {
    set.seed(seed)
    as.integer(find.clusters(g, n.pca = n_pca_all, n.clust = k, n.start = n_start,
                             choose.n.clust = FALSE)$grp)
  }), as.character(k_focus))

  ### snmf (one project folder per draw)
  d <- file.path(dir_snmf, "replicates", draw)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

  geno_r <- matrix(9L, nrow = nrow(gt), ncol = ncol(gt))
  geno_r[which(gt == "0")] <- 0L
  geno_r[which(gt == "1")] <- 1L

  f_g <- file.path(d, "input.geno")
  write.table(geno_r, f_g, col.names = FALSE, row.names = FALSE, sep = "",
              quote = FALSE)

  proj <- snmf(f_g, K = k_range, entropy = TRUE, repetitions = n_reps,
               project = "new", ploidy = 1, seed = seed, CPU = 4)

  snmf_k <- set_names(map(k_focus, function(k) {
    best <- which.min(cross.entropy(proj, K = k))
    apply(Q(proj, K = k, run = best), 1, which.max)
  }), as.character(k_focus))

  curves <- tibble(draw = draw, K = k_range) %>%
    left_join(tibble(K = seq_along(scan$Kstat), BIC = as.numeric(scan$Kstat)),
              by = "K") %>%
    mutate(mean_ce = map_dbl(K, ~ mean(cross.entropy(proj, K = .x))))

  list(draw = draw, curves = curves, dapc = dapc_k, snmf = snmf_k)
}


### per sample, the share of its primary cluster it stays with in a draw
keep_company <- function(ref, others) {
  a <- split(seq_along(ref), ref)[as.character(ref)]
  rowMeans(sapply(others, function(h) {
    b <- split(seq_along(h), h)[as.character(h)]
    map2_dbl(a, b, ~ length(intersect(.x, .y)) / length(union(.x, .y)))
  }))
}


if (run_replicates && length(f_rep) > 0) {

  reps <- map(f_rep, cluster_draw)

  ### the curves of every draw against the primary
  rep_curves <- bind_rows(
    tibble(draw = "primary", K = k_range) %>%
      left_join(bic, by = "K") %>%
      left_join(ce_mean, by = "K"),
    map_dfr(reps, "curves"))

  write_tsv(rep_curves, file.path(dir_sum, "replicate_curves.tsv"))

  rep_long <- rep_curves %>%
    pivot_longer(c(BIC, mean_ce), names_to = "metric", values_to = "value") %>%
    mutate(metric = if_else(metric == "mean_ce", "cross-entropy", metric))

  p_rep <- ggplot(rep_long, aes(x = K, y = value, group = draw)) +
    geom_vline(xintercept = k_range, colour = "grey85", linetype = "dashed",
               linewidth = 0.25) +
    geom_line(data = function(d) filter(d, draw != "primary"),
              colour = "grey70", linewidth = 0.4) +
    geom_line(data = function(d) filter(d, draw == "primary"),
              colour = "black", linewidth = 0.8) +
    geom_point(data = function(d) filter(d, draw == "primary"),
               shape = 21, colour = "black", fill = "black", size = 1.4) +
    facet_wrap(~ metric, scales = "free_y", ncol = 1, axes = "all_x") +
    scale_x_continuous(breaks = k_range) +
    labs(x = "Number of clusters (K)", y = NULL)

  save_fig(file.path(dir_sum, "replicate_curves.pdf"), p_rep, 7.1, 5.5)


  if (length(k_focus) > 0) {

    ### each draw against the primary per method
    ### (plus the two methods within the draw)
    rep_ari <- map_dfr(reps, function(r) {
      map_dfr(as.character(k_focus), function(kk) {
        tibble(draw            = r$draw,
               K               = as.integer(kk),
               dapc_vs_primary = mclust::adjustedRandIndex(r$dapc[[kk]], grp_dapc[[kk]]),
               snmf_vs_primary = mclust::adjustedRandIndex(r$snmf[[kk]], grp_snmf[[kk]]),
               dapc_vs_snmf    = mclust::adjustedRandIndex(r$dapc[[kk]], r$snmf[[kk]]))
      })
    })

    write_tsv(rep_ari, file.path(dir_sum, "replicate_ari.tsv"))

    rep_ari_long <- rep_ari %>%
      pivot_longer(-c(draw, K), names_to = "comparison", values_to = "ARI") %>%
      mutate(comparison = factor(comparison,
                                 levels = c("dapc_vs_primary", "snmf_vs_primary",
                                            "dapc_vs_snmf"),
                                 labels = c("DAPC vs primary", "sNMF vs primary",
                                            "DAPC vs sNMF")))

    cat("\n=== replicate draws against the primary ===\n")
    rep_ari_long %>%
      group_by(K, comparison) %>%
      summarise(n_draws = n(), min = min(ARI), median = median(ARI),
                max = max(ARI), .groups = "drop") %>%
      mutate(across(where(is.double), ~ round(.x, 3))) %>%
      print(n = Inf)

    p_ari <- ggplot(rep_ari_long, aes(x = comparison, y = ARI)) +
      geom_jitter(width = 0.12, height = 0, shape = 21, colour = "black",
                  fill = "grey55", size = 1.3, alpha = 0.7) +
      stat_summary(fun = median, fun.min = median, fun.max = median,
                   geom = "crossbar", width = 0.4, linewidth = 0.3,
                   colour = "#D73027") +
      facet_wrap(~ K, labeller = k_lab) +
      labs(x = NULL, y = "Adjusted Rand Index")

    save_fig(file.path(dir_sum, "replicate_ari.pdf"), p_ari, 7.1, 3.6)


    ### against the DAPC partition carried forward as lineages (per sample)
    rep_stability <- map_dfr(as.character(k_focus), function(kk) {
      tibble(K          = as.integer(kk),
             sample     = sample_ids,
             dapc       = grp_dapc[[kk]],
             stability  = keep_company(grp_dapc[[kk]], map(reps, ~ .x$dapc[[kk]])))
    }) %>%
      left_join(select(sample_info, sample, project, country, site), by = "sample")

    write_tsv(rep_stability, file.path(dir_sum, "replicate_sample_stability.tsv"))

    cat("\n=== samples that change company across the draws ===\n")
    rep_stability %>%
      group_by(K) %>%
      summarise(median_stability = median(stability),
                below_90 = sum(stability < 0.9),
                below_50 = sum(stability < 0.5),
                .groups = "drop") %>%
      mutate(across(where(is.double), ~ round(.x, 3))) %>%
      print(n = Inf)
  }

} else {
  cat("\nno replicate draws read (run_replicates is FALSE or",
      rep_dir, "holds none)\n")
}



# missing data -----------------------------------------------------------------

### per-sample missingness (as sNMF reads it), with run and depth
missing_df <- tibble(sample = sample_ids, missing = colMeans(geno == 9L)) %>%
  left_join(select(sample_info, sample, project), by = "sample") %>%
  left_join(select(meta, sample, run_id, mean_seq_depth, primary_mapped_pct),
            by = "sample")

write_tsv(missing_df, file.path(dir_sum, "sample_missingness.tsv"))

cat("\n=== per-sample missingness by run ===\n")
missing_df %>%
  summarise(n = n(), median = median(missing), max = max(missing),
            above_miss_max = sum(missing > miss_max),
            median_depth = median(mean_seq_depth), .by = c(project, run_id)) %>%
  arrange(project, run_id) %>%
  mutate(across(where(is.double), ~ round(.x, 3))) %>%
  print(n = Inf)

cat(sprintf("\nSpearman rho, missingness vs depth: %.3f   vs mapped reads: %.3f\n",
            cor(missing_df$missing, missing_df$mean_seq_depth,
                method = "spearman", use = "complete.obs"),
            cor(missing_df$missing, missing_df$primary_mapped_pct,
                method = "spearman", use = "complete.obs")))

p_depth <- ggplot(missing_df, aes(x = mean_seq_depth, y = missing, fill = project)) +
  geom_hline(yintercept = miss_max, colour = "grey60", linetype = "dashed",
             linewidth = 0.3) +
  geom_point(shape = 21, colour = "black", size = 1, alpha = 0.7) +
  scale_x_log10() +
  scale_fill_manual(values = pal_project, labels = lab_project, name = NULL,
                    na.value = "grey60") +
  guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "Mean sequencing depth", y = "Missing genotypes per sample")

save_fig(file.path(dir_sum, "missingness_vs_depth.pdf"), p_depth, 7.1, 4.5)


if (length(k_focus) > 0) {

  ### assignment measures at k_focus (DAPC-sNMF agreement as the overlap of the
  ### two co-member sets)
  miss_quality <- quality %>%
    filter(K %in% k_focus) %>%
    select(K, sample, dapc, snmf_max_q, snmf_consistency) %>%
    left_join(map_dfr(as.character(k_focus), ~ tibble(
      K                   = as.integer(.x),
      sample              = sample_ids,
      dapc_snmf_agreement = keep_company(grp_dapc[[.x]], list(grp_snmf[[.x]])))),
      by = c("K", "sample")) %>%
    left_join(missing_df, by = "sample")

  ### stability across the draws
  if (exists("rep_stability")) {
    miss_quality <- left_join(miss_quality, select(rep_stability, K, sample, stability),
                              by = c("K", "sample"))
  }

  write_tsv(miss_quality, file.path(dir_sum, "missingness_vs_assignment.tsv"))

  measure_labels <- c(snmf_max_q          = "sNMF max ancestry",
                      snmf_consistency    = "sNMF run-to-run consistency",
                      dapc_snmf_agreement = "DAPC-sNMF agreement",
                      stability           = "stability across draws")

  miss_long <- miss_quality %>%
    pivot_longer(any_of(names(measure_labels)), names_to = "measure",
                 values_to = "value")

  ### rank correlation and medians above and below miss_max
  cat("\n=== missingness against assignment ===\n")
  miss_long %>%
    summarise(rho         = suppressWarnings(cor(missing, value, method = "spearman")),
              median_high = median(value[missing > miss_max]),
              median_rest = median(value[missing <= miss_max]),
              .by = c(K, measure)) %>%
    mutate(across(where(is.double), ~ round(.x, 3))) %>%
    print(n = Inf)

  ### per DAPC cluster (a cluster made of high-missing samples is an artefact)
  miss_cluster <- miss_quality %>%
    summarise(n              = n(),
              median_missing = median(missing),
              max_missing    = max(missing),
              n_above        = sum(missing > miss_max),
              median_depth   = median(mean_seq_depth),
              .by = c(K, dapc)) %>%
    mutate(share_above = n_above / n) %>%
    arrange(K, dapc)

  write_tsv(miss_cluster, file.path(dir_sum, "missingness_by_cluster.tsv"))

  cat("\n=== missingness per DAPC cluster ===\n")
  miss_cluster %>%
    mutate(across(where(is.double), ~ round(.x, 3))) %>%
    print(n = Inf)

  cat("\n=== samples above", miss_max, "missing ===\n")
  miss_quality %>%
    filter(missing > miss_max) %>%
    select(K, sample, project, run_id, mean_seq_depth, missing, dapc) %>%
    pivot_wider(names_from = K, values_from = dapc, names_prefix = "dapc_K") %>%
    arrange(desc(missing)) %>%
    mutate(missing = round(missing, 3)) %>%
    print(n = Inf)

  p_miss <- miss_long %>%
    mutate(measure = factor(measure_labels[measure], levels = unname(measure_labels))) %>%
    ggplot(aes(x = missing, y = value, fill = project)) +
    geom_vline(xintercept = miss_max, colour = "grey60", linetype = "dashed",
               linewidth = 0.3) +
    geom_point(shape = 21, colour = "black", size = 0.8, stroke = 0.2,
               alpha = 0.6) +
    scale_fill_manual(values = pal_project, labels = lab_project, name = NULL,
                      na.value = "grey60") +
    guides(fill = guide_legend(override.aes = list(size = 2, alpha = 1))) +
    facet_grid(measure ~ K, labeller = k_lab) +
    labs(x = "Missing genotypes per sample", y = NULL)

  save_fig(file.path(dir_sum, "missingness_vs_assignment.pdf"), p_miss, 7.1, 9)


  ### sensitivity run: DAPC and sNMF at k_focus without the samples above
  ### miss_max (against the full run on the shared samples)
  keep_lo <- missing_df$missing <= miss_max
  n_lo    <- sum(keep_lo)

  ## check if any sample is above miss_max
  if (n_lo == n_samples) {
    cat("\nno sample above", miss_max, "missing, sensitivity run skipped\n")
  } else {

    gi_lo <- gi[keep_lo, , drop = TRUE]

    sens_dapc <- set_names(map(k_focus, function(k) {
      set.seed(seed)
      as.integer(find.clusters(gi_lo, n.pca = n_lo - 1L, n.clust = k,
                               n.start = n_start, choose.n.clust = FALSE)$grp)
    }), as.character(k_focus))

    d_lo <- file.path(dir_snmf, "missingness_sensitivity")
    dir.create(d_lo, recursive = TRUE, showWarnings = FALSE)

    f_lo <- file.path(d_lo, "input.geno")
    write.table(geno[, keep_lo], f_lo, col.names = FALSE, row.names = FALSE,
                sep = "", quote = FALSE)

    proj_lo <- snmf(f_lo, K = k_focus, entropy = TRUE, repetitions = n_reps,
                    project = "new", ploidy = 1, seed = seed, CPU = 4)

    sens_snmf <- set_names(map(k_focus, function(k) {
      best <- which.min(cross.entropy(proj_lo, K = k))
      apply(Q(proj_lo, K = k, run = best), 1, which.max)
    }), as.character(k_focus))

    ### ARI and reassigned samples (those not in the majority of their cluster)
    sens <- map_dfr(as.character(k_focus), function(kk) {
      tibble(K           = as.integer(kk),
             n_dropped   = n_samples - n_lo,
             n_shared    = n_lo,
             dapc_ARI    = mclust::adjustedRandIndex(sens_dapc[[kk]], grp_dapc[[kk]][keep_lo]),
             snmf_ARI    = mclust::adjustedRandIndex(sens_snmf[[kk]], grp_snmf[[kk]][keep_lo]),
             dapc_moved  = n_crossing(grp_dapc[[kk]][keep_lo], sens_dapc[[kk]]),
             snmf_moved  = n_crossing(grp_snmf[[kk]][keep_lo], sens_snmf[[kk]]))
    })

    write_tsv(sens, file.path(dir_sum, "missingness_sensitivity.tsv"))

    cat("\n=== without the", n_samples - n_lo, "samples above", miss_max, "missing ===\n")
    sens %>%
      mutate(across(where(is.double), ~ round(.x, 3))) %>%
      print(n = Inf)
  }
}



# focal K ----------------------------------------------------------------------

### the DAPC scatter and the per-sample table
write_focus <- function(k) {

  kk <- as.character(k)
  grp <- factor(grp_dapc[[kk]])

  ### number of PCs by the alpha-score (penalises overfitting)
  set.seed(seed)
  a_res <- optim.a.score(dapc(gi, grp, n.pca = n_pca_max, n.da = k - 1),
                         plot = FALSE)

  ### discriminant axes at most the retained PCs
  n_pca_opt <- as.integer(a_res$best)
  n_da      <- min(k - 1L, n_pca_opt)

  set.seed(seed)
  dapc_k <- dapc(gi, grp, n.pca = n_pca_opt, n.da = n_da)

  scores <- as_tibble(dapc_k$ind.coord) %>%
    mutate(cluster = as.character(grp))

  p_scatter <- if ("LD2" %in% names(scores)) {
    ggplot(scores, aes(x = LD1, y = LD2, fill = cluster)) +
      geom_point(shape = 21, colour = "black", size = 1.6, alpha = 0.85) +
      labs(x = "Discriminant axis 1", y = "Discriminant axis 2")
  } else {
    ggplot(scores, aes(x = LD1, fill = cluster)) +
      geom_density(alpha = 0.6, colour = NA) +
      labs(x = "Discriminant axis 1", y = "Density")
  }

  cl_cols <- set_names(pal_cluster(k), as.character(seq_len(k)))

  p_scatter <- p_scatter +
    scale_fill_manual(values = cl_cols, name = NULL)

  cat(sprintf("DAPC K = %d   PCs retained %d   discriminant axes %d\n",
              k, n_pca_opt, n_da))

  save_fig(file.path(dir_dapc, sprintf("dapc_scatter_K%d.pdf", k)), p_scatter, 7.1, 5)

  quality %>%
    filter(K == k) %>%
    select(sample, dapc, snmf, snmf_max_q, snmf_consistency,
           project, country, site) %>%
    write_tsv(file.path(dir_sum, sprintf("clusters_K%d.tsv", k)))
}

if (length(k_focus) > 0) {
  walk(k_focus, write_focus)
} else {
  cat("\nk_focus is empty, no per-level output written.\n")
}



# fineRADstructure at focal K --------------------------------------------------

### fineRADstructure populations per DAPC cluster
### (a population split across clusters does not support that split)
if (length(k_focus) > 0) {

  frs_focus <- map_dfr(as.character(k_focus), function(kk) {
    tibble(K = as.integer(kk), dapc = grp_dapc[[kk]], population = frs_pop) %>%
      count(K, dapc, population)
  })

  write_tsv(frs_focus, file.path(dir_sum, "fineradstructure_by_cluster.tsv"))

  cat("\n=== fineRADstructure populations per DAPC cluster ===\n")
  frs_focus %>%
    summarise(n = sum(n), n_populations = n_distinct(population), .by = c(K, dapc)) %>%
    arrange(K, dapc) %>%
    print(n = Inf)

  cat("\n=== fineRADstructure populations split across DAPC clusters ===\n")
  frs_focus %>%
    filter(n_distinct(dapc) > 1, .by = c(K, population)) %>%
    arrange(K, population, desc(n)) %>%
    print(n = Inf)
}



# lineages ---------------------------------------------------------------------

### lineage labels per level (the first also as lineage_assignment.tsv, the
### table read downstream)
write_lineages <- function(k, primary) {

  clustered <- tibble(sample  = sample_ids,
                      lineage = as.character(grp_dapc[[as.character(k)]]))

  ### clonemates inherit the label of their MLG (clustered: samples with their
  ### own label)
  mlg_lineage <- meta %>%
    select(sample, mlg) %>%
    inner_join(clustered, by = "sample") %>%
    filter(!is.na(mlg)) %>%
    distinct(mlg, .keep_all = TRUE) %>%
    select(mlg, lineage_mlg = lineage)

  lineages_all <- meta %>%
    filter(in_ingroup %in% TRUE) %>%
    select(sample, mlg, clone_status, project, country, site) %>%
    left_join(clustered, by = "sample") %>%
    left_join(mlg_lineage, by = "mlg") %>%
    mutate(clustered = !is.na(lineage),
           lineage   = coalesce(lineage, lineage_mlg)) %>%
    select(sample, lineage, clustered, mlg, clone_status, project, country, site) %>%
    arrange(sample)

  ## check if every ingroup sample got a label
  n_unlabelled <- sum(is.na(lineages_all$lineage))
  if (n_unlabelled > 0) {
    warning(n_unlabelled, " ingroup samples at K = ", k, " have no lineage, their ",
            "MLG has no clustered representative")
  }

  write_tsv(lineages_all, file.path(dir_sum, sprintf("lineage_assignment_K%d.tsv", k)))
  if (primary) write_tsv(lineages_all, file.path(dir_sum, "lineage_assignment.tsv"))

  cat("\n=== lineage labels at K =", k, "===\n")
  cat("clustered:", sum(lineages_all$clustered),
      "  inherited from their MLG:", sum(!lineages_all$clustered & !is.na(lineages_all$lineage)),
      "  unlabelled:", n_unlabelled, "\n")

  ### the tables below on the clustered samples (one per MLG)
  lineages <- filter(lineages_all, clustered)

  ### lineages per site (structure against sampling design)
  count(lineages, lineage, site) %>%
    pivot_wider(names_from = site, values_from = n, values_fill = 0) %>%
    write_tsv(file.path(dir_sum, sprintf("lineage_by_site_K%d.tsv", k)))

  ### lineages per project (a lineage confined to one project may track the
  ### collection)
  lineage_project <- count(lineages, lineage, project) %>%
    pivot_wider(names_from = project, values_from = n, values_fill = 0)

  write_tsv(lineage_project,
            file.path(dir_sum, sprintf("lineage_by_project_K%d.tsv", k)))

  cat("\n=== lineages at K =", k, if (primary) "(primary)" else "", "===\n")
  print(lineage_project, n = Inf)

  cat("\nsites per lineage\n")
  lineages %>%
    summarise(n_sites = n_distinct(site), n = n(), .by = lineage) %>%
    arrange(lineage) %>%
    print(n = Inf)

  cat("\nlineages per site\n")
  lineages %>%
    summarise(n_lineages = n_distinct(lineage), n = n(), .by = site) %>%
    arrange(desc(n_lineages)) %>%
    print(n = Inf)
}

if (length(k_lineage) > 0) {
  walk2(k_lineage, seq_along(k_lineage) == 1L, write_lineages)
} else {
  cat("\nk_lineage is empty, no lineage assignment written.\n")
}




# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  summary/sessionInfo.txt\n")
cat("  summary/assignment_quality.tsv          assignment per sample and K\n")
cat("  summary/k_diagnostics.tsv               every diagnostic across K\n")
cat("  summary/k_diagnostics.pdf\n")
cat("  snmf/input.geno, input.snmfProject      sNMF input and runs\n")
cat("  snmf/snmf_barplots.pdf                  ancestry across K\n")
if (run_replicates && length(f_rep) > 0) {
  cat("  summary/replicate_curves.tsv            BIC and cross-entropy of every draw\n")
  cat("  summary/replicate_curves.pdf\n")
  cat("  snmf/replicates/rep*/                   sNMF runs per draw\n")
  if (length(k_focus) > 0) {
    cat("  summary/replicate_ari.tsv               each draw against the primary at k_focus\n")
    cat("  summary/replicate_ari.pdf\n")
    cat("  summary/replicate_sample_stability.tsv  per sample, across the draws\n")
  }
}
cat("  summary/sample_missingness.tsv          missing genotypes per sample, run, depth\n")
cat("  summary/missingness_vs_depth.pdf\n")
if (length(k_focus) > 0) {
  cat("  summary/missingness_vs_assignment.tsv\n")
  cat("  summary/missingness_vs_assignment.pdf\n")
  cat("  summary/missingness_by_cluster.tsv\n")
  cat("  summary/missingness_sensitivity.tsv     partitions without the samples above", miss_max, "\n")
  cat("  snmf/missingness_sensitivity/\n")
  cat("  dapc/dapc_scatter_K*.pdf\n")
  cat("  summary/clusters_K*.tsv\n")
  cat("  summary/fineradstructure_by_cluster.tsv fineRADstructure populations per DAPC cluster\n")
}
if (length(k_lineage) > 0) {
  cat("  summary/lineage_assignment_K*.tsv       every ingroup sample, clonemates with their MLG label\n")
  cat("  summary/lineage_by_site_K*.tsv\n")
  cat("  summary/lineage_by_project_K*.tsv\n")
  cat("  summary/lineage_assignment.tsv          K =", k_lineage[1], "(read downstream)\n")
}
