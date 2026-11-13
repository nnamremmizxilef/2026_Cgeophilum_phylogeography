## ============================================================================
## 100_cluster_monophyly.R
## Monophyly of the candidate clusters of 014 in the four phylogenies
##
## Files needed:
##   clusters_K4.tsv                  cluster per sample at K = 4
##   clusters_K5.tsv                  the same at K = 5
##   clusters_K6.tsv                  the same at K = 6
##                                    (all three from local/02_structure_phylogeny/014_clustering.R,
##                                     results/014_results_clustering/summary)
##   aster_collapsed_tree_rooted.tre  rooted coalescent tree, collapsed gene trees,
##                                    from local/02_structure_phylogeny/080_aster_collapsed.R
##   aster_tree_rooted.tre            rooted coalescent tree, from local/02_structure_phylogeny/081_aster.R
##   species_tree_grafted_rooted.tre  rooted concatenation tree, UFBoot of the decloned tree, clonemates grafted,
##                                    from local/02_structure_phylogeny/052_species_tree_decloned.R
##   caster_site_tree_rooted.tre      rooted CASTER-site tree, from local/02_structure_phylogeny/082_caster_site.R
##   site_concordance.cf.stat, site_concordance.cf.branch
##                                    sCF per branch of the concatenation tree,
##                                    from hpc/02_structure_phylogeny/060_site_concordance.sh
##                                    (copy both from ${ROOT}/02_structure_phylogeny/results/060_results_site_concordance
##                                     to local/02_structure_phylogeny/data/060_data_site_concordance)
##   gene_concordance_<tree>.cf.stat, gene_concordance_<tree>.cf.branch
##                                    gCF per branch of each tree (concat, aster_collapsed, aster, caster),
##                                    from hpc/02_structure_phylogeny/090_gene_concordance.sh
##                                    (copy the eight files from ${ROOT}/02_structure_phylogeny/results/090_results_gene_concordance
##                                     to local/02_structure_phylogeny/data/090_data_gene_concordance)
##   sample_metadata_ALL_FIXED.tsv    project, country and site per sample,
##                                    from local/01_variant_calling/103_mating_type.R
## ============================================================================

library(ape)
library(ggtree)
library(tidyverse)
library(ggpubr)
library(patchwork)



# files ------------------------------------------------------------------------

res_root <- "02_structure_phylogeny/results"
out_root <- file.path(res_root, "100_results_cluster_monophyly")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

f_meta <- "01_variant_calling/results/103_results_mating_type_collect/sample_metadata_ALL_FIXED.tsv"

### rooted trees of 052, 080, 081 and 082
f_trees <- c(aster_collapsed = file.path(res_root, "080_results_aster_collapsed", "aster_collapsed_tree_rooted.tre"),
             aster           = file.path(res_root, "081_results_aster", "aster_tree_rooted.tre"),
             concat          = file.path(res_root, "052_results_species_tree_decloned", "species_tree_grafted_rooted.tre"),
             caster          = file.path(res_root, "082_results_caster_site", "caster_site_tree_rooted.tre"))

### concordance factors of 060 and 090 (file prefixes of .cf.stat and .cf.branch)
f_scf <- file.path("02_structure_phylogeny/data/060_data_site_concordance", "site_concordance")
f_gcf <- setNames(file.path("02_structure_phylogeny/data/090_data_gene_concordance",
                            paste0("gene_concordance_", names(f_trees))), names(f_trees))

### clusters from 014
cluster_dir <- file.path(res_root, "014_results_clustering", "summary")



# parameters -------------------------------------------------------------------

### outgroup assemblies (mapped outgroup samples dropped by prefix)
outgroup <- c("PsefloM402_1", "PsefloM405_1", "PsefloM406_1")

### assemblies (tip name pattern)
genome_pattern <- "^(Cenge|Pseflo)"

### K levels carried forward from 014
k_focus <- c(4, 5, 6)

### clustering tested (DAPC; sNMF agrees at these K)
cluster_col <- "dapc"

### purity of the smallest containing clade for nearly monophyletic
purity_near <- 0.95

### tips printed (those breaking the most clusters)
n_show <- 25



# style ------------------------------------------------------------------------

### status colours
col_status <- c("monophyletic" = "#1A9850",
                "split by the root" = "#4575B4",
                "nearly monophyletic" = "#FDAE61",
                "not monophyletic" = "#D73027")

pal_discrete <- viridisLite::viridis

### qualitative cluster palette (viridis above 12 clusters)
pal_cluster_base <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                      "#A65628", "#F781BF", "#00868B", "#B8860B", "#666666",
                      "#7570B3", "#8C564B")

pal_cluster <- function(n) {
  if (n <= length(pal_cluster_base)) pal_cluster_base[seq_len(n)] else pal_discrete(n)
}

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf, limitsize = FALSE)
}



# read -------------------------------------------------------------------------

meta <- read_tsv(f_meta, show_col_types = FALSE)

clusters <- map_dfr(k_focus, function(k) {
  f <- file.path(cluster_dir, sprintf("clusters_K%d.tsv", k))
  read_tsv(f, show_col_types = FALSE) %>%
    transmute(K = k, sample, cluster = as.character(.data[[cluster_col]]))
})

### rooted trees of 052, 080, 081 and 082 (outgroup kept, for the figures)
trees_fig <- f_trees %>% map(read.tree)

### the same trees without the outgroup, for the test (no cluster; it would
### inflate the clade sizes)
trees <- map(trees_fig, function(tr) {
  drop <- c(intersect(outgroup, tr$tip.label),
            grep("^M4-", tr$tip.label, value = TRUE))
  ladderize(drop.tip(tr, drop))
})

### concordance factors: per-branch table and the tree carrying the branch IDs
read_cf <- function(prefix) {
  list(stat   = read_tsv(paste0(prefix, ".cf.stat"), comment = "#", show_col_types = FALSE,
                         col_types = cols(.default = "c")),
       branch = read.tree(paste0(prefix, ".cf.branch")))
}

scf <- read_cf(f_scf)
gcf <- map(f_gcf, read_cf)

message("trees: ", paste(names(trees), collapse = ", "))
message("gCF: ", paste(names(gcf), collapse = ", "))
message("clusters: ", paste(sprintf("K%d n=%d", k_focus,
                                    map_int(k_focus, ~ sum(clusters$K == .x))), collapse = "  "))



# checks -----------------------------------------------------------------------

## check if the metadata hold project, country and site
if (!all(c("sample", "project", "country", "site") %in% names(meta))) {
  stop("sample_metadata_ALL_FIXED.tsv lacks sample, project, country or site")
}

## check if the cluster tables cover the same samples at every K
n_per_k <- clusters %>% count(K)

if (n_distinct(n_per_k$n) > 1) {
  warning("the cluster tables cover different samples at different K")
}

## check if the clustered samples are in the trees
for (nm in names(trees)) {
  missing <- setdiff(unique(clusters$sample), trees[[nm]]$tip.label)
  if (length(missing) > 0) {
    message(nm, ": ", length(missing), " clustered samples not in the tree")
  }
}



# clade statistics -------------------------------------------------------------

### tips of the clade below a node
clade_tips <- function(tr, node) {
  if (node <= Ntip(tr)) return(tr$tip.label[node])
  extract.clade(tr, node)$tip.label
}

### tips in plotting order (to count the runs of a cluster)
plot_order <- function(tr) {
  e <- reorder(tr, "cladewise")$edge
  tips <- e[e[, 2] <= Ntip(tr), 2]
  tr$tip.label[tips]
}

### support on a node (local posterior for ASTER, bootstrap for the
### concatenation, local bootstrap for CASTER-site)
node_support <- function(tr, node) {
  if (is.null(tr$node.label)) return(NA_real_)
  lab <- tr$node.label[node - Ntip(tr)]
  if (is.na(lab) || lab == "") return(NA_real_)
  m <- str_match(lab, "pp1=(-?[0-9.]+(?:[eE][+-]?[0-9]+)?)")[, 2]
  if (!is.na(m)) return(as.numeric(m))
  suppressWarnings(as.numeric(str_extract(lab, "[0-9.]+")))
}

### every clade of a tree as a tip set (once per tree)
clade_sets <- function(tr) {
  pp   <- prop.part(tr)
  labs <- attr(pp, "labels")
  lapply(pp, function(i) labs[i])
}

### best matching clade (smallest symmetric difference: tips to move)
best_clade <- function(sets, members, assigned) {
  m <- members
  d <- vapply(sets, function(s) {
    s <- intersect(s, assigned)
    length(setdiff(m, s)) + length(setdiff(s, m))
  }, numeric(1))
  i <- which.min(d)
  s <- intersect(sets[[i]], assigned)
  list(size    = length(s),
       missing = length(setdiff(m, s)),
       extra   = length(setdiff(s, m)),
       diff    = d[i],
       tips    = s,
       missing_tips = setdiff(m, s),
       extra_tips   = setdiff(s, m))
}

### tips on either side of the root
root_children <- function(tr) {
  r    <- Ntip(tr) + 1
  kids <- tr$edge[tr$edge[, 1] == r, 2]
  lapply(kids, function(k) clade_tips(tr, k))
}

### unrooted splits as keys (a cluster that is a split but not a clade is
### divided by the root)
split_keys <- function(tr, assigned) {
  pp   <- prop.part(tr)
  labs <- attr(pp, "labels")
  keep <- intersect(labs, assigned)
  vapply(pp, function(i) {
    a <- sort(intersect(labs[i], keep))
    b <- sort(setdiff(keep, a))
    if (length(a) < 1 || length(b) < 1) return(NA_character_)
    ka <- paste(a, collapse = "|")
    kb <- paste(b, collapse = "|")
    if (length(a) < length(b)) ka else if (length(b) < length(a)) kb else min(ka, kb)
  }, character(1))
}

### split key of a tip set (smaller side)
set_key <- function(x, keep) {
  a <- sort(intersect(x, keep))
  b <- sort(setdiff(keep, a))
  ka <- paste(a, collapse = "|")
  kb <- paste(b, collapse = "|")
  if (length(a) < length(b)) ka else if (length(b) < length(a)) kb else min(ka, kb)
}

### concordance per bipartition (branch IDs are the node labels of the .cf.branch tree), keyed on
### the tips both trees hold (IQ-TREE drops identical sequences from the sCF tree)
cf_keys <- function(cf, tips) {
  ctr    <- cf$branch
  common <- intersect(tips, ctr$tip.label)
  inner  <- Ntip(ctr) + seq_len(ctr$Nnode)
  tbl <- tibble(ID = ctr$node.label, node = inner) %>%
    filter(ID != "") %>%
    mutate(key = map_chr(node, ~ set_key(clade_tips(ctr, .x), common))) %>%
    select(ID, key) %>%
    left_join(cf$stat, by = "ID")
  list(tbl = tbl, common = common)
}

### concordance row of the branch above a clade
cf_row <- function(cf, below) {
  if (is.null(cf)) return(tibble())
  key <- set_key(below, cf$common)
  filter(cf$tbl, key == !!key)
}


### clade statistics of one cluster (intruders from other clusters and
### unassigned tips counted apart; the verdict rests on the clustered tips)
cluster_stats <- function(tr, members, assigned, sets, skeys, cf_s = NULL, cf_g = NULL) {

  members  <- intersect(members, tr$tip.label)
  assigned <- intersect(assigned, tr$tip.label)

  if (length(members) < 2) {
    return(tibble(n_tips = length(members), monophyletic = NA,
                  mrca_size = NA_integer_, intruders_clustered = NA_integer_,
                  intruders_unassigned = NA_integer_, purity = NA_real_,
                  purity_all = NA_real_, n_runs = NA_integer_,
                  best_clade_size = NA_integer_, tips_missing = NA_integer_,
                  tips_extra = NA_integer_, tips_to_move = NA_integer_,
                  unrooted_split = NA, support = NA_real_,
                  sCF = NA_real_, gCF = NA_real_, gN = NA_real_,
                  intruder_tips = NA_character_, missing_tips = NA_character_,
                  extra_tips = NA_character_, best_tips = NA_character_))
  }

  bc <- best_clade(sets, members, assigned)

  nd    <- getMRCA(tr, members)
  below <- clade_tips(tr, nd)
  intr  <- setdiff(below, members)
  intr_cl <- intersect(intr, assigned)
  intr_un <- setdiff(intr, assigned)

  ### runs over the clustered tips only
  ord  <- plot_order(tr)
  ord  <- ord[ord %in% assigned]
  inc  <- ord %in% members
  runs <- sum(inc & !c(FALSE, head(inc, -1)))

  ### site and gene concordance of the stem branch (monophyletic clusters only)
  mono  <- length(intr_cl) == 0
  s_row <- if (mono) cf_row(cf_s, below) else tibble()
  g_row <- if (mono) cf_row(cf_g, below) else tibble()

  tibble(n_tips               = length(members),
         monophyletic         = length(intr_cl) == 0,
         mrca_size            = length(below),
         intruders_clustered  = length(intr_cl),
         intruders_unassigned = length(intr_un),
         purity               = length(members) / (length(members) + length(intr_cl)),
         purity_all           = length(members) / length(below),
         n_runs               = runs,
         best_clade_size      = bc$size,
         tips_missing         = bc$missing,
         tips_extra           = bc$extra,
         tips_to_move         = bc$diff,
         unrooted_split       = set_key(members, assigned) %in% skeys,
         support              = node_support(tr, nd),
         sCF                  = if (nrow(s_row) == 1) as.numeric(s_row$sCF) else NA_real_,
         gCF                  = if (nrow(g_row) == 1) as.numeric(g_row$gCF) else NA_real_,
         gN                   = if (nrow(g_row) == 1) as.numeric(g_row$gN) else NA_real_,
         intruder_tips        = paste(head(sort(intr_cl), 20), collapse = ","),
         missing_tips         = paste(sort(bc$missing_tips), collapse = ","),
         extra_tips           = paste(sort(bc$extra_tips), collapse = ","),
         best_tips            = paste(sort(bc$tips), collapse = ","))
}


grid <- expand_grid(tree = names(trees), K = k_focus) %>%
  mutate(clusters = map2(tree, K, function(tn, k) {
    cl    <- filter(clusters, K == k)
    tr    <- trees[[tn]]
    ### clade sets and split keys once per tree and K
    sets  <- clade_sets(tr)
    skeys <- split_keys(tr, cl$sample)
    cf_s  <- if (tn == "concat") cf_keys(scf, tr$tip.label) else NULL
    cf_g  <- cf_keys(gcf[[tn]], tr$tip.label)
    map_dfr(sort(unique(cl$cluster)), function(cc) {
      members <- cl$sample[cl$cluster == cc]
      bind_cols(tibble(cluster = cc),
                cluster_stats(tr, members, cl$sample, sets, skeys, cf_s, cf_g))
    })
  })) %>%
  unnest(clusters) %>%
  mutate(status = case_when(monophyletic          ~ "monophyletic",
                            unrooted_split        ~ "split by the root",
                            purity >= purity_near ~ "nearly monophyletic",
                            TRUE                  ~ "not monophyletic"),
         status = factor(status, levels = names(col_status)))

write_tsv(select(grid, -intruder_tips, -best_tips),
          file.path(out_root, "cluster_monophyly.tsv"))



# the tips that do not fit -----------------------------------------------------

### cluster members outside the best matching clade and clade members outside
### the cluster
misplaced <- grid %>%
  select(tree, K, cluster, missing_tips, extra_tips) %>%
  pivot_longer(c(missing_tips, extra_tips), names_to = "type", values_to = "tip") %>%
  mutate(type = str_remove(type, "_tips")) %>%
  separate_rows(tip, sep = ",") %>%
  filter(!is.na(tip), tip != "") %>%
  left_join(select(meta, tip = sample, project, country, site), by = "tip")

write_tsv(misplaced, file.path(out_root, "cluster_misplaced_tips.tsv"))

cat("\n=== the tips that do not fit their best matching clade ===\n")
misplaced %>%
  count(tree, K, cluster, type, project) %>%
  print(n = Inf)

cat("\n=== clusters off by at most five tips, named ===\n")
grid %>%
  filter(tips_to_move > 0, tips_to_move <= 5) %>%
  select(tree, K, cluster, tips_to_move, missing_tips, extra_tips) %>%
  print(n = Inf)



# placement of the assemblies --------------------------------------------------

### smallest clade holding an assembly and at least one clustered sample
place_tip <- function(tr, tip, cl) {
  root <- Ntip(tr) + 1
  node <- tr$edge[tr$edge[, 2] == match(tip, tr$tip.label), 1]
  hit  <- intersect(clade_tips(tr, node), names(cl))
  while (length(hit) == 0 && node != root) {
    node <- tr$edge[tr$edge[, 2] == node, 1]
    hit  <- intersect(clade_tips(tr, node), names(cl))
  }
  n_cl <- sort(table(cl[hit]), decreasing = TRUE)
  tibble(cluster     = names(n_cl)[1],
         share       = as.numeric(n_cl[1]) / length(hit),
         n_clustered = length(hit),
         n_clade     = length(clade_tips(tr, node)),
         clusters    = paste(sprintf("%s:%d", names(n_cl), n_cl), collapse = ","),
         support     = node_support(tr, node))
}

### one row per tree, K and assembly
placement <- expand_grid(tree = names(trees), K = k_focus) %>%
  mutate(rows = map2(tree, K, function(tn, k) {
    tr      <- trees[[tn]]
    cl      <- filter(clusters, K == k) %>% select(sample, cluster) %>% deframe()
    genomes <- grep(genome_pattern, tr$tip.label, value = TRUE)
    map_dfr(genomes, ~ bind_cols(tibble(genome = .x), place_tip(tr, .x, cl)))
  })) %>%
  unnest(rows)

write_tsv(placement, file.path(out_root, "cluster_genome_placement.tsv"))

cat("\n=== placement of the assemblies (cluster, share if the clade is mixed) ===\n")
placement %>%
  mutate(cell = if_else(share < 1, sprintf("%s (%.2f)", cluster, share), cluster)) %>%
  select(genome, K, tree, cell) %>%
  pivot_wider(names_from = tree, values_from = cell) %>%
  arrange(K, genome) %>%
  print(n = Inf)



# the root ---------------------------------------------------------------------

### cluster members on each side of the root
root_tbl <- expand_grid(tree = names(trees), K = k_focus) %>%
  mutate(rows = map2(tree, K, function(tn, k) {
    tr   <- trees[[tn]]
    kids <- root_children(tr)
    cl   <- filter(clusters, K == k)
    map_dfr(sort(unique(cl$cluster)), function(cc) {
      m  <- intersect(cl$sample[cl$cluster == cc], tr$tip.label)
      a  <- intersect(m, kids[[1]])
      b  <- intersect(m, kids[[2]])
      mn <- if (length(a) <= length(b)) a else b
      tibble(cluster       = cc,
             n_side_1      = length(a),
             n_side_2      = length(b),
             straddles     = length(a) > 0 && length(b) > 0,
             n_minority    = min(length(a), length(b)),
             minority_tips = paste(sort(mn), collapse = ","))
    })
  })) %>%
  unnest(rows)

write_tsv(root_tbl, file.path(out_root, "cluster_root_partition.tsv"))

cat("\n=== root children ===\n")
for (tn in names(trees)) {
  kids <- root_children(trees[[tn]])
  cat(sprintf("  %-7s %d and %d tips\n", tn, length(kids[[1]]), length(kids[[2]])))
}

cat("\n=== clusters the root divides ===\n")
root_tbl %>%
  filter(straddles) %>%
  select(tree, K, cluster, n_side_1, n_side_2, n_minority) %>%
  print(n = Inf)



# project composition ----------------------------------------------------------

### project composition of the best matching clades (project confound)
clade_members <- grid %>%
  select(tree, K, cluster, best_tips) %>%
  separate_rows(best_tips, sep = ",") %>%
  filter(!is.na(best_tips), best_tips != "") %>%
  left_join(select(meta, best_tips = sample, project, country, site),
            by = "best_tips") %>%
  rename(tip = best_tips)

clade_project <- clade_members %>%
  count(tree, K, cluster, project) %>%
  pivot_wider(names_from = project, values_from = n, values_fill = 0)

write_tsv(clade_project, file.path(out_root, "clade_project.tsv"))

### country composition of the same clades (MAGIC is Swiss only)
clade_country <- clade_members %>%
  count(tree, K, cluster, country) %>%
  pivot_wider(names_from = country, values_from = n, values_fill = 0)

write_tsv(clade_country, file.path(out_root, "clade_country.tsv"))

### project composition of the clusters
cluster_project <- clusters %>%
  left_join(select(meta, sample, project), by = "sample") %>%
  count(K, cluster, project) %>%
  pivot_wider(names_from = project, values_from = n, values_fill = 0)

write_tsv(cluster_project, file.path(out_root, "cluster_project.tsv"))

### Swiss test (MAGIC against Swiss and other PhytOakmeter samples)
ch_test <- clade_members %>%
  mutate(origin = case_when(project == "magic"                     ~ "MAGIC (CH)",
                            country == "CH"                        ~ "PhytOakmeter CH",
                            TRUE                                   ~ "PhytOakmeter other")) %>%
  count(tree, K, cluster, origin) %>%
  pivot_wider(names_from = origin, values_from = n, values_fill = 0)

write_tsv(ch_test, file.path(out_root, "clade_swiss_test.tsv"))

cat("\n=== project composition of the best matching clades ===\n")
print(clade_project, n = Inf)

cat("\n=== the Swiss test, MAGIC against PhytOakmeter from Switzerland ===\n")
print(ch_test, n = Inf)

cat("\n=== where the Swiss PhytOakmeter samples sit ===\n")
clade_members %>%
  filter(project == "phytoakmeter", country == "CH") %>%
  count(tree, K, cluster, site) %>%
  print(n = Inf)



# summary ----------------------------------------------------------------------

summary_k <- grid %>%
  group_by(tree, K) %>%
  summarise(clusters          = n(),
            monophyletic      = sum(status == "monophyletic"),
            root_split        = sum(status == "split by the root"),
            nearly            = sum(status == "nearly monophyletic"),
            broken            = sum(status == "not monophyletic"),
            median_to_move    = median(tips_to_move, na.rm = TRUE),
            max_to_move       = max(tips_to_move, na.rm = TRUE),
            median_purity     = median(purity, na.rm = TRUE),
            intruders_clustered  = sum(intruders_clustered, na.rm = TRUE),
            intruders_unassigned = sum(intruders_unassigned, na.rm = TRUE),
            .groups = "drop")

write_tsv(summary_k, file.path(out_root, "monophyly_summary.tsv"))

cat("\n=== clusters against the topologies ===\n")
print(summary_k, n = Inf)

cat("\n=== every cluster ===\n")
grid %>%
  select(tree, K, cluster, n_tips, mrca_size, intruders_clustered,
         purity, n_runs, best_clade_size, tips_missing, tips_extra,
         tips_to_move, unrooted_split, support, sCF, gCF, status) %>%
  arrange(tree, K, cluster) %>%
  print(n = Inf)

cat("\n=== the", n_show, "tips that break the most clusters ===\n")
grid %>%
  filter(intruders_clustered > 0, !is.na(intruder_tips)) %>%
  separate_rows(intruder_tips, sep = ",") %>%
  filter(intruder_tips != "") %>%
  count(tip = intruder_tips, sort = TRUE) %>%
  left_join(select(meta, tip = sample, project, country), by = "tip") %>%
  print(n = n_show)



# figures ----------------------------------------------------------------------

### purity of the smallest containing clade (clustered tips only)
p_purity <- ggplot(grid, aes(x = cluster, y = purity, fill = status)) +
  geom_col() +
  geom_hline(yintercept = purity_near, colour = "grey40", linetype = "dashed",
             linewidth = 0.4) +
  facet_grid(tree ~ K, scales = "free_x", space = "free_x",
             labeller = labeller(K = function(x) paste("K =", x))) +
  scale_fill_manual(values = col_status, name = NULL, drop = FALSE) +
  labs(x = "Cluster", y = "Purity of the smallest containing clade")

### unassigned tips (assemblies, reference) inside each clade
p_unassigned <- ggplot(grid, aes(x = cluster, y = intruders_unassigned)) +
  geom_col(fill = "grey50") +
  facet_grid(tree ~ K, scales = "free_x", space = "free_x",
             labeller = labeller(K = function(x) paste("K =", x))) +
  labs(x = "Cluster", y = "N unassigned tips inside the clade")

save_fig("cluster_unassigned.pdf", p_unassigned, 7.1, 5.2)


### tips to move for an exact clade
p_move <- ggplot(grid, aes(x = cluster, y = tips_to_move, fill = status)) +
  geom_col() +
  facet_grid(tree ~ K, scales = "free_x", space = "free_x",
             labeller = labeller(K = function(x) paste("K =", x))) +
  scale_fill_manual(values = col_status, name = NULL, drop = FALSE) +
  labs(x = "Cluster", y = "N tips to move for an exact clade")

save_fig("cluster_tips_to_move.pdf", p_move, 7.1, 5.2)

save_fig("cluster_purity.pdf", p_purity, 7.1, 5.2)


### separate positions of each cluster in the tree
p_runs <- ggplot(grid, aes(x = cluster, y = n_runs, fill = status)) +
  geom_col() +
  facet_grid(tree ~ K, scales = "free_x", space = "free_x",
             labeller = labeller(K = function(x) paste("K =", x))) +
  scale_fill_manual(values = col_status, name = NULL, drop = FALSE) +
  labs(x = "Cluster", y = "N separate positions in the tree")

save_fig("cluster_runs.pdf", p_runs, 7.1, 5.2)


### stem support of the monophyletic clusters: support/sCF/gCF for the concatenation,
### support/gCF for the others (UFBoot, local posterior, CASTER local bootstrap); intruders otherwise
fmt_num <- function(x, d) if_else(is.na(x), "-", formatC(x, format = "f", digits = d))

p_stem <- grid %>%
  mutate(label = case_when(
    !monophyletic                ~ sprintf("+%d", intruders_clustered),
    tree == "concat"             ~ sprintf("%s/%s/%s", fmt_num(support, 0), fmt_num(sCF, 0), fmt_num(gCF, 0)),
    str_starts(tree, "aster")    ~ sprintf("%s/%s", fmt_num(support, 2), fmt_num(gCF, 0)),
    TRUE                         ~ sprintf("%s/%s", fmt_num(support, 0), fmt_num(gCF, 0))),
    cluster = fct_rev(factor(cluster, levels = str_sort(unique(cluster), numeric = TRUE)))) %>%
  ggplot(aes(x = tree, y = cluster, fill = status)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = label), size = 2.4) +
  facet_grid(K ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(K = function(x) paste("K =", x))) +
  scale_fill_manual(values = col_status, name = NULL, drop = FALSE) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = "Cluster") +
  theme(axis.line = element_blank(), axis.ticks = element_blank())

save_fig("cluster_stem_support.pdf", p_stem, 7.1, 6)


### trees with the tips coloured by cluster, one figure per tree and K
for (nm in names(trees_fig)) {
  for (k in k_focus) {

    tr <- trees_fig[[nm]]

    tip_cl <- tibble(label = tr$tip.label) %>%
      left_join(filter(clusters, K == k) %>% select(label = sample, cluster),
                by = "label")

    ### black tip labels, points filled by cluster
    p <- ggtree(tr, linewidth = 0.2) %<+% tip_cl +
      geom_tiplab(size = 1.1, align = TRUE, linesize = 0.1, linetype = "dotted") +
      geom_tippoint(aes(fill = cluster), shape = 21, colour = "black",
                    stroke = 0.2, size = 0.9, na.rm = TRUE) +
      scale_fill_manual(values = pal_cluster(n_distinct(na.omit(tip_cl$cluster))),
                        name = "Cluster", na.value = "grey75") +
      hexpand(0.06)

    save_fig(sprintf("tree_%s_K%d.pdf", nm, k), p, 24, 60)
  }
}



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
cat("  cluster_monophyly.tsv         one row per tree, K and cluster (incl. sCF and gCF of the stem)\n")
cat("  cluster_misplaced_tips.tsv    the tips missing from or extra to the best clade\n")
cat("  cluster_genome_placement.tsv  cluster of the smallest clade around each assembly\n")
cat("  cluster_root_partition.tsv    how each cluster divides over the two root children\n")
cat("  clade_project.tsv             project composition of each best matching clade\n")
cat("  clade_country.tsv             country composition of each best matching clade\n")
cat("  clade_swiss_test.tsv          MAGIC against Swiss and non-Swiss PhytOakmeter\n")
cat("  cluster_project.tsv           project composition of each cluster\n")
cat("  monophyly_summary.tsv         clusters per tree and K by status\n")
cat("  cluster_purity.pdf            purity of the smallest containing clade\n")
cat("  cluster_unassigned.pdf        genomes and reference inside each clade\n")
cat("  cluster_tips_to_move.pdf      tips to move for the cluster to be a clade\n")
cat("  cluster_runs.pdf              how many places in the tree each cluster occupies\n")
cat("  cluster_stem_support.pdf      support, sCF and gCF of the stem of each monophyletic cluster\n")
cat("  tree_<tree>_K<k>.pdf          the trees as 052, 080, 081 and 082 root them, outgroup\n")
cat("                                included, tips coloured by cluster\n")
