## ============================================================================
## 103_mating_type.R
## Mating type ratio per site and lineage, agreement with the PCR, clonal groups
##
## Files needed:
##   mating_type.tsv                call per sample
##   mat_variants.tsv               variable sites per sample
##   mat_callable.tsv               callable sites per sample and idiomorph
##   mat_primer_genotypes.tsv       genotypes at the primer sites
##                                  (all four from hpc/01_variant_calling/103_mating_type_collect.sh;
##                                   copy them from ${ROOT}/01_variant_calling/results/103_results_mating_type_collect
##                                   to local/01_variant_calling/data/103_data_mating_type_collect)
##   mat_primer_sites.tsv           primer sites in the reference, from hpc/01_variant_calling/101_mat_idiomorphs.sh
##                                  (copy ${ROOT}/01_variant_calling/results/101_results_mat_idiomorphs/mat_primer_sites.tsv
##                                   to local/01_variant_calling/data/103_data_mating_type_collect)
##   sample_metadata_edit.tsv       site, project, MLG, genome-wide depth and PCR call,
##                                  from hpc/01_variant_calling/100_metadata.sh
##                                  (copy ${ROOT}/01_variant_calling/results/100_results_metadata/sample_metadata.tsv
##                                   to local/01_variant_calling/data/103_data_mating_type_collect
##                                   as sample_metadata_edit.tsv and add the PCR mating type column mat by hand)
##   clone_groups.tsv               multilocus genotypes, from local/01_variant_calling/090_clone_id.R
##   lineage_labels.tsv             lineages at K = 4, 5 and 6, from local/03_delimitation/001_lineage_labels.R
##   samples_removed_missingness.tsv
##                                  samples removed by missingness, from local/01_variant_calling/080_filtering.R
## ============================================================================

library(tidyverse)
library(ggpubr)



# files ------------------------------------------------------------------------

data_root <- "01_variant_calling/data/103_data_mating_type_collect"
out_root  <- "01_variant_calling/results/103_results_mating_type_collect"

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

f_mat           <- file.path(data_root, "mating_type.tsv")
f_variants      <- file.path(data_root, "mat_variants.tsv")
f_callable      <- file.path(data_root, "mat_callable.tsv")
f_primer_gt     <- file.path(data_root, "mat_primer_genotypes.tsv")
f_primers       <- file.path(data_root, "mat_primer_sites.tsv")
f_meta          <- file.path(data_root, "sample_metadata_edit.tsv")
f_clones        <- "01_variant_calling/results/090_results_clone_id/clone_groups.tsv"
f_lineage       <- "03_delimitation/results/001_results_lineage_labels/lineage_labels.tsv"
f_removed       <- "01_variant_calling/results/080_results_filtering/samples_removed_missingness.tsv"



# parameters -------------------------------------------------------------------

### min. genotypes per site and lineage for the test
site_min <- 8

### significance level (BH-adjusted)
alpha <- 0.05

### breadth thresholds of the call (as in hpc/01_variant_calling/102_mating_type.sh)
breadth_present <- 0.75
breadth_absent  <- 0.40

### min. depth of a substitution
depth_min <- 3

### primer bases counted as the 3' end
prime3 <- 5

### read file names that differ from the metadata
sample_rename <- c("CH-GOX-W2" = "CH-GOX-W2C", "FR-BOX-E23C" = "FR-BOX-E3C")

### lineage levels (the first for the site tables)
lineage_k <- c("K4", "K5", "K6")

### genome-wide depth column of the metadata
col_depth <- "mean_seq_depth"



# style ------------------------------------------------------------------------

### call levels and labels
lab_mat <- c("MAT1"      = "MAT1-1",
             "MAT2"      = "MAT1-2",
             "MAT12"     = "Both idiomorphs",
             "ambiguous" = "Ambiguous",
             "none"      = "No idiomorph")

### project labels
lab_set <- c(all = "All samples", magic = "MAGIC", phytoakmeter = "PhytOakmeter")

### lineage colours
pal_unit <- c(L1  = "#E69F00",
              L2  = "#CFE8FF", L2a = "#CFE8FF", L2b = "#1E5AA8",
              L3  = "#D9F0D3", L3a = "#D9F0D3", L3b = "#007225",
              L4  = "#CC79A7")

theme_set(theme_pubr(base_size = 10, legend = "bottom"))
theme_update(plot.tag   = element_text(size = 12, face = "bold"),
             axis.title = element_text(size = 10))

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(out_root, name), plot, width = width, height = height,
         device = cairo_pdf, limitsize = FALSE)
}



# functions --------------------------------------------------------------------

### read names to metadata names
fix_names <- function(d) {
  mutate(d, sample = recode(sample, !!!sample_rename))
}

### binomial test against 1:1 with 95% CI
binom_ci <- function(d) {
  d %>%
    mutate(bt      = map2(mat1, n, ~ binom.test(.x, .y, p = 0.5)),
           lo      = map_dbl(bt, ~ .x$conf.int[1]),
           hi      = map_dbl(bt, ~ .x$conf.int[2]),
           p_value = map_dbl(bt, "p.value")) %>%
    select(-bt)
}

### reverse complement
revcomp <- function(x) {
  chartr("ACGT", "TGCA", stringi::stri_reverse(toupper(x)))
}

### allele of a haploid genotype
gt_allele <- function(ref, alt, gt) {
  alleles <- str_split(paste(ref, alt, sep = ","), ",")
  map2_chr(alleles, suppressWarnings(as.integer(gt)),
           ~ if (is.na(.y)) NA_character_ else .x[.y + 1])
}

### rows of the called idiomorph
own_idiomorph <- function(d) {
  filter(d, (mat_called == "MAT1" & idiomorph == "MAT1_1") |
           (mat_called == "MAT2" & idiomorph == "MAT1_2"))
}



# read -------------------------------------------------------------------------

mat <- read_tsv(f_mat, show_col_types = FALSE) %>%
  fix_names()

meta <- read_tsv(f_meta, show_col_types = FALSE)

clones <- read_tsv(f_clones, show_col_types = FALSE) %>%
  mutate(mlg = as.character(mlg))

removed <- read_tsv(f_removed, show_col_types = FALSE)

### lineages, one column per K (clustered samples; clonemates follow their MLG)
lineages <- read_tsv(f_lineage, show_col_types = FALSE) %>%
  filter(clustered) %>%
  transmute(sample, K4 = k4, K5 = k5, K6 = k6) %>%
  select(sample, all_of(lineage_k))

### locus tables (optional)
variants <- if (file.exists(f_variants)) {
  read_tsv(f_variants, show_col_types = FALSE,
           col_types = cols(gt = col_character(), ad = col_character())) %>%
    fix_names()
} else {
  message("no variant table, the locus section is skipped")
  NULL
}

callable <- if (file.exists(f_callable)) {
  read_tsv(f_callable, show_col_types = FALSE) %>%
    fix_names()
} else {
  NULL
}

### primer sites (primer in reference orientation)
primers <- if (file.exists(f_primers)) {
  read_tsv(f_primers, show_col_types = FALSE) %>%
    transmute(idiomorph = seqID, primer = patternName, strand, start, end,
              primer_plus = if_else(strand == "+", toupper(pattern), revcomp(pattern)))
} else {
  NULL
}

primer_gt <- if (file.exists(f_primer_gt)) {
  read_tsv(f_primer_gt, show_col_types = FALSE, col_types = cols(gt = col_character())) %>%
    fix_names()
} else {
  NULL
}

### analysis set (ingroup passing missingness, single samples as own MLG,
### relative depth)
dat <- mat %>%
  left_join(select(meta, sample, project, country, site, in_ingroup,
                   mlg, clone_status, mat_pcr = mat,
                   depth_genome = all_of(col_depth)),
            by = "sample") %>%
  filter(in_ingroup, project != "outgroup", !sample %in% removed$sample) %>%
  mutate(mat_called = factor(mat_called, levels = names(lab_mat)),
         mlg = coalesce(as.character(mlg), sample),
         rel_mat1 = depth_mat1 / depth_genome,
         rel_mat2 = depth_mat2 / depth_genome)

### lineage per MLG (clonemates inherit it)
lineage_mlg <- lineages %>%
  inner_join(select(dat, sample, mlg), by = "sample") %>%
  select(mlg, all_of(lineage_k)) %>%
  distinct()

dat <- left_join(dat, lineage_mlg, by = "mlg", relationship = "many-to-one")



# checks -----------------------------------------------------------------------

## check if every called sample is in the metadata
missing <- anti_join(mat, meta, by = "sample")

if (nrow(missing) > 0) {
  warning(nrow(missing), " samples called but not in the metadata: ",
          paste(head(missing$sample, 10), collapse = ", "))
}

## check if every metadata sample has a call
no_call <- anti_join(meta, mat, by = "sample")

if (nrow(no_call) > 0) {
  warning(nrow(no_call), " samples in the metadata without a call: ",
          paste(head(no_call$sample, 10), collapse = ", "))
}

### called samples removed by missingness
cat("\n=== samples removed by missingness ===\n")
cat(sum(mat$sample %in% removed$sample), "of", nrow(removed), "\n")

## check if every sample has a lineage
no_lineage <- filter(dat, is.na(.data[[lineage_k[1]]]))

if (nrow(no_lineage) > 0) {
  warning(nrow(no_lineage), " samples without a lineage: ",
          paste(head(no_lineage$sample, 10), collapse = ", "))
}

## check if every lineage has a colour
no_colour <- setdiff(unlist(select(lineages, -sample)), names(pal_unit))

if (length(no_colour) > 0) {
  warning("lineages without a colour in pal_unit: ", paste(no_colour, collapse = ", "))
}

## check if every sample has a genome-wide depth
if (any(is.na(dat$depth_genome))) {
  warning(sum(is.na(dat$depth_genome)), " samples without a genome-wide depth")
}

### calls per category
cat("\n=== calls ===\n")
count(dat, mat_called) %>% print(n = Inf)

### samples between the breadth thresholds
between <- filter(dat, (breadth_mat1 > breadth_absent & breadth_mat1 < breadth_present) |
                    (breadth_mat2 > breadth_absent & breadth_mat2 < breadth_present))

cat("\n=== samples between the breadth thresholds ===\n")
cat(nrow(between), "of", nrow(dat), "\n")

if (nrow(between) > 0) {
  between %>%
    select(sample, breadth_mat1, breadth_mat2, mat_called) %>%
    arrange(desc(pmax(breadth_mat1, breadth_mat2))) %>%
    print(n = 25)
}



# genome-based mating type -----------------------------------------------------

### metadata with the read-based call next to the PCR call (all samples)
meta_all <- meta %>%
  left_join(select(mat, sample, mat_genome = mat_called), by = "sample") %>%
  relocate(mat_genome, .after = mat)

### analysis set (ingroup passing missingness)
meta_analysis <- filter(meta_all, in_ingroup, project != "outgroup",
                        !sample %in% removed$sample)

### one sample per MLG (unique genotypes and clone representatives)
meta_declone <- filter(meta_analysis, clone_status %in% c("Unique", "Representative"))

write_tsv(meta_all,      file.path(out_root, "sample_metadata_ALL_FIXED.tsv"))
write_tsv(meta_analysis, file.path(out_root, "sample_metadata_ANALYSIS_FIXED.tsv"))
write_tsv(meta_declone,  file.path(out_root, "sample_metadata_ANALYSIS_DECLONED_FIXED.tsv"))

## check if the declone set matches the clustering samples
if (!setequal(meta_declone$sample, lineages$sample)) {
  warning("declone set and clustering samples differ: ",
          paste(head(symdiff(meta_declone$sample, lineages$sample), 10), collapse = ", "))
}

cat("\n=== genome-based mating type per metadata file ===\n")
list(ALL = meta_all, ANALYSIS = meta_analysis, ANALYSIS_DECLONED = meta_declone) %>%
  imap_dfr(~ count(.x, mat_genome) %>% mutate(file = .y)) %>%
  pivot_wider(names_from = file, values_from = n, values_fill = 0) %>%
  print(n = Inf)



# agreement with the PCR -------------------------------------------------------

### read call against the PCR call (MAT1 or MAT2 by both)
agree <- dat %>%
  filter(mat_called %in% c("MAT1", "MAT2"),
         mat_pcr %in% c("MAT1", "MAT2")) %>%
  mutate(agreement = if_else(as.character(mat_called) == mat_pcr, "agree", "differ"))

write_tsv(agree, file.path(out_root, "mating_type_vs_pcr.tsv"))

cat("\n=== reads against the PCR ===\n")
count(agree, agreement) %>% mutate(f = n / sum(n)) %>% print(n = Inf)

cat("\n=== where they differ, by site ===\n")
agree %>%
  filter(agreement == "differ") %>%
  count(site, sort = TRUE) %>%
  print(n = Inf)

### read calls against every PCR category
pcr_tab <- dat %>%
  mutate(mat_pcr = factor(coalesce(lab_mat[mat_pcr], "No PCR call"),
                          levels = c(lab_mat[c("MAT1", "MAT2", "MAT12")], "No PCR call")),
         mat_called = fct_drop(factor(lab_mat[as.character(mat_called)], levels = lab_mat))) %>%
  count(mat_pcr, mat_called) %>%
  complete(mat_pcr, mat_called, fill = list(n = 0)) %>%
  group_by(mat_pcr) %>%
  mutate(f = n / sum(n)) %>%
  ungroup()

write_tsv(pivot_wider(pcr_tab, id_cols = mat_pcr, names_from = mat_called, values_from = n),
          file.path(out_root, "mating_type_vs_pcr_categories.tsv"))

cat("\n=== read calls against every PCR category ===\n")
pcr_tab %>%
  pivot_wider(id_cols = mat_pcr, names_from = mat_called, values_from = n) %>%
  print(n = Inf)



# relative depth ---------------------------------------------------------------

### idiomorph depth relative to genome-wide depth
rel_depth <- dat %>%
  select(sample, project, site, mat_called, depth_genome,
         depth_mat1, depth_mat2, rel_mat1, rel_mat2)

write_tsv(rel_depth, file.path(out_root, "mating_type_relative_depth.tsv"))

cat("\n=== relative depth by call ===\n")
rel_depth %>%
  group_by(mat_called) %>%
  summarise(n = n(),
            depth_genome = median(depth_genome, na.rm = TRUE),
            rel_mat1     = median(rel_mat1, na.rm = TRUE),
            rel_mat2     = median(rel_mat2, na.rm = TRUE),
            .groups = "drop") %>%
  print(n = Inf)



# ratio per site ---------------------------------------------------------------

### typed samples
typed <- filter(dat, mat_called %in% c("MAT1", "MAT2"))

### one sample per MLG and site (MLGs with both types left out)
by_genotype <- typed %>%
  group_by(mlg) %>%
  filter(n_distinct(mat_called) == 1) %>%
  ungroup() %>%
  filter(!is.na(site)) %>%
  distinct(site, mlg, .keep_all = TRUE)

### MAT1 and MAT2 counts per group
ratio_tbl <- function(d, label, by = "site") {
  d %>%
    mutate(mat_called = droplevels(mat_called)) %>%
    count(across(all_of(by)), mat_called, .drop = FALSE) %>%
    pivot_wider(names_from = mat_called, values_from = n, values_fill = 0) %>%
    transmute(across(all_of(by)),
              set = label,
              mat1 = MAT1, mat2 = MAT2,
              n = MAT1 + MAT2,
              f_mat1 = MAT1 / (MAT1 + MAT2))
}

### ratio per site, all samples and one per genotype
site_ratio <- bind_rows(ratio_tbl(filter(typed, !is.na(site)), "all samples"),
                        ratio_tbl(by_genotype, "one per genotype")) %>%
  binom_ci()

### ratio per site and lineage
site_lineage <- by_genotype %>%
  mutate(lineage = .data[[lineage_k[1]]]) %>%
  filter(!is.na(lineage)) %>%
  ratio_tbl("one per genotype", by = c("site", "lineage")) %>%
  mutate(both = mat1 > 0 & mat2 > 0)

### binomial test against 1:1 per site and lineage
site_test <- site_lineage %>%
  filter(n >= site_min) %>%
  binom_ci() %>%
  mutate(p_adj  = p.adjust(p_value, method = "BH"),
         skewed = p_adj < alpha)

write_tsv(site_ratio,   file.path(out_root, "mating_type_by_site.tsv"))
write_tsv(site_lineage, file.path(out_root, "mating_type_by_site_lineage.tsv"))
write_tsv(site_test,    file.path(out_root, "mating_type_site_test.tsv"))

cat("\n=== sites holding both types within a lineage (two or more genotypes) ===\n")
site_lineage %>%
  filter(n >= 2) %>%
  group_by(lineage) %>%
  summarise(sites = n(), both = sum(both), .groups = "drop") %>%
  print(n = Inf)

cat("\n=== sites departing from 1:1 over genotypes, per lineage ===\n")
site_test %>%
  arrange(p_adj) %>%
  select(site, lineage, mat1, mat2, n, f_mat1, p_value, p_adj, skewed) %>%
  print(n = Inf)



# clonal groups ----------------------------------------------------------------

### mating types within each MLG
clone_mat <- typed %>%
  group_by(mlg) %>%
  summarise(n_typed = n(),
            n_mat1 = sum(mat_called == "MAT1"),
            n_mat2 = sum(mat_called == "MAT2"),
            members = paste(sample, collapse = ","),
            .groups = "drop") %>%
  filter(n_typed >= 2) %>%
  mutate(consistent = n_mat1 == 0 | n_mat2 == 0) %>%
  left_join(select(clones, mlg, n_members, max_n_diff, min_n_compared),
            by = "mlg")

write_tsv(clone_mat, file.path(out_root, "mating_type_by_clone.tsv"))

cat("\n=== clonal groups with more than one mating type ===\n")
clone_mat %>%
  filter(!consistent) %>%
  arrange(desc(max_n_diff)) %>%
  select(mlg, n_members, n_mat1, n_mat2, max_n_diff, members) %>%
  print(n = Inf)

cat("\n=== consistency against how far apart the members are ===\n")
clone_mat %>%
  group_by(consistent) %>%
  summarise(groups = n(),
            median_max_diff = median(max_n_diff, na.rm = TRUE),
            max_max_diff    = max(max_n_diff, na.rm = TRUE),
            .groups = "drop") %>%
  print(n = Inf)



# the locus --------------------------------------------------------------------

### substitution rate per sample (filled below)
sub_sample <- tibble(sample = character(), sub_rate = numeric())

if (!is.null(variants)) {

  ### substitutions in the called idiomorph (SNPs at callable depth)
  var_own <- variants %>%
    inner_join(select(dat, sample, site, country, mat_called), by = "sample") %>%
    own_idiomorph() %>%
    filter(!gt %in% c("0", "."), dp >= depth_min) %>%
    mutate(alt = gt_allele(ref, alt, gt)) %>%
    filter(nchar(ref) == 1, nchar(alt) == 1)

  ### substitutions under a primer
  if (!is.null(primers) && !is.null(primer_gt)) {

    var_own <- var_own %>%
      left_join(select(primers, idiomorph, primer, start, end),
                by = join_by(idiomorph, between(pos, start, end))) %>%
      mutate(region = if_else(is.na(primer), "idiomorph", "primer")) %>%
      select(-start, -end)

    cat("\n=== substitutions by position in the locus ===\n")
    var_own %>%
      group_by(idiomorph, region) %>%
      summarise(sites = n_distinct(pos), samples = n_distinct(sample),
                .groups = "drop") %>%
      print(n = Inf)

    ### sample bases against the primer sequence, with the PCR call
    primer_mm <- primer_gt %>%
      inner_join(select(dat, sample, site, mat_called, mat_pcr), by = "sample") %>%
      own_idiomorph() %>%
      filter(gt != ".", dp >= depth_min) %>%
      mutate(base = gt_allele(ref, alt, gt),
             base = if_else(nchar(base) == 1, toupper(base), NA_character_)) %>%
      filter(!is.na(base)) %>%
      distinct(sample, idiomorph, pos, .keep_all = TRUE) %>%
      inner_join(primers, by = join_by(idiomorph, between(pos, start, end))) %>%
      mutate(primer_base = str_sub(primer_plus, pos - start + 1, pos - start + 1),
             from_3prime = if_else(strand == "+", end - pos, pos - start),
             mismatch    = base != primer_base) %>%
      group_by(sample, site, mat_called, mat_pcr, primer) %>%
      summarise(n_sites   = n(),
                mm        = sum(mismatch),
                mm_3prime = sum(mismatch & from_3prime < prime3),
                .groups = "drop") %>%
      mutate(pcr = case_when(!mat_pcr %in% c("MAT1", "MAT2")    ~ "no PCR call",
                             mat_pcr == as.character(mat_called) ~ "agree",
                             TRUE                                ~ "differ"))

    write_tsv(primer_mm, file.path(out_root, "primer_mismatches.tsv"))

    cat("\n=== samples with a mismatch under a primer, by site ===\n")
    primer_mm %>%
      filter(mm > 0) %>%
      count(site, primer, name = "n_samples") %>%
      arrange(desc(n_samples)) %>% print(n = Inf)

    cat("\n=== 3' primer mismatches against the PCR ===\n")
    primer_mm %>%
      group_by(sample, pcr) %>%
      summarise(mm_3prime = sum(mm_3prime), .groups = "drop") %>%
      count(pcr, mismatch_3prime = mm_3prime > 0) %>%
      print(n = Inf)
  }

  ### substitution spectrum over unique sites
  spectrum <- var_own %>%
    distinct(idiomorph, pos, ref, alt) %>%
    mutate(change = paste0(ref, ">", alt),
           class  = if_else(change %in% c("A>G", "G>A", "C>T", "T>C"),
                            "transition", "transversion")) %>%
    count(idiomorph, change, class, name = "n") %>%
    group_by(idiomorph) %>%
    mutate(f = n / sum(n)) %>%
    ungroup()

  ### transitions and transversions per idiomorph
  ti_tv <- spectrum %>%
    group_by(idiomorph) %>%
    summarise(n_variants    = sum(n),
              transitions   = sum(n[class == "transition"]),
              transversions = sum(n[class == "transversion"]),
              .groups = "drop") %>%
    mutate(ti_tv = transitions / transversions)

  write_tsv(spectrum, file.path(out_root, "substitution_spectrum.tsv"))
  write_tsv(ti_tv,    file.path(out_root, "substitution_summary.tsv"))

  cat("\n=== substitution spectrum per idiomorph ===\n")
  print(ti_tv, n = Inf)

  ### substitutions per callable site
  if (!is.null(callable)) {

    sub_sample <- callable %>%
      inner_join(select(dat, sample, mat_called), by = "sample") %>%
      own_idiomorph() %>%
      left_join(count(var_own, sample, name = "n_sub"), by = "sample") %>%
      mutate(n_sub = replace_na(n_sub, 0L),
             sub_rate = n_sub / n_callable)

    write_tsv(sub_sample, file.path(out_root, "substitutions_by_sample.tsv"))
  }
}



# ratio per lineage ------------------------------------------------------------

### one sample per MLG (MLGs with both types left out)
mlg_typed <- typed %>%
  group_by(mlg) %>%
  filter(n_distinct(mat_called) == 1) %>%
  ungroup()

### ratio and binomial test per lineage, per set and K
lineage_ratio <- bind_rows(mutate(mlg_typed, set = "all"),
                           mutate(mlg_typed, set = project)) %>%
  distinct(set, mlg, .keep_all = TRUE) %>%
  pivot_longer(all_of(lineage_k), names_to = "k", values_to = "lineage") %>%
  filter(!is.na(lineage)) %>%
  mutate(mat_called = droplevels(mat_called)) %>%
  count(set, k, lineage, mat_called, .drop = FALSE) %>%
  pivot_wider(names_from = mat_called, values_from = n, values_fill = 0) %>%
  transmute(set, k, lineage,
            mat1 = MAT1, mat2 = MAT2,
            n = MAT1 + MAT2,
            f_mat1 = MAT1 / n) %>%
  filter(n > 0) %>%
  binom_ci() %>%
  group_by(set, k) %>%
  mutate(p_adj  = p.adjust(p_value, method = "BH"),
         skewed = p_adj < alpha) %>%
  ungroup()

### ratio differing among lineages (Fisher, simulated p-value)
set.seed(42)

lineage_test <- lineage_ratio %>%
  group_by(set, k) %>%
  filter(n() > 1) %>%
  summarise(lineages  = n(),
            genotypes = sum(n),
            p_value   = fisher.test(cbind(mat1, mat2), simulate.p.value = TRUE,
                                    B = 1e5)$p.value,
            .groups = "drop")

write_tsv(lineage_ratio, file.path(out_root, "mating_type_by_lineage.tsv"))
write_tsv(lineage_test,  file.path(out_root, "mating_type_lineage_test.tsv"))

cat("\n=== ratio per lineage over genotypes ===\n")
lineage_ratio %>%
  select(set, k, lineage, mat1, mat2, n, f_mat1, p_value, p_adj, skewed) %>%
  print(n = Inf)

cat("\n=== ratio differing among lineages (Fisher) ===\n")
print(lineage_test, n = Inf)

### call outcome, coverage and divergence by lineage
lineage_qc <- dat %>%
  mutate(breadth_own = case_when(mat_called == "MAT1" ~ breadth_mat1,
                                 mat_called == "MAT2" ~ breadth_mat2),
         rel_own     = case_when(mat_called == "MAT1" ~ rel_mat1,
                                 mat_called == "MAT2" ~ rel_mat2)) %>%
  left_join(select(sub_sample, sample, sub_rate), by = "sample") %>%
  pivot_longer(all_of(lineage_k), names_to = "k", values_to = "lineage") %>%
  filter(!is.na(lineage)) %>%
  group_by(k, lineage) %>%
  summarise(n_samples     = n(),
            f_typed       = mean(mat_called %in% c("MAT1", "MAT2")),
            f_mat12       = mean(mat_called == "MAT12"),
            f_ambiguous   = mean(mat_called == "ambiguous"),
            f_none        = mean(mat_called == "none"),
            breadth_own   = median(breadth_own, na.rm = TRUE),
            rel_depth_own = median(rel_own, na.rm = TRUE),
            sub_rate      = median(sub_rate, na.rm = TRUE),
            .groups = "drop")

write_tsv(lineage_qc, file.path(out_root, "mating_type_lineage_qc.tsv"))

cat("\n=== call outcome, coverage and divergence by lineage ===\n")
print(lineage_qc, n = Inf, width = Inf)



# figures ----------------------------------------------------------------------

### ratio per lineage, one sample per MLG
lin_dat <- lineage_ratio %>%
  filter(set == "all") %>%
  mutate(k = str_replace(k, "^K", "K = "))

p_lineage <- ggplot(lin_dat, aes(x = lineage, y = f_mat1)) +
  geom_hline(aes(yintercept = 0.5, linetype = "1:1 ratio"),
             colour = "grey50", linewidth = 0.3) +
  geom_linerange(aes(ymin = lo, ymax = hi), colour = "grey30", linewidth = 0.5) +
  geom_point(aes(fill = lineage), shape = 21, size = 2.4, colour = "black",
             stroke = 0.3) +
  geom_text(aes(y = 1.06, label = paste0("n = ", n)), size = 2.4, colour = "black") +
  facet_grid(. ~ k, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal_unit, guide = "none") +
  scale_linetype_manual(values = c("1:1 ratio" = "dashed"), name = NULL) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Lineage", y = "Proportion MAT1-1 (95% CI)")

save_fig("mating_type_by_lineage.pdf", p_lineage, 7.1, 3.5)

### ratio per site, one sample per MLG
site_dat <- site_ratio %>%
  filter(set == "one per genotype") %>%
  left_join(distinct(dat, site, .keep_all = TRUE) %>% select(site, project),
            by = "site") %>%
  mutate(project = recode(project, !!!lab_set),
         site    = fct_rev(factor(site)))

p_site <- ggplot(site_dat, aes(x = f_mat1, y = site)) +
  geom_vline(aes(xintercept = 0.5, linetype = "1:1 ratio"),
             colour = "grey50", linewidth = 0.3) +
  geom_linerange(aes(xmin = lo, xmax = hi), colour = "grey30", linewidth = 0.5) +
  geom_point(shape = 21, size = 1.8, fill = "grey10", colour = "black") +
  geom_text(aes(x = 1.04, label = paste0("n = ", n)), hjust = 0, size = 2.4,
            colour = "black") +
  facet_grid(project ~ ., scales = "free_y", space = "free_y") +
  scale_linetype_manual(values = c("1:1 ratio" = "dashed"), name = NULL) +
  scale_x_continuous(limits = c(0, 1.16), breaks = seq(0, 1, 0.25)) +
  labs(x = "Proportion MAT1-1 (95% CI)", y = "Site") +
  theme(axis.text.y = element_text(size = 7))

save_fig("mating_type_by_site.pdf", p_site, 7.1, 1.5 + 0.16 * nrow(site_dat))

### read calls against the PCR calls
p_pcr <- ggplot(pcr_tab, aes(x = mat_called, y = mat_pcr, fill = f)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = n, colour = f > 0.5), size = 3) +
  scale_fill_gradient(low = "grey95", high = "#4575B4", limits = c(0, 1),
                      name = "Proportion per PCR call") +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "grey20"), guide = "none") +
  scale_y_discrete(limits = rev) +
  labs(x = "Mating type from reads", y = "Mating type from PCR") +
  theme(axis.line = element_blank(), axis.ticks = element_blank())

save_fig("mating_type_vs_pcr.pdf", p_pcr, 7.1, 4)



# output -----------------------------------------------------------------------

cat("\n\nOutput written to", out_root, "\n")
### copy the three sample_metadata_*_FIXED.tsv to ${ROOT}/03_delimitation/supplementary
cat("  sample_metadata_ALL_FIXED.tsv                metadata with mat_genome, all samples\n")
cat("  sample_metadata_ANALYSIS_FIXED.tsv           the same, ingroup passing missingness\n")
cat("  sample_metadata_ANALYSIS_DECLONED_FIXED.tsv  the same, one sample per MLG\n")
cat("  mating_type_by_site.tsv           ratio, CI and test per site, both sets\n")
cat("  mating_type_by_site_lineage.tsv   counts per site and lineage over genotypes\n")
cat("  mating_type_site_test.tsv         binomial test against 1:1, per site and lineage\n")
cat("  mating_type_by_lineage.tsv        ratio, CI and test per lineage, per set and K\n")
cat("  mating_type_lineage_test.tsv      Fisher test of the ratio among lineages\n")
cat("  mating_type_lineage_qc.tsv        call outcome, coverage, divergence by lineage\n")
cat("  mating_type_relative_depth.tsv    idiomorph depth over genome-wide depth\n")
cat("  mating_type_by_clone.tsv          mating types per clonal group\n")
cat("  mating_type_vs_pcr.tsv            the read call against the PCR call\n")
cat("  mating_type_vs_pcr_categories.tsv read calls against every PCR category\n")
cat("  primer_mismatches.tsv             sample bases against the primers, with the PCR\n")
cat("  substitutions_by_sample.tsv       substitutions per callable site\n")
cat("  substitution_spectrum.tsv         the changes seen, per idiomorph (unique sites)\n")
cat("  substitution_summary.tsv          counts and the transition to transversion ratio\n")
cat("  mating_type_by_lineage.pdf        ratio per lineage with CI, K = 4 to 6\n")
cat("  mating_type_by_site.pdf           ratio per site with CI\n")
cat("  mating_type_vs_pcr.pdf            read calls against the PCR calls\n")
