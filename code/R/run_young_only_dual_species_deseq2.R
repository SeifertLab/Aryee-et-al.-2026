get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hit <- grep(file_arg, args, value = TRUE)
  if (length(hit) == 0) return(NULL)
  sub(file_arg, "", hit[[1]], fixed = TRUE)
}

script_path <- get_script_path()
project_root <- if (!is.null(script_path)) {
  normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(".", winslash = "/", mustWork = TRUE)
}

.libPaths(unique(c(
  "C:/Users/theaw/AppData/Local/R/win-library/4.5",
  "C:/Users/theaw/AppData/Local/R/win-library/4.4",
  .libPaths()
)))

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

shared_config_path <- file.path(project_root, "aging_project_config.R")
if (file.exists(shared_config_path)) {
  source(shared_config_path, local = environment())
}

if (!exists("get_aging_project_paths")) {
  stop("Could not load get_aging_project_paths() from aging_project_config.R")
}

cfg <- get_aging_project_paths(project_root)

paths <- list(
  dds_dual = cfg$rnaseq$dds_dual,
  ortholog_annotation = file.path(project_root, "rnaseq_publication_pipeline_dual_rebuild_clean", "tables", "ortholog_annotation.csv"),
  master_table = file.path(project_root, "rnaseq_publication_pipeline_dual_rebuild_clean", "results", "tables", "master_table_canonical.csv"),
  mito_long = file.path(project_root, "rnaseq_publication_pipeline", "input", "mitochondrial_gene_sets", "mitochondrial_gene_sets_long.csv"),
  out_dir = file.path(project_root, "rnaseq_publication_pipeline", "results", "tables", "young_only_dual_species_deseq2")
)

dir.create(paths$out_dir, recursive = TRUE, showWarnings = FALSE)

expected_dual_gene_count <- 13101L
species_levels <- c("Acomys dimidiatus", "Mus musculus")
young_label <- "young"

normalize_species <- function(x) {
  x_chr <- as.character(x)
  dplyr::case_when(
    x_chr %in% c("Acomys", "Acomys dimidiatus") ~ "Acomys dimidiatus",
    x_chr %in% c("Mus", "Mus musculus") ~ "Mus musculus",
    TRUE ~ x_chr
  )
}

drop_class_columns <- function(tbl) {
  tbl %>% select(-any_of(c("class_main", "class_sub")))
}

make_results_df <- function(res_obj, shrink_obj, annot_tbl, master_tbl, mito_tbl) {
  raw_df <- as.data.frame(res_obj) %>%
    tibble::rownames_to_column("gene_id") %>%
    as_tibble() %>%
    rename(
      log2FoldChange_raw = log2FoldChange,
      lfcSE_raw = lfcSE,
      stat_raw = stat,
      pvalue_raw = pvalue,
      padj_raw = padj
    )

  shr_df <- as.data.frame(shrink_obj) %>%
    tibble::rownames_to_column("gene_id") %>%
    as_tibble() %>%
    transmute(
      gene_id,
      log2FoldChange = log2FoldChange,
      lfcSE = lfcSE
    )

  raw_df %>%
    left_join(shr_df, by = "gene_id") %>%
    left_join(annot_tbl, by = "gene_id") %>%
    left_join(master_tbl, by = "gene_id") %>%
    left_join(mito_tbl, by = "gene_id") %>%
    mutate(
      species_reference = "Acomys dimidiatus",
      species_comparison = "Mus musculus",
      contrast = "young_Mus_vs_Acomys",
      direction = case_when(
        is.na(log2FoldChange) ~ NA_character_,
        log2FoldChange > 0 ~ "higher_in_Mus",
        log2FoldChange < 0 ~ "higher_in_Acomys",
        TRUE ~ "no_difference"
      ),
      significant_fdr_0_05 = !is.na(padj_raw) & padj_raw < 0.05,
      significant_fdr_0_05_lfc_1 = !is.na(padj_raw) & padj_raw < 0.05 & abs(log2FoldChange) >= 1
    ) %>%
    arrange(padj_raw, desc(abs(log2FoldChange)), desc(baseMean))
}

dds_dual <- readRDS(paths$dds_dual)
if (nrow(dds_dual) != expected_dual_gene_count) {
  stop(
    "dds_dual has ", nrow(dds_dual), " rows; expected ", expected_dual_gene_count,
    " in the strict ortholog universe."
  )
}

cd <- as.data.frame(SummarizedExperiment::colData(dds_dual))
if (!all(c("ageGroup", "species") %in% colnames(cd))) {
  stop("dds_dual must contain ageGroup and species in colData.")
}

cd$ageGroup <- as.character(cd$ageGroup)
cd$species <- normalize_species(cd$species)

young_keep <- cd$ageGroup == young_label
if (sum(young_keep) != 6L) {
  stop("Expected exactly 6 young samples, found ", sum(young_keep), ".")
}

dds_young <- dds_dual[, young_keep]
colData(dds_young)$ageGroup <- droplevels(factor(as.character(colData(dds_young)$ageGroup)))
colData(dds_young)$species <- factor(
  normalize_species(colData(dds_young)$species),
  levels = species_levels
)

young_cd <- as.data.frame(SummarizedExperiment::colData(dds_young)) %>%
  tibble::rownames_to_column("sample_id") %>%
  mutate(
    species_short = recode(species, "Acomys dimidiatus" = "Acomys", "Mus musculus" = "Mus")
  )

species_counts <- table(as.character(colData(dds_young)$species))
if (!identical(as.integer(species_counts[species_levels]), c(3L, 3L))) {
  stop("Expected 3 Acomys and 3 Mus young samples.")
}

design(dds_young) <- ~ species
dds_young <- DESeq(dds_young, test = "Wald")

coef_name <- resultsNames(dds_young)[grepl("^species_", resultsNames(dds_young))]
if (length(coef_name) != 1L) {
  stop("Expected exactly one species coefficient, found: ", paste(coef_name, collapse = ", "))
}

res_young_raw <- results(dds_young, name = coef_name)
res_young_shr <- lfcShrink(dds_young, coef = coef_name, type = "apeglm")

ortholog_annotation <- readr::read_csv(paths$ortholog_annotation, show_col_types = FALSE) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  select(
    gene_id,
    acomys_gene_id,
    mus_gene_id,
    symbol,
    description,
    strict_1to1_ortholog,
    n_mus_per_acomys,
    n_acomys_per_mus
  )

master_tbl <- readr::read_csv(paths$master_table, show_col_types = FALSE, name_repair = "unique") %>%
  transmute(
    gene_id,
    species_young_log2FC_full_model = species_young_log2FC_dual,
    species_young_padj_full_model = species_young_padj_dual,
    interaction_padj_full_model = interaction_padj_dual
  )

mito_tbl <- readr::read_csv(paths$mito_long, show_col_types = FALSE) %>%
  filter(in_project_universe_any %in% TRUE) %>%
  group_by(gene_id) %>%
  summarise(
    mito_focus_categories = str_c(sort(unique(focus_category)), collapse = ";"),
    mito_set_names = str_c(sort(unique(set_name)), collapse = ";"),
    mito_set_count = n_distinct(set_name),
    mito_preferred_symbol = {
      source_symbols <- sort(unique(na.omit(source_gene_symbol)))
      source_symbols <- source_symbols[source_symbols != ""]
      source_symbols <- source_symbols[!str_detect(source_symbols, "^(ENSMUSG|ENSADMG)")]
      if (length(source_symbols) == 1L) source_symbols else NA_character_
    },
    is_mito_target = TRUE,
    .groups = "drop"
  )

results_tbl <- make_results_df(
  res_obj = res_young_raw,
  shrink_obj = res_young_shr,
  annot_tbl = ortholog_annotation,
  master_tbl = master_tbl,
  mito_tbl = mito_tbl
) %>%
  mutate(
    symbol = case_when(
      is_mito_target %in% TRUE &
        !is.na(mito_preferred_symbol) &
        str_detect(coalesce(symbol, ""), "^(Gm|RIK|Rik|LOC|Fam|AU0)") ~ mito_preferred_symbol,
      TRUE ~ symbol
    ),
    is_mito_target = coalesce(is_mito_target, FALSE),
    mito_focus_categories = na_if(mito_focus_categories, ""),
    mito_set_names = na_if(mito_set_names, "")
  )

results_tbl_export <- drop_class_columns(results_tbl)

summary_tbl <- tibble(
  analysis = "young_only_dual_species_deseq2",
  gene_universe_n = nrow(results_tbl),
  sample_n = ncol(dds_young),
  sample_n_acomys = unname(species_counts["Acomys dimidiatus"]),
  sample_n_mus = unname(species_counts["Mus musculus"]),
  reference_species = "Acomys dimidiatus",
  comparison_species = "Mus musculus",
  coefficient = coef_name,
  sig_fdr_0_05_n = sum(results_tbl$significant_fdr_0_05, na.rm = TRUE),
  sig_fdr_0_05_lfc_1_n = sum(results_tbl$significant_fdr_0_05_lfc_1, na.rm = TRUE),
  mito_target_n = sum(results_tbl$is_mito_target, na.rm = TRUE),
  mito_sig_fdr_0_05_n = sum(results_tbl$is_mito_target & results_tbl$significant_fdr_0_05, na.rm = TRUE),
  mito_sig_fdr_0_05_lfc_1_n = sum(results_tbl$is_mito_target & results_tbl$significant_fdr_0_05_lfc_1, na.rm = TRUE)
)

readr::write_csv(young_cd, file.path(paths$out_dir, "young_sample_manifest.csv"))
readr::write_csv(summary_tbl, file.path(paths$out_dir, "young_only_dual_species_deseq2_summary.csv"))
readr::write_csv(results_tbl_export, file.path(paths$out_dir, "young_only_dual_species_deseq2_full_universe.csv"))
readr::write_csv(
  results_tbl_export %>% filter(is_mito_target),
  file.path(paths$out_dir, "young_only_dual_species_deseq2_mito_targets.csv")
)
saveRDS(dds_young, file.path(paths$out_dir, "young_only_dual_species_dds.rds"))
saveRDS(res_young_raw, file.path(paths$out_dir, "young_only_dual_species_results_raw.rds"))
saveRDS(res_young_shr, file.path(paths$out_dir, "young_only_dual_species_results_shrunk.rds"))

message("Young-only dual-species DESeq2 run complete.")
message("Output directory: ", paths$out_dir)
print(summary_tbl)
