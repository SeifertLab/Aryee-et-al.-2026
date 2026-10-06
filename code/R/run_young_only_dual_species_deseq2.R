#!/usr/bin/env Rscript

# Reproduce the young Acomys-versus-Mus analysis from the processed prepDE
# gene-count matrices and the OrthoFinder table.
#
# The analysis has three stages:
#   1. construct the reciprocal one-to-one cross-species count matrix;
#   2. apply the paper's count filter across all 18 samples;
#   3. run the six-sample young-only DESeq2 Wald comparison.

get_script_path <- function() {
  hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(hit)) return(NULL)
  sub("^--file=", "", hit[[1]])
}

script_path <- get_script_path()
project_root <- if (is.null(script_path)) {
  normalizePath(".", winslash = "/", mustWork = TRUE)
} else {
  normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
}

data_root <- Sys.getenv("ARYEE_DATA_ROOT", unset = file.path(project_root, "data"))
out_dir <- Sys.getenv(
  "ARYEE_DE_RESULTS_DIR",
  unset = file.path(project_root, "results", "young_only_deseq2")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  acomys_counts = Sys.getenv(
    "ARYEE_ACOMYS_COUNTS",
    unset = file.path(data_root, "processed_counts", "gene_count_matrix_GTF_acomys.csv")
  ),
  mus_counts = Sys.getenv(
    "ARYEE_MUS_COUNTS",
    unset = file.path(data_root, "processed_counts", "gene_count_matrix_GTF_mus.csv")
  ),
  sample_sheet = Sys.getenv(
    "ARYEE_SAMPLE_SHEET",
    unset = file.path(data_root, "metadata", "sample_sheet_dual.csv")
  ),
  orthofinder_pairs = Sys.getenv(
    "ARYEE_ORTHOFINDER_PAIRS",
    unset = file.path(
      data_root, "orthology",
      "Acomys_dimidiatus__v__Mus_musculus_Ensembl102.tsv"
    )
  )
)

missing_paths <- unlist(paths, use.names = FALSE)
missing_paths <- missing_paths[!file.exists(missing_paths)]
if (length(missing_paths)) {
  stop("Missing required input(s):\n", paste(missing_paths, collapse = "\n"))
}

required_packages <- c("DESeq2", "apeglm", "dplyr", "readr", "stringr", "tidyr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

EXPECTED_STRICT_PAIRS <- 16314L
EXPECTED_STRICT_PRESENT <- 16310L
EXPECTED_FILTERED_PAIRS <- 12944L
EXPECTED_YOUNG_SAMPLES <- 6L

strip_version <- function(x) sub("[.][0-9]+$", "", as.character(x))

split_ids <- function(x) {
  stringr::str_split(as.character(x), ",[[:space:]]*")
}

read_count_matrix <- function(path) {
  tbl <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  raw_labels <- as.character(tbl[[1]])
  stable_ids <- strip_version(sub("[|].*$", "", raw_labels))
  symbols <- ifelse(
    grepl("[|]", raw_labels),
    sub("^[^|]*[|]", "", raw_labels),
    stable_ids
  )
  symbols[is.na(symbols) | symbols == ""] <- stable_ids[is.na(symbols) | symbols == ""]

  mat <- as.matrix(tbl[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- stable_ids
  if (anyDuplicated(rownames(mat))) {
    stop("Duplicated stable gene IDs in ", path)
  }

  list(
    counts = mat,
    annotation = tibble(gene_id = stable_ids, symbol = symbols) %>%
      distinct(gene_id, .keep_all = TRUE)
  )
}

normalize_species <- function(x) {
  dplyr::case_when(
    as.character(x) %in% c("Acomys", "Acomys dimidiatus") ~ "Acomys dimidiatus",
    as.character(x) %in% c("Mus", "Mus musculus") ~ "Mus musculus",
    TRUE ~ as.character(x)
  )
}

# 1. Construct the strict reciprocal one-to-one ortholog map.
pair_raw <- readr::read_tsv(
  paths$orthofinder_pairs,
  show_col_types = FALSE,
  progress = FALSE
)
required_pair_cols <- c("Orthogroup", "Acomys_dimidiatus", "Mus_musculus")
if (!all(required_pair_cols %in% names(pair_raw))) {
  stop("Unexpected OrthoFinder table columns in ", paths$orthofinder_pairs)
}

pair_edges <- pair_raw %>%
  transmute(
    orthogroup = Orthogroup,
    acomys_gene_id = split_ids(Acomys_dimidiatus),
    mus_gene_id = split_ids(Mus_musculus)
  ) %>%
  unnest_longer(acomys_gene_id) %>%
  unnest_longer(mus_gene_id) %>%
  transmute(
    orthogroup,
    acomys_gene_id = strip_version(trimws(acomys_gene_id)),
    mus_gene_id = strip_version(trimws(mus_gene_id))
  ) %>%
  filter(acomys_gene_id != "", mus_gene_id != "") %>%
  distinct()

pair_support <- pair_edges %>% distinct(acomys_gene_id, mus_gene_id)
acomys_degree <- pair_support %>% count(acomys_gene_id, name = "n_mus_per_acomys")
mus_degree <- pair_support %>% count(mus_gene_id, name = "n_acomys_per_mus")

strict_pairs <- pair_support %>%
  left_join(acomys_degree, by = "acomys_gene_id") %>%
  left_join(mus_degree, by = "mus_gene_id") %>%
  filter(n_mus_per_acomys == 1L, n_acomys_per_mus == 1L) %>%
  arrange(mus_gene_id, acomys_gene_id)

if (nrow(strict_pairs) != EXPECTED_STRICT_PAIRS) {
  stop(
    "Strict reciprocal map has ", nrow(strict_pairs),
    " pairs; expected ", EXPECTED_STRICT_PAIRS
  )
}

# 2. Align the two prepDE count matrices by ortholog pair.
acomys_input <- read_count_matrix(paths$acomys_counts)
mus_input <- read_count_matrix(paths$mus_counts)

strict_presence <- strict_pairs %>%
  mutate(
    present_acomys = acomys_gene_id %in% rownames(acomys_input$counts),
    present_mus = mus_gene_id %in% rownames(mus_input$counts),
    present_both = present_acomys & present_mus
  )

strict_present <- strict_presence %>% filter(present_both)
if (nrow(strict_present) != EXPECTED_STRICT_PRESENT) {
  stop(
    "Strict map has ", nrow(strict_present),
    " pairs represented in both count matrices; expected ", EXPECTED_STRICT_PRESENT
  )
}

acomys_sub <- acomys_input$counts[strict_present$acomys_gene_id, , drop = FALSE]
mus_sub <- mus_input$counts[strict_present$mus_gene_id, , drop = FALSE]
rownames(acomys_sub) <- strict_present$mus_gene_id
rownames(mus_sub) <- strict_present$mus_gene_id
combined_counts <- cbind(acomys_sub, mus_sub)

sample_meta <- read.csv(
  paths$sample_sheet,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
missing_count_columns <- setdiff(rownames(sample_meta), colnames(combined_counts))
unexpected_count_columns <- setdiff(colnames(combined_counts), rownames(sample_meta))
if (length(missing_count_columns) || length(unexpected_count_columns)) {
  stop(
    "Count/sample mismatch. Missing count columns: ",
    paste(missing_count_columns, collapse = ", "),
    "; unexpected count columns: ",
    paste(unexpected_count_columns, collapse = ", ")
  )
}
combined_counts <- combined_counts[, rownames(sample_meta), drop = FALSE]

sample_meta$ageGroup <- factor(
  as.character(sample_meta$ageGroup),
  levels = c("embryonic", "young", "old")
)
sample_meta$species <- factor(
  normalize_species(sample_meta$species),
  levels = c("Acomys dimidiatus", "Mus musculus")
)

dds_all <- DESeqDataSetFromMatrix(
  countData = round(combined_counts),
  colData = sample_meta,
  design = ~ species + ageGroup + species:ageGroup
)
count_filter <- rowSums(counts(dds_all) >= 10) >= 3L
dds_all <- dds_all[count_filter, ]
if (nrow(dds_all) != EXPECTED_FILTERED_PAIRS) {
  stop(
    "Count-filtered object has ", nrow(dds_all),
    " rows; expected ", EXPECTED_FILTERED_PAIRS
  )
}

tested_manifest <- strict_present %>%
  filter(mus_gene_id %in% rownames(dds_all)) %>%
  mutate(
    gene_id = mus_gene_id,
    strict_1to1_ortholog = TRUE
  ) %>%
  left_join(
    mus_input$annotation %>% rename(mus_gene_id = gene_id),
    by = "mus_gene_id"
  ) %>%
  mutate(symbol = coalesce(na_if(symbol, ""), gene_id)) %>%
  select(
    gene_id, symbol, acomys_gene_id, mus_gene_id,
    strict_1to1_ortholog, n_mus_per_acomys, n_acomys_per_mus
  )

if (nrow(tested_manifest) != EXPECTED_FILTERED_PAIRS || anyDuplicated(tested_manifest$gene_id)) {
  stop("Tested ortholog manifest is incomplete or duplicated.")
}

# 3. Run the six-sample young-only species comparison.
young_keep <- as.character(colData(dds_all)$ageGroup) == "young"
if (sum(young_keep) != EXPECTED_YOUNG_SAMPLES) {
  stop("Expected six young samples; found ", sum(young_keep))
}

dds_young <- dds_all[, young_keep]
colData(dds_young)$species <- factor(
  normalize_species(colData(dds_young)$species),
  levels = c("Acomys dimidiatus", "Mus musculus")
)
design(dds_young) <- ~ species
dds_young <- DESeq(dds_young, test = "Wald", quiet = TRUE)

coefficient <- grep("^species_", resultsNames(dds_young), value = TRUE)
if (length(coefficient) != 1L) {
  stop("Expected one species coefficient; found: ", paste(coefficient, collapse = ", "))
}

raw_results <- results(dds_young, name = coefficient)
shrunk_results <- lfcShrink(dds_young, coef = coefficient, type = "apeglm")

results_tbl <- as.data.frame(raw_results) %>%
  rownames_to_column("gene_id") %>%
  as_tibble() %>%
  rename(
    log2FoldChange_raw = log2FoldChange,
    lfcSE_raw = lfcSE,
    stat_raw = stat,
    pvalue_raw = pvalue,
    padj_raw = padj
  ) %>%
  left_join(
    as.data.frame(shrunk_results) %>%
      rownames_to_column("gene_id") %>%
      transmute(
        gene_id,
        log2FoldChange = log2FoldChange,
        lfcSE = lfcSE
      ),
    by = "gene_id"
  ) %>%
  left_join(tested_manifest, by = "gene_id") %>%
  mutate(
    analysis = "young_Mus_vs_Acomys_strict_1to1",
    species_reference = "Acomys dimidiatus",
    species_comparison = "Mus musculus",
    direction = case_when(
      log2FoldChange > 0 ~ "higher_in_Mus",
      log2FoldChange < 0 ~ "higher_in_Acomys",
      TRUE ~ "no_difference"
    ),
    significant_fdr_0_05 = !is.na(padj_raw) & padj_raw < 0.05
  ) %>%
  arrange(padj_raw, desc(abs(log2FoldChange)))

normalized_counts_tbl <- as.data.frame(
  counts(dds_young, normalized = TRUE),
  check.names = FALSE
) %>%
  rownames_to_column("gene_id") %>%
  left_join(tested_manifest %>% select(gene_id, symbol), by = "gene_id") %>%
  relocate(gene_id, symbol)

aligned_raw_counts_tbl <- as.data.frame(
  counts(dds_all, normalized = FALSE),
  check.names = FALSE
) %>%
  rownames_to_column("gene_id") %>%
  left_join(tested_manifest, by = "gene_id") %>%
  relocate(
    gene_id, symbol, acomys_gene_id, mus_gene_id,
    strict_1to1_ortholog, n_mus_per_acomys, n_acomys_per_mus
  )

summary_tbl <- tibble(
  analysis = "young_Mus_vs_Acomys_strict_1to1",
  orthofinder_rows = nrow(pair_raw),
  expanded_unique_gene_pairs = nrow(pair_support),
  strict_reciprocal_pairs = nrow(strict_pairs),
  strict_pairs_present_in_both_count_matrices = nrow(strict_present),
  tested_gene_pairs = nrow(results_tbl),
  count_filter = "count >= 10 in at least 3 of 18 samples",
  young_samples = ncol(dds_young),
  young_acomys_samples = sum(as.character(colData(dds_young)$species) == "Acomys dimidiatus"),
  young_mus_samples = sum(as.character(colData(dds_young)$species) == "Mus musculus"),
  coefficient = coefficient,
  positive_log2fc = "higher in Mus musculus",
  significant_fdr_0_05 = sum(results_tbl$significant_fdr_0_05, na.rm = TRUE)
)

readr::write_csv(summary_tbl, file.path(out_dir, "analysis_summary.csv"))
readr::write_csv(strict_presence, file.path(out_dir, "ortholog_pair_presence_audit.csv"))
readr::write_csv(tested_manifest, file.path(out_dir, "tested_ortholog_manifest.csv"))
readr::write_tsv(
  aligned_raw_counts_tbl,
  file.path(out_dir, "aligned_filtered_raw_counts.tsv.gz")
)
readr::write_csv(results_tbl, file.path(out_dir, "young_only_dual_species_deseq2.csv"))
readr::write_tsv(
  normalized_counts_tbl,
  file.path(out_dir, "young_only_normalized_counts.tsv.gz")
)
saveRDS(dds_young, file.path(out_dir, "young_only_dds.rds"))
saveRDS(raw_results, file.path(out_dir, "young_only_results_raw.rds"))
saveRDS(shrunk_results, file.path(out_dir, "young_only_results_shrunk.rds"))

message("Young-only dual-species analysis complete: ", out_dir)
print(summary_tbl, width = Inf)
