#!/usr/bin/env Rscript

get_script_path <- function() {
  hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (!length(hit)) return(NULL)
  sub("^--file=", "", hit[[1]])
}

script_path <- get_script_path()
fork_root <- if (is.null(script_path)) normalizePath(".", winslash = "/") else
  normalizePath(dirname(script_path), winslash = "/", mustWork = TRUE)
aging_root <- normalizePath(file.path(fork_root, "..", ".."), winslash = "/", mustWork = TRUE)
clean_root <- normalizePath(
  Sys.getenv("AGING_CLEAN_ROOT", unset = "D:/Github/aging-clean"),
  winslash = "/", mustWork = TRUE
)

source(file.path(clean_root, "R", "setup_env.R"))
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

out_dir <- file.path(fork_root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  corrected_dds = file.path(clean_root, "rnaseq_publication_pipeline", "input",
                            "dds_dual_ortholog_mm10_v102.rds"),
  corrected_annotation = file.path(clean_root,
    "rnaseq_publication_pipeline_dual_rebuild_mm10_v102", "tables", "ortholog_annotation.csv"),
  corrected_master = file.path(clean_root,
    "rnaseq_publication_pipeline_dual_rebuild_mm10_v102", "results", "tables",
    "master_table_canonical.csv"),
  legacy_results = file.path(aging_root, "rnaseq_publication_pipeline", "results", "tables",
    "young_only_dual_species_deseq2", "young_only_dual_species_deseq2_full_universe.csv"),
  pair_tsv = file.path(clean_root, "rnaseq_publication_pipeline", "input",
    "orthofinder_2species_Ensembl102_AcomysRapidRelease_Acomys_dimidiatus__v__Mus_musculus.tsv"),
  aco_counts = file.path(clean_root, "rnaseq_publication_pipeline", "input",
                         "gene_count_matrix_GTF_acomys.csv"),
  mus_counts = file.path(clean_root, "rnaseq_publication_pipeline", "input",
                         "gene_count_matrix_GTF_mus.csv"),
  sample_sheet = file.path(clean_root, "rnaseq_publication_pipeline", "input",
                           "sample_sheet_dual.csv")
)
missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing)) stop("Missing required inputs: ", paste(missing, collapse = ", "))

strip_version <- function(x) sub("[.][0-9]+$", "", as.character(x))
normalize_species <- function(x) case_when(
  as.character(x) %in% c("Acomys", "Acomys dimidiatus") ~ "Acomys dimidiatus",
  as.character(x) %in% c("Mus", "Mus musculus") ~ "Mus musculus",
  TRUE ~ as.character(x)
)
read_count_matrix <- function(path) {
  x <- read_csv(path, show_col_types = FALSE, progress = FALSE)
  ids <- strip_version(sub("[|].*$", "", x[[1]]))
  m <- as.matrix(x[, -1, drop = FALSE]); storage.mode(m) <- "numeric"
  rownames(m) <- ids
  if (anyDuplicated(rownames(m))) stop("Duplicated stable IDs in ", path)
  m
}
stars <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "***",
  ifelse(p < .01, "**", ifelse(p < .05, "*", ""))))

run_young_wald <- function(dds, analysis_label) {
  cd <- as.data.frame(SummarizedExperiment::colData(dds))
  if (!all(c("ageGroup", "species") %in% names(cd)))
    stop(analysis_label, ": ageGroup/species missing from colData")
  keep <- as.character(cd$ageGroup) == "young"
  if (sum(keep) != 6L) stop(analysis_label, ": expected six young samples")
  y <- dds[, keep]
  colData(y)$species <- factor(normalize_species(colData(y)$species),
                               levels = c("Acomys dimidiatus", "Mus musculus"))
  design(y) <- ~ species
  y <- DESeq(y, test = "Wald", quiet = TRUE)
  coef_name <- grep("^species_", resultsNames(y), value = TRUE)
  if (length(coef_name) != 1L) stop("Ambiguous species coefficient: ", paste(coef_name, collapse = ", "))
  raw <- results(y, name = coef_name)
  shr <- lfcShrink(y, coef = coef_name, type = "apeglm")
  tbl <- as.data.frame(raw) %>% rownames_to_column("gene_id") %>% as_tibble() %>%
    rename(log2FoldChange_raw = log2FoldChange, lfcSE_raw = lfcSE,
           stat_raw = stat, pvalue_raw = pvalue, padj_raw = padj) %>%
    left_join(as.data.frame(shr) %>% rownames_to_column("gene_id") %>%
      transmute(gene_id, log2FoldChange = log2FoldChange, lfcSE = lfcSE), by = "gene_id") %>%
    mutate(analysis = analysis_label,
           direction = case_when(log2FoldChange > 0 ~ "higher_in_Mus",
                                 log2FoldChange < 0 ~ "higher_in_Acomys",
                                 TRUE ~ "no_difference"),
           significant_fdr_0_05 = !is.na(padj_raw) & padj_raw < .05)
  list(dds = y, raw = raw, shrunk = shr, table = tbl, coefficient = coef_name)
}

# ---- 1. Corrected strict reciprocal 1:1 analysis --------------------------
dds_corrected <- readRDS(paths$corrected_dds)
if (nrow(dds_corrected) != 12944L)
  stop("Corrected dds has ", nrow(dds_corrected), " rows; expected 12,944")
strict_fit <- run_young_wald(dds_corrected, "mm10_v102_strict_1to1")

annot <- read_csv(paths$corrected_annotation, show_col_types = FALSE) %>%
  distinct(gene_id, .keep_all = TRUE)
strict_tbl <- strict_fit$table %>% left_join(annot, by = "gene_id") %>%
  arrange(padj_raw, desc(abs(log2FoldChange)))
write_csv(strict_tbl, file.path(out_dir, "01_corrected_strict_young_deseq2.csv"))
saveRDS(strict_fit$dds, file.path(out_dir, "01_corrected_strict_young_dds.rds"))
saveRDS(strict_fit$raw, file.path(out_dir, "01_corrected_strict_young_results_raw.rds"))
saveRDS(strict_fit$shrunk, file.path(out_dir, "01_corrected_strict_young_results_shrunk.rds"))

# ---- 2. Direct comparison to the submitted 13,101-gene analysis ----------
legacy <- read_csv(paths$legacy_results, show_col_types = FALSE, name_repair = "unique") %>%
  transmute(gene_id, legacy_acomys_gene_id = acomys_gene_id,
            legacy_symbol = symbol, legacy_lfc = log2FoldChange,
            legacy_lfc_raw = log2FoldChange_raw, legacy_stat = stat_raw,
            legacy_pvalue = pvalue_raw, legacy_padj = padj_raw,
            legacy_sig = significant_fdr_0_05, in_legacy = TRUE)
comparison <- full_join(legacy, strict_tbl %>%
  transmute(gene_id, corrected_acomys_gene_id = acomys_gene_id,
            corrected_symbol = symbol, corrected_lfc = log2FoldChange,
            corrected_lfc_raw = log2FoldChange_raw, corrected_stat = stat_raw,
            corrected_pvalue = pvalue_raw, corrected_padj = padj_raw,
            corrected_sig = significant_fdr_0_05, in_corrected = TRUE), by = "gene_id") %>%
  mutate(
    membership = case_when(in_legacy %in% TRUE & in_corrected %in% TRUE ~ "shared",
                           in_legacy %in% TRUE ~ "legacy_only",
                           TRUE ~ "corrected_only"),
    partner_changed = membership == "shared" &
      legacy_acomys_gene_id != corrected_acomys_gene_id,
    call_transition = case_when(
      membership == "corrected_only" & corrected_sig ~ "rescued_significant",
      membership == "corrected_only" ~ "rescued_not_significant",
      membership == "legacy_only" & legacy_sig ~ "legacy_significant_not_admissible",
      membership == "legacy_only" ~ "legacy_not_admissible",
      legacy_sig & corrected_sig ~ "significant_both",
      !legacy_sig & corrected_sig ~ "gained_significance",
      legacy_sig & !corrected_sig ~ "lost_significance",
      TRUE ~ "not_significant_either"),
    direction_flip = membership == "shared" & sign(legacy_lfc) != sign(corrected_lfc),
    delta_lfc = corrected_lfc - legacy_lfc,
    corrected_stars = stars(corrected_padj)
  )
write_csv(comparison, file.path(out_dir, "02_legacy_vs_corrected_gene_audit.csv"))
write_csv(comparison %>% filter(membership != "shared" | partner_changed |
  call_transition %in% c("gained_significance", "lost_significance") | direction_flip),
  file.path(out_dir, "03_changed_or_at_risk_genes.csv"))

# ---- 3. Family-aware sensitivity analysis --------------------------------
# Parse the supplied pairwise graph. Connected components are calculated from
# explicit Acomys--Mus edges, not whole OrthoFinder orthogroups, so independent
# paralog pairs in a broad orthogroup are not unnecessarily merged.
pair_raw <- read_tsv(paths$pair_tsv, show_col_types = FALSE, progress = FALSE)
edges <- pair_raw %>% transmute(
  orthogroup = Orthogroup,
  aco = str_split(as.character(Acomys_dimidiatus), ",[[:space:]]*"),
  mus = str_split(as.character(Mus_musculus), ",[[:space:]]*")
) %>% unnest_longer(aco) %>% unnest_longer(mus) %>%
  transmute(orthogroup, aco = strip_version(trimws(aco)), mus = strip_version(trimws(mus))) %>%
  filter(aco != "", mus != "") %>% distinct()

nodes <- unique(c(paste0("A:", edges$aco), paste0("M:", edges$mus)))
parent <- setNames(nodes, nodes)
find_root <- function(x) {
  while (parent[[x]] != x) x <- parent[[x]]
  x
}
union_nodes <- function(a, b) {
  ra <- find_root(a); rb <- find_root(b)
  if (ra != rb) parent[[rb]] <<- ra
}
for (i in seq_len(nrow(edges))) union_nodes(paste0("A:", edges$aco[i]), paste0("M:", edges$mus[i]))
roots <- vapply(nodes, find_root, character(1))
component_key <- match(roots, unique(roots))
node_tbl <- tibble(node = nodes, component_n = component_key,
                   species = substr(node, 1, 1), gene_id = substring(node, 3))
components <- node_tbl %>% group_by(component_n) %>% summarise(
  aco_gene_ids = paste(sort(gene_id[species == "A"]), collapse = ";"),
  mus_gene_ids = paste(sort(gene_id[species == "M"]), collapse = ";"),
  n_aco = sum(species == "A"), n_mus = sum(species == "M"), .groups = "drop") %>%
  mutate(component_id = sprintf("ORTHOCOMP_%05d", component_n),
         relationship = paste0(n_aco, ":", n_mus))

aco_counts <- read_count_matrix(paths$aco_counts)
mus_counts <- read_count_matrix(paths$mus_counts)
sample_meta <- read.csv(paths$sample_sheet, row.names = 1, check.names = FALSE)
young_samples <- rownames(sample_meta)[as.character(sample_meta$ageGroup) == "young"]
sample_meta <- sample_meta[young_samples, , drop = FALSE]
sample_meta$species <- factor(normalize_species(sample_meta$species),
                              levels = c("Acomys dimidiatus", "Mus musculus"))

sum_present <- function(ids, mat, samples) {
  ids <- ids[ids %in% rownames(mat)]
  if (!length(ids)) return(rep(NA_real_, length(samples)))
  colSums(mat[ids, samples, drop = FALSE])
}
family_counts <- matrix(NA_real_, nrow(components), length(young_samples),
                        dimnames = list(components$component_id, young_samples))
for (i in seq_len(nrow(components))) {
  a_ids <- strsplit(components$aco_gene_ids[i], ";", fixed = TRUE)[[1]]
  m_ids <- strsplit(components$mus_gene_ids[i], ";", fixed = TRUE)[[1]]
  a_samp <- rownames(sample_meta)[sample_meta$species == "Acomys dimidiatus"]
  m_samp <- rownames(sample_meta)[sample_meta$species == "Mus musculus"]
  family_counts[i, a_samp] <- sum_present(a_ids, aco_counts, a_samp)
  family_counts[i, m_samp] <- sum_present(m_ids, mus_counts, m_samp)
}
complete <- rowSums(is.na(family_counts)) == 0
family_counts <- round(family_counts[complete, , drop = FALSE])
components_tested <- components[complete, , drop = FALSE]
keep <- rowSums(family_counts >= 10) >= 3L
family_counts <- family_counts[keep, , drop = FALSE]
components_tested <- components_tested[keep, , drop = FALSE]

dds_family <- DESeqDataSetFromMatrix(family_counts, sample_meta[colnames(family_counts), , drop = FALSE],
                                     design = ~ species)
family_fit <- run_young_wald(dds_family, "mm10_v102_connected_component_sensitivity")
family_tbl <- family_fit$table %>% left_join(components_tested, by = c("gene_id" = "component_id")) %>%
  arrange(padj_raw, desc(abs(log2FoldChange)))
write_csv(family_tbl, file.path(out_dir, "04_family_component_young_deseq2.csv"))
write_csv(components_tested, file.path(out_dir, "05_family_component_membership.csv"))
saveRDS(family_fit$dds, file.path(out_dir, "04_family_component_young_dds.rds"))

# ENO1/ENO3 audit: stable IDs are authoritative; Acomys source symbols alone
# cannot distinguish the two loci safely.
eno_ids <- c("ENSADMG00000016109", "ENSADMG00000022897",
             "ENSMUSG00000060600", "ENSMUSG00000063524", "ENSMUSG00000059040")
eno_components <- components_tested %>% filter(
  str_detect(aco_gene_ids, paste(eno_ids[1:2], collapse = "|")) |
  str_detect(mus_gene_ids, paste(eno_ids[3:5], collapse = "|")))
eno_audit <- eno_components %>% left_join(family_tbl, by = c("component_id" = "gene_id"))
write_csv(eno_audit, file.path(out_dir, "06_ENO1_ENO3_component_audit.csv"))

summary_tbl <- bind_rows(
  tibble(layer = "legacy_strict", universe_n = nrow(legacy),
         significant_n = sum(legacy$legacy_sig, na.rm = TRUE)),
  tibble(layer = "corrected_strict", universe_n = nrow(strict_tbl),
         significant_n = sum(strict_tbl$significant_fdr_0_05, na.rm = TRUE)),
  tibble(layer = "family_component_sensitivity", universe_n = nrow(family_tbl),
         significant_n = sum(family_tbl$significant_fdr_0_05, na.rm = TRUE))
)
write_csv(summary_tbl, file.path(out_dir, "00_analysis_summary.csv"))
write_csv(comparison %>% count(membership, call_transition, name = "n_genes"),
          file.path(out_dir, "00_call_transition_summary.csv"))

message("Young-only orthology refresh complete: ", out_dir)
print(summary_tbl)
