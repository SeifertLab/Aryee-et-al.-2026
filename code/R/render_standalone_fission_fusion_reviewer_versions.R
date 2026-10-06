#!/usr/bin/env Rscript

# Standalone mitochondrial fission/fusion heatmaps requested during review.
#
# This renderer uses the mm10/v102 strict reciprocal 1:1 analysis
# and the same row-Z/inferno/ComplexHeatmap style in the main heatmap figure
# It creates two transparent variants:
#   1) every curated fission/fusion gene, whether significant or not;
#   2) the significant-only baseline plus Opa1, Dnm1l, Oma1, and Gdap1 
#      as per reviewer request
# Significant genes are marked with an asterisk; reviewer-requested genes are
# not promoted to significance.

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
results_root <- Sys.getenv("ARYEE_RESULTS_ROOT", unset = file.path(project_root, "results"))
out_dir <- Sys.getenv(
  "ARYEE_SUPPLEMENT_FIGURE_OUT_DIR",
  unset = file.path(project_root, "figures", "supplement")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("ComplexHeatmap", "circlize")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Install required package(s) before running: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

paths <- list(
  panel_csv = Sys.getenv(
    "ARYEE_PANEL_CSV",
    unset = file.path(data_root, "figure_inputs", "young_mito2_redox_panel_groups.csv")
  ),
  results_csv = Sys.getenv(
    "ARYEE_DE_RESULTS_CSV",
    unset = file.path(results_root, "young_only_deseq2", "young_only_dual_species_deseq2.csv")
  ),
  normalized_counts_tsv = Sys.getenv(
    "ARYEE_NORMALIZED_COUNTS",
    unset = file.path(results_root, "young_only_deseq2", "young_only_normalized_counts.tsv.gz")
  )
)

missing_paths <- unlist(paths, use.names = FALSE)
missing_paths <- missing_paths[!file.exists(missing_paths)]
if (length(missing_paths) > 0L) {
  stop("Missing required input(s): ", paste(missing_paths, collapse = ", "))
}

if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}

font_family <- "Arial"
sample_cols <- c("SMA1", "SMA5", "SMA6", "MUSY1", "MUSY2", "MUSY3")
sample_labels <- c("Aco 1", "Aco 2", "Aco 3", "Mus 1", "Mus 2", "Mus 3")
reviewer_genes <- c("Opa1", "Dnm1l", "Oma1", "Gdap1")
block_levels <- c("Mitochondrial fission", "Mitochondrial fusion")
z_limits <- c(-2, 2)
cell_mm <- 6.2
section_gap_mm <- 5

read_csv_base <- function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

read_tsv_gz_base <- function(path) {
  utils::read.delim(gzfile(path), sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
}

make_row_z <- function(mat) {
  log_norm <- log2(mat + 1)
  row_z <- t(scale(t(log_norm)))
  row_z[is.na(row_z)] <- 0
  row_z
}

panel_tbl <- read_csv_base(paths$panel_csv)
results_tbl <- read_csv_base(paths$results_csv)
counts_tbl <- read_tsv_gz_base(paths$normalized_counts_tsv)

for (nm in c("panel_group", "panel_subgroup", "gene_symbol")) {
  panel_tbl[[nm]] <- trimws(panel_tbl[[nm]])
}
results_tbl$symbol <- trimws(results_tbl$symbol)
counts_tbl$symbol <- trimws(counts_tbl$symbol)

results_tbl <- results_tbl[!duplicated(results_tbl$symbol), , drop = FALSE]
counts_tbl <- counts_tbl[!duplicated(counts_tbl$symbol), , drop = FALSE]

missing_count_cols <- setdiff(sample_cols, names(counts_tbl))
if (length(missing_count_cols) > 0L) {
  stop("Normalized-count table is missing sample columns: ", paste(missing_count_cols, collapse = ", "))
}


dynamics_tbl <- panel_tbl[
  panel_tbl$panel_group %in% c("Mitochondria fission", "Mitochondria fusion"),
  c("panel_group", "panel_subgroup", "gene_symbol", "gene_order", "source_note"),
  drop = FALSE
]
dynamics_tbl$source_block_duplicate <- duplicated(dynamics_tbl$gene_symbol) |
  duplicated(dynamics_tbl$gene_symbol, fromLast = TRUE)
dynamics_tbl <- dynamics_tbl[
  dynamics_tbl$gene_symbol != "Oma1" | dynamics_tbl$panel_group == "Mitochondria fusion",
  ,
  drop = FALSE
]

gdap_row <- data.frame(
  panel_group = "Mitochondria fission",
  panel_subgroup = "Mitochondria fission",
  gene_symbol = "Gdap1",
  gene_order = max(suppressWarnings(as.numeric(
    dynamics_tbl$gene_order[dynamics_tbl$panel_group == "Mitochondria fission"]
  )), na.rm = TRUE) + 1,
  source_note = "Reviewer-requested addition; not in the pre-review Mito2 panel",
  source_block_duplicate = FALSE,
  stringsAsFactors = FALSE
)
dynamics_tbl <- rbind(dynamics_tbl, gdap_row)
dynamics_tbl$block_label <- ifelse(
  dynamics_tbl$panel_group == "Mitochondria fission",
  "Mitochondrial fission",
  "Mitochondrial fusion"
)
dynamics_tbl$gene_order <- suppressWarnings(as.numeric(dynamics_tbl$gene_order))
dynamics_tbl$reviewer_requested <- dynamics_tbl$gene_symbol %in% reviewer_genes
dynamics_tbl$pre_review_panel_member <- dynamics_tbl$gene_symbol != "Gdap1"

res_idx <- match(dynamics_tbl$gene_symbol, results_tbl$symbol)
dynamics_tbl$gene_id <- results_tbl$gene_id[res_idx]
dynamics_tbl$acomys_gene_id <- results_tbl$acomys_gene_id[res_idx]
dynamics_tbl$mus_gene_id <- results_tbl$mus_gene_id[res_idx]
dynamics_tbl$baseMean <- suppressWarnings(as.numeric(results_tbl$baseMean[res_idx]))
dynamics_tbl$log2FoldChange_raw <- suppressWarnings(as.numeric(results_tbl$log2FoldChange_raw[res_idx]))
dynamics_tbl$log2FoldChange_shrunken <- suppressWarnings(as.numeric(results_tbl$log2FoldChange[res_idx]))
dynamics_tbl$pvalue_raw <- suppressWarnings(as.numeric(results_tbl$pvalue_raw[res_idx]))
dynamics_tbl$padj_raw <- suppressWarnings(as.numeric(results_tbl$padj_raw[res_idx]))
dynamics_tbl$significant_fdr_0_05 <- !is.na(dynamics_tbl$padj_raw) & dynamics_tbl$padj_raw < 0.05
dynamics_tbl$direction <- results_tbl$direction[res_idx]
dynamics_tbl$strict_1to1_ortholog <- toupper(as.character(
  results_tbl$strict_1to1_ortholog[res_idx]
)) %in% "TRUE"

count_idx <- match(dynamics_tbl$gene_symbol, counts_tbl$symbol)
dynamics_tbl$in_analysis_universe <- !is.na(res_idx)
dynamics_tbl$has_normalized_expression <- !is.na(count_idx)

if (!all(dynamics_tbl$in_analysis_universe)) {
  stop("Curated genes absent from the DESeq2 universe: ", paste(
    dynamics_tbl$gene_symbol[!dynamics_tbl$in_analysis_universe], collapse = ", "
  ))
}
if (!all(dynamics_tbl$has_normalized_expression)) {
  stop("Curated genes absent from normalized-count input: ", paste(
    dynamics_tbl$gene_symbol[!dynamics_tbl$has_normalized_expression], collapse = ", "
  ))
}
if (!setequal(reviewer_genes, dynamics_tbl$gene_symbol[dynamics_tbl$reviewer_requested])) {
  stop("Reviewer-requested gene membership is incomplete.")
}
if (any(dynamics_tbl$significant_fdr_0_05[dynamics_tbl$reviewer_requested])) {
  stop("A reviewer-requested gene is now significant; update the figure wording and audit.")
}

raw_count_mat <- as.matrix(counts_tbl[count_idx, sample_cols, drop = FALSE])
storage.mode(raw_count_mat) <- "numeric"
rownames(raw_count_mat) <- dynamics_tbl$gene_symbol
row_z_mat <- make_row_z(raw_count_mat)
row_z_mat <- pmax(pmin(row_z_mat, z_limits[2]), z_limits[1])

for (sample_name in sample_cols) {
  dynamics_tbl[[paste0("normalized_", sample_name)]] <- raw_count_mat[, sample_name]
  dynamics_tbl[[paste0("row_z_", sample_name)]] <- row_z_mat[, sample_name]
}

dynamics_tbl$baseline_significant_only <- dynamics_tbl$significant_fdr_0_05
dynamics_tbl$targeted_four_added <- dynamics_tbl$significant_fdr_0_05 |
  dynamics_tbl$reviewer_requested
dynamics_tbl$all_curated_genes <- TRUE
dynamics_tbl$change_from_pre_review_significant_only <- ifelse(
  dynamics_tbl$baseline_significant_only,
  "retained_significant_baseline",
  ifelse(
    dynamics_tbl$reviewer_requested,
    "reviewer_requested_nonsignificant_addition",
    "other_curated_nonsignificant_addition_full_version_only"
  )
)
dynamics_tbl$asterisk_displayed <- dynamics_tbl$significant_fdr_0_05
dynamics_tbl$asterisk_definition <- "BH FDR < 0.05, young Mus versus Acomys"

cluster_within_blocks <- function(tbl) {
  ordered_parts <- lapply(block_levels, function(block_name) {
    part <- tbl[tbl$block_label == block_name, , drop = FALSE]
    if (nrow(part) <= 1L) return(part)
    mat <- row_z_mat[part$gene_symbol, sample_cols, drop = FALSE]
    hc <- stats::hclust(stats::dist(mat), method = "complete")
    part[hc$order, , drop = FALSE]
  })
  out <- do.call(rbind, ordered_parts)
  rownames(out) <- NULL
  out$display_order <- seq_len(nrow(out))
  out
}

z_fun <- circlize::colorRamp2(
  c(z_limits[1], 0, z_limits[2]),
  c("#2D0B59", "#B63679", "#FCA50A")
)

suppress_windows_arial_device_warning <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("font family 'Arial' not found", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

render_variant <- function(membership_col, stub, title_text) {
  display_tbl <- dynamics_tbl[dynamics_tbl[[membership_col]], , drop = FALSE]
  display_tbl <- cluster_within_blocks(display_tbl)

  expr_mat <- t(row_z_mat[display_tbl$gene_symbol, sample_cols, drop = FALSE])
  rownames(expr_mat) <- sample_labels
  colnames(expr_mat) <- display_tbl$gene_symbol

  column_split <- factor(display_tbl$block_label, levels = block_levels)
  column_labels <- paste0(
    display_tbl$gene_symbol,
    ifelse(display_tbl$significant_fdr_0_05, "*", "")
  )

  ht <- Heatmap(
    expr_mat,
    name = "Row Z Score",
    col = z_fun,
    rect_gp = gpar(col = NA),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_split = column_split,
    cluster_column_slices = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    column_gap = unit(section_gap_mm, "mm"),
    row_names_side = "left",
    row_names_gp = gpar(fontfamily = font_family, fontface = "bold", fontsize = 11),
    column_names_side = "bottom",
    column_labels = column_labels,
    column_names_gp = gpar(fontfamily = font_family, fontface = "bold.italic", fontsize = 10.5),
    column_names_rot = 45,
    column_names_max_height = unit(27, "mm"),
    row_title = NULL,
    column_title = block_levels,
    column_title_rot = 0,
    column_title_gp = gpar(fontfamily = font_family, fontface = "bold", fontsize = 11),
    show_heatmap_legend = FALSE,
    width = unit(ncol(expr_mat) * cell_mm, "mm"),
    height = unit(nrow(expr_mat) * cell_mm, "mm")
  )

  legend_obj <- Legend(
    title = "Row Z Score",
    at = c(-2, -1, 0, 1, 2),
    col_fun = z_fun,
    direction = "horizontal",
    legend_width = unit(42, "mm"),
    grid_height = unit(3.5, "mm"),
    grid_width = unit(7.5, "mm"),
    title_gp = gpar(fontfamily = font_family, fontface = "bold", fontsize = 9),
    labels_gp = gpar(fontfamily = font_family, fontsize = 8)
  )

  draw_panel <- function() {
    draw(
      ht,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "bottom",
      merge_legends = TRUE,
      padding = unit(c(7, 8, 14, 10), "mm"),
      annotation_legend_list = list(legend_obj)
    )
    grid.text(
      title_text,
      x = unit(0.5, "npc"), y = unit(1, "npc") - unit(2.2, "mm"),
      just = c("center", "top"),
      gp = gpar(fontfamily = font_family, fontface = "bold", fontsize = 12)
    )
    grid.text(
      "* BH FDR < 0.05 (young Mus vs Acomys)",
      x = unit(0.5, "npc"), y = unit(1.8, "mm"),
      just = c("center", "bottom"),
      gp = gpar(fontfamily = font_family, fontsize = 8.5)
    )
  }

  pdf_width <- max(9.2, ncol(expr_mat) * 0.255 + 2.7)
  pdf_height <- 4.9
  base_stub <- file.path(out_dir, stub)

  grDevices::svg(
    filename = paste0(base_stub, ".svg"),
    width = pdf_width,
    height = pdf_height,
    family = font_family
  )

  suppress_windows_arial_device_warning(draw_panel())
  grDevices::dev.off()

  if (capabilities("cairo")) {
    grDevices::cairo_pdf(
      filename = paste0(base_stub, ".pdf"),
      width = pdf_width,
      height = pdf_height,
      family = font_family
    )
  } else {
    grDevices::pdf(file = paste0(base_stub, ".pdf"), width = pdf_width, height = pdf_height)
  }
  suppress_windows_arial_device_warning(draw_panel())
  grDevices::dev.off()

  grDevices::png(
    filename = paste0(base_stub, ".png"),
    width = round(pdf_width * 300),
    height = round(pdf_height * 300),
    res = 300
  )
  suppress_windows_arial_device_warning(draw_panel())
  grDevices::dev.off()

  display_tbl$variant <- membership_col
  display_tbl$display_label <- column_labels
  utils::write.csv(
    display_tbl,
    paste0(base_stub, "_displayed_genes.csv"),
    row.names = FALSE
  )

  invisible(display_tbl)
}

targeted_display <- render_variant(
  "targeted_four_added",
  "standalone_fission_fusion_significant_plus_reviewer_four",
  "Mitochondrial fission and fusion - significant genes plus reviewer-requested genes"
)
full_display <- render_variant(
  "all_curated_genes",
  "standalone_fission_fusion_all_curated_genes_with_significance",
  "Mitochondrial fission and fusion - all curated genes"
)

targeted_added <- setdiff(
  targeted_display$gene_symbol,
  dynamics_tbl$gene_symbol[dynamics_tbl$baseline_significant_only]
)
if (!setequal(targeted_added, reviewer_genes)) {
  stop("Targeted variant does not add exactly the four reviewer-requested genes.")
}

audit_cols <- c(
  "block_label", "gene_symbol", "gene_id", "acomys_gene_id", "mus_gene_id",
  "pre_review_panel_member", "source_block_duplicate", "reviewer_requested",
  "in_analysis_universe", "has_normalized_expression", "strict_1to1_ortholog",
  "baseMean", "log2FoldChange_raw", "log2FoldChange_shrunken", "pvalue_raw",
  "padj_raw", "significant_fdr_0_05", "direction", "baseline_significant_only",
  "targeted_four_added", "all_curated_genes", "asterisk_displayed",
  "asterisk_definition", "change_from_pre_review_significant_only", "source_note",
  paste0("normalized_", sample_cols), paste0("row_z_", sample_cols)
)
utils::write.csv(
  dynamics_tbl[, audit_cols, drop = FALSE],
  file.path(out_dir, "standalone_fission_fusion_full_gene_audit.csv"),
  row.names = FALSE
)

summary_tbl <- data.frame(
  variant = c("pre_review_significant_only", "targeted_four_added", "all_curated_genes"),
  n_unique_genes = c(
    sum(dynamics_tbl$baseline_significant_only),
    sum(dynamics_tbl$targeted_four_added),
    sum(dynamics_tbl$all_curated_genes)
  ),
  n_significant = c(
    sum(dynamics_tbl$significant_fdr_0_05[dynamics_tbl$baseline_significant_only]),
    sum(dynamics_tbl$significant_fdr_0_05[dynamics_tbl$targeted_four_added]),
    sum(dynamics_tbl$significant_fdr_0_05)
  ),
  n_nonsignificant = c(
    0,
    sum(!dynamics_tbl$significant_fdr_0_05[dynamics_tbl$targeted_four_added]),
    sum(!dynamics_tbl$significant_fdr_0_05)
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  summary_tbl,
  file.path(out_dir, "standalone_fission_fusion_variant_summary.csv"),
  row.names = FALSE
)

reviewer_stats_lines <- vapply(reviewer_genes, function(gene_name) {
  row <- dynamics_tbl[dynamics_tbl$gene_symbol == gene_name, , drop = FALSE]
  sprintf(
    "- %s: shrunken log2FC = %.3f; BH FDR = %.4f; significant = %s.",
    gene_name,
    row$log2FoldChange_shrunken,
    row$padj_raw,
    ifelse(row$significant_fdr_0_05, "yes", "no")
  )
}, character(1))



message("Wrote standalone fission/fusion reviewer variants to: ", out_dir)
