get_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hit <- grep(file_arg, args, value = TRUE)
  if (length(hit) == 0) {
    return(NULL)
  }
  sub(file_arg, "", hit[[1]], fixed = TRUE)
}

script_path <- get_script_path()
project_root <- if (!is.null(script_path)) {
  normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
} else {
  normalizePath(".", winslash = "/", mustWork = TRUE)
}

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

data_root <- Sys.getenv("ARYEE_DATA_ROOT", unset = file.path(project_root, "data"))
results_root <- Sys.getenv("ARYEE_RESULTS_ROOT", unset = file.path(project_root, "results"))
out_dir <- Sys.getenv(
  "ARYEE_MAIN_FIGURE_OUT_DIR",
  unset = file.path(project_root, "figures", "main")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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

required_paths <- unlist(paths, use.names = FALSE)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop("Missing required input(s): ", paste(missing_paths, collapse = ", "))
}

if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
}

font_family <- "Arial"
sample_cols <- c("SMA1", "SMA5", "SMA6", "MUSY1", "MUSY2", "MUSY3")
sample_labels <- c("Aco 1", "Aco 2", "Aco 3", "Mus 1", "Mus 2", "Mus 3")
z_limits <- c(-2, 2)
section_gap_mm <- 4
cell_mm <- 6.2

selected_blocks <- data.frame(
  panel_group = c(
    "Glycolysis",
    "Pyruvate metabolism",
    "Krebs cycle",
    "Glutathione synthesis/maintenance",
    "Glutathione synthesis/maintenance"
  ),
  panel_subgroup = c(
    "Glycolysis",
    "Pyruvate metabolism",
    "Krebs cycle",
    "GSH synthesis/import",
    "GSH recycling/NADPH"
  ),
  block_label = c(
    "Glycolysis",
    "Pyruvate metabolism",
    "Krebs cycle",
    "Glutathione synthesis/import",
    "Glutathione / NADPH recycling"
  ),
  block_order = c(1, 2, 3, 4, 5),
  stringsAsFactors = FALSE
)

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

clamp_mat <- function(mat, limits) {
  pmax(pmin(mat, limits[2]), limits[1])
}

panel_tbl <- read_csv_base(paths$panel_csv)
results_tbl <- read_csv_base(paths$results_csv)
normalized_counts_tbl <- read_tsv_gz_base(paths$normalized_counts_tsv)

panel_tbl$gene_symbol <- trimws(panel_tbl$gene_symbol)
panel_tbl$panel_group <- trimws(panel_tbl$panel_group)
panel_tbl$panel_subgroup <- trimws(panel_tbl$panel_subgroup)
results_tbl$symbol <- trimws(results_tbl$symbol)
normalized_counts_tbl$symbol <- trimws(normalized_counts_tbl$symbol)

results_tbl <- results_tbl[!duplicated(results_tbl$symbol), , drop = FALSE]
normalized_counts_tbl <- normalized_counts_tbl[!duplicated(normalized_counts_tbl$symbol), , drop = FALSE]

norm_missing <- setdiff(sample_cols, names(normalized_counts_tbl))
if (length(norm_missing) > 0) {
  stop("Normalized-count table is missing sample columns: ", paste(norm_missing, collapse = ", "))
}

normalized_mat <- as.matrix(normalized_counts_tbl[, sample_cols, drop = FALSE])
storage.mode(normalized_mat) <- "numeric"
row_z_mat <- make_row_z(normalized_mat)
colnames(row_z_mat) <- sample_cols

expr_tbl <- cbind(
  normalized_counts_tbl["symbol"],
  as.data.frame(row_z_mat, stringsAsFactors = FALSE, check.names = FALSE)
)

selected_key <- paste(selected_blocks$panel_group, selected_blocks$panel_subgroup, sep = "||")
panel_key <- paste(panel_tbl$panel_group, panel_tbl$panel_subgroup, sep = "||")
panel_sub_tbl <- panel_tbl[panel_key %in% selected_key, , drop = FALSE]

if (nrow(panel_sub_tbl) == 0) {
  stop("No panel rows matched the selected rotated blocks.")
}

block_label_map <- stats::setNames(selected_blocks$block_label, selected_key)
block_order_map <- stats::setNames(selected_blocks$block_order, selected_key)

panel_sub_tbl$block_key <- paste(panel_sub_tbl$panel_group, panel_sub_tbl$panel_subgroup, sep = "||")
panel_sub_tbl$block_label <- unname(block_label_map[panel_sub_tbl$block_key])
panel_sub_tbl$block_order <- unname(block_order_map[panel_sub_tbl$block_key])

res_idx <- match(panel_sub_tbl$gene_symbol, results_tbl$symbol)
panel_sub_tbl$gene_id <- results_tbl$gene_id[res_idx]
panel_sub_tbl$log2FoldChange <- suppressWarnings(as.numeric(results_tbl$log2FoldChange[res_idx]))
panel_sub_tbl$padj_raw <- suppressWarnings(as.numeric(results_tbl$padj_raw[res_idx]))
panel_sub_tbl$significant_fdr_0_05 <- toupper(as.character(results_tbl$significant_fdr_0_05[res_idx])) %in% "TRUE"

expr_idx <- match(panel_sub_tbl$gene_symbol, expr_tbl$symbol)
for (sample_name in sample_cols) {
  panel_sub_tbl[[sample_name]] <- suppressWarnings(as.numeric(expr_tbl[[sample_name]][expr_idx]))
}

panel_sub_tbl$in_universe <- !is.na(panel_sub_tbl$gene_id)
panel_sub_tbl$has_expression <- !is.na(expr_idx)

display_tbl <- panel_sub_tbl[
  panel_sub_tbl$significant_fdr_0_05 & panel_sub_tbl$in_universe & panel_sub_tbl$has_expression,
  ,
  drop = FALSE
]

if (nrow(display_tbl) == 0) {
  stop("No significant genes remained in the selected rotated blocks.")
}

display_tbl <- display_tbl[order(display_tbl$block_order, display_tbl$gene_order, display_tbl$gene_symbol), , drop = FALSE]
display_tbl$block_label <- factor(display_tbl$block_label, levels = selected_blocks$block_label)

expr_mat <- t(as.matrix(display_tbl[, sample_cols, drop = FALSE]))
storage.mode(expr_mat) <- "numeric"
expr_mat <- clamp_mat(expr_mat, z_limits)
rownames(expr_mat) <- sample_labels
colnames(expr_mat) <- display_tbl$gene_symbol

column_split <- droplevels(factor(
  as.character(display_tbl$block_label), levels = selected_blocks$block_label
))
column_title_map <- c(
  "Glycolysis" = "Glycolysis",
  "Pyruvate metabolism" = "Pyruvate metabolism",
  "Krebs cycle" = "Krebs cycle",
  "Glutathione synthesis/import" = "GSH Synthesis +\nImport",
  "Glutathione / NADPH recycling" = "GSH / NADPH\nRecycling",
  "Mitochondrial fission" = "Mitochondrial fission",
  "Mitochondrial fusion" = "Mitochondrial fusion"
)
column_title_labels <- unname(column_title_map[levels(column_split)])

gene_list_export <- display_tbl[, c(
  "block_label",
  "panel_group",
  "panel_subgroup",
  "gene_symbol",
  "gene_id",
  "log2FoldChange",
  "padj_raw",
  "significant_fdr_0_05"
), drop = FALSE]

audit_tbl <- panel_sub_tbl[, c(
  "block_label",
  "panel_group",
  "panel_subgroup",
  "gene_symbol",
  "gene_id",
  "in_universe",
  "has_expression",
  "significant_fdr_0_05",
  "log2FoldChange",
  "padj_raw"
), drop = FALSE]
audit_tbl$current_rotated_member <- audit_tbl$gene_symbol %in% display_tbl$gene_symbol

z_fun <- circlize::colorRamp2(
  c(z_limits[1], 0, z_limits[2]),
  c("#2D0B59", "#B63679", "#FCA50A")
)

ht <- Heatmap(
  expr_mat,
  name = "Row Z Score",
  col = z_fun,
  rect_gp = gpar(col = NA),
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  column_split = column_split,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  column_gap = unit(section_gap_mm, "mm"),
  row_names_side = "left",
  row_names_gp = gpar(fontfamily = font_family, fontface = "bold", fontsize = 12),
  column_names_side = "bottom",
  column_names_gp = gpar(fontfamily = font_family, fontface = "italic", fontsize = 12),
  column_names_rot = 45,
  column_names_max_height = unit(28, "mm"),
  row_title = NULL,
  column_title = column_title_labels,
  column_title_rot = 0,
  column_title_gp = gpar(fontfamily = font_family, fontface = "bold", fontsize = 10),
  top_annotation = NULL,
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
    padding = unit(c(2, 8, 7, 10), "mm"),
    annotation_legend_list = list(legend_obj)
  )
}

pdf_width <- max(11.5, ncol(expr_mat) * 0.205 + 2.2)
pdf_height <- 4.9
png_width <- round(pdf_width * 300)
png_height <- round(pdf_height * 300)

base_stub <- file.path(out_dir, "rotated_sig_only_blocks_inferno_complexheatmap")

grDevices::svg(
  filename = paste0(base_stub, ".svg"),
  width = pdf_width,
  height = pdf_height,
  family = font_family
)
draw_panel()
grDevices::dev.off()

if (capabilities("cairo")) {
  grDevices::cairo_pdf(
    filename = paste0(base_stub, ".pdf"),
    width = pdf_width,
    height = pdf_height,
    family = font_family
  )
} else {
  grDevices::pdf(
    file = paste0(base_stub, ".pdf"),
    width = pdf_width,
    height = pdf_height
  )
}
draw_panel()
grDevices::dev.off()

grDevices::png(
  filename = paste0(base_stub, ".png"),
  width = png_width,
  height = png_height,
  res = 300
)
draw_panel()
grDevices::dev.off()

utils::write.csv(
  gene_list_export,
  file.path(out_dir, "rotated_sig_only_blocks_gene_list.csv"),
  row.names = FALSE
)

utils::write.csv(
  audit_tbl,
  file.path(out_dir, "rotated_sig_only_blocks_audit.csv"),
  row.names = FALSE
)

writeLines(
  c(
    "Rotated significant-only inferno panel",
    "Renderer: ComplexHeatmap",
    "Blocks included:",
    "- Glycolysis",
    "- Pyruvate metabolism",
    "- Krebs cycle",
    "- Glutathione synthesis/import",
    "- Glutathione / NADPH recycling",
    "",
    "Figure rules:",
    "- Existing panel genes only;",
    "- Only significant genes (young Mus vs Acomys, FDR < 0.05) are displayed.",
    "- Samples are on the Y axis in fixed order: Aco 1, Aco 2, Aco 3, Mus 1, Mus 2, Mus 3.",
    "- Genes are on the X axis and clustered hierarchically within each block.",
    "- No dendrogram is displayed.",
    "- Gene symbols are bold italic at 45 degrees.",
    "- Sample labels are bold Arial.",
    "- A single shared horizontal Row Z Score legend applies to all blocks."
  ),
  file.path(out_dir, "README_rotated_sig_only_blocks.txt")
)

message("Wrote ComplexHeatmap rotated inferno panel to: ", out_dir)
