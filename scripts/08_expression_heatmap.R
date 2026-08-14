##############################################################
# Project:
# Transcriptomic Characterization of Sepsis Using RNA-seq
#
# Script:
# 08_expression_heatmap.R
#
# Purpose:
# Reproduce the manuscript heatmap showing the top 30
# differentially expressed genes from the Cluster 2 vs
# Cluster 1 comparison across all 348 sepsis samples.
#
# Recovered manuscript workflow:
#   1. Use validated Cluster 2 vs Cluster 1 DESeq2 results
#   2. Keep genes with adjusted p < 0.05 and |log2FC| >= 1
#   3. Rank genes by:
#        |log2FC| * log10(baseMean + 1)
#   4. Select the top 30 genes
#   5. Extract VST expression values
#   6. Standardize expression within each gene (row Z-score)
#   7. Clip displayed Z-scores to [-2.5, 2.5]
#   8. Order samples by validated transcriptomic cluster
#   9. Cluster genes (rows) only; do not recluster samples
#  10. Display one continuous gene heatmap without numbered
#      row splits
#
# Inputs:
#   results/r_objects/vst_data.rds
#   results/r_objects/dds_processed.rds
#   results/tables/differential_expression/
#     Cluster_2_vs_Cluster_1_all_genes.csv
#
# Outputs:
#   results/figures/expression_heatmap/
#   results/tables/expression_heatmap/
#
# Author:
# Patrick Culver
#
# Last Updated:
# August 2026
##############################################################

rm(list = ls())
graphics.off()
set.seed(12345)

#------------------------------------------------------------
# Required packages
#------------------------------------------------------------

required_packages <- c(
  "DESeq2",
  "ComplexHeatmap",
  "circlize",
  "grid",
  "svglite"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    ". Run scripts/00_setup.R first."
  )
}

library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(grid)

#------------------------------------------------------------
# Paths
#------------------------------------------------------------

vst_file <- "results/r_objects/vst_data.rds"
dds_file <- "results/r_objects/dds_processed.rds"

de_file <- file.path(
  "results",
  "tables",
  "differential_expression",
  "Cluster_2_vs_Cluster_1_all_genes.csv"
)

figure_directory <- file.path(
  "results",
  "figures",
  "expression_heatmap"
)

table_directory <- file.path(
  "results",
  "tables",
  "expression_heatmap"
)

log_directory <- file.path(
  "results",
  "logs"
)

for (directory in c(
  figure_directory,
  table_directory,
  log_directory
)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

#------------------------------------------------------------
# Validate inputs
#------------------------------------------------------------

for (file_path in c(
  vst_file,
  dds_file,
  de_file
)) {
  if (!file.exists(file_path)) {
    stop(
      "Required input file not found:\n",
      file_path
    )
  }
}

#------------------------------------------------------------
# Load project objects
#------------------------------------------------------------

vst_data <- readRDS(vst_file)
dds <- readRDS(dds_file)

if (!is(vst_data, "DESeqTransform")) {
  stop(
    "vst_data.rds is not a DESeqTransform object."
  )
}

if (!is(dds, "DESeqDataSet")) {
  stop(
    "dds_processed.rds is not a DESeqDataSet object."
  )
}

if (!identical(
  colnames(vst_data),
  colnames(dds)
)) {
  stop(
    "Sample order differs between vst_data.rds and ",
    "dds_processed.rds."
  )
}

#------------------------------------------------------------
# Load validated Cluster 2 vs Cluster 1 DE results
#------------------------------------------------------------

deg_table <- read.csv(
  de_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_de_columns <- c(
  "Ensembl",
  "SYMBOL",
  "baseMean",
  "log2FoldChange",
  "padj"
)

missing_de_columns <- setdiff(
  required_de_columns,
  colnames(deg_table)
)

if (length(missing_de_columns) > 0) {
  stop(
    "Required DE columns are missing: ",
    paste(
      missing_de_columns,
      collapse = ", "
    )
  )
}

#------------------------------------------------------------
# Heatmap settings
#------------------------------------------------------------

adjusted_p_threshold <- 0.05
absolute_lfc_threshold <- 1
number_heatmap_genes <- 30
z_score_limit <- 2.5

expected_samples <- 348

expected_cluster_sizes <- c(
  "1" = 37,
  "2" = 157,
  "3" = 134,
  "4" = 20
)

#------------------------------------------------------------
# Prepare validated cluster metadata
#------------------------------------------------------------

metadata <- as.data.frame(
  colData(dds)
)

cluster_values <- as.character(
  metadata$cluster
)

cluster_values <- gsub(
  "^Cluster\\s*",
  "",
  cluster_values,
  ignore.case = TRUE
)

metadata$cluster <- factor(
  cluster_values,
  levels = c(
    "1",
    "2",
    "3",
    "4"
  )
)

if (any(is.na(metadata$cluster))) {
  stop(
    "One or more cluster assignments could not be ",
    "interpreted as clusters 1-4."
  )
}

if (nrow(metadata) != expected_samples) {
  stop(
    "Expected ",
    expected_samples,
    " samples but found ",
    nrow(metadata),
    "."
  )
}

observed_cluster_sizes <- table(
  metadata$cluster
)

if (!identical(
  as.integer(observed_cluster_sizes),
  as.integer(expected_cluster_sizes)
)) {
  stop(
    "Cluster sizes do not match validated values.\n",
    "Expected: ",
    paste(
      expected_cluster_sizes,
      collapse = ", "
    ),
    "\nObserved: ",
    paste(
      observed_cluster_sizes,
      collapse = ", "
    )
  )
}

#------------------------------------------------------------
# Select significant genes
#------------------------------------------------------------

heatmap_genes <- deg_table[
  !is.na(deg_table$padj) &
    deg_table$padj < adjusted_p_threshold &
    !is.na(deg_table$log2FoldChange) &
    abs(deg_table$log2FoldChange) >=
      absolute_lfc_threshold,
  ,
  drop = FALSE
]

if (nrow(heatmap_genes) < number_heatmap_genes) {
  stop(
    "Fewer than ",
    number_heatmap_genes,
    " significant genes were available."
  )
}

#------------------------------------------------------------
# Rank genes using recovered manuscript score
#
# score = |log2FC| * log10(baseMean + 1)
#------------------------------------------------------------

heatmap_genes$HeatmapScore <-
  abs(heatmap_genes$log2FoldChange) *
  log10(heatmap_genes$baseMean + 1)

heatmap_genes <- heatmap_genes[
  order(
    heatmap_genes$HeatmapScore,
    decreasing = TRUE,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

selected_genes <- head(
  heatmap_genes,
  number_heatmap_genes
)

#------------------------------------------------------------
# Validate selected gene IDs against VST matrix
#------------------------------------------------------------

vst_matrix <- assay(vst_data)

selected_ensembl <- selected_genes$Ensembl

missing_expression_genes <- setdiff(
  selected_ensembl,
  rownames(vst_matrix)
)

if (length(missing_expression_genes) > 0) {
  stop(
    "Selected DE genes were not found in the VST matrix:\n",
    paste(
      missing_expression_genes,
      collapse = "\n"
    )
  )
}

#------------------------------------------------------------
# Build expression matrix
#------------------------------------------------------------

heatmap_matrix <- vst_matrix[
  selected_ensembl,
  ,
  drop = FALSE
]

gene_labels <- selected_genes$SYMBOL

missing_symbol <- is.na(gene_labels) |
  trimws(gene_labels) == ""

gene_labels[missing_symbol] <-
  selected_ensembl[missing_symbol]

rownames(heatmap_matrix) <- make.unique(
  gene_labels
)

#------------------------------------------------------------
# Standardize each gene across samples
#------------------------------------------------------------

heatmap_matrix_scaled <- t(
  scale(
    t(heatmap_matrix)
  )
)

heatmap_matrix_scaled[
  !is.finite(heatmap_matrix_scaled)
] <- 0

# Clip extreme values only for display.
heatmap_matrix_scaled[
  heatmap_matrix_scaled > z_score_limit
] <- z_score_limit

heatmap_matrix_scaled[
  heatmap_matrix_scaled < -z_score_limit
] <- -z_score_limit

#------------------------------------------------------------
# Order samples by validated cluster
#
# Samples are not reclustered. This preserves the manuscript
# cluster structure and makes the heatmap a visualization of
# the validated cluster solution rather than a new analysis.
#------------------------------------------------------------

metadata <- metadata[
  colnames(heatmap_matrix_scaled),
  ,
  drop = FALSE
]

sample_order <- order(
  metadata$cluster
)

heatmap_matrix_scaled <- heatmap_matrix_scaled[
  ,
  sample_order,
  drop = FALSE
]

cluster_ordered <- metadata$cluster[
  sample_order
]

#------------------------------------------------------------
# Save selected-gene table
#------------------------------------------------------------

selected_gene_output <- selected_genes[
  ,
  c(
    "Ensembl",
    "SYMBOL",
    "GENENAME",
    "baseMean",
    "log2FoldChange",
    "padj",
    "HeatmapScore"
  )[
    c(
      "Ensembl",
      "SYMBOL",
      "GENENAME",
      "baseMean",
      "log2FoldChange",
      "padj",
      "HeatmapScore"
    ) %in%
      colnames(selected_genes)
  ],
  drop = FALSE
]

selected_gene_output$HeatmapRank <-
  seq_len(
    nrow(selected_gene_output)
  )

selected_gene_output <- selected_gene_output[
  ,
  c(
    "HeatmapRank",
    setdiff(
      colnames(selected_gene_output),
      "HeatmapRank"
    )
  ),
  drop = FALSE
]

write.csv(
  selected_gene_output,
  file.path(
    table_directory,
    "Cluster_2_vs_Cluster_1_top30_heatmap_genes.csv"
  ),
  row.names = FALSE
)

#------------------------------------------------------------
# Cluster annotation
#------------------------------------------------------------

cluster_colors <- c(
  "1" = "#1F77B4",
  "2" = "#D62728",
  "3" = "#2CA02C",
  "4" = "#9467BD"
)

top_annotation <- HeatmapAnnotation(
  Cluster = cluster_ordered,
  col = list(
    Cluster = cluster_colors
  ),
  annotation_name_gp = gpar(
    fontsize = 11,
    fontface = "bold"
  ),
  simple_anno_size = unit(
    5,
    "mm"
  )
)

#------------------------------------------------------------
# Expression color scale
#------------------------------------------------------------

expression_colors <- circlize::colorRamp2(
  c(
    -z_score_limit,
    0,
    z_score_limit
  ),
  c(
    "#2166AC",
    "white",
    "#B2182B"
  )
)

#------------------------------------------------------------
# Create heatmap
#------------------------------------------------------------

heatmap_object <- Heatmap(
  heatmap_matrix_scaled,

  name = "Expression\nZ-score",

  col = expression_colors,

  top_annotation = top_annotation,

  show_column_names = FALSE,

  show_row_names = TRUE,

  row_names_side = "left",

  row_names_gp = gpar(
    fontsize = 10
  ),

  column_title =
    "Top 30 Differentially Expressed Genes: Cluster 2 vs Cluster 1",

  column_title_gp = gpar(
    fontsize = 15,
    fontface = "bold"
  ),

  cluster_rows = TRUE,

  cluster_columns = FALSE,

  column_split = cluster_ordered,

  column_gap = unit(
    1.5,
    "mm"
  ),

  row_dend_width = unit(
    2,
    "cm"
  ),

  border = FALSE,

  heatmap_legend_param = list(
    title = "Expression\nZ-score",
    at = c(
      -z_score_limit,
      0,
      z_score_limit
    ),
    labels = c(
      paste0("-", z_score_limit),
      "0",
      paste0("+", z_score_limit)
    ),
    title_gp = gpar(
      fontsize = 11,
      fontface = "bold"
    ),
    labels_gp = gpar(
      fontsize = 10
    )
  )
)

#------------------------------------------------------------
# Helper for drawing the heatmap consistently
#------------------------------------------------------------

draw_heatmap <- function() {

  draw(
    heatmap_object,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    merge_legends = TRUE
  )
}

#------------------------------------------------------------
# Save PNG
#------------------------------------------------------------

png_file <- file.path(
  figure_directory,
  "Cluster_2_vs_Cluster_1_top30_heatmap.png"
)

png(
  filename = png_file,
  width = 4400,
  height = 3200,
  res = 300
)

draw_heatmap()

dev.off()

#------------------------------------------------------------
# Save PDF
#------------------------------------------------------------

pdf_file <- file.path(
  figure_directory,
  "Cluster_2_vs_Cluster_1_top30_heatmap.pdf"
)

cairo_pdf(
  filename = pdf_file,
  width = 14.67,
  height = 10.67
)

draw_heatmap()

dev.off()

#------------------------------------------------------------
# Save SVG
#------------------------------------------------------------

svg_file <- file.path(
  figure_directory,
  "Cluster_2_vs_Cluster_1_top30_heatmap.svg"
)

svglite::svglite(
  file = svg_file,
  width = 14.67,
  height = 10.67
)

draw_heatmap()

dev.off()

#------------------------------------------------------------
# Save analysis settings
#------------------------------------------------------------

settings_table <- data.frame(
  Setting = c(
    "DE comparison",
    "Adjusted p-value threshold",
    "Absolute log2 fold-change threshold",
    "Gene ranking score",
    "Number of displayed genes",
    "Expression transformation",
    "Row scaling",
    "Displayed Z-score limits",
    "Sample ordering",
    "Column clustering",
    "Row clustering",
    "Heatmap legend"
  ),
  Value = c(
    "Cluster 2 vs Cluster 1",
    adjusted_p_threshold,
    absolute_lfc_threshold,
    "|log2FC| * log10(baseMean + 1)",
    number_heatmap_genes,
    "DESeq2 variance stabilizing transformation",
    "Z-score within gene",
    paste0(
      "-",
      z_score_limit,
      " to +",
      z_score_limit
    ),
    "Validated transcriptomic cluster 1-4",
    "No",
    "Yes",
    "Expression Z-score; -2.5, 0, +2.5"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  settings_table,
  file.path(
    table_directory,
    "expression_heatmap_settings.csv"
  ),
  row.names = FALSE
)

#------------------------------------------------------------
# Reproducibility log
#------------------------------------------------------------

log_lines <- c(
  "Expression Heatmap Summary",
  "==================================================",
  paste("Date:", Sys.Date()),
  paste("R version:", R.version.string),
  paste(
    "ComplexHeatmap version:",
    as.character(
      packageVersion("ComplexHeatmap")
    )
  ),
  "",
  paste("VST input:", vst_file),
  paste("Cluster input:", dds_file),
  paste("DE input:", de_file),
  "",
  "Comparison: Cluster 2 vs Cluster 1",
  paste(
    "Significance threshold: adjusted p <",
    adjusted_p_threshold
  ),
  paste(
    "Absolute log2FC threshold:",
    absolute_lfc_threshold
  ),
  paste(
    "Significant DEGs available:",
    nrow(heatmap_genes)
  ),
  paste(
    "Genes displayed:",
    number_heatmap_genes
  ),
  "Ranking score: |log2FC| * log10(baseMean + 1)",
  "Expression: VST",
  "Row scaling: gene-wise Z-score",
  paste(
    "Displayed Z-score range:",
    paste0(
      "-",
      z_score_limit,
      " to +",
      z_score_limit
    )
  ),
  "Samples ordered by validated cluster",
  "Samples not reclustered",
  "Genes clustered",
  "",
  paste(
    "Validated cluster sizes:",
    paste(
      as.integer(observed_cluster_sizes),
      collapse = ", "
    )
  )
)

writeLines(
  log_lines,
  file.path(
    log_directory,
    "expression_heatmap_summary.txt"
  )
)

#------------------------------------------------------------
# Completion summary
#------------------------------------------------------------

cat(
  "\n",
  "============================================================\n",
  "08_expression_heatmap.R COMPLETE\n",
  "============================================================\n",
  "Comparison: Cluster 2 vs Cluster 1\n",
  "Significant DEGs available: ",
  nrow(heatmap_genes),
  "\n",
  "Top genes displayed: ",
  number_heatmap_genes,
  "\n",
  "Samples displayed: ",
  ncol(heatmap_matrix_scaled),
  "\n",
  "Cluster sizes: ",
  paste(
    as.integer(observed_cluster_sizes),
    collapse = ", "
  ),
  "\n\n",
  "Figures saved:\n",
  "  ",
  png_file,
  "\n",
  "  ",
  pdf_file,
  "\n",
  "  ",
  svg_file,
  "\n\n",
  "Selected-gene table:\n",
  "  ",
  file.path(
    table_directory,
    "Cluster_2_vs_Cluster_1_top30_heatmap_genes.csv"
  ),
  "\n",
  "============================================================\n",
  sep = ""
)
