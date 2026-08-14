##############################################################
# Project:
# Transcriptomic Characterization of Sepsis Using RNA-seq
#
# Script:
# 05_differential_expression.R
#
# Purpose:
# Perform pairwise differential expression analysis among all
# four validated transcriptomic clusters using DESeq2.
#
# Six pairwise comparisons are evaluated:
#   Cluster 2 vs Cluster 1
#   Cluster 3 vs Cluster 1
#   Cluster 4 vs Cluster 1
#   Cluster 3 vs Cluster 2
#   Cluster 4 vs Cluster 2
#   Cluster 4 vs Cluster 3
#
# Log2 fold changes are shrunken using apeglm.
#
# Author:
# Patrick Culver
#
# Last Updated:
# July 2026
##############################################################

rm(list = ls())
graphics.off()
set.seed(12345)

#------------------------------------------------------------
# Load required packages
#------------------------------------------------------------

required_packages <- c(
  "DESeq2",
  "apeglm",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "ggplot2"
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
library(apeglm)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)

#------------------------------------------------------------
# Analysis settings
#------------------------------------------------------------

adjusted_p_threshold <- 0.05
absolute_lfc_threshold <- 1

expected_genes <- 19208
expected_samples <- 348

expected_cluster_sizes <- c(
  "1" = 37,
  "2" = 157,
  "3" = 134,
  "4" = 20
)

comparisons <- data.frame(
  Numerator = c(
    "2",
    "3",
    "4",
    "3",
    "4",
    "4"
  ),
  Denominator = c(
    "1",
    "1",
    "1",
    "2",
    "2",
    "3"
  ),
  stringsAsFactors = FALSE
)

#------------------------------------------------------------
# Create output directories
#------------------------------------------------------------

output_directories <- c(
  "results/tables/differential_expression",
  "results/r_objects/differential_expression",
  "figures/differential_expression"
)

for (directory in output_directories) {
  
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

#------------------------------------------------------------
# Load processed DESeq2 object
#------------------------------------------------------------

dds_file <- "results/r_objects/dds_processed.rds"

if (!file.exists(dds_file)) {
  stop(
    "Processed DESeq2 object not found:\n",
    dds_file
  )
}

dds <- readRDS(dds_file)

if (!is(dds, "DESeqDataSet")) {
  stop(
    "dds_processed.rds is not a DESeqDataSet object."
  )
}

#------------------------------------------------------------
# Validate dimensions
#------------------------------------------------------------

if (nrow(dds) != expected_genes) {
  stop(
    "Unexpected number of genes. Expected ",
    expected_genes,
    " but found ",
    nrow(dds),
    "."
  )
}

if (ncol(dds) != expected_samples) {
  stop(
    "Unexpected number of samples. Expected ",
    expected_samples,
    " but found ",
    ncol(dds),
    "."
  )
}

if (!"cluster" %in% colnames(colData(dds))) {
  stop(
    "The DESeq2 object does not contain a cluster metadata column."
  )
}

#------------------------------------------------------------
# Prepare cluster variable
#------------------------------------------------------------

cluster_values <- as.character(
  colData(dds)$cluster
)

# This also handles labels such as "Cluster 1"
cluster_values <- gsub(
  pattern = "^Cluster\\s*",
  replacement = "",
  x = cluster_values,
  ignore.case = TRUE
)

dds$cluster <- factor(
  cluster_values,
  levels = c(
    "1",
    "2",
    "3",
    "4"
  )
)

if (any(is.na(dds$cluster))) {
  stop(
    "One or more samples have missing or invalid cluster assignments."
  )
}

#------------------------------------------------------------
# Validate cluster sizes
#------------------------------------------------------------

observed_cluster_sizes <- table(dds$cluster)

if (!identical(
  as.integer(observed_cluster_sizes),
  as.integer(expected_cluster_sizes)
)) {
  stop(
    "Cluster sizes do not match the validated manuscript values.\n",
    "Expected: ",
    paste(
      names(expected_cluster_sizes),
      expected_cluster_sizes,
      sep = "=",
      collapse = ", "
    ),
    "\nObserved: ",
    paste(
      names(observed_cluster_sizes),
      observed_cluster_sizes,
      sep = "=",
      collapse = ", "
    )
  )
}

#------------------------------------------------------------
# Define DESeq2 design and fit the model
#------------------------------------------------------------

design(dds) <- ~ cluster

cat(
  "\nFitting the DESeq2 negative-binomial model...\n"
)

dds <- DESeq(
  dds,
  quiet = FALSE
)

saveRDS(
  dds,
  "results/r_objects/dds_deseq_all_clusters.rds"
)

#------------------------------------------------------------
# Helper function: annotate Ensembl gene IDs
#------------------------------------------------------------

annotate_results <- function(results_table) {
  
  results_table$Ensembl <- rownames(
    results_table
  )
  
  # Remove Ensembl version suffix, when present
  results_table$EnsemblBase <- sub(
    pattern = "\\..*$",
    replacement = "",
    x = results_table$Ensembl
  )
  
  gene_annotation <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(results_table$EnsemblBase),
    keytype = "ENSEMBL",
    columns = c(
      "SYMBOL",
      "GENENAME",
      "ENTREZID"
    )
  )
  
  # Some Ensembl IDs can map to more than one annotation.
  # Retain the first mapping to prevent duplicate DE rows.
  gene_annotation <- gene_annotation[
    !duplicated(gene_annotation$ENSEMBL),
    ,
    drop = FALSE
  ]
  
  annotated_table <- merge(
    results_table,
    gene_annotation,
    by.x = "EnsemblBase",
    by.y = "ENSEMBL",
    all.x = TRUE,
    sort = FALSE
  )
  
  original_order <- match(
    results_table$EnsemblBase,
    annotated_table$EnsemblBase
  )
  
  annotated_table <- annotated_table[
    original_order,
    ,
    drop = FALSE
  ]
  
  rownames(annotated_table) <- NULL
  
  desired_columns <- c(
    "Ensembl",
    "EnsemblBase",
    "SYMBOL",
    "GENENAME",
    "ENTREZID",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "pvalue",
    "padj"
  )
  
  remaining_columns <- setdiff(
    colnames(annotated_table),
    desired_columns
  )
  
  annotated_table[
    ,
    c(
      desired_columns[
        desired_columns %in%
          colnames(annotated_table)
      ],
      remaining_columns
    ),
    drop = FALSE
  ]
}

#------------------------------------------------------------
# Helper function: significance classification
#------------------------------------------------------------

classify_significance <- function(
    padj,
    log2_fold_change
) {
  
  classification <- rep(
    "Not significant",
    length(padj)
  )
  
  classification[
    !is.na(padj) &
      padj < adjusted_p_threshold &
      log2_fold_change >= absolute_lfc_threshold
  ] <- "Upregulated"
  
  classification[
    !is.na(padj) &
      padj < adjusted_p_threshold &
      log2_fold_change <= -absolute_lfc_threshold
  ] <- "Downregulated"
  
  factor(
    classification,
    levels = c(
      "Downregulated",
      "Not significant",
      "Upregulated"
    )
  )
}

#------------------------------------------------------------
# Helper function: readable p-value labels
#------------------------------------------------------------

format_p_value <- function(p_value) {
  
  if (is.na(p_value)) {
    return(NA_character_)
  }
  
  if (p_value < 0.0001) {
    return("<0.0001")
  }
  
  formatC(
    p_value,
    format = "f",
    digits = 4
  )
}

#------------------------------------------------------------
# Figure colors
#------------------------------------------------------------

volcano_colors <- c(
  "Downregulated" = "#0072B2",
  "Not significant" = "grey75",
  "Upregulated" = "#D55E00"
)

#------------------------------------------------------------
# Storage objects
#------------------------------------------------------------

summary_list <- list()
comparison_objects <- list()

#------------------------------------------------------------
# Run all pairwise comparisons
#------------------------------------------------------------

for (comparison_index in seq_len(nrow(comparisons))) {
  
  numerator <- comparisons$Numerator[
    comparison_index
  ]
  
  denominator <- comparisons$Denominator[
    comparison_index
  ]
  
  comparison_label <- paste0(
    "Cluster_",
    numerator,
    "_vs_Cluster_",
    denominator
  )
  
  comparison_title <- paste0(
    "Cluster ",
    numerator,
    " vs Cluster ",
    denominator
  )
  
  cat(
    "\n",
    "----------------------------------------\n",
    "Analyzing ",
    comparison_title,
    "\n",
    "----------------------------------------\n",
    sep = ""
  )
  
  #----------------------------------------------------------
  # Relevel cluster factor
  #
  # apeglm requires a named model coefficient rather than a
  # general contrast. Releveling produces the required
  # numerator-versus-denominator coefficient.
  #----------------------------------------------------------
  
  dds_comparison <- dds
  
  dds_comparison$cluster <- relevel(
    dds_comparison$cluster,
    ref = denominator
  )
  
  design(dds_comparison) <- ~ cluster
  
  # Re-estimate model coefficients using the chosen reference.
  # Existing dispersion estimates are retained.
  dds_comparison <- nbinomWaldTest(
    dds_comparison,
    quiet = TRUE
  )
  
  available_coefficients <- resultsNames(
    dds_comparison
  )
  
  expected_coefficient <- paste0(
    "cluster_",
    numerator,
    "_vs_",
    denominator
  )
  
  if (!expected_coefficient %in%
      available_coefficients) {
    
    stop(
      "Expected coefficient was not found for ",
      comparison_title,
      ".\nExpected: ",
      expected_coefficient,
      "\nAvailable coefficients: ",
      paste(
        available_coefficients,
        collapse = ", "
      )
    )
  }
  
  #----------------------------------------------------------
  # Obtain unshrunken DESeq2 results
  #----------------------------------------------------------
  
  raw_results <- results(
    dds_comparison,
    name = expected_coefficient,
    alpha = adjusted_p_threshold,
    independentFiltering = TRUE
  )
  
  #----------------------------------------------------------
  # Shrink log2 fold changes using apeglm
  #----------------------------------------------------------
  
  shrunken_results <- lfcShrink(
    dds = dds_comparison,
    coef = expected_coefficient,
    type = "apeglm"
  )
  
  # Preserve p-values and adjusted p-values from the standard
  # DESeq2 results object.
  shrunken_results$pvalue <- raw_results$pvalue
  shrunken_results$padj <- raw_results$padj
  
  #----------------------------------------------------------
  # Convert and annotate results
  #----------------------------------------------------------
  
  results_table <- as.data.frame(
    shrunken_results
  )
  
  results_table <- annotate_results(
    results_table
  )
  
  results_table$Comparison <- comparison_title
  
  results_table$Direction <- classify_significance(
    padj = results_table$padj,
    log2_fold_change =
      results_table$log2FoldChange
  )
  
  # Sort by adjusted p-value, followed by absolute LFC
  results_table$AbsoluteLFC <- abs(
    results_table$log2FoldChange
  )
  
  results_table <- results_table[
    order(
      results_table$padj,
      -results_table$AbsoluteLFC,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  #----------------------------------------------------------
  # Identify significant genes
  #----------------------------------------------------------
  
  significant_results <- results_table[
    !is.na(results_table$padj) &
      results_table$padj <
      adjusted_p_threshold &
      abs(results_table$log2FoldChange) >=
      absolute_lfc_threshold,
    ,
    drop = FALSE
  ]
  
  upregulated_results <- significant_results[
    significant_results$log2FoldChange >=
      absolute_lfc_threshold,
    ,
    drop = FALSE
  ]
  
  downregulated_results <- significant_results[
    significant_results$log2FoldChange <=
      -absolute_lfc_threshold,
    ,
    drop = FALSE
  ]
  
  #----------------------------------------------------------
  # Export result tables
  #----------------------------------------------------------
  
  write.csv(
    results_table,
    file.path(
      "results/tables/differential_expression",
      paste0(
        comparison_label,
        "_all_genes.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    significant_results,
    file.path(
      "results/tables/differential_expression",
      paste0(
        comparison_label,
        "_significant_DEGs.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    upregulated_results,
    file.path(
      "results/tables/differential_expression",
      paste0(
        comparison_label,
        "_upregulated.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    downregulated_results,
    file.path(
      "results/tables/differential_expression",
      paste0(
        comparison_label,
        "_downregulated.csv"
      )
    ),
    row.names = FALSE
  )
  
  #----------------------------------------------------------
  # Save comparison-specific R objects
  #----------------------------------------------------------
  
  comparison_object <- list(
    numerator = numerator,
    denominator = denominator,
    comparison = comparison_title,
    coefficient = expected_coefficient,
    raw_results = raw_results,
    shrunken_results = shrunken_results,
    annotated_results = results_table,
    significant_results = significant_results,
    upregulated_results = upregulated_results,
    downregulated_results = downregulated_results
  )
  
  saveRDS(
    comparison_object,
    file.path(
      "results/r_objects/differential_expression",
      paste0(
        comparison_label,
        "_results.rds"
      )
    )
  )
  
  comparison_objects[[comparison_label]] <- comparison_object
  
  #----------------------------------------------------------
  # Create volcano plot data
  #----------------------------------------------------------
  
  volcano_data <- results_table
  
  volcano_data$NegativeLog10Padj <- -log10(
    volcano_data$padj
  )
  
  finite_y_values <- volcano_data$NegativeLog10Padj[
    is.finite(volcano_data$NegativeLog10Padj)
  ]
  
  if (length(finite_y_values) == 0) {
    
    y_axis_cap <- 10
    
  } else {
    
    y_axis_cap <- max(
      finite_y_values,
      na.rm = TRUE
    ) * 1.05
  }
  
  volcano_data$NegativeLog10Padj[
    is.infinite(
      volcano_data$NegativeLog10Padj
    )
  ] <- y_axis_cap
  
  # Select the most statistically significant genes for labels
  label_candidates <- volcano_data[
    volcano_data$Direction !=
      "Not significant" &
      !is.na(volcano_data$SYMBOL) &
      volcano_data$SYMBOL != "",
    ,
    drop = FALSE
  ]
  
  label_candidates <- label_candidates[
    order(
      label_candidates$padj,
      -label_candidates$AbsoluteLFC,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  label_candidates <- head(
    label_candidates,
    10
  )
  
  #----------------------------------------------------------
  # Create volcano plot
  #----------------------------------------------------------
  
  volcano_plot <- ggplot(
    volcano_data,
    aes(
      x = log2FoldChange,
      y = NegativeLog10Padj,
      color = Direction
    )
  ) +
    geom_point(
      alpha = 0.65,
      size = 1.5,
      na.rm = TRUE
    ) +
    geom_vline(
      xintercept = c(
        -absolute_lfc_threshold,
        absolute_lfc_threshold
      ),
      linetype = 2,
      linewidth = 0.6,
      color = "grey40"
    ) +
    geom_hline(
      yintercept = -log10(
        adjusted_p_threshold
      ),
      linetype = 2,
      linewidth = 0.6,
      color = "grey40"
    ) +
    geom_text(
      data = label_candidates,
      aes(
        label = SYMBOL
      ),
      size = 3,
      check_overlap = TRUE,
      vjust = -0.7,
      show.legend = FALSE,
      na.rm = TRUE
    ) +
    scale_color_manual(
      values = volcano_colors,
      drop = FALSE
    ) +
    labs(
      title = paste0(
        "Differential Expression: ",
        comparison_title
      ),
      subtitle = paste0(
        nrow(significant_results),
        " significant genes; adjusted P < ",
        adjusted_p_threshold,
        " and |log2FC| \u2265 ",
        absolute_lfc_threshold
      ),
      x = paste0(
        "Shrunken log2 fold change\n",
        "(",
        numerator,
        " relative to ",
        denominator,
        ")"
      ),
      y = expression(
        -log[10](
          adjusted~italic(P)
        )
      ),
      color = "Expression"
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      plot.subtitle = element_text(
        size = 10.5
      ),
      axis.title = element_text(
        face = "bold"
      ),
      legend.title = element_text(
        face = "bold"
      ),
      legend.position = "right",
      plot.margin = margin(
        t = 10,
        r = 15,
        b = 10,
        l = 10
      )
    )
  
  ggsave(
    filename = file.path(
      "figures/differential_expression",
      paste0(
        comparison_label,
        "_volcano.png"
      )
    ),
    plot = volcano_plot,
    width = 7.5,
    height = 6,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(
      "figures/differential_expression",
      paste0(
        comparison_label,
        "_volcano.pdf"
      )
    ),
    plot = volcano_plot,
    width = 7.5,
    height = 6,
    units = "in",
    device = cairo_pdf
  )
  
  #----------------------------------------------------------
  # Store manuscript summary statistics
  #----------------------------------------------------------
  
  summary_list[[comparison_index]] <- data.frame(
    Comparison = comparison_title,
    NumeratorCluster = numerator,
    DenominatorCluster = denominator,
    NumeratorSamples = sum(
      dds$cluster == numerator
    ),
    DenominatorSamples = sum(
      dds$cluster == denominator
    ),
    GenesTested = sum(
      !is.na(results_table$padj)
    ),
    SignificantDEGs = nrow(
      significant_results
    ),
    Upregulated = nrow(
      upregulated_results
    ),
    Downregulated = nrow(
      downregulated_results
    ),
    AdjustedPThreshold =
      adjusted_p_threshold,
    AbsoluteLFCThreshold =
      absolute_lfc_threshold,
    stringsAsFactors = FALSE
  )
  
  cat(
    "Genes with adjusted p-values: ",
    sum(!is.na(results_table$padj)),
    "\n",
    "Significant DEGs: ",
    nrow(significant_results),
    "\n",
    "Upregulated in Cluster ",
    numerator,
    ": ",
    nrow(upregulated_results),
    "\n",
    "Downregulated in Cluster ",
    numerator,
    ": ",
    nrow(downregulated_results),
    "\n",
    sep = ""
  )
}

#------------------------------------------------------------
# Combine and export comparison summary
#------------------------------------------------------------

differential_expression_summary <- do.call(
  rbind,
  summary_list
)

write.csv(
  differential_expression_summary,
  "results/tables/differential_expression_summary.csv",
  row.names = FALSE
)

#------------------------------------------------------------
# Create manuscript-formatted summary table
#------------------------------------------------------------

manuscript_summary <- data.frame(
  Comparison =
    differential_expression_summary$Comparison,
  Samples = paste0(
    differential_expression_summary$NumeratorSamples,
    " vs ",
    differential_expression_summary$DenominatorSamples
  ),
  GenesTested =
    differential_expression_summary$GenesTested,
  SignificantDEGs =
    differential_expression_summary$SignificantDEGs,
  Upregulated =
    differential_expression_summary$Upregulated,
  Downregulated =
    differential_expression_summary$Downregulated,
  stringsAsFactors = FALSE
)

write.csv(
  manuscript_summary,
  "results/tables/differential_expression_manuscript_table.csv",
  row.names = FALSE
)

#------------------------------------------------------------
# Save combined differential-expression object
#------------------------------------------------------------

combined_de_results <- list(
  settings = list(
    adjusted_p_threshold =
      adjusted_p_threshold,
    absolute_lfc_threshold =
      absolute_lfc_threshold,
    shrinkage_method = "apeglm",
    design = "~ cluster"
  ),
  comparisons = comparisons,
  summary =
    differential_expression_summary,
  results = comparison_objects
)

saveRDS(
  combined_de_results,
  "results/r_objects/differential_expression_all_comparisons.rds"
)

#------------------------------------------------------------
# Create DEG-count summary plot
#------------------------------------------------------------

deg_count_plot_data <- rbind(
  data.frame(
    Comparison =
      differential_expression_summary$Comparison,
    Direction = "Upregulated",
    Count =
      differential_expression_summary$Upregulated
  ),
  data.frame(
    Comparison =
      differential_expression_summary$Comparison,
    Direction = "Downregulated",
    Count =
      differential_expression_summary$Downregulated
  )
)

deg_count_plot_data$Comparison <- factor(
  deg_count_plot_data$Comparison,
  levels = rev(
    differential_expression_summary$Comparison
  )
)

deg_count_plot_data$Direction <- factor(
  deg_count_plot_data$Direction,
  levels = c(
    "Downregulated",
    "Upregulated"
  )
)

deg_summary_plot <- ggplot(
  deg_count_plot_data,
  aes(
    x = Comparison,
    y = Count,
    fill = Direction
  )
) +
  geom_col(
    position = "dodge",
    width = 0.7
  ) +
  geom_text(
    aes(
      label = Count
    ),
    position = position_dodge(
      width = 0.7
    ),
    hjust = -0.15,
    size = 3
  ) +
  coord_flip(
    clip = "off"
  ) +
  scale_fill_manual(
    values = c(
      "Downregulated" = "#0072B2",
      "Upregulated" = "#D55E00"
    )
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.12
      )
    )
  ) +
  labs(
    title = "Differentially Expressed Genes Across Cluster Comparisons",
    subtitle = paste0(
      "Adjusted P < ",
      adjusted_p_threshold,
      " and |log2FC| \u2265 ",
      absolute_lfc_threshold
    ),
    x = NULL,
    y = "Number of significant genes",
    fill = "Direction"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 10.5
    ),
    axis.title.x = element_text(
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    ),
    plot.margin = margin(
      t = 10,
      r = 30,
      b = 10,
      l = 10
    )
  )

ggsave(
  filename =
    "figures/differential_expression/DEG_counts_all_comparisons.png",
  plot = deg_summary_plot,
  width = 8,
  height = 5.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename =
    "figures/differential_expression/DEG_counts_all_comparisons.pdf",
  plot = deg_summary_plot,
  width = 8,
  height = 5.5,
  units = "in",
  device = cairo_pdf
)

#------------------------------------------------------------
# Print completion summary
#------------------------------------------------------------

cat(
  "\n",
  "=================================================\n",
  "Differential-expression analysis complete\n",
  "=================================================\n",
  "Model: DESeq2 negative-binomial generalized ",
  "linear model\n",
  "Design: ~ cluster\n",
  "LFC shrinkage: apeglm\n",
  "Adjusted p-value threshold: ",
  adjusted_p_threshold,
  "\n",
  "Absolute log2 fold-change threshold: ",
  absolute_lfc_threshold,
  "\n",
  "Comparisons completed: ",
  nrow(comparisons),
  "\n\n",
  sep = ""
)

print(
  manuscript_summary,
  row.names = FALSE
)

cat(
  "\nMain summary tables:\n",
  "results/tables/differential_expression_summary.csv\n",
  "results/tables/differential_expression_manuscript_table.csv\n",
  "\nComparison tables:\n",
  "results/tables/differential_expression/\n",
  "\nVolcano plots:\n",
  "figures/differential_expression/\n",
  "\nCombined R object:\n",
  "results/r_objects/",
  "differential_expression_all_comparisons.rds\n",
  "\nDESeq2 fitted object:\n",
  "results/r_objects/dds_deseq_all_clusters.rds\n",
  "=================================================\n",
  sep = ""
)