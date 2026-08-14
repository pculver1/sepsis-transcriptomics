##############################################################
# Project:
# Transcriptomic Characterization of Sepsis Using RNA-seq
#
# Script:
# 03_pca_clustering.R
#
# Purpose:
# Reproduce principal component analysis and hierarchical
# clustering from variance-stabilized RNA-seq expression data.
#
# Recovered/validated workflow:
#   PCA:
#     - VST-transformed expression
#     - 500 genes with highest variance
#     - prcomp(center = TRUE, scale. = FALSE)
#
#   Hierarchical clustering:
#     - VST-transformed expression for ALL 19,208 retained genes
#     - Euclidean sample distance
#     - Ward's minimum-variance method (ward.D2)
#     - dendrogram cut into k = 4 clusters
#
# This clustering workflow was recovered by comparison with the
# original validated manuscript assignments and reproduced the
# partition exactly (ARI = 1.0; 100% sample agreement).
#
# Author:
# Patrick Culver
# Last Updated:
# August 2026
##############################################################

rm(list = ls())
graphics.off()
set.seed(12345)

#------------------------------------------------------------
# Packages
#------------------------------------------------------------

required_packages <- c("DESeq2", "ggplot2")

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
library(ggplot2)

#------------------------------------------------------------
# Output directories
#------------------------------------------------------------

for (directory in c(
  "results/tables",
  "results/r_objects",
  "results/logs",
  "figures"
)) {
  dir.create(
    directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

#------------------------------------------------------------
# Input/output files
#------------------------------------------------------------

raw_reconstructed_dds_file <-
  "results/r_objects/dds_from_raw.rds"

legacy_validated_dds_file <-
  "results/r_objects/dds_processed.rds"

vst_file <-
  "results/r_objects/vst_data.rds"

processed_dds_output <-
  "results/r_objects/dds_processed.rds"

processed_vst_output <-
  "results/r_objects/vst_data.rds"

#------------------------------------------------------------
# Load reconstructed DESeq2 object
#------------------------------------------------------------

if (file.exists(raw_reconstructed_dds_file)) {

  dds <- readRDS(raw_reconstructed_dds_file)

  message(
    "Using DESeq2 object reconstructed from raw data: ",
    raw_reconstructed_dds_file
  )

} else {

  stop(
    "Reconstructed DESeq2 object not found:\n",
    raw_reconstructed_dds_file,
    "\nRun scripts/01_preprocessing.R first."
  )
}

if (!is(dds, "DESeqDataSet")) {
  stop("dds_from_raw.rds is not a DESeqDataSet.")
}

#------------------------------------------------------------
# Load VST object
#------------------------------------------------------------

if (!file.exists(vst_file)) {
  stop(
    "VST object not found:\n",
    vst_file,
    "\nRun scripts/02_exploratory_analysis.R first."
  )
}

vst_data <- readRDS(vst_file)

if (!is(vst_data, "DESeqTransform")) {
  stop("vst_data.rds is not a DESeqTransform.")
}

#------------------------------------------------------------
# Validate reconstructed dimensions
#------------------------------------------------------------

expected_genes <- 19208
expected_samples <- 348

if (nrow(dds) != expected_genes) {
  stop(
    "Expected ",
    expected_genes,
    " genes but found ",
    nrow(dds),
    "."
  )
}

if (ncol(dds) != expected_samples) {
  stop(
    "Expected ",
    expected_samples,
    " samples but found ",
    ncol(dds),
    "."
  )
}

if (!setequal(colnames(dds), colnames(vst_data))) {
  stop(
    "Sample identifiers differ between dds_from_raw.rds ",
    "and vst_data.rds."
  )
}

if (!identical(rownames(dds), rownames(vst_data))) {
  stop(
    "Gene identifiers differ between dds_from_raw.rds ",
    "and vst_data.rds."
  )
}

# Match sample order exactly.
vst_data <- vst_data[, colnames(dds)]

stopifnot(
  identical(colnames(dds), colnames(vst_data))
)

#------------------------------------------------------------
# Extract VST matrix
#------------------------------------------------------------

vst_matrix <- assay(vst_data)

if (any(!is.finite(vst_matrix))) {
  stop("VST matrix contains non-finite values.")
}

#============================================================
# PART 1: PCA
#============================================================

#------------------------------------------------------------
# Select 500 genes with highest variance
#------------------------------------------------------------

number_variable_genes <- 500

gene_variances <- apply(
  vst_matrix,
  1,
  var
)

variable_gene_order <- order(
  gene_variances,
  decreasing = TRUE
)

top_variable_genes <- names(
  gene_variances[
    variable_gene_order[
      seq_len(number_variable_genes)
    ]
  ]
)

top_variable_gene_table <- data.frame(
  Rank = seq_len(number_variable_genes),
  Ensembl = top_variable_genes,
  Variance = gene_variances[top_variable_genes],
  stringsAsFactors = FALSE
)

write.csv(
  top_variable_gene_table,
  "results/tables/PCA_top500_variable_genes.csv",
  row.names = FALSE
)

#------------------------------------------------------------
# Perform PCA
#------------------------------------------------------------

cat(
  "\nPerforming PCA using the 500 most variable VST genes...\n"
)

pca_input <- t(
  vst_matrix[
    top_variable_genes,
    ,
    drop = FALSE
  ]
)

pca <- prcomp(
  pca_input,
  center = TRUE,
  scale. = FALSE
)

variance_proportion <- pca$sdev^2 /
  sum(pca$sdev^2)

variance_table <- data.frame(
  PrincipalComponent = paste0(
    "PC",
    seq_along(variance_proportion)
  ),
  VarianceProportion = variance_proportion,
  VariancePercent = variance_proportion * 100,
  CumulativePercent = cumsum(
    variance_proportion * 100
  )
)

write.csv(
  variance_table,
  "results/tables/PCA_variance_explained.csv",
  row.names = FALSE
)

saveRDS(
  pca,
  "results/r_objects/pca_top500_genes.rds"
)

#============================================================
# PART 2: HIERARCHICAL CLUSTERING
#============================================================

#------------------------------------------------------------
# IMPORTANT:
# Clustering uses ALL 19,208 retained VST genes.
# It does NOT use only the 500 genes selected for PCA.
#------------------------------------------------------------

cat(
  "Performing hierarchical clustering using ALL ",
  nrow(vst_matrix),
  " retained VST genes...\n",
  sep = ""
)

clustering_input <- t(vst_matrix)

sample_distances <- dist(
  clustering_input,
  method = "euclidean"
)

hierarchical_clustering <- hclust(
  sample_distances,
  method = "ward.D2"
)

raw_cluster_assignments <- cutree(
  hierarchical_clustering,
  k = 4
)

#------------------------------------------------------------
# Validate cluster-size structure
#------------------------------------------------------------

raw_cluster_sizes <- table(
  raw_cluster_assignments
)

cat(
  "\nRaw cutree cluster sizes:\n"
)

print(raw_cluster_sizes)

# Cluster numbers produced by cutree() are arbitrary.
# The manuscript clusters have unique sizes, allowing the
# recovered raw groups to be mapped deterministically.
expected_cluster_sizes <- c(
  "Cluster 1" = 37,
  "Cluster 2" = 157,
  "Cluster 3" = 134,
  "Cluster 4" = 20
)

if (!setequal(
  as.integer(raw_cluster_sizes),
  as.integer(expected_cluster_sizes)
)) {

  diagnostic_table <- data.frame(
    RawCluster = names(raw_cluster_sizes),
    Samples = as.integer(raw_cluster_sizes),
    stringsAsFactors = FALSE
  )

  write.csv(
    diagnostic_table,
    "results/tables/hierarchical_clustering_size_diagnostic.csv",
    row.names = FALSE
  )

  stop(
    "Recovered clustering does not reproduce the expected ",
    "manuscript cluster-size structure.\nExpected: ",
    paste(expected_cluster_sizes, collapse = ", "),
    "\nObserved: ",
    paste(raw_cluster_sizes, collapse = ", "),
    "\nNo downstream object was overwritten."
  )
}

#------------------------------------------------------------
# Map arbitrary raw labels to manuscript labels by size
#------------------------------------------------------------

size_to_manuscript_label <- setNames(
  names(expected_cluster_sizes),
  as.character(expected_cluster_sizes)
)

raw_to_manuscript_label <- setNames(
  size_to_manuscript_label[
    as.character(
      as.integer(raw_cluster_sizes)
    )
  ],
  names(raw_cluster_sizes)
)

manuscript_cluster_labels <-
  raw_to_manuscript_label[
    as.character(raw_cluster_assignments)
  ]

manuscript_cluster_labels <- factor(
  manuscript_cluster_labels,
  levels = c(
    "Cluster 1",
    "Cluster 2",
    "Cluster 3",
    "Cluster 4"
  )
)

names(manuscript_cluster_labels) <-
  names(raw_cluster_assignments)

if (any(is.na(manuscript_cluster_labels))) {
  stop(
    "Recovered raw clusters could not be mapped to ",
    "manuscript cluster labels."
  )
}

#------------------------------------------------------------
# Optional validation against legacy manuscript assignments
#
# If the old dds_processed.rds exists, compare the newly
# reconstructed partition at the SAMPLE level before replacing
# it. This is a validation step only; the new clusters are
# generated independently from raw/VST data.
#------------------------------------------------------------

legacy_ari <- NA_real_
legacy_accuracy <- NA_real_

adjusted_rand_index <- function(a, b) {

  tab <- table(a, b)

  choose2 <- function(x) {
    x * (x - 1) / 2
  }

  n <- sum(tab)

  cell_sum <- sum(choose2(tab))
  row_sum <- sum(choose2(rowSums(tab)))
  col_sum <- sum(choose2(colSums(tab)))
  total <- choose2(n)

  expected <- (row_sum * col_sum) / total
  maximum <- 0.5 * (row_sum + col_sum)

  denominator <- maximum - expected

  if (denominator == 0) {
    return(NA_real_)
  }

  (cell_sum - expected) / denominator
}

if (file.exists(legacy_validated_dds_file)) {

  legacy_dds <- readRDS(
    legacy_validated_dds_file
  )

  if (
    is(legacy_dds, "DESeqDataSet") &&
    "cluster" %in% colnames(colData(legacy_dds)) &&
    setequal(colnames(legacy_dds), colnames(dds))
  ) {

    legacy_cluster <- as.character(
      colData(legacy_dds)$cluster
    )

    legacy_cluster <- gsub(
      "^Cluster\\s*",
      "",
      legacy_cluster,
      ignore.case = TRUE
    )

    names(legacy_cluster) <-
      colnames(legacy_dds)

    legacy_cluster <- legacy_cluster[
      colnames(dds)
    ]

    new_numeric_cluster <- sub(
      "^Cluster\\s*",
      "",
      as.character(
        manuscript_cluster_labels[
          colnames(dds)
        ]
      )
    )

    legacy_ari <- adjusted_rand_index(
      legacy_cluster,
      new_numeric_cluster
    )

    legacy_accuracy <- mean(
      legacy_cluster == new_numeric_cluster
    )

    cat(
      "\nValidation against existing manuscript assignments:\n",
      "  ARI: ",
      round(legacy_ari, 6),
      "\n",
      "  Exact sample agreement: ",
      round(100 * legacy_accuracy, 2),
      "%\n",
      sep = ""
    )

    if (
      !isTRUE(all.equal(legacy_ari, 1)) ||
      !isTRUE(all.equal(legacy_accuracy, 1))
    ) {

      stop(
        "Newly generated clusters do not exactly reproduce ",
        "the validated manuscript sample assignments.\n",
        "dds_processed.rds was NOT overwritten."
      )
    }
  }
}

#------------------------------------------------------------
# Add recovered clusters to objects
#------------------------------------------------------------

cluster_numeric <- as.integer(
  sub(
    "^Cluster\\s*",
    "",
    as.character(manuscript_cluster_labels)
  )
)

names(cluster_numeric) <-
  names(manuscript_cluster_labels)

dds$cluster <- factor(
  cluster_numeric[colnames(dds)],
  levels = c(1, 2, 3, 4)
)

colData(vst_data)$cluster <- factor(
  cluster_numeric[colnames(vst_data)],
  levels = c(1, 2, 3, 4)
)

#------------------------------------------------------------
# Save cluster tables and objects
#------------------------------------------------------------

cluster_assignment_table <- data.frame(
  Sample = colnames(dds),
  ClusterNumber = as.integer(dds$cluster),
  Cluster = paste0(
    "Cluster ",
    as.integer(dds$cluster)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  cluster_assignment_table,
  "results/tables/cluster_assignments.csv",
  row.names = FALSE
)

cluster_size_table <- data.frame(
  Cluster = names(expected_cluster_sizes),
  Samples = as.integer(expected_cluster_sizes),
  stringsAsFactors = FALSE
)

write.csv(
  cluster_size_table,
  "results/tables/cluster_sizes.csv",
  row.names = FALSE
)

saveRDS(
  sample_distances,
  "results/r_objects/sample_distances_all_vst_genes.rds"
)

saveRDS(
  hierarchical_clustering,
  "results/r_objects/hierarchical_clustering_all_vst_genes.rds"
)

# The newly reproduced object becomes the canonical downstream
# object used by scripts 04 onward.
saveRDS(
  dds,
  processed_dds_output
)

saveRDS(
  vst_data,
  processed_vst_output
)

#============================================================
# PART 3: PCA OUTPUT WITH RECOVERED CLUSTERS
#============================================================

metadata <- as.data.frame(
  colData(dds)
)

metadata <- metadata[
  colnames(dds),
  ,
  drop = FALSE
]

pca_coordinates <- as.data.frame(
  pca$x
)

pca_coordinates$Sample <- rownames(
  pca_coordinates
)

pca_df <- data.frame(
  Sample = pca_coordinates$Sample,
  PC1 = pca_coordinates$PC1,
  PC2 = pca_coordinates$PC2,
  PC3 = pca_coordinates$PC3,
  Cluster = factor(
    paste0(
      "Cluster ",
      metadata[
        pca_coordinates$Sample,
        "cluster"
      ]
    ),
    levels = c(
      "Cluster 1",
      "Cluster 2",
      "Cluster 3",
      "Cluster 4"
    )
  ),
  Cohort = if (
    "cohort" %in% colnames(metadata)
  ) {
    metadata[
      pca_coordinates$Sample,
      "cohort"
    ]
  } else {
    NA
  },
  Age = if (
    "age" %in% colnames(metadata)
  ) {
    metadata[
      pca_coordinates$Sample,
      "age"
    ]
  } else {
    NA
  },
  SOFA = if (
    "sofa" %in% colnames(metadata)
  ) {
    metadata[
      pca_coordinates$Sample,
      "sofa"
    ]
  } else {
    NA
  },
  stringsAsFactors = FALSE
)

write.csv(
  pca_df,
  "results/tables/PCA_coordinates.csv",
  row.names = FALSE
)

#------------------------------------------------------------
# Figure appearance
#------------------------------------------------------------

cluster_colors <- c(
  "Cluster 1" = "#0072B2",
  "Cluster 2" = "#E69F00",
  "Cluster 3" = "#009E73",
  "Cluster 4" = "#CC79A7"
)

cluster_shapes <- c(
  "Cluster 1" = 16,
  "Cluster 2" = 17,
  "Cluster 3" = 15,
  "Cluster 4" = 18
)

pc1_percent <- round(
  variance_table$VariancePercent[1],
  1
)

pc2_percent <- round(
  variance_table$VariancePercent[2],
  1
)

#------------------------------------------------------------
# PCA plot
#------------------------------------------------------------

pca_plot <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = Cluster,
    shape = Cluster
  )
) +
  stat_ellipse(
    aes(
      group = Cluster,
      color = Cluster
    ),
    type = "norm",
    level = 0.95,
    linewidth = 0.8,
    linetype = 2,
    show.legend = FALSE
  ) +
  geom_point(
    size = 3,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = cluster_colors,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = cluster_shapes,
    drop = FALSE
  ) +
  labs(
    title = "Principal Component Analysis of Sepsis Samples",
    subtitle = paste0(
      "PCA: top 500 variable genes; ",
      "clusters: all 19,208 VST genes"
    ),
    x = paste0(
      "PC1 (",
      pc1_percent,
      "% variance)"
    ),
    y = paste0(
      "PC2 (",
      pc2_percent,
      "% variance)"
    ),
    color = "Transcriptomic cluster",
    shape = "Transcriptomic cluster"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 11
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    ),
    legend.position = "right",
    plot.margin = margin(
      10, 10, 10, 10
    )
  )

ggsave(
  "figures/PCA_clusters.png",
  pca_plot,
  width = 7.5,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  "figures/PCA_clusters.pdf",
  pca_plot,
  width = 7.5,
  height = 6,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

#------------------------------------------------------------
# Scree plot
#------------------------------------------------------------

number_of_pcs <- min(
  10,
  nrow(variance_table)
)

scree_df <- variance_table[
  seq_len(number_of_pcs),
  ,
  drop = FALSE
]

scree_df$PrincipalComponent <- factor(
  scree_df$PrincipalComponent,
  levels = scree_df$PrincipalComponent
)

scree_plot <- ggplot(
  scree_df,
  aes(
    x = PrincipalComponent,
    y = VariancePercent
  )
) +
  geom_col(width = 0.75) +
  geom_text(
    aes(
      label = paste0(
        round(VariancePercent, 1),
        "%"
      )
    ),
    vjust = -0.4,
    size = 3.2
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.12)
    )
  ) +
  labs(
    title = "Variance Explained by Principal Components",
    subtitle = "PCA based on the 500 most variable genes",
    x = "Principal component",
    y = "Variance explained (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 11
    ),
    axis.title = element_text(
      face = "bold"
    )
  )

ggsave(
  "figures/PCA_variance.png",
  scree_plot,
  width = 7,
  height = 4.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  "figures/PCA_variance.pdf",
  scree_plot,
  width = 7,
  height = 4.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

#============================================================
# PART 4: REPRODUCIBILITY LOG
#============================================================

final_sizes <- table(dds$cluster)

log_lines <- c(
  "PCA and Hierarchical Clustering Summary",
  "==================================================",
  paste("Date:", Sys.Date()),
  paste("R version:", R.version.string),
  paste(
    "DESeq2 version:",
    as.character(packageVersion("DESeq2"))
  ),
  "",
  "INPUT",
  paste(
    "DESeq2 object:",
    raw_reconstructed_dds_file
  ),
  paste("VST object:", vst_file),
  paste("Samples:", ncol(dds)),
  paste("Retained genes:", nrow(vst_matrix)),
  "",
  "PCA",
  "Expression: VST",
  "Genes: 500 genes with highest variance",
  "Function: prcomp",
  "center = TRUE",
  "scale = FALSE",
  paste("PC1 variance (%):", pc1_percent),
  paste("PC2 variance (%):", pc2_percent),
  "",
  "HIERARCHICAL CLUSTERING",
  "Expression: VST",
  paste(
    "Genes:",
    nrow(vst_matrix),
    "(all retained genes)"
  ),
  "Distance: Euclidean",
  "Linkage: ward.D2",
  "k = 4",
  "",
  paste(
    "Final cluster sizes:",
    paste(
      paste0(
        "Cluster ",
        names(final_sizes),
        "=",
        as.integer(final_sizes)
      ),
      collapse = ", "
    )
  ),
  "",
  paste(
    "Legacy validation ARI:",
    ifelse(
      is.na(legacy_ari),
      "not performed",
      round(legacy_ari, 6)
    )
  ),
  paste(
    "Legacy exact sample agreement:",
    ifelse(
      is.na(legacy_accuracy),
      "not performed",
      paste0(
        round(100 * legacy_accuracy, 2),
        "%"
      )
    )
  ),
  "",
  paste(
    "Canonical downstream object:",
    processed_dds_output
  )
)

writeLines(
  log_lines,
  "results/logs/pca_clustering_summary.txt"
)

#------------------------------------------------------------
# Completion message
#------------------------------------------------------------

cat(
  "\n",
  "============================================================\n",
  "03_pca_clustering.R COMPLETE\n",
  "============================================================\n",
  "PCA:\n",
  "  Top 500 variable VST genes\n",
  "  PC1 variance: ",
  pc1_percent,
  "%\n",
  "  PC2 variance: ",
  pc2_percent,
  "%\n\n",
  "Hierarchical clustering:\n",
  "  All ",
  nrow(vst_matrix),
  " retained VST genes\n",
  "  Euclidean distance\n",
  "  ward.D2 linkage\n",
  "  k = 4\n\n",
  "Final cluster sizes:\n",
  paste(
    paste0(
      "  Cluster ",
      names(final_sizes),
      ": ",
      as.integer(final_sizes)
    ),
    collapse = "\n"
  ),
  "\n\n",
  "Canonical downstream object saved:\n",
  "  ",
  processed_dds_output,
  "\n",
  "============================================================\n",
  sep = ""
)
