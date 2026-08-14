##############################################################
# Project:
# Transcriptomic Characterization of Sepsis Using RNA-seq
#
# Script:
# 02_exploratory_analysis.R
#
# Purpose:
# Perform exploratory data analysis and quality-control
# summaries on the reconstructed sepsis RNA-seq dataset and
# create the variance-stabilized expression object used by
# downstream PCA and hierarchical clustering.
#
# Input:
#   results/r_objects/dds_from_raw.rds
#
# Output:
#   results/r_objects/vst_data.rds
#   exploratory summary tables and QC figures
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
# Load required packages
#------------------------------------------------------------

required_packages <- c(
  "DESeq2",
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
library(ggplot2)

#------------------------------------------------------------
# Create output directories
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
# Input file
#------------------------------------------------------------

dds_file <- "results/r_objects/dds_from_raw.rds"

if (!file.exists(dds_file)) {
  stop(
    "Reconstructed DESeq2 object not found:\n",
    dds_file,
    "\nRun scripts/01_preprocessing.R first."
  )
}

#------------------------------------------------------------
# Load reconstructed DESeq2 object
#------------------------------------------------------------

dds <- readRDS(dds_file)

if (!is(dds, "DESeqDataSet")) {
  stop(
    "dds_from_raw.rds is not a DESeqDataSet."
  )
}

#------------------------------------------------------------
# Validate expected dimensions
#------------------------------------------------------------

expected_samples <- 348
expected_genes <- 19208

if (ncol(dds) != expected_samples) {
  stop(
    "Unexpected number of samples. Expected ",
    expected_samples,
    " but found ",
    ncol(dds),
    "."
  )
}

if (nrow(dds) != expected_genes) {
  stop(
    "Unexpected number of genes. Expected ",
    expected_genes,
    " but found ",
    nrow(dds),
    "."
  )
}

#------------------------------------------------------------
# Prepare metadata
#------------------------------------------------------------

metadata <- as.data.frame(
  colData(dds)
)

metadata <- metadata[
  colnames(dds),
  ,
  drop = FALSE
]

# Add a stable sample identifier column if it is absent.
if (!"sample_id" %in% colnames(metadata)) {
  metadata$sample_id <- rownames(metadata)
}

# Convert clinical variables to numeric when present.
if ("age" %in% colnames(metadata)) {
  metadata$age <- suppressWarnings(
    as.numeric(
      as.character(metadata$age)
    )
  )
}

if ("sofa" %in% colnames(metadata)) {
  metadata$sofa <- suppressWarnings(
    as.numeric(
      as.character(metadata$sofa)
    )
  )
}

#------------------------------------------------------------
# Calculate library sizes
#------------------------------------------------------------

library_size <- colSums(
  counts(dds)
)

metadata$LibrarySize <- library_size

#------------------------------------------------------------
# Variance stabilizing transformation
#
# This object is used by 03_pca_clustering.R.
# VST is applied before PCA and clustering and therefore does
# not require cluster assignments.
#------------------------------------------------------------

cat(
  "Running variance stabilizing transformation...\n"
)

vst_data <- vst(
  dds,
  blind = TRUE
)

saveRDS(
  vst_data,
  "results/r_objects/vst_data.rds"
)

cat(
  "VST object saved to results/r_objects/vst_data.rds\n"
)

#------------------------------------------------------------
# Exploratory summary statistics
#------------------------------------------------------------

metric_names <- c(
  "Samples",
  "Genes",
  "Median library size",
  "Minimum library size",
  "Maximum library size"
)

metric_values <- c(
  ncol(dds),
  nrow(dds),
  median(library_size, na.rm = TRUE),
  min(library_size, na.rm = TRUE),
  max(library_size, na.rm = TRUE)
)

if ("age" %in% colnames(metadata)) {
  metric_names <- c(
    metric_names,
    "Mean age",
    "Median age",
    "Age missing"
  )

  metric_values <- c(
    metric_values,
    round(
      mean(metadata$age, na.rm = TRUE),
      1
    ),
    median(
      metadata$age,
      na.rm = TRUE
    ),
    sum(
      is.na(metadata$age)
    )
  )
}

if ("sofa" %in% colnames(metadata)) {
  metric_names <- c(
    metric_names,
    "Mean SOFA",
    "Median SOFA",
    "SOFA missing"
  )

  metric_values <- c(
    metric_values,
    round(
      mean(metadata$sofa, na.rm = TRUE),
      1
    ),
    median(
      metadata$sofa,
      na.rm = TRUE
    ),
    sum(
      is.na(metadata$sofa)
    )
  )
}

summary_table <- data.frame(
  Metric = metric_names,
  Value = metric_values,
  stringsAsFactors = FALSE
)

write.csv(
  summary_table,
  "results/tables/exploratory_summary.csv",
  row.names = FALSE
)

#------------------------------------------------------------
# Missing-value summary
#------------------------------------------------------------

missing_table <- data.frame(
  Variable = names(metadata),
  Missing = vapply(
    metadata,
    function(x) sum(is.na(x)),
    FUN.VALUE = integer(1)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  missing_table,
  "results/tables/missing_values.csv",
  row.names = FALSE
)

#------------------------------------------------------------
# Save metadata snapshot used for exploratory analysis
#------------------------------------------------------------

write.csv(
  metadata,
  "results/tables/exploratory_metadata.csv",
  row.names = FALSE
)

#============================================================
# FIGURES
#============================================================

#------------------------------------------------------------
# Library-size distribution
#------------------------------------------------------------

library_df <- data.frame(
  Sample = metadata$sample_id,
  LibrarySize = library_size,
  stringsAsFactors = FALSE
)

library_size_distribution <- ggplot(
  library_df,
  aes(
    x = log10(LibrarySize)
  )
) +
  geom_boxplot(
    width = 0.45,
    outlier.alpha = 0.65
  ) +
  labs(
    title = "Distribution of Sequencing Library Sizes",
    x = expression(log[10]~"Library Size"),
    y = NULL
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.title.x = element_text(
      face = "bold"
    )
  )

ggsave(
  "figures/library_size_distribution.png",
  library_size_distribution,
  width = 6,
  height = 4.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  "figures/library_size_distribution.pdf",
  library_size_distribution,
  width = 6,
  height = 4.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

#------------------------------------------------------------
# Library sizes by sample
#------------------------------------------------------------

library_plot <- ggplot(
  library_df,
  aes(
    x = reorder(
      Sample,
      LibrarySize
    ),
    y = LibrarySize
  )
) +
  geom_col() +
  labs(
    title = "Library Sizes Across Sepsis Samples",
    x = "Samples",
    y = "Total counts"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(
      face = "bold"
    )
  )

ggsave(
  "figures/library_sizes.png",
  library_plot,
  width = 9,
  height = 4,
  units = "in",
  dpi = 300,
  bg = "white"
)

#------------------------------------------------------------
# Age distribution, when available
#------------------------------------------------------------

if ("age" %in% colnames(metadata)) {

  age_plot <- ggplot(
    metadata,
    aes(x = age)
  ) +
    geom_histogram(
      bins = 20,
      color = "black",
      fill = "steelblue"
    ) +
    labs(
      title = "Age Distribution",
      x = "Age (years)",
      y = "Number of samples"
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      axis.title = element_text(
        face = "bold"
      )
    )

  ggsave(
    "figures/age_distribution.png",
    age_plot,
    width = 6,
    height = 4,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

#------------------------------------------------------------
# SOFA distribution, when available
#------------------------------------------------------------

if ("sofa" %in% colnames(metadata)) {

  sofa_plot <- ggplot(
    metadata,
    aes(x = sofa)
  ) +
    geom_histogram(
      bins = 15,
      color = "black",
      fill = "darkred"
    ) +
    labs(
      title = "SOFA Score Distribution",
      x = "SOFA score",
      y = "Number of samples"
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      axis.title = element_text(
        face = "bold"
      )
    )

  ggsave(
    "figures/sofa_distribution.png",
    sofa_plot,
    width = 6,
    height = 4,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

#------------------------------------------------------------
# Cohort distribution, when available
#------------------------------------------------------------

if ("cohort" %in% colnames(metadata)) {

  cohort_counts <- as.data.frame(
    table(metadata$cohort),
    stringsAsFactors = FALSE
  )

  colnames(cohort_counts) <- c(
    "Cohort",
    "Samples"
  )

  cohort_counts <- cohort_counts[
    !is.na(cohort_counts$Cohort) &
      cohort_counts$Cohort != "",
    ,
    drop = FALSE
  ]

  write.csv(
    cohort_counts,
    "results/tables/cohort_sample_counts.csv",
    row.names = FALSE
  )

  cohort_plot <- ggplot(
    cohort_counts,
    aes(
      x = Cohort,
      y = Samples
    )
  ) +
    geom_col() +
    geom_text(
      aes(label = Samples),
      vjust = -0.4,
      size = 3.5
    ) +
    scale_y_continuous(
      expand = expansion(
        mult = c(0, 0.12)
      )
    ) +
    labs(
      title = "Distribution of Samples by Cohort",
      x = "Cohort",
      y = "Number of samples"
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      axis.title = element_text(
        face = "bold"
      )
    )

  ggsave(
    "figures/cohort_distribution.png",
    cohort_plot,
    width = 7,
    height = 5,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  ggsave(
    "figures/cohort_distribution.pdf",
    cohort_plot,
    width = 7,
    height = 5,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
}

#------------------------------------------------------------
# Reproducibility log
#------------------------------------------------------------

log_lines <- c(
  "Exploratory Analysis Summary",
  "==================================================",
  paste("Date:", Sys.Date()),
  paste("R version:", R.version.string),
  paste(
    "DESeq2 version:",
    as.character(
      packageVersion("DESeq2")
    )
  ),
  "",
  paste("Input object:", dds_file),
  paste("Samples:", ncol(dds)),
  paste("Genes:", nrow(dds)),
  "",
  "Variance-stabilizing transformation:",
  "  DESeq2::vst",
  "  blind = TRUE",
  "",
  paste(
    "VST output:",
    "results/r_objects/vst_data.rds"
  )
)

writeLines(
  log_lines,
  "results/logs/exploratory_analysis_summary.txt"
)

#------------------------------------------------------------
# Completion summary
#------------------------------------------------------------

cat(
  "\n",
  "============================================================\n",
  "02_exploratory_analysis.R COMPLETE\n",
  "============================================================\n",
  "Input object:\n",
  "  ",
  dds_file,
  "\n\n",
  "Samples: ",
  ncol(dds),
  "\n",
  "Genes:   ",
  nrow(dds),
  "\n\n",
  "VST object saved:\n",
  "  results/r_objects/vst_data.rds\n\n",
  "Tables saved to:\n",
  "  results/tables/\n\n",
  "QC figures saved to:\n",
  "  figures/\n",
  "============================================================\n",
  sep = ""
)
