##############################################################
# Project:
# Transcriptomic Characterization of Sepsis Using RNA-seq
#
# Script:
# 04_cluster_statistics.R
#
# Purpose:
# Compare age and SOFA score distributions among the four
# reproducibly generated transcriptomic clusters.
#
# Statistical workflow:
#   - Descriptive statistics by cluster
#   - Kruskal-Wallis tests for age and SOFA
#   - Epsilon-squared effect sizes
#   - Dunn post-hoc tests with Benjamini-Hochberg correction
#     for any variable with a significant overall
#     Kruskal-Wallis test
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
  "ggplot2",
  "FSA"
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
library(FSA)

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
# Load processed DESeq2 object
#------------------------------------------------------------

dds_file <- "results/r_objects/dds_processed.rds"

if (!file.exists(dds_file)) {
  stop(
    "Processed DESeq2 object not found:\n",
    dds_file,
    "\nRun scripts/03_pca_clustering.R first."
  )
}

dds <- readRDS(dds_file)

if (!is(dds, "DESeqDataSet")) {
  stop(
    "dds_processed.rds is not a DESeqDataSet object."
  )
}

#------------------------------------------------------------
# Validate required metadata
#------------------------------------------------------------

required_metadata <- c(
  "age",
  "sofa",
  "cluster"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(colData(dds))
)

if (length(missing_metadata) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_metadata, collapse = ", ")
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

if (!"sample_id" %in% colnames(metadata)) {
  metadata$sample_id <- rownames(metadata)
}

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
  ),
  labels = c(
    "Cluster 1",
    "Cluster 2",
    "Cluster 3",
    "Cluster 4"
  )
)

metadata$age <- suppressWarnings(
  as.numeric(
    as.character(metadata$age)
  )
)

metadata$sofa <- suppressWarnings(
  as.numeric(
    as.character(metadata$sofa)
  )
)

if (any(is.na(metadata$cluster))) {
  stop(
    "One or more samples have missing or invalid cluster assignments."
  )
}

#------------------------------------------------------------
# Validate expected cluster sizes
#------------------------------------------------------------

expected_cluster_sizes <- c(
  "Cluster 1" = 37,
  "Cluster 2" = 157,
  "Cluster 3" = 134,
  "Cluster 4" = 20
)

observed_cluster_sizes <- table(
  metadata$cluster
)

if (!identical(
  as.integer(observed_cluster_sizes),
  as.integer(expected_cluster_sizes)
)) {
  stop(
    "Observed cluster sizes do not match the validated values.\n",
    "Expected: ",
    paste(expected_cluster_sizes, collapse = ", "),
    "\nObserved: ",
    paste(observed_cluster_sizes, collapse = ", ")
  )
}

#------------------------------------------------------------
# Helper for descriptive statistics
#------------------------------------------------------------

summarize_variable <- function(
    data,
    variable_name
) {

  split_values <- split(
    data[[variable_name]],
    data$cluster
  )

  summary_rows <- lapply(
    names(split_values),
    function(cluster_name) {

      values <- split_values[[cluster_name]]

      nonmissing_values <- values[
        !is.na(values)
      ]

      data.frame(
        Variable = variable_name,
        Cluster = cluster_name,
        TotalSamples = length(values),
        Nonmissing = length(nonmissing_values),
        Missing = sum(is.na(values)),
        Mean = mean(
          nonmissing_values,
          na.rm = TRUE
        ),
        SD = sd(
          nonmissing_values,
          na.rm = TRUE
        ),
        Median = median(
          nonmissing_values,
          na.rm = TRUE
        ),
        Q1 = unname(
          quantile(
            nonmissing_values,
            probs = 0.25,
            na.rm = TRUE
          )
        ),
        Q3 = unname(
          quantile(
            nonmissing_values,
            probs = 0.75,
            na.rm = TRUE
          )
        ),
        Minimum = min(
          nonmissing_values,
          na.rm = TRUE
        ),
        Maximum = max(
          nonmissing_values,
          na.rm = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  )

  do.call(
    rbind,
    summary_rows
  )
}

#------------------------------------------------------------
# Descriptive-statistics tables
#------------------------------------------------------------

age_summary <- summarize_variable(
  metadata,
  "age"
)

sofa_summary <- summarize_variable(
  metadata,
  "sofa"
)

clinical_summary <- rbind(
  age_summary,
  sofa_summary
)

numeric_columns <- c(
  "Mean",
  "SD",
  "Median",
  "Q1",
  "Q3",
  "Minimum",
  "Maximum"
)

clinical_summary[numeric_columns] <- lapply(
  clinical_summary[numeric_columns],
  round,
  digits = 2
)

write.csv(
  clinical_summary,
  "results/tables/cluster_clinical_summary.csv",
  row.names = FALSE
)

#------------------------------------------------------------
# Manuscript-style summary table
#------------------------------------------------------------

format_median_iqr <- function(
    median_value,
    q1_value,
    q3_value
) {

  paste0(
    round(median_value, 1),
    " (",
    round(q1_value, 1),
    "-",
    round(q3_value, 1),
    ")"
  )
}

age_manuscript <- age_summary
age_manuscript$Summary <- mapply(
  format_median_iqr,
  age_manuscript$Median,
  age_manuscript$Q1,
  age_manuscript$Q3
)

sofa_manuscript <- sofa_summary
sofa_manuscript$Summary <- mapply(
  format_median_iqr,
  sofa_manuscript$Median,
  sofa_manuscript$Q1,
  sofa_manuscript$Q3
)

manuscript_table <- merge(
  age_manuscript[
    ,
    c(
      "Cluster",
      "TotalSamples",
      "Summary"
    )
  ],
  sofa_manuscript[
    ,
    c(
      "Cluster",
      "Summary"
    )
  ],
  by = "Cluster",
  suffixes = c(
    "_Age",
    "_SOFA"
  ),
  sort = FALSE
)

colnames(manuscript_table) <- c(
  "Cluster",
  "SampleSize",
  "AgeMedianIQR",
  "SOFAMedianIQR"
)

write.csv(
  manuscript_table,
  "results/tables/cluster_clinical_manuscript_table.csv",
  row.names = FALSE
)

#------------------------------------------------------------
# Complete-case data
#------------------------------------------------------------

age_complete <- metadata[
  !is.na(metadata$age),
  ,
  drop = FALSE
]

sofa_complete <- metadata[
  !is.na(metadata$sofa),
  ,
  drop = FALSE
]

#------------------------------------------------------------
# Kruskal-Wallis tests
#------------------------------------------------------------

age_kruskal <- kruskal.test(
  age ~ cluster,
  data = age_complete
)

sofa_kruskal <- kruskal.test(
  sofa ~ cluster,
  data = sofa_complete
)

#------------------------------------------------------------
# Effect sizes
#
# Epsilon-squared:
# epsilon^2 = (H - k + 1) / (n - k)
#------------------------------------------------------------

kruskal_epsilon_squared <- function(
    test_result,
    sample_size,
    number_groups
) {

  effect_size <- (
    as.numeric(test_result$statistic) -
      number_groups +
      1
  ) / (
    sample_size -
      number_groups
  )

  max(
    0,
    effect_size
  )
}

age_effect_size <- kruskal_epsilon_squared(
  age_kruskal,
  nrow(age_complete),
  nlevels(age_complete$cluster)
)

sofa_effect_size <- kruskal_epsilon_squared(
  sofa_kruskal,
  nrow(sofa_complete),
  nlevels(sofa_complete$cluster)
)

overall_tests <- data.frame(
  Variable = c(
    "Age",
    "SOFA"
  ),
  Test = c(
    "Kruskal-Wallis",
    "Kruskal-Wallis"
  ),
  CompleteSamples = c(
    nrow(age_complete),
    nrow(sofa_complete)
  ),
  Statistic = c(
    as.numeric(age_kruskal$statistic),
    as.numeric(sofa_kruskal$statistic)
  ),
  DegreesOfFreedom = c(
    as.numeric(age_kruskal$parameter),
    as.numeric(sofa_kruskal$parameter)
  ),
  PValue = c(
    age_kruskal$p.value,
    sofa_kruskal$p.value
  ),
  EpsilonSquared = c(
    age_effect_size,
    sofa_effect_size
  )
)

write.csv(
  overall_tests,
  "results/tables/cluster_kruskal_wallis_tests.csv",
  row.names = FALSE
)

#------------------------------------------------------------
# Dunn-test helper
#------------------------------------------------------------

run_dunn_bh <- function(
    formula_object,
    data_object,
    variable_label
) {

  dunn_result <- FSA::dunnTest(
    formula_object,
    data = data_object,
    method = "bh"
  )

  dunn_table <- dunn_result$res

  data.frame(
    Variable = variable_label,
    Comparison = dunn_table$Comparison,
    Z = dunn_table$Z,
    UnadjustedPValue = dunn_table$P.unadj,
    AdjustedPValue = dunn_table$P.adj,
    stringsAsFactors = FALSE
  )
}

significance_label <- function(p_value) {

  if (is.na(p_value)) {
    return(NA_character_)
  }

  if (p_value < 0.0001) {
    return("****")
  }

  if (p_value < 0.001) {
    return("***")
  }

  if (p_value < 0.01) {
    return("**")
  }

  if (p_value < 0.05) {
    return("*")
  }

  "ns"
}

#------------------------------------------------------------
# Dunn post-hoc test for age
#------------------------------------------------------------

if (age_kruskal$p.value < 0.05) {

  age_dunn_output <- run_dunn_bh(
    age ~ cluster,
    age_complete,
    "Age"
  )

  age_dunn_output$Significance <- vapply(
    age_dunn_output$AdjustedPValue,
    significance_label,
    FUN.VALUE = character(1)
  )

} else {

  age_dunn_output <- data.frame(
    Variable = character(0),
    Comparison = character(0),
    Z = numeric(0),
    UnadjustedPValue = numeric(0),
    AdjustedPValue = numeric(0),
    Significance = character(0),
    stringsAsFactors = FALSE
  )
}

write.csv(
  age_dunn_output,
  "results/tables/cluster_age_dunn_tests.csv",
  row.names = FALSE
)

#------------------------------------------------------------
# Dunn post-hoc test for SOFA
#------------------------------------------------------------

if (sofa_kruskal$p.value < 0.05) {

  sofa_dunn_output <- run_dunn_bh(
    sofa ~ cluster,
    sofa_complete,
    "SOFA"
  )

  sofa_dunn_output$Significance <- vapply(
    sofa_dunn_output$AdjustedPValue,
    significance_label,
    FUN.VALUE = character(1)
  )

} else {

  sofa_dunn_output <- data.frame(
    Variable = character(0),
    Comparison = character(0),
    Z = numeric(0),
    UnadjustedPValue = numeric(0),
    AdjustedPValue = numeric(0),
    Significance = character(0),
    stringsAsFactors = FALSE
  )
}

write.csv(
  sofa_dunn_output,
  "results/tables/cluster_sofa_dunn_tests.csv",
  row.names = FALSE
)

#------------------------------------------------------------
# Combined post-hoc table
#------------------------------------------------------------

combined_posthoc <- rbind(
  age_dunn_output,
  sofa_dunn_output
)

write.csv(
  combined_posthoc,
  "results/tables/cluster_dunn_tests_all.csv",
  row.names = FALSE
)

#------------------------------------------------------------
# Save statistical objects
#------------------------------------------------------------

statistical_results <- list(
  age_kruskal = age_kruskal,
  sofa_kruskal = sofa_kruskal,
  age_effect_size = age_effect_size,
  sofa_effect_size = sofa_effect_size,
  age_dunn = age_dunn_output,
  sofa_dunn = sofa_dunn_output
)

saveRDS(
  statistical_results,
  "results/r_objects/cluster_statistics.rds"
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

format_p_value <- function(p_value) {

  if (p_value < 0.001) {
    return("P < 0.001")
  }

  paste0(
    "P = ",
    formatC(
      p_value,
      format = "f",
      digits = 3
    )
  )
}

age_p_label <- paste0(
  "Kruskal-Wallis ",
  format_p_value(
    age_kruskal$p.value
  )
)

sofa_p_label <- paste0(
  "Kruskal-Wallis ",
  format_p_value(
    sofa_kruskal$p.value
  )
)

#------------------------------------------------------------
# Age boxplot
#------------------------------------------------------------

age_plot <- ggplot(
  age_complete,
  aes(
    x = cluster,
    y = age,
    fill = cluster
  )
) +
  geom_boxplot(
    width = 0.65,
    alpha = 0.75,
    outlier.shape = NA,
    linewidth = 0.7
  ) +
  geom_jitter(
    aes(
      color = cluster
    ),
    width = 0.16,
    height = 0,
    size = 1.7,
    alpha = 0.55,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = cluster_colors,
    drop = FALSE
  ) +
  scale_color_manual(
    values = cluster_colors,
    drop = FALSE
  ) +
  labs(
    title = "Age Distribution by Transcriptomic Cluster",
    subtitle = age_p_label,
    x = "Transcriptomic cluster",
    y = "Age (years)"
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
      size = 11
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.position = "none",
    plot.margin = margin(
      10, 10, 10, 10
    )
  )

ggsave(
  "figures/age_by_cluster.png",
  age_plot,
  width = 7,
  height = 5.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  "figures/age_by_cluster.pdf",
  age_plot,
  width = 7,
  height = 5.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

#------------------------------------------------------------
# SOFA boxplot
#------------------------------------------------------------

sofa_plot <- ggplot(
  sofa_complete,
  aes(
    x = cluster,
    y = sofa,
    fill = cluster
  )
) +
  geom_boxplot(
    width = 0.65,
    alpha = 0.75,
    outlier.shape = NA,
    linewidth = 0.7
  ) +
  geom_jitter(
    aes(
      color = cluster
    ),
    width = 0.16,
    height = 0,
    size = 1.7,
    alpha = 0.55,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = cluster_colors,
    drop = FALSE
  ) +
  scale_color_manual(
    values = cluster_colors,
    drop = FALSE
  ) +
  labs(
    title = "SOFA Score Distribution by Transcriptomic Cluster",
    subtitle = sofa_p_label,
    x = "Transcriptomic cluster",
    y = "SOFA score"
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
      size = 11
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.position = "none",
    plot.margin = margin(
      10, 10, 10, 10
    )
  )

ggsave(
  "figures/sofa_by_cluster.png",
  sofa_plot,
  width = 7,
  height = 5.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  "figures/sofa_by_cluster.pdf",
  sofa_plot,
  width = 7,
  height = 5.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

#------------------------------------------------------------
# Reproducibility log
#------------------------------------------------------------

log_lines <- c(
  "Cluster Clinical Statistics Summary",
  "==================================================",
  paste("Date:", Sys.Date()),
  paste("R version:", R.version.string),
  paste(
    "FSA version:",
    as.character(
      packageVersion("FSA")
    )
  ),
  "",
  paste("Input object:", dds_file),
  paste("Samples:", nrow(metadata)),
  "",
  "AGE",
  paste(
    "Complete observations:",
    nrow(age_complete)
  ),
  paste(
    "Kruskal-Wallis statistic:",
    round(
      as.numeric(age_kruskal$statistic),
      6
    )
  ),
  paste(
    "Kruskal-Wallis df:",
    as.numeric(age_kruskal$parameter)
  ),
  paste(
    "Kruskal-Wallis p-value:",
    age_kruskal$p.value
  ),
  paste(
    "Epsilon-squared:",
    age_effect_size
  ),
  paste(
    "Post-hoc performed:",
    age_kruskal$p.value < 0.05
  ),
  "Post-hoc method: Dunn test",
  "P-value adjustment: Benjamini-Hochberg",
  "",
  "SOFA",
  paste(
    "Complete observations:",
    nrow(sofa_complete)
  ),
  paste(
    "Kruskal-Wallis statistic:",
    round(
      as.numeric(sofa_kruskal$statistic),
      6
    )
  ),
  paste(
    "Kruskal-Wallis df:",
    as.numeric(sofa_kruskal$parameter)
  ),
  paste(
    "Kruskal-Wallis p-value:",
    sofa_kruskal$p.value
  ),
  paste(
    "Epsilon-squared:",
    sofa_effect_size
  ),
  paste(
    "Post-hoc performed:",
    sofa_kruskal$p.value < 0.05
  ),
  "Post-hoc method: Dunn test",
  "P-value adjustment: Benjamini-Hochberg"
)

writeLines(
  log_lines,
  "results/logs/cluster_statistics_summary.txt"
)

#------------------------------------------------------------
# Completion summary
#------------------------------------------------------------

cat(
  "\n",
  "============================================================\n",
  "04_cluster_statistics.R COMPLETE\n",
  "============================================================\n",
  "Samples: ",
  nrow(metadata),
  "\n\n",
  "Age:\n",
  "  Complete observations: ",
  nrow(age_complete),
  "\n",
  "  Kruskal-Wallis statistic: ",
  round(
    as.numeric(age_kruskal$statistic),
    3
  ),
  "\n",
  "  p-value: ",
  format.pval(
    age_kruskal$p.value,
    digits = 5,
    eps = 0.00001
  ),
  "\n",
  "  Epsilon-squared: ",
  round(
    age_effect_size,
    4
  ),
  "\n",
  "  Post-hoc: ",
  ifelse(
    age_kruskal$p.value < 0.05,
    "Dunn + BH",
    "not performed"
  ),
  "\n\n",
  "Age Dunn tests:\n",
  sep = ""
)

if (nrow(age_dunn_output) > 0) {
  print(
    age_dunn_output[
      ,
      c(
        "Comparison",
        "AdjustedPValue",
        "Significance"
      )
    ],
    row.names = FALSE
  )
} else {
  cat("  None\n")
}

cat(
  "\nSOFA:\n",
  "  Complete observations: ",
  nrow(sofa_complete),
  "\n",
  "  Kruskal-Wallis statistic: ",
  round(
    as.numeric(sofa_kruskal$statistic),
    3
  ),
  "\n",
  "  p-value: ",
  format.pval(
    sofa_kruskal$p.value,
    digits = 5,
    eps = 0.00001
  ),
  "\n",
  "  Epsilon-squared: ",
  round(
    sofa_effect_size,
    4
  ),
  "\n",
  "  Post-hoc: ",
  ifelse(
    sofa_kruskal$p.value < 0.05,
    "Dunn + BH",
    "not performed"
  ),
  "\n\n",
  "SOFA Dunn tests:\n",
  sep = ""
)

if (nrow(sofa_dunn_output) > 0) {
  print(
    sofa_dunn_output[
      ,
      c(
        "Comparison",
        "AdjustedPValue",
        "Significance"
      )
    ],
    row.names = FALSE
  )
} else {
  cat("  None\n")
}

cat(
  "\nTables saved:\n",
  "  results/tables/cluster_clinical_summary.csv\n",
  "  results/tables/cluster_clinical_manuscript_table.csv\n",
  "  results/tables/cluster_kruskal_wallis_tests.csv\n",
  "  results/tables/cluster_age_dunn_tests.csv\n",
  "  results/tables/cluster_sofa_dunn_tests.csv\n",
  "  results/tables/cluster_dunn_tests_all.csv\n\n",
  "Figures saved:\n",
  "  figures/age_by_cluster.png/.pdf\n",
  "  figures/sofa_by_cluster.png/.pdf\n\n",
  "R object saved:\n",
  "  results/r_objects/cluster_statistics.rds\n",
  "============================================================\n",
  sep = ""
)
