############################################################
# 06_go_enrichment.R
#
# Gene Ontology Biological Process enrichment analysis for
# all pairwise transcriptomic-cluster comparisons.
#
# Input:
#   results/tables/differential_expression/
#     *_all_genes.csv
#     *_upregulated.csv
#     *_downregulated.csv
#
# Output:
#   results/tables/go_enrichment/
#   results/r_objects/go_enrichment/
############################################################


############################################################
# 1. Install required packages if needed
############################################################

cran_packages <- c(
  "dplyr",
  "readr",
  "tibble"
)

bioconductor_packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db"
)

for (package_name in cran_packages) {

  if (!requireNamespace(package_name, quietly = TRUE)) {

    install.packages(
      package_name,
      dependencies = TRUE
    )
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {

  install.packages("BiocManager")
}

for (package_name in bioconductor_packages) {

  if (!requireNamespace(package_name, quietly = TRUE)) {

    BiocManager::install(
      package_name,
      ask = FALSE,
      update = FALSE
    )
  }
}


############################################################
# 2. Load packages
############################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(readr)
library(tibble)


############################################################
# 3. Define project directories
############################################################

differential_expression_directory <- file.path(
  "results",
  "tables",
  "differential_expression"
)

go_table_directory <- file.path(
  "results",
  "tables",
  "go_enrichment"
)

go_object_directory <- file.path(
  "results",
  "r_objects",
  "go_enrichment"
)

go_summary_directory <- file.path(
  go_table_directory,
  "summary"
)


############################################################
# 4. Validate the differential-expression directory
############################################################

if (!dir.exists(differential_expression_directory)) {

  stop(
    paste0(
      "The differential-expression directory was not found:\n",
      normalizePath(
        differential_expression_directory,
        winslash = "/",
        mustWork = FALSE
      ),
      "\n\nRun this script from the Sepsis_RNAseq_project root directory."
    )
  )
}


############################################################
# 5. Create GO enrichment output directories
############################################################

dir.create(
  go_table_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  go_object_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  go_summary_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 6. Analysis settings
############################################################

go_ontology <- "BP"
p_adjust_method <- "BH"
p_value_cutoff <- 0.05
q_value_cutoff <- 0.05
minimum_gene_set_size <- 10
maximum_gene_set_size <- 500


############################################################
# 7. Helper function: validate required columns
############################################################

validate_required_columns <- function(
    data,
    required_columns,
    file_path
) {

  missing_columns <- setdiff(
    required_columns,
    colnames(data)
  )

  if (length(missing_columns) > 0) {

    stop(
      paste0(
        "Required column(s) missing from:\n",
        file_path,
        "\n\nMissing column(s): ",
        paste(missing_columns, collapse = ", "),
        "\n\nAvailable columns: ",
        paste(colnames(data), collapse = ", ")
      )
    )
  }

  invisible(TRUE)
}


############################################################
# 8. Helper function: standardize Entrez identifiers
############################################################

extract_entrez_ids <- function(data) {

  entrez_ids <- data$ENTREZID

  entrez_ids <- as.character(entrez_ids)

  entrez_ids <- trimws(entrez_ids)

  entrez_ids <- entrez_ids[
    !is.na(entrez_ids) &
      entrez_ids != "" &
      entrez_ids != "NA"
  ]

  unique(entrez_ids)
}


############################################################
# 9. Helper function: run GO enrichment safely
############################################################

run_go_enrichment <- function(
    gene_ids,
    universe_ids,
    comparison_label,
    direction_label
) {

  gene_ids <- unique(
    intersect(
      gene_ids,
      universe_ids
    )
  )

  if (length(gene_ids) == 0) {

    warning(
      paste0(
        comparison_label,
        " contained no usable ",
        direction_label,
        " Entrez IDs after filtering."
      )
    )

    return(NULL)
  }

  enrichment_result <- tryCatch(

    enrichGO(
      gene = gene_ids,
      universe = universe_ids,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = go_ontology,
      pAdjustMethod = p_adjust_method,
      pvalueCutoff = p_value_cutoff,
      qvalueCutoff = q_value_cutoff,
      minGSSize = minimum_gene_set_size,
      maxGSSize = maximum_gene_set_size,
      readable = TRUE
    ),

    error = function(error_message) {

      warning(
        paste0(
          "GO enrichment failed for ",
          comparison_label,
          " (",
          direction_label,
          "): ",
          conditionMessage(error_message)
        )
      )

      NULL
    }
  )

  enrichment_result
}


############################################################
# 10. Helper function: convert GO result to a table
############################################################

go_result_to_table <- function(go_result) {

  if (is.null(go_result)) {

    return(
      tibble()
    )
  }

  result_table <- as.data.frame(go_result)

  if (nrow(result_table) == 0) {

    return(
      tibble()
    )
  }

  as_tibble(result_table)
}


############################################################
# 11. Locate all differential-expression comparisons
############################################################

all_gene_files <- list.files(
  differential_expression_directory,
  pattern = "_all_genes\\.csv$",
  full.names = TRUE
)

if (length(all_gene_files) == 0) {

  stop(
    paste0(
      "No files ending in '_all_genes.csv' were found in:\n",
      differential_expression_directory
    )
  )
}

all_gene_files <- sort(all_gene_files)

comparison_labels <- sub(
  "_all_genes\\.csv$",
  "",
  basename(all_gene_files)
)

message(
  "\nDetected ",
  length(comparison_labels),
  " differential-expression comparisons:"
)

message(
  paste0(
    "  - ",
    comparison_labels,
    collapse = "\n"
  )
)


############################################################
# 12. Initialize summary storage
############################################################

go_summary_list <- vector(
  mode = "list",
  length = length(all_gene_files)
)

processing_log_list <- vector(
  mode = "list",
  length = length(all_gene_files)
)


############################################################
# 13. Run GO enrichment for every comparison
############################################################

for (comparison_index in seq_along(all_gene_files)) {

  all_genes_file <- all_gene_files[comparison_index]

  comparison_label <- comparison_labels[comparison_index]

  upregulated_file <- file.path(
    differential_expression_directory,
    paste0(
      comparison_label,
      "_upregulated.csv"
    )
  )

  downregulated_file <- file.path(
    differential_expression_directory,
    paste0(
      comparison_label,
      "_downregulated.csv"
    )
  )

  message(
    "\n============================================================"
  )

  message(
    "Running GO enrichment for: ",
    comparison_label
  )

  message(
    "============================================================"
  )


  ##########################################################
  # Validate comparison files
  ##########################################################

  missing_input_files <- c(
    upregulated_file,
    downregulated_file
  )[
    !file.exists(
      c(
        upregulated_file,
        downregulated_file
      )
    )
  ]

  if (length(missing_input_files) > 0) {

    warning(
      paste0(
        "Skipping ",
        comparison_label,
        " because the following file(s) are missing:\n",
        paste(missing_input_files, collapse = "\n")
      )
    )

    processing_log_list[[comparison_index]] <- tibble(
      Comparison = comparison_label,
      Status = "Skipped",
      Message = paste(
        basename(missing_input_files),
        collapse = "; "
      )
    )

    next
  }


  ##########################################################
  # Read differential-expression tables
  ##########################################################

  all_genes <- read_csv(
    all_genes_file,
    show_col_types = FALSE
  )

  upregulated_genes <- read_csv(
    upregulated_file,
    show_col_types = FALSE
  )

  downregulated_genes <- read_csv(
    downregulated_file,
    show_col_types = FALSE
  )


  ##########################################################
  # Validate required columns
  ##########################################################

  validate_required_columns(
    all_genes,
    "ENTREZID",
    all_genes_file
  )

  validate_required_columns(
    upregulated_genes,
    "ENTREZID",
    upregulated_file
  )

  validate_required_columns(
    downregulated_genes,
    "ENTREZID",
    downregulated_file
  )


  ##########################################################
  # Extract tested-gene universe and DEG identifiers
  ##########################################################

  universe_entrez_ids <- extract_entrez_ids(
    all_genes
  )

  upregulated_entrez_ids <- extract_entrez_ids(
    upregulated_genes
  )

  downregulated_entrez_ids <- extract_entrez_ids(
    downregulated_genes
  )

  upregulated_entrez_ids <- intersect(
    upregulated_entrez_ids,
    universe_entrez_ids
  )

  downregulated_entrez_ids <- intersect(
    downregulated_entrez_ids,
    universe_entrez_ids
  )

  if (length(universe_entrez_ids) == 0) {

    stop(
      paste0(
        "No usable Entrez IDs were found in the tested-gene universe for ",
        comparison_label,
        "."
      )
    )
  }


  ##########################################################
  # Report gene counts
  ##########################################################

  message(
    "All tested rows: ",
    nrow(all_genes)
  )

  message(
    "Unique tested Entrez IDs: ",
    length(universe_entrez_ids)
  )

  message(
    "Upregulated rows: ",
    nrow(upregulated_genes)
  )

  message(
    "Unique upregulated Entrez IDs: ",
    length(upregulated_entrez_ids)
  )

  message(
    "Downregulated rows: ",
    nrow(downregulated_genes)
  )

  message(
    "Unique downregulated Entrez IDs: ",
    length(downregulated_entrez_ids)
  )


  ##########################################################
  # Run GO Biological Process enrichment
  ##########################################################

  go_upregulated <- run_go_enrichment(
    gene_ids = upregulated_entrez_ids,
    universe_ids = universe_entrez_ids,
    comparison_label = comparison_label,
    direction_label = "upregulated"
  )

  go_downregulated <- run_go_enrichment(
    gene_ids = downregulated_entrez_ids,
    universe_ids = universe_entrez_ids,
    comparison_label = comparison_label,
    direction_label = "downregulated"
  )


  ##########################################################
  # Convert GO objects to tables
  ##########################################################

  go_upregulated_table <- go_result_to_table(
    go_upregulated
  )

  go_downregulated_table <- go_result_to_table(
    go_downregulated
  )


  ##########################################################
  # Create comparison-specific output directories
  ##########################################################

  comparison_table_directory <- file.path(
    go_table_directory,
    comparison_label
  )

  comparison_object_directory <- file.path(
    go_object_directory,
    comparison_label
  )

  dir.create(
    comparison_table_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    comparison_object_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )


  ##########################################################
  # Save GO enrichment tables
  ##########################################################

  write_csv(
    go_upregulated_table,
    file.path(
      comparison_table_directory,
      paste0(
        comparison_label,
        "_GO_BP_upregulated.csv"
      )
    )
  )

  write_csv(
    go_downregulated_table,
    file.path(
      comparison_table_directory,
      paste0(
        comparison_label,
        "_GO_BP_downregulated.csv"
      )
    )
  )


  ##########################################################
  # Save GO enrichment R objects
  ##########################################################

  saveRDS(
    go_upregulated,
    file.path(
      comparison_object_directory,
      paste0(
        comparison_label,
        "_GO_BP_upregulated.rds"
      )
    )
  )

  saveRDS(
    go_downregulated,
    file.path(
      comparison_object_directory,
      paste0(
        comparison_label,
        "_GO_BP_downregulated.rds"
      )
    )
  )


  ##########################################################
  # Save Entrez ID lists used for enrichment
  ##########################################################

  write_csv(
    tibble(
      ENTREZID = universe_entrez_ids
    ),
    file.path(
      comparison_table_directory,
      paste0(
        comparison_label,
        "_GO_universe_ENTREZIDs.csv"
      )
    )
  )

  write_csv(
    tibble(
      ENTREZID = upregulated_entrez_ids
    ),
    file.path(
      comparison_table_directory,
      paste0(
        comparison_label,
        "_GO_upregulated_ENTREZIDs.csv"
      )
    )
  )

  write_csv(
    tibble(
      ENTREZID = downregulated_entrez_ids
    ),
    file.path(
      comparison_table_directory,
      paste0(
        comparison_label,
        "_GO_downregulated_ENTREZIDs.csv"
      )
    )
  )


  ##########################################################
  # Record summary information
  ##########################################################

  go_summary_list[[comparison_index]] <- tibble(
    Comparison = comparison_label,
    Tested_Rows = nrow(all_genes),
    Unique_Tested_ENTREZIDs =
      length(universe_entrez_ids),
    Upregulated_Rows =
      nrow(upregulated_genes),
    Unique_Upregulated_ENTREZIDs =
      length(upregulated_entrez_ids),
    Downregulated_Rows =
      nrow(downregulated_genes),
    Unique_Downregulated_ENTREZIDs =
      length(downregulated_entrez_ids),
    Significant_Upregulated_GO_Terms =
      nrow(go_upregulated_table),
    Significant_Downregulated_GO_Terms =
      nrow(go_downregulated_table)
  )

  processing_log_list[[comparison_index]] <- tibble(
    Comparison = comparison_label,
    Status = "Completed",
    Message = "GO enrichment completed successfully."
  )


  ##########################################################
  # Report completion
  ##########################################################

  message(
    "Significant upregulated GO terms: ",
    nrow(go_upregulated_table)
  )

  message(
    "Significant downregulated GO terms: ",
    nrow(go_downregulated_table)
  )

  message(
    "Completed: ",
    comparison_label
  )
}


############################################################
# 14. Combine summary and processing-log tables
############################################################

go_summary <- bind_rows(
  go_summary_list
)

processing_log <- bind_rows(
  processing_log_list
)


############################################################
# 15. Save summary and processing-log tables
############################################################

write_csv(
  go_summary,
  file.path(
    go_summary_directory,
    "GO_BP_enrichment_summary.csv"
  )
)

write_csv(
  processing_log,
  file.path(
    go_summary_directory,
    "GO_BP_processing_log.csv"
  )
)

saveRDS(
  go_summary,
  file.path(
    go_object_directory,
    "GO_BP_enrichment_summary.rds"
  )
)


############################################################
# 16. Save analysis settings
############################################################

analysis_settings <- tibble(
  Setting = c(
    "Ontology",
    "P-value adjustment method",
    "P-value cutoff",
    "Q-value cutoff",
    "Minimum gene-set size",
    "Maximum gene-set size",
    "Gene identifier",
    "Background universe"
  ),
  Value = c(
    go_ontology,
    p_adjust_method,
    as.character(p_value_cutoff),
    as.character(q_value_cutoff),
    as.character(minimum_gene_set_size),
    as.character(maximum_gene_set_size),
    "ENTREZID",
    "Unique mapped genes in each comparison's all_genes table"
  )
)

write_csv(
  analysis_settings,
  file.path(
    go_summary_directory,
    "GO_BP_analysis_settings.csv"
  )
)


############################################################
# 17. Save session information
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    go_summary_directory,
    "GO_BP_session_info.txt"
  )
)


############################################################
# 18. Display final results
############################################################

message(
  "\n============================================================"
)

message(
  "GO Biological Process enrichment analysis complete."
)

message(
  "============================================================"
)

message(
  "GO tables: ",
  normalizePath(
    go_table_directory,
    winslash = "/",
    mustWork = FALSE
  )
)

message(
  "GO R objects: ",
  normalizePath(
    go_object_directory,
    winslash = "/",
    mustWork = FALSE
  )
)

message(
  "Summary files: ",
  normalizePath(
    go_summary_directory,
    winslash = "/",
    mustWork = FALSE
  )
)

print(go_summary)

print(processing_log)
