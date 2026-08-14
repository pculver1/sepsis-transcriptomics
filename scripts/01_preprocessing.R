##############################################################
# Project:
# Transcriptomic Characterization of Sepsis Using RNA-seq
#
# Script:
# 01_preprocessing.R
#
# Purpose:
# Import raw RNA-seq counts and GEO metadata, filter
# low-expression genes across all samples, retain sepsis
# samples, normalize counts, and create processed objects.
#
# Author:
# Patrick Culver
#
# Last Updated:
# July 2026
##############################################################

#------------------------------------------------------------
# Reset environment
#------------------------------------------------------------

rm(list = ls())
graphics.off()
set.seed(12345)

#------------------------------------------------------------
# Load required packages
#------------------------------------------------------------

required_packages <- c(
  "DESeq2",
  "data.table"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Required package is not installed: ",
      pkg,
      "\nRun scripts/00_setup.R first."
    )
  }
  
  library(pkg, character.only = TRUE)
}

#------------------------------------------------------------
# File paths
#------------------------------------------------------------

counts_file <- "data/GSE185263_raw_counts.csv"

series_matrix_file <- "data/GSE185263_series_matrix.txt"

validated_dds_file <- "results/r_objects/dds_processed.rds"

recreated_dds_file <- "results/r_objects/dds_from_raw.rds"

#------------------------------------------------------------
# Confirm required files exist
#------------------------------------------------------------

required_files <- c(
  counts_file,
  series_matrix_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  
  stop(
    "The following required files were not found:\n",
    paste(missing_files, collapse = "\n")
  )
}

#------------------------------------------------------------
# Create output directories
#------------------------------------------------------------

output_directories <- c(
  "results/processed_data",
  "results/r_objects",
  "results/logs"
)

for (directory in output_directories) {
  
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
}

#------------------------------------------------------------
# Function for parsing GEO metadata rows
#------------------------------------------------------------

parse_geo_line <- function(line) {
  
  values <- strsplit(line, "\t", fixed = TRUE)[[1]]
  
  values <- values[-1]
  
  values <- gsub('^"|"$', "", values)
  
  values
}

remove_characteristic_label <- function(values) {
  
  sub("^[^:]+:\\s*", "", values)
}

#------------------------------------------------------------
# Read GEO series matrix
#------------------------------------------------------------

cat("Reading GEO series matrix...\n")

series_lines <- readLines(
  series_matrix_file,
  warn = FALSE
)

# Sample titles
title_line <- grep(
  "^!Sample_title",
  series_lines,
  value = TRUE
)

if (length(title_line) != 1) {
  stop("Could not uniquely identify the !Sample_title row.")
}

sample_titles <- parse_geo_line(title_line)

# GSM accession numbers
gsm_line <- grep(
  "^!Sample_geo_accession",
  series_lines,
  value = TRUE
)

if (length(gsm_line) != 1) {
  stop("Could not uniquely identify the !Sample_geo_accession row.")
}

gsm_accessions <- parse_geo_line(gsm_line)

# Repeated characteristics rows
characteristic_lines <- grep(
  "^!Sample_characteristics_ch1",
  series_lines,
  value = TRUE
)

if (length(characteristic_lines) < 8) {
  
  stop(
    "Fewer than eight sample-characteristic rows were found. ",
    "Inspect the GEO series matrix format."
  )
}

characteristics <- lapply(
  characteristic_lines,
  parse_geo_line
)

characteristics <- lapply(
  characteristics,
  remove_characteristic_label
)

#------------------------------------------------------------
# Construct metadata table
#------------------------------------------------------------

metadata_all <- data.frame(
  sample_title = sample_titles,
  gsm_accession = gsm_accessions,
  disease_state = characteristics[[1]],
  age = characteristics[[2]],
  sex = characteristics[[3]],
  collection_location = characteristics[[4]],
  collection_site = characteristics[[5]],
  sofa = characteristics[[6]],
  mortality = characteristics[[7]],
  tissue = characteristics[[8]],
  stringsAsFactors = FALSE
)

# Convert missing-value labels to NA
missing_labels <- c(
  "",
  "NA",
  "N/A",
  "not available",
  "unknown",
  "not reported"
)

metadata_all[metadata_all %in% missing_labels] <- NA

# Convert numeric variables
metadata_all$age <- suppressWarnings(
  as.numeric(metadata_all$age)
)

metadata_all$sofa <- suppressWarnings(
  as.numeric(metadata_all$sofa)
)

# Derive cohort from sample title
metadata_all$cohort <- sub(
  "[0-9].*$",
  "",
  metadata_all$sample_title
)

# Save complete metadata
write.csv(
  metadata_all,
  "results/processed_data/metadata_all_samples.csv",
  row.names = FALSE
)

cat(
  "Metadata records created:",
  nrow(metadata_all),
  "\n"
)

#------------------------------------------------------------
# Read raw count matrix
#------------------------------------------------------------

cat("Reading raw count matrix...\n")

raw_counts <- data.table::fread(
  counts_file,
  data.table = FALSE,
  check.names = FALSE
)

if (ncol(raw_counts) < 2) {
  stop("The raw count file does not contain count columns.")
}

gene_ids <- raw_counts[[1]]

count_data <- raw_counts[, -1, drop = FALSE]

rownames(count_data) <- gene_ids

# Force count columns to numeric
count_data[] <- lapply(
  count_data,
  function(x) as.numeric(as.character(x))
)

count_matrix <- as.matrix(count_data)

storage.mode(count_matrix) <- "integer"

# Remove duplicated gene identifiers, if present
duplicate_gene_ids <- duplicated(rownames(count_matrix))

if (any(duplicate_gene_ids)) {
  
  warning(
    sum(duplicate_gene_ids),
    " duplicated gene identifiers were removed."
  )
  
  count_matrix <- count_matrix[
    !duplicate_gene_ids,
    ,
    drop = FALSE
  ]
}

#------------------------------------------------------------
# Match count columns to metadata
#------------------------------------------------------------

count_sample_ids <- colnames(count_matrix)

title_matches <- sum(
  count_sample_ids %in% metadata_all$sample_title
)

gsm_matches <- sum(
  count_sample_ids %in% metadata_all$gsm_accession
)

if (title_matches == ncol(count_matrix)) {
  
  metadata_all$count_sample_id <- metadata_all$sample_title
  
  matching_identifier <- "Sample title"
  
} else if (gsm_matches == ncol(count_matrix)) {
  
  metadata_all$count_sample_id <- metadata_all$gsm_accession
  
  matching_identifier <- "GSM accession"
  
} else {
  
  stop(
    "The count-matrix column names do not completely match ",
    "either GEO sample titles or GSM accession numbers.\n",
    "Sample-title matches: ",
    title_matches,
    "\nGSM matches: ",
    gsm_matches,
    "\nCount columns: ",
    ncol(count_matrix)
  )
}

cat(
  "Count columns matched using:",
  matching_identifier,
  "\n"
)

# Reorder metadata to match the count matrix
metadata_all <- metadata_all[
  match(
    count_sample_ids,
    metadata_all$count_sample_id
  ),
  ,
  drop = FALSE
]

if (!identical(
  count_sample_ids,
  metadata_all$count_sample_id
)) {
  stop("Metadata could not be aligned to count-matrix columns.")
}

rownames(metadata_all) <- metadata_all$count_sample_id

#------------------------------------------------------------
# Validate raw dimensions
#------------------------------------------------------------

original_gene_count <- nrow(count_matrix)
original_sample_count <- ncol(count_matrix)

cat(
  "Original genes:",
  original_gene_count,
  "\n"
)

cat(
  "Original samples:",
  original_sample_count,
  "\n"
)

#------------------------------------------------------------
# Filter genes across ALL samples
#------------------------------------------------------------

# Retain genes with at least 10 reads in at least 5 samples.
#
# Important:
# Filtering is performed before healthy controls are removed.

minimum_count <- 10
minimum_samples <- 5

keep_gene <- rowSums(
  count_matrix >= minimum_count
) >= minimum_samples

filtered_counts_all <- count_matrix[
  keep_gene,
  ,
  drop = FALSE
]

genes_retained_all <- nrow(filtered_counts_all)

cat(
  "Genes retained after filtering all samples:",
  genes_retained_all,
  "\n"
)

#------------------------------------------------------------
# Identify sepsis samples
#------------------------------------------------------------

disease_clean <- tolower(
  trimws(metadata_all$disease_state)
)

sepsis_samples <- grepl(
  "sepsis",
  disease_clean
)

healthy_samples <- grepl(
  "healthy|control",
  disease_clean
)

cat(
  "Sepsis labels identified:",
  sum(sepsis_samples),
  "\n"
)

cat(
  "Healthy/control labels identified:",
  sum(healthy_samples),
  "\n"
)

if (sum(sepsis_samples) != 348) {
  
  warning(
    "Expected 348 sepsis samples, but identified ",
    sum(sepsis_samples),
    ". Inspect metadata_all_samples.csv."
  )
}

#------------------------------------------------------------
# Retain sepsis samples
#------------------------------------------------------------

counts_sepsis <- filtered_counts_all[
  ,
  sepsis_samples,
  drop = FALSE
]

metadata_sepsis <- metadata_all[
  sepsis_samples,
  ,
  drop = FALSE
]

if (!identical(
  colnames(counts_sepsis),
  rownames(metadata_sepsis)
)) {
  stop("Sepsis metadata and count columns are not aligned.")
}

# Use sample titles as final sample identifiers
colnames(counts_sepsis) <- metadata_sepsis$sample_title
rownames(metadata_sepsis) <- metadata_sepsis$sample_title

metadata_sepsis$sample_id <- metadata_sepsis$sample_title

# Keep primary analysis variables
metadata_deseq <- metadata_sepsis[
  ,
  c(
    "sample_id",
    "age",
    "sofa",
    "cohort",
    "disease_state",
    "sex",
    "mortality",
    "collection_location",
    "collection_site",
    "tissue"
  ),
  drop = FALSE
]

#------------------------------------------------------------
# Save processed metadata and counts
#------------------------------------------------------------

write.csv(
  metadata_sepsis,
  "results/processed_data/metadata_sepsis.csv",
  row.names = FALSE
)

write.csv(
  counts_sepsis,
  "results/processed_data/sepsis_counts_filtered.csv",
  row.names = TRUE
)

#------------------------------------------------------------
# Create DESeq2 object
#------------------------------------------------------------

dds_from_raw <- DESeqDataSetFromMatrix(
  countData = counts_sepsis,
  colData = metadata_deseq,
  design = ~1
)

# Estimate normalization size factors
dds_from_raw <- estimateSizeFactors(
  dds_from_raw
)

#------------------------------------------------------------
# Save normalized counts
#------------------------------------------------------------

normalized_counts <- counts(
  dds_from_raw,
  normalized = TRUE
)

write.csv(
  normalized_counts,
  "results/processed_data/sepsis_counts_normalized.csv",
  row.names = TRUE
)

#------------------------------------------------------------
# Save recreated DESeq2 object
#------------------------------------------------------------

saveRDS(
  dds_from_raw,
  recreated_dds_file
)

#------------------------------------------------------------
# Compare against validated processed object
#------------------------------------------------------------

expected_genes <- 19208
expected_samples <- 348

gene_difference <- nrow(dds_from_raw) - expected_genes
sample_difference <- ncol(dds_from_raw) - expected_samples

validated_object_exists <- file.exists(
  validated_dds_file
)

gene_set_matches <- NA
sample_set_matches <- NA

if (validated_object_exists) {
  
  validated_dds <- readRDS(
    validated_dds_file
  )
  
  gene_set_matches <- setequal(
    rownames(dds_from_raw),
    rownames(validated_dds)
  )
  
  sample_set_matches <- setequal(
    colnames(dds_from_raw),
    colnames(validated_dds)
  )
}

#------------------------------------------------------------
# Write preprocessing summary log
#------------------------------------------------------------

summary_lines <- c(
  "Preprocessing Summary",
  "==================================================",
  paste("Date:", Sys.Date()),
  paste("R version:", R.version.string),
  paste(
    "DESeq2 version:",
    as.character(packageVersion("DESeq2"))
  ),
  "",
  paste("Original genes:", original_gene_count),
  paste("Original samples:", original_sample_count),
  "",
  paste("Minimum count threshold:", minimum_count),
  paste("Minimum number of samples:", minimum_samples),
  paste(
    "Filtering performed before sepsis subsetting:",
    "Yes"
  ),
  "",
  paste(
    "Genes retained after filtering:",
    nrow(dds_from_raw)
  ),
  paste(
    "Sepsis samples retained:",
    ncol(dds_from_raw)
  ),
  "",
  paste("Expected genes:", expected_genes),
  paste("Expected samples:", expected_samples),
  paste("Gene difference:", gene_difference),
  paste("Sample difference:", sample_difference),
  "",
  paste(
    "Validated object found:",
    validated_object_exists
  ),
  paste(
    "Gene identifiers match validated object:",
    gene_set_matches
  ),
  paste(
    "Sample identifiers match validated object:",
    sample_set_matches
  ),
  "",
  paste(
    "Recreated object:",
    recreated_dds_file
  ),
  paste(
    "Validated downstream object:",
    validated_dds_file
  )
)

writeLines(
  summary_lines,
  "results/logs/preprocessing_summary.txt"
)

#------------------------------------------------------------
# Print final summary
#------------------------------------------------------------

cat(
  "\n",
  "----------------------------------------\n",
  "Preprocessing complete\n",
  "----------------------------------------\n",
  "Original samples:      ",
  original_sample_count,
  "\n",
  "Sepsis samples:        ",
  ncol(dds_from_raw),
  "\n\n",
  "Original genes:        ",
  original_gene_count,
  "\n",
  "Genes retained:        ",
  nrow(dds_from_raw),
  "\n",
  "Expected genes:        ",
  expected_genes,
  "\n",
  "Difference:            ",
  gene_difference,
  "\n\n",
  "Recreated object:\n",
  recreated_dds_file,
  "\n\n",
  "Validated downstream object:\n",
  validated_dds_file,
  "\n",
  "----------------------------------------\n",
  sep = ""
)

if (validated_object_exists) {
  
  cat(
    "Gene identifiers match validated object: ",
    gene_set_matches,
    "\n",
    "Sample identifiers match validated object: ",
    sample_set_matches,
    "\n",
    sep = ""
  )
}
