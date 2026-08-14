##############################################################
# Project:
# Transcriptomic Characterization of Sepsis Using RNA-seq
#
# Script:
# 00_setup.R
#
# Purpose:
# Install and load required packages, create project folders,
# and initialize the analysis environment.
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
# Package lists
#------------------------------------------------------------

cran_packages <- c(
  "tidyverse",
  "data.table",
  "matrixStats",
  "pheatmap",
  "RColorBrewer",
  "ggrepel",
  "factoextra",
  "openxlsx"
)

bioconductor_packages <- c(
  "DESeq2",
  "clusterProfiler",
  "org.Hs.eg.db"
)

#------------------------------------------------------------
# Install CRAN packages
#------------------------------------------------------------

for (pkg in cran_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  
  library(pkg, character.only = TRUE)
}

#------------------------------------------------------------
# Install Bioconductor packages
#------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in bioconductor_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  
  library(pkg, character.only = TRUE)
}

#------------------------------------------------------------
# Create project folders
#------------------------------------------------------------

project_directories <- c(
  "data",
  "scripts",
  "results",
  "results/processed_data",
  "results/r_objects",
  "results/tables",
  "results/logs",
  "figures",
  "manuscript"
)

for (directory in project_directories) {
  
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
}

cat(
  "\n",
  "Setup complete.\n",
  "R version: ", R.version.string, "\n",
  "DESeq2 version: ", as.character(packageVersion("DESeq2")), "\n",
  sep = ""
)