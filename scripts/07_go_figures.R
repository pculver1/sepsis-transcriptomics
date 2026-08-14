############################################################
# 07_go_figures.R
#
# Publication-quality Gene Ontology Biological Process
# enrichment figures for all transcriptomic-cluster
# comparisons analyzed in 06_go_enrichment.R.
#
# Input:
#   results/r_objects/go_enrichment/
#     <comparison>/
#       *_GO_BP_upregulated.rds
#       *_GO_BP_downregulated.rds
#
# Output:
#   results/figures/go_enrichment/
#     PNG, PDF, and SVG figures for every comparison
#
# Figure design:
#   - Plain black GO-term labels
#   - GeneRatio shown as decimal values
#   - Shared gene-count and adjusted-p-value scales across
#     upregulated/downregulated panels
############################################################


############################################################
# 1. Install required packages if needed
############################################################

cran_packages <- c(
  "dplyr",
  "ggplot2",
  "patchwork",
  "readr",
  "scales",
  "stringr",
  "tibble"
)

for (package_name in cran_packages) {

  if (!requireNamespace(package_name, quietly = TRUE)) {

    install.packages(
      package_name,
      dependencies = TRUE
    )
  }
}


############################################################
# 2. Load packages
############################################################

library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(scales)
library(stringr)
library(tibble)


############################################################
# 3. Define project directories
############################################################

go_object_directory <- file.path(
  "results",
  "r_objects",
  "go_enrichment"
)

go_figure_directory <- file.path(
  "results",
  "figures",
  "go_enrichment"
)

go_figure_summary_directory <- file.path(
  go_figure_directory,
  "summary"
)


############################################################
# 4. Validate input directory
############################################################

if (!dir.exists(go_object_directory)) {

  stop(
    paste0(
      "The GO enrichment object directory was not found:\n",
      normalizePath(
        go_object_directory,
        winslash = "/",
        mustWork = FALSE
      ),
      "\n\nRun this script from the Sepsis_RNAseq_project root ",
      "directory after completing 06_go_enrichment.R."
    )
  )
}


############################################################
# 5. Create output directories
############################################################

dir.create(
  go_figure_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  go_figure_summary_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 6. Figure settings
############################################################

top_terms <- 15

figure_width_two_panel <- 15
figure_width_single_panel <- 8.5
figure_height <- 9

png_resolution <- 600

base_font_size <- 11
axis_label_font_size <- 9.5
title_font_size <- 13
subtitle_font_size <- 10.5
legend_font_size <- 9

go_label_wrap_width <- 42

dot_size_range <- c(2.5, 8)

############################################################
# 7. Helper function: parse GeneRatio values
############################################################

parse_gene_ratio <- function(gene_ratio) {

  gene_ratio <- as.character(gene_ratio)

  vapply(
    gene_ratio,
    FUN.VALUE = numeric(1),
    FUN = function(ratio_value) {

      ratio_parts <- strsplit(
        ratio_value,
        split = "/",
        fixed = TRUE
      )[[1]]

      if (length(ratio_parts) != 2) {

        return(NA_real_)
      }

      numerator <- suppressWarnings(
        as.numeric(ratio_parts[1])
      )

      denominator <- suppressWarnings(
        as.numeric(ratio_parts[2])
      )

      if (
        is.na(numerator) ||
        is.na(denominator) ||
        denominator == 0
      ) {

        return(NA_real_)
      }

      numerator / denominator
    }
  )
}


############################################################
# 8. Helper function: clean and shorten GO labels
############################################################

shorten_go_terms <- function(go_terms) {

  cleaned_terms <- as.character(go_terms)

  cleaned_terms <- str_replace_all(
    cleaned_terms,
    regex(
      "^positive regulation of ",
      ignore_case = TRUE
    ),
    "Positive regulation of "
  )

  cleaned_terms <- str_replace_all(
    cleaned_terms,
    regex(
      "^negative regulation of ",
      ignore_case = TRUE
    ),
    "Negative regulation of "
  )

  cleaned_terms <- str_replace_all(
    cleaned_terms,
    regex(
      "^regulation of ",
      ignore_case = TRUE
    ),
    "Regulation of "
  )

  cleaned_terms <- str_replace_all(
    cleaned_terms,
    regex(
      "leukocyte migration involved in inflammatory response",
      ignore_case = TRUE
    ),
    "Leukocyte migration (inflammatory response)"
  )

  cleaned_terms <- str_replace_all(
    cleaned_terms,
    regex(
      "cellular response to molecule of bacterial origin",
      ignore_case = TRUE
    ),
    "Cellular response to bacterial molecules"
  )

  cleaned_terms <- str_replace_all(
    cleaned_terms,
    regex(
      "response to molecule of bacterial origin",
      ignore_case = TRUE
    ),
    "Response to bacterial molecules"
  )

  cleaned_terms <- str_replace_all(
    cleaned_terms,
    regex(
      "regulation of immune effector process",
      ignore_case = TRUE
    ),
    "Regulation of immune effector processes"
  )

  cleaned_terms <- str_replace_all(
    cleaned_terms,
    regex(
      "production of molecular mediator of immune response",
      ignore_case = TRUE
    ),
    "Production of immune mediators"
  )

  cleaned_terms <- str_replace_all(
    cleaned_terms,
    regex(
      "antigen processing and presentation of peptide antigen",
      ignore_case = TRUE
    ),
    "Peptide antigen processing and presentation"
  )

  cleaned_terms <- str_replace_all(
    cleaned_terms,
    regex(
      "antigen processing and presentation of exogenous peptide antigen",
      ignore_case = TRUE
    ),
    "Exogenous peptide antigen presentation"
  )

  cleaned_terms <- str_replace_all(
    cleaned_terms,
    regex(
      "mitochondrial respiratory chain complex assembly",
      ignore_case = TRUE
    ),
    "Mitochondrial respiratory-chain assembly"
  )

  cleaned_terms <- str_replace_all(
    cleaned_terms,
    regex(
      "aerobic electron transport chain",
      ignore_case = TRUE
    ),
    "Aerobic electron-transport chain"
  )

  cleaned_terms <- str_replace_all(
    cleaned_terms,
    regex(
      "oxidative phosphorylation",
      ignore_case = TRUE
    ),
    "Oxidative phosphorylation"
  )

  cleaned_terms <- str_squish(
    cleaned_terms
  )

  str_wrap(
    cleaned_terms,
    width = go_label_wrap_width
  )
}


############################################################
# 9. Helper function: convert enrichment object to table
############################################################

prepare_go_table <- function(
    go_result,
    top_n = top_terms
) {

  if (is.null(go_result)) {

    return(
      tibble()
    )
  }

  go_table <- as.data.frame(
    go_result
  )

  if (nrow(go_table) == 0) {

    return(
      tibble()
    )
  }

  required_columns <- c(
    "Description",
    "GeneRatio",
    "Count",
    "p.adjust"
  )

  missing_columns <- setdiff(
    required_columns,
    colnames(go_table)
  )

  if (length(missing_columns) > 0) {

    stop(
      paste0(
        "The GO result is missing required column(s): ",
        paste(missing_columns, collapse = ", ")
      )
    )
  }

  go_table <- go_table %>%
    mutate(
      GeneRatioNumeric = parse_gene_ratio(
        GeneRatio
      ),
      AdjustedPValue = as.numeric(
        p.adjust
      ),
      Count = as.numeric(
        Count
      )
    ) %>%
    filter(
      !is.na(AdjustedPValue),
      !is.na(GeneRatioNumeric),
      !is.na(Count)
    ) %>%
    arrange(
      AdjustedPValue,
      desc(GeneRatioNumeric),
      desc(Count)
    ) %>%
    slice_head(
      n = top_n
    ) %>%
    mutate(
      DisplayLabel = shorten_go_terms(
        Description
      )
    )

  if (nrow(go_table) == 0) {

    return(
      tibble()
    )
  }

  go_table <- go_table %>%
    arrange(
      desc(GeneRatioNumeric),
      AdjustedPValue
    ) %>%
    mutate(
      DisplayLabel = factor(
        DisplayLabel,
        levels = rev(
          unique(DisplayLabel)
        )
      )
    )

  go_table
}


############################################################
# 10. Helper function: create one GO dot-plot panel
############################################################

create_go_panel <- function(
    go_table,
    panel_title,
    count_limits,
    p_limits
) {

  if (nrow(go_table) == 0) {

    return(NULL)
  }

  ggplot(
    go_table,
    aes(
      x = GeneRatioNumeric,
      y = DisplayLabel,
      size = Count,
      color = AdjustedPValue
    )
  ) +
    geom_point(
      alpha = 0.9
    ) +
    scale_size_continuous(
      range = dot_size_range,
      limits = count_limits,
      breaks = pretty_breaks(
        n = 4
      ),
      name = "Gene count"
    ) +
    scale_color_viridis_c(
      option = "D",
      direction = -1,
      trans = "reverse",
      limits = p_limits,
      name = "Adjusted\np-value"
    ) +
    scale_x_continuous(
      labels = label_number(
        accuracy = 0.01
      ),
      expand = expansion(
        mult = c(0.02, 0.12)
      )
    ) +
    labs(
      title = panel_title,
      x = "Gene ratio",
      y = NULL
    ) +
    theme_minimal(
      base_size = base_font_size
    ) +
    theme(
      plot.title = element_text(
        size = title_font_size,
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = 10
        )
      ),
      axis.text.y = element_text(
        size = axis_label_font_size,
        lineheight = 1.05,
        hjust = 1,
        color = "black"
      ),
      axis.text.x = element_text(
        size = base_font_size
      ),
      axis.title.x = element_text(
        size = base_font_size,
        face = "bold",
        margin = margin(
          t = 8
        )
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(
        linewidth = 0.35
      ),
      legend.title = element_text(
        size = legend_font_size,
        face = "bold"
      ),
      legend.text = element_text(
        size = legend_font_size
      ),
      legend.position = "right",
      plot.margin = margin(
        t = 12,
        r = 14,
        b = 12,
        l = 12
      )
    )
}


############################################################
# 11. Helper function: format comparison labels
############################################################

format_comparison_label <- function(comparison_label) {

  comparison_label %>%
    str_replace_all(
      "_",
      " "
    ) %>%
    str_replace(
      " vs ",
      " versus "
    )
}


############################################################
# 12. Helper function: create complete GO figure
############################################################

create_go_plot <- function(
    go_up,
    go_down,
    comparison,
    top_n = top_terms
) {

  up_table <- prepare_go_table(
    go_up,
    top_n = top_n
  )

  down_table <- prepare_go_table(
    go_down,
    top_n = top_n
  )

  combined_tables <- bind_rows(
    up_table,
    down_table
  )

  if (nrow(combined_tables) > 0) {
    count_limits <- range(
      combined_tables$Count,
      na.rm = TRUE
    )

    p_limits <- range(
      combined_tables$AdjustedPValue,
      na.rm = TRUE
    )

    if (count_limits[1] == count_limits[2]) {
      count_limits <- c(
        max(0, count_limits[1] - 1),
        count_limits[2] + 1
      )
    }

    if (p_limits[1] == p_limits[2]) {
      p_limits <- c(
        max(.Machine$double.xmin, p_limits[1] * 0.9),
        p_limits[2] * 1.1
      )
    }
  } else {
    count_limits <- c(1, 2)
    p_limits <- c(0.001, 0.05)
  }

  up_panel <- create_go_panel(
    up_table,
    "A. Upregulated genes",
    count_limits = count_limits,
    p_limits = p_limits
  )

  down_panel <- create_go_panel(
    down_table,
    "B. Downregulated genes",
    count_limits = count_limits,
    p_limits = p_limits
  )

  comparison_title <- paste0(
    "GO Biological Process enrichment: ",
    format_comparison_label(
      comparison
    )
  )

  figure_subtitle <- paste0(
    "Top ",
    top_n,
    " enriched terms ranked by adjusted p-value"
  )

  if (
    !is.null(up_panel) &&
    !is.null(down_panel)
  ) {

    combined_figure <- (
      up_panel |
        down_panel
    ) +
      plot_layout(
        guides = "collect"
      ) +
      plot_annotation(
        title = comparison_title,
        subtitle = figure_subtitle,
        theme = theme(
          plot.title = element_text(
            size = 15,
            face = "bold",
            hjust = 0.5
          ),
          plot.subtitle = element_text(
            size = subtitle_font_size,
            hjust = 0.5,
            margin = margin(
              b = 12
            )
          )
        )
      ) &
      theme(
        legend.position = "right"
      )

    return(
      list(
        figure = combined_figure,
        panel_type = "Two-panel",
        up_terms = nrow(up_table),
        down_terms = nrow(down_table),
        width = figure_width_two_panel
      )
    )
  }

  if (!is.null(up_panel)) {

    single_figure <- up_panel +
      plot_annotation(
        title = comparison_title,
        subtitle = figure_subtitle,
        theme = theme(
          plot.title = element_text(
            size = 15,
            face = "bold",
            hjust = 0.5
          ),
          plot.subtitle = element_text(
            size = subtitle_font_size,
            hjust = 0.5,
            margin = margin(
              b = 12
            )
          )
        )
      )

    return(
      list(
        figure = single_figure,
        panel_type = "Upregulated only",
        up_terms = nrow(up_table),
        down_terms = 0,
        width = figure_width_single_panel
      )
    )
  }

  if (!is.null(down_panel)) {

    single_figure <- down_panel +
      plot_annotation(
        title = comparison_title,
        subtitle = figure_subtitle,
        theme = theme(
          plot.title = element_text(
            size = 15,
            face = "bold",
            hjust = 0.5
          ),
          plot.subtitle = element_text(
            size = subtitle_font_size,
            hjust = 0.5,
            margin = margin(
              b = 12
            )
          )
        )
      )

    return(
      list(
        figure = single_figure,
        panel_type = "Downregulated only",
        up_terms = 0,
        down_terms = nrow(down_table),
        width = figure_width_single_panel
      )
    )
  }

  list(
    figure = NULL,
    panel_type = "No significant terms",
    up_terms = 0,
    down_terms = 0,
    width = figure_width_single_panel
  )
}


############################################################
# 13. Locate comparison directories automatically
############################################################

comparison_directories <- list.dirs(
  go_object_directory,
  recursive = FALSE,
  full.names = TRUE
)

comparison_directories <- comparison_directories[
  basename(comparison_directories) != ""
]

comparison_directories <- sort(
  comparison_directories
)

if (length(comparison_directories) == 0) {

  stop(
    paste0(
      "No comparison directories were found in:\n",
      go_object_directory
    )
  )
}

comparison_labels <- basename(
  comparison_directories
)

message(
  "\nDetected ",
  length(comparison_labels),
  " GO enrichment comparisons:"
)

message(
  paste0(
    "  - ",
    comparison_labels,
    collapse = "\n"
  )
)


############################################################
# 14. Initialize figure summary and processing log
############################################################

figure_summary_list <- vector(
  mode = "list",
  length = length(comparison_directories)
)

processing_log_list <- vector(
  mode = "list",
  length = length(comparison_directories)
)


############################################################
# 15. Generate figures for every comparison
############################################################

for (
  comparison_index in seq_along(
    comparison_directories
  )
) {

  comparison_directory <- comparison_directories[
    comparison_index
  ]

  comparison_label <- comparison_labels[
    comparison_index
  ]

  upregulated_rds_file <- file.path(
    comparison_directory,
    paste0(
      comparison_label,
      "_GO_BP_upregulated.rds"
    )
  )

  downregulated_rds_file <- file.path(
    comparison_directory,
    paste0(
      comparison_label,
      "_GO_BP_downregulated.rds"
    )
  )

  message(
    "\n============================================================"
  )

  message(
    "Creating GO figure for: ",
    comparison_label
  )

  message(
    "============================================================"
  )


  ##########################################################
  # Validate RDS files
  ##########################################################

  missing_rds_files <- c(
    upregulated_rds_file,
    downregulated_rds_file
  )[
    !file.exists(
      c(
        upregulated_rds_file,
        downregulated_rds_file
      )
    )
  ]

  if (length(missing_rds_files) > 0) {

    warning(
      paste0(
        "Skipping ",
        comparison_label,
        " because the following RDS file(s) are missing:\n",
        paste(
          missing_rds_files,
          collapse = "\n"
        )
      )
    )

    processing_log_list[[comparison_index]] <- tibble(
      Comparison = comparison_label,
      Status = "Skipped",
      Message = paste(
        basename(missing_rds_files),
        collapse = "; "
      )
    )

    next
  }


  ##########################################################
  # Read enrichment objects
  ##########################################################

  go_upregulated <- readRDS(
    upregulated_rds_file
  )

  go_downregulated <- readRDS(
    downregulated_rds_file
  )


  ##########################################################
  # Create the complete figure
  ##########################################################

  figure_result <- create_go_plot(
    go_up = go_upregulated,
    go_down = go_downregulated,
    comparison = comparison_label,
    top_n = top_terms
  )

  if (is.null(figure_result$figure)) {

    warning(
      paste0(
        "No significant GO terms were available for ",
        comparison_label,
        ". No figure was created."
      )
    )

    figure_summary_list[[comparison_index]] <- tibble(
      Comparison = comparison_label,
      Panel_Type = figure_result$panel_type,
      Displayed_Upregulated_Terms = 0,
      Displayed_Downregulated_Terms = 0,
      PNG_File = NA_character_,
      PDF_File = NA_character_,
      SVG_File = NA_character_
    )

    processing_log_list[[comparison_index]] <- tibble(
      Comparison = comparison_label,
      Status = "No figure",
      Message = "No significant GO terms were available."
    )

    next
  }


  ##########################################################
  # Define output file names
  ##########################################################

  png_file <- file.path(
    go_figure_directory,
    paste0(
      comparison_label,
      "_GO_BP.png"
    )
  )

  pdf_file <- file.path(
    go_figure_directory,
    paste0(
      comparison_label,
      "_GO_BP.pdf"
    )
  )

  svg_file <- file.path(
    go_figure_directory,
    paste0(
      comparison_label,
      "_GO_BP.svg"
    )
  )


  ##########################################################
  # Save PNG figure
  ##########################################################

  ggsave(
    filename = png_file,
    plot = figure_result$figure,
    width = figure_result$width,
    height = figure_height,
    units = "in",
    dpi = png_resolution,
    bg = "white"
  )


  ##########################################################
  # Save PDF figure
  ##########################################################

  ggsave(
    filename = pdf_file,
    plot = figure_result$figure,
    width = figure_result$width,
    height = figure_height,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )


  ##########################################################
  # Save SVG figure
  ##########################################################

  ggsave(
    filename = svg_file,
    plot = figure_result$figure,
    width = figure_result$width,
    height = figure_height,
    units = "in",
    device = "svg",
    bg = "white"
  )


  ##########################################################
  # Record figure summary
  ##########################################################

  figure_summary_list[[comparison_index]] <- tibble(
    Comparison = comparison_label,
    Panel_Type = figure_result$panel_type,
    Displayed_Upregulated_Terms =
      figure_result$up_terms,
    Displayed_Downregulated_Terms =
      figure_result$down_terms,
    PNG_File = normalizePath(
      png_file,
      winslash = "/",
      mustWork = FALSE
    ),
    PDF_File = normalizePath(
      pdf_file,
      winslash = "/",
      mustWork = FALSE
    ),
    SVG_File = normalizePath(
      svg_file,
      winslash = "/",
      mustWork = FALSE
    )
  )

  processing_log_list[[comparison_index]] <- tibble(
    Comparison = comparison_label,
    Status = "Completed",
    Message = paste0(
      figure_result$panel_type,
      " GO figure created successfully."
    )
  )

  message(
    "Figure type: ",
    figure_result$panel_type
  )

  message(
    "Displayed upregulated terms: ",
    figure_result$up_terms
  )

  message(
    "Displayed downregulated terms: ",
    figure_result$down_terms
  )

  message(
    "Saved PNG, PDF, and SVG files."
  )
}


############################################################
# 16. Combine summary and processing-log tables
############################################################

figure_summary <- bind_rows(
  figure_summary_list
)

processing_log <- bind_rows(
  processing_log_list
)


############################################################
# 17. Save figure summary and processing log
############################################################

write_csv(
  figure_summary,
  file.path(
    go_figure_summary_directory,
    "GO_BP_figure_summary.csv"
  )
)

write_csv(
  processing_log,
  file.path(
    go_figure_summary_directory,
    "GO_BP_figure_processing_log.csv"
  )
)


############################################################
# 18. Save figure settings
############################################################

figure_settings <- tibble(
  Setting = c(
    "Top terms per direction",
    "GO label wrap width",
    "PNG resolution",
    "Two-panel width",
    "Single-panel width",
    "Figure height",
    "Dot-size minimum",
    "Dot-size maximum",
    "Adjusted p-value color scale",
    "Figure formats"
  ),
  Value = c(
    as.character(top_terms),
    as.character(go_label_wrap_width),
    as.character(png_resolution),
    as.character(figure_width_two_panel),
    as.character(figure_width_single_panel),
    as.character(figure_height),
    as.character(dot_size_range[1]),
    as.character(dot_size_range[2]),
    "viridis option D, reversed",
    "PNG; PDF; SVG"
  )
)

write_csv(
  figure_settings,
  file.path(
    go_figure_summary_directory,
    "GO_BP_figure_settings.csv"
  )
)


############################################################
# 19. Save session information
############################################################

capture.output(
  sessionInfo(),
  file = file.path(
    go_figure_summary_directory,
    "GO_BP_figure_session_info.txt"
  )
)


############################################################
# 20. Display final output
############################################################

message(
  "\n============================================================"
)

message(
  "GO enrichment figure generation complete."
)

message(
  "============================================================"
)

message(
  "Figures saved to: ",
  normalizePath(
    go_figure_directory,
    winslash = "/",
    mustWork = FALSE
  )
)

message(
  "Figure summaries saved to: ",
  normalizePath(
    go_figure_summary_directory,
    winslash = "/",
    mustWork = FALSE
  )
)

print(
  figure_summary
)

print(
  processing_log
)
