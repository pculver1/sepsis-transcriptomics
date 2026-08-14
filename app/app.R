############################################################
# Sepsis Transcriptomics Dashboard
############################################################

library(shiny)
library(bslib)
library(plotly)
library(DT)
library(ggplot2)

# ==========================================================
# LOCATE PROJECT ROOT
# ==========================================================

# This lets the app run whether R's working directory is:
#   Sepsis_RNAseq_project/
# or
#   Sepsis_RNAseq_project/app/

if (
  file.exists(
    "results/tables/PCA_coordinates.csv"
  )
) {

  project_root <- "."

} else if (
  file.exists(
    "../results/tables/PCA_coordinates.csv"
  )
) {

  project_root <- ".."

} else {

  stop(
    paste0(
      "Could not locate the Sepsis_RNAseq_project results folder.\n",
      "Open the project in RStudio and run the app from either the ",
      "project root or the app folder."
    )
  )
}

# ==========================================================
# LOAD FINALIZED PROJECT OUTPUTS
# ==========================================================

pca <- read.csv(
  file.path(
    project_root,
    "results",
    "tables",
    "PCA_coordinates.csv"
  ),
  stringsAsFactors = FALSE
)

clinical <- read.csv(
  file.path(
    project_root,
    "results",
    "tables",
    "cluster_clinical_summary.csv"
  ),
  stringsAsFactors = FALSE
)

kw <- read.csv(
  file.path(
    project_root,
    "results",
    "tables",
    "cluster_kruskal_wallis_tests.csv"
  ),
  stringsAsFactors = FALSE
)

age_dunn <- read.csv(
  file.path(
    project_root,
    "results",
    "tables",
    "cluster_age_dunn_tests.csv"
  ),
  stringsAsFactors = FALSE
)

sofa_dunn <- read.csv(
  file.path(
    project_root,
    "results",
    "tables",
    "cluster_sofa_dunn_tests.csv"
  ),
  stringsAsFactors = FALSE
)


# ==========================================================
# PREPARE DATA
# ==========================================================

cluster_levels <- c(
  "Cluster 1",
  "Cluster 2",
  "Cluster 3",
  "Cluster 4"
)

pca$Cluster <- factor(
  pca$Cluster,
  levels = cluster_levels
)

clinical$Cluster <- factor(
  clinical$Cluster,
  levels = cluster_levels
)

cluster_sizes <- table(pca$Cluster)

# Differential-expression comparison definitions
de_comparisons <- c(
  "Cluster 2 vs Cluster 1" = "Cluster_2_vs_Cluster_1",
  "Cluster 3 vs Cluster 1" = "Cluster_3_vs_Cluster_1",
  "Cluster 4 vs Cluster 1" = "Cluster_4_vs_Cluster_1",
  "Cluster 3 vs Cluster 2" = "Cluster_3_vs_Cluster_2",
  "Cluster 4 vs Cluster 2" = "Cluster_4_vs_Cluster_2",
  "Cluster 4 vs Cluster 3" = "Cluster_4_vs_Cluster_3"
)

# ==========================================================
# HELPER FUNCTIONS
# ==========================================================

stat_card <- function(title, value, subtitle = NULL) {
  div(
    class = "stat-card",
    div(class = "stat-card-title", title),
    div(class = "stat-card-value", value),
    if (!is.null(subtitle)) {
      div(class = "stat-card-subtitle", subtitle)
    }
  )
}

load_de_table <- function(comparison_key) {

  file_path <- file.path(
    project_root,
    "results",
    "tables",
    "differential_expression",
    paste0(
      comparison_key,
      "_all_genes.csv"
    )
  )

  if (!file.exists(file_path)) {
    stop(
      "Differential-expression file not found: ",
      file_path
    )
  }

  read.csv(
    file_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

format_p <- function(x) {

  ifelse(
    is.na(x),
    NA_character_,
    ifelse(
      x < 0.001,
      format(
        x,
        scientific = TRUE,
        digits = 3
      ),
      formatC(
        x,
        format = "f",
        digits = 4
      )
    )
  )
}


parse_ratio <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)

  vapply(
    parts,
    function(z) {
      if (length(z) != 2) return(NA_real_)
      numerator <- suppressWarnings(as.numeric(z[1]))
      denominator <- suppressWarnings(as.numeric(z[2]))

      if (is.na(numerator) || is.na(denominator) || denominator == 0) {
        return(NA_real_)
      }

      numerator / denominator
    },
    FUN.VALUE = numeric(1)
  )
}

ratio_denominator <- function(x) {
  if (length(x) == 0) return(NA_integer_)

  parts <- strsplit(
    as.character(x[1]),
    "/",
    fixed = TRUE
  )[[1]]

  if (length(parts) != 2) return(NA_integer_)

  suppressWarnings(as.integer(parts[2]))
}

wrap_go_term <- function(x, width = 42) {
  vapply(
    as.character(x),
    function(term) {
      paste(
        strwrap(term, width = width),
        collapse = "\n"
      )
    },
    FUN.VALUE = character(1)
  )
}

load_go_table <- function(comparison_key, direction) {

  direction_file <- if (
    direction == "Upregulated"
  ) {
    "upregulated"
  } else {
    "downregulated"
  }

  file_path <- file.path(
    project_root,
    "results",
    "tables",
    "go_enrichment",
    comparison_key,
    paste0(
      comparison_key,
      "_GO_BP_",
      direction_file,
      ".csv"
    )
  )

  if (!file.exists(file_path)) {
    stop(
      "GO enrichment file not found: ",
      file_path
    )
  }

  read.csv(
    file_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

# ==========================================================
# USER INTERFACE
# ==========================================================

ui <- page_navbar(

  title = "Sepsis Transcriptomics Dashboard",

  theme = bs_theme(
    version = 5,
    bootswatch = "flatly"
  ),

  header = tags$head(
    tags$style(
      HTML("
        body {
          background-color: #f7f9fb;
        }

        .navbar {
          margin-bottom: 0;
        }

        .page-title {
          margin-top: 1rem;
          margin-bottom: 0.25rem;
        }

        .page-subtitle {
          margin-bottom: 1.25rem;
          color: #4f5b66;
        }

        .stat-grid {
          display: grid;
          grid-template-columns: repeat(4, minmax(0, 1fr));
          gap: 1rem;
          margin-bottom: 1.25rem;
        }

        .stat-grid-three {
          display: grid;
          grid-template-columns: repeat(3, minmax(0, 1fr));
          gap: 1rem;
          margin-bottom: 1.25rem;
        }

        .stat-card {
          background: white;
          border: 1px solid #d9e0e6;
          border-radius: 10px;
          padding: 1rem 1.1rem;
          box-shadow: 0 1px 2px rgba(0,0,0,0.04);
        }

        .stat-card-title {
          font-size: 0.95rem;
          color: #5f6b75;
          margin-bottom: 0.25rem;
        }

        .stat-card-value {
          font-size: 1.75rem;
          font-weight: 700;
          line-height: 1.1;
        }

        .stat-card-subtitle {
          font-size: 0.8rem;
          color: #74808a;
          margin-top: 0.3rem;
        }

        .dashboard-card {
          background: white;
          border: 1px solid #d9e0e6;
          border-radius: 10px;
          padding: 1rem;
          box-shadow: 0 1px 2px rgba(0,0,0,0.04);
          margin-bottom: 1.25rem;
          overflow: visible;
        }

        .dashboard-card h3 {
          margin-top: 0;
          margin-bottom: 0.5rem;
          font-size: 1.15rem;
        }

        .pca-wrap {
          width: 100%;
          min-height: 720px;
        }

        .clinical-grid {
          display: grid;
          grid-template-columns: minmax(0, 2fr) minmax(280px, 1fr);
          gap: 1rem;
          align-items: stretch;
        }

        .clinical-plot-wrap {
          min-height: 560px;
        }

        .clinical-stat-card {
          background: #fbfcfd;
          border: 1px solid #dfe5ea;
          border-radius: 8px;
          padding: 1rem;
        }

        .clinical-grid > .clinical-stat-card {
          height: 100%;
        }

        .volcano-wrap {
          width: 100%;
          min-height: 650px;
        }

        .controls-row {
          display: flex;
          gap: 1rem;
          flex-wrap: wrap;
          align-items: end;
          margin-bottom: 1rem;
        }

        .controls-row .form-group {
          margin-bottom: 0;
        }

        .author-affiliation {
          margin-top: -0.5rem;
          margin-bottom: 1.5rem;
          color: #5f6b76;
        }

        .author-name {
          color: #2f3942;
          font-size: 1.05rem;
          font-weight: 600;
        }

        @media (max-width: 1000px) {

          .stat-grid,
          .stat-grid-three {
            grid-template-columns: repeat(2, minmax(0, 1fr));
          }

          .clinical-grid {
            grid-template-columns: 1fr;
          }

          .workflow-grid {
            grid-template-columns: repeat(2, minmax(0, 1fr)) !important;
          }
        }

        @media (max-width: 650px) {

          .stat-grid,
          .stat-grid-three,
          .workflow-grid {
            grid-template-columns: 1fr !important;
          }
        }
      ")
    )
  ),

  # --------------------------------------------------------
  # OVERVIEW
  # --------------------------------------------------------

  nav_panel(
    "Overview",

    div(
      class = "page-title",
      h2("Transcriptomic Characterization of Sepsis")
    ),

    div(
      class = "page-subtitle",
      p(
        "Interactive exploration of transcriptomic heterogeneity, clinical ",
        "characteristics, differential gene expression, and biological ",
        "process enrichment among patients with sepsis."
      )
    ),

    div(
      class = "author-affiliation",
      tags$div(class = "author-name", "Patrick Culver"),
      tags$div("Arizona State University"),
      tags$div("New College of Interdisciplinary Arts and Sciences")
    ),

    div(
      class = "stat-grid",

      stat_card(
        "Sepsis patients",
        format(nrow(pca), big.mark = ","),
        "Whole-blood RNA-seq samples"
      ),

      stat_card(
        "Retained genes",
        "19,208",
        "After low-expression filtering"
      ),

      stat_card(
        "Transcriptomic clusters",
        length(unique(na.omit(pca$Cluster))),
        "Ward hierarchical clustering"
      ),

      stat_card(
        "Pairwise comparisons",
        length(de_comparisons),
        "Differential expression and GO analysis"
      )
    ),

    div(
      class = "dashboard-card",

      h3("Study Overview"),

      p(
        "Whole-blood RNA-seq data from GEO dataset GSE185263 were used ",
        "to investigate transcriptomic heterogeneity among 348 patients ",
        "with sepsis. After low-expression filtering, 19,208 genes were ",
        "retained for downstream analysis."
      ),

      p(
        "Variance-stabilized expression values were used for principal ",
        "component analysis and hierarchical clustering. Four transcriptomic ",
        "clusters were identified and then compared using patient age, SOFA ",
        "score, differential gene expression, and Gene Ontology Biological ",
        "Process enrichment."
      )
    ),

    div(
      class = "dashboard-card",

      h3("Analysis Workflow"),

      div(
        class = "workflow-grid",
        style = paste0(
          "display:grid;grid-template-columns:repeat(5,minmax(0,1fr));",
          "gap:0.75rem;text-align:center;"
        ),

        div(class = "clinical-stat-card",
            strong("1. RNA-seq Data"),
            p("GSE185263 whole-blood counts and clinical metadata")),

        div(class = "clinical-stat-card",
            strong("2. Preprocessing"),
            p("Low-expression filtering and variance-stabilizing transformation")),

        div(class = "clinical-stat-card",
            strong("3. Clustering"),
            p("PCA visualization and Ward hierarchical clustering")),

        div(class = "clinical-stat-card",
            strong("4. Differential Expression"),
            p("DESeq2 with apeglm log2 fold-change shrinkage")),

        div(class = "clinical-stat-card",
            strong("5. Functional Enrichment"),
            p("GO Biological Process enrichment of significant DEGs"))
      )
    ),

    div(
      class = "dashboard-card",

      h3("Key Findings"),

      div(
        class = "clinical-grid",

        div(
          class = "clinical-stat-card",
          h4("Clinical Characteristics"),
          p(
            strong("Age differed among clusters. "),
            "Cluster 4 was substantially younger than Clusters 1, 2, and 3, ",
            "while no significant age differences were detected among ",
            "Clusters 1, 2, and 3."
          ),
          p(
            strong("SOFA scores also differed overall among clusters. "),
            "Cluster 2 had significantly higher SOFA scores than Clusters ",
            "1 and 3 after Dunn testing with Benjamini-Hochberg correction. ",
            "The overall SOFA effect size was small."
          )
        ),

        div(
          class = "clinical-stat-card",
          h4("Transcriptomic Differences"),
          p(
            "The four clusters showed substantial transcriptomic separation. ",
            "For Cluster 2 versus Cluster 1, 3,916 genes met the differential ",
            "expression criteria of adjusted p-value < 0.05 and absolute ",
            "apeglm-shrunken log2 fold change >= 1."
          ),
          p(
            "The dashboard provides all six pairwise cluster comparisons so ",
            "gene-level differences can be explored beyond the primary ",
            "Cluster 2 versus Cluster 1 comparison highlighted in the manuscript."
          )
        )
      ),

      br(),

      div(
        class = "clinical-stat-card",
        h4("Biological Interpretation"),
        p(
          "Cluster 2 showed a distinct immune-response profile relative to Cluster 1. ",
          "Genes with higher expression in Cluster 2 were enriched for inflammatory ",
          "and innate immune processes, including regulation of inflammatory response, ",
          "chemotaxis, and defense response to bacterium. Genes with lower expression ",
          "were enriched for adaptive immune processes, including antigen receptor ",
          "signaling, lymphocyte activation, and T-cell differentiation. Together, ",
          "these findings suggest an imbalance between innate inflammatory activity ",
          "and adaptive immune function rather than a uniform increase in immune activity."
        ),
        p(
          "This molecular pattern was accompanied by significantly higher SOFA scores ",
          "in Cluster 2 than in Clusters 1 and 3, although the overall SOFA effect size ",
          "was small. Cluster 4 was the most age-distinct group and was substantially ",
          "younger than the other clusters, indicating that age may contribute to some ",
          "of its transcriptomic separation. Overall, the results support biologically ",
          "meaningful heterogeneity among patients with sepsis while emphasizing that ",
          "these clusters should not be interpreted as established clinical endotypes ",
          "without external validation."
        )
      )
    )
  ),

  # --------------------------------------------------------
  # CLUSTERS
  # --------------------------------------------------------

  nav_panel(
    "Clusters",

    div(
      class = "page-title",
      h2("Transcriptomic Clusters")
    ),

    div(
      class = "page-subtitle",
      p(
        "Explore the four transcriptomic clusters identified ",
        "by hierarchical clustering of variance-stabilized ",
        "RNA-seq expression data."
      )
    ),

    div(
      class = "stat-grid",

      stat_card(
        "Cluster 1",
        as.integer(cluster_sizes[["Cluster 1"]])
      ),

      stat_card(
        "Cluster 2",
        as.integer(cluster_sizes[["Cluster 2"]])
      ),

      stat_card(
        "Cluster 3",
        as.integer(cluster_sizes[["Cluster 3"]])
      ),

      stat_card(
        "Cluster 4",
        as.integer(cluster_sizes[["Cluster 4"]])
      )
    ),

    div(
      class = "dashboard-card",

      h3("Principal Component Analysis"),

      div(
        class = "pca-wrap",

        plotlyOutput(
          "pca_plot",
          width = "100%",
          height = "720px"
        )
      )
    ),

    div(
      class = "dashboard-card",

      h3("Clinical Characteristics by Cluster"),

      selectInput(
        inputId = "clinical_variable",
        label = "Clinical variable",
        choices = c("Age", "SOFA"),
        selected = "Age",
        width = "260px"
      ),

      div(
        class = "clinical-grid",

        div(
          class = "clinical-plot-wrap",

          plotlyOutput(
            "clinical_plot",
            width = "100%",
            height = "560px"
          )
        ),

        div(
          class = "clinical-stat-card",

          h4("Overall Statistical Test"),

          uiOutput("kw_result")
        )
      ),

      h4("Summary Statistics"),

      DTOutput("clinical_summary_table"),

      br(),

      h4("Pairwise Dunn Tests"),

      DTOutput("dunn_table")
    )
  ),

  # --------------------------------------------------------
  # DIFFERENTIAL EXPRESSION
  # --------------------------------------------------------

  nav_panel(
    "Differential Expression",

    div(
      class = "page-title",
      h2("Differential Expression")
    ),

    div(
      class = "page-subtitle",
      p(
        "Explore the six pairwise transcriptomic-cluster ",
        "comparisons generated with DESeq2 and apeglm shrinkage."
      )
    ),

    div(
      class = "dashboard-card",

      div(
        class = "controls-row",

        selectInput(
          inputId = "de_comparison",
          label = "Cluster comparison",
          choices = de_comparisons,
          selected = "Cluster_2_vs_Cluster_1",
          width = "320px"
        ),

        selectInput(
          inputId = "de_table_filter",
          label = "Genes shown in table",
          choices = c(
            "Significant DEGs" = "significant",
            "All tested genes" = "all",
            "Upregulated only" = "up",
            "Downregulated only" = "down"
          ),
          selected = "significant",
          width = "220px"
        )
      ),

      uiOutput("de_comparison_title")
    ),

    uiOutput("de_stat_cards"),

    div(
      class = "dashboard-card",

      h3("Interactive Volcano Plot"),

      p(
        "Significant DEGs are defined as adjusted p-value < 0.05 ",
        "and absolute apeglm-shrunken log2 fold change >= 1."
      ),

      div(
        class = "volcano-wrap",

        plotlyOutput(
          "volcano_plot",
          width = "100%",
          height = "650px"
        )
      )
    ),

    div(
      class = "dashboard-card",

      h3("Gene-Level Results"),

      p(
        "Use the search box to locate genes by symbol, Ensembl ID, ",
        "gene name, or Entrez ID."
      ),

      DTOutput("de_table")
    )
  ),

  # --------------------------------------------------------
  # GO ENRICHMENT
  # --------------------------------------------------------

  nav_panel(
    "GO Enrichment",

    div(
      class = "page-title",
      h2("Gene Ontology Enrichment")
    ),

    div(
      class = "page-subtitle",
      p(
        "Explore Gene Ontology Biological Process enrichment ",
        "for upregulated and downregulated genes from each ",
        "pairwise transcriptomic-cluster comparison."
      )
    ),

    div(
      class = "dashboard-card",

      div(
        class = "controls-row",

        selectInput(
          inputId = "go_comparison",
          label = "Cluster comparison",
          choices = de_comparisons,
          selected = "Cluster_2_vs_Cluster_1",
          width = "320px"
        ),

        selectInput(
          inputId = "go_direction",
          label = "Gene direction",
          choices = c(
            "Upregulated",
            "Downregulated"
          ),
          selected = "Upregulated",
          width = "220px"
        ),

        selectInput(
          inputId = "go_top_n",
          label = "Terms displayed",
          choices = c(
            5,
            10,
            15,
            20,
            25
          ),
          selected = 15,
          width = "180px"
        )
      ),

      uiOutput("go_selection_title")
    ),

    uiOutput("go_stat_cards"),

    div(
      class = "dashboard-card",

      h3("Interactive GO Biological Process Plot"),

      p(
        "Terms are ranked by Benjamini-Hochberg adjusted ",
        "p-value. Dot position represents GeneRatio, dot size ",
        "represents gene count, and color represents adjusted p-value."
      ),

      uiOutput("go_plot_message"),

      plotlyOutput(
        "go_plot",
        width = "100%",
        height = "720px"
      )
    ),

    div(
      class = "dashboard-card",

      h3("GO Enrichment Results"),

      p(
        "Search the complete significant enrichment results ",
        "for the selected comparison and direction."
      ),

      DTOutput("go_table")
    )
  )
  ,

  nav_panel(
    "About",

    div(
      class = "page-title",
      h2("About This Project")
    ),

    div(
      class = "page-subtitle",
      p(
        "Methods, data sources, software, and references supporting ",
        "the interactive sepsis transcriptomics dashboard."
      )
    ),

    div(
      class = "dashboard-card",

      h3("Author & Affiliation"),

      p(
        strong("Patrick Culver"),
        tags$br(),
        "Arizona State University",
        tags$br(),
        "New College of Interdisciplinary Arts and Sciences"
      )
    ),

    div(
      class = "dashboard-card",

      h3("Methods Summary"),

      p(
        "Whole-blood RNA-seq count data and associated clinical metadata ",
        "were obtained from GEO accession GSE185263. Healthy control samples ",
        "were excluded, leaving 348 sepsis samples for analysis."
      ),

      p(
        "Genes were retained if they had a count of at least 10 in at least ",
        "five samples, resulting in 19,208 genes. Variance-stabilizing ",
        "transformation was applied before PCA and hierarchical clustering. ",
        "PCA used the 500 most variable VST-transformed genes, while ",
        "hierarchical clustering used all 19,208 retained VST-transformed genes ",
        "with Euclidean distance, ward.D2 linkage, and k = 4 clusters."
      ),

      p(
        "Age and SOFA score were compared among transcriptomic clusters using ",
        "Kruskal-Wallis tests. Significant overall tests were followed by ",
        "pairwise Dunn tests with Benjamini-Hochberg correction."
      ),

      p(
        "Differential expression analysis was performed with DESeq2 using ",
        "cluster as the design variable. Log2 fold changes were shrunk using ",
        "apeglm. Genes were considered differentially expressed when the ",
        "Benjamini-Hochberg adjusted p-value was < 0.05 and the absolute ",
        "shrunken log2 fold change was >= 1."
      ),

      p(
        "Gene Ontology Biological Process enrichment was performed separately ",
        "for upregulated and downregulated genes using clusterProfiler. ",
        "Unique mapped genes from the corresponding tested-gene table were ",
        "used as the background universe, with Benjamini-Hochberg adjusted ",
        "p-value and q-value thresholds of 0.05."
      )
    ),

    div(
      class = "dashboard-card",

      h3("Data and Software"),

      p(
        strong("Primary dataset: "),
        "GSE185263, Predicting Sepsis Severity at First Clinical Presentation: ",
        "The Role of Endotypes and Mechanistic Signatures."
      ),

      p(
        strong("Core software: "),
        "R; DESeq2; apeglm; clusterProfiler; ggplot2; plotly; shiny; bslib; DT."
      ),

      p(
        strong("Interactive application: "),
        "Built in R with Shiny and Plotly. The dashboard reads finalized ",
        "analysis outputs rather than rerunning the full RNA-seq pipeline ",
        "during use."
      )
    ),

    div(
      class = "dashboard-card",

      h3("References"),

      tags$ol(

        tags$li(
          "Baghela AS, et al. Predicting sepsis severity at first clinical ",
          "presentation: the role of endotypes and mechanistic signatures ",
          "(GSE185263) [dataset]. Gene Expression Omnibus; 2022."
        ),

        tags$li(
          "Love MI, Huber W, Anders S. Moderated estimation of fold change ",
          "and dispersion for RNA-seq data with DESeq2. Genome Biology. ",
          "2014;15(12):550."
        ),

        tags$li(
          "Zhu A, Ibrahim JG, Love MI. Heavy-tailed prior distributions for ",
          "sequence count data: removing the noise and preserving large ",
          "differences. Bioinformatics. 2019;35:2084-2092."
        ),

        tags$li(
          "Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for ",
          "comparing biological themes among gene clusters. OMICS. ",
          "2012;16(5):284-287."
        ),

        tags$li(
          "Chenoweth JG, et al. Sepsis endotypes identified by host gene ",
          "expression across global cohorts. Communications Medicine. ",
          "2024;4(1):120."
        ),

        tags$li(
          "Yang JO, et al. Whole blood transcriptomics identifies subclasses ",
          "of pediatric septic shock. Critical Care. 2023;27(1):486."
        )
      )
    ),

    div(
      class = "dashboard-card",

      h3("Project Scope"),

      p(
        "This project is an exploratory transcriptomic analysis intended to ",
        "characterize molecular heterogeneity within sepsis. The identified ",
        "clusters are transcriptomic subgroups and should not be interpreted ",
        "as clinically validated sepsis endotypes."
      ),

      p(
        "The analysis is based on a public multicohort dataset. Cohort was not ",
        "included as a covariate in the clustering or differential expression ",
        "models, so some observed transcriptomic differences may reflect ",
        "cohort-associated variation."
      )
    )
  )

)

# ==========================================================
# SERVER
# ==========================================================

server <- function(input, output, session) {

  # --------------------------------------------------------
  # PCA plot
  # --------------------------------------------------------

  output$pca_plot <- renderPlotly({

    pca_colors <- c(
      "Cluster 1" = "#F8766D",
      "Cluster 2" = "#7CAE00",
      "Cluster 3" = "#00BFC4",
      "Cluster 4" = "#C77CFF"
    )

    p <- ggplot(
      pca,
      aes(
        x = PC1,
        y = PC2,
        color = Cluster,
        text = paste0(
          "Sample: ", Sample,
          "<br>Cluster: ", Cluster,
          "<br>Cohort: ", Cohort,
          "<br>Age: ", Age,
          "<br>SOFA: ",
          ifelse(
            is.na(SOFA),
            "Missing",
            SOFA
          )
        )
      )
    ) +

      stat_ellipse(
        aes(
          x = PC1,
          y = PC2,
          group = Cluster,
          color = Cluster
        ),
        type = "norm",
        level = 0.95,
        linewidth = 0.8,
        linetype = "dashed",
        show.legend = FALSE,
        inherit.aes = FALSE
      ) +

      geom_point(
        size = 3,
        alpha = 0.85,
        shape = 16
      ) +

      scale_color_manual(
        values = pca_colors,
        name = "Transcriptomic cluster"
      ) +

      labs(
        title = "Principal Component Analysis of Sepsis Samples",
        subtitle = paste0(
          "PCA: top 500 variable genes; ",
          "clusters: all 19,208 VST genes"
        ),
        x = "PC1 (23.7% variance)",
        y = "PC2 (12.5% variance)"
      ) +

      theme_classic(
        base_size = 14
      ) +

      theme(
        plot.title = element_text(
          face = "bold",
          size = 17
        ),
        plot.subtitle = element_text(
          size = 12
        ),
        axis.title = element_text(
          face = "bold"
        ),
        legend.title = element_text(
          face = "bold"
        ),
        legend.position = "right"
      )

    ggplotly(
      p,
      tooltip = "text"
    ) |>
      layout(
        autosize = TRUE,
        margin = list(
          l = 80,
          r = 40,
          b = 80,
          t = 90
        )
      ) |>
      config(
        responsive = TRUE,
        displaylogo = FALSE
      )
  })

  # --------------------------------------------------------
  # Clinical variable selector
  # --------------------------------------------------------

  selected_variable <- reactive({

    if (input$clinical_variable == "Age") {

      list(
        column = "Age",
        summary_name = "age",
        y_label = "Age (years)"
      )

    } else {

      list(
        column = "SOFA",
        summary_name = "sofa",
        y_label = "SOFA score"
      )
    }
  })

  # --------------------------------------------------------
  # Clinical plot
  # --------------------------------------------------------

  output$clinical_plot <- renderPlotly({

    selected <- selected_variable()

    plot_data <- pca[
      !is.na(pca[[selected$column]]),
      ,
      drop = FALSE
    ]

    p <- ggplot(
      plot_data,
      aes(
        x = Cluster,
        y = .data[[selected$column]],
        fill = Cluster
      )
    ) +
      geom_boxplot(
        alpha = 0.72,
        outlier.shape = NA
      ) +
      geom_jitter(
        width = 0.14,
        alpha = 0.50,
        size = 1.6
      ) +
      labs(
        x = "Transcriptomic cluster",
        y = selected$y_label
      ) +
      theme_minimal(
        base_size = 13
      ) +
      theme(
        legend.position = "none"
      )

    ggplotly(
      p,
      tooltip = c("x", "y")
    ) |>
      layout(
        autosize = TRUE,
        margin = list(
          l = 75,
          r = 25,
          b = 70,
          t = 20
        )
      ) |>
      config(
        responsive = TRUE,
        displaylogo = FALSE
      )
  })

  # --------------------------------------------------------
  # Kruskal-Wallis result
  # --------------------------------------------------------

  output$kw_result <- renderUI({

    selected <- selected_variable()

    result <- kw[
      kw$Variable == selected$column,
      ,
      drop = FALSE
    ]

    if (nrow(result) == 0) {
      return(
        p("Statistical result not found.")
      )
    }

    p_value <- result$PValue[1]
    effect_size <- result$EpsilonSquared[1]

    p_display <- if (p_value < 0.001) {
      format(
        p_value,
        scientific = TRUE,
        digits = 3
      )
    } else {
      formatC(
        p_value,
        format = "f",
        digits = 4
      )
    }

    tagList(

      p(
        strong("Test: "),
        "Kruskal-Wallis"
      ),

      p(
        strong("\u03C7\u00B2("),
        result$DegreesOfFreedom[1],
        strong(") = "),
        round(result$Statistic[1], 3)
      ),

      p(
        strong("p = "),
        p_display
      ),

      p(
        strong("Epsilon-squared = "),
        round(effect_size, 4)
      ),

      hr(),

      p(
        if (p_value < 0.05) {
          "The overall difference among clusters is statistically significant."
        } else {
          "The overall difference among clusters is not statistically significant."
        }
      )
    )
  })

  # --------------------------------------------------------
  # Clinical summary table
  # --------------------------------------------------------

  output$clinical_summary_table <- renderDT({

    selected <- selected_variable()

    table_data <- clinical[
      tolower(clinical$Variable) == selected$summary_name,
      ,
      drop = FALSE
    ]

    table_data <- table_data[
      ,
      c(
        "Cluster",
        "TotalSamples",
        "Nonmissing",
        "Missing",
        "Mean",
        "SD",
        "Median",
        "Q1",
        "Q3",
        "Minimum",
        "Maximum"
      )
    ]

    datatable(
      table_data,
      rownames = FALSE,
      options = list(
        pageLength = 4,
        searching = FALSE,
        paging = FALSE,
        info = FALSE,
        autoWidth = TRUE,
        scrollX = TRUE
      ),
      colnames = c(
        "Cluster",
        "Total N",
        "Available",
        "Missing",
        "Mean",
        "SD",
        "Median",
        "Q1",
        "Q3",
        "Minimum",
        "Maximum"
      )
    )
  })

  # --------------------------------------------------------
  # Dunn post-hoc table
  # --------------------------------------------------------

  output$dunn_table <- renderDT({

    if (input$clinical_variable == "Age") {
      table_data <- age_dunn
    } else {
      table_data <- sofa_dunn
    }

    table_data <- table_data[
      ,
      c(
        "Comparison",
        "Z",
        "AdjustedPValue",
        "Significance"
      )
    ]

    table_data$Z <- round(
      table_data$Z,
      3
    )

    table_data$AdjustedPValue <- signif(
      table_data$AdjustedPValue,
      4
    )

    datatable(
      table_data,
      rownames = FALSE,
      options = list(
        pageLength = 6,
        searching = FALSE,
        paging = FALSE,
        info = FALSE,
        autoWidth = TRUE,
        scrollX = TRUE
      ),
      colnames = c(
        "Comparison",
        "Z",
        "Adjusted p-value",
        "Significance"
      )
    )
  })

  # ========================================================
  # DIFFERENTIAL EXPRESSION
  # ========================================================

  de_data <- reactive({

    req(input$de_comparison)

    data <- load_de_table(
      input$de_comparison
    )

    data$DisplayGene <- ifelse(
      is.na(data$SYMBOL) |
        trimws(data$SYMBOL) == "",
      data$Ensembl,
      data$SYMBOL
    )

    data$VolcanoStatus <- "Not significant"

    significant <- !is.na(data$padj) &
      data$padj < 0.05 &
      !is.na(data$log2FoldChange) &
      abs(data$log2FoldChange) >= 1

    data$VolcanoStatus[
      significant &
        data$log2FoldChange > 0
    ] <- "Upregulated"

    data$VolcanoStatus[
      significant &
        data$log2FoldChange < 0
    ] <- "Downregulated"

    data$MinusLog10Padj <- -log10(
      pmax(
        data$padj,
        .Machine$double.xmin,
        na.rm = FALSE
      )
    )

    data
  })

  output$de_comparison_title <- renderUI({

    data <- de_data()

    comparison_name <- unique(
      data$Comparison
    )

    h4(
      comparison_name[1]
    )
  })

  output$de_stat_cards <- renderUI({

    data <- de_data()

    sig <- data[
      data$VolcanoStatus != "Not significant",
      ,
      drop = FALSE
    ]

    up_count <- sum(
      sig$VolcanoStatus == "Upregulated"
    )

    down_count <- sum(
      sig$VolcanoStatus == "Downregulated"
    )

    div(
      class = "stat-grid-three",

      stat_card(
        "Significant DEGs",
        format(
          nrow(sig),
          big.mark = ","
        ),
        "padj < 0.05 and |log2FC| >= 1"
      ),

      stat_card(
        "Upregulated",
        format(
          up_count,
          big.mark = ","
        ),
        "Higher in the first cluster"
      ),

      stat_card(
        "Downregulated",
        format(
          down_count,
          big.mark = ","
        ),
        "Lower in the first cluster"
      )
    )
  })

  output$volcano_plot <- renderPlotly({

    data <- de_data()

    volcano_colors <- c(
      "Downregulated" = "#2C7FB8",
      "Not significant" = "#B8B8B8",
      "Upregulated" = "#D95F02"
    )

    p <- ggplot(
      data,
      aes(
        x = log2FoldChange,
        y = MinusLog10Padj,
        color = VolcanoStatus,
        text = paste0(
          "Gene: ", DisplayGene,
          "<br>Ensembl: ", Ensembl,
          "<br>Gene name: ",
          ifelse(
            is.na(GENENAME),
            "Not available",
            GENENAME
          ),
          "<br>log2FC: ",
          round(
            log2FoldChange,
            3
          ),
          "<br>Adjusted p: ",
          format_p(padj),
          "<br>Base mean: ",
          round(
            baseMean,
            2
          ),
          "<br>Direction: ",
          VolcanoStatus
        )
      )
    ) +

      geom_point(
        alpha = 0.65,
        size = 1.7
      ) +

      geom_vline(
        xintercept = c(
          -1,
          1
        ),
        linetype = "dashed",
        linewidth = 0.6
      ) +

      geom_hline(
        yintercept = -log10(0.05),
        linetype = "dashed",
        linewidth = 0.6
      ) +

      scale_color_manual(
        values = volcano_colors,
        name = "Direction"
      ) +

      labs(
        x = "Shrunken log2 fold change",
        y = "-log10 adjusted p-value"
      ) +

      theme_classic(
        base_size = 14
      ) +

      theme(
        axis.title = element_text(
          face = "bold"
        ),
        legend.position = "right"
      )

    ggplotly(
      p,
      tooltip = "text"
    ) |>
      layout(
        autosize = TRUE,
        margin = list(
          l = 80,
          r = 40,
          b = 80,
          t = 25
        )
      ) |>
      config(
        responsive = TRUE,
        displaylogo = FALSE
      )
  })

  output$de_table <- renderDT({

    data <- de_data()

    table_data <- switch(
      input$de_table_filter,

      "significant" = data[
        data$VolcanoStatus !=
          "Not significant",
        ,
        drop = FALSE
      ],

      "up" = data[
        data$VolcanoStatus ==
          "Upregulated",
        ,
        drop = FALSE
      ],

      "down" = data[
        data$VolcanoStatus ==
          "Downregulated",
        ,
        drop = FALSE
      ],

      data
    )

    table_data <- table_data[
      ,
      c(
        "DisplayGene",
        "Ensembl",
        "GENENAME",
        "ENTREZID",
        "baseMean",
        "log2FoldChange",
        "padj",
        "VolcanoStatus"
      )
    ]

    names(table_data) <- c(
      "Gene",
      "Ensembl",
      "Gene name",
      "Entrez ID",
      "Base mean",
      "log2FC",
      "Adjusted p-value",
      "Direction"
    )

    table_data$`Base mean` <- round(
      table_data$`Base mean`,
      2
    )

    table_data$log2FC <- round(
      table_data$log2FC,
      3
    )

    table_data$`Adjusted p-value` <- signif(
      table_data$`Adjusted p-value`,
      4
    )

    datatable(
      table_data,
      rownames = FALSE,
      filter = "top",
      extensions = "Buttons",
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        autoWidth = TRUE,
        dom = "Bfrtip",
        buttons = c(
          "copy",
          "csv"
        )
      )
    )
  })

  # ========================================================
  # GO ENRICHMENT
  # ========================================================

  go_data <- reactive({

    req(
      input$go_comparison,
      input$go_direction
    )

    data <- load_go_table(
      input$go_comparison,
      input$go_direction
    )

    if (nrow(data) == 0) {
      return(data)
    }

    data$GeneRatioNumeric <- parse_ratio(
      data$GeneRatio
    )

    data$DescriptionWrapped <- wrap_go_term(
      data$Description,
      width = 42
    )

    data
  })

  output$go_selection_title <- renderUI({

    label <- names(de_comparisons)[
      de_comparisons == input$go_comparison
    ]

    h4(
      paste0(
        label,
        " | ",
        input$go_direction,
        " genes"
      )
    )
  })

  output$go_stat_cards <- renderUI({

    data <- go_data()

    term_count <- nrow(data)

    gene_set_size <- if (term_count > 0) {
      ratio_denominator(data$GeneRatio)
    } else {
      NA_integer_
    }

    background_size <- if (term_count > 0) {
      ratio_denominator(data$BgRatio)
    } else {
      NA_integer_
    }

    div(
      class = "stat-grid-three",

      stat_card(
        "Significant GO terms",
        format(term_count, big.mark = ","),
        "GO Biological Process"
      ),

      stat_card(
        "GO gene-set size",
        ifelse(
          is.na(gene_set_size),
          "N/A",
          format(gene_set_size, big.mark = ",")
        ),
        "Annotated genes used by enrichGO"
      ),

      stat_card(
        "GO background size",
        ifelse(
          is.na(background_size),
          "N/A",
          format(background_size, big.mark = ",")
        ),
        "Annotated background genes"
      )
    )
  })

  output$go_plot_message <- renderUI({

    data <- go_data()

    if (nrow(data) == 0) {
      div(
        class = "alert alert-secondary",
        paste0(
          "No significantly enriched GO Biological Process terms ",
          "were identified for ",
          input$go_direction,
          " genes in this comparison."
        )
      )
    }
  })

  output$go_plot <- renderPlotly({

    data <- go_data()

    validate(
      need(
        nrow(data) > 0,
        ""
      )
    )

    top_n <- min(
      as.integer(input$go_top_n),
      nrow(data)
    )

    plot_data <- data[
      order(
        data$p.adjust,
        data$qvalue,
        -data$Count,
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]

    plot_data <- head(
      plot_data,
      top_n
    )

    plot_data$DescriptionWrapped <- factor(
      plot_data$DescriptionWrapped,
      levels = rev(plot_data$DescriptionWrapped)
    )

    p <- ggplot(
      plot_data,
      aes(
        x = GeneRatioNumeric,
        y = DescriptionWrapped,
        size = Count,
        color = p.adjust,
        text = paste0(
          "GO term: ", Description,
          "<br>GO ID: ", ID,
          "<br>GeneRatio: ", GeneRatio,
          "<br>Background ratio: ", BgRatio,
          "<br>Gene count: ", Count,
          "<br>Adjusted p: ", format_p(p.adjust),
          "<br>q-value: ", format_p(qvalue)
        )
      )
    ) +

      geom_point(alpha = 0.9) +

      scale_size_continuous(
        name = "Gene count",
        range = c(5, 14)
      ) +

      scale_color_viridis_c(
        name = "Adjusted\np-value",
        direction = -1,
        trans = "reverse"
      ) +

      labs(
        x = "Gene ratio",
        y = NULL
      ) +

      theme_classic(
        base_size = 13
      ) +

      theme(
        axis.title.x = element_text(
          face = "bold"
        ),
        axis.text.y = element_text(
          size = 10
        ),
        legend.title = element_text(
          face = "bold"
        ),
        legend.position = "right"
      )

    ggplotly(
      p,
      tooltip = "text"
    ) |>
      layout(
        autosize = TRUE,
        margin = list(
          l = 320,
          r = 50,
          b = 80,
          t = 20
        )
      ) |>
      config(
        responsive = TRUE,
        displaylogo = FALSE
      )
  })

  output$go_table <- renderDT({

    data <- go_data()

    if (nrow(data) == 0) {

      empty_table <- data.frame(
        Message = paste0(
          "No significant GO Biological Process terms for ",
          input$go_direction,
          " genes in this comparison."
        )
      )

      return(
        datatable(
          empty_table,
          rownames = FALSE,
          options = list(
            searching = FALSE,
            paging = FALSE,
            info = FALSE
          )
        )
      )
    }

    table_data <- data[
      ,
      c(
        "ID",
        "Description",
        "GeneRatio",
        "BgRatio",
        "Count",
        "p.adjust",
        "qvalue"
      )
    ]

    names(table_data) <- c(
      "GO ID",
      "Biological Process",
      "Gene Ratio",
      "Background Ratio",
      "Gene Count",
      "Adjusted p-value",
      "q-value"
    )

    table_data$`Adjusted p-value` <- signif(
      table_data$`Adjusted p-value`,
      4
    )

    table_data$`q-value` <- signif(
      table_data$`q-value`,
      4
    )

    datatable(
      table_data,
      rownames = FALSE,
      filter = "top",
      extensions = "Buttons",
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        autoWidth = TRUE,
        dom = "Bfrtip",
        buttons = c(
          "copy",
          "csv"
        )
      )
    )
  })

}

# ==========================================================
# RUN APPLICATION
# ==========================================================

shinyApp(
  ui = ui,
  server = server
)
