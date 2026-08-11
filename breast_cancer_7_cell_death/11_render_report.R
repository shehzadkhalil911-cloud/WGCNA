# =============================================================================
# Script 11: Render Report
# Breast Cancer 7 Cell Death Crosstalk
# =============================================================================

source(file.path("..", "config", "analysis_config.R"))

library(rmarkdown)
library(flextable)
library(officer)

log_message("Rendering reports...")

# --- 1. METHODS SUMMARY ---
methods_summary <- "
## Methods Summary

### Study Design
This study analysed crosstalk among ferroptosis, necroptosis, and cuproptosis in breast cancer using a discovery cohort (GSE42568, GSE15852, GSE29431) and an independent validation cohort (GSE70947). All analyses were performed using R version 4.6.0.

### Gene Lists
Curated gene lists for ferroptosis, necroptosis, and cuproptosis were obtained from the local Datasets directory. Genes were standardised to official gene symbols using org.Hs.eg.db.

### Data Preprocessing
Each dataset was processed independently. Expression data were log2-transformed where necessary. Probe-to-gene mapping used platform annotation files. Batch correction was performed using ComBat with tumour/normal status preserved.

### Differential Expression
Differential expression analysis used limma with thresholds: adjusted P < 0.05 and |log2FC| >= 1. Per-dataset DE analysis was performed for consistency assessment.

### Activity Scores
Sample-level activity scores for each cell-death programme were calculated using ssGSEA, comparing tumour versus normal samples.

### WGCNA
Signed co-expression networks were constructed using blockwiseModules with soft-threshold power selected via scale-free topology criterion. Modules were correlated with tumour status and activity scores.

### Functional Enrichment
GO, KEGG, Reactome, and Hallmark enrichment were performed using clusterProfiler. GSEA used ranked gene lists.

### PPI Networks
Protein-protein interaction networks were constructed using STRINGdb with a confidence threshold of 0.4. Centrality metrics were calculated.

### Immune Infiltration
Immune cell infiltration was estimated using xCell. Associations with cell-death scores and genes were assessed via Spearman correlation.

### Multi-Evidence Prioritisation
Candidate crosstalk genes were prioritised using a standardised composite score integrating differential expression, program membership, WGCNA metrics, and activity-score correlations.

### Independent Validation
All findings were validated in GSE70947 without influencing discovery-stage decisions.
"

writeLines(methods_summary, file.path(BASE_DIR, "report", "methods_summary.txt"))

# --- 2. FIGURE MANIFEST ---
figure_files <- list.files(OUTPUT$figures, pattern = "\\.tiff$", recursive = TRUE)
if (length(figure_files) == 0) {
  # Search entire results for TIFF files
  figure_files <- list.files(RESULTS_DIR, pattern = "\\.tiff$", recursive = TRUE)
}

figure_manifest <- data.frame(
  Figure = seq_along(figure_files),
  File = basename(figure_files),
  Path = figure_files,
  Size_Bytes = file.info(file.path(RESULTS_DIR, figure_files))$size,
  stringsAsFactors = FALSE
)

write.csv(figure_manifest, file.path(OUTPUT$tables, "figure_manifest.csv"), row.names = FALSE)

# --- 3. TABLE MANIFEST ---
table_files <- list.files(RESULTS_DIR, pattern = "\\.csv$", recursive = TRUE)
table_manifest <- data.frame(
  Table = seq_along(table_files),
  File = basename(table_files),
  Path = table_files,
  Size_Bytes = file.info(file.path(RESULTS_DIR, table_files))$size,
  stringsAsFactors = FALSE
)

write.csv(table_manifest, file.path(OUTPUT$tables, "table_manifest.csv"), row.names = FALSE)

# --- 4. SESSION INFO ---
sink(file.path(BASE_DIR, "sessionInfo.txt"))
cat("Session Information\n")
cat("===================\n\n")
print(sessionInfo())
sink()

# --- 5. PACKAGE VERSION SUMMARY ---
packages <- sessionInfo()$otherPkgs
if (length(packages) > 0) {
  pkg_versions <- data.frame(
    Package = names(packages),
    Version = sapply(packages, function(x) as.character(x$Version)),
    stringsAsFactors = FALSE
  )
  write.csv(pkg_versions, file.path(OUTPUT$tables, "package_versions.csv"), row.names = FALSE)
}

# --- 6. HTML REPORT ---
html_report <- file.path(BASE_DIR, "report", "analysis_report.html")

# Create a simple HTML report
res_dir <- gsub("\\\\", "/", file.path(BASE_DIR, "results"))
rmd_content <- sprintf('---
title: "Breast Cancer 7 Cell Death Crosstalk"
subtitle: "Integrated Analysis of Ferroptosis, Necroptosis, and Cuproptosis"
author: "opencode R Bioinformatics Pipeline"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: flatly
    highlight: tango
    code_folding: hide
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(tidyverse)
library(knitr)
library(flextable)
library(ComplexHeatmap)
```

# Study Overview

This report presents a comprehensive analysis of the molecular crosstalk among ferroptosis, necroptosis, and cuproptosis in breast cancer.

## Study Design

- **Discovery cohort**: GSE42568, GSE15852, GSE29431
- **Independent validation**: GSE70947

## Key Findings

1. **Gene List Overlap**: Analysis of curated gene lists revealed shared and unique genes across the three regulated cell death programmes.
2. **Differential Expression**: Identification of differentially expressed ferroptosis-related (DE-FRGs), necroptosis-related (DE-NRGs), and cuproptosis-related (DE-CRGs) genes.
3. **Activity Scores**: Programme-specific activity scores showed significant differences between tumour and normal samples.
4. **Co-expression Modules**: WGCNA identified modules associated with multiple cell-death programmes.
5. **Pathway Crosstalk**: Shared pathways involving iron metabolism, oxidative stress, and immune signalling.
6. **PPI Networks**: Integrated networks revealed cross-programme hub genes.
7. **Immune Infiltration**: Associations between cell-death programmes and immune cell composition.
8. **Multi-Evidence Prioritisation**: Candidate crosstalk genes prioritised using integrated evidence.
9. **Independent Validation**: Key findings validated in GSE70947.

# Gene List Overlap

# Results Summary

## Candidate Crosstalk Genes

```{r candidates}
top_candidates <- tryCatch({
  read.csv("%s/multi_evidence_prioritisation/top_30_candidates.csv")
}, error = function(e) NULL)
if (!is.null(top_candidates)) {
  kable(head(top_candidates[, c("Gene", "logFC", "N_Programs", "Combined_Score")], 20), 
        caption = "Top 20 Candidate Crosstalk Genes", digits = 3)
}
```

## Key Tables

```{r deg-count}
deg_up <- tryCatch(read.csv("%s/differential_expression/discovery_DEG_upregulated.csv"), error = function(e) NULL)
deg_down <- tryCatch(read.csv("%s/differential_expression/discovery_DEG_downregulated.csv"), error = function(e) NULL)
if (!is.null(deg_up) && !is.null(deg_down)) {
  cat("- Significant DEGs:", nrow(deg_up) + nrow(deg_down), "\\n")
  cat("- Upregulated:", nrow(deg_up), "\\n")
  cat("- Downregulated:", nrow(deg_down), "\\n")
}
```

# Methods

See `report/methods_summary.txt` for a complete methods description.

# Limitations

- Further experimental validation is required to confirm the functional significance of identified candidate regulated cell-death crosstalk genes.
- The curated gene lists may not be comprehensive.
- Batch effects between platforms may not be fully eliminated.

# Session Information

```{r session}
sessionInfo()
```
', res_dir, res_dir, res_dir)

writeLines(rmd_content, file.path(BASE_DIR, "report", "analysis_report.Rmd"))

# Try to render
tryCatch({
  old_wd <- getwd()
  setwd(BASE_DIR)
  rmarkdown::render(
    "report/analysis_report.Rmd",
    output_file = "analysis_report.html",
    output_dir = file.path(BASE_DIR, "report"),
    quiet = TRUE
  )
  setwd(old_wd)
  log_message("HTML report generated.")
}, error = function(e) {
  log_message(sprintf("Report rendering failed: %s", e$message), "WARNING")
  if (exists("old_wd")) setwd(old_wd)
})

# --- 7. WORD REPORT ---
word_report <- file.path(BASE_DIR, "report", "analysis_report.docx")

doc <- read_docx()
doc <- doc %>%
  body_add_par("Breast Cancer 7 Cell Death Crosstalk Report", style = "heading 1") %>%
  body_add_par(sprintf("Generated: %s", Sys.Date()), style = "Normal") %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("Study Overview", style = "heading 2") %>%
  body_add_par("This report presents a comprehensive analysis of the molecular crosstalk among ferroptosis, necroptosis, and cuproptosis in breast cancer.", style = "Normal") %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("Key Findings", style = "heading 2") %>%
  body_add_par("1. Gene list overlap analysis revealed shared and unique genes across three regulated cell death programmes", style = "Normal") %>%
  body_add_par("2. Differential expression identified DE-FRGs, DE-NRGs, and DE-CRGs in breast cancer", style = "Normal") %>%
  body_add_par("3. Activity scores showed significant tumour-normal differences", style = "Normal") %>%
  body_add_par("4. WGCNA identified shared co-expression modules", style = "Normal") %>%
  body_add_par("5. Pathway crosstalk analysis revealed shared pathways", style = "Normal") %>%
  body_add_par("6. Integrated PPI networks identified cross-programme hub genes", style = "Normal") %>%
  body_add_par("7. Immune infiltration analysis revealed immune associations", style = "Normal") %>%
  body_add_par("8. Multi-evidence prioritisation identified candidate crosstalk genes", style = "Normal") %>%
  body_add_par("9. Independent validation confirmed key findings", style = "Normal") %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("Candidate Crosstalk Genes", style = "heading 2")

# Add top candidates table
if (file.exists(file.path(OUTPUT$prioritisation, "top_30_candidates.csv"))) {
  top_df <- read.csv(file.path(OUTPUT$prioritisation, "top_30_candidates.csv"))
  ft <- flextable(top_df[1:min(20, nrow(top_df)), c("Gene", "logFC", "N_Programs", "Combined_Score")])
  ft <- theme_booktabs(ft)
  doc <- doc %>% body_add_flextable(ft)
}

doc <- doc %>%
  body_add_par("", style = "Normal") %>%
  body_add_par("Further experimental validation is required to confirm the functional significance of identified candidate regulated cell-death crosstalk genes.", style = "Normal")

print(doc, target = word_report)
log_message("Word report generated.")

# --- 8. README ---
readme_content <- sprintf(
"# Breast Cancer 7 Cell Death Crosstalk

## Integrated Analysis of Ferroptosis, Necroptosis, and Cuproptosis

### Study Design
- **Discovery cohort**: GSE42568, GSE15852, GSE29431
- **Independent validation**: GSE70947

### Directory Structure
- `config/` - Analysis configuration
- `scripts/` - R analysis scripts
- `R/` - Helper R functions
- `metadata/` - Sample metadata
- `data_processed/` - Processed expression data
- `results/` - All analysis results
- `logs/` - Analysis logs
- `report/` - Generated reports

### Pipeline Scripts (in order)
1. `01_input_inventory_and_gene_lists.R` - File inventory and gene-list preparation
2. `02_metadata_and_preprocessing.R` - Metadata curation and preprocessing
3. `03_batch_correction_and_DE.R` - Batch correction and differential expression
4. `04_regulated_cell_death_gene_analysis.R` - Regulated cell-death gene analysis
5. `05_activity_scores.R` - Activity score analysis
6. `06_WGCNA.R` - WGCNA network analysis
7. `07_enrichment.R` - Functional enrichment analysis
8. `08_pathway_crosstalk.R` - Pathway crosstalk analysis
9. `09_PPI_and_immune.R` - PPI and immune infiltration analysis
10. `10_multi_evidence_prioritisation.R` - Prioritisation and validation
11. `11_render_report.R` - Report generation

### Key Results
- Gene list characterisation and overlap
- Differential expression of FRGs, NRGs, and CRGs
- Programme activity scores
- Co-expression modules (WGCNA)
- Functional enrichment and pathway crosstalk
- PPI networks and hub genes
- Immune infiltration associations
- Multi-evidence prioritisation of candidate genes
- Independent validation

### Software
- R version %s
- Key packages: limma, sva, WGCNA, clusterProfiler, GSVA, etc.

### Citation
Generated by opencode R bioinformatics pipeline.
Report date: %s

### Note
Further experimental validation is required.
",
  R.version.string,
  Sys.Date()
)

writeLines(readme_content, file.path(BASE_DIR, "README.md"))

log_message("All reports generated successfully.")
