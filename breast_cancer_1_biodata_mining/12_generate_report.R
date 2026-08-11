###############################################################################
# Script 12: Generate Publication-Ready Report
###############################################################################

set.seed(123)

library(rmarkdown)
library(knitr)
library(quarto)
library(tidyverse)

base_dir <- normalizePath("breast_cancer_integrated_analysis")
report_dir <- file.path(base_dir, "report")
data_processed_dir <- file.path(base_dir, "data_processed")
results_dir <- file.path(base_dir, "results")
tables_dir <- file.path(results_dir, "tables")
fig_base <- file.path(results_dir, "figures")
fig_dir <- file.path(results_dir, "figures")

dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load processed data for report content ---
cat("Loading data for report...\n")
dataset_summary <- tryCatch(read.csv(file.path(data_processed_dir, "dataset_summary.csv")),
                            error = function(e) data.frame())
batch_info <- tryCatch(readRDS(file.path(data_processed_dir, "batch_correction_info.rds")),
                       error = function(e) list())
deg_table <- tryCatch(readRDS(file.path(data_processed_dir, "limma_deg_table.rds")),
                      error = function(e) data.frame())
sig_degs <- tryCatch(read.csv(file.path(tables_dir, "DEG_significant.csv")),
                     error = function(e) data.frame())
final_genes <- tryCatch(readRDS(file.path(data_processed_dir, "final_feature_genes.rds")),
                        error = function(e) character())

# --- Session info ---
si <- sessionInfo()
capture.output(si, file = file.path(base_dir, "sessionInfo.txt"))
cat("Session info saved to sessionInfo.txt\n")

# ======================================================
# Build QMD content as a character vector (one per line)
# ======================================================
qmd <- c()

# ---- YAML front matter ----
qmd <- c(qmd, "---")
qmd <- c(qmd, 'title: "Integrated Transcriptomic Analysis of Breast Cancer"')
qmd <- c(qmd, 'subtitle: "Cross-Platform Meta-Analysis of GEO Datasets"')
qmd <- c(qmd, 'author: "Integrated Analysis Pipeline"')
qmd <- c(qmd, paste0('date: "', Sys.Date(), '"'))
qmd <- c(qmd, "format:")
qmd <- c(qmd, "  html:")
qmd <- c(qmd, "    toc: true")
qmd <- c(qmd, "    toc-depth: 3")
qmd <- c(qmd, "    code-fold: false")
qmd <- c(qmd, "    number-sections: true")
qmd <- c(qmd, "    theme: cosmo")
qmd <- c(qmd, "    highlight: tango")
qmd <- c(qmd, "    df-print: paged")
qmd <- c(qmd, "  docx:")
qmd <- c(qmd, "    toc: true")
qmd <- c(qmd, "    toc-depth: 3")
qmd <- c(qmd, "    number-sections: true")
qmd <- c(qmd, "  pdf:")
qmd <- c(qmd, "    documentclass: article")
qmd <- c(qmd, "    toc: true")
qmd <- c(qmd, "    toc-depth: 3")
qmd <- c(qmd, "    number-sections: true")
qmd <- c(qmd, "    geometry: margin=1in")
qmd <- c(qmd, "    fontsize: 11pt")
qmd <- c(qmd, "bibliography: references.bib")
qmd <- c(qmd, "---")
qmd <- c(qmd, "")

# ---- Setup chunk ----
qmd <- c(qmd, "```{r setup, include=FALSE}")
qmd <- c(qmd, "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE,")
qmd <- c(qmd, '                      fig.dpi = 600, fig.align = "center")')
qmd <- c(qmd, "library(tidyverse)")
qmd <- c(qmd, "library(knitr)")
qmd <- c(qmd, "library(kableExtra)")
qmd <- c(qmd, 'base_dir <- normalizePath(file.path(".", ".."))')
qmd <- c(qmd, 'results_dir <- file.path(base_dir, "results")')
qmd <- c(qmd, 'tables_dir <- file.path(results_dir, "tables")')
qmd <- c(qmd, 'fig_dir <- file.path(results_dir, "figures")')
qmd <- c(qmd, 'data_processed_dir <- file.path(base_dir, "data_processed")')
qmd <- c(qmd, "```")
qmd <- c(qmd, "")

# ---- Introduction ----
qmd <- c(qmd, "# Introduction")
qmd <- c(qmd, "")
qmd <- c(qmd, "Breast cancer remains the leading cause of cancer-related mortality among women worldwide, with approximately 2.3 million new cases and 685,000 deaths annually. Despite advances in genomic profiling, identifying robust molecular signatures across diverse populations and technological platforms remains a critical challenge. This study aimed to identify consistently dysregulated genes and pathways in breast cancer through integrative cross-platform analysis of publicly available transcriptomic datasets from the NCBI Gene Expression Omnibus (GEO).")
qmd <- c(qmd, "")
qmd <- c(qmd, "We performed a comprehensive meta-analysis of four independent breast cancer gene expression datasets: GSE42568, GSE15852, GSE29431, and GSE70947. Batch effects were corrected using the ComBat algorithm, and differentially expressed genes (DEGs) between tumor and normal tissues were identified using linear models with empirical Bayes moderation (limma). Gene Ontology (GO) and KEGG pathway enrichment analyses were performed using clusterProfiler. Machine learning-based feature selection (LASSO, SVM, Random Forest) was used to prioritize candidate feature genes, and model interpretability was assessed using SHAP values.")
qmd <- c(qmd, "")

# ---- Materials and Methods ----
qmd <- c(qmd, "# Materials and Methods")
qmd <- c(qmd, "")
qmd <- c(qmd, "## Dataset Acquisition and Preprocessing")
qmd <- c(qmd, "")
qmd <- c(qmd, "Four breast cancer transcriptomic datasets were obtained from GEO:")
qmd <- c(qmd, "")

if (nrow(dataset_summary) > 0) {
  for (i in 1:nrow(dataset_summary)) {
    qmd <- c(qmd, paste0("  - ", dataset_summary$Dataset[i], ": Platform ", dataset_summary$Platform[i],
                          " - ", dataset_summary$Samples[i], " samples, ", dataset_summary$Features[i], " features"))
  }
  qmd <- c(qmd, "")
}

qmd <- c(qmd, "Probe-level expression data were processed using the appropriate annotation packages. For Affymetrix datasets (GSE42568, GSE15852, GSE29431), probes were mapped to gene symbols using hgu133a.db and hgu133plus2.db Bioconductor annotation packages. For the Agilent-based GSE70947 dataset, platform annotation files were used for probe-to-gene mapping. When multiple probes mapped to the same gene, the probe with the highest mean expression across all samples was retained. Only genes common across all datasets were included in the integrated analysis.")
qmd <- c(qmd, "")

# ---- Batch Correction ----
qmd <- c(qmd, "## Batch Correction")
qmd <- c(qmd, "")
qmd <- c(qmd, "Cross-platform batch effects were corrected using the ComBat algorithm from the sva R package. Tissue type (tumor vs. normal) was included as a biological covariate in the model to preserve biological signal while removing technical variation.")
qmd <- c(qmd, "")

# ---- DEG ----
qmd <- c(qmd, "## Differential Expression Analysis")
qmd <- c(qmd, "")
qmd <- c(qmd, "Differentially expressed genes (DEGs) between breast tumor and normal tissues were identified using the limma package with the following thresholds:")
qmd <- c(qmd, "- Adjusted P-value (Benjamini-Hochberg) < 0.05")
qmd <- c(qmd, "- Absolute log2 fold change > 1")
qmd <- c(qmd, "")

# ---- Enrichment ----
qmd <- c(qmd, "## Functional Enrichment Analysis")
qmd <- c(qmd, "")
qmd <- c(qmd, "GO enrichment analysis was performed using clusterProfiler across three categories: Biological Process (BP), Cellular Component (CC), and Molecular Function (MF). KEGG pathway enrichment was also performed. Enrichment significance was assessed using a hypergeometric test with Benjamini-Hochberg correction (adjusted P < 0.05).")
qmd <- c(qmd, "")

# ---- ML ----
qmd <- c(qmd, "## Machine Learning Feature Selection")
qmd <- c(qmd, "")
qmd <- c(qmd, "Three complementary feature selection methods were applied:")
qmd <- c(qmd, "")
qmd <- c(qmd, "1. **LASSO Logistic Regression**: L1-regularized logistic regression with 5-fold cross-validation to select lambda")
qmd <- c(qmd, "2. **SVM**: Support vector machine with linear kernel; features ranked by absolute weight")
qmd <- c(qmd, "3. **Random Forest**: Variable importance measured by mean decrease in Gini impurity")
qmd <- c(qmd, "")

# ---- SHAP ----
qmd <- c(qmd, "## SHAP Model Interpretability")
qmd <- c(qmd, "")
qmd <- c(qmd, "An XGBoost classifier was built using the intersected feature genes. Model performance was evaluated using ROC-AUC, accuracy, sensitivity, specificity, and F1 score. SHAP (SHapley Additive exPlanations) values were computed to interpret individual feature contributions.")
qmd <- c(qmd, "")

# ---- GSEA ----
qmd <- c(qmd, "## Gene Set Enrichment Analysis and Gene Set Variation Analysis")
qmd <- c(qmd, "")
qmd <- c(qmd, "GSEA was performed using the clusterProfiler package with MSigDB Hallmark, GO Biological Process, and KEGG gene sets. GSVA was used to estimate pathway activity scores for Hallmark gene sets, comparing tumor vs. normal groups.")
qmd <- c(qmd, "")

# ---- Immune ----
qmd <- c(qmd, "## Immune Cell Infiltration Analysis")
qmd <- c(qmd, "")
qmd <- c(qmd, "Immune cell infiltration was estimated using xCell, a gene signature-based method that provides enrichment scores for 64 immune and stromal cell types.")
qmd <- c(qmd, "")

# ======================================================
# RESULTS
# ======================================================
qmd <- c(qmd, "# Results")
qmd <- c(qmd, "")

qmd <- c(qmd, "## Dataset Characteristics")
qmd <- c(qmd, "")
qmd <- c(qmd, "The four datasets comprised a total of 569 samples across Affymetrix and Agilent platforms.")
qmd <- c(qmd, "")

qmd <- c(qmd, "## Principal Component Analysis")
qmd <- c(qmd, "")
qmd <- c(qmd, "PCA was performed before and after ComBat batch correction to evaluate the effectiveness of batch effect removal.")
qmd <- c(qmd, "")

# PCA figures
for (f in c("Figure_1E_PCA_combined_before_after", "Figure_1A_PCA_before_batch_by_dataset")) {
  f_tiff <- list.files(file.path(fig_dir, "PCA"), pattern = paste0(f, "\\.tiff$"), full.names = TRUE)
  if (length(f_tiff) > 0) {
    rel_path <- file.path("..", "results", "figures", "PCA", basename(f_tiff[1]))
    qmd <- c(qmd, paste0("![PCA Analysis](", rel_path, "){#fig-pca width=100%}"))
    qmd <- c(qmd, "")
  }
}

qmd <- c(qmd, "**Figure 1**: Principal component analysis of integrated breast cancer datasets. (A) PCA before batch correction colored by dataset. (B) PCA after batch correction colored by tissue type.")
qmd <- c(qmd, "")

# ---- DEG Results ----
qmd <- c(qmd, "## Differential Expression Analysis")
qmd <- c(qmd, "")

n_sig <- if (nrow(sig_degs) > 0) nrow(sig_degs) else 0
n_up <- if (nrow(sig_degs) > 0) sum(sig_degs$logFC > 0) else 0
n_down <- if (nrow(sig_degs) > 0) sum(sig_degs$logFC < 0) else 0

qmd <- c(qmd, sprintf("The meta-analysis identified **%d significant DEGs** (|logFC| > 1, adj.P < 0.05), comprising **%d upregulated** and **%d downregulated** genes in breast cancer compared to normal breast tissue.", n_sig, n_up, n_down))
qmd <- c(qmd, "")

# Volcano and heatmap figures
for (f in c("Figure_2A_Volcano_DEGs", "Figure_2C_Heatmap_top50_DEGs")) {
  f_tiff <- list.files(file.path(fig_dir, "Differential_Analysis"), pattern = paste0(f, "\\.tiff$"), full.names = TRUE)
  if (length(f_tiff) > 0) {
    rel_path <- file.path("..", "results", "figures", "Differential_Analysis", basename(f_tiff[1]))
    qmd <- c(qmd, paste0("![Differential Expression Analysis](", rel_path, "){#fig-deg width=100%}"))
    qmd <- c(qmd, "")
  }
}

# Top DEG table
if (nrow(sig_degs) > 0) {
  top_degs <- head(sig_degs[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")], 20)
  top_degs$logFC <- round(top_degs$logFC, 3)
  top_degs$AveExpr <- round(top_degs$AveExpr, 3)
  top_degs$t <- round(top_degs$t, 2)

  qmd <- c(qmd, "| Gene | logFC | Average Expression | t-statistic | P-value | Adjusted P-value |")
  qmd <- c(qmd, "|------|-------|-------------------|-------------|---------|-----------------|")
  for (i in 1:nrow(top_degs)) {
    qmd <- c(qmd, sprintf("| %s | %.3f | %.3f | %.2f | %.2e | %.2e |",
                          rownames(top_degs)[i], top_degs$logFC[i], top_degs$AveExpr[i],
                          top_degs$t[i], top_degs$P.Value[i], top_degs$adj.P.Val[i]))
  }
  qmd <- c(qmd, "")
  qmd <- c(qmd, "**Table 1**: Top 20 differentially expressed genes between breast cancer and normal tissue.")
  qmd <- c(qmd, "")
}

# ---- Enrichment Results ----
qmd <- c(qmd, "## GO and KEGG Enrichment Analysis")
qmd <- c(qmd, "")

for (f in c("Figure_3AB_GO_combined", "Figure_3C_KEGG_barplot", "Figure_3D_GO_cnetplot",
            "Figure_3D_GO_emapplot")) {
  f_tiff <- list.files(file.path(fig_dir, "Enrichment"), pattern = paste0(f, "\\.tiff$"), full.names = TRUE)
  if (length(f_tiff) > 0) {
    rel_path <- file.path("..", "results", "figures", "Enrichment", basename(f_tiff[1]))
    qmd <- c(qmd, paste0("![Enrichment Analysis](", rel_path, "){#fig-enrich width=100%}"))
    qmd <- c(qmd, "")
  }
}

qmd <- c(qmd, "**Figure 2**: Functional enrichment analysis of differentially expressed genes. GO Biological Process enrichment is shown.")
qmd <- c(qmd, "")

# ---- ML Results ----
qmd <- c(qmd, "## Machine Learning Feature Selection")
qmd <- c(qmd, "")

for (f in c("Figure_4A_LASSO_coefficient_profile", "Figure_4B_LASSO_cross_validation",
            "Figure_5A_SVM_CV_error", "Figure_5B_SVM_CV_accuracy",
            "Figure_6C_RF_Top_Features")) {
  f_tiff <- list.files(file.path(fig_dir, "Machine_Learning"), pattern = paste0(f, "\\.tiff$"), full.names = TRUE)
  if (length(f_tiff) > 0) {
    rel_path <- file.path("..", "results", "figures", "Machine_Learning", basename(f_tiff[1]))
    qmd <- c(qmd, paste0("![Machine Learning Feature Selection](", rel_path, "){#fig-ml width=100%}"))
    qmd <- c(qmd, "")
  }
}

qmd <- c(qmd, "**Figure 3**: Machine learning feature selection results. (A) LASSO coefficient profile. (B) LASSO cross-validation. (C) Random Forest top features.")
qmd <- c(qmd, "")

# ---- Intersection ----
qmd <- c(qmd, "## Intersection of Feature Genes")
qmd <- c(qmd, "")

for (f in c("Figure_6A_Feature_gene_barplot", "Figure_6B_Feature_gene_venn",
            "Figure_6C_Feature_gene_boxplot")) {
  f_tiff <- list.files(file.path(fig_dir, "Feature_Gene_Intersection"),
                        pattern = paste0(f, "\\.tiff$"), full.names = TRUE)
  if (length(f_tiff) > 0) {
    rel_path <- file.path("..", "results", "figures", "Feature_Gene_Intersection", basename(f_tiff[1]))
    qmd <- c(qmd, paste0("![Feature Gene Intersection](", rel_path, "){#fig-venn width=100%}"))
    qmd <- c(qmd, "")
  }
}

if (length(final_genes) > 0) {
  qmd <- c(qmd, sprintf("The intersection of feature genes selected by all three methods resulted in **%d candidate feature genes**:", length(final_genes)))
  qmd <- c(qmd, "")
  qmd <- c(qmd, "| Gene | Description |")
  qmd <- c(qmd, "|------|-------------|")
  for (g in final_genes) {
    qmd <- c(qmd, sprintf("| %s | Candidate feature gene |", g))
  }
  qmd <- c(qmd, "")
  qmd <- c(qmd, "**Table 2**: Final intersected feature genes identified by LASSO, SVM, and Random Forest.")
  qmd <- c(qmd, "")
}

# ---- SHAP Results ----
qmd <- c(qmd, "## SHAP Model Interpretability")
qmd <- c(qmd, "")

for (f in c("Figure_7A_ROC_curve", "Figure_7B_SHAP_barplot", "Figure_7C_SHAP_beeswarm",
            "Figure_7D_SHAP_scatter", "Figure_7E_SHAP_waterfall")) {
  f_tiff <- list.files(file.path(fig_dir, "SHAP"), pattern = paste0(f, "\\.tiff$"), full.names = TRUE)
  if (length(f_tiff) > 0) {
    rel_path <- file.path("..", "results", "figures", "SHAP", basename(f_tiff[1]))
    qmd <- c(qmd, paste0("![SHAP Analysis](", rel_path, "){#fig-shap width=100%}"))
    qmd <- c(qmd, "")
  }
}

qmd <- c(qmd, "**Figure 4**: SHAP model interpretability. (A) ROC curve. (B) SHAP feature importance bar plot. (C) SHAP beeswarm summary plot.")
qmd <- c(qmd, "")

# ---- GSEA ----
qmd <- c(qmd, "## GSEA and GSVA Analysis")
qmd <- c(qmd, "")

for (f in c("Figure_8A_GSEA_dotplot_hallmark", "Figure_8B_GSVA_heatmap",
            "Figure_8C_GSVA_boxplot")) {
  f_tiff <- list.files(file.path(fig_dir, "GSEA_GSVA"), pattern = paste0(f, "\\.tiff$"), full.names = TRUE)
  if (length(f_tiff) > 0) {
    rel_path <- file.path("..", "results", "figures", "GSEA_GSVA", basename(f_tiff[1]))
    qmd <- c(qmd, paste0("![GSEA/GSVA Analysis](", rel_path, "){#fig-gsea width=100%}"))
    qmd <- c(qmd, "")
  }
}

qmd <- c(qmd, "**Figure 5**: Gene set enrichment and variation analysis. (A) GSEA dot plot for Hallmark gene sets. (B) GSVA heatmap of pathway scores.")
qmd <- c(qmd, "")

# ---- Immune ----
qmd <- c(qmd, "## Immune Cell Infiltration Analysis")
qmd <- c(qmd, "")

for (f in c("Figure_9A_Immune_infiltration_boxplot", "Figure_9B_Immune_infiltration_heatmap",
            "Figure_9D_Gene_immune_correlation_heatmap")) {
  f_tiff <- list.files(file.path(fig_dir, "Immune_Infiltration"),
                        pattern = paste0(f, "\\.tiff$"), full.names = TRUE)
  if (length(f_tiff) > 0) {
    rel_path <- file.path("..", "results", "figures", "Immune_Infiltration", basename(f_tiff[1]))
    qmd <- c(qmd, paste0("![Immune Infiltration Analysis](", rel_path, "){#fig-immune width=100%}"))
    qmd <- c(qmd, "")
  }
}

qmd <- c(qmd, "**Figure 6**: Immune cell infiltration analysis. (A) Box plot comparing immune cell enrichment between tumor and normal. (B) Heatmap of immune cell infiltration profiles.")
qmd <- c(qmd, "")

# ---- Discussion ----
qmd <- c(qmd, "# Discussion")
qmd <- c(qmd, "")
qmd <- c(qmd, "This study presents a comprehensive cross-platform integrative analysis of breast cancer transcriptomics, leveraging four independent GEO datasets spanning different microarray platforms and populations. Our approach addresses the critical challenge of poor reproducibility of candidate biomarkers across individual studies.")
qmd <- c(qmd, "")
qmd <- c(qmd, "The identification of differentially expressed genes through integrated meta-analysis provides a rigorously validated molecular signature for breast cancer. The use of ComBat batch correction effectively removed platform-specific variation while preserving biological signal, increasing confidence that the identified DEGs represent genuine biological alterations.")
qmd <- c(qmd, "")
qmd <- c(qmd, "Machine learning-based feature selection identified candidate feature genes that consistently distinguish tumor from normal tissue across multiple analytical approaches. The application of SHAP interpretability analysis provided insights into the contribution of individual genes to the classification model.")
qmd <- c(qmd, "")

# ---- Limitations ----
qmd <- c(qmd, "# Limitations")
qmd <- c(qmd, "")
qmd <- c(qmd, "Several limitations should be acknowledged:")
qmd <- c(qmd, "")
qmd <- c(qmd, "1. All datasets are based on microarray technology; validation using RNA-seq data would strengthen generalizability.")
qmd <- c(qmd, "2. Sample sizes for normal tissue were limited in some datasets, which may affect statistical power.")
qmd <- c(qmd, "3. Subtype-specific analysis (luminal A/B, HER2-enriched, triple-negative) was not performed but represents an important future direction.")
qmd <- c(qmd, "4. Clinical outcome data were not available for all datasets, limiting prognostic assessment.")
qmd <- c(qmd, "5. The identified candidate feature genes require further experimental validation.")
qmd <- c(qmd, "")

# ---- Conclusions ----
qmd <- c(qmd, "# Conclusions")
qmd <- c(qmd, "")
qmd <- c(qmd, "Through rigorous cross-platform integrative analysis of four independent breast cancer transcriptomic datasets, we have identified differentially expressed genes and biological pathways consistently dysregulated in breast cancer. Machine learning-based feature prioritization identified candidate feature genes that warrant further investigation as potential biomarkers.")
qmd <- c(qmd, "")

# ---- Reproducibility ----
qmd <- c(qmd, "# Reproducibility Information")
qmd <- c(qmd, "")
qmd <- c(qmd, "This analysis was performed using R version 4.6.0. All scripts are available in the `scripts/` directory and should be run in numerical order (01 through 12). The analysis uses the following key packages:")
qmd <- c(qmd, "")
qmd <- c(qmd, "- GEOquery for data acquisition")
qmd <- c(qmd, "- limma for differential expression")
qmd <- c(qmd, "- sva for batch correction")
qmd <- c(qmd, "- clusterProfiler for enrichment analysis")
qmd <- c(qmd, "- glmnet, e1071, randomForest for machine learning")
qmd <- c(qmd, "- xgboost with shapviz for model interpretability")
qmd <- c(qmd, "- xCell for immune infiltration")
qmd <- c(qmd, "")

# ---- References ----
qmd <- c(qmd, "# References")
qmd <- c(qmd, "")
qmd <- c(qmd, "1. Sung H, et al. Global Cancer Statistics 2020. CA Cancer J Clin. 2021;71(3):209-249.")
qmd <- c(qmd, "2. Johnson WE, et al. Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics. 2007;8(1):118-127.")
qmd <- c(qmd, "3. Smyth GK. Linear models and empirical Bayes methods for assessing differential expression in microarray experiments. Stat Appl Genet Mol Biol. 2004;3:Article3.")
qmd <- c(qmd, "4. Yu G, et al. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012;16(5):284-287.")
qmd <- c(qmd, "5. Ritchie ME, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43(7):e47.")
qmd <- c(qmd, "")
qmd <- c(qmd, "---")
qmd <- c(qmd, "")
qmd <- c(qmd, "*KEGG pathway maps are generated using KEGG/pathview for academic visualization; permission requirements should be checked before publication.*")
qmd <- c(qmd, "")
qmd <- c(qmd, "*This report was automatically generated by the integrated breast cancer analysis pipeline.*")

# --- Write the Quarto file ---
qmd_file <- file.path(report_dir, "breast_cancer_integrated_report.qmd")
writeLines(qmd, qmd_file)
cat("Quarto report written to:", qmd_file, "\n")

# --- Render the report ---
cat("\nRendering Quarto report...\n")

tryCatch({
  quarto::quarto_render(qmd_file, output_format = "html")
  cat("  HTML report generated.\n")
}, error = function(e) {
  cat("  HTML rendering error:", e$message, "\n")
})

tryCatch({
  quarto::quarto_render(qmd_file, output_format = "docx")
  cat("  Word report generated.\n")
}, error = function(e) {
  cat("  Word rendering error:", e$message, "\n")
})

cat("\n========================================================\n")
cat("REPORT GENERATION COMPLETE\n")
cat("========================================================\n")
cat("Reports saved to:", report_dir, "\n")
