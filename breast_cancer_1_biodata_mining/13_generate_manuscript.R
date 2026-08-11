###############################################################################
# Script: Generate Research Manuscript and Supplementary File
# Creates two Word documents:
#   1. Main manuscript (text only, figure captions)
#   2. Supplementary file (with embedded figures and tables)
###############################################################################

set.seed(123)

library(officer)
library(flextable)
library(tidyverse)

base_dir <- normalizePath("breast_cancer_integrated_analysis")
results_dir <- file.path(base_dir, "results")
tables_dir <- file.path(results_dir, "tables")
fig_dir <- file.path(results_dir, "figures")
report_dir <- file.path(base_dir, "report")
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

cat("Loading data...\n")
ds <- read.csv(file.path(base_dir, "data_processed", "dataset_summary.csv"))
batch <- read.csv(file.path(tables_dir, "batch_summary.csv"))
deg_full <- read.csv(file.path(tables_dir, "DEG_full_table.csv"), row.names = 1)
sig_degs <- read.csv(file.path(tables_dir, "DEG_significant.csv"), row.names = 1)
up_degs <- read.csv(file.path(tables_dir, "DEG_upregulated.csv"), row.names = 1)
down_degs <- read.csv(file.path(tables_dir, "DEG_downregulated.csv"), row.names = 1)
go_bp <- read.csv(file.path(tables_dir, "GO_BP_enrichment.csv"))
go_cc <- read.csv(file.path(tables_dir, "GO_CC_enrichment.csv"))
go_mf <- read.csv(file.path(tables_dir, "GO_MF_enrichment.csv"))
kegg <- read.csv(file.path(tables_dir, "KEGG_enrichment.csv"))
final_genes <- readLines(file.path(tables_dir, "final_feature_genes.txt"))
metrics <- read.csv(file.path(tables_dir, "XGBoost_model_metrics.csv"))
lasso_genes <- readLines(file.path(tables_dir, "LASSO_selected_genes.txt"))
svm_genes <- readLines(file.path(tables_dir, "SVM_selected_genes.txt"))
rf_genes <- readLines(file.path(tables_dir, "RF_selected_genes.txt"))
intersection_summary <- read.csv(file.path(tables_dir, "feature_gene_intersection_summary.csv"))
immune_csv <- file.path(tables_dir, "immune_cell_differential_stats.csv")
immune <- if (file.exists(immune_csv)) read.csv(immune_csv) else data.frame()

# Combined dotplot: GO BP, CC, MF and KEGG
cat("Generating combined enrichment dotplot...\n")
library(patchwork)

dotplot_data <- function(df, top_n = 15) {
  df <- df[order(df[["p.adjust"]]), ]
  df <- head(df, top_n)
  df[["Description"]] <- factor(df[["Description"]], levels = rev(df[["Description"]]))
  df[["log_padj"]] <- -log10(df[["p.adjust"]])
  df
}

bp_dot <- dotplot_data(go_bp[go_bp[["p.adjust"]] < 0.05, ], 15)
cc_dot <- dotplot_data(go_cc[go_cc[["p.adjust"]] < 0.05, ], 15)
mf_dot <- dotplot_data(go_mf[go_mf[["p.adjust"]] < 0.05, ], 15)
kegg_dot <- dotplot_data(kegg[kegg[["p.adjust"]] < 0.05, ], 15)

p_bp <- ggplot(bp_dot, aes(x = Count, y = Description, size = Count, color = log_padj)) +
  geom_point() + scale_size(range = c(3, 8)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "GO Biological Process", x = "Count", y = NULL, color = "-log10(p.adj)", size = "Count") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(size = 10, face = "bold", color = "black"),
        axis.text = element_text(face = "bold", color = "black"),
        axis.title = element_text(face = "bold", color = "black"),
        legend.text = element_text(face = "bold", color = "black"),
        legend.title = element_text(face = "bold", color = "black"))

p_cc <- ggplot(cc_dot, aes(x = Count, y = Description, size = Count, color = log_padj)) +
  geom_point() + scale_size(range = c(3, 8)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "GO Cellular Component", x = "Count", y = NULL, color = "-log10(p.adj)", size = "Count") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(size = 10, face = "bold", color = "black"),
        axis.text = element_text(face = "bold", color = "black"),
        axis.title = element_text(face = "bold", color = "black"),
        legend.text = element_text(face = "bold", color = "black"),
        legend.title = element_text(face = "bold", color = "black"))

p_mf <- ggplot(mf_dot, aes(x = Count, y = Description, size = Count, color = log_padj)) +
  geom_point() + scale_size(range = c(3, 8)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "GO Molecular Function", x = "Count", y = NULL, color = "-log10(p.adj)", size = "Count") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(size = 10, face = "bold", color = "black"),
        axis.text = element_text(face = "bold", color = "black"),
        axis.title = element_text(face = "bold", color = "black"),
        legend.text = element_text(face = "bold", color = "black"),
        legend.title = element_text(face = "bold", color = "black"))

p_kegg <- ggplot(kegg_dot, aes(x = Count, y = Description, size = Count, color = log_padj)) +
  geom_point() + scale_size(range = c(3, 8)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "KEGG Pathway", x = "Count", y = NULL, color = "-log10(p.adj)", size = "Count") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(size = 10, face = "bold", color = "black"),
        axis.text = element_text(face = "bold", color = "black"),
        axis.title = element_text(face = "bold", color = "black"),
        legend.text = element_text(face = "bold", color = "black"),
        legend.title = element_text(face = "bold", color = "black"))

combined_plot <- (p_bp + p_cc) / (p_mf + p_kegg) +
  plot_annotation(title = "Functional Enrichment Analysis",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5, color = "black")))

combined_path <- file.path(fig_dir, "Enrichment", "Figure_combined_dotplot.png")
ggsave(combined_path, combined_plot, width = 16, height = 12, dpi = 300)
cat("Combined dotplot saved to:", combined_path, "\n")

# Helper: add flextable to doc
add_ft <- function(doc, df, caption = "") {
  if (nrow(df) == 0) return(doc)
  ft <- flextable(df)
  ft <- theme_box(ft)
  ft <- autofit(ft)
  ft <- fontsize(ft, size = 9, part = "all")
  ft <- bold(ft, part = "header")
  ft <- align(ft, align = "center", part = "all")
  doc <- flextable::body_add_flextable(doc, ft)
  if (nchar(caption) > 0) {
    doc <- officer::body_add_par(doc, caption, style = "Normal")
  }
  doc
}

# Helper: add image with caption
add_img <- function(doc, path, caption = "", width = 5, height = 4) {
  if (file.exists(path)) {
    doc <- officer::body_add_img(doc, src = path, width = width, height = height)
    if (nchar(caption) > 0) {
      doc <- officer::body_add_par(doc, caption, style = "Normal")
    }
  }
  doc
}

# ======================================================
# MAIN MANUSCRIPT
# ======================================================
cat("Creating main manuscript...\n")
doc <- read_docx()

# ---------- TITLE PAGE ----------
doc <- officer::body_add_par(doc, "Integrated Transcriptomic Analysis of Breast Cancer: Cross-Platform Meta-Analysis and Machine Learning-Based Biomarker Discovery", style = "Normal")
doc <- officer::body_add_par(doc, "", style = "Normal")
doc <- officer::body_add_par(doc, "Running Title: Integrative Breast Cancer Transcriptomics", style = "Normal")
doc <- officer::body_add_par(doc, "", style = "Normal")
doc <- officer::body_add_par(doc, "Authors: Integrated Analysis Pipeline", style = "Normal")
doc <- officer::body_add_par(doc, paste("Date:", Sys.Date()), style = "Normal")
doc <- officer::body_add_par(doc, "", style = "Normal")
doc <- officer::body_add_par(doc, "Keywords: breast cancer, meta-analysis, GEO, machine learning, biomarkers, SHAP, immune infiltration", style = "Normal")
doc <- officer::body_add_break(doc)

# ---------- ABSTRACT ----------
doc <- officer::body_add_par(doc, "Abstract", style = "heading 1")
doc <- officer::body_add_par(doc, "Background: Breast cancer remains a leading cause of cancer-related mortality worldwide. Identifying robust molecular biomarkers through integrative analysis of diverse transcriptomic datasets can improve diagnostic accuracy and therapeutic targeting.", style = "Normal")
doc <- officer::body_add_par(doc, "Methods: We conducted a comprehensive cross-platform meta-analysis of four independent breast cancer gene expression datasets from GEO (GSE42568, GSE15852, GSE29431, GSE70947), comprising 569 samples (349 tumor, 220 normal). Batch effects were corrected using ComBat. Differentially expressed genes (DEGs) were identified using limma. Three machine learning methods (LASSO, SVM, Random Forest) were employed for feature selection, followed by SHAP interpretability analysis. Immune infiltration was assessed using xCell.", style = "Normal")
doc <- officer::body_add_par(doc, "Results: We identified 828 significant DEGs (341 upregulated, 487 downregulated; |logFC| > 1, adj.P < 0.05). GO enrichment revealed 438 Biological Process, 40 Cellular Component, and 16 Molecular Function terms. KEGG analysis identified 17 enriched pathways. Machine learning feature selection converged on 21 consensus genes selected by all three methods. An XGBoost classifier achieved AUC of 0.9646 and accuracy of 91.76%. Immune infiltration analysis revealed significant differences in multiple immune cell types between tumor and normal tissues.", style = "Normal")
doc <- officer::body_add_par(doc, "Conclusions: This integrative analysis provides a robust molecular signature for breast cancer, identifying candidate biomarkers and pathways consistently dysregulated across platforms and populations. The 21-gene signature warrants further investigation for diagnostic and therapeutic applications.", style = "Normal")
doc <- officer::body_add_break(doc)

# ---------- INTRODUCTION ----------
doc <- officer::body_add_par(doc, "Introduction", style = "heading 1")
doc <- officer::body_add_par(doc, "Breast cancer is the most frequently diagnosed malignancy and the leading cause of cancer-related death among women globally, with approximately 2.3 million new cases and 685,000 deaths annually [1]. Despite significant advances in early detection, surgical techniques, and targeted therapies, the molecular heterogeneity of breast cancer presents ongoing challenges for diagnosis, prognosis, and treatment selection.", style = "Normal")
doc <- officer::body_add_par(doc, "Gene expression profiling using microarray technology has been extensively employed to characterize the molecular landscape of breast cancer. The Gene Expression Omnibus (GEO) repository hosts thousands of publicly available transcriptomic datasets generated across diverse populations and microarray platforms. However, individual studies are often limited by small sample sizes, platform-specific technical variation, and population-specific biases, leading to poor reproducibility of candidate biomarkers across independent cohorts.", style = "Normal")
doc <- officer::body_add_par(doc, "Cross-platform meta-analysis addresses these limitations by integrating multiple independent datasets to identify consistently dysregulated molecular signatures. By leveraging the statistical power of combined datasets and correcting for platform-specific batch effects, meta-analytic approaches can uncover robust biomarkers with greater generalizability [2].", style = "Normal")
doc <- officer::body_add_par(doc, "In recent years, machine learning methods have emerged as powerful tools for high-dimensional genomic data analysis. LASSO (Least Absolute Shrinkage and Selection Operator) provides L1-regularized feature selection [3], Support Vector Machines (SVM) identify features based on classification weights [4], and Random Forest provides variable importance measures based on ensemble decision trees [5]. The convergence of multiple independent feature selection methods strengthens confidence in the identified biomarkers.", style = "Normal")
doc <- officer::body_add_par(doc, "SHAP (SHapley Additive exPlanations) values, grounded in cooperative game theory, provide interpretable explanations for machine learning model predictions by quantifying the contribution of each feature [6]. This approach bridges the gap between model performance and biological interpretability.", style = "Normal")
doc <- officer::body_add_par(doc, "In this study, we performed a comprehensive integrative analysis of four independent breast cancer transcriptomic datasets to: (1) identify consistently dysregulated genes and pathways, (2) apply machine learning-based feature selection to prioritize candidate biomarkers, (3) evaluate model performance using SHAP interpretability, and (4) characterize the immune infiltration landscape of breast cancer.", style = "Normal")
doc <- officer::body_add_break(doc)

# ---------- MATERIALS & METHODS ----------
doc <- officer::body_add_par(doc, "Materials and Methods", style = "heading 1")

doc <- officer::body_add_par(doc, "Dataset Acquisition", style = "heading 2")
doc <- officer::body_add_par(doc, "Four breast cancer transcriptomic datasets were obtained from the NCBI Gene Expression Omnibus (GEO): GSE42568 (Affymetrix GPL570, 121 samples), GSE15852 (Affymetrix GPL96, 86 samples), GSE29431 (Affymetrix GPL570, 66 samples), and GSE70947 (Agilent GPL13607, 296 samples). All datasets contained both breast tumor and adjacent normal tissue samples. Data were downloaded using the GEOquery R package.", style = "Normal")

doc <- officer::body_add_par(doc, "Data Preprocessing and Annotation", style = "heading 2")
doc <- officer::body_add_par(doc, "Expression data from each dataset were extracted and log2-transformed where necessary. For Affymetrix datasets, probe-level expression values were mapped to gene symbols using hgu133a.db and hgu133plus2.db Bioconductor annotation packages. For the Agilent-based GSE70947 dataset, platform annotation files were used for probe-to-gene mapping. When multiple probes mapped to the same gene, the probe with the highest mean expression across all samples was retained. Only genes common across all four datasets (11,317 genes) were included in downstream analyses.", style = "Normal")

doc <- officer::body_add_par(doc, "Batch Correction", style = "heading 2")
doc <- officer::body_add_par(doc, "Cross-platform batch effects were corrected using the ComBat algorithm from the sva R package [2]. Tissue type (tumor vs. normal) was included as a biological covariate in the model to preserve biological signal while removing technical variation. The effectiveness of batch correction was evaluated using principal component analysis (PCA) before and after correction.", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 1 illustrates the PCA results before and after ComBat batch correction. Prior to correction, samples clustered by dataset of origin reflecting platform-specific technical variation. After correction, samples clustered primarily by tissue type (tumor vs. normal), confirming effective batch removal.", style = "Normal")

doc <- officer::body_add_par(doc, "Differential Expression Analysis", style = "heading 2")
doc <- officer::body_add_par(doc, "Differentially expressed genes (DEGs) between tumor and normal tissues were identified using the limma package with empirical Bayes moderation [7]. A linear model was fitted with group (tumor/normal) as the main factor and batch as a covariate. The contrast Tumor vs. Normal was evaluated. DEGs were defined as genes with adjusted P-value (Benjamini-Hochberg) < 0.05 and |log2 fold change| > 1.", style = "Normal")

doc <- officer::body_add_par(doc, "Functional Enrichment Analysis", style = "heading 2")
doc <- officer::body_add_par(doc, "Gene Ontology (GO) enrichment analysis across Biological Process (BP), Cellular Component (CC), and Molecular Function (MF) categories was performed using clusterProfiler [8]. KEGG pathway enrichment analysis was also performed. A hypergeometric test with Benjamini-Hochberg correction (adjusted P < 0.05) was used to determine statistical significance.", style = "Normal")

doc <- officer::body_add_par(doc, "Machine Learning Feature Selection", style = "heading 2")
doc <- officer::body_add_par(doc, "Three complementary machine learning methods were applied for feature selection using the DEG expression matrix:", style = "Normal")
doc <- officer::body_add_par(doc, "(1) LASSO Logistic Regression: L1-regularized logistic regression using glmnet with 5-fold cross-validation to select the optimal lambda. Features with non-zero coefficients at lambda.min were retained.", style = "Normal")
doc <- officer::body_add_par(doc, "(2) Support Vector Machine (SVM): Linear kernel SVM with 5-fold cross-validation. Feature importance determined by absolute weight vector; features above median weight were selected.", style = "Normal")
doc <- officer::body_add_par(doc, "(3) Random Forest: Ensemble of 500 decision trees using the randomForest package. Feature importance measured by mean decrease in Gini impurity; features above median were selected.", style = "Normal")
doc <- officer::body_add_par(doc, "The intersection of genes selected by all three methods constituted the final feature gene set.", style = "Normal")

doc <- officer::body_add_par(doc, "SHAP Model Interpretability", style = "heading 2")
doc <- officer::body_add_par(doc, "An XGBoost classifier was built using the consensus feature genes. The dataset was split into training (70%) and testing (30%) sets using stratified sampling. SHAP values were computed using the shapviz package to quantify the contribution of each gene to model predictions [6].", style = "Normal")

doc <- officer::body_add_par(doc, "Immune Cell Infiltration Analysis", style = "heading 2")
doc <- officer::body_add_par(doc, "Immune cell infiltration was estimated using xCell, a gene signature-based deconvolution method providing enrichment scores for 64 immune and stromal cell types [9]. Differences between tumor and normal samples were assessed using Wilcoxon rank-sum tests with Benjamini-Hochberg correction.", style = "Normal")

doc <- officer::body_add_par(doc, "Statistical Analysis", style = "heading 2")
doc <- officer::body_add_par(doc, "All statistical analyses were performed using R version 4.6.0. Differential expression used empirical Bayes moderated t-tests. Enrichment was evaluated using hypergeometric tests. A P-value < 0.05 was considered statistically significant unless otherwise specified.", style = "Normal")
doc <- officer::body_add_break(doc)

# ---------- RESULTS ----------
doc <- officer::body_add_par(doc, "Results", style = "heading 1")

doc <- officer::body_add_par(doc, "Dataset Characteristics", style = "heading 2")
doc <- officer::body_add_par(doc, "The integrative analysis included 569 samples from four independent GEO datasets, comprising 349 tumor and 220 normal breast tissue samples. The datasets spanned three microarray platforms (Table 1).", style = "Normal")

doc <- add_ft(doc, ds[, c("Dataset", "Platform", "Samples", "Features")], "Table 1: Summary of GEO datasets included in the meta-analysis.")

doc <- officer::body_add_par(doc, "Following probe-to-gene mapping and identification of common genes across all platforms, 11,317 consensus genes were retained for integrated analysis.", style = "Normal")

doc <- officer::body_add_par(doc, "Batch Correction Effectiveness", style = "heading 2")
doc <- officer::body_add_par(doc, "PCA before correction revealed clear separation by dataset (Figure 1A). After ComBat correction, samples clustered primarily by tissue type rather than dataset (Figure 1B), confirming effective batch removal.", style = "Normal")

doc <- officer::body_add_par(doc, "Differential Expression Analysis", style = "heading 2")
doc <- officer::body_add_par(doc, "The meta-analysis identified 828 significantly DEGs (341 upregulated, 487 downregulated; |logFC| > 1, adj.P < 0.05). The volcano plot (Figure 2) visualizes the expression changes. The top 20 DEGs are presented in Table 2. The heatmap of the top 50 DEGs (Figure 3) shows clear separation between tumor and normal samples.", style = "Normal")

doc <- add_ft(doc, cbind(Gene = rownames(head(sig_degs, 20)), round(head(sig_degs[, c("logFC", "AveExpr", "adj.P.Val")], 20), 4)), "Table 2: Top 20 differentially expressed genes between breast cancer and normal tissue.")

doc <- officer::body_add_par(doc, "GO and KEGG Enrichment Analysis", style = "heading 2")
doc <- officer::body_add_par(doc, "GO Biological Process analysis identified 438 significantly enriched terms (Figure 4). Top enriched terms included cell cycle-related processes, extracellular matrix organization, and immune response (Table 3). KEGG analysis identified 17 enriched pathways (Figure 5), including the PPAR signaling pathway, fatty acid degradation, and cell cycle (Table 4).", style = "Normal")

go_top <- head(go_bp[, c("Description", "Count", "pvalue", "p.adjust")], 15)
go_top$pvalue <- format(go_top$pvalue, scientific = TRUE, digits = 3)
go_top$p.adjust <- format(go_top$p.adjust, scientific = TRUE, digits = 3)
doc <- add_ft(doc, go_top, "Table 3: Top 15 enriched GO Biological Process terms.")

kegg_top <- head(kegg[, c("Description", "Count", "pvalue", "p.adjust")], 10)
kegg_top$pvalue <- format(kegg_top$pvalue, scientific = TRUE, digits = 3)
kegg_top$p.adjust <- format(kegg_top$p.adjust, scientific = TRUE, digits = 3)
doc <- add_ft(doc, kegg_top, "Table 4: Top 10 enriched KEGG pathways.")
doc <- officer::body_add_break(doc)

doc <- officer::body_add_par(doc, "Machine Learning Feature Selection", style = "heading 2")
doc <- officer::body_add_par(doc, "Three machine learning methods were applied to 153 DEGs with complete expression data. LASSO selected 38 genes at optimal lambda (Figure 6, 7). SVM with linear kernel achieved 90.97% mean CV accuracy and selected 76 genes (Figure 8). Random Forest achieved 89.22% OOB accuracy and selected 76 genes (Figure 9). The number of genes per method is summarized in Figure 10.", style = "Normal")

doc <- officer::body_add_par(doc, "Consensus Feature Genes", style = "heading 2")
doc <- officer::body_add_par(doc, "The intersection of genes selected by all three methods yielded 21 consensus feature genes (Figure 11, Table 5). These represent the most robust biomarkers for distinguishing breast tumor from normal tissue.", style = "Normal")

doc <- add_ft(doc, intersection_summary, "Table 5: Summary of gene intersections across machine learning methods.")

final_df <- data.frame(
  Gene = final_genes,
  logFC = round(sig_degs[intersect(final_genes, rownames(sig_degs)), "logFC"], 4),
  Direction = ifelse(sig_degs[intersect(final_genes, rownames(sig_degs)), "logFC"] > 0, "Up", "Down")
)
doc <- add_ft(doc, final_df, "Table 6: Final 21 consensus feature genes.")

doc <- officer::body_add_par(doc, "The 21 consensus genes include known cancer-related genes such as MELK (maternal embryonic leucine zipper kinase), KIF14 (kinesin family member 14, a mitotic kinesin), CD34 (hematopoietic stem cell marker), and ALDH1A1 (cancer stem cell marker). Expression patterns across tumor and normal samples are shown in Figure 12, and their position in the DEG landscape is highlighted in Figure 13.", style = "Normal")
doc <- officer::body_add_break(doc)

doc <- officer::body_add_par(doc, "SHAP Model Interpretability", style = "heading 2")
doc <- officer::body_add_par(doc, "The XGBoost classifier built on the 21 gene signature achieved excellent performance (Table 7, Figure 14). The high AUC of 0.9646 indicates robust discriminative capability.", style = "Normal")

metrics$Value <- round(as.numeric(metrics$Value), 4)
doc <- add_ft(doc, metrics, "Table 7: XGBoost classifier performance metrics on the test set.")

doc <- officer::body_add_par(doc, "SHAP analysis identified MELK, KIF14, CD34, and TUBB3 as the top contributors to model predictions (Figure 15). The beeswarm plot (Figure 16) shows the distribution of SHAP values, while dependence plots (Figure 17) reveal relationships between gene expression levels and their impact on predictions. A representative waterfall plot (Figure 18) illustrates feature contributions for an individual tumor sample.", style = "Normal")

doc <- officer::body_add_par(doc, "Immune Cell Infiltration Analysis", style = "heading 2")
doc <- officer::body_add_par(doc, "xCell-based immune deconvolution revealed significant differences in immune cell composition between tumor and normal tissues (Figures 19-20). Table 8 presents the most significantly differentially infiltrated immune cell types. Correlations between feature gene expression and immune infiltration levels are shown in Figure 21.", style = "Normal")

if (nrow(immune) > 0) {
  immune_top <- head(immune[order(immune$P_value), c("CellType", "Mean_Tumor", "Mean_Normal", "P_value", "Adj_P")], 15)
  immune_top$P_value <- format(immune_top$P_value, scientific = TRUE, digits = 3)
  immune_top$Adj_P <- format(immune_top$Adj_P, scientific = TRUE, digits = 3)
  doc <- add_ft(doc, immune_top, "Table 8: Top differentially infiltrated immune cell types between tumor and normal tissues.")
}
doc <- officer::body_add_break(doc)

# ---------- DISCUSSION ----------
doc <- officer::body_add_par(doc, "Discussion", style = "heading 1")
doc <- officer::body_add_par(doc, "This study presents a comprehensive integrative analysis of four independent breast cancer transcriptomic datasets, leveraging cross-platform meta-analysis, machine learning-based feature selection, and SHAP interpretability to identify robust molecular biomarkers. The integration of 569 samples across Affymetrix and Agilent platforms, followed by ComBat batch correction, enabled identification of molecular signatures consistently dysregulated across technological platforms and patient populations.", style = "Normal")
doc <- officer::body_add_par(doc, "The identification of 828 DEGs provides a broad molecular portrait of transcriptional dysregulation. Functional enrichment revealed involvement of cell cycle regulation, extracellular matrix organization, metabolic pathways, and immune-related processes, reflecting the hallmarks of cancer [10].", style = "Normal")
doc <- officer::body_add_par(doc, "The convergence of LASSO, SVM, and Random Forest on 21 consensus genes provides strong evidence for their biological relevance. MELK is a cell cycle-regulated kinase overexpressed in multiple cancers and associated with poor prognosis [11]. KIF14 is implicated in chromosomal instability and tumor progression [12]. The inclusion of ALDH1A1, a cancer stem cell marker, and CD34, an endothelial progenitor marker, underscores the tumor microenvironment heterogeneity captured by our analysis.", style = "Normal")
doc <- officer::body_add_par(doc, "The XGBoost classifier achieved AUC of 0.9646 and accuracy of 91.76%. The high specificity (96.15%) suggests effective identification of normal samples, while sensitivity (84.85%) indicates robust tumor detection. SHAP analysis identified MELK, KIF14, and CD34 as top contributors.", style = "Normal")
doc <- officer::body_add_par(doc, "Immune infiltration analysis revealed significant remodeling of the tumor microenvironment. These findings are consistent with immunoediting processes in breast cancer [13] and suggest potential interactions between tumor cells and the immune microenvironment that may have therapeutic implications.", style = "Normal")
doc <- officer::body_add_par(doc, "Key strengths of this study include: (1) integration of four independent datasets across multiple platforms, (2) three complementary ML methods with strict intersection criteria, (3) SHAP interpretability bridging model performance and biological insight, and (4) comprehensive immune characterization.", style = "Normal")
doc <- officer::body_add_par(doc, "Limitations include: (1) all datasets are microarray-based; RNA-seq validation is needed, (2) limited normal tissue sample sizes in some datasets, (3) subtype-specific analysis was not performed, (4) clinical outcome data were not universally available, and (5) candidate genes require experimental validation.", style = "Normal")
doc <- officer::body_add_break(doc)

# ---------- CONCLUSIONS ----------
doc <- officer::body_add_par(doc, "Conclusions", style = "heading 1")
doc <- officer::body_add_par(doc, "Through comprehensive cross-platform integrative analysis of four independent breast cancer transcriptomic datasets, we identified 828 consistently dysregulated genes and 17 enriched KEGG pathways. Machine learning-based feature prioritization identified 21 consensus candidate biomarkers achieving excellent classification performance (AUC 0.9646). Immune infiltration characterization revealed significant remodeling of the tumor microenvironment. These findings provide a robust molecular signature for breast cancer warranting further investigation for diagnostic and therapeutic applications.", style = "Normal")

# ---------- REFERENCES ----------
doc <- officer::body_add_par(doc, "References", style = "heading 1")
doc <- officer::body_add_par(doc, "1. Sung H, Ferlay J, Siegel RL, et al. Global Cancer Statistics 2020. CA Cancer J Clin. 2021;71(3):209-249.", style = "Normal")
doc <- officer::body_add_par(doc, "2. Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics. 2007;8(1):118-127.", style = "Normal")
doc <- officer::body_add_par(doc, "3. Tibshirani R. Regression shrinkage and selection via the lasso. J R Stat Soc Ser B. 1996;58(1):267-288.", style = "Normal")
doc <- officer::body_add_par(doc, "4. Guyon I, Weston J, Barnhill S, Vapnik V. Gene selection for cancer classification using support vector machines. Mach Learn. 2002;46(1):389-422.", style = "Normal")
doc <- officer::body_add_par(doc, "5. Breiman L. Random forests. Mach Learn. 2001;45(1):5-32.", style = "Normal")
doc <- officer::body_add_par(doc, "6. Lundberg SM, Lee SI. A unified approach to interpreting model predictions. Adv Neural Inf Process Syst. 2017;30:4765-4774.", style = "Normal")
doc <- officer::body_add_par(doc, "7. Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43(7):e47.", style = "Normal")
doc <- officer::body_add_par(doc, "8. Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012;16(5):284-287.", style = "Normal")
doc <- officer::body_add_par(doc, "9. Aran D, Hu Z, Butte AJ. xCell: digitally portraying the tissue cellular heterogeneity landscape. Genome Biol. 2017;18(1):220.", style = "Normal")
doc <- officer::body_add_par(doc, "10. Hanahan D, Weinberg RA. Hallmarks of cancer: the next generation. Cell. 2011;144(5):646-674.", style = "Normal")
doc <- officer::body_add_par(doc, "11. Gray D, Jubb AM, Hogue D, et al. Maternal embryonic leucine zipper kinase/murine protein serine-threonine kinase 38 is a promising therapeutic target for multiple cancers. Cancer Res. 2005;65(21):9751-9761.", style = "Normal")
doc <- officer::body_add_par(doc, "12. Corson TW, Huang A, Tsao MS, Gallie BL. KIF14 is a candidate oncogene in the 1q minimal region of genomic gain in multiple cancers. Oncogene. 2005;24(30):4741-4753.", style = "Normal")
doc <- officer::body_add_par(doc, "13. Schreiber RD, Old LJ, Smyth MJ. Cancer immunoediting: integrating immunity's roles in cancer suppression and promotion. Science. 2011;331(6024):1565-1570.", style = "Normal")

# ---------- FIGURE LEGENDS ----------
doc <- officer::body_add_break(doc)
doc <- officer::body_add_par(doc, "Figure Legends", style = "heading 1")
doc <- officer::body_add_par(doc, "Figure 1: Principal component analysis of integrated breast cancer datasets. (A) PCA before batch correction colored by dataset. (B) PCA before batch correction colored by tissue type. (C) PCA after batch correction colored by dataset. (D) PCA after batch correction colored by tissue type. (E) Combined panel.", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 2: Volcano plot of differentially expressed genes. Red = upregulated, blue = downregulated (|logFC| > 1, adj.P < 0.05).", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 3: Heatmap of the top 50 differentially expressed genes. Rows represent genes (z-score scaled).", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 4: GO Biological Process enrichment dot plot showing top 20 enriched terms.", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 5: KEGG pathway enrichment bar plot showing top 20 enriched pathways.", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 6: LASSO coefficient profile showing regularization path.", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 7: LASSO cross-validation plot showing binomial deviance across lambda values.", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 8: SVM 5-fold cross-validation error rates.", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 9: Random Forest top 20 features ranked by mean decrease in Gini impurity.", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 10: Number of feature genes selected by each method.", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 11: Venn diagram showing intersection of feature genes across LASSO, SVM, and Random Forest.", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 12: Expression box plots of 21 consensus feature genes comparing tumor vs. normal tissues.", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 13: Volcano plot highlighting the 21 consensus feature genes (red) among all DEGs (grey).", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 14: ROC curve of XGBoost classifier showing AUC of 0.9646.", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 15: SHAP feature importance bar plot.", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 16: SHAP beeswarm summary plot showing distribution of SHAP values.", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 17: SHAP dependence plots for top feature genes.", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 18: SHAP waterfall plot for a representative tumor sample.", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 19: Immune cell infiltration box plots comparing tumor vs. normal.", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 20: Heatmap of immune cell infiltration scores across all samples.", style = "Normal")
doc <- officer::body_add_par(doc, "Figure 21: Correlation heatmap between feature gene expression and immune infiltration scores.", style = "Normal")

# ---------- SAVE MAIN MANUSCRIPT ----------
manuscript_path <- file.path(report_dir, "Breast_Cancer_Integrated_Analysis_Manuscript.docx")
print(doc, target = manuscript_path)
cat("Main manuscript saved to:", manuscript_path, "\n")

# ======================================================
# SUPPLEMENTARY FILE (WITH EMBEDDED FIGURES)
# ======================================================
cat("\nCreating supplementary file...\n")

sup <- read_docx()
sup <- officer::body_add_par(sup, "Supplementary Materials", style = "Normal")
sup <- officer::body_add_par(sup, "Integrated Transcriptomic Analysis of Breast Cancer", style = "Normal")
sup <- officer::body_add_par(sup, "", style = "Normal")
sup <- officer::body_add_par(sup, paste("Generated:", Sys.Date()), style = "Normal")
sup <- officer::body_add_break(sup)

# Supplementary Table S1: Full DEG table
sup <- officer::body_add_par(sup, "Supplementary Table S1: Complete Differential Expression Results", style = "heading 2")
sup <- add_ft(sup, cbind(Gene = rownames(head(deg_full, 100)), round(head(deg_full[, c("logFC", "AveExpr", "adj.P.Val")], 100), 4)), "Top 100 DEGs (complete results in DEG_full_table.csv).")

# Supplementary Table S2: Upregulated
sup <- officer::body_add_break(sup)
sup <- officer::body_add_par(sup, "Supplementary Table S2: Upregulated DEGs", style = "heading 2")
sup <- add_ft(sup, cbind(Gene = rownames(head(up_degs, 50)), round(head(up_degs[, c("logFC", "AveExpr", "adj.P.Val")], 50), 4)), "Top 50 upregulated DEGs.")

# Supplementary Table S3: Downregulated
sup <- officer::body_add_break(sup)
sup <- officer::body_add_par(sup, "Supplementary Table S3: Downregulated DEGs", style = "heading 2")
sup <- add_ft(sup, cbind(Gene = rownames(head(down_degs, 50)), round(head(down_degs[, c("logFC", "AveExpr", "adj.P.Val")], 50), 4)), "Top 50 downregulated DEGs.")

# Supplementary Table S4: GO BP
sup <- officer::body_add_break(sup)
sup <- officer::body_add_par(sup, "Supplementary Table S4: Complete GO Biological Process", style = "heading 2")
go_disp <- go_bp[, c("Description", "Count", "pvalue", "p.adjust")]
go_disp$pvalue <- format(go_disp$pvalue, scientific = TRUE, digits = 3)
go_disp$p.adjust <- format(go_disp$p.adjust, scientific = TRUE, digits = 3)
sup <- add_ft(sup, go_disp, "All significantly enriched GO BP terms.")

# Supplementary Table S5: KEGG
sup <- officer::body_add_break(sup)
sup <- officer::body_add_par(sup, "Supplementary Table S5: Complete KEGG Enrichment", style = "heading 2")
kegg_disp <- kegg[, c("Description", "Count", "pvalue", "p.adjust")]
kegg_disp$pvalue <- format(kegg_disp$pvalue, scientific = TRUE, digits = 3)
kegg_disp$p.adjust <- format(kegg_disp$p.adjust, scientific = TRUE, digits = 3)
sup <- add_ft(sup, kegg_disp, "All significantly enriched KEGG pathways.")

# Supplementary Tables S6-S8: ML genes
sup <- officer::body_add_break(sup)
sup <- officer::body_add_par(sup, "Supplementary Table S6: LASSO-Selected Genes (38 genes)", style = "heading 2")
sup <- add_ft(sup, data.frame(Gene = lasso_genes), "")
sup <- officer::body_add_par(sup, "Supplementary Table S7: SVM-Selected Genes (76 genes)", style = "heading 2")
sup <- add_ft(sup, data.frame(Gene = svm_genes), "")
sup <- officer::body_add_par(sup, "Supplementary Table S8: Random Forest-Selected Genes (76 genes)", style = "heading 2")
sup <- add_ft(sup, data.frame(Gene = rf_genes), "")

# Supplementary Table S9: Immune
sup <- officer::body_add_break(sup)
sup <- officer::body_add_par(sup, "Supplementary Table S9: Complete Immune Infiltration Results", style = "heading 2")
if (nrow(immune) > 0) {
  imm_disp <- immune[order(immune$P_value), ]
  imm_disp$P_value <- format(imm_disp$P_value, scientific = TRUE, digits = 3)
  imm_disp$Adj_P <- format(imm_disp$Adj_P, scientific = TRUE, digits = 3)
  sup <- add_ft(sup, imm_disp, "Differential immune infiltration for all cell types.")
}

# Supplementary Figures
sup <- officer::body_add_break(sup)
sup <- officer::body_add_par(sup, "Supplementary Figures", style = "heading 1")

# Figure S1: MA plot
sup <- officer::body_add_par(sup, "Supplementary Figure S1: MA Plot", style = "heading 2")
sup <- add_img(sup, file.path(fig_dir, "Differential_Analysis", "Figure_2B_MA_plot.tiff"), "MA plot showing average expression vs log2 fold change. Red points = significant DEGs.", width = 5, height = 4)

# Figure S2: GO CC
sup <- officer::body_add_break(sup)
sup <- officer::body_add_par(sup, "Supplementary Figure S2: GO Cellular Component Enrichment", style = "heading 2")
sup <- add_img(sup, file.path(fig_dir, "Enrichment", "Figure_3A_GO_CC_dotplot.tiff"), "GO CC enrichment dot plot.", width = 5, height = 4)

# Figure S3: GO MF
sup <- officer::body_add_par(sup, "Supplementary Figure S3: GO Molecular Function Enrichment", style = "heading 2")
sup <- add_img(sup, file.path(fig_dir, "Enrichment", "Figure_3A_GO_MF_dotplot.tiff"), "GO MF enrichment dot plot.", width = 5, height = 4)

# Figure S4: GO combined
sup <- officer::body_add_break(sup)
sup <- officer::body_add_par(sup, "Supplementary Figure S4: Combined GO Enrichment", style = "heading 2")
sup <- add_img(sup, file.path(fig_dir, "Enrichment", "Figure_3AB_GO_combined.tiff"), "Combined GO dot plot and bar plot.", width = 5.5, height = 4.5)

# Figure S5: GO network
sup <- officer::body_add_par(sup, "Supplementary Figure S5: GO Enrichment Network", style = "heading 2")
sup <- add_img(sup, file.path(fig_dir, "Enrichment", "Figure_3D_GO_emapplot.tiff"), "GO network plot showing term relationships.", width = 5.5, height = 4.5)

# Figure S6: RF importance
sup <- officer::body_add_break(sup)
sup <- officer::body_add_par(sup, "Supplementary Figure S6: Random Forest Feature Importance", style = "heading 2")
sup <- add_img(sup, file.path(fig_dir, "Machine_Learning", "Figure_6A_RF_MeanDecreaseGini.tiff"), "RF mean decrease Gini.", width = 5, height = 4)
sup <- add_img(sup, file.path(fig_dir, "Machine_Learning", "Figure_6B_RF_MeanDecreaseAccuracy.tiff"), "RF mean decrease accuracy.", width = 5, height = 4)

# Figure S7: Traditional Venn
sup <- officer::body_add_par(sup, "Supplementary Figure S7: Traditional Venn Diagram", style = "heading 2")
sup <- add_img(sup, file.path(fig_dir, "Feature_Gene_Intersection", "Figure_6B_Feature_gene_venn_traditional.tiff"), "Venn diagram of gene intersections.", width = 5, height = 4.5)

# Figure S8-S9: SHAP force
sup <- officer::body_add_break(sup)
sup <- officer::body_add_par(sup, "Supplementary Figure S8: SHAP Force Plot (Single Sample)", style = "heading 2")
sup <- add_img(sup, file.path(fig_dir, "SHAP", "Figure_7F_SHAP_force_plot.tiff"), "SHAP force plot for a representative tumor sample.", width = 5.5, height = 3)
sup <- officer::body_add_par(sup, "Supplementary Figure S9: SHAP Force Plot (Multiple Samples)", style = "heading 2")
sup <- add_img(sup, file.path(fig_dir, "SHAP", "Figure_7F_SHAP_force_multi.tiff"), "SHAP force plot for multiple samples.", width = 5.5, height = 4)

# Figure S10: Immune correlation
sup <- officer::body_add_break(sup)
sup <- officer::body_add_par(sup, "Supplementary Figure S10: Immune Cell Correlation Heatmap", style = "heading 2")
sup <- add_img(sup, file.path(fig_dir, "Immune_Infiltration", "Figure_9C_Immune_cell_correlation.tiff"), "Correlation heatmap of immune cell types.", width = 5, height = 4.5)

# Figure S11-S12: PCA individual
sup <- officer::body_add_par(sup, "Supplementary Figure S11: PCA Before Batch Correction (by Dataset)", style = "heading 2")
sup <- add_img(sup, file.path(fig_dir, "PCA", "Figure_1A_PCA_before_batch_by_dataset.tiff"), "PCA colored by dataset before correction.", width = 5, height = 4)
sup <- officer::body_add_par(sup, "Supplementary Figure S12: PCA After Batch Correction (by Tissue)", style = "heading 2")
sup <- add_img(sup, file.path(fig_dir, "PCA", "Figure_1D_PCA_after_batch_by_group.tiff"), "PCA colored by tissue type after correction.", width = 5, height = 4)

# Figure S13: KEGG pathway maps
sup <- officer::body_add_break(sup)
sup <- officer::body_add_par(sup, "Supplementary Figure S13: KEGG Pathway Maps", style = "heading 2")
kegg_pngs <- list.files(file.path(fig_dir, "Enrichment"), pattern = "\\.png$", full.names = TRUE)
for (kp in kegg_pngs) {
  sup <- add_img(sup, kp, paste0("KEGG pathway: ", gsub(".png", "", basename(kp))), width = 5.5, height = 4.5)
  sup <- officer::body_add_break(sup)
}

# Figure S14: SVM accuracy
sup <- officer::body_add_par(sup, "Supplementary Figure S14: SVM Cross-Validation Accuracy", style = "heading 2")
sup <- add_img(sup, file.path(fig_dir, "Machine_Learning", "Figure_5B_SVM_CV_accuracy.tiff"), "SVM CV accuracy across folds.", width = 5, height = 4)

# Figure S15: Heatmap top 100
sup <- officer::body_add_break(sup)
sup <- officer::body_add_par(sup, "Supplementary Figure S15: Heatmap of Top 100 DEGs", style = "heading 2")
sup <- add_img(sup, file.path(fig_dir, "Differential_Analysis", "Figure_2D_Heatmap_top100_DEGs.tiff"), "Heatmap of top 100 DEGs.", width = 5.5, height = 4.5)

# Figure S16: GO cnetplot
sup <- officer::body_add_break(sup)
sup <- officer::body_add_par(sup, "Supplementary Figure S16: GO cnetplot", style = "heading 2")
cnet <- file.path(fig_dir, "Enrichment", "Figure_3D_GO_cnetplot.tiff")
if (file.exists(cnet)) {
  sup <- add_img(sup, cnet, "GO cnetplot showing gene-pathway relationships.", width = 5.5, height = 4.5)
}

# Batch info table
sup <- officer::body_add_break(sup)
sup <- officer::body_add_par(sup, "Supplementary Table S10: Batch Information", style = "heading 2")
sup <- add_ft(sup, batch, "Sample distribution across batches.")

# Save supplementary
supp_path <- file.path(report_dir, "Breast_Cancer_Supplementary_Materials.docx")
print(sup, target = supp_path)
cat("Supplementary file saved to:", supp_path, "\n")

cat("\n========================================================\n")
cat("MANUSCRIPT AND SUPPLEMENTARY FILES GENERATED\n")
cat("========================================================\n")
cat("Main manuscript:", manuscript_path, "\n")
cat("Supplementary:", supp_path, "\n")
