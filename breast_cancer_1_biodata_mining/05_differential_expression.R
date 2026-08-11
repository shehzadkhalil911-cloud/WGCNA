###############################################################################
# Script 05: Differential Expression Analysis using limma
###############################################################################

set.seed(123)

library(limma)
library(ggplot2)
library(tidyverse)
library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
library(viridis)

base_dir <- normalizePath("breast_cancer_integrated_analysis")
data_processed_dir <- file.path(base_dir, "data_processed")
fig_dir <- file.path(base_dir, "results", "figures", "Differential_Analysis")
tables_dir <- file.path(base_dir, "results", "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helper: Save 600 dpi figures ---
save_publication_figure <- function(plot, filename, width = 8, height = 6) {
  filepath_tiff <- file.path(fig_dir, paste0(filename, ".tiff"))
  filepath_pdf <- file.path(fig_dir, paste0(filename, ".pdf"))
  ggsave(filepath_tiff, plot = plot, width = width, height = height, 
         dpi = 600, compression = "lzw", device = "tiff")
  ggsave(filepath_pdf, plot = plot, width = width, height = height, 
         dpi = 600, device = "pdf")
  cat("  Saved figures\n")
}

# --- Publication theme ---
theme_publication <- theme_classic() +
  theme(
    text = element_text(size = 12, family = "sans"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

# --- Load data ---
cat("Loading batch-corrected expression data...\n")
exprs_after <- readRDS(file.path(data_processed_dir, "merged_exprs_after_batch.rds"))
pheno <- readRDS(file.path(data_processed_dir, "merged_pheno_with_batch.rds"))

cat("Expression matrix:", dim(exprs_after), "\n")
cat("Samples:", ncol(exprs_after), "\n")
cat("Genes:", nrow(exprs_after), "\n")

# --- Prepare design matrix ---
group <- factor(pheno$group, levels = c("Normal", "Tumor"))
batch <- factor(pheno$batch)

# Design with group as main factor, batch as covariate
design <- model.matrix(~ 0 + group + batch)
colnames(design) <- make.names(colnames(design))

cat("\nDesign matrix columns:", colnames(design), "\n")
print(table(group, batch))

# --- Fit linear model ---
cat("\nFitting linear model...\n")
fit <- lmFit(exprs_after, design)

# Contrast: Tumor vs Normal
contrast_matrix <- makeContrasts(Tumor_vs_Normal = groupTumor - groupNormal, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, robust = TRUE)

# --- Extract results ---
deg_table <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "P")
cat("\nTotal DEGs tested:", nrow(deg_table), "\n")

# Add direction
deg_table$Direction <- ifelse(deg_table$logFC > 0, "Up", "Down")

# Save full DEG table
write.csv(deg_table, file.path(tables_dir, "DEG_full_table.csv"))
cat("Full DEG table saved:", nrow(deg_table), "genes\n")

# --- Significant DEGs (adj.P < 0.05, |logFC| > 1) ---
pval_cutoff <- 0.05
lfc_cutoff <- 1

sig_degs <- deg_table[deg_table$adj.P.Val < pval_cutoff & 
                        abs(deg_table$logFC) > lfc_cutoff, ]
cat("\nSignificant DEGs (adj.P <", pval_cutoff, ", |logFC| >", lfc_cutoff, "):", 
    nrow(sig_degs), "\n")

if (nrow(sig_degs) > 0) {
  sig_degs <- sig_degs[order(sig_degs$adj.P.Val), ]
  
  # Up and down regulated
  up_degs <- sig_degs[sig_degs$logFC > 0, ]
  down_degs <- sig_degs[sig_degs$logFC < 0, ]
  
  cat("  Upregulated:", nrow(up_degs), "\n")
  cat("  Downregulated:", nrow(down_degs), "\n")
  
  # Save tables
  write.csv(sig_degs, file.path(tables_dir, "DEG_significant.csv"))
  write.csv(up_degs, file.path(tables_dir, "DEG_upregulated.csv"))
  write.csv(down_degs, file.path(tables_dir, "DEG_downregulated.csv"))
  
  # Save DEG gene list for downstream analysis
  deg_genes <- rownames(sig_degs)
  saveRDS(deg_genes, file.path(data_processed_dir, "deg_genes.rds"))
  writeLines(deg_genes, file.path(tables_dir, "DEG_gene_list.txt"))
} else {
  cat("WARNING: No significant DEGs found at current thresholds.\n")
  cat("Trying relaxed threshold (|logFC| > 0.5)...\n")
  
  sig_degs <- deg_table[deg_table$adj.P.Val < pval_cutoff & 
                          abs(deg_table$logFC) > 0.5, ]
  cat("Significant DEGs (adj.P <", pval_cutoff, ", |logFC| > 0.5):", nrow(sig_degs), "\n")
  
  if (nrow(sig_degs) > 0) {
    up_degs <- sig_degs[sig_degs$logFC > 0, ]
    down_degs <- sig_degs[sig_degs$logFC < 0, ]
    write.csv(sig_degs, file.path(tables_dir, "DEG_significant_relaxed.csv"))
  }
}

# --- Volcano plot ---
cat("\nGenerating volcano plot...\n")
key_genes <- c()
if (nrow(sig_degs) > 10) {
  key_genes <- rownames(head(sig_degs[order(abs(sig_degs$logFC), decreasing = TRUE), ], 10))
} else if (nrow(sig_degs) > 0) {
  key_genes <- rownames(sig_degs)
}

volcano_plot <- EnhancedVolcano(
  deg_table,
  lab = rownames(deg_table),
  x = "logFC",
  y = "adj.P.Val",
  selectLab = key_genes,
  pCutoff = pval_cutoff,
  FCcutoff = lfc_cutoff,
  title = "Differential Expression: Tumor vs Normal",
  subtitle = paste0("Cutoffs: |log2FC| > ", lfc_cutoff, ", adj.P < ", pval_cutoff),
  caption = paste0("Total: ", nrow(sig_degs), " DEGs"),
  pointSize = 2.0,
  labSize = 3.5,
  colAlpha = 0.6,
  legendPosition = "right",
  legendLabSize = 10,
  legendIconSize = 4.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  max.overlaps = 15
) + theme_publication

ggsave(file.path(fig_dir, "Figure_2A_Volcano_DEGs.tiff"), 
       plot = volcano_plot, width = 10, height = 8, dpi = 600, compression = "lzw")
ggsave(file.path(fig_dir, "Figure_2A_Volcano_DEGs.pdf"), 
       plot = volcano_plot, width = 10, height = 8, dpi = 600)
cat("  Volcano plot saved\n")

# --- MA plot ---
cat("Generating MA plot...\n")
deg_table$AveExpr <- deg_table$AveExpr %||% rowMeans(exprs_after[rownames(deg_table), ])
ma_data <- deg_table
ma_data$Significant <- ifelse(ma_data$adj.P.Val < pval_cutoff & 
                                abs(ma_data$logFC) > lfc_cutoff, "Yes", "No")

ma_plot <- ggplot(ma_data, aes(x = AveExpr, y = logFC, color = Significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", 
             color = "red", linewidth = 0.5) +
  scale_color_manual(values = c("Yes" = "red", "No" = "grey50")) +
  labs(title = "MA Plot: Tumor vs Normal",
       x = "Average Expression",
       y = "Log2 Fold Change") +
  theme_publication

ggsave(file.path(fig_dir, "Figure_2B_MA_plot.tiff"), 
       plot = ma_plot, width = 8, height = 6, dpi = 600, compression = "lzw")
ggsave(file.path(fig_dir, "Figure_2B_MA_plot.pdf"), 
       plot = ma_plot, width = 8, height = 6, dpi = 600)

# --- Heatmap of top DEGs ---
if (nrow(sig_degs) > 0) {
  cat("Generating heatmaps...\n")
  
  n_top <- min(50, nrow(sig_degs))
  top_genes <- rownames(sig_degs)[1:n_top]
  
  # Prepare expression data for heatmap
  heatmap_data <- exprs_after[top_genes, ]
  
  # Annotation for samples
  annotation_col <- data.frame(
    Dataset = pheno$batch,
    Group = pheno$group,
    row.names = colnames(exprs_after)
  )
  
  annotation_colors <- list(
    Dataset = c(
      GSE42568 = "#1f78b4",
      GSE15852 = "#33a02c",
      GSE29431 = "#e31a1c",
      GSE70947 = "#ff7f00"
    ),
    Group = c(
      Normal = "#4daf4a",
      Tumor = "#e41a1c"
    )
  )
  
  # Heatmap of top 50 DEGs
  htmap <- pheatmap(
    heatmap_data,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    scale = "row",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 6,
    fontsize_col = 6,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    main = paste0("Top ", n_top, " Differentially Expressed Genes"),
    filename = file.path(fig_dir, "Figure_2C_Heatmap_top50_DEGs.tiff"),
    width = 10,
    height = 10,
    resolution = 600
  )
  
  # Also save as PDF
  pheatmap(
    heatmap_data,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    scale = "row",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 6,
    fontsize_col = 6,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    main = paste0("Top ", n_top, " Differentially Expressed Genes"),
    filename = file.path(fig_dir, "Figure_2C_Heatmap_top50_DEGs.pdf"),
    width = 10,
    height = 10
  )
  
  # Heatmap of top 100 if available
  if (nrow(sig_degs) >= 100) {
    n_top100 <- min(100, nrow(sig_degs))
    top_genes100 <- rownames(sig_degs)[1:n_top100]
    
    pheatmap(
      exprs_after[top_genes100, ],
      annotation_col = annotation_col,
      annotation_colors = annotation_colors,
      scale = "row",
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      clustering_method = "ward.D2",
      show_rownames = FALSE,
      show_colnames = FALSE,
      fontsize_row = 4,
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      main = paste0("Top ", n_top100, " Differentially Expressed Genes"),
      filename = file.path(fig_dir, "Figure_2D_Heatmap_top100_DEGs.tiff"),
      width = 12,
      height = 10,
      resolution = 600
    )
  }
} else {
  cat("No DEGs found - skipping heatmaps.\n")
}

# --- Save limma results for downstream ---
saveRDS(deg_table, file.path(data_processed_dir, "limma_deg_table.rds"))
saveRDS(fit2, file.path(data_processed_dir, "limma_fit.rds"))

cat("\n========================================================\n")
cat("DIFFERENTIAL EXPRESSION ANALYSIS COMPLETE\n")
cat("========================================================\n")
cat("Total DEGs tested:", nrow(deg_table), "\n")
if (exists("sig_degs")) cat("Significant DEGs:", nrow(sig_degs), "\n")
cat("Results saved to:\n")
cat("  - results/tables/DEG_full_table.csv\n")
cat("  - results/tables/DEG_significant.csv\n")
cat("  - results/figures/Differential_Analysis/\n")
