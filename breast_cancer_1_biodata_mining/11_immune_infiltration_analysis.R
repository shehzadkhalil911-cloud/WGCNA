###############################################################################
# Script 11: Immune Cell Infiltration Analysis
# Using xCell for immune cell deconvolution
###############################################################################

set.seed(123)

library(xCell)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(corrplot)

base_dir <- normalizePath("breast_cancer_integrated_analysis")
data_processed_dir <- file.path(base_dir, "data_processed")
fig_dir <- file.path(base_dir, "results", "figures", "Immune_Infiltration")
tables_dir <- file.path(base_dir, "results", "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helper: Save 600 dpi ---
save_publication_figure <- function(plot, filename, width = 9, height = 7) {
  ggsave(file.path(fig_dir, paste0(filename, ".tiff")), plot = plot,
         width = width, height = height, dpi = 600, compression = "lzw")
  ggsave(file.path(fig_dir, paste0(filename, ".pdf")), plot = plot,
         width = width, height = height, dpi = 600)
}

theme_publication <- theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  )

# --- Load data ---
cat("Loading expression data...\n")
exprs_after <- readRDS(file.path(data_processed_dir, "merged_exprs_after_batch.rds"))
pheno <- readRDS(file.path(data_processed_dir, "merged_pheno_with_batch.rds"))
final_genes <- readRDS(file.path(data_processed_dir, "final_feature_genes.rds"))

cat("Expression matrix:", dim(exprs_after), "\n")
cat("Samples:", ncol(exprs_after), "\n")

# ================================================================
# Immune Cell Infiltration using xCell
# ================================================================
cat("\n========================================================\n")
cat("xCell Immune Cell Deconvolution\n")
cat("========================================================\n")

# xCell requires gene-level expression matrix
# Need to ensure gene symbols are valid
exprs_xcell <- exprs_after
rownames(exprs_xcell) <- make.names(rownames(exprs_xcell), unique = TRUE)

# Run xCell
cat("Running xCell (this may take a few minutes)...\n")
xcell_results <- tryCatch({
  xCellAnalysis(exprs_xcell, 
                rnaseq = FALSE,
                file.name = NULL,
                scale = TRUE,
                alpha = 0.5,
                save.raw = FALSE)
}, error = function(e) {
  cat("xCell error:", e$message, "\n")
  cat("Trying with adjusted parameters...\n")
  xCellAnalysis(exprs_xcell, 
                rnaseq = FALSE,
                scale = FALSE,
                alpha = 0.5)
})

cat("xCell result dimensions:", dim(xcell_results), "\n")
cat("Cell types analyzed:", nrow(xcell_results), "\n")

# Save xCell results
saveRDS(xcell_results, file.path(data_processed_dir, "xcell_results.rds"))
write.csv(xcell_results, file.path(tables_dir, "immune_cell_infiltration_xCell.csv"))

# ================================================================
# A. Box plot: Immune cell infiltration comparison
# ================================================================
cat("\nGenerating immune cell infiltration box plots...\n")

# Prepare data for plotting
immune_scores <- as.data.frame(t(xcell_results))
immune_scores$Sample <- rownames(immune_scores)
immune_scores$Group <- pheno$group[match(immune_scores$Sample, colnames(exprs_after))]

# Select immune cell types (exclude stroma and scores)
cell_types <- rownames(xcell_results)
# Filter for actual immune cell types (exclude scores like "ImmuneScore", "StromaScore")
immune_cells <- cell_types[!grepl("Score$", cell_types)]

if (length(immune_cells) > 0) {
  # Long format for ggplot
  plot_data <- immune_scores %>%
    pivot_longer(cols = all_of(immune_cells), 
                 names_to = "CellType", values_to = "Score") %>%
    filter(!is.na(Group))
  
  # Box plot
  p_immune_box <- ggplot(plot_data, aes(x = CellType, y = Score, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5, notch = FALSE) +
    scale_fill_manual(values = c("Normal" = "#4daf4a", "Tumor" = "#e41a1c")) +
    labs(title = "Immune Cell Infiltration: Tumor vs Normal",
         x = "Cell Type", y = "xCell Score") +
    theme_publication +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  
  save_publication_figure(p_immune_box, "Figure_9A_Immune_infiltration_boxplot",
                          width = 12, height = 7)
  
  # Perform statistical tests
  immune_stats <- plot_data %>%
    group_by(CellType) %>%
    summarise(
      Mean_Tumor = mean(Score[Group == "Tumor"], na.rm = TRUE),
      Mean_Normal = mean(Score[Group == "Normal"], na.rm = TRUE),
      P_value = tryCatch(
        t.test(Score[Group == "Tumor"], Score[Group == "Normal"])$p.value,
        error = function(e) NA
      ),
      .groups = "drop"
    )
  
  if (nrow(immune_stats) > 0) {
    immune_stats$Adj_P <- p.adjust(immune_stats$P_value, method = "BH")
    immune_stats <- immune_stats[order(immune_stats$P_value), ]
    write.csv(immune_stats, 
              file.path(tables_dir, "immune_cell_differential_stats.csv"), 
              row.names = FALSE)
    cat("Differential immune infiltration statistics saved.\n")
  }
}

# ================================================================
# B. Heatmap of immune cell infiltration
# ================================================================
cat("Generating immune infiltration heatmap...\n")

# Annotations
annotation_col <- data.frame(
  Group = pheno$group,
  row.names = colnames(exprs_after)
)
annotation_colors <- list(
  Group = c(Normal = "#4daf4a", Tumor = "#e41a1c")
)

# Filter to immune cells only
if (length(immune_cells) > 1) {
  pheatmap(
    xcell_results[immune_cells, ],
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    scale = "row",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 8,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    main = "Immune Cell Infiltration Profiles",
    filename = file.path(fig_dir, "Figure_9B_Immune_infiltration_heatmap.tiff"),
    width = 10,
    height = 7,
    resolution = 600
  )
  
  pheatmap(
    xcell_results[immune_cells, ],
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    scale = "row",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "ward.D2",
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 8,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    main = "Immune Cell Infiltration Profiles",
    filename = file.path(fig_dir, "Figure_9B_Immune_infiltration_heatmap.pdf"),
    width = 10,
    height = 7
  )
}

# ================================================================
# Immune Cell Correlation Heatmap
# ================================================================
cat("Generating immune cell correlation heatmap...\n")

if (length(immune_cells) > 2) {
  immune_corr <- cor(t(xcell_results[immune_cells, ]), method = "spearman",
                     use = "pairwise.complete.obs")
  
  pdf(file.path(fig_dir, "Figure_9C_Immune_cell_correlation.pdf"),
      width = 8, height = 7)
  corrplot(immune_corr, method = "color", type = "upper",
           order = "hclust", tl.col = "black", tl.cex = 0.6,
           col = colorRampPalette(c("#4daf4a", "white", "#e41a1c"))(100),
           title = "Immune Cell Type Correlations",
           mar = c(0, 0, 2, 0))
  dev.off()
  
  tiff(file.path(fig_dir, "Figure_9C_Immune_cell_correlation.tiff"),
       width = 8, height = 7, units = "in", res = 600, compression = "lzw")
  corrplot(immune_corr, method = "color", type = "upper",
           order = "hclust", tl.col = "black", tl.cex = 0.6,
           col = colorRampPalette(c("#4daf4a", "white", "#e41a1c"))(100),
           title = "Immune Cell Type Correlations",
           mar = c(0, 0, 2, 0))
  dev.off()
}

# ================================================================
# Gene-Immune Cell Associations (using final feature genes)
# ================================================================
cat("Generating gene-immune cell correlation heatmap...\n")

if (length(final_genes) > 0 && length(immune_cells) > 0) {
  # Get expression of final feature genes
  final_expr <- exprs_after[intersect(final_genes, rownames(exprs_after)), ]
  
  if (nrow(final_expr) > 0) {
    # Compute correlation between genes and immune cell scores
    gene_immune_corr <- cor(t(final_expr), t(xcell_results[immune_cells, ]),
                            method = "spearman", use = "pairwise.complete.obs")
    
    # Save correlation matrix
    write.csv(gene_immune_corr, 
              file.path(tables_dir, "gene_immune_correlation.csv"))
    
    # Heatmap
    pheatmap(
      gene_immune_corr,
      scale = "none",
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      clustering_method = "ward.D2",
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 8,
      fontsize_col = 8,
      color = colorRampPalette(c("#4daf4a", "white", "#e41a1c"))(100),
      main = "Feature Gene - Immune Cell Correlations",
      filename = file.path(fig_dir, "Figure_9D_Gene_immune_correlation_heatmap.tiff"),
      width = 8,
      height = max(6, nrow(gene_immune_corr) * 0.4),
      resolution = 600
    )
    
    pheatmap(
      gene_immune_corr,
      scale = "none",
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      clustering_method = "ward.D2",
      show_rownames = TRUE,
      show_colnames = TRUE,
      fontsize_row = 8,
      fontsize_col = 8,
      color = colorRampPalette(c("#4daf4a", "white", "#e41a1c"))(100),
      main = "Feature Gene - Immune Cell Correlations",
      filename = file.path(fig_dir, "Figure_9D_Gene_immune_correlation_heatmap.pdf"),
      width = 8,
      height = max(6, nrow(gene_immune_corr) * 0.4)
    )
    
    cat("  Gene-immune correlation heatmap saved.\n")
    
    # Create long-format correlation table
    gene_immune_long <- as.data.frame(as.table(gene_immune_corr))
    colnames(gene_immune_long) <- c("Gene", "Immune_Cell", "Correlation")
    gene_immune_long <- gene_immune_long[order(-abs(gene_immune_long$Correlation)), ]
    write.csv(gene_immune_long, 
              file.path(tables_dir, "gene_immune_correlation_long.csv"),
              row.names = FALSE)
  }
}

cat("\n========================================================\n")
cat("IMMUNE INFILTRATION ANALYSIS COMPLETE\n")
cat("========================================================\n")
cat("Results saved to results/figures/Immune_Infiltration/ and results/tables/\n")
