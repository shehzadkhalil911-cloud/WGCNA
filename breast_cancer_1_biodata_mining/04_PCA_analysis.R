###############################################################################
# Script 04: PCA Analysis - Before and After Batch Correction
###############################################################################

set.seed(123)

library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(patchwork)

base_dir <- normalizePath("breast_cancer_integrated_analysis")
data_processed_dir <- file.path(base_dir, "data_processed")
fig_dir <- file.path(base_dir, "results", "figures", "PCA")
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
  
  cat("  Saved:", filepath_tiff, "\n")
  cat("  Saved:", filepath_pdf, "\n")
}

# --- Helper: theme for publication ---
theme_publication <- theme_classic() +
  theme(
    text = element_text(size = 12, family = "sans"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(size = 11, face = "bold")
  )

# --- Helper: Run PCA ---
run_pca <- function(exprs_mat, scale = TRUE) {
  # Transpose: samples in rows, genes in columns
  mat <- t(exprs_mat)
  # Remove any columns with NA or zero variance
  keep_cols <- apply(mat, 2, function(x) {
    !any(is.na(x)) && var(x, na.rm = TRUE) > 0
  })
  mat <- mat[, keep_cols]
  
  pca_result <- prcomp(mat, center = TRUE, scale. = scale)
  var_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 2)
  
  list(
    pca = pca_result,
    var_explained = var_explained,
    x = as.data.frame(pca_result$x),
    samples_used = nrow(mat),
    genes_used = ncol(mat)
  )
}

# --- Helper: PCA plot ---
plot_pca <- function(pca_obj, metadata, color_var, color_name, 
                     title = "", ellipse = TRUE) {
  
  pca_df <- pca_obj$x
  pca_df[[color_var]] <- as.factor(metadata[[color_var]])
  
  var_ex <- pca_obj$var_explained
  
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[color_var]])) +
    geom_point(size = 2.5, alpha = 0.8) +
    labs(
      title = title,
      x = paste0("PC1 (", var_ex[1], "%)"),
      y = paste0("PC2 (", var_ex[2], "%)"),
      color = color_name
    ) +
    theme_publication +
    scale_color_viridis_d(option = "D", end = 0.9) +
    coord_fixed(ratio = 1, xlim = NULL, ylim = NULL)
  
  if (ellipse && length(unique(pca_df[[color_var]])) > 1) {
    p <- p + stat_ellipse(aes(fill = .data[[color_var]]), 
                          alpha = 0.1, level = 0.95, show.legend = FALSE) +
      scale_fill_viridis_d(option = "D", end = 0.9, guide = "none")
  }
  
  p
}

# --- Load data ---
cat("Loading data...\n")
exprs_before <- readRDS(file.path(data_processed_dir, "merged_exprs_before_batch.rds"))
exprs_after <- readRDS(file.path(data_processed_dir, "merged_exprs_after_batch.rds"))
pheno <- readRDS(file.path(data_processed_dir, "merged_pheno_with_batch.rds"))

cat("Expression matrix before batch:", dim(exprs_before), "\n")
cat("Expression matrix after batch:", dim(exprs_after), "\n")

# --- Run PCA before batch correction ---
cat("\nRunning PCA before batch correction...\n")
pca_before <- run_pca(exprs_before)
cat("  Samples used:", pca_before$samples_used, "\n")
cat("  Genes used:", pca_before$genes_used, "\n")

pca_before_df <- cbind(
  pca_before$x,
  Batch = pheno$batch,
  Group = pheno$group
)

# Save PCA coordinates
write.csv(pca_before_df, 
          file.path(tables_dir, "PCA_before_batch_coordinates.csv"), row.names = TRUE)

# --- Run PCA after batch correction ---
cat("Running PCA after batch correction...\n")
pca_after <- run_pca(exprs_after)
cat("  Samples used:", pca_after$samples_used, "\n")
cat("  Genes used:", pca_after$genes_used, "\n")

pca_after_df <- cbind(
  pca_after$x,
  Batch = pheno$batch,
  Group = pheno$group
)
write.csv(pca_after_df, 
          file.path(tables_dir, "PCA_after_batch_coordinates.csv"), row.names = TRUE)

# --- Generate PCA plots ---

# 1. PCA before batch - colored by dataset/batch
cat("\nGenerating PCA plots...\n")
p1 <- plot_pca(pca_before, pheno, "batch", "Dataset",
               title = "PCA Before Batch Correction\n(Colored by Dataset)",
               ellipse = TRUE)

# 2. PCA before batch - colored by group
p2 <- plot_pca(pca_before, pheno, "group", "Tissue Type",
               title = "PCA Before Batch Correction\n(Colored by Tissue Type)",
               ellipse = TRUE)

# 3. PCA after batch - colored by dataset/batch
p3 <- plot_pca(pca_after, pheno, "batch", "Dataset",
               title = "PCA After Batch Correction\n(Colored by Dataset)",
               ellipse = TRUE)

# 4. PCA after batch - colored by group
p4 <- plot_pca(pca_after, pheno, "group", "Tissue Type",
               title = "PCA After Batch Correction\n(Colored by Tissue Type)",
               ellipse = TRUE)

# 5. Combined before-after figure
p5 <- (p1 + p2) / (p3 + p4) +
  plot_annotation(
    title = "Principal Component Analysis of Integrated Breast Cancer Datasets",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

# --- Save all figures ---
cat("\nSaving figures...\n")

save_publication_figure(p1, "Figure_1A_PCA_before_batch_by_dataset", 
                        width = 8, height = 6)
save_publication_figure(p2, "Figure_1B_PCA_before_batch_by_group", 
                        width = 8, height = 6)
save_publication_figure(p3, "Figure_1C_PCA_after_batch_by_dataset", 
                        width = 8, height = 6)
save_publication_figure(p4, "Figure_1D_PCA_after_batch_by_group", 
                        width = 8, height = 6)
save_publication_figure(p5, "Figure_1E_PCA_combined_before_after", 
                        width = 14, height = 12)

cat("\nAll PCA figures saved to:", fig_dir, "\n")

# --- Save PCA objects for downstream use ---
saveRDS(pca_before, file.path(data_processed_dir, "pca_before_batch.rds"))
saveRDS(pca_after, file.path(data_processed_dir, "pca_after_batch.rds"))

cat("\nPCA analysis complete.\n")
