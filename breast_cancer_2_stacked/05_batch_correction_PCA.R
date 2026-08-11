setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 05_batch_correction_PCA.R
# Discovery-only batch correction using ComBat and PCA visualization
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")
source("R/plotting_functions.R")
load_required_packages()

disc <- readRDS(file.path(DISCOVERY_DIR, "discovery_cohort.rds"))
expr <- disc$expression
meta <- disc$metadata

cat("Discovery expression matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")
cat("Batch (dataset) distribution:\n")
print(table(meta$dataset))
cat("Class distribution:\n")
print(table(meta$class))

# --- PCA before batch correction ---
pca_before <- prcomp(t(expr), center = TRUE, scale. = TRUE)
var_exp <- round(100 * pca_before$sdev^2 / sum(pca_before$sdev^2), 2)
cat("\nVariance explained by PC1:", var_exp[1], "%\n")
cat("Variance explained by PC2:", var_exp[2], "%\n")

pca_df <- data.frame(
    PC1 = pca_before$x[, 1],
    PC2 = pca_before$x[, 2],
    Dataset = meta$dataset,
    Class = meta$class
)

p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Dataset)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(title = "PCA Before Batch Correction (by Dataset)",
         x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)")) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Class)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = c(Tumor = TUMOR_COLOR, Normal = NORMAL_COLOR)) +
    labs(title = "PCA Before Batch Correction (by Class)",
         x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)")) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

save_publication_plot(p1, "PCA_before_correction_by_dataset", width = 8, height = 6)
save_publication_plot(p2, "PCA_before_correction_by_class", width = 8, height = 6)

# --- ComBat batch correction ---
cat("\nApplying ComBat batch correction...\n")
design <- model.matrix(~ class, data = meta)
combat_expr <- sva::ComBat(
    dat = as.matrix(expr),
    batch = meta$dataset,
    mod = design,
    par.prior = TRUE,
    prior.plots = FALSE
)

cat("Batch correction complete.\n")

# --- PCA after batch correction ---
pca_after <- prcomp(t(combat_expr), center = TRUE, scale. = TRUE)
var_exp_after <- round(100 * pca_after$sdev^2 / sum(pca_after$sdev^2), 2)

pca_df_after <- data.frame(
    PC1 = pca_after$x[, 1],
    PC2 = pca_after$x[, 2],
    Dataset = meta$dataset,
    Class = meta$class
)

p3 <- ggplot(pca_df_after, aes(x = PC1, y = PC2, color = Dataset)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(title = "PCA After Batch Correction (by Dataset)",
         x = paste0("PC1 (", var_exp_after[1], "%)"),
         y = paste0("PC2 (", var_exp_after[2], "%)")) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p4 <- ggplot(pca_df_after, aes(x = PC1, y = PC2, color = Class)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = c(Tumor = TUMOR_COLOR, Normal = NORMAL_COLOR)) +
    labs(title = "PCA After Batch Correction (by Class)",
         x = paste0("PC1 (", var_exp_after[1], "%)"),
         y = paste0("PC2 (", var_exp_after[2], "%)")) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

save_publication_plot(p3, "PCA_after_correction_by_dataset", width = 8, height = 6)
save_publication_plot(p4, "PCA_after_correction_by_class", width = 8, height = 6)

# Combined before/after
p5 <- p1 + p3 + plot_layout(ncol = 2) + 
    plot_annotation(title = "PCA Before vs After Batch Correction")
save_publication_plot(p5, "PCA_before_vs_after_combined", width = 16, height = 7)

# Save batch-corrected data
saveRDS(combat_expr, file.path(DISCOVERY_DIR, "discovery_combat.rds"))
saveRDS(pca_before, file.path(PCA_DIR, "pca_before.rds"))
saveRDS(pca_after, file.path(PCA_DIR, "pca_after.rds"))

cat("\nBatch correction and PCA complete.\n")

