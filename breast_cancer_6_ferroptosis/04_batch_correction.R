source("R/utils.R")
load_libraries()
set.seed(123)

BASE_DIR <- "F:/independent projects/cancer/breast cancer"
OUT_DIR <- file.path(BASE_DIR, "Breast cancer 6 Ferroptosis")

log_msg("=== Step 4: Discovery Cohort Integration and Batch Correction ===")

discovery_list <- readRDS(file.path(OUT_DIR, "data_processed/discovery/discovery_list.rds"))
ferr_genes <- readRDS(file.path(OUT_DIR, "data_processed/ferroptosis_gene_list/ferroptosis_genes_final.rds"))

# Get common genes across discovery datasets
gene_sets <- lapply(discovery_list, function(d) rownames(d$expr))
common_genes <- Reduce(intersect, gene_sets)
log_msg(sprintf("Common genes across discovery: %d", length(common_genes)))

# Merge expression matrices
merged_expr <- NULL
merged_metadata <- NULL
batch_vec <- c()

for (ds in names(discovery_list)) {
  d <- discovery_list[[ds]]
  expr_sub <- d$expr[common_genes, , drop = FALSE]
  merged_expr <- cbind(merged_expr, expr_sub)
  batch_vec <- c(batch_vec, rep(ds, ncol(expr_sub)))
  meta_df <- data.frame(
    sample_id = colnames(expr_sub),
    dataset = ds,
    group = d$group,
    stringsAsFactors = FALSE
  )
  merged_metadata <- rbind(merged_metadata, meta_df)
}

log_msg(sprintf("Merged expression: %d genes x %d samples", nrow(merged_expr), nrow(merged_metadata)))
log_msg(sprintf("Batch distribution: %s", paste(table(batch_vec), collapse = ", ")))

# Check confounding between batch and phenotype
confusion <- table(merged_metadata$dataset, merged_metadata$group)
log_msg(sprintf("Batch x Group contingency:\n%s", paste(capture.output(print(confusion)), collapse = "\n")))

# ComBat batch correction
library(sva)

batch <- as.factor(batch_vec)
mod <- model.matrix(~ group, data = merged_metadata)

combat_expr <- tryCatch({
  ComBat(dat = merged_expr, batch = batch, mod = mod, par.prior = TRUE)
}, error = function(e) {
  log_msg(sprintf("ComBat parametric failed: %s", e$message))
  log_msg("Trying non-parametric ComBat...")
  ComBat(dat = merged_expr, batch = batch, mod = mod, par.prior = FALSE)
})

log_msg("ComBat batch correction completed")

# PCA before and after
pca_before <- prcomp(t(merged_expr), scale. = TRUE, center = TRUE)
pca_after <- prcomp(t(combat_expr), scale. = TRUE, center = TRUE)
pct_var_before <- summary(pca_before)$importance[2, 1:6] * 100
pct_var_after <- summary(pca_after)$importance[2, 1:6] * 100

pca_df_before <- data.frame(
  PC1 = pca_before$x[, 1], PC2 = pca_before$x[, 2],
  dataset = merged_metadata$dataset,
  group = merged_metadata$group
)
pca_df_after <- data.frame(
  PC1 = pca_after$x[, 1], PC2 = pca_after$x[, 2],
  dataset = merged_metadata$dataset,
  group = merged_metadata$group
)

# PCA before correction by dataset
p1 <- ggplot(pca_df_before, aes(PC1, PC2, color = dataset)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "PCA Before Correction", subtitle = "Colored by Dataset",
       x = sprintf("PC1 (%.1f%%)", pct_var_before[1]),
       y = sprintf("PC2 (%.1f%%)", pct_var_before[2])) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
save_tiff_600dpi(p1, file.path(OUT_DIR, "results/figures/Figure2A_PCA_before_dataset.tiff"))

p2 <- ggplot(pca_df_after, aes(PC1, PC2, color = dataset)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "PCA After Correction", subtitle = "Colored by Dataset",
       x = sprintf("PC1 (%.1f%%)", pct_var_after[1]),
       y = sprintf("PC2 (%.1f%%)", pct_var_after[2])) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
save_tiff_600dpi(p2, file.path(OUT_DIR, "results/figures/Figure2B_PCA_after_dataset.tiff"))

# PCA by group after correction
p3 <- ggplot(pca_df_after, aes(PC1, PC2, color = group)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(title = "PCA After Correction", subtitle = "Colored by Group",
       x = sprintf("PC1 (%.1f%%)", pct_var_after[1]),
       y = sprintf("PC2 (%.1f%%)", pct_var_after[2])) +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
save_tiff_600dpi(p3, file.path(OUT_DIR, "results/figures/Figure2C_PCA_after_group.tiff"))

# Density plots
expr_df_before <- data.frame(
  expression = c(merged_expr[, sample(1:ncol(merged_expr), min(20, ncol(merged_expr)))]),
  type = "Before correction"
)
expr_df_after <- data.frame(
  expression = c(combat_expr[, sample(1:ncol(combat_expr), min(20, ncol(combat_expr)))]),
  type = "After correction"
)
expr_density <- rbind(expr_df_before, expr_df_after)

p4 <- ggplot(expr_density, aes(expression, color = type)) +
  geom_density(alpha = 0.5) +
  labs(title = "Expression Density: Before vs After Correction") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
save_tiff_600dpi(p4, file.path(OUT_DIR, "results/figures/Figure2D_Density_comparison.tiff"))

# Sample correlation heatmaps
cor_before <- cor(merged_expr, method = "spearman")
cor_after <- cor(combat_expr, method = "spearman")

anno_col <- data.frame(
  Dataset = merged_metadata$dataset,
  Group = merged_metadata$group,
  stringsAsFactors = FALSE
)
rownames(anno_col) <- colnames(merged_expr)
anno_colors <- list(
  Dataset = c(GSE42568 = "#E41A1C", GSE15852 = "#377EB8", GSE29431 = "#4DAF4A"),
  Group = c(Tumor = "black", Normal = "white")
)

open_tiff(file.path(OUT_DIR, "results/figures/Figure2E_Correlation_before.tiff"), width = 8, height = 7)
pheatmap::pheatmap(cor_before, annotation_col = anno_col, annotation_colors = anno_colors,
                   main = "Sample Correlation Before Correction",
                   show_rownames = FALSE, show_colnames = FALSE,
                   fontsize = 6, cutree_cols = 3)
dev.off()

open_tiff(file.path(OUT_DIR, "results/figures/Figure2F_Correlation_after.tiff"), width = 8, height = 7)
pheatmap::pheatmap(cor_after, annotation_col = anno_col, annotation_colors = anno_colors,
                   main = "Sample Correlation After Correction",
                   show_rownames = FALSE, show_colnames = FALSE,
                   fontsize = 6, cutree_cols = 3)
dev.off()

# Save batch-corrected data
saveRDS(combat_expr, file.path(OUT_DIR, "data_processed/discovery/discovery_combat.rds"))
saveRDS(merged_expr, file.path(OUT_DIR, "data_processed/discovery/discovery_raw_merged.rds"))
write.csv(merged_metadata, file.path(OUT_DIR, "data_processed/discovery/discovery_metadata.csv"), row.names = FALSE)
write.csv(combat_expr, file.path(OUT_DIR, "data_processed/discovery/discovery_combat_expr.csv"))

# PCA data
pca_data <- list(
  before = pca_df_before, after = pca_df_after,
  var_before = pct_var_before, var_after = pct_var_after
)
saveRDS(pca_data, file.path(OUT_DIR, "results/batch_correction/pca_results.rds"))
write.csv(pca_df_before, file.path(OUT_DIR, "results/batch_correction/pca_before.csv"), row.names = FALSE)
write.csv(pca_df_after, file.path(OUT_DIR, "results/batch_correction/pca_after.csv"), row.names = FALSE)

# Batch effect assessment via ANOVA
batch_effect_before <- sapply(1:min(1000, nrow(merged_expr)), function(i) {
  summary(aov(merged_expr[i, ] ~ batch))[[1]]$`Pr(>F)`[1]
})
batch_effect_after <- sapply(1:min(1000, nrow(combat_expr)), function(i) {
  summary(aov(combat_expr[i, ] ~ batch))[[1]]$`Pr(>F)`[1]
})
batch_assessment <- data.frame(
  metric = c("Median P-value (batch effect)", "% significant batch (p<0.05)"),
  before = c(median(batch_effect_before, na.rm = TRUE),
             mean(batch_effect_before < 0.05, na.rm = TRUE) * 100),
  after = c(median(batch_effect_after, na.rm = TRUE),
            mean(batch_effect_after < 0.05, na.rm = TRUE) * 100)
)
write.csv(batch_assessment, file.path(OUT_DIR, "results/batch_correction/batch_effect_assessment.csv"), row.names = FALSE)
log_msg(sprintf("Batch effect assessment done: before median p=%.4f, after median p=%.4f",
                batch_assessment$before[1], batch_assessment$after[1]))

log_msg("Step 4 complete: Batch correction and integration finished")
