# =============================================================================
# Script 03: Discovery Cohort Integration, Batch Correction, and DE Analysis
# Breast Cancer 7 Cell Death Crosstalk
# =============================================================================

source(file.path("..", "config", "analysis_config.R"))

library(tidyverse)
library(sva)
library(limma)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(circlize)

log_message("Starting batch correction and differential expression analysis.")

# --- 1. LOAD PROCESSED DISCOVERY DATA ---
processed_list <- readRDS(file.path(OUTPUT$preprocessing, "all_processed_datasets.rds"))
discovery_datasets <- DISCOVERY_DATASETS

# --- 2. FIND COMMON GENES ---
common_genes <- NULL
for (ds in discovery_datasets) {
  common_genes <- if (is.null(common_genes)) {
    rownames(processed_list[[ds]]$expression)
  } else {
    intersect(common_genes, rownames(processed_list[[ds]]$expression))
  }
}
log_message(sprintf("Common genes across discovery datasets: %d", length(common_genes)))

# --- 3. MERGE DISCOVERY COHORT ---
discovery_expr_list <- list()
discovery_meta_list <- list()

for (i in seq_along(discovery_datasets)) {
  ds <- discovery_datasets[i]
  mat <- processed_list[[ds]]$expression[common_genes, , drop = FALSE]
  colnames(mat) <- paste0(ds, "_", colnames(mat))
  discovery_expr_list[[ds]] <- mat
  
  meta <- processed_list[[ds]]$metadata
  meta$merged_id <- paste0(ds, "_", meta$sample_id)
  meta$dataset <- ds
  discovery_meta_list[[ds]] <- meta
}

discovery_matrix <- do.call(cbind, discovery_expr_list)
discovery_metadata <- bind_rows(discovery_meta_list)

# Create group variable
discovery_metadata$group <- factor(discovery_metadata$class, levels = c("Normal", "Tumor"))
discovery_metadata$dataset <- factor(discovery_metadata$dataset, levels = discovery_datasets)

# Fix column order
discovery_matrix <- discovery_matrix[, discovery_metadata$merged_id]

log_message(sprintf("Discovery matrix: %d genes x %d samples", nrow(discovery_matrix), ncol(discovery_matrix)))

# --- 4. PCA BEFORE BATCH CORRECTION ---
pca_before <- prcomp(t(discovery_matrix), scale = TRUE, center = TRUE)
pca_before_df <- data.frame(
  PC1 = pca_before$x[, 1],
  PC2 = pca_before$x[, 2],
  Dataset = discovery_metadata$dataset,
  Group = discovery_metadata$group
)
pca_before_var <- summary(pca_before)$importance[2, 1:2] * 100

# By dataset
p_pca_before_ds <- ggplot(pca_before_df, aes(x = PC1, y = PC2, color = Dataset)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = CB_PALETTE[1:3]) +
  labs(title = "PCA Before Batch Correction (by Dataset)",
       x = sprintf("PC1 (%.1f%%)", pca_before_var[1]),
       y = sprintf("PC2 (%.1f%%)", pca_before_var[2])) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_tiff_600dpi(p_pca_before_ds, file.path(OUTPUT$batch_correction, "PCA_before_dataset.tiff"))

# By group
p_pca_before_grp <- ggplot(pca_before_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = c("Normal" = "#0072B2", "Tumor" = "#D55E00")) +
  labs(title = "PCA Before Batch Correction (by Group)",
       x = sprintf("PC1 (%.1f%%)", pca_before_var[1]),
       y = sprintf("PC2 (%.1f%%)", pca_before_var[2])) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_tiff_600dpi(p_pca_before_grp, file.path(OUTPUT$batch_correction, "PCA_before_group.tiff"))

# --- 5. BATCH CORRECTION WITH ComBat ---
design <- model.matrix(~ group, data = discovery_metadata)

combat_matrix <- sva::ComBat(
  dat = discovery_matrix,
  batch = discovery_metadata$dataset,
  mod = design,
  par.prior = TRUE,
  prior.plots = FALSE
)

log_message("ComBat batch correction completed.")

# --- 6. PCA AFTER BATCH CORRECTION ---
pca_after <- prcomp(t(combat_matrix), scale = TRUE, center = TRUE)
pca_after_df <- data.frame(
  PC1 = pca_after$x[, 1],
  PC2 = pca_after$x[, 2],
  Dataset = discovery_metadata$dataset,
  Group = discovery_metadata$group
)
pca_after_var <- summary(pca_after)$importance[2, 1:2] * 100

# By dataset
p_pca_after_ds <- ggplot(pca_after_df, aes(x = PC1, y = PC2, color = Dataset)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = CB_PALETTE[1:3]) +
  labs(title = "PCA After Batch Correction (by Dataset)",
       x = sprintf("PC1 (%.1f%%)", pca_after_var[1]),
       y = sprintf("PC2 (%.1f%%)", pca_after_var[2])) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_tiff_600dpi(p_pca_after_ds, file.path(OUTPUT$batch_correction, "PCA_after_dataset.tiff"))

# By group
p_pca_after_grp <- ggplot(pca_after_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  scale_color_manual(values = c("Normal" = "#0072B2", "Tumor" = "#D55E00")) +
  labs(title = "PCA After Batch Correction (by Group)",
       x = sprintf("PC1 (%.1f%%)", pca_after_var[1]),
       y = sprintf("PC2 (%.1f%%)", pca_after_var[2])) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_tiff_600dpi(p_pca_after_grp, file.path(OUTPUT$batch_correction, "PCA_after_group.tiff"))

# Combined before/after
pca_combined <- rbind(
  data.frame(pca_before_df, Batch = "Before"),
  data.frame(pca_after_df, Batch = "After")
)
p_combined <- ggplot(pca_combined, aes(x = PC1, y = PC2, color = Dataset, shape = Batch)) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_color_manual(values = CB_PALETTE[1:3]) +
  labs(title = "PCA Before vs After Batch Correction",
       x = sprintf("PC1 (%.1f%%)", mean(c(pca_before_var[1], pca_after_var[1]))),
       y = sprintf("PC2 (%.1f%%)", mean(c(pca_before_var[2], pca_after_var[2])))) +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_tiff_600dpi(p_combined, file.path(OUTPUT$batch_correction, "PCA_combined.tiff"))

# --- 7. SAVE COMBAT-CORRECTED DATA ---
saveRDS(combat_matrix, file.path(OUTPUT$batch_correction, "combat_corrected_matrix.rds"))
write.csv(combat_matrix, file.path(OUTPUT$batch_correction, "combat_corrected_matrix.csv"))
saveRDS(discovery_metadata, file.path(OUTPUT$batch_correction, "discovery_metadata.rds"))

# --- 8. DIFFERENTIAL EXPRESSION ANALYSIS (tumour vs normal) ---
design_de <- model.matrix(~ group, data = discovery_metadata)

fit <- lmFit(combat_matrix, design_de)
fit <- eBayes(fit)
de_results <- topTable(fit, coef = "groupTumor", number = Inf, sort.by = "P")

log_message(sprintf("DE analysis complete. %d genes tested.", nrow(de_results)))

# Save full results
write.csv(de_results, file.path(OUTPUT$de, "discovery_DEG_full.csv"))
saveRDS(de_results, file.path(OUTPUT$de, "discovery_DEG_full.rds"))

# Significant DEGs
deg_sig <- de_results[de_results$adj.P.Val < ADJ_PVAL_THRESHOLD & 
                        abs(de_results$logFC) >= LOG2FC_THRESHOLD, ]
deg_up <- deg_sig[deg_sig$logFC > 0, ]
deg_down <- deg_sig[deg_sig$logFC < 0, ]

log_message(sprintf("Significant DEGs: %d (Up: %d, Down: %d)", 
                    nrow(deg_sig), nrow(deg_up), nrow(deg_down)))

write.csv(deg_sig, file.path(OUTPUT$de, "discovery_DEG_significant.csv"))
write.csv(deg_up, file.path(OUTPUT$de, "discovery_DEG_upregulated.csv"))
write.csv(deg_down, file.path(OUTPUT$de, "discovery_DEG_downregulated.csv"))

# Ranked gene list (for GSEA)
de_results$rank_metric <- -log10(de_results$P.Value) * sign(de_results$logFC)
ranked_genes <- de_results[order(de_results$rank_metric, decreasing = TRUE), "rank_metric", drop = FALSE]
write.csv(ranked_genes, file.path(OUTPUT$de, "ranked_gene_list.csv"))

# --- 9. VOLCANO PLOT ---
p_volcano <- EnhancedVolcano(de_results,
                             lab = rownames(de_results),
                             x = "logFC",
                             y = "adj.P.Val",
                             title = "Discovery Cohort: Tumour vs Normal",
                             pCutoff = ADJ_PVAL_THRESHOLD,
                             FCcutoff = LOG2FC_THRESHOLD,
                             pointSize = 1.5,
                             labSize = 2.5,
                             colAlpha = 0.6,
                             legendPosition = "right",
                             legendLabSize = 10,
                             legendIconSize = 3.0,
                             col = c("grey30", "skyblue", "orange", "red3"),
                             cutoffLineType = "dashed",
                             cutoffLineCol = "black",
                             cutoffLineWidth = 0.5,
                             gridlines.major = FALSE,
                             gridlines.minor = FALSE,
                             border = "partial",
                             borderWidth = 1,
                             borderColour = "black")
save_tiff_600dpi(p_volcano, file.path(OUTPUT$de, "volcano.tiff"))

# --- 10. MA PLOT ---
de_results$AveExpr <- de_results$AveExpr
p_ma <- ggplot(de_results, aes(x = AveExpr, y = logFC, color = adj.P.Val < ADJ_PVAL_THRESHOLD & abs(logFC) >= LOG2FC_THRESHOLD)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_manual(values = c("grey60", "red3"), name = "Significant") +
  geom_hline(yintercept = c(-LOG2FC_THRESHOLD, 0, LOG2FC_THRESHOLD), 
             linetype = "dashed", color = "blue", linewidth = 0.3) +
  labs(title = "MA Plot: Discovery Cohort", x = "Average Expression", y = "log2 Fold Change") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_tiff_600dpi(p_ma, file.path(OUTPUT$de, "MA_plot.tiff"))

# --- 11. TOP DEG HEATMAPS ---
top_genes <- function(n, direction = "both") {
  if (direction == "both") {
    head(rownames(deg_sig[order(abs(deg_sig$logFC), decreasing = TRUE), ]), n)
  } else if (direction == "up") {
    head(rownames(deg_up[order(deg_up$logFC, decreasing = TRUE), ]), min(n, nrow(deg_up)))
  } else {
    head(rownames(deg_down[order(deg_down$logFC, decreasing = FALSE), ]), min(n, nrow(deg_down)))
  }
}

for (n in c(25, 50, 100)) {
  topn <- top_genes(n)
  if (length(topn) > 2) {
    hm_data <- combat_matrix[topn, ]
    hm_data <- t(scale(t(hm_data)))
    
    col_ann <- HeatmapAnnotation(
      Group = discovery_metadata$group,
      Dataset = discovery_metadata$dataset,
      col = list(
        Group = c(Normal = "#0072B2", Tumor = "#D55E00"),
        Dataset = c(GSE42568 = unname(CB_PALETTE[1]), GSE15852 = unname(CB_PALETTE[2]), GSE29431 = unname(CB_PALETTE[3]))
      )
    )
    
    ht <- Heatmap(hm_data,
                  name = "Z-score",
                  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
                  top_annotation = col_ann,
                  show_row_names = TRUE,
                  show_column_names = FALSE,
                  row_names_gp = gpar(fontsize = 6),
                  column_title = sprintf("Top %d DEGs", n),
                  clustering_distance_columns = "euclidean",
                  clustering_method_columns = "ward.D2")
    
    open_tiff(file.path(OUTPUT$de, sprintf("heatmap_top%d.tiff", n)))
    draw(ht)
    dev.off()
  }
}

# --- 12. DEG COUNT BAR PLOT ---
deg_counts <- data.frame(
  Direction = c("Upregulated", "Downregulated"),
  Count = c(nrow(deg_up), nrow(deg_down)),
  stringsAsFactors = FALSE
)

p_deg_counts <- ggplot(deg_counts, aes(x = Direction, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.5, color = "black") +
  geom_text(aes(label = Count), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("Upregulated" = "#D55E00", "Downregulated" = "#0072B2")) +
  labs(title = "Differentially Expressed Genes", y = "Count") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")
save_tiff_600dpi(p_deg_counts, file.path(OUTPUT$de, "DEG_counts.tiff"))

# --- 13. PER-DATASET DE ANALYSIS ---
per_dataset_deg <- list()
for (ds in discovery_datasets) {
  log_message(sprintf("Running per-dataset DE for %s...", ds))
  ds_samples <- discovery_metadata$merged_id[discovery_metadata$dataset == ds]
  ds_expr <- combat_matrix[, ds_samples]
  ds_meta <- discovery_metadata[discovery_metadata$dataset == ds, ]
  
  if (length(unique(ds_meta$group)) < 2) {
    log_message(sprintf("WARNING: %s has only one group, skipping.", ds), "WARNING")
    next
  }
  
  ds_design <- model.matrix(~ group, data = ds_meta)
  ds_fit <- lmFit(ds_expr, ds_design)
  ds_fit <- eBayes(ds_fit)
  ds_de <- topTable(ds_fit, coef = "groupTumor", number = Inf, sort.by = "P")
  ds_de$Gene <- rownames(ds_de)
  per_dataset_deg[[ds]] <- ds_de
  
  write.csv(ds_de, file.path(OUTPUT$de, sprintf("DEG_%s.csv", ds)))
}

# --- 14. FOLD CHANGE CONSISTENCY ---
if (length(per_dataset_deg) >= 2) {
  fc_consistency <- data.frame(
    Gene = common_genes,
    GSE42568_FC = per_dataset_deg$GSE42568$logFC[match(common_genes, per_dataset_deg$GSE42568$Gene)],
    GSE15852_FC = per_dataset_deg$GSE15852$logFC[match(common_genes, per_dataset_deg$GSE15852$Gene)],
    GSE29431_FC = per_dataset_deg$GSE29431$logFC[match(common_genes, per_dataset_deg$GSE29431$Gene)],
    Discovery_FC = de_results$logFC[match(common_genes, rownames(de_results))],
    stringsAsFactors = FALSE
  )
  
  fc_long <- fc_consistency %>%
    pivot_longer(cols = c(GSE42568_FC, GSE15852_FC, GSE29431_FC),
                 names_to = "Dataset", values_to = "FC") %>%
    mutate(Dataset = gsub("_FC", "", Dataset))
  
  p_fc_consistency <- ggplot(fc_long, aes(x = FC, color = Dataset)) +
    geom_density(alpha = 0.5, linewidth = 1) +
    scale_color_manual(values = CB_PALETTE[1:3]) +
    labs(title = "Fold-Change Distribution Across Datasets",
         x = "log2 Fold Change", y = "Density") +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_tiff_600dpi(p_fc_consistency, file.path(OUTPUT$de, "FC_consistency_by_dataset.tiff"))
  
  # Direction consistency
  for (ds in names(per_dataset_deg)) {
    fc_consistency[[paste0(ds, "_Direction")]] <- sign(fc_consistency[[paste0(ds, "_FC")]])
  }
  
  write.csv(fc_consistency, file.path(OUTPUT$de, "fold_change_consistency.csv"), row.names = FALSE)
}

log_message("Batch correction and differential expression analysis complete.")
