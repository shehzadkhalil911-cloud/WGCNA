source("R/utils.R")
load_libraries()
set.seed(123)

BASE_DIR <- "F:/independent projects/cancer/breast cancer"
OUT_DIR <- file.path(BASE_DIR, "Breast cancer 6 Ferroptosis")

log_msg("=== Step 7: Ferroptosis Activity Scoring ===")

ferr_genes <- readRDS(file.path(OUT_DIR, "data_processed/ferroptosis_gene_list/ferroptosis_genes_final.rds"))
combat_expr <- readRDS(file.path(OUT_DIR, "data_processed/discovery/discovery_combat.rds"))
metadata <- read.csv(file.path(OUT_DIR, "data_processed/discovery/discovery_metadata.csv"))

# GSVA for ferroptosis scoring
library(GSVA)

ferr_present <- intersect(ferr_genes, rownames(combat_expr))
log_msg(sprintf("Ferroptosis genes present in expression: %d", length(ferr_present)))

gsva_params <- gsvaParam(as.matrix(combat_expr), list(ferr_present),
                          kcdf = "Gaussian", minSize = 5, maxSize = 500)
ferr_score <- gsva(gsva_params)
ferr_score_vec <- as.numeric(ferr_score)
names(ferr_score_vec) <- colnames(combat_expr)

# Add to metadata
metadata$ferroptosis_score <- ferr_score_vec[metadata$sample_id]

# Compare tumor vs normal
score_comparison <- wilcox.test(ferroptosis_score ~ group, data = metadata)
log_msg(sprintf("Ferroptosis score: Tumor vs Normal p = %.6f", score_comparison$p.value))

# Score box plot
p_score_box <- ggplot(metadata, aes(group, ferroptosis_score, fill = group)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  scale_fill_manual(values = c(Tumor = "#E41A1C", Normal = "#377EB8")) +
  labs(title = "Ferroptosis Activity Score",
       subtitle = sprintf("Wilcoxon p = %.4e", score_comparison$p.value),
       x = "", y = "Ferroptosis Score") +
  theme_minimal() + theme(legend.position = "none")
save_tiff_600dpi(p_score_box, file.path(OUT_DIR, "results/figures/Figure4A_ferroptosis_score_box.tiff"))

# Score density plot
p_score_density <- ggplot(metadata, aes(ferroptosis_score, fill = group, color = group)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c(Tumor = "#E41A1C", Normal = "#377EB8")) +
  scale_color_manual(values = c(Tumor = "#E41A1C", Normal = "#377EB8")) +
  labs(title = "Ferroptosis Score Density", x = "Ferroptosis Score", y = "Density") +
  theme_minimal()
save_tiff_600dpi(p_score_density, file.path(OUT_DIR, "results/figures/Figure4B_score_density.tiff"))

# Score by dataset
p_score_dataset <- ggplot(metadata, aes(dataset, ferroptosis_score, fill = group)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = c(Tumor = "#E41A1C", Normal = "#377EB8")) +
  labs(title = "Ferroptosis Score by Dataset", x = "Dataset", y = "Ferroptosis Score") +
  theme_minimal()
save_tiff_600dpi(p_score_dataset, file.path(OUT_DIR, "results/figures/Figure4C_score_by_dataset.tiff"))

# Score correlation with DE-FRG expression
de_frg <- readRDS(file.path(OUT_DIR, "results/ferroptosis_genes/DE_FRG_results.rds"))
if (nrow(de_frg) > 0) {
  de_frg_present <- intersect(rownames(de_frg), rownames(combat_expr))
  if (length(de_frg_present) > 0) {
    score_gene_cor <- apply(combat_expr[de_frg_present, ], 1, function(x) {
      cor(x, ferr_score_vec, method = "spearman")
    })
    score_gene_p <- apply(combat_expr[de_frg_present, ], 1, function(x) {
      cor.test(x, ferr_score_vec, method = "spearman")$p.value
    })

    score_gene_df <- data.frame(
      Gene = names(score_gene_cor),
      Correlation = score_gene_cor,
      P_value = score_gene_p,
      Adj_P = p.adjust(score_gene_p, method = "BH"),
      stringsAsFactors = FALSE
    )
    write.csv(score_gene_df, file.path(OUT_DIR, "results/ferroptosis_score/score_gene_correlation.csv"), row.names = FALSE)

    # Score-gene correlation heatmap
    score_gene_df <- score_gene_df[order(abs(score_gene_df$Correlation), decreasing = TRUE), ]
    if (nrow(score_gene_df) > 1) {
      open_tiff(file.path(OUT_DIR, "results/figures/Figure4D_score_gene_cor_heatmap.tiff"), width = 7, height = 9)
      pheatmap::pheatmap(as.matrix(score_gene_df$Correlation),
                         labels_row = score_gene_df$Gene,
                         cluster_cols = FALSE,
                         color = colorRampPalette(c("blue", "white", "red"))(100),
                         main = "DE-FRG Correlation with Ferroptosis Score",
                         fontsize = 7)
      dev.off()
    }
  }
}

# Score heatmap (samples)
sample_order <- order(ferr_score_vec)
score_heatmap_data <- data.frame(
  Score = ferr_score_vec[sample_order],
  Group = metadata$group[match(names(ferr_score_vec)[sample_order], metadata$sample_id)],
  Dataset = metadata$dataset[match(names(ferr_score_vec)[sample_order], metadata$sample_id)],
  Sample = names(ferr_score_vec)[sample_order]
)

p_score_heatmap <- ggplot(score_heatmap_data, aes(Sample, 1, fill = Score)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white") +
  labs(title = "Ferroptosis Score Across Samples", x = "Samples", y = "") +
  theme_minimal() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
save_tiff_600dpi(p_score_heatmap, file.path(OUT_DIR, "results/figures/Figure4E_score_heatmap.tiff"),
                 width = 10, height = 3)

# Paired score plot for GSE15852
discovery_list <- readRDS(file.path(OUT_DIR, "data_processed/discovery/discovery_list.rds"))
if (discovery_list[["GSE15852"]]$is_paired) {
  pair_info <- discovery_list[["GSE15852"]]$pair_info
  if (!is.null(pair_info)) {
    pair_info$ferroptosis_score <- ferr_score_vec[pair_info$sample_id]
    pair_plot_data <- pair_info[!is.na(pair_info$ferroptosis_score), ]
    p_paired <- ggplot(pair_plot_data, aes(class, ferroptosis_score, group = pair_id, color = class)) +
      geom_line(alpha = 0.3) + geom_point(size = 2) +
      scale_color_manual(values = c(Tumor = "#E41A1C", Normal = "#377EB8")) +
      labs(title = "Paired Ferroptosis Score (GSE15852)", x = "", y = "Ferroptosis Score") +
      theme_minimal()
    save_tiff_600dpi(p_paired, file.path(OUT_DIR, "results/figures/Figure4F_paired_score.tiff"))
  }
}

# Save score data
metadata_out <- metadata[, c("sample_id", "dataset", "group", "ferroptosis_score")]
write.csv(metadata_out, file.path(OUT_DIR, "results/ferroptosis_score/ferroptosis_scores.csv"), row.names = FALSE)
saveRDS(ferr_score_vec, file.path(OUT_DIR, "results/ferroptosis_score/ferroptosis_score_vector.rds"))

log_msg("Step 7 complete: Ferroptosis activity scoring finished")
