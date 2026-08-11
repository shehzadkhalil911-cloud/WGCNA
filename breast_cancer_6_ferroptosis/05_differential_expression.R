source("R/utils.R")
load_libraries()
set.seed(123)

BASE_DIR <- "F:/independent projects/cancer/breast cancer"
OUT_DIR <- file.path(BASE_DIR, "Breast cancer 6 Ferroptosis")

log_msg("=== Step 5: Differential Expression Analysis ===")

combat_expr <- readRDS(file.path(OUT_DIR, "data_processed/discovery/discovery_combat.rds"))
metadata <- read.csv(file.path(OUT_DIR, "data_processed/discovery/discovery_metadata.csv"))
ferr_genes <- readRDS(file.path(OUT_DIR, "data_processed/ferroptosis_gene_list/ferroptosis_genes_final.rds"))

# --- Combined limma analysis ---
design <- model.matrix(~ 0 + group, data = metadata)
colnames(design) <- gsub("group", "", colnames(design))
contrast <- makeContrasts(Tumor - Normal, levels = design)

fit <- lmFit(combat_expr, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2, trend = TRUE)
deg_table <- topTable(fit2, number = Inf, adjust.method = "BH")

log_msg(sprintf("Total DEGs tested: %d", nrow(deg_table)))
log_msg(sprintf("Significant DEGs (adj.P<0.05, |logFC|>=1): %d",
                sum(deg_table$adj.P.Val < 0.05 & abs(deg_table$logFC) >= 1)))

# Add direction
deg_table$Direction <- ifelse(deg_table$logFC > 0, "Up", "Down")
deg_table$Significant <- deg_table$adj.P.Val < 0.05 & abs(deg_table$logFC) >= 1

sig_deg <- deg_table[deg_table$adj.P.Val < 0.05 & abs(deg_table$logFC) >= 1, ]
up_deg <- sig_deg[sig_deg$logFC > 0, ]
down_deg <- sig_deg[sig_deg$logFC < 0, ]

# Save DEG results
write.csv(deg_table, file.path(OUT_DIR, "results/differential_expression/complete_DEG_table.csv"))
write.csv(deg_table[deg_table$Significant, ],
          file.path(OUT_DIR, "results/differential_expression/significant_DEGs.csv"), row.names = FALSE)
write.csv(data.frame(Gene = rownames(sig_deg), t = sig_deg$t,
                     stringsAsFactors = FALSE),
          file.path(OUT_DIR, "results/differential_expression/ranked_gene_list.csv"), row.names = FALSE)

# Volcano plot
deg_table$Gene <- rownames(deg_table)
deg_table$label <- ifelse(deg_table$Significant & abs(deg_table$logFC) > 2, deg_table$Gene, "")

p_volc <- ggplot(deg_table, aes(logFC, -log10(P.Value), color = Significant)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50")) +
  geom_text_repel(aes(label = label), size = 2, max.overlaps = 15) +
  labs(title = "Volcano Plot: Tumor vs Normal",
       x = expression(log[2]~(Fold~Change)),
       y = expression(-log[10]~(P-value))) +
  theme_minimal() + theme(legend.position = "none")
save_tiff_600dpi(p_volc, file.path(OUT_DIR, "results/figures/Figure3A_volcano.tiff"))

# MA plot
deg_table$AveExpr <- deg_table$AveExpr
p_ma <- ggplot(deg_table, aes(AveExpr, logFC, color = Significant)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50")) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "blue", alpha = 0.5) +
  labs(title = "MA Plot", x = "Average Expression", y = expression(log[2]~(Fold~Change))) +
  theme_minimal()
save_tiff_600dpi(p_ma, file.path(OUT_DIR, "results/figures/Figure3B_MA_plot.tiff"))

# Top DEG heatmaps
for (n in c(25, 50, 100)) {
  top_genes <- rownames(sig_deg)[1:min(n, nrow(sig_deg))]
  if (length(top_genes) > 0) {
    mat <- combat_expr[top_genes, , drop = FALSE]
    mat_scaled <- t(scale(t(mat)))
    anno_col <- data.frame(Group = metadata$group, Dataset = metadata$dataset,
                           row.names = colnames(combat_expr))
    anno_colors <- list(Group = c(Tumor = "black", Normal = "white"),
                        Dataset = c(GSE42568 = "#E41A1C", GSE15852 = "#377EB8", GSE29431 = "#4DAF4A"))

    open_tiff(file.path(OUT_DIR, sprintf("results/figures/Figure3C_top%d_DEG_heatmap.tiff", n)),
              width = 9, height = 7)
    pheatmap::pheatmap(mat_scaled, annotation_col = anno_col,
                       annotation_colors = anno_colors,
                       show_rownames = TRUE, show_colnames = FALSE,
                       fontsize_row = 5, cutree_cols = 2, scale = "none",
                       main = paste("Top", n, "DEGs"))
    dev.off()
  }
}

# Up/Down count plot
count_data <- data.frame(
  Direction = c("Upregulated", "Downregulated"),
  Count = c(nrow(up_deg), nrow(down_deg))
)
p_counts <- ggplot(count_data, aes(Direction, Count, fill = Direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  labs(title = "Differentially Expressed Genes", x = "", y = "Count") +
  theme_minimal() + theme(legend.position = "none")
save_tiff_600dpi(p_counts, file.path(OUT_DIR, "results/figures/Figure3D_DEG_counts.tiff"))

# --- Per-dataset differential expression ---
discovery_list <- readRDS(file.path(OUT_DIR, "data_processed/discovery/discovery_list.rds"))

per_dataset_deg <- list()
fc_matrix <- NULL
for (ds in names(discovery_list)) {
  d <- discovery_list[[ds]]
  expr <- d$expr
  # Only use common genes
  common <- intersect(rownames(expr), common <- rownames(combat_expr))
  expr <- expr[common, ]
  grp <- factor(d$group, levels = c("Normal", "Tumor"))

  design_ds <- model.matrix(~ grp)
  fit_ds <- lmFit(expr, design_ds)
  fit_ds <- eBayes(fit_ds, trend = TRUE)
  deg_ds <- topTable(fit_ds, coef = 2, number = Inf, adjust.method = "BH")
  deg_ds$Direction <- ifelse(deg_ds$logFC > 0, "Up", "Down")
  per_dataset_deg[[ds]] <- deg_ds

  write.csv(deg_ds, file.path(OUT_DIR, sprintf("results/differential_expression/%s_DEG_table.csv", ds)))

  if (is.null(fc_matrix)) {
    fc_matrix <- data.frame(Gene = rownames(deg_ds), stringsAsFactors = FALSE)
  }
  fc_df <- data.frame(Gene = rownames(deg_ds), fc = deg_ds$logFC, stringsAsFactors = FALSE)
  colnames(fc_df)[2] <- ds
  fc_matrix <- merge(fc_matrix, fc_df, by = "Gene", all = TRUE)
}

rownames(fc_matrix) <- fc_matrix$Gene
fc_matrix$Gene <- NULL
fc_matrix <- fc_matrix[complete.cases(fc_matrix), ]

# FC correlation across datasets
fc_cor <- cor(fc_matrix, method = "spearman", use = "pairwise")
write.csv(fc_cor, file.path(OUT_DIR, "results/differential_expression/fold_change_correlation.csv"))
write.csv(fc_matrix, file.path(OUT_DIR, "results/differential_expression/per_dataset_fold_changes.csv"))

# Direction consistency
direction_mat <- sign(fc_matrix)
direction_consistency <- apply(direction_mat, 1, function(x) all(x == x[1]))
log_msg(sprintf("Direction-consistent genes across datasets: %d/%d", sum(direction_consistency), nrow(fc_matrix)))

p_consistency <- ggplot(data.frame(Consistent = c("Consistent", "Not Consistent"),
                                   Count = c(sum(direction_consistency), sum(!direction_consistency))),
                        aes(Consistent, Count, fill = Consistent)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Consistent" = "#4DAF4A", "Not Consistent" = "#E41A1C")) +
  labs(title = "Fold-Change Direction Consistency Across Discovery Datasets") +
  theme_minimal() + theme(legend.position = "none")
save_tiff_600dpi(p_consistency, file.path(OUT_DIR, "results/figures/Figure3E_FC_consistency.tiff"))

saveRDS(deg_table, file.path(OUT_DIR, "results/differential_expression/combined_DEG_results.rds"))
saveRDS(per_dataset_deg, file.path(OUT_DIR, "results/differential_expression/per_dataset_DEG_results.rds"))
saveRDS(list(up = up_deg, down = down_deg, sig = sig_deg), file.path(OUT_DIR, "results/differential_expression/DEG_subsets.rds"))

log_msg(sprintf("DEG summary: %d up, %d down, %d total significant",
                nrow(up_deg), nrow(down_deg), nrow(sig_deg)))
log_msg("Step 5 complete: Differential expression analysis finished")
