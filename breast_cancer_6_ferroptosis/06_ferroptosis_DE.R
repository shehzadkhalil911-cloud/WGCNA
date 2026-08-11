source("R/utils.R")
load_libraries()
set.seed(123)

BASE_DIR <- "F:/independent projects/cancer/breast cancer"
OUT_DIR <- file.path(BASE_DIR, "Breast cancer 6 Ferroptosis")

log_msg("=== Step 6: Ferroptosis-Related Differential Analysis ===")

ferr_genes <- readRDS(file.path(OUT_DIR, "data_processed/ferroptosis_gene_list/ferroptosis_genes_final.rds"))
deg_table <- readRDS(file.path(OUT_DIR, "results/differential_expression/combined_DEG_results.rds"))
combat_expr <- readRDS(file.path(OUT_DIR, "data_processed/discovery/discovery_combat.rds"))
metadata <- read.csv(file.path(OUT_DIR, "data_processed/discovery/discovery_metadata.csv"))

# Intersect ferroptosis genes with DEG results
ferr_detected <- intersect(ferr_genes, rownames(deg_table))
ferr_undetected <- setdiff(ferr_genes, rownames(deg_table))
ferr_in_deg <- intersect(ferr_genes, rownames(deg_table)[deg_table$Significant])

# DE-FRGs
de_frg <- deg_table[ferr_in_deg, ]
de_frg_up <- de_frg[de_frg$logFC > 0, ]
de_frg_down <- de_frg[de_frg$logFC < 0, ]
ferr_non_sig <- setdiff(ferr_detected, ferr_in_deg)

log_msg(sprintf("Ferroptosis genes detected: %d/%d", length(ferr_detected), length(ferr_genes)))
log_msg(sprintf("Significant DE-FRGs: %d (%d up, %d down)", nrow(de_frg), nrow(de_frg_up), nrow(de_frg_down)))
log_msg(sprintf("Undetected ferroptosis genes: %d", length(ferr_undetected)))

# Save tables
write.csv(data.frame(Gene = ferr_detected), file.path(OUT_DIR, "results/ferroptosis_genes/detected_ferroptosis_genes.csv"), row.names = FALSE)
write.csv(data.frame(Gene = ferr_undetected), file.path(OUT_DIR, "results/ferroptosis_genes/undetected_ferroptosis_genes.csv"), row.names = FALSE)
write.csv(de_frg, file.path(OUT_DIR, "results/ferroptosis_genes/DE_FRGs.csv"))
write.csv(de_frg_up, file.path(OUT_DIR, "results/ferroptosis_genes/upregulated_DE_FRGs.csv"))
write.csv(de_frg_down, file.path(OUT_DIR, "results/ferroptosis_genes/downregulated_DE_FRGs.csv"))
write.csv(data.frame(Gene = ferr_non_sig), file.path(OUT_DIR, "results/ferroptosis_genes/non_significant_ferroptosis_genes.csv"), row.names = FALSE)

# Venn diagram: DEGs vs Ferroptosis genes
library(VennDiagram)
venn_list <- list(
  DEGs = rownames(deg_table)[deg_table$Significant],
  Ferroptosis = ferr_genes
)
open_tiff(file.path(OUT_DIR, "results/figures/Figure3F_DEG_ferroptosis_venn.tiff"))
venn.plot <- venn.diagram(venn_list, filename = NULL,
                           fill = c("red", "steelblue"), alpha = 0.5,
                           cat.cex = 1.2, cex = 1.5,
                           main = "DEGs vs Ferroptosis Genes")
grid.draw(venn.plot)
dev.off()

# DE-FRG heatmap
if (nrow(de_frg) > 0) {
  de_frg_expr <- combat_expr[rownames(de_frg), , drop = FALSE]
  de_frg_scaled <- t(scale(t(de_frg_expr)))
  anno_col <- data.frame(Group = metadata$group, row.names = colnames(combat_expr))

  open_tiff(file.path(OUT_DIR, "results/figures/Figure3G_DE_FRG_heatmap.tiff"), width = 9, height = 7)
  pheatmap::pheatmap(de_frg_scaled, annotation_col = anno_col,
                     show_rownames = TRUE, show_colnames = FALSE,
                     fontsize_row = 7, cutree_cols = 2,
                     main = "DE-FRG Expression Heatmap")
  dev.off()
}

# DE-FRG volcano plot
deg_table$is_ferroptosis <- rownames(deg_table) %in% ferr_genes
deg_table$is_de_frg <- rownames(deg_table) %in% ferr_in_deg
deg_table$label <- ifelse(deg_table$is_de_frg, rownames(deg_table), "")

p_volc_ferr <- ggplot(deg_table, aes(logFC, -log10(P.Value),
                                      color = is_de_frg)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70")) +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20, color = "red") +
  labs(title = "Ferroptosis Genes in DEG Volcano Plot",
       x = expression(log[2]~(Fold~Change)),
       y = expression(-log[10]~(P-value))) +
  theme_minimal() + theme(legend.position = "none")
save_tiff_600dpi(p_volc_ferr, file.path(OUT_DIR, "results/figures/Figure3H_ferroptosis_volcano.tiff"))

# FC bar plot for DE-FRGs
if (nrow(de_frg) > 0) {
  de_frg_ordered <- de_frg[order(de_frg$logFC), ]
  de_frg_ordered$Gene <- factor(rownames(de_frg_ordered), levels = rownames(de_frg_ordered))

  p_fc_bar <- ggplot(de_frg_ordered, aes(Gene, logFC, fill = Direction)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c(Up = "red", Down = "blue")) +
    labs(title = "DE-FRG Fold Changes", x = "", y = expression(log[2]~(FC))) +
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
  save_tiff_600dpi(p_fc_bar, file.path(OUT_DIR, "results/figures/Figure3I_DE_FRG_FC_bar.tiff"),
                   width = 8, height = 5)
}

# Effect direction heatmap across datasets
discovery_list <- readRDS(file.path(OUT_DIR, "data_processed/discovery/discovery_list.rds"))

fc_across <- NULL
for (ds in names(discovery_list)) {
  deg_path <- file.path(OUT_DIR, "results/differential_expression", paste0(ds, "_DEG_table.csv"))
  if (file.exists(deg_path)) {
    deg_ds <- read.csv(deg_path, row.names = 1)
    fc_ds <- deg_ds[intersect(rownames(de_frg), rownames(deg_ds)), "logFC", drop = FALSE]
    colnames(fc_ds) <- ds
    if (is.null(fc_across)) fc_across <- fc_ds else fc_across[[ds]] <- fc_ds[[ds]]
  }
}

if (!is.null(fc_across) && nrow(fc_across) > 1) {
  fc_across <- fc_across[complete.cases(fc_across), ]
  direction_mat <- sign(fc_across)

  open_tiff(file.path(OUT_DIR, "results/figures/Figure3J_DE_FRG_direction_heatmap.tiff"), width = 7, height = 8)
  pheatmap::pheatmap(direction_mat, display_numbers = TRUE,
                     color = c("blue", "white", "red"),
                     breaks = c(-1.5, -0.5, 0.5, 1.5),
                     main = "DE-FRG Direction Across Datasets",
                     fontsize = 8)
  dev.off()
}

saveRDS(de_frg, file.path(OUT_DIR, "results/ferroptosis_genes/DE_FRG_results.rds"))
log_msg("Step 6 complete: Ferroptosis-related differential analysis finished")
