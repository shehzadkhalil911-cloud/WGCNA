source("R/utils.R")
load_libraries()
set.seed(123)

BASE_DIR <- "F:/independent projects/cancer/breast cancer"
OUT_DIR <- file.path(BASE_DIR, "Breast cancer 6 Ferroptosis")

log_msg("=== Step 13: Immune Infiltration Analysis ===")

combat_expr <- readRDS(file.path(OUT_DIR, "data_processed/discovery/discovery_combat.rds"))
metadata <- read.csv(file.path(OUT_DIR, "data_processed/discovery/discovery_metadata.csv"))
ferr_score <- readRDS(file.path(OUT_DIR, "results/ferroptosis_score/ferroptosis_score_vector.rds"))
de_frg <- readRDS(file.path(OUT_DIR, "results/ferroptosis_genes/DE_FRG_results.rds"))
module_info <- read.csv(file.path(OUT_DIR, "results/WGCNA/module_assignments.csv"), row.names = 1)
MEs <- readRDS(file.path(OUT_DIR, "results/WGCNA/module_eigengenes.rds"))
top_hubs <- readRDS(file.path(OUT_DIR, "results/WGCNA/top_hub_genes.rds"))

# Use ssGSEA with immune gene signatures
immune_signatures <- list(
  B_cells = c("CD19", "MS4A1", "CD79A", "CD79B", "PAX5", "BLK", "BANK1"),
  CD4_T_cells = c("CD4", "CD3D", "CD3E", "CD3G", "IL7R", "LCK", "ZAP70"),
  CD8_T_cells = c("CD8A", "CD8B", "CD3D", "CD3E", "GZMK", "GZMA", "PRF1"),
  NK_cells = c("KLRD1", "NKG7", "KIR2DL1", "KIR3DL1", "NCR1", "FCGR3A"),
  Macrophages = c("CD68", "CD163", "CD14", "FCGR2A", "ITGAM", "CSF1R", "MRC1"),
  Monocytes = c("CD14", "FCGR3A", "CSF1R", "LYZ", "S100A8", "S100A9"),
  Neutrophils = c("FCGR3B", "CEACAM3", "S100A12", "CXCR1", "CXCR2", "CSF3R"),
  Dendritic_cells = c("ITGAX", "FCER1A", "CLEC4C", "CD1C", "THBD", "NRP1"),
  Mast_cells = c("KIT", "TPSAB1", "TPSB2", "CPA3", "FCER1A", "MS4A2"),
  Myeloid_Derived_Suppressor = c("CD33", "CD11B", "ARG1", "NOS2", "S100A8", "S100A9", "IL10"),
  Regulatory_T = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "IKZF2"),
  Stromal = c("COL1A1", "COL1A2", "FAP", "PDPN", "ACTA2", "TAGLN", "VIM"),
  Angiogenesis = c("PECAM1", "VWF", "CDH5", "KDR", "FLT1", "ANGPT1", "TEK")
)

library(GSVA)
common_immune_genes <- intersect(unique(unlist(immune_signatures)), rownames(combat_expr))
log_msg(sprintf("Common immune genes: %d/%d", length(common_immune_genes), length(unique(unlist(immune_signatures)))))

# Filter signatures to include only present genes
immune_sigs <- lapply(immune_signatures, function(sig) intersect(sig, rownames(combat_expr)))
immune_sigs <- immune_sigs[sapply(immune_sigs, length) >= 3]
log_msg(sprintf("Immune signatures with >= 3 genes: %d", length(immune_sigs)))

if (length(immune_sigs) >= 2) {
  gsva_params <- gsvaParam(as.matrix(combat_expr), immune_sigs, kcdf = "Gaussian", minSize = 3, maxSize = 500)
  immune_scores <- gsva(gsva_params)
  log_msg(sprintf("Immune scores computed: %d cell types x %d samples", nrow(immune_scores), ncol(immune_scores)))

  rownames(immune_scores) <- gsub("\\.", "_", rownames(immune_scores))
  immune_df <- as.data.frame(t(immune_scores))
  immune_df$Group <- metadata$group[match(rownames(immune_df), metadata$sample_id)]
  immune_df$Dataset <- metadata$dataset[match(rownames(immune_df), metadata$sample_id)]

  # Save immune scores
  write.csv(immune_df, file.path(OUT_DIR, "results/immune_infiltration/immune_scores.csv"))
  saveRDS(immune_scores, file.path(OUT_DIR, "results/immune_infiltration/immune_scores.rds"))

  # Immune cell box plots
  imm_long <- reshape2::melt(immune_df, id.vars = c("Group", "Dataset"),
                               variable.name = "Cell_Type", value.name = "Score")

  p_immune_box <- ggplot(imm_long, aes(Cell_Type, Score, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    scale_fill_manual(values = c(Tumor = "#E41A1C", Normal = "#377EB8")) +
    labs(title = "Immune Cell Infiltration Scores", x = "", y = "ssGSEA Score") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
  save_tiff_600dpi(p_immune_box, file.path(OUT_DIR, "results/figures/Figure11A_immune_boxplots.tiff"),
                   width = 10, height = 6)

  # Immune infiltration heatmap (Tumor vs Normal comparison)
  imm_means <- aggregate(. ~ Group, data = immune_df[, -ncol(immune_df)], FUN = mean)
  rownames(imm_means) <- imm_means$Group
  imm_means$Group <- NULL
  imm_means_scaled <- t(scale(t(imm_means)))

  open_tiff(file.path(OUT_DIR, "results/figures/Figure11B_immune_heatmap.tiff"), width = 6, height = 5)
  pheatmap::pheatmap(imm_means_scaled, display_numbers = TRUE,
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     main = "Immune Infiltration (Tumor vs Normal)",
                     fontsize = 8)
  dev.off()

  # Immune cell correlation heatmap
  imm_cor <- cor(immune_scores, method = "spearman")
  open_tiff(file.path(OUT_DIR, "results/figures/Figure11C_immune_correlation_heatmap.tiff"), width = 8, height = 7)
  pheatmap::pheatmap(imm_cor, display_numbers = TRUE, fontsize = 7,
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     main = "Immune Cell Correlation")
  dev.off()

  # Differential immune infiltration (Wilcoxon)
  immune_diff <- data.frame(Cell_Type = colnames(immune_df)[1:(ncol(immune_df)-2)])
  immune_diff$Wilcoxon_P <- sapply(immune_diff$Cell_Type, function(ct) {
    wilcox.test(immune_df[[ct]] ~ immune_df$Group)$p.value
  })
  immune_diff$Adj_P <- p.adjust(immune_diff$Wilcoxon_P, method = "BH")
  immune_diff$Mean_Tumor <- sapply(immune_diff$Cell_Type, function(ct) {
    mean(immune_df[[ct]][immune_df$Group == "Tumor"], na.rm = TRUE)
  })
  immune_diff$Mean_Normal <- sapply(immune_diff$Cell_Type, function(ct) {
    mean(immune_df[[ct]][immune_df$Group == "Normal"], na.rm = TRUE)
  })
  immune_diff$Effect_Size <- immune_diff$Mean_Tumor - immune_diff$Mean_Normal
  write.csv(immune_diff, file.path(OUT_DIR, "results/immune_infiltration/differential_immune.csv"), row.names = FALSE)

  # --- Ferroptosis-Immune Crosstalk ---
  # Ferroptosis score vs Immune cells
  imm_for_cor <- immune_df[, !colnames(immune_df) %in% c("Group", "Dataset")]
  ferr_immune_cor <- sapply(imm_for_cor, function(x) {
    cor(x, ferr_score[rownames(imm_for_cor)], method = "spearman", use = "pairwise")
  })
  ferr_immune_p <- sapply(imm_for_cor, function(x) {
    cor.test(x, ferr_score[rownames(imm_for_cor)], method = "spearman")$p.value
  })

  ferr_immune_df <- data.frame(
    Cell_Type = names(ferr_immune_cor),
    Correlation = ferr_immune_cor,
    P_value = ferr_immune_p,
    Adj_P = p.adjust(ferr_immune_p, method = "BH"),
    stringsAsFactors = FALSE
  )
  write.csv(ferr_immune_df, file.path(OUT_DIR, "results/immune_infiltration/ferroptosis_immune_correlation.csv"), row.names = FALSE)

  # Ferroptosis score-immune heatmap
  ferr_imm_mat <- as.matrix(ferr_immune_df$Correlation)
  rownames(ferr_imm_mat) <- ferr_immune_df$Cell_Type
  colnames(ferr_imm_mat) <- "Ferroptosis Score"

  open_tiff(file.path(OUT_DIR, "results/figures/Figure11D_ferroptosis_immune_heatmap.tiff"), width = 5, height = 6)
  pheatmap::pheatmap(ferr_imm_mat, display_numbers = TRUE, fontsize = 8,
                     color = colorRampPalette(c("blue", "white", "red"))(100),
                     main = "Ferroptosis Score vs Immune Cells",
                     cluster_rows = FALSE, cluster_cols = FALSE)
  dev.off()

  # DE-FRG vs Immune correlation
  de_frg_present <- intersect(rownames(de_frg), rownames(combat_expr))
  if (length(de_frg_present) > 5) {
    gene_immune_cor <- cor(t(combat_expr[de_frg_present, ]), t(immune_scores),
                           method = "spearman")
    write.csv(gene_immune_cor, file.path(OUT_DIR, "results/immune_infiltration/gene_immune_correlation.csv"))

    open_tiff(file.path(OUT_DIR, "results/figures/Figure11E_gene_immune_heatmap.tiff"), width = 9, height = 8)
    pheatmap::pheatmap(gene_immune_cor, display_numbers = FALSE, fontsize = 6,
                       color = colorRampPalette(c("blue", "white", "red"))(100),
                       main = "DE-FRG vs Immune Cell Correlation",
                       )
    dev.off()
  }

  # Module-Immune relationship
  common_samples <- intersect(rownames(MEs), colnames(immune_scores))
  if (length(common_samples) > 10 && ncol(MEs) > 2) {
    me_immune_cor <- cor(MEs[common_samples, ], t(immune_scores[, common_samples]),
                         method = "spearman")

    open_tiff(file.path(OUT_DIR, "results/figures/Figure11F_module_immune_heatmap.tiff"), width = 9, height = 7)
    pheatmap::pheatmap(me_immune_cor, display_numbers = TRUE, fontsize = 7,
                       color = colorRampPalette(c("blue", "white", "red"))(100),
                       main = "Module-Immune Relationships",
                       )
    dev.off()
    write.csv(me_immune_cor, file.path(OUT_DIR, "results/immune_infiltration/module_immune_correlations.csv"))
  }

  # Hub gene-immune network
  hub_genes_present <- intersect(top_hubs$Gene[1:min(20, nrow(top_hubs))], rownames(combat_expr))
  if (length(hub_genes_present) > 3) {
    hub_immune_cor <- cor(t(combat_expr[hub_genes_present, ]), t(immune_scores), method = "spearman")
    hub_x <- t(combat_expr[hub_genes_present, ])
    hub_y <- t(immune_scores)
    hub_immune_pval <- apply(hub_x, 2, function(x) {
      apply(hub_y, 2, function(y) {
        cor.test(x, y, method = "spearman")$p.value
      })
    })

    # Significant hub-immune pairs
    sig_pairs <- which(hub_immune_pval < 0.05, arr.ind = TRUE)
    if (nrow(sig_pairs) > 0) {
      sig_network <- data.frame(
        Gene = colnames(hub_immune_pval)[sig_pairs[, 2]],
        Immune_Cell = rownames(hub_immune_pval)[sig_pairs[, 1]],
        Correlation = hub_immune_cor[cbind(sig_pairs[, 2], sig_pairs[, 1])],
        P_Value = hub_immune_pval[sig_pairs],
        stringsAsFactors = FALSE
      )
      write.csv(sig_network, file.path(OUT_DIR, "results/immune_infiltration/significant_gene_immune_network.csv"), row.names = FALSE)
    }
  }

  saveRDS(immune_scores, file.path(OUT_DIR, "results/immune_infiltration/immune_scores.rds"))
} else {
  log_msg("WARNING: Too few immune signatures with sufficient genes")
}

log_msg("Step 13 complete: Immune infiltration analysis finished")
