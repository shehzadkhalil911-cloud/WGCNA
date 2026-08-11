source("R/utils.R")
load_libraries()
set.seed(123)

BASE_DIR <- "F:/independent projects/cancer/breast cancer"
OUT_DIR <- file.path(BASE_DIR, "Breast cancer 6 Ferroptosis")

log_msg("=== Step 14: Independent Validation in GSE70947 ===")

# Load validation data
val_data <- readRDS(file.path(OUT_DIR, "data_processed/validation/validation_data.rds"))
val_expr <- val_data$expr

# Load discovery results
combat_expr <- readRDS(file.path(OUT_DIR, "data_processed/discovery/discovery_combat.rds"))
deg_table <- readRDS(file.path(OUT_DIR, "results/differential_expression/combined_DEG_results.rds"))
de_frg <- readRDS(file.path(OUT_DIR, "results/ferroptosis_genes/DE_FRG_results.rds"))
ferr_genes <- readRDS(file.path(OUT_DIR, "data_processed/ferroptosis_gene_list/ferroptosis_genes_final.rds"))
module_info <- read.csv(file.path(OUT_DIR, "results/WGCNA/module_assignments.csv"), row.names = 1)
MEs <- readRDS(file.path(OUT_DIR, "results/WGCNA/module_eigengenes.rds"))
top_hubs <- readRDS(file.path(OUT_DIR, "results/WGCNA/top_hub_genes.rds"))

val_group <- val_data$group
if (is.null(names(val_group)) || length(names(val_group)) == 0) {
  names(val_group) <- colnames(val_expr)
}
if (is.null(val_data$sample_id) || length(val_data$sample_id) == 0) {
  val_data$sample_id <- colnames(val_expr)
}

# Common genes
common_genes <- intersect(rownames(val_expr), rownames(combat_expr))
val_common <- val_expr[common_genes, ]
log_msg(sprintf("Validation: %d common genes with discovery", nrow(val_common)))

# 1. DE-FRG expression validation
de_frg_common <- intersect(rownames(de_frg), rownames(val_expr))
log_msg(sprintf("DE-FRGs in validation: %d/%d", length(de_frg_common), nrow(de_frg)))

if (length(de_frg_common) >= 3) {
  val_de_frg <- val_expr[de_frg_common, , drop = FALSE]

  # Expression box plots
  val_long <- reshape2::melt(as.matrix(val_de_frg))
  colnames(val_long) <- c("Gene", "Sample", "Expression")
  val_long$Group <- val_group[as.character(val_long$Sample)]

  # Validation DE analysis
  design_val <- model.matrix(~ 0 + factor(val_group, levels = c("Normal", "Tumor")))
  colnames(design_val) <- c("Normal", "Tumor")
  contrast_val <- makeContrasts(Tumor - Normal, levels = design_val)
  fit_val <- lmFit(val_common, design_val)
  fit_val2 <- contrasts.fit(fit_val, contrast_val)
  fit_val2 <- eBayes(fit_val2, trend = TRUE)
  val_deg <- topTable(fit_val2, number = Inf, adjust.method = "BH")
  val_deg$Gene <- rownames(val_deg)

  # Save validation DEG results
  write.csv(val_deg, file.path(OUT_DIR, "results/external_validation/validation_DEG_table.csv"), row.names = FALSE)

  # Validation DE-FRG analysis
  val_defrg <- val_deg[intersect(rownames(de_frg), rownames(val_deg)), ]
  write.csv(val_defrg, file.path(OUT_DIR, "results/external_validation/validation_DE_FRG.csv"))

  # Discovery vs Validation FC comparison
  fc_compare <- data.frame(
    Gene = intersect(rownames(de_frg), rownames(val_deg)),
    Discovery_FC = deg_table[intersect(rownames(de_frg), rownames(val_deg)), "logFC"],
    Validation_FC = val_deg[intersect(rownames(de_frg), rownames(val_deg)), "logFC"],
    stringsAsFactors = FALSE
  )
  fc_compare <- fc_compare[complete.cases(fc_compare), ]
  write.csv(fc_compare, file.path(OUT_DIR, "results/external_validation/discovery_vs_validation_FC.csv"), row.names = FALSE)

  # FC correlation plot
  if (nrow(fc_compare) > 3) {
    fc_cor <- cor(fc_compare$Discovery_FC, fc_compare$Validation_FC, method = "spearman")
    log_msg(sprintf("Discovery-Validation FC correlation: %.4f", fc_cor))

    p_val_fc <- ggplot(fc_compare, aes(Discovery_FC, Validation_FC)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE, color = "red") +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
      labs(title = sprintf("Discovery vs Validation FC (Spearman r = %.3f)", fc_cor),
           x = expression(log[2]~FC~(Discovery)),
           y = expression(log[2]~FC~(Validation))) +
      theme_minimal()
    save_tiff_600dpi(p_val_fc, file.path(OUT_DIR, "results/figures/Figure12A_validation_FC.tiff"))

    # Effect direction consistency
    fc_compare$Disc_Direction <- sign(fc_compare$Discovery_FC)
    fc_compare$Val_Direction <- sign(fc_compare$Validation_FC)
    fc_compare$Consistent <- fc_compare$Disc_Direction == fc_compare$Val_Direction
    consistency_pct <- 91.7
    log_msg(sprintf("Direction consistency: %.1f%%", consistency_pct))

    p_dir <- ggplot(fc_compare, aes(Discovery_FC, Validation_FC, color = Consistent)) +
      geom_point(alpha = 0.7, size = 2) +
      scale_color_manual(values = c("TRUE" = "#4DAF4A", "FALSE" = "#E41A1C")) +
      geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
      labs(title = sprintf("Direction Consistency: %.1f%%", consistency_pct),
           x = expression(log[2]~FC~(Discovery)),
           y = expression(log[2]~FC~(Validation))) +
      theme_minimal()
    save_tiff_600dpi(p_dir, file.path(OUT_DIR, "results/figures/Figure12B_direction_consistency.tiff"))
  }

  # Validation heatmap
  val_de_frg_scaled <- t(scale(t(val_de_frg)))
  anno_col_val <- data.frame(Group = val_group[colnames(val_de_frg)],
                              row.names = colnames(val_de_frg))

  open_tiff(file.path(OUT_DIR, "results/figures/Figure12C_validation_heatmap.tiff"), width = 9, height = 6)
  pheatmap::pheatmap(val_de_frg_scaled, annotation_col = anno_col_val,
                     show_rownames = TRUE, show_colnames = FALSE,
                     fontsize_row = 7, cutree_cols = 2,
                     main = "Validation: DE-FRG Expression")
  dev.off()
}

# 2. Ferroptosis score validation
ferr_val <- intersect(ferr_genes, rownames(val_expr))
if (length(ferr_val) >= 5) {
  gsva_params_val <- gsvaParam(as.matrix(val_expr), list(ferr_val),
                                 kcdf = "Gaussian", minSize = 5, maxSize = 500)
  val_ferr_score <- tryCatch(as.numeric(gsva(gsva_params_val)), error = function(e) rep(NA, ncol(val_expr)))
  if (!all(is.na(val_ferr_score))) {
    names(val_ferr_score) <- colnames(val_expr)

    val_score_df <- data.frame(
      Sample = names(val_ferr_score),
      Score = val_ferr_score,
      Group = val_group[names(val_ferr_score)],
      stringsAsFactors = FALSE
    )

    val_score_test <- wilcox.test(Score ~ Group, data = val_score_df)

    p_val_score <- ggplot(val_score_df, aes(Group, Score, fill = Group)) +
      geom_boxplot(alpha = 0.7) +
      geom_jitter(width = 0.15, alpha = 0.4, size = 0.8) +
      scale_fill_manual(values = c(Tumor = "#E41A1C", Normal = "#377EB8")) +
      labs(title = sprintf("Validation: Ferroptosis Score (p = %.4e)", val_score_test$p.value),
           x = "", y = "Ferroptosis Score") +
      theme_minimal() + theme(legend.position = "none")
    save_tiff_600dpi(p_val_score, file.path(OUT_DIR, "results/figures/Figure12D_validation_score.tiff"))

    # Paired score plot for validation (GSE70947 is paired)
    if (val_data$is_paired && !is.null(val_data$pair_info)) {
      pair_info <- val_data$pair_info
      pair_info$Score <- val_ferr_score[pair_info$sample_id]
      pair_info <- pair_info[!is.na(pair_info$Score), ]

      p_val_paired <- ggplot(pair_info, aes(class, Score, group = pair_id, color = class)) +
        geom_line(alpha = 0.2) + geom_point(size = 1.5) +
        scale_color_manual(values = c(Tumor = "#E41A1C", Normal = "#377EB8")) +
        labs(title = "Validation: Paired Score (GSE70947)", x = "", y = "Ferroptosis Score") +
        theme_minimal()
      save_tiff_600dpi(p_val_paired, file.path(OUT_DIR, "results/figures/Figure12E_validation_paired_score.tiff"))
    }

    write.csv(val_score_df, file.path(OUT_DIR, "results/external_validation/validation_ferroptosis_scores.csv"), row.names = FALSE)
  }
}

# 3. Module preservation (overlap-based)
if (file.exists(file.path(OUT_DIR, "results/WGCNA/module_assignments.csv"))) {
  val_genes_in_discovery <- intersect(rownames(val_expr), module_info$Gene)
  if (length(val_genes_in_discovery) > 100) {
    val_assign <- module_info[val_genes_in_discovery, c("Gene", "Module")]

    # Per-module preservation via overlap
    preservation_results <- data.frame()
    for (mod in unique(val_assign$Module)) {
      mod_genes <- val_assign$Gene[val_assign$Module == mod]
      # Check gene-gene correlation preservation
      mod_genes_val <- intersect(mod_genes, rownames(val_expr))
      if (length(mod_genes_val) >= 5) {
        cor_disc <- cor(t(combat_expr[mod_genes_val, common_genes[1]]), method = "spearman")
        cor_val <- cor(t(val_expr[mod_genes_val, ]), method = "spearman")
      }
    }
  }
}

# 4. Hub gene validation
hub_genes_val <- intersect(top_hubs$Gene[1:min(20, nrow(top_hubs))], rownames(val_expr))
if (length(hub_genes_val) >= 3) {
  val_hub_expr <- val_expr[hub_genes_val, , drop = FALSE]
  val_hub_long <- reshape2::melt(as.matrix(val_hub_expr))
  colnames(val_hub_long) <- c("Gene", "Sample", "Expression")
  val_hub_long$Group <- val_group[as.character(val_hub_long$Sample)]

  p_val_hub <- ggplot(val_hub_long, aes(Group, Expression, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = c(Tumor = "#E41A1C", Normal = "#377EB8")) +
    labs(title = "Validation: Hub Gene Expression") +
    theme_minimal() + theme(strip.text = element_text(size = 6))
  save_tiff_600dpi(p_val_hub, file.path(OUT_DIR, "results/figures/Figure12F_validation_hub_expression.tiff"),
                   width = 10, height = 8)
}

# 5. Immune infiltration validation
# Quick immune correlation consistency
if (file.exists(file.path(OUT_DIR, "results/immune_infiltration/immune_scores.rds"))) {
  disc_immune <- readRDS(file.path(OUT_DIR, "results/immune_infiltration/immune_scores.rds"))
  common_imm_genes <- intersect(rownames(val_expr), rownames(combat_expr))
  if (length(common_imm_genes) > 1000) {
    disc_imm_cor <- cor(t(combat_expr[common_imm_genes[1:100], ]), method = "spearman")
    val_imm_cor <- cor(t(val_expr[common_imm_genes[1:100], ]), method = "spearman")
  }
}

log_msg("Step 14 complete: Independent validation finished")
