# =============================================================================
# Script 10: Multi-Evidence Prioritisation and Independent Validation
# Breast Cancer 7 Cell Death Crosstalk
# =============================================================================

source(file.path("..", "config", "analysis_config.R"))

library(tidyverse)
library(ggrepel)
library(ggpubr)
library(GSVA)
library(ComplexHeatmap)
library(circlize)

log_message("Starting multi-evidence prioritisation.")

# --- 0. LOAD GENE LISTS ---
gene_lists <- readRDS(file.path(OUTPUT$inventory, "gene_lists_raw.rds"))
ferr_genes <- gene_lists$ferroptosis
necr_genes <- gene_lists$necroptosis
crg_genes <- gene_lists$cuproptosis

# --- 1. LOAD DATA ---
shared_info <- readRDS(file.path(OUTPUT$shared, "shared_gene_classification.rds"))
de_results <- readRDS(file.path(OUTPUT$de, "discovery_DEG_full.rds"))
deg_sig <- read.csv(file.path(OUTPUT$de, "discovery_DEG_significant.csv"), row.names = 1)
module_df <- read.csv(file.path(OUTPUT$wgcna, "module_assignments.csv"))
module_trait <- readRDS(file.path(OUTPUT$wgcna, "module_trait_correlations.rds"))
intra_conn <- read.csv(file.path(OUTPUT$wgcna, "intramodular_connectivity.csv"))
scores <- readRDS(file.path(OUTPUT$activity, "program_scores.rds"))
combat_matrix <- readRDS(file.path(OUTPUT$batch_correction, "combat_corrected_matrix.rds"))
discovery_metadata <- readRDS(file.path(OUTPUT$batch_correction, "discovery_metadata.rds"))

# --- 2. FOCUS ON ALL REGULATED CELL DEATH GENES IN DEG SPACE ---
all_de_rcd <- unique(c(
  shared_info$Gene[shared_info$DE_FRG],
  shared_info$Gene[shared_info$DE_NRG],
  shared_info$Gene[shared_info$DE_CRG]
))

log_message(sprintf("Building evidence matrix for %d candidate genes", length(all_de_rcd)))

# --- 3. BUILD EVIDENCE MATRIX ---
evidence <- data.frame(
  Gene = all_de_rcd,
  stringsAsFactors = FALSE
)

# 3a. Differential expression
evidence$logFC <- de_results[all_de_rcd, "logFC"]
evidence$adj.P.Val <- de_results[all_de_rcd, "adj.P.Val"]
evidence$abs_logFC <- abs(evidence$logFC)
evidence$DE_Significant <- ifelse(evidence$adj.P.Val < ADJ_PVAL_THRESHOLD & 
                                    evidence$abs_logFC >= LOG2FC_THRESHOLD, 1, 0)

# 3b. Program membership
evidence$Ferroptosis <- ifelse(all_de_rcd %in% shared_info$Gene[shared_info$DE_FRG], 1, 0)
evidence$Necroptosis <- ifelse(all_de_rcd %in% shared_info$Gene[shared_info$DE_NRG], 1, 0)
evidence$Cuproptosis <- ifelse(all_de_rcd %in% shared_info$Gene[shared_info$DE_CRG], 1, 0)
evidence$N_Programs <- evidence$Ferroptosis + evidence$Necroptosis + evidence$Cuproptosis

# 3c. Shared gene status
evidence$Shared <- ifelse(all_de_rcd %in% shared_info$Gene[
  shared_info$Category %in% c("All Three", "Ferroptosis + Necroptosis",
                               "Ferroptosis + Cuproptosis", "Necroptosis + Cuproptosis")], 1, 0)

# 3d. WGCNA module membership
gene_mod <- module_df$Module[match(all_de_rcd, module_df$Gene)]
gene_mod_color <- module_df$ModuleColor[match(all_de_rcd, module_df$Gene)]
evidence$Module <- gene_mod

# Gene significance (correlation with group)
gs <- as.numeric(cor(t(combat_matrix[all_de_rcd, ]), 
                     as.numeric(discovery_metadata$group) - 1, 
                     method = "spearman"))
evidence$Gene_Significance <- abs(gs)

# Intramodular connectivity
ic <- intra_conn$kWithin[match(all_de_rcd, intra_conn$Gene)]
evidence$Intramodular_Connectivity <- ifelse(is.na(ic), 0, ic)

# 3e. Activity score correlation
for (prog in c("Ferroptosis", "Necroptosis", "Cuproptosis")) {
  prog_score <- scores[prog, ]
  gene_cor <- apply(combat_matrix[intersect(all_de_rcd, rownames(combat_matrix)), ], 1, 
                    function(x) cor(x, prog_score, method = "spearman"))
  evidence[[paste0("Cor_", prog)]] <- gene_cor[match(all_de_rcd, names(gene_cor))]
}

# 3f. Module-trait association
module_color_map <- setNames(module_df$ModuleColor, module_df$Gene)
evidence$Module_Tumor_Association <- NA
module_cor_mat <- module_trait$cor
available_modules <- intersect(unique(module_df$ModuleColor), rownames(module_cor_mat))
for (gene in all_de_rcd) {
  mod_color <- module_color_map[gene]
  if (!is.na(mod_color) && mod_color %in% available_modules) {
    evidence$Module_Tumor_Association[evidence$Gene == gene] <- module_cor_mat[mod_color, "Group"]
  }
}

# --- 4. STANDARDISED EVIDENCE SCORE ---
# Standardise each numeric evidence dimension to [0, 1]
standardise_01 <- function(x) {
  if (length(unique(x)) == 1) return(rep(0.5, length(x)))
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}

evidence$Score_DE <- standardise_01(-log10(evidence$adj.P.Val + 1e-300)) * 0.5 + 
  standardise_01(evidence$abs_logFC) * 0.5

evidence$Score_Program <- evidence$N_Programs / 3

evidence$Score_WGCNA <- standardise_01(evidence$Gene_Significance) * 0.5 +
  standardise_01(evidence$Intramodular_Connectivity) * 0.5

evidence$Score_Activity <- rowMeans(
  cbind(standardise_01(abs(evidence$Cor_Ferroptosis)),
        standardise_01(abs(evidence$Cor_Necroptosis)),
        standardise_01(abs(evidence$Cor_Cuproptosis))),
  na.rm = TRUE
)

# Combined evidence score
evidence$Combined_Score <- rowMeans(
  cbind(evidence$Score_DE, 
        evidence$Score_Program,
        evidence$Score_WGCNA,
        evidence$Score_Activity),
  na.rm = TRUE
)

# Handle NAs
evidence$Combined_Score[is.na(evidence$Combined_Score)] <- 0

# Rank
evidence <- evidence[order(evidence$Combined_Score, decreasing = TRUE), ]
evidence$Rank <- 1:nrow(evidence)

write.csv(evidence, file.path(OUTPUT$prioritisation, "candidate_evidence_matrix.csv"), row.names = FALSE)
saveRDS(evidence, file.path(OUTPUT$prioritisation, "candidate_evidence_matrix.rds"))

log_message(sprintf("Evidence matrix: %d candidates, top 3: %s", 
                    nrow(evidence), paste(head(evidence$Gene, 3), collapse = ", ")))

# --- 5. EVIDENCE HEATMAP ---
evidence_cols <- c("Score_DE", "Score_Program", "Score_WGCNA", "Score_Activity", "Combined_Score")
# Top 50 candidates
top_n <- min(50, nrow(evidence))
ev_hm <- as.matrix(evidence[1:top_n, evidence_cols])
rownames(ev_hm) <- evidence$Gene[1:top_n]

ht_ev <- Heatmap(ev_hm,
                 name = "Score",
                 col = colorRamp2(c(0, 0.5, 1), c("white", "orange", "red3")),
                 rect_gp = gpar(col = "grey60", lwd = 0.5),
                 row_names_gp = gpar(fontsize = 8),
                 column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                 column_title = "Multi-Evidence Scores",
                 row_title = "Candidate Genes",
                 cluster_rows = FALSE,
                 cluster_columns = FALSE)

open_tiff(file.path(OUTPUT$prioritisation, "evidence_heatmap.tiff"), width = 8, height = 10)
draw(ht_ev)
dev.off()

# --- 6. BUBBLE PLOT ---
top_evidence <- evidence[1:top_n, ]
p_bubble <- ggplot(top_evidence, aes(x = Score_DE, y = Score_WGCNA, 
                                      size = Combined_Score, color = Score_Program)) +
  geom_point(alpha = 0.7) +
  geom_text_repel(aes(label = Gene), size = 2.5, max.overlaps = 15) +
  scale_size_continuous(range = c(2, 10)) +
  scale_color_gradient(low = "blue", high = "red") +
  labs(title = "Evidence Dimensions for Top Candidates",
       x = "DE Score", y = "WGCNA Score",
       size = "Combined", color = "Program Score") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_tiff_600dpi(p_bubble, file.path(OUTPUT$prioritisation, "evidence_bubble.tiff"))

# --- 7. SAVE TOP CANDIDATE LIST ---
write.csv(evidence[1:min(30, nrow(evidence)), ], 
          file.path(OUTPUT$prioritisation, "top_30_candidates.csv"), row.names = FALSE)

# --- 8. INDEPENDENT VALIDATION IN GSE70947 ---
log_message("Starting independent validation in GSE70947...")

# Load GSE70947
validation_expr <- read.csv(EXPRESSION_FILES$GSE70947, row.names = 1, check.names = FALSE)
validation_classes <- read.csv(SAMPLE_CLASS_FILES$GSE70947, stringsAsFactors = FALSE)

# Match samples
common_valid <- intersect(colnames(validation_expr), validation_classes$sample_id)
validation_expr <- as.matrix(validation_expr[, common_valid])
validation_classes <- validation_classes[match(common_valid, validation_classes$sample_id), ]

group_valid <- factor(validation_classes$class, levels = c("Normal", "Tumor"))

log_message(sprintf("Validation cohort: %d genes x %d samples", nrow(validation_expr), ncol(validation_expr)))

# 8a. Validate expression of candidate genes
candidate_genes <- evidence$Gene[evidence$Combined_Score > 0.3]
candidate_genes <- intersect(candidate_genes, rownames(validation_expr))

log_message(sprintf("Validating %d candidate genes in GSE70947", length(candidate_genes)))

if (length(candidate_genes) > 3) {
  # Per-gene validation
  validation_de <- data.frame(Gene = candidate_genes, stringsAsFactors = FALSE)
  for (gene in candidate_genes) {
    tumor_expr <- validation_expr[gene, group_valid == "Tumor"]
    normal_expr <- validation_expr[gene, group_valid == "Normal"]
    if (length(tumor_expr) > 2 && length(normal_expr) > 2) {
      fc <- mean(tumor_expr) - mean(normal_expr)
      wt <- wilcox.test(tumor_expr, normal_expr)
      validation_de$FC[validation_de$Gene == gene] <- fc
      validation_de$P_value[validation_de$Gene == gene] <- wt$p.value
    }
  }
  
  validation_de$FDR <- p.adjust(validation_de$P_value, method = "BH")
  validation_de$Discovery_FC <- de_results[candidate_genes, "logFC"]
  validation_de$FC_Consistent <- sign(validation_de$FC) == sign(validation_de$Discovery_FC)
  
  write.csv(validation_de, file.path(OUTPUT$validation, "candidate_validation.csv"), row.names = FALSE)
  
  # Discovery vs validation FC plot
  p_fc_valid <- ggplot(validation_de, aes(x = Discovery_FC, y = FC, color = FC_Consistent)) +
    geom_point(size = 2.5, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
    scale_color_manual(values = c("grey60", "#009E73")) +
    labs(title = "Discovery vs Validation Fold Change",
         x = "Discovery log2FC", y = "Validation log2FC") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_tiff_600dpi(p_fc_valid, file.path(OUTPUT$validation, "FC_discovery_vs_validation.tiff"))
  
  # Expression heatmap in validation
  if (length(candidate_genes) > 2) {
    valid_mat <- validation_expr[candidate_genes, ]
    # Remove rows with zero variance
    row_var <- apply(valid_mat, 1, var, na.rm = TRUE)
    valid_mat <- valid_mat[row_var > 0, , drop = FALSE]
    if (nrow(valid_mat) > 2) {
      valid_hm <- t(scale(t(valid_mat)))
      
      col_ann <- HeatmapAnnotation(
        Group = group_valid,
        col = list(Group = c(Normal = "#0072B2", Tumor = "#D55E00"))
      )
      
      ht_valid <- tryCatch({
        Heatmap(valid_hm,
                            name = "Z-score",
                            col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
                            top_annotation = col_ann,
                            show_row_names = TRUE,
                            show_column_names = FALSE,
                            row_names_gp = gpar(fontsize = 7),
                            column_title = "Validation: Candidate Gene Expression",
                            cluster_rows = FALSE,
                            cluster_columns = FALSE,
                            show_row_dend = FALSE,
                            show_column_dend = FALSE,
                            na_col = "grey90")
      }, error = function(e) {
        log_message(sprintf("Validation heatmap error: %s", e$message), "WARNING")
        NULL
      })
      if (!is.null(ht_valid)) {
        open_tiff(file.path(OUTPUT$validation, "validation_expression_heatmap.tiff"), width = 10, height = 8)
        draw(ht_valid)
        dev.off()
      }
    }
  }
}

# 8b. Validate activity scores
# Map to common gene symbols
common_gsva_genes <- intersect(rownames(validation_expr), 
                                unique(c(ferr_genes, necr_genes, crg_genes)))
valid_gs_list <- list(
  Ferroptosis = intersect(ferr_genes, rownames(validation_expr)),
  Necroptosis = intersect(necr_genes, rownames(validation_expr)),
  Cuproptosis = intersect(crg_genes, rownames(validation_expr))
)

valid_gs_list <- valid_gs_list[sapply(valid_gs_list, length) >= 5]

if (length(valid_gs_list) >= 2) {
  gsva_par_valid <- ssgseaParam(validation_expr, valid_gs_list, minSize = 5, maxSize = 500)
  valid_scores <- gsva(gsva_par_valid)
  
  saveRDS(valid_scores, file.path(OUTPUT$validation, "validation_scores.rds"))
  
  # Score boxplots
  valid_score_df <- as.data.frame(t(valid_scores))
  valid_score_df$Group <- group_valid
  
  valid_long <- valid_score_df %>%
    pivot_longer(cols = -Group, names_to = "Program", values_to = "Score")
  
  p_valid_scores <- ggplot(valid_long, aes(x = Program, y = Score, fill = Group)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
    scale_fill_manual(values = c("Normal" = "#0072B2", "Tumor" = "#D55E00")) +
    stat_compare_means(aes(label = ..p.signif..), method = "wilcox.test",
                       label.x = 1.5, size = 4) +
    labs(title = "Validation: Activity Scores", y = "ssGSEA Score") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_tiff_600dpi(p_valid_scores, file.path(OUTPUT$validation, "validation_scores_boxplot.tiff"))
  
  # Score correlation validation
  valid_cor <- tryCatch(cor(t(valid_scores), method = "spearman", use = "complete.obs"),
                        error = function(e) NULL)
  
  if (!is.null(valid_cor) && nrow(valid_cor) > 1) {
    open_tiff(file.path(OUTPUT$validation, "validation_score_correlation.tiff"), width = 5, height = 4)
    ht_vc <- Heatmap(valid_cor,
                     name = "Spearman r",
                     col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(sprintf("%.3f", valid_cor[i, j]), x, y, gp = gpar(fontsize = 12))
                     },
                     rect_gp = gpar(col = "grey60", lwd = 1.5),
                     column_title = "Validation Score Correlation",
                     cluster_rows = FALSE,
                     cluster_columns = FALSE,
                     show_row_dend = FALSE,
                     show_column_dend = FALSE)
    draw(ht_vc)
    dev.off()
  }
}

log_message("Multi-evidence prioritisation and validation complete.")
