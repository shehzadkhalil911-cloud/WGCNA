# =============================================================================
# Script 05: Regulated Cell Death Activity Scores
# Breast Cancer 7 Cell Death Crosstalk
# =============================================================================

source(file.path("..", "config", "analysis_config.R"))

library(tidyverse)
library(GSVA)
library(GSEABase)
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
library(patchwork)

log_message("Starting activity score analysis.")

# --- 1. LOAD DATA ---
gene_lists <- readRDS(file.path(OUTPUT$inventory, "gene_lists_raw.rds"))
combat_matrix <- readRDS(file.path(OUTPUT$batch_correction, "combat_corrected_matrix.rds"))
discovery_metadata <- readRDS(file.path(OUTPUT$batch_correction, "discovery_metadata.rds"))

ferr_genes <- gene_lists$ferroptosis
necr_genes <- gene_lists$necroptosis
crg_genes <- gene_lists$cuproptosis

# --- 2. BUILD GENE SETS FOR GSVA ---
gs_list <- list(
  Ferroptosis = intersect(ferr_genes, rownames(combat_matrix)),
  Necroptosis = intersect(necr_genes, rownames(combat_matrix)),
  Cuproptosis = intersect(crg_genes, rownames(combat_matrix))
)

log_message(sprintf("Genes in expression matrix - Ferroptosis: %d, Necroptosis: %d, Cuproptosis: %d",
                    length(gs_list$Ferroptosis), length(gs_list$Necroptosis), length(gs_list$Cuproptosis)))

# --- 3. CALCULATE SSGSEA SCORES ---
gsva_par <- ssgseaParam(combat_matrix, gs_list)
scores <- gsva(gsva_par, verbose = FALSE)

log_message(sprintf("Activity scores calculated: %d programs x %d samples", nrow(scores), ncol(scores)))

saveRDS(scores, file.path(OUTPUT$activity, "program_scores.rds"))
write.csv(t(scores), file.path(OUTPUT$activity, "program_scores.csv"))

# --- 4. SCORE DISTRIBUTIONS ---
score_df <- as.data.frame(t(scores))
score_df$Sample <- rownames(score_df)
score_df$Group <- discovery_metadata$group[match(score_df$Sample, discovery_metadata$merged_id)]
score_df$Dataset <- discovery_metadata$dataset[match(score_df$Sample, discovery_metadata$merged_id)]

score_long <- score_df %>%
  pivot_longer(cols = c(Ferroptosis, Necroptosis, Cuproptosis),
               names_to = "Program", values_to = "Score")

# Boxplot: Tumour vs Normal for each program
p_boxes <- ggplot(score_long, aes(x = Program, y = Score, fill = Group)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7, width = 0.6) +
  scale_fill_manual(values = c("Normal" = "#0072B2", "Tumor" = "#D55E00")) +
  stat_compare_means(aes(label = ..p.signif..), method = "wilcox.test", 
                     label.x = 1.5, label.y = max(score_long$Score) * 1.05, size = 4) +
  labs(title = "Regulated Cell Death Activity Scores",
       subtitle = "Tumour vs Normal", y = "ssGSEA Score") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
save_tiff_600dpi(p_boxes, file.path(OUTPUT$activity, "score_boxplots.tiff"))

# Density plot
p_density <- ggplot(score_long, aes(x = Score, color = Group, linetype = Program)) +
  geom_density(linewidth = 0.8) +
  scale_color_manual(values = c("Normal" = "#0072B2", "Tumor" = "#D55E00")) +
  labs(title = "Activity Score Densities", x = "ssGSEA Score", y = "Density") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_tiff_600dpi(p_density, file.path(OUTPUT$activity, "score_density.tiff"))

# --- 5. PAIRED ANALYSIS ---
# GSE15852 has paired samples: even-indexed are Normal, odd-indexed are Tumor
gse15852_meta <- discovery_metadata[discovery_metadata$dataset == "GSE15852", ]
gse15852_normals <- gse15852_meta$merged_id[gse15852_meta$group == "Normal"]
gse15852_tumors <- gse15852_meta$merged_id[gse15852_meta$group == "Tumor"]

# Try to match pairs by extracting GSM IDs
normal_gsm <- gsub(".*_(GSM\\d+)", "\\1", gse15852_normals)
tumor_gsm <- gsub(".*_(GSM\\d+)", "\\1", gse15852_tumors)
normal_num <- as.numeric(gsub("GSM", "", normal_gsm))
tumor_num <- as.numeric(gsub("GSM", "", tumor_gsm))

# Find pairs where normal is even and tumor is normal+1
paired_normal_idx <- which(normal_num %% 2 == 0 & normal_num + 1 %in% tumor_num)
paired_tumor_idx <- match(normal_num[paired_normal_idx] + 1, tumor_num)

if (length(paired_normal_idx) > 0 && length(paired_tumor_idx) > 0 &&
    length(paired_normal_idx) == length(paired_tumor_idx)) {
  paired_df <- data.frame(
    Normal_Score = scores["Ferroptosis", gse15852_normals[paired_normal_idx]],
    Tumor_Score = scores["Ferroptosis", gse15852_tumors[paired_tumor_idx]],
    Program = "Ferroptosis"
  )
  
  if (nrow(paired_df) > 3 && var(paired_df$Normal_Score) > 0 && var(paired_df$Tumor_Score) > 0) {
    p_paired <- ggplot(paired_df, aes(x = Normal_Score, y = Tumor_Score)) +
      geom_point(size = 3, alpha = 0.7, color = "#0072B2") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      stat_cor(method = "spearman") +
      labs(title = "Paired Ferroptosis Scores (GSE15852)",
           x = "Normal Score", y = "Tumor Score") +
      theme_classic(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    save_tiff_600dpi(p_paired, file.path(OUTPUT$activity, "paired_ferroptosis.tiff"))
  }
}

# --- 6. SCORE CORRELATION ---
cor_scores <- cor(scores, method = "spearman")
write.csv(cor_scores, file.path(OUTPUT$activity, "score_correlation.csv"))

p_cor_heat <- Heatmap(cor_scores,
                      name = "Spearman r",
                      col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
                      cell_fun = function(j, i, x, y, width, height, fill) {
                        grid.text(sprintf("%.3f", cor_scores[i, j]), x, y, 
                                  gp = gpar(fontsize = 14, fontface = "bold"))
                      },
                      rect_gp = gpar(col = "grey60", lwd = 1.5),
                      column_title = "Score Correlation",
                      row_title = "Program")

open_tiff(file.path(OUTPUT$activity, "score_correlation_heatmap.tiff"), width = 5, height = 4)
draw(p_cor_heat)
dev.off()

# --- 7. PAIRWISE SCATTER PLOTS ---
score_df_long <- score_df
pairs_list <- list(
  c("Ferroptosis", "Necroptosis"),
  c("Ferroptosis", "Cuproptosis"),
  c("Necroptosis", "Cuproptosis")
)

for (prog_pair in pairs_list) {
  p <- ggplot(score_df_long, aes_string(x = prog_pair[1], y = prog_pair[2], color = "Group")) +
    geom_point(size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
    scale_color_manual(values = c("Normal" = "#0072B2", "Tumor" = "#D55E00")) +
    stat_cor(aes(label = ..r.label..), method = "spearman") +
    labs(title = sprintf("%s vs %s", prog_pair[1], prog_pair[2])) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_tiff_600dpi(p, file.path(OUTPUT$activity, sprintf("scatter_%s_vs_%s.tiff", prog_pair[1], prog_pair[2])))
}

# --- 8. SAMPLE-BY-SCORE HEATMAP ---
ht_score <- Heatmap(scores,
                    name = "Score",
                    col = colorRamp2(c(-0.3, 0, 0.3), c("#2166AC", "white", "#B2182B")),
                    top_annotation = HeatmapAnnotation(
                      Group = discovery_metadata$group,
                      Dataset = discovery_metadata$dataset,
                      col = list(
                        Group = c(Normal = "#0072B2", Tumor = "#D55E00"),
                        Dataset = c(GSE42568 = unname(CB_PALETTE[1]), GSE15852 = unname(CB_PALETTE[2]), GSE29431 = unname(CB_PALETTE[3]))
                      )
                    ),
                    show_column_names = FALSE,
                    column_title = "Sample-Level Activity Scores",
                    row_names_gp = gpar(fontsize = 12, fontface = "bold"))

open_tiff(file.path(OUTPUT$activity, "sample_score_heatmap.tiff"), width = 10, height = 4)
draw(ht_score)
dev.off()

# --- 9. SCORE BY DATASET ---
p_score_by_ds <- ggplot(score_long, aes(x = Dataset, y = Score, fill = Group)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_wrap(~ Program, scales = "free_y") +
  scale_fill_manual(values = c("Normal" = "#0072B2", "Tumor" = "#D55E00")) +
  labs(title = "Activity Scores by Dataset", y = "ssGSEA Score") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))
save_tiff_600dpi(p_score_by_ds, file.path(OUTPUT$activity, "score_by_dataset.tiff"))

# --- 10. DIFFERENTIAL ACTIVITY ANALYSIS ---
# Test score differences between tumour and normal
diff_results <- data.frame()
for (prog in c("Ferroptosis", "Necroptosis", "Cuproptosis")) {
  tumor_scores <- score_df[[prog]][score_df$Group == "Tumor"]
  normal_scores <- score_df[[prog]][score_df$Group == "Normal"]
  
  wt <- wilcox.test(tumor_scores, normal_scores)
  cohen_d <- (mean(tumor_scores) - mean(normal_scores)) / 
    sqrt(((length(tumor_scores)-1)*var(tumor_scores) + (length(normal_scores)-1)*var(normal_scores)) / 
           (length(tumor_scores) + length(normal_scores) - 2))
  
  diff_results <- rbind(diff_results, data.frame(
    Program = prog,
    Tumor_Mean = mean(tumor_scores),
    Normal_Mean = mean(normal_scores),
    Difference = mean(tumor_scores) - mean(normal_scores),
    P_value = wt$p.value,
    Cohen_d = cohen_d,
    stringsAsFactors = FALSE
  ))
}

diff_results$FDR <- p.adjust(diff_results$P_value, method = "BH")
write.csv(diff_results, file.path(OUTPUT$activity, "differential_activity_results.csv"), row.names = FALSE)
log_message("Differential activity results:")
print(diff_results)

log_message("Activity score analysis complete.")
