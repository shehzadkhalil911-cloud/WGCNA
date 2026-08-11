###############################################################################
# SCRIPT 7: Independent Validation in GSE70947
# Multi-Cohort WGCNA: Breast Cancer Hub Gene Discovery
###############################################################################
library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(pROC)
library(dplyr)
library(tidyr)
library(gridExtra)
library(cowplot)

options(stringsAsFactors = FALSE)

OUT_DIR <- "F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)"
FIG_DIR <- file.path(OUT_DIR, "figures/validation")
TAB_DIR <- file.path(OUT_DIR, "tables")
RDA_DIR <- file.path(OUT_DIR, "RData")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

load(file.path(RDA_DIR, "01_processed_data_common.RData"))
load(file.path(RDA_DIR, "03_WGCNA_GSE42568.RData"))
wgcna_42568 <- wgcna_result; rm(wgcna_result)
load(file.path(RDA_DIR, "03_WGCNA_GSE15852.RData"))
wgcna_15852 <- wgcna_result; rm(wgcna_result)
load(file.path(RDA_DIR, "03_WGCNA_GSE29431.RData"))
wgcna_29431 <- wgcna_result; rm(wgcna_result)
load(file.path(RDA_DIR, "04_module_matching.RData"))
load(file.path(RDA_DIR, "05_hub_genes.RData"))

cat("========== Independent Validation in GSE70947 ==========\n")

# Expression validation dataset
val_expr <- comm_70947$expr
val_meta <- comm_70947$meta

cat(sprintf("Validation dataset: %d genes, %d samples (%d Normal, %d Tumor)\n", 
            nrow(val_expr), ncol(val_expr), sum(val_meta$class=="Normal"), sum(val_meta$class=="Tumor")))

# Find overlapping hub genes
hub_genes_to_validate <- intersect(final_hubs, rownames(val_expr))
cat(sprintf("Hub genes found in validation dataset: %d/%d\n", length(hub_genes_to_validate), length(final_hubs)))

###############################################################################
# A. Expression Validation
###############################################################################
cat("\n--- Expression Validation ---\n")

val_results <- data.frame()
for (gene in hub_genes_to_validate) {
  normal_expr <- val_expr[gene, val_meta$class == "Normal"]
  tumor_expr <- val_expr[gene, val_meta$class == "Tumor"]
  
  # Welch t-test
  tt <- t.test(tumor_expr, normal_expr)
  fc <- mean(tumor_expr) - mean(normal_expr)
  
  val_results <- rbind(val_results, data.frame(
    Gene = gene,
    Normal_mean = mean(normal_expr),
    Tumor_mean = mean(tumor_expr),
    logFC = fc,
    P_value = tt$p.value,
    stringsAsFactors = FALSE
  ))
}
val_results$FDR <- p.adjust(val_results$P_value, method = "BH")
val_results$Direction_validated <- ifelse(val_results$logFC > 0, "Up", "Down")
val_results$Significant <- val_results$FDR < 0.05

# Check direction consistency with discovery
direction_consistent <- sapply(hub_genes_to_validate, function(g) {
  disc_direction <- sign(mean(c(
    hub_summary$logFC_42568[hub_summary$Gene == g],
    hub_summary$logFC_15852[hub_summary$Gene == g],
    hub_summary$logFC_29431[hub_summary$Gene == g]
  ), na.rm = TRUE))
  val_direction <- sign(val_results$logFC[val_results$Gene == g])
  disc_direction == val_direction
})
val_results$Direction_consistent <- direction_consistent

cat(sprintf("Validated (FDR<0.05): %d/%d\n", sum(val_results$Significant), nrow(val_results)))
cat(sprintf("Direction consistent: %d/%d\n", sum(val_results$Direction_consistent), nrow(val_results)))

write.csv(val_results, file.path(TAB_DIR, "validation_expression_results.csv"), row.names = FALSE)

# Box plots for top validated genes
top_valid <- val_results %>% filter(Significant) %>% arrange(P_value) %>% head(20)
for (gene in top_valid$Gene) {
  plot_df <- data.frame(
    Expression = val_expr[gene, ],
    Class = val_meta$class
  )
  p <- ggplot(plot_df, aes(x = Class, y = Expression, fill = Class)) +
    geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
    scale_fill_manual(values = c("Normal" = "blue", "Tumor" = "red")) +
    labs(title = paste(gene, "- Validation in GSE70947"),
         subtitle = sprintf("logFC=%.2f, FDR=%.2e", 
                           val_results$logFC[val_results$Gene == gene],
                           val_results$FDR[val_results$Gene == gene])) +
    theme_bw() + theme(legend.position = "none")
  ggsave(file.path(FIG_DIR, paste0("boxplot_", gene, "_validation.tiff")), p, width = 6, height = 6, dpi = 600, compression = "lzw")
}

# Validation heatmap
if (nrow(top_valid) > 5) {
  heat_data <- val_expr[top_valid$Gene, ]
  rownames(heat_data) <- top_valid$Gene
  annotation_col <- data.frame(Class = val_meta$class, row.names = colnames(val_expr))
  annotation_colors <- list(Class = c(Normal = "blue", Tumor = "red"))
  
  heat_data_scaled <- t(scale(t(heat_data)))
  
  tiff(file.path(FIG_DIR, "validation_heatmap.tiff"), width = 10, height = 8, units = "in", res = 600, compression = "lzw")
  pheatmap::pheatmap(heat_data_scaled,
           annotation_col = annotation_col,
           annotation_colors = annotation_colors,
           show_rownames = TRUE, show_colnames = FALSE,
           main = "Hub Gene Expression Validation (GSE70947)")
  dev.off()
}

# Discovery vs Validation fold-change plot
fc_comparison <- data.frame(
  Gene = hub_genes_to_validate,
  Discovery_logFC = sapply(hub_genes_to_validate, function(g) {
    mean(c(hub_summary$logFC_42568[hub_summary$Gene == g],
           hub_summary$logFC_15852[hub_summary$Gene == g],
           hub_summary$logFC_29431[hub_summary$Gene == g]), na.rm = TRUE)
  }),
  Validation_logFC = val_results$logFC
)

p_fc <- ggplot(fc_comparison, aes(x = Discovery_logFC, y = Validation_logFC)) +
  geom_point(aes(color = val_results$Significant), alpha = 0.7, size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = TRUE, color = "blue", alpha = 0.3) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50")) +
  labs(title = "Discovery vs Validation Fold Change",
       x = "Mean Discovery logFC", y = "Validation logFC (GSE70947)") +
  theme_bw() + theme(legend.position = "none")

ggsave(file.path(FIG_DIR, "discovery_vs_validation_foldchange.tiff"), p_fc, width = 8, height = 7, dpi = 600, compression = "lzw")

###############################################################################
# B. Diagnostic Validation (ROC)
###############################################################################
cat("\n--- Diagnostic Validation ---\n")

# Individual gene ROC
roc_results <- data.frame()
for (gene in hub_genes_to_validate) {
  roc_obj <- roc(val_meta$class_binary, val_expr[gene, ], quiet = TRUE)
  ci <- ci.auc(roc_obj)
  roc_results <- rbind(roc_results, data.frame(
    Gene = gene,
    AUC = as.numeric(roc_obj$auc),
    CI_lower = ci[1],
    CI_upper = ci[3],
    stringsAsFactors = FALSE
  ))
}
roc_results <- roc_results[order(-roc_results$AUC), ]
write.csv(roc_results, file.path(TAB_DIR, "validation_roc_individual.csv"), row.names = FALSE)

cat(sprintf("Top 10 hub genes by AUC:\n"))
print(head(roc_results, 10))

# Plot top 5 ROC curves
top_roc <- head(roc_results, 5)
roc_list <- lapply(top_roc$Gene, function(g) {
  roc(val_meta$class_binary, val_expr[g, ], quiet = TRUE)
})
names(roc_list) <- top_roc$Gene

roc_plot <- ggroc(roc_list, legacy.axes = TRUE) + 
  labs(title = "Hub Gene ROC Curves (GSE70947)",
       x = "1 - Specificity", y = "Sensitivity") +
  scale_color_brewer(palette = "Set1") +
  theme_bw()
ggsave(file.path(FIG_DIR, "ROC_curves_top_hubs.tiff"), roc_plot, width = 8, height = 7, dpi = 600, compression = "lzw")

# Combined signature ROC (using mean expression of top hubs)
top_signature_genes <- head(roc_results$Gene, 10)
signature_score <- colMeans(val_expr[top_signature_genes, , drop = FALSE], na.rm = TRUE)
roc_signature <- roc(val_meta$class_binary, signature_score, quiet = TRUE)
ci_sig <- ci.auc(roc_signature)

cat(sprintf("Combined signature (top 10 hubs) AUC: %.3f (%.3f-%.3f)\n", 
            roc_signature$auc, ci_sig[1], ci_sig[3]))

# Plot combined ROC
tiff(file.path(FIG_DIR, "ROC_combined_signature.tiff"), width = 8, height = 7, units = "in", res = 600, compression = "lzw")
plot(roc_signature, main = paste("Combined Hub Gene Signature ROC (GSE70947)\nAUC =", 
                                 round(roc_signature$auc, 3)),
     col = "red", lwd = 2, legacy.axes = TRUE)
dev.off()

###############################################################################
# C. Module Eigengene Validation  
###############################################################################
cat("\n--- Module Eigengene Validation ---\n")

# Get consensus module gene sets
for (cm in consensus_modules) {
  cm_name <- paste(cm$GSE42568, cm$GSE15852, cm$GSE29431, sep="_")
  
  # Get genes in this module that exist in validation set
  mod1_genes <- names(wgcna_42568$moduleColors[wgcna_42568$moduleColors == cm$GSE42568])
  val_genes <- intersect(mod1_genes, rownames(val_expr))
  
  if (length(val_genes) >= 5) {
    # Calculate eigengene in validation set
    val_dat <- t(val_expr[val_genes, , drop = FALSE])
    val_ME <- svd(scale(val_dat), nu = 1, nv = 1)$u[,1] * sqrt(ncol(val_dat))
    
    # Test association with tumor status
    tt <- t.test(val_ME ~ val_meta$class)
    
    cat(sprintf("Module %s: ME validation P=%.2e, diff=%.3f (%d genes)\n", 
                cm_name, tt$p.value, 
                mean(val_ME[val_meta$class == "Tumor"]) - mean(val_ME[val_meta$class == "Normal"]),
                length(val_genes)))
    
    # Plot
    plot_df <- data.frame(ME = val_ME, Class = val_meta$class)
    p <- ggplot(plot_df, aes(x = Class, y = ME, fill = Class)) +
      geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.2, alpha = 0.4) +
      scale_fill_manual(values = c("Normal" = "blue", "Tumor" = "red")) +
      labs(title = paste("Module Eigengene Validation:", cm_name),
           subtitle = sprintf("P=%.2e", tt$p.value)) +
      theme_bw() + theme(legend.position = "none")
    ggsave(file.path(FIG_DIR, paste0("ME_validation_", cm_name, ".tiff")), p, width = 6, height = 6, dpi = 600, compression = "lzw")
  }
}

save(val_results, roc_results, roc_signature, fc_comparison,
     file = file.path(RDA_DIR, "07_validation.RData"))

cat("\nScript 7 complete. Validation results saved.\n")
