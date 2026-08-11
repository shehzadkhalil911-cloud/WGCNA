###############################################################################
# SCRIPT 2: Differential Expression Analysis
# Multi-Cohort WGCNA: Breast Cancer Hub Gene Discovery
###############################################################################
load("F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)/RData/01_processed_data_common.RData")
load("F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)/RData/01_processed_data_full.RData")

OUT_DIR <- "F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)"
FIG_DIR <- file.path(OUT_DIR, "figures/DEA")
TAB_DIR <- file.path(OUT_DIR, "tables")
RDA_DIR <- file.path(OUT_DIR, "RData")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(dplyr)
library(tidyr)

perform_DE_analysis <- function(expr, meta, dataset_name, output_prefix) {
  cat(sprintf("\n========== DE Analysis: %s ==========\n", dataset_name))
  
  design <- model.matrix(~ class, data = meta)
  colnames(design) <- gsub("class", "", colnames(design))
  
  fit <- lmFit(expr, design)
  fit <- eBayes(fit)
  
  # Get results for Tumor vs Normal
  coef_name <- setdiff(colnames(design), "(Intercept)")
  results <- topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH")
  results$Gene <- rownames(results)
  
  # Add direction and significance
  results$Direction <- ifelse(results$logFC > 0, "Up", "Down")
  results$Significance <- "NS"
  results$Significance[results$adj.P.Val < 0.05 & abs(results$logFC) > 1] <- 
    ifelse(results$logFC[results$adj.P.Val < 0.05 & abs(results$logFC) > 1] > 0, "Up", "Down")
  
  n_up <- sum(results$adj.P.Val < 0.05 & results$logFC > 0)
  n_down <- sum(results$adj.P.Val < 0.05 & results$logFC < 0)
  n_up_fc <- sum(results$adj.P.Val < 0.05 & results$logFC > 1)
  n_down_fc <- sum(results$adj.P.Val < 0.05 & results$logFC < -1)
  
  cat(sprintf("  Upregulated (p<0.05): %d\n", n_up))
  cat(sprintf("  Downregulated (p<0.05): %d\n", n_down))
  cat(sprintf("  Upregulated (p<0.05 & |logFC|>1): %d\n", n_up_fc))
  cat(sprintf("  Downregulated (p<0.05 & |logFC|>1): %d\n", n_down_fc))
  
  # Save DEG table
  write.csv(results, file.path(TAB_DIR, paste0(output_prefix, "_DEG_table.csv")), row.names = FALSE)
  
  # --- Volcano Plot ---
  p <- ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey50")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    labs(title = paste(dataset_name, "- Volcano Plot"),
         subtitle = paste0("Up: ", n_up_fc, " | Down: ", n_down_fc, " | NS: ", nrow(results) - n_up_fc - n_down_fc),
         x = "log2(Fold Change)", y = "-log10(Adjusted P-value)") +
    theme_bw() + theme(legend.position = "right")
  
  # Label top genes
  top_genes <- results %>% filter(Significance != "NS") %>% arrange(adj.P.Val) %>% head(15)
  p <- p + geom_text_repel(data = top_genes, aes(label = Gene), size = 2.5, max.overlaps = 15)
  
  ggsave(file.path(FIG_DIR, paste0(output_prefix, "_volcano.tiff")), p, width = 10, height = 8, dpi = 600, compression = "lzw")
  
  # --- MA Plot ---
  ma_data <- results
  ma_data$AveExpr <- rowMeans(expr, na.rm = TRUE)
  ma_data$MA_sig <- ifelse(ma_data$adj.P.Val < 0.05, 
                           ifelse(ma_data$logFC > 0, "Up", "Down"), "NS")
  
  p_ma <- ggplot(ma_data, aes(x = AveExpr, y = logFC, color = MA_sig)) +
    geom_point(alpha = 0.4, size = 1) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey60")) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    labs(title = paste(dataset_name, "- MA Plot"),
         x = "Average Expression", y = "log2(Fold Change)") +
    theme_bw() + theme(legend.position = "none")
  
  ggsave(file.path(FIG_DIR, paste0(output_prefix, "_MA_plot.tiff")), p_ma, width = 8, height = 6, dpi = 600, compression = "lzw")
  
  # --- DEG Heatmap (top 100 DEGs) ---
  sig_genes <- results %>% filter(adj.P.Val < 0.05) %>% arrange(adj.P.Val) %>% head(100)
  if (nrow(sig_genes) > 10) {
    heat_data <- expr[rownames(expr) %in% sig_genes$Gene, , drop = FALSE]
    heat_data <- heat_data[match(intersect(sig_genes$Gene, rownames(heat_data)), rownames(heat_data)), ]
    
    annotation_col <- data.frame(Class = meta$class, row.names = colnames(expr))
    annotation_colors <- list(Class = c(Normal = "blue", Tumor = "red"))
    
    # Scale rows for heatmap
    heat_data_scaled <- t(scale(t(heat_data)))
    
    tiff(file.path(FIG_DIR, paste0(output_prefix, "_DEG_heatmap.tiff")),
         width = 10, height = 8, units = "in", res = 600, compression = "lzw")
    pheatmap::pheatmap(heat_data_scaled,
             annotation_col = annotation_col,
             annotation_colors = annotation_colors,
             show_rownames = FALSE, show_colnames = FALSE,
             main = paste(dataset_name, "- Top 100 DEGs"))
    dev.off()
  }
  
  list(results = results, n_up = n_up, n_down = n_down,
       n_up_fc = n_up_fc, n_down_fc = n_down_fc, fit = fit)
}

# Run DE analysis on full filtered discovery datasets
cat("Running DE analysis on full discovery datasets...\n")
de_42568 <- perform_DE_analysis(processed_42568$expr, processed_42568$meta, 
                                 "GSE42568", "GSE42568")
de_15852 <- perform_DE_analysis(processed_15852$expr, processed_15852$meta,
                                 "GSE15852", "GSE15852")
de_29431 <- perform_DE_analysis(processed_29431$expr, processed_29431$meta,
                                 "GSE29431", "GSE29431")

# Also run DE on common-gene datasets for meta-analysis
cat("\nRunning DE on common-gene datasets for meta-analysis...\n")
de_comm_42568 <- perform_DE_analysis(comm_42568$expr, comm_42568$meta,
                                      "GSE42568_common", "GSE42568_common")
de_comm_15852 <- perform_DE_analysis(comm_15852$expr, comm_15852$meta,
                                      "GSE15852_common", "GSE15852_common")
de_comm_29431 <- perform_DE_analysis(comm_29431$expr, comm_29431$meta,
                                      "GSE29431_common", "GSE29431_common")

###############################################################################
# Cross-dataset DEG comparison and meta-analysis
###############################################################################
cat("\n========== Cross-dataset DEG Comparison ==========\n")

# Get DEG lists (p < 0.05)
deg_42568 <- de_42568$results %>% filter(adj.P.Val < 0.05) %>% pull(Gene)
deg_15852 <- de_15852$results %>% filter(adj.P.Val < 0.05) %>% pull(Gene)
deg_29431 <- de_29431$results %>% filter(adj.P.Val < 0.05) %>% pull(Gene)

# Consistent DEGs (in at least 2 datasets)
all_deg <- list(GSE42568 = deg_42568, GSE15852 = deg_15852, GSE29431 = deg_29431)
consistent_deg <- names(which(table(unlist(all_deg)) >= 2))
common_deg_3 <- Reduce(intersect, all_deg)

cat(sprintf("DEGs in GSE42568: %d\n", length(deg_42568)))
cat(sprintf("DEGs in GSE15852: %d\n", length(deg_15852)))
cat(sprintf("DEGs in GSE29431: %d\n", length(deg_29431)))
cat(sprintf("DEGs in >=2 datasets: %d\n", length(consistent_deg)))
cat(sprintf("DEGs in all 3 datasets: %d\n", length(common_deg_3)))

# Venn diagram
library(ggvenn)
venn_data <- list(GSE42568 = deg_42568, GSE15852 = deg_15852, GSE29431 = deg_29431)
p_venn <- ggvenn(venn_data, fill_color = c("red", "blue", "green"), 
                 stroke_size = 0.5, set_name_size = 4)
ggsave(file.path(FIG_DIR, "DEG_venn_3_datasets.tiff"), p_venn, width = 8, height = 8, dpi = 600, compression = "lzw")

# Direction consistency
cat("\nDirection consistency analysis...\n")
deg_42568_dir <- de_42568$results %>% filter(adj.P.Val < 0.05) %>% dplyr::select(Gene, logFC)
deg_15852_dir <- de_15852$results %>% filter(adj.P.Val < 0.05) %>% dplyr::select(Gene, logFC)
deg_29431_dir <- de_29431$results %>% filter(adj.P.Val < 0.05) %>% dplyr::select(Gene, logFC)

common_sig <- Reduce(intersect, list(deg_42568_dir$Gene, deg_15852_dir$Gene, deg_29431_dir$Gene))
cat(sprintf("Genes significant in all 3: %d\n", length(common_sig)))

# Check direction consistency
direction_check <- data.frame(
  Gene = common_sig,
  logFC_42568 = deg_42568_dir$logFC[match(common_sig, deg_42568_dir$Gene)],
  logFC_15852 = deg_15852_dir$logFC[match(common_sig, deg_15852_dir$Gene)],
  logFC_29431 = deg_29431_dir$logFC[match(common_sig, deg_29431_dir$Gene)]
)
direction_check$consistent <- sign(direction_check$logFC_42568) == sign(direction_check$logFC_15852) &
  sign(direction_check$logFC_15852) == sign(direction_check$logFC_29431)
cat(sprintf("Genes with consistent direction across all 3: %d\n", sum(direction_check$consistent)))

# Meta-analysis using Stouffer's method on common genes
cat("\nMeta-analysis (common genes)...\n")
common_genes <- common_genes_discovery
meta_p <- data.frame(Gene = common_genes)

for (ds_name in c("comm_42568", "comm_15852", "comm_29431")) {
  res <- get(paste0("de_", ds_name))$results
  meta_p[[paste0("P_", ds_name)]] <- res$P.Value[match(meta_p$Gene, res$Gene)]
  meta_p[[paste0("logFC_", ds_name)]] <- res$logFC[match(meta_p$Gene, res$Gene)]
}

# Combined effect size (average logFC)
meta_p$meta_logFC <- rowMeans(meta_p[, grep("logFC", colnames(meta_p))], na.rm = TRUE)

# Stouffer's combined p-value
meta_p$meta_P <- apply(meta_p[, grep("^P_", colnames(meta_p))], 1, function(x) {
  z <- qnorm(1 - x/2)
  z_sum <- sum(z, na.rm = TRUE)
  k <- sum(!is.na(x))
  p_combined <- 2 * (1 - pnorm(abs(z_sum) / sqrt(k)))
  return(p_combined)
})
meta_p$meta_FDR <- p.adjust(meta_p$meta_P, method = "BH")
meta_p$Direction <- ifelse(meta_p$meta_logFC > 0, "Up", "Down")

# Consistent direction across all datasets
fc_cols <- grep("logFC", colnames(meta_p), value = TRUE)
meta_p$direction_consistent <- apply(meta_p[, fc_cols], 1, function(x) {
  all(sign(x) == sign(x[1]), na.rm = TRUE)
})

cat(sprintf("Meta-analysis significant (FDR<0.05): %d\n", sum(meta_p$meta_FDR < 0.05)))
cat(sprintf("Consistent direction + significant: %d\n", 
    sum(meta_p$meta_FDR < 0.05 & meta_p$direction_consistent)))

meta_p <- meta_p[order(meta_p$meta_P), ]
write.csv(meta_p, file.path(TAB_DIR, "meta_analysis_DEGs.csv"), row.names = FALSE)

# Save all DE results
save(de_42568, de_15852, de_29431, de_comm_42568, de_comm_15852, de_comm_29431,
     meta_p, direction_check,
     file = file.path(RDA_DIR, "02_DE_analysis.RData"))

cat("\nScript 2 complete.\n")
