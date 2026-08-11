###############################################################################
# SCRIPT 8: Immune Infiltration Analysis (xCell)
# Multi-Cohort WGCNA: Breast Cancer Hub Gene Discovery
###############################################################################
library(xCell)
library(GSVA)
library(ggplot2)
library(dplyr)
library(tidyr)
library(corrplot)

options(stringsAsFactors = FALSE)

OUT_DIR <- "F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)"
FIG_DIR <- file.path(OUT_DIR, "figures/immune")
TAB_DIR <- file.path(OUT_DIR, "tables")
RDA_DIR <- file.path(OUT_DIR, "RData")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

load(file.path(RDA_DIR, "01_processed_data_common.RData"))
load(file.path(RDA_DIR, "01_processed_data_full.RData"))
load(file.path(RDA_DIR, "05_hub_genes.RData"))

datasets <- list(
  GSE42568 = processed_42568,
  GSE15852 = processed_15852,
  GSE29431 = processed_29431
)

immune_results <- list()

for (ds_name in names(datasets)) {
  cat(sprintf("\n========== %s: Immune Infiltration ==========\n", ds_name))
  
  expr <- datasets[[ds_name]]$expr
  meta <- datasets[[ds_name]]$meta
  
  # Run xCell
  xcell_scores <- xCellAnalysis(expr, rnaseq = FALSE)
  
  # Transpose: samples x cell types
  xcell_df <- as.data.frame(t(xcell_scores))
  
  # Add class info
  xcell_df$Class <- meta$class
  
  # Differential immune infiltration
  cell_types <- setdiff(colnames(xcell_df), "Class")
  diff_results <- data.frame()
  for (ct in cell_types) {
    normal_vals <- xcell_df[xcell_df$Class == "Normal", ct]
    tumor_vals <- xcell_df[xcell_df$Class == "Tumor", ct]
    tt <- tryCatch(t.test(tumor_vals, normal_vals), error = function(e) NULL)
    if (!is.null(tt)) {
      diff_results <- rbind(diff_results, data.frame(
        CellType = ct,
        Normal_mean = mean(normal_vals),
        Tumor_mean = mean(tumor_vals),
        logFC = tt$estimate[1] - tt$estimate[2],
        P_value = tt$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
  diff_results$FDR <- p.adjust(diff_results$P_value, method = "BH")
  diff_results <- diff_results[order(diff_results$P_value), ]
  
  cat(sprintf("  Significant cell types (FDR<0.05): %d/%d\n",
              sum(diff_results$FDR < 0.05), nrow(diff_results)))
  
  write.csv(diff_results, file.path(TAB_DIR, paste0("immune_diff_", ds_name, ".csv")), row.names = FALSE)
  
  # Box plots for top 8 significant cell types
  top_ct <- diff_results %>% filter(FDR < 0.05) %>% head(8)
  if (nrow(top_ct) > 0) {
    for (ct in top_ct$CellType) {
      plot_df <- data.frame(Score = xcell_df[[ct]], Class = xcell_df$Class)
      p <- ggplot(plot_df, aes(x = Class, y = Score, fill = Class)) +
        geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.2, alpha = 0.4, size = 1) +
        scale_fill_manual(values = c("Normal" = "blue", "Tumor" = "red")) +
        labs(title = paste(ds_name, "-", ct),
             subtitle = sprintf("P=%.2e", diff_results$P_value[diff_results$CellType == ct])) +
        theme_bw() + theme(legend.position = "none")
      ggsave(file.path(FIG_DIR, paste0(ds_name, "_", gsub("[ /]", "_", ct), "_boxplot.tiff")),
             p, width = 6, height = 6, dpi = 300, compression = "lzw")
    }
    
    # Heatmap of top cell types
    heat_data <- t(xcell_scores[top_ct$CellType, , drop = FALSE])
    colnames(heat_data) <- top_ct$CellType
    heat_df <- scale(heat_data)
    annotation_row <- data.frame(Class = meta$class, row.names = rownames(heat_df))
    
    tiff(file.path(FIG_DIR, paste0(ds_name, "_immune_heatmap.tiff")),
         width = 8, height = 6, units = "in", res = 300, compression = "lzw")
    pheatmap::pheatmap(t(heat_df), annotation_col = annotation_row,
                        annotation_colors = list(Class = c(Normal = "blue", Tumor = "red")),
                        main = paste(ds_name, "- Immune Infiltration"),
                        show_colnames = FALSE, fontsize_row = 8)
    dev.off()
  }
  
  # Correlation between hub genes and immune scores
  hub_expr <- expr[rownames(expr) %in% final_hubs, , drop = FALSE]
  common_hubs <- intersect(final_hubs, rownames(hub_expr))
  
  if (length(common_hubs) > 3) {
    hub_expr <- hub_expr[common_hubs, ]
    immune_mat <- xcell_scores
    
    # Correlation matrix
    cor_matrix <- cor(t(hub_expr), t(immune_mat), method = "spearman")
    
    # Save correlation table
    cor_df <- as.data.frame(as.table(cor_matrix))
    colnames(cor_df) <- c("HubGene", "CellType", "SpearmanR")
    write.csv(cor_df, file.path(TAB_DIR, paste0("immune_hub_correlation_", ds_name, ".csv")), row.names = FALSE)
    
    # Heatmap of hub-immune correlations
    if (nrow(top_ct) > 0) {
      tiff(file.path(FIG_DIR, paste0(ds_name, "_hub_immune_corheatmap.tiff")),
           width = 10, height = 8, units = "in", res = 300, compression = "lzw")
      pheatmap::pheatmap(cor_matrix[, intersect(top_ct$CellType, colnames(cor_matrix)), drop = FALSE],
                          main = paste(ds_name, "- Hub Gene vs Immune Cell Correlation"),
                          fontsize_row = 7, fontsize_col = 7,
                          color = colorRampPalette(c("blue", "white", "red"))(100))
      dev.off()
    }
  }
  
  immune_results[[ds_name]] <- list(
    scores = xcell_scores,
    diff = diff_results
  )
}

save(immune_results, file = file.path(RDA_DIR, "08_immune_infiltration.RData"))
cat("\nScript 8 complete.\n")
