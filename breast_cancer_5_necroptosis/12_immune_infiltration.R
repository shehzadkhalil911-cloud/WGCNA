# ============================================================
# Step 12: Immune Cell Infiltration Analysis
# ============================================================
set.seed(123)
library(xCell)
library(ggplot2)
library(pheatmap)
library(corrplot)

base_dir <- normalizePath(file.path(".", "Breast Cancer 5 Necroptosis"))
output_dir <- base_dir

sink(file.path(output_dir, "logs", "12_immune_infiltration.log"), split=TRUE)
cat("IMMUNE INFILTRATION ANALYSIS\n")
cat("============================\n\n")

save_tiff_600dpi <- function(plot, filename, width = 7, height = 6) {
  ggsave(filename=filename, plot=plot, width=width, height=height,
         units="in", dpi=600, compression="lzw", bg="white", device="tiff")
}

combat_expr <- readRDS(file.path(output_dir, "data_processed", "discovery_combat_corrected.rds"))
merged_meta <- readRDS(file.path(output_dir, "data_processed", "discovery_metadata.rds"))

# xCell immune infiltration
cat("Running xCell...\n")
xcell_results <- xCellAnalysis(as.matrix(combat_expr), 
                                 rnaseq=FALSE, 
                                 parallel.sz=1)
saveRDS(xcell_results, file.path(output_dir, "results", "immune_infiltration", "xcell_results.rds"))
write.csv(xcell_results, file.path(output_dir, "results", "immune_infiltration", "xcell_scores.csv"))

cat("xCell results:", nrow(xcell_results), "cell types x", ncol(xcell_results), "samples\n")

# Remove rows with all zeros
xcell_filt <- xcell_results[rowSums(xcell_results) > 0, ]
cat("Non-zero cell types:", nrow(xcell_filt), "\n")

# Tumor vs normal immune cell box plots
immune_cell_types <- rownames(xcell_filt)
plot_list <- list()
for(ct in immune_cell_types[1:min(12, length(immune_cell_types))]) {
  df <- data.frame(
    Score = as.numeric(xcell_filt[ct, ]),
    Group = merged_meta$class
  )
  p <- ggplot(df, aes(x=Group, y=Score, fill=Group)) +
    geom_boxplot() + theme_bw() +
    scale_fill_manual(values=c("Normal"="blue", "Tumor"="red")) +
    labs(title=ct, x="", y="xCell Score") +
    theme(legend.position="none", plot.title=element_text(hjust=0.5, size=10))
  plot_list[[ct]] <- p
}

if(length(plot_list) > 0) {
  library(gridExtra)
  n_col <- min(4, length(plot_list))
  n_row <- ceiling(length(plot_list)/n_col)
  p_grid <- grid.arrange(grobs=plot_list, ncol=n_col,
                          top="Immune Cell Infiltration: Tumor vs Normal")
  ggsave(file.path(output_dir, "results", "immune_infiltration", "immune_boxplots.tiff"),
         plot=p_grid, width=4*n_col, height=4*n_row, dpi=600, compression="lzw", bg="white")
}

# Immune infiltration heatmap
heat_annot <- data.frame(Group=merged_meta$class, row.names=colnames(combat_expr))
heat_colors <- list(Group=c(Normal="blue", Tumor="red"))

tiff(file.path(output_dir, "results", "immune_infiltration", "immune_heatmap.tiff"),
     width=12, height=8, units="in", res=600, compression="lzw", bg="white")
pheatmap(xcell_filt, annotation_col=heat_annot, annotation_colors=heat_colors,
         scale="row", main="Immune Cell Infiltration Heatmap",
         show_rownames=TRUE, show_colnames=FALSE,
         color=colorRampPalette(c("navy", "white", "firebrick3"))(100),
         fontsize_row=8, cutree_cols=2)
dev.off()

# Immune cell correlation heatmap
immune_cor <- cor(t(xcell_filt), method="spearman")
tiff(file.path(output_dir, "results", "immune_infiltration", "immune_correlation_heatmap.tiff"),
     width=8, height=7, units="in", res=600, compression="lzw", bg="white")
corrplot(immune_cor, method="color", type="upper", order="hclust",
         tl.cex=0.7, tl.col="black",
         col=colorRampPalette(c("navy", "white", "firebrick3"))(100),
         title="Immune Cell Correlation", mar=c(0,0,2,0))
dev.off()

# Differential immune cell analysis
diff_immune <- data.frame()
for(ct in rownames(xcell_filt)) {
  tumor_scores <- as.numeric(xcell_filt[ct, merged_meta$class == "Tumor"])
  normal_scores <- as.numeric(xcell_filt[ct, merged_meta$class == "Normal"])
  wt <- tryCatch(wilcox.test(tumor_scores, normal_scores), error=function(e) NULL)
  if(!is.null(wt)) {
    diff_immune <- rbind(diff_immune, data.frame(
      CellType = ct,
      Tumor_Mean = mean(tumor_scores),
      Normal_Mean = mean(normal_scores),
      p_value = wt$p.value,
      FDR = p.adjust(wt$p.value, "BH")
    ))
  }
}
diff_immune <- diff_immune[order(diff_immune$FDR), ]
write.csv(diff_immune, file.path(output_dir, "results", "immune_infiltration", "differential_immune_cells.csv"), row.names=FALSE)
cat("Differential immune cells:\n")
print(head(diff_immune, 10))

# ESTIMATE scores (stromal/immune/estimate) - skip if not available
if(require(estimate, quietly=TRUE)) {
  cat("Running ESTIMATE...\n")
  write.table(combat_expr, file.path(output_dir, "temp_estimate_input.txt"), sep="\t", quote=FALSE)
  estimateScore(file.path(output_dir, "temp_estimate_input.txt"),
                file.path(output_dir, "temp_estimate_output.gct"),
                platform="illumina")
  estimate_scores <- read.table(file.path(output_dir, "temp_estimate_output.gct"), skip=2, header=TRUE)
  unlink(file.path(output_dir, "temp_estimate_input.txt"))
  unlink(file.path(output_dir, "temp_estimate_output.gct"))
  
  if(ncol(estimate_scores) > 3) {
    estimate_df <- data.frame(
      Sample = colnames(estimate_scores)[-c(1:2)],
      StromalScore = as.numeric(estimate_scores[1, -c(1:2)]),
      ImmuneScore = as.numeric(estimate_scores[2, -c(1:2)]),
      ESTIMATEScore = as.numeric(estimate_scores[3, -c(1:2)])
    )
    write.csv(estimate_df, file.path(output_dir, "results", "immune_infiltration", "ESTIMATE_scores.csv"), row.names=FALSE)
    
    estimate_df$Group <- merged_meta$class[match(estimate_df$Sample, merged_meta$sample_id)]
    
    p_est <- ggplot(estimate_df, aes(x=Group, y=ESTIMATEScore, fill=Group)) +
      geom_boxplot() + theme_bw() +
      scale_fill_manual(values=c("Normal"="blue", "Tumor"="red")) +
      labs(title="ESTIMATE Score", x="", y="Score")
    save_tiff_600dpi(p_est, file.path(output_dir, "results", "immune_infiltration", "ESTIMATE_boxplot.tiff"), width=6, height=5)
  }
}

cat("\n\nIMMUNE INFILTRATION COMPLETE\n")
sink()
