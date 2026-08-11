# ============================================================
# Step 2: Data Preprocessing and Normalization
# Each dataset processed independently
# ============================================================
set.seed(123)
library(limma)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

base_dir <- normalizePath(file.path(".", "Breast Cancer 5 Necroptosis"))
data_dir <- normalizePath(file.path(".", "Datasets", "data"))
output_dir <- base_dir

sink(file.path(output_dir, "logs", "02_preprocessing.log"), split=TRUE)
cat("DATA PREPROCESSING\n")
cat("==================\n\n")

save_tiff_600dpi <- function(plot, filename, width = 7, height = 6) {
  ggsave(filename=filename, plot=plot, width=width, height=height,
         units="in", dpi=600, compression="lzw", bg="white", device="tiff")
}

process_dataset <- function(ds_name) {
  cat("\n--- Processing", ds_name, "---\n")
  
  csv_file <- file.path(data_dir, "processed", paste0(ds_name, "_expression.csv"))
  meta_file <- file.path(data_dir, "metadata", paste0(ds_name, "_sample_classes.csv"))
  
  expr <- read.csv(csv_file, row.names=1, check.names=FALSE)
  meta <- read.csv(meta_file)
  
  # Ensure alignment
  common_samples <- intersect(colnames(expr), meta$sample_id)
  expr <- expr[, common_samples, drop=FALSE]
  meta <- meta[match(common_samples, meta$sample_id), , drop=FALSE]
  rownames(meta) <- meta$sample_id
  
  cat("  Expression matrix:", nrow(expr), "genes x", ncol(expr), "samples\n")
  cat("  Expression range:", round(range(as.matrix(expr), na.rm=TRUE), 3), "\n")
  cat("  Expression mean:", round(mean(as.matrix(expr), na.rm=TRUE), 3), "\n")
  
  # Check if log2 transformation needed
  expr_mat <- as.matrix(expr)
  q75 <- apply(expr_mat, 2, function(x) quantile(x, 0.75, na.rm=TRUE))
  q25 <- apply(expr_mat, 2, function(x) quantile(x, 0.25, na.rm=TRUE))
  iqr_ratio <- median(q75 - q25)
  cat("  Median IQR:", round(iqr_ratio, 3), "\n")
  
  # Usually expression values > 20 suggest not log2 transformed
  max_val <- max(expr_mat, na.rm=TRUE)
  needs_log2 <- max_val > 30
  cat("  Max value:", round(max_val, 3), "-> Need log2 transform:", needs_log2, "\n")
  
  if(needs_log2 && min(expr_mat, na.rm=TRUE) >= 0) {
    expr_mat <- log2(expr_mat + 1)
    cat("  Applied log2(x+1) transformation\n")
  }
  
  # Remove genes with low variance
  gene_vars <- apply(expr_mat, 1, var, na.rm=TRUE)
  keep <- gene_vars > quantile(gene_vars, 0.05, na.rm=TRUE) & !is.na(gene_vars)
  expr_mat <- expr_mat[keep, ]
  cat("  After low-variance filter:", nrow(expr_mat), "genes retained\n")
  
  # Quantile normalization for platform consistency
  expr_norm <- normalizeBetweenArrays(expr_mat, method="quantile")
  cat("  Applied quantile normalization\n")
  
  # Sample QC: check for outliers
  sample_cors <- cor(expr_norm, method="spearman")
  mean_cor <- apply(sample_cors, 2, function(x) mean(x, na.rm=TRUE))
  outlier_samples <- names(which(mean_cor < quantile(mean_cor, 0.01)))
  cat("  Potential outlier samples (cor < 1st percentile):", length(outlier_samples), "\n")
  if(length(outlier_samples) > 0) cat("  ", paste(outlier_samples, collapse=", "), "\n")
  
  # Gene-level QC: remove genes with excessive missing
  na_frac <- apply(expr_norm, 1, function(x) mean(is.na(x)))
  keep_genes <- na_frac < 0.5
  expr_norm <- expr_norm[keep_genes, ]
  cat("  After gene-level QC:", nrow(expr_norm), "genes retained\n")
  
  # Pre-normalization boxplot
  pdf(NULL)
  expr_df <- as.data.frame(expr_mat[, 1:min(50, ncol(expr_mat))])
  expr_long <- stack(expr_df)
  colnames(expr_long) <- c("expression", "sample")
  p_before <- ggplot(expr_long, aes(x=sample, y=expression)) +
    geom_boxplot(fill="steelblue", outlier.size=0.5) +
    theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1, size=4)) +
    labs(title=paste(ds_name, "- Before Normalization"), x="", y="Expression")
  save_tiff_600dpi(p_before, file.path(output_dir, "results", "preprocessing",
                   paste0(ds_name, "_boxplot_before.tiff")), width=10, height=6)
  
  expr_df2 <- as.data.frame(expr_norm[, 1:min(50, ncol(expr_norm))])
  expr_long2 <- stack(expr_df2)
  colnames(expr_long2) <- c("expression", "sample")
  p_after <- ggplot(expr_long2, aes(x=sample, y=expression)) +
    geom_boxplot(fill="darkgreen", outlier.size=0.5) +
    theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1, size=4)) +
    labs(title=paste(ds_name, "- After Normalization"), x="", y="Expression")
  save_tiff_600dpi(p_after, file.path(output_dir, "results", "preprocessing",
                   paste0(ds_name, "_boxplot_after.tiff")), width=10, height=6)
  
  # Density plots
  dens_before <- ggplot(expr_long, aes(x=expression, color=sample)) +
    geom_density(alpha=0.3) + theme_bw() +
    theme(legend.position="none") +
    labs(title=paste(ds_name, "- Density Before"), x="Expression", y="Density")
  save_tiff_600dpi(dens_before, file.path(output_dir, "results", "preprocessing",
                   paste0(ds_name, "_density_before.tiff")), width=8, height=6)
  
  dens_after <- ggplot(expr_long2, aes(x=expression, color=sample)) +
    geom_density(alpha=0.3) + theme_bw() +
    theme(legend.position="none") +
    labs(title=paste(ds_name, "- Density After"), x="Expression", y="Density")
  save_tiff_600dpi(dens_after, file.path(output_dir, "results", "preprocessing",
                   paste0(ds_name, "_density_after.tiff")), width=8, height=6)
  
  # PCA with robust handling
  expr_norm_clean <- expr_norm
  expr_norm_clean[is.na(expr_norm_clean) | is.infinite(expr_norm_clean)] <- 0
  pca_try <- tryCatch({
    pca_res <- prcomp(t(expr_norm_clean), scale.=TRUE, center=TRUE)
    pca_df <- data.frame(PC1=pca_res$x[,1], PC2=pca_res$x[,2], 
                         sample=colnames(expr_norm))
    pca_df$group <- meta[colnames(expr_norm), "class"]
    var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)
    
    p_pca <- ggplot(pca_df, aes(x=PC1, y=PC2, color=group)) +
      geom_point(size=2.5, alpha=0.8) + theme_bw() +
      scale_color_manual(values=c("Normal"="blue", "Tumor"="red")) +
      labs(title=paste(ds_name, "- PCA"), 
           x=paste0("PC1 (", var_exp[1], "%)"), y=paste0("PC2 (", var_exp[2], "%)"))
    save_tiff_600dpi(p_pca, file.path(output_dir, "results", "preprocessing",
                     paste0(ds_name, "_PCA.tiff")), width=8, height=6)
    pca_res
  }, error=function(e) {
    cat("  PCA failed:", e$message, "- skipping PCA plot\n")
    NULL
  })
  
  # Sample correlation heatmap
  expr_norm_clean <- expr_norm
  expr_norm_clean[is.na(expr_norm_clean) | is.infinite(expr_norm_clean)] <- 0
  cor_mat <- cor(expr_norm_clean, method="spearman")
  annotation_col <- data.frame(Group=meta[colnames(cor_mat), "class"],
                               row.names=colnames(cor_mat))
  ann_colors <- list(Group=c(Normal="blue", Tumor="red"))
  
  tiff(file.path(output_dir, "results", "preprocessing", paste0(ds_name, "_correlation_heatmap.tiff")),
       width=10, height=9, units="in", res=600, compression="lzw", bg="white")
  pheatmap(cor_mat, annotation_col=annotation_col, annotation_colors=ann_colors,
           main=paste(ds_name, "- Sample Correlation"),
           show_rownames=FALSE, show_colnames=FALSE,
           color=colorRampPalette(c("navy", "white", "firebrick3"))(100),
           fontsize=8)
  dev.off()
  
  # Dendrogram
  hc_try <- tryCatch(hclust(dist(t(expr_norm_clean)), method="ward.D2"), error=function(e) NULL)
  if(!is.null(hc_try)) {
    tiff(file.path(output_dir, "results", "preprocessing", paste0(ds_name, "_dendrogram.tiff")),
         width=12, height=6, units="in", res=600, compression="lzw", bg="white")
    plot(hc_try, main=paste(ds_name, "- Sample Dendrogram"), 
         xlab="", sub="", cex=0.5)
    dev.off()
  }
  
  # Annotation retention plot
  probe_retention <- data.frame(
    Step = c("Raw", "LowVarFilter", "GeneQC"),
    Genes = c(nrow(expr), nrow(expr_mat), nrow(expr_norm))
  )
  p_ret <- ggplot(probe_retention, aes(x=Step, y=Genes, group=1)) +
    geom_line(color="blue", size=1) + geom_point(size=3, color="red") +
    theme_bw() + labs(title=paste(ds_name, "- Gene Retention"), x="", y="Number of Genes") +
    ylim(0, max(probe_retention$Genes)*1.1)
  save_tiff_600dpi(p_ret, file.path(output_dir, "results", "preprocessing",
                   paste0(ds_name, "_retention.tiff")), width=6, height=5)
  
  # Save processed
  write.csv(expr_norm, file.path(output_dir, "data_processed", paste0(ds_name, "_normalized.csv")))
  saveRDS(list(expression=expr_norm, metadata=meta, pca=pca_res),
          file.path(output_dir, "data_processed", paste0(ds_name, "_processed.rds")))
  
  cat("  Saved processed data\n")
  
  list(expression=expr_norm, metadata=meta, pca=pca_try)
}

# Process each dataset
processed <- list()
for(ds in c("GSE42568", "GSE15852", "GSE29431", "GSE70947")) {
  processed[[ds]] <- process_dataset(ds)
}

# Save all processed data references
saveRDS(processed, file.path(output_dir, "data_processed", "all_processed.rds"))

cat("\n\nPREPROCESSING COMPLETE\n")
sink()
