# ============================================================
# Step 3: Discovery Cohort Integration and Batch Correction
# ============================================================
set.seed(123)
library(limma)
library(sva)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

base_dir <- normalizePath(file.path(".", "Breast Cancer 5 Necroptosis"))
data_dir <- normalizePath(file.path(".", "Datasets", "data"))
output_dir <- base_dir

sink(file.path(output_dir, "logs", "03_batch_correction.log"), split=TRUE)
cat("BATCH CORRECTION\n")
cat("================\n\n")

save_tiff_600dpi <- function(plot, filename, width = 7, height = 6) {
  ggsave(filename=filename, plot=plot, width=width, height=height,
         units="in", dpi=600, compression="lzw", bg="white", device="tiff")
}

discovery_sets <- c("GSE42568", "GSE15852", "GSE29431")

# Load processed discovery datasets
processed <- readRDS(file.path(output_dir, "data_processed", "all_processed.rds"))

# Find common genes across discovery datasets
gene_lists <- lapply(discovery_sets, function(ds) rownames(processed[[ds]]$expression))
common_genes <- Reduce(intersect, gene_lists)
cat("Common genes across discovery datasets:", length(common_genes), "\n")

saveRDS(common_genes, file.path(output_dir, "metadata", "common_genes_discovery.rds"))

# Build merged expression matrix
merged_expr <- NULL
merged_meta <- NULL
merged_batch <- c()

for(ds in discovery_sets) {
  e <- processed[[ds]]$expression[common_genes, , drop=FALSE]
  m <- processed[[ds]]$metadata
  
  merged_expr <- cbind(merged_expr, e)
  merged_meta <- rbind(merged_meta, m)
  merged_batch <- c(merged_batch, rep(ds, ncol(e)))
}

colnames(merged_expr) <- merged_meta$sample_id
rownames(merged_meta) <- merged_meta$sample_id
merged_meta$batch <- merged_batch
merged_meta$class <- factor(merged_meta$class, levels=c("Normal", "Tumor"))

cat("Merged discovery matrix:", nrow(merged_expr), "genes x", ncol(merged_expr), "samples\n")
cat("Batch distribution:\n")
print(table(merged_batch))
cat("Class distribution:\n")
print(table(merged_meta$class))

# Check batch-phenotype confounding
cat("\n--- Batch-Phenotype Confounding Check ---\n")
conf_tbl <- table(merged_meta$batch, merged_meta$class)
print(conf_tbl)
write.csv(conf_tbl, file.path(output_dir, "results", "tables", "batch_class_contingency.csv"))

# PCA before batch correction
pca_before <- prcomp(t(merged_expr), scale.=TRUE, center=TRUE)
pca_before_df <- data.frame(PC1=pca_before$x[,1], PC2=pca_before$x[,2],
                            PC3=pca_before$x[,3], sample=colnames(merged_expr))
pca_before_df$batch <- merged_meta[colnames(merged_expr), "batch"]
pca_before_df$class <- merged_meta[colnames(merged_expr), "class"]
var_exp_before <- round(summary(pca_before)$importance[2, 1:3] * 100, 1)

p_pca_batch_before <- ggplot(pca_before_df, aes(x=PC1, y=PC2, color=batch)) +
  geom_point(size=2.5, alpha=0.8) + theme_bw() +
  labs(title="Before Batch Correction - Colored by Dataset",
       x=paste0("PC1 (", var_exp_before[1], "%)"), y=paste0("PC2 (", var_exp_before[2], "%)"))
save_tiff_600dpi(p_pca_batch_before, file.path(output_dir, "results", "batch_correction",
                 "PCA_before_correction_by_dataset.tiff"), width=8, height=6)

p_pca_class_before <- ggplot(pca_before_df, aes(x=PC1, y=PC2, color=class)) +
  geom_point(size=2.5, alpha=0.8) + theme_bw() +
  scale_color_manual(values=c("Normal"="blue", "Tumor"="red")) +
  labs(title="Before Batch Correction - Colored by Group",
       x=paste0("PC1 (", var_exp_before[1], "%)"), y=paste0("PC2 (", var_exp_before[2], "%)"))
save_tiff_600dpi(p_pca_class_before, file.path(output_dir, "results", "batch_correction",
                 "PCA_before_correction_by_group.tiff"), width=8, height=6)

# Correlation heatmap before correction
cor_before <- cor(merged_expr, method="spearman")
annotation_before <- data.frame(Dataset=merged_meta$batch, Group=merged_meta$class,
                                row.names=colnames(merged_expr))
ann_colors_before <- list(Dataset=c(GSE42568="#E41A1C", GSE15852="#377EB8", GSE29431="#4DAF4A"),
                          Group=c(Normal="blue", Tumor="red"))

tiff(file.path(output_dir, "results", "batch_correction", "correlation_heatmap_before.tiff"),
     width=10, height=9, units="in", res=600, compression="lzw", bg="white")
pheatmap(cor_before, annotation_col=annotation_before, annotation_colors=ann_colors_before,
         main="Sample Correlation Before Batch Correction",
         show_rownames=FALSE, show_colnames=FALSE,
         color=colorRampPalette(c("navy", "white", "firebrick3"))(100), fontsize=8)
dev.off()

# Dendrogram before correction
hc_before <- hclust(dist(t(merged_expr)), method="ward.D2")
tiff(file.path(output_dir, "results", "batch_correction", "dendrogram_before.tiff"),
     width=14, height=6, units="in", res=600, compression="lzw", bg="white")
plot(hc_before, main="Sample Dendrogram Before Batch Correction", xlab="", sub="", cex=0.4)
dev.off()

# ComBat batch correction
cat("\n--- Running ComBat ---\n")
mod <- model.matrix(~class, data=merged_meta)
combat_expr <- ComBat(dat=merged_expr, batch=merged_meta$batch, mod=mod, par.prior=TRUE, prior.plots=FALSE)
cat("ComBat completed. Corrected matrix:", nrow(combat_expr), "x", ncol(combat_expr), "\n")

# PCA after batch correction
pca_after <- prcomp(t(combat_expr), scale.=TRUE, center=TRUE)
pca_after_df <- data.frame(PC1=pca_after$x[,1], PC2=pca_after$x[,2],
                           PC3=pca_after$x[,3], sample=colnames(combat_expr))
pca_after_df$batch <- merged_meta[colnames(combat_expr), "batch"]
pca_after_df$class <- merged_meta[colnames(combat_expr), "class"]
var_exp_after <- round(summary(pca_after)$importance[2, 1:3] * 100, 1)

p_pca_batch_after <- ggplot(pca_after_df, aes(x=PC1, y=PC2, color=batch)) +
  geom_point(size=2.5, alpha=0.8) + theme_bw() +
  labs(title="After Batch Correction - Colored by Dataset",
       x=paste0("PC1 (", var_exp_after[1], "%)"), y=paste0("PC2 (", var_exp_after[2], "%)"))
save_tiff_600dpi(p_pca_batch_after, file.path(output_dir, "results", "batch_correction",
                 "PCA_after_correction_by_dataset.tiff"), width=8, height=6)

p_pca_class_after <- ggplot(pca_after_df, aes(x=PC1, y=PC2, color=class)) +
  geom_point(size=2.5, alpha=0.8) + theme_bw() +
  scale_color_manual(values=c("Normal"="blue", "Tumor"="red")) +
  labs(title="After Batch Correction - Colored by Group",
       x=paste0("PC1 (", var_exp_after[1], "%)"), y=paste0("PC2 (", var_exp_after[2], "%)"))
save_tiff_600dpi(p_pca_class_after, file.path(output_dir, "results", "batch_correction",
                 "PCA_after_correction_by_group.tiff"), width=8, height=6)

# Correlation heatmap after correction
cor_after <- cor(combat_expr, method="spearman")
tiff(file.path(output_dir, "results", "batch_correction", "correlation_heatmap_after.tiff"),
     width=10, height=9, units="in", res=600, compression="lzw", bg="white")
pheatmap(cor_after, annotation_col=annotation_before, annotation_colors=ann_colors_before,
         main="Sample Correlation After Batch Correction",
         show_rownames=FALSE, show_colnames=FALSE,
         color=colorRampPalette(c("navy", "white", "firebrick3"))(100), fontsize=8)
dev.off()

# Dendrogram after correction
hc_after <- hclust(dist(t(combat_expr)), method="ward.D2")
tiff(file.path(output_dir, "results", "batch_correction", "dendrogram_after.tiff"),
     width=14, height=6, units="in", res=600, compression="lzw", bg="white")
plot(hc_after, main="Sample Dendrogram After Batch Correction", xlab="", sub="", cex=0.4)
dev.off()

# Batch effect summary: variance explained by batch before vs after
batch_var_before <- apply(merged_expr, 1, function(x) {
  summary(lm(x ~ merged_meta$batch))$r.squared
})
batch_var_after <- apply(combat_expr, 1, function(x) {
  summary(lm(x ~ merged_meta$batch))$r.squared
})

batch_summary <- data.frame(
  Metric = c("Median R2 by batch", "Mean R2 by batch"),
  Before = c(median(batch_var_before), mean(batch_var_before)),
  After = c(median(batch_var_after), mean(batch_var_after))
)
cat("\n--- Batch Effect Summary ---\n")
print(batch_summary)
write.csv(batch_summary, file.path(output_dir, "results", "tables", "batch_effect_summary.csv"), row.names=FALSE)

# Batch effect boxplot
batch_df <- data.frame(
  R2 = c(batch_var_before, batch_var_after),
  Phase = rep(c("Before", "After"), each=length(batch_var_before))
)
p_batch <- ggplot(batch_df, aes(x=Phase, y=R2, fill=Phase)) +
  geom_boxplot() + theme_bw() +
  scale_fill_manual(values=c("Before"="#E41A1C", "After"="#377EB8")) +
  labs(title="Batch Effect: Variance Explained by Batch", x="", y="R-squared")
save_tiff_600dpi(p_batch, file.path(output_dir, "results", "batch_correction",
                 "batch_effect_summary.tiff"), width=7, height=6)

# Save combat-corrected data
saveRDS(combat_expr, file.path(output_dir, "data_processed", "discovery_combat_corrected.rds"))
saveRDS(merged_meta, file.path(output_dir, "data_processed", "discovery_metadata.rds"))
write.csv(combat_expr, file.path(output_dir, "data_processed", "discovery_combat_expression.csv"))
write.csv(merged_meta, file.path(output_dir, "data_processed", "discovery_metadata.csv"), row.names=FALSE)

cat("\n\nBATCH CORRECTION COMPLETE\n")
sink()
