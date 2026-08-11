# ============================================================
# Step 6: Necroptosis Activity Scoring
# ============================================================
set.seed(123)
library(GSVA)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

base_dir <- normalizePath(file.path(".", "Breast Cancer 5 Necroptosis"))
output_dir <- base_dir

sink(file.path(output_dir, "logs", "06_necroptosis_score.log"), split=TRUE)
cat("NECROPTOSIS ACTIVITY SCORING\n")
cat("============================\n\n")

save_tiff_600dpi <- function(plot, filename, width = 7, height = 6) {
  ggsave(filename=filename, plot=plot, width=width, height=height,
         units="in", dpi=600, compression="lzw", bg="white", device="tiff")
}

combat_expr <- readRDS(file.path(output_dir, "data_processed", "discovery_combat_corrected.rds"))
merged_meta <- readRDS(file.path(output_dir, "data_processed", "discovery_metadata.rds"))
ng_official <- readRDS(file.path(output_dir, "metadata", "necroptosis_genes_official.rds"))
ng_detected <- intersect(ng_official, rownames(combat_expr))

cat("Using", length(ng_detected), "necroptosis genes for GSVA scoring\n")

# GSVA necroptosis score
gsva_params <- gsvaParam(as.matrix(combat_expr), list(necroptosis=ng_detected),
                         maxDiff=TRUE, kcdf="Gaussian")
gsva_score <- gsva(gsva_params)
necroptosis_score <- as.numeric(gsva_score[1, ])
names(necroptosis_score) <- colnames(combat_expr)

saveRDS(necroptosis_score, file.path(output_dir, "results", "necroptosis_score", "necroptosis_score.rds"))
write.csv(data.frame(Sample=names(necroptosis_score), Score=necroptosis_score),
          file.path(output_dir, "results", "necroptosis_score", "necroptosis_score.csv"), row.names=FALSE)

# Box plot: tumor vs normal
score_df <- data.frame(
  Score = necroptosis_score,
  Group = merged_meta$class
)
library(ggpubr)
p_score_box <- ggplot(score_df, aes(x=Group, y=Score, fill=Group)) +
  geom_boxplot() + theme_bw() +
  scale_fill_manual(values=c("Normal"="blue", "Tumor"="red")) +
  labs(title="Necroptosis Activity Score", x="", y="GSVA Score") +
  stat_compare_means(method="wilcox.test", label="p.format")
save_tiff_600dpi(p_score_box, file.path(output_dir, "results", "necroptosis_score",
                 "score_boxplot.tiff"), width=6, height=5)

# Score density plot
p_density <- ggplot(score_df, aes(x=Score, fill=Group)) +
  geom_density(alpha=0.5) + theme_bw() +
  scale_fill_manual(values=c("Normal"="blue", "Tumor"="red")) +
  labs(title="Necroptosis Score Density", x="GSVA Score", y="Density")
save_tiff_600dpi(p_density, file.path(output_dir, "results", "necroptosis_score",
                 "score_density.tiff"), width=7, height=5)

# Score heatmap
score_mat <- matrix(necroptosis_score, nrow=1)
rownames(score_mat) <- "Necroptosis_Score"
colnames(score_mat) <- names(necroptosis_score)
heat_annot <- data.frame(Group=merged_meta$class, row.names=colnames(combat_expr))
heat_colors <- list(Group=c(Normal="blue", Tumor="red"))

tiff(file.path(output_dir, "results", "necroptosis_score", "score_heatmap.tiff"),
     width=12, height=3, units="in", res=600, compression="lzw", bg="white")
pheatmap(score_mat, annotation_col=heat_annot, annotation_colors=heat_colors,
         main="Necroptosis Score Across Samples",
         show_rownames=TRUE, show_colnames=FALSE,
         color=colorRampPalette(c("navy", "white", "firebrick3"))(100),
         fontsize=8, cluster_cols=TRUE, cluster_rows=FALSE)
dev.off()

# Correlation with DE-NRG expression
de_nrgs <- readLines(file.path(output_dir, "results", "necroptosis_genes", "DE_NRGs.txt"))
if(length(de_nrgs) > 0) {
  de_nrg_expr <- combat_expr[intersect(de_nrgs, rownames(combat_expr)), , drop=FALSE]
  
  score_cor <- apply(de_nrg_expr, 1, function(x) cor(x, necroptosis_score, method="spearman"))
  score_pval <- apply(de_nrg_expr, 1, function(x) cor.test(x, necroptosis_score, method="spearman")$p.value)
  
  cor_df <- data.frame(Gene=names(score_cor), Correlation=score_cor, Pvalue=score_pval)
  cor_df$FDR <- p.adjust(cor_df$Pvalue, "BH")
  cor_df <- cor_df[order(abs(cor_df$Correlation), decreasing=TRUE), ]
  write.csv(cor_df, file.path(output_dir, "results", "necroptosis_score", "DE_NRG_score_correlation.csv"), row.names=FALSE)
  
  # Plot
  p_cor <- ggplot(cor_df, aes(x=Correlation, y=-log10(FDR), label=Gene)) +
    geom_point(size=2, alpha=0.7) + theme_bw() +
    ggrepel::geom_text_repel(size=3, max.overlaps=15) +
    labs(title="DE-NRG vs Necroptosis Score Correlation", x="Spearman Correlation", y="-log10(FDR)")
  save_tiff_600dpi(p_cor, file.path(output_dir, "results", "necroptosis_score",
                   "DE_NRG_score_correlation.tiff"), width=9, height=7)
}

cat("\n\nNECROPTOSIS SCORING COMPLETE\n")
sink()
