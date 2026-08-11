# ============================================================
# Step 4: Differential Expression Analysis (limma)
# ============================================================
set.seed(123)
library(limma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

base_dir <- normalizePath(file.path(".", "Breast Cancer 5 Necroptosis"))
output_dir <- base_dir

sink(file.path(output_dir, "logs", "04_differential_expression.log"), split=TRUE)
cat("DIFFERENTIAL EXPRESSION ANALYSIS\n")
cat("================================\n\n")

save_tiff_600dpi <- function(plot, filename, width = 7, height = 6) {
  ggsave(filename=filename, plot=plot, width=width, height=height,
         units="in", dpi=600, compression="lzw", bg="white", device="tiff")
}

# Load batch-corrected discovery data
combat_expr <- readRDS(file.path(output_dir, "data_processed", "discovery_combat_corrected.rds"))
merged_meta <- readRDS(file.path(output_dir, "data_processed", "discovery_metadata.rds"))

# Configurable thresholds
pval_threshold <- 0.05
lfc_threshold <- 1.0

cat("Thresholds: adj.P.Val <", pval_threshold, "|log2FC| >=", lfc_threshold, "\n\n")

# Design matrix
group <- factor(merged_meta$class, levels=c("Normal", "Tumor"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Check sample alignment
stopifnot(all(colnames(combat_expr) == rownames(merged_meta)))

# Contrast
contrast_matrix <- makeContrasts(Tumor - Normal, levels=design)

# Fit linear model
fit <- lmFit(combat_expr, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Full results
deg_all <- topTable(fit2, number=Inf, sort.by="P")
cat("Total genes tested:", nrow(deg_all), "\n")
cat("DEGs (adj.P <", pval_threshold, "and |logFC| >=", lfc_threshold, "):",
    sum(deg_all$adj.P.Val < pval_threshold & abs(deg_all$logFC) >= lfc_threshold), "\n")

# Significant DEGs
deg_sig <- deg_all[deg_all$adj.P.Val < pval_threshold & abs(deg_all$logFC) >= lfc_threshold, ]
deg_up <- deg_sig[deg_sig$logFC > 0, ]
deg_down <- deg_sig[deg_sig$logFC < 0, ]

cat("Upregulated:", nrow(deg_up), "\n")
cat("Downregulated:", nrow(deg_down), "\n")

# Ranked gene list for GSEA
deg_ranked <- deg_all[order(deg_all$t, decreasing=TRUE), ]
ranked_list <- deg_ranked$t
names(ranked_list) <- rownames(deg_ranked)

write.csv(deg_all, file.path(output_dir, "results", "tables", "DEG_all.csv"))
write.csv(deg_sig, file.path(output_dir, "results", "tables", "DEG_significant.csv"))
write.csv(deg_up, file.path(output_dir, "results", "tables", "DEG_upregulated.csv"))
write.csv(deg_down, file.path(output_dir, "results", "tables", "DEG_downregulated.csv"))
saveRDS(ranked_list, file.path(output_dir, "results", "differential_expression", "ranked_gene_list.rds"))
write.csv(data.frame(Gene=names(ranked_list), Score=ranked_list),
          file.path(output_dir, "results", "differential_expression", "ranked_gene_list.csv"), row.names=FALSE)

cat("\nSaved DEG tables to results/tables/\n")

# Volcano plot
deg_all$sig <- "NS"
deg_all$sig[deg_all$adj.P.Val < pval_threshold & deg_all$logFC >= lfc_threshold] <- "Up"
deg_all$sig[deg_all$adj.P.Val < pval_threshold & deg_all$logFC <= -lfc_threshold] <- "Down"
deg_all$sig_label <- ""
top_genes <- rownames(deg_sig[order(abs(deg_sig$logFC), decreasing=TRUE), ])[1:15]
deg_all$sig_label[rownames(deg_all) %in% top_genes] <- rownames(deg_all)[rownames(deg_all) %in% top_genes]

p_volcano <- ggplot(deg_all, aes(x=logFC, y=-log10(adj.P.Val), color=sig)) +
  geom_point(alpha=0.6, size=0.8) +
  scale_color_manual(values=c("NS"="grey60", "Up"="red", "Down"="blue")) +
  theme_bw() +
  geom_vline(xintercept=c(-lfc_threshold, lfc_threshold), linetype="dashed", color="grey40") +
  geom_hline(yintercept=-log10(pval_threshold), linetype="dashed", color="grey40") +
  labs(title="Volcano Plot: Tumor vs Normal", x="log2 Fold Change", y="-log10(adjusted p-value)") +
  ggrepel::geom_text_repel(aes(label=sig_label), size=2.5, max.overlaps=10) +
  theme(legend.position="right")
save_tiff_600dpi(p_volcano, file.path(output_dir, "results", "differential_expression",
                 "volcano_plot.tiff"), width=9, height=7)

# MA plot
deg_all$AveExpr <- deg_all$AveExpr
p_ma <- ggplot(deg_all, aes(x=AveExpr, y=logFC, color=sig)) +
  geom_point(alpha=0.6, size=0.8) +
  scale_color_manual(values=c("NS"="grey60", "Up"="red", "Down"="blue")) +
  theme_bw() +
  geom_hline(yintercept=0, linetype="solid", color="black") +
  labs(title="MA Plot: Tumor vs Normal", x="Average Expression", y="log2 Fold Change")
save_tiff_600dpi(p_ma, file.path(output_dir, "results", "differential_expression",
                 "MA_plot.tiff"), width=9, height=7)

# Top DEG heatmap
n_top <- min(50, nrow(deg_sig))
if(n_top > 0) {
  top_genes_heat <- rownames(deg_sig)[order(abs(deg_sig$logFC), decreasing=TRUE)][1:n_top]
  heat_expr <- combat_expr[top_genes_heat, ]
  
  heat_annot <- data.frame(Group=merged_meta$class, row.names=colnames(combat_expr))
  heat_colors <- list(Group=c(Normal="blue", Tumor="red"))
  
  tiff(file.path(output_dir, "results", "differential_expression", "top_DEG_heatmap.tiff"),
       width=10, height=12, units="in", res=600, compression="lzw", bg="white")
  pheatmap(heat_expr, annotation_col=heat_annot, annotation_colors=heat_colors,
           scale="row", main="Top 50 Differentially Expressed Genes",
           show_rownames=TRUE, show_colnames=FALSE,
           color=colorRampPalette(c("navy", "white", "firebrick3"))(100),
           fontsize_row=6, cutree_cols=2)
  dev.off()
}

# DEG box plots for top genes
if(n_top > 0) {
  top5 <- rownames(deg_sig)[order(abs(deg_sig$logFC), decreasing=TRUE)][1:5]
  for(gene in top5) {
    df <- data.frame(
      Expression = combat_expr[gene, ],
      Group = merged_meta$class
    )
    p_box <- ggplot(df, aes(x=Group, y=Expression, fill=Group)) +
      geom_boxplot() + theme_bw() +
      scale_fill_manual(values=c("Normal"="blue", "Tumor"="red")) +
      labs(title=paste(gene, "- Expression"), x="", y="Normalized Expression") +
      theme(legend.position="none")
    save_tiff_600dpi(p_box, file.path(output_dir, "results", "differential_expression",
                     paste0("boxplot_", gene, ".tiff")), width=5, height=5)
  }
}

# Up/Down count plot
count_df <- data.frame(
  Direction = c("Upregulated", "Downregulated"),
  Count = c(nrow(deg_up), nrow(deg_down))
)
p_count <- ggplot(count_df, aes(x=Direction, y=Count, fill=Direction)) +
  geom_bar(stat="identity") + theme_bw() +
  scale_fill_manual(values=c("Upregulated"="red", "Downregulated"="blue")) +
  labs(title="Differentially Expressed Genes", x="", y="Count") +
  geom_text(aes(label=Count), vjust=-0.3) +
  theme(legend.position="none")
save_tiff_600dpi(p_count, file.path(output_dir, "results", "differential_expression",
                 "DEG_count.tiff"), width=5, height=5)

cat("\n\nDIFFERENTIAL EXPRESSION ANALYSIS COMPLETE\n")
sink()
