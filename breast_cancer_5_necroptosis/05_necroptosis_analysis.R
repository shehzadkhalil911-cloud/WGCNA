# ============================================================
# Step 5: Necroptosis Gene Analysis
# ============================================================
set.seed(123)
library(limma)
library(ggplot2)
library(pheatmap)
library(VennDiagram)
library(RColorBrewer)
library(corrplot)

base_dir <- normalizePath(file.path(".", "Breast Cancer 5 Necroptosis"))
output_dir <- base_dir

sink(file.path(output_dir, "logs", "05_necroptosis_analysis.log"), split=TRUE)
cat("NECROPTOSIS GENE ANALYSIS\n")
cat("=========================\n\n")

save_tiff_600dpi <- function(plot, filename, width = 7, height = 6) {
  ggsave(filename=filename, plot=plot, width=width, height=height,
         units="in", dpi=600, compression="lzw", bg="white", device="tiff")
}

# Load data
combat_expr <- readRDS(file.path(output_dir, "data_processed", "discovery_combat_corrected.rds"))
merged_meta <- readRDS(file.path(output_dir, "data_processed", "discovery_metadata.rds"))
deg_sig <- read.csv(file.path(output_dir, "results", "tables", "DEG_significant.csv"), row.names=1)
deg_all <- read.csv(file.path(output_dir, "results", "tables", "DEG_all.csv"), row.names=1)
deg_up <- read.csv(file.path(output_dir, "results", "tables", "DEG_upregulated.csv"), row.names=1)
deg_down <- read.csv(file.path(output_dir, "results", "tables", "DEG_downregulated.csv"), row.names=1)
ng_official <- readRDS(file.path(output_dir, "metadata", "necroptosis_genes_official.rds"))

# Check necroptosis genes in expression data
ng_detected <- intersect(ng_official, rownames(combat_expr))
ng_not_detected <- setdiff(ng_official, rownames(combat_expr))
cat("Necroptosis genes detected in expression data:", length(ng_detected), "\n")
cat("Not detected:", length(ng_not_detected), "\n")
if(length(ng_not_detected) > 0) {
  cat("Non-detected genes:", paste(ng_not_detected, collapse=", "), "\n")
}

writeLines(ng_detected, file.path(output_dir, "results", "necroptosis_genes", "detected_genes.txt"))
writeLines(ng_not_detected, file.path(output_dir, "results", "necroptosis_genes", "non_detected_genes.txt"))

# DE-NRGs: intersect DEGs with necroptosis gene list
deg_sig_genes <- rownames(deg_sig)
de_nrgs <- intersect(deg_sig_genes, ng_detected)
cat("DE-NRGs (DEGs that are also necroptosis genes):", length(de_nrgs), "\n")

# Check up/down
de_nrg_up <- intersect(rownames(deg_up), de_nrgs)
de_nrg_down <- intersect(rownames(deg_down), de_nrgs)
cat("  Upregulated DE-NRGs:", length(de_nrg_up), "\n")
cat("  Downregulated DE-NRGs:", length(de_nrg_down), "\n")

writeLines(de_nrgs, file.path(output_dir, "results", "necroptosis_genes", "DE_NRGs.txt"))
writeLines(de_nrg_up, file.path(output_dir, "results", "necroptosis_genes", "DE_NRGs_up.txt"))
writeLines(de_nrg_down, file.path(output_dir, "results", "necroptosis_genes", "DE_NRGs_down.txt"))

# Save DE-NRG expression data
de_nrg_expr <- combat_expr[intersect(de_nrgs, rownames(combat_expr)), ]
write.csv(de_nrg_expr, file.path(output_dir, "results", "necroptosis_genes", "DE_NRG_expression.csv"))

# Save DE-NRG detailed table
de_nrg_table <- deg_all[intersect(de_nrgs, rownames(deg_all)), , drop=FALSE]
write.csv(de_nrg_table, file.path(output_dir, "results", "necroptosis_genes", "DE_NRG_table.csv"))

# DEG-NRG Venn diagram
tiff(file.path(output_dir, "results", "necroptosis_genes", "DEG_NRG_venn.tiff"),
     width=7, height=7, units="in", res=600, compression="lzw", bg="white")
g <- venn.diagram(
  x = list(DEGs = deg_sig_genes, Necroptosis = ng_detected),
  category.names = c("DEGs", "Necroptosis Genes"),
  filename = NULL,
  col = c("red", "blue"),
  fill = c(rgb(1,0,0,0.3), rgb(0,0,1,0.3)),
  cat.col = c("red", "blue"),
  cat.cex = 1.5,
  cex = 1.5,
  main = "DEG \u2229 Necroptosis Genes"
)
grid.draw(g)
dev.off()

# DE-NRG volcano plot
de_nrg_volcano <- deg_all[ng_detected, ]
de_nrg_volcano$status <- "NS"
de_nrg_volcano$status[de_nrg_volcano$adj.P.Val < 0.05 & de_nrg_volcano$logFC >= 1] <- "Up"
de_nrg_volcano$status[de_nrg_volcano$adj.P.Val < 0.05 & de_nrg_volcano$logFC <= -1] <- "Down"
de_nrg_volcano$label <- rownames(de_nrg_volcano)

p_de_nrg_volcano <- ggplot(de_nrg_volcano, aes(x=logFC, y=-log10(adj.P.Val), color=status)) +
  geom_point(size=2, alpha=0.8) +
  scale_color_manual(values=c("NS"="grey60", "Up"="red", "Down"="blue")) +
  theme_bw() +
  geom_vline(xintercept=c(-1,1), linetype="dashed", color="grey40") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="grey40") +
  labs(title="DE-NRG Volcano Plot", x="log2 Fold Change", y="-log10(adj.P.Val)") +
  ggrepel::geom_text_repel(aes(label=label), size=3, max.overlaps=20)
save_tiff_600dpi(p_de_nrg_volcano, file.path(output_dir, "results", "necroptosis_genes",
                 "DE_NRG_volcano.tiff"), width=9, height=7)

# DE-NRG heatmap
if(length(de_nrgs) > 1) {
  de_nrg_heat_expr <- combat_expr[de_nrgs, ]
  heat_annot <- data.frame(Group=merged_meta$class, row.names=colnames(combat_expr))
  heat_colors <- list(Group=c(Normal="blue", Tumor="red"))
  
  tiff(file.path(output_dir, "results", "necroptosis_genes", "DE_NRG_heatmap.tiff"),
       width=10, height=max(6, length(de_nrgs)*0.3), units="in", res=600, compression="lzw", bg="white")
  pheatmap(de_nrg_heat_expr, annotation_col=heat_annot, annotation_colors=heat_colors,
           scale="row", main="DE-NRG Expression Heatmap",
           show_rownames=TRUE, show_colnames=FALSE,
           color=colorRampPalette(c("navy", "white", "firebrick3"))(100),
           fontsize_row=8, cutree_cols=2)
  dev.off()
}

# DE-NRG expression box plots
if(length(de_nrgs) > 0) {
  n_plot <- min(9, length(de_nrgs))
  plot_genes <- de_nrgs[order(abs(deg_all[de_nrgs, "logFC"]), decreasing=TRUE)][1:n_plot]
  plot_list <- list()
  for(gene in plot_genes) {
    df <- data.frame(
      Expression = combat_expr[gene, ],
      Group = merged_meta$class
    )
    p <- ggplot(df, aes(x=Group, y=Expression, fill=Group)) +
      geom_boxplot() + theme_bw() +
      scale_fill_manual(values=c("Normal"="blue", "Tumor"="red")) +
      labs(title=gene, x="", y="Expression") +
      theme(legend.position="none", plot.title=element_text(hjust=0.5))
    plot_list[[gene]] <- p
  }
  if(length(plot_list) > 0) {
    library(gridExtra)
    p_grid <- grid.arrange(grobs=plot_list, ncol=3, top="DE-NRG Expression Boxplots")
    ggsave(file.path(output_dir, "results", "necroptosis_genes", "DE_NRG_boxplots.tiff"),
           plot=p_grid, width=10, height=3*ceiling(n_plot/3), dpi=600, compression="lzw", bg="white")
  }
}

# Fold change bar plot
if(length(de_nrgs) > 0) {
  fc_df <- data.frame(
    Gene = de_nrgs,
    logFC = deg_all[de_nrgs, "logFC"],
    Direction = ifelse(deg_all[de_nrgs, "logFC"] > 0, "Up", "Down")
  )
  fc_df <- fc_df[order(fc_df$logFC), ]
  fc_df$Gene <- factor(fc_df$Gene, levels=fc_df$Gene)
  
  p_fc <- ggplot(fc_df, aes(x=Gene, y=logFC, fill=Direction)) +
    geom_bar(stat="identity") + theme_bw() +
    scale_fill_manual(values=c("Up"="red", "Down"="blue")) +
    labs(title="DE-NRG Fold Changes", x="", y="log2 Fold Change") +
    theme(axis.text.x=element_text(angle=90, hjust=1, size=8)) +
    coord_flip()
  save_tiff_600dpi(p_fc, file.path(output_dir, "results", "necroptosis_genes",
                   "DE_NRG_foldchange.tiff"), width=8, height=max(5, length(de_nrgs)*0.25))
}

# Correlation heatmap among DE-NRGs
if(length(de_nrgs) > 1) {
  de_nrg_cor <- cor(t(de_nrg_heat_expr), method="spearman")
  tiff(file.path(output_dir, "results", "necroptosis_genes", "DE_NRG_correlation_heatmap.tiff"),
       width=8, height=7, units="in", res=600, compression="lzw", bg="white")
  corrplot(de_nrg_cor, method="color", type="upper", order="hclust",
           tl.cex=0.6, tl.col="black", 
           col=colorRampPalette(c("navy", "white", "firebrick3"))(100),
           title="DE-NRG Correlation", mar=c(0,0,2,0))
  dev.off()
}

cat("\n\nNECROPTOSIS GENE ANALYSIS COMPLETE\n")
sink()
