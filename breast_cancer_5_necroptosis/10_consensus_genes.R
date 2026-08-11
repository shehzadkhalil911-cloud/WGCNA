# ============================================================
# Step 10: Consensus Necroptosis Genes
# ============================================================
set.seed(123)
library(ggplot2)
library(VennDiagram)
library(UpSetR)

base_dir <- normalizePath(file.path(".", "Breast Cancer 5 Necroptosis"))
output_dir <- base_dir

sink(file.path(output_dir, "logs", "10_consensus_genes.log"), split=TRUE)
cat("CONSENSUS GENES\n")
cat("===============\n\n")

save_tiff_600dpi <- function(plot, filename, width = 7, height = 6) {
  ggsave(filename=filename, plot=plot, width=width, height=height,
         units="in", dpi=600, compression="lzw", bg="white", device="tiff")
}

# Load gene lists from different methods
deg_sig <- read.csv(file.path(output_dir, "results", "tables", "DEG_significant.csv"), row.names=1)
lasso_genes <- readLines(file.path(output_dir, "results", "machine_learning", "lasso_selected_genes.txt"))
svm_genes <- readLines(file.path(output_dir, "results", "machine_learning", "svm_selected_genes.txt"))
rf_genes <- readLines(file.path(output_dir, "results", "machine_learning", "rf_selected_genes.txt"))
de_nrgs <- readLines(file.path(output_dir, "results", "necroptosis_genes", "DE_NRGs.txt"))

# PPI centrality genes
cent_file <- file.path(output_dir, "results", "networks", "centrality_measures.csv")
if(file.exists(cent_file)) {
  cent_df <- read.csv(cent_file)
  ppi_genes <- cent_df$Gene[cent_df$Degree > quantile(cent_df$Degree, 0.75, na.rm=TRUE)]
} else {
  ppi_genes <- character(0)
}

# Correlation analysis genes
cor_file <- file.path(output_dir, "results", "networks", "gene_correlation_matrix.csv")
if(file.exists(cor_file)) {
  cor_mat <- as.matrix(read.csv(cor_file, row.names=1))
  # Genes with high average correlation to other genes
  mean_cor <- apply(abs(cor_mat), 1, mean, na.rm=TRUE)
  cor_genes <- names(which(mean_cor > quantile(mean_cor, 0.75, na.rm=TRUE)))
} else {
  cor_genes <- character(0)
}

gene_lists <- list(
  LASSO = lasso_genes,
  SVM = svm_genes,
  RandomForest = rf_genes,
  DifferentialExpression = de_nrgs,
  PPICentrality = ppi_genes,
  Correlation = cor_genes
)

# Feature count plot
count_df <- data.frame(
  Method = names(gene_lists),
  Count = sapply(gene_lists, length)
)

p_count <- ggplot(count_df, aes(x=reorder(Method, Count), y=Count, fill=Method)) +
  geom_bar(stat="identity") + theme_bw() +
  coord_flip() + labs(title="Features Selected by Each Method", x="", y="Number of Genes") +
  geom_text(aes(label=Count), hjust=-0.2) +
  theme(legend.position="none")
save_tiff_600dpi(p_count, file.path(output_dir, "results", "consensus_genes", "feature_count.tiff"), width=8, height=5)

# Strict 3-method intersection (ML methods only)
ml_intersection <- Reduce(intersect, list(lasso_genes, svm_genes, rf_genes))
cat("Strict 3 ML method intersection:", length(ml_intersection), "\n")
if(length(ml_intersection) > 0) cat("  ", paste(ml_intersection, collapse=", "), "\n")

# Intersection of at least 2 ML methods
ml_union <- unique(c(lasso_genes, svm_genes, rf_genes))
ml_2_intersection <- ml_union[sapply(ml_union, function(g) {
  sum(g %in% lasso_genes, g %in% svm_genes, g %in% rf_genes) >= 2
})]
cat("Selected by >= 2 ML methods:", length(ml_2_intersection), "\n")

# Weighted consensus score
all_genes <- unique(unlist(gene_lists))
consensus_score <- sapply(all_genes, function(g) {
  score <- 0
  score <- score + (g %in% de_nrgs) * 3  # DE-NRG is primary
  score <- score + (g %in% lasso_genes) * 2
  score <- score + (g %in% svm_genes) * 2
  score <- score + (g %in% rf_genes) * 2
  score <- score + (g %in% ppi_genes) * 1
  score <- score + (g %in% cor_genes) * 1
  score
})
consensus_df <- data.frame(Gene = names(consensus_score), Score = consensus_score)
consensus_df <- consensus_df[order(consensus_df$Score, decreasing=TRUE), ]
write.csv(consensus_df, file.path(output_dir, "results", "consensus_genes", "consensus_scores.csv"), row.names=FALSE)

# Final locked signature: genes with consensus score >= threshold
top_quartile <- quantile(consensus_df$Score, 0.75)
final_genes <- consensus_df$Gene[consensus_df$Score >= top_quartile]
cat("Top quartile consensus genes:", length(final_genes), "\n")
writeLines(final_genes, file.path(output_dir, "results", "consensus_genes", "final_consensus_genes.txt"))
cat("  ", paste(final_genes, collapse=", "), "\n")

# Venn diagram for ML methods
if(length(lasso_genes) > 0 || length(svm_genes) > 0 || length(rf_genes) > 0) {
  tiff(file.path(output_dir, "results", "consensus_genes", "ML_venn.tiff"),
       width=7, height=7, units="in", res=600, compression="lzw", bg="white")
  g <- venn.diagram(
    x = list(LASSO = lasso_genes, SVM = svm_genes, RF = rf_genes),
    category.names = c("LASSO", "SVM", "RF"),
    filename = NULL,
    col = c("#E41A1C", "#377EB8", "#4DAF4A"),
    fill = c(rgb(228,26,28, max=255, alpha=80), 
             rgb(55,126,184, max=255, alpha=80),
             rgb(77,175,74, max=255, alpha=80)),
    cat.col = c("#E41A1C", "#377EB8", "#4DAF4A"),
    cat.cex = 1.5, cex = 1.5,
    main = "ML Method Intersection"
  )
  grid.draw(g)
  dev.off()
}

# UpSet plot
tiff(file.path(output_dir, "results", "consensus_genes", "upset_plot.tiff"),
     width=12, height=8, units="in", res=600, compression="lzw", bg="white")
upset(fromList(gene_lists), nsets=length(gene_lists), nintersects=30,
      order.by="freq", decreasing=TRUE,
      mainbar.y.label="Intersection Size", sets.x.label="Set Size",
      text.scale=c(1.3, 1.3, 1, 1, 1.5, 1.2))
dev.off()

# Method-by-gene heatmap
method_gene_mat <- sapply(gene_lists, function(gl) as.numeric(all_genes %in% gl))
rownames(method_gene_mat) <- all_genes
colnames(method_gene_mat) <- names(gene_lists)

# Show top consensus genes
top_n <- min(50, length(all_genes))
top_consensus <- consensus_df$Gene[1:top_n]
method_gene_top <- method_gene_mat[top_consensus, ]

tiff(file.path(output_dir, "results", "consensus_genes", "method_gene_heatmap.tiff"),
     width=8, height=10, units="in", res=600, compression="lzw", bg="white")
pheatmap::pheatmap(method_gene_top, cluster_rows=TRUE, cluster_cols=FALSE,
                   color=c("white", "steelblue"),
                   main="Method-by-Gene Selection Matrix",
                   display_numbers=TRUE, fontsize_number=8,
                   fontsize_row=8, fontsize_col=10)
dev.off()

# Circos plot for consensus genes
if(length(final_genes) >= 3 && length(final_genes) <= 50) {
  library(circlize)
  tiff(file.path(output_dir, "results", "consensus_genes", "consensus_circos.tiff"),
       width=10, height=10, units="in", res=600, compression="lzw", bg="white")
  circos.clear()
  
  # Prepare chord diagram data
  chord_data <- data.frame()
  for(m1 in names(gene_lists)) {
    for(m2 in names(gene_lists)) {
      if(m1 < m2) {
        overlap <- length(intersect(gene_lists[[m1]], gene_lists[[m2]]))
        if(overlap > 0) {
          chord_data <- rbind(chord_data, data.frame(from=m1, to=m2, value=overlap))
        }
      }
    }
  }
  
  if(nrow(chord_data) > 0) {
    chordDiagram(chord_data, transparency=0.3,
                 annotationTrack=c("grid", "axis"),
                 annotationTrackHeight=c(0.03, 0.01))
    title("Method Overlap Circos Plot")
  }
  dev.off()
  circos.clear()
}

cat("\n\nCONSENSUS GENES COMPLETE\n")
sink()
