###############################################################################
# Script 08: Feature Gene Intersection and Visualization
###############################################################################

set.seed(123)

library(ggplot2)
library(tidyverse)
library(VennDiagram)
library(ggvenn)
library(viridis)
library(RColorBrewer)
library(circlize)
library(grid)
library(gridExtra)

base_dir <- normalizePath("breast_cancer_integrated_analysis")
data_processed_dir <- file.path(base_dir, "data_processed")
fig_dir <- file.path(base_dir, "results", "figures", "Feature_Gene_Intersection")
tables_dir <- file.path(base_dir, "results", "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helper: Save 600 dpi ---
save_publication_figure <- function(plot, filename, width = 8, height = 7) {
  if (inherits(plot, "grob")) {
    tiff(file.path(fig_dir, paste0(filename, ".tiff")), 
         width = width, height = height, units = "in", res = 600, compression = "lzw")
    grid.draw(plot)
    dev.off()
    
    pdf(file.path(fig_dir, paste0(filename, ".pdf")), 
        width = width, height = height)
    grid.draw(plot)
    dev.off()
  } else {
    ggsave(file.path(fig_dir, paste0(filename, ".tiff")), plot = plot,
           width = width, height = height, dpi = 600, compression = "lzw")
    ggsave(file.path(fig_dir, paste0(filename, ".pdf")), plot = plot,
           width = width, height = height, dpi = 600)
  }
  cat("  Saved:", filename, "\n")
}

theme_publication <- theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  )

# --- Load ML results ---
cat("Loading ML feature selection results...\n")
ml_results <- readRDS(file.path(data_processed_dir, "ml_selection_results.rds"))
lasso_genes <- ml_results$lasso_genes
svm_genes <- ml_results$svm_genes
rf_genes <- ml_results$rf_genes

cat("LASSO genes:", length(lasso_genes), "\n")
cat("SVM genes:", length(svm_genes), "\n")
cat("RF genes:", length(rf_genes), "\n")

# --- Load expression data ---
exprs_after <- readRDS(file.path(data_processed_dir, "merged_exprs_after_batch.rds"))
pheno <- readRDS(file.path(data_processed_dir, "merged_pheno_with_batch.rds"))
deg_table <- readRDS(file.path(data_processed_dir, "limma_deg_table.rds"))

# ======================================
# 1. Bar plot: Number of selected genes per method
# ======================================
cat("\nGenerating feature gene bar plot...\n")
n_genes_df <- data.frame(
  Method = c("LASSO", "SVM", "Random Forest"),
  Count = c(length(lasso_genes), length(svm_genes), length(rf_genes))
)

p_bar <- ggplot(n_genes_df, aes(x = Method, y = Count, fill = Method)) +
  geom_bar(stat = "identity", width = 0.6, alpha = 0.85) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("LASSO" = "#e41a1c", "SVM" = "#377eb8", 
                                "Random Forest" = "#4daf4a")) +
  labs(title = "Number of Selected Feature Genes by Method",
       x = "Feature Selection Method", y = "Number of Genes") +
  theme_publication +
  ylim(0, max(n_genes_df$Count) * 1.2) +
  theme(legend.position = "none")

save_publication_figure(p_bar, "Figure_6A_Feature_gene_barplot")

# ======================================
# 2. Venn diagram: Intersection of feature genes
# ======================================
cat("Generating Venn diagram...\n")

# Clean gene lists (remove any NAs or empty)
lasso_genes <- lasso_genes[!is.na(lasso_genes) & nchar(lasso_genes) > 0]
svm_genes <- svm_genes[!is.na(svm_genes) & nchar(svm_genes) > 0]
rf_genes <- rf_genes[!is.na(rf_genes) & nchar(rf_genes) > 0]

# Create Venn diagram
venn_data <- list(
  LASSO = lasso_genes,
  SVM = svm_genes,
  RF = rf_genes
)

p_venn <- ggvenn(venn_data, 
                  fill_color = c("#e41a1c", "#377eb8", "#4daf4a"),
                  stroke_size = 0.5,
                  set_name_size = 4,
                  text_size = 3) +
  labs(title = "Intersection of Feature Genes Across Methods") +
  theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))

save_publication_figure(p_venn, "Figure_6B_Feature_gene_venn", width = 7, height = 6)

# Also create a traditional VennDiagram
pdf(file.path(fig_dir, "Figure_6B_Feature_gene_venn_traditional.pdf"),
    width = 7, height = 7)
venn.plot <- venn.diagram(
  x = venn_data,
  category.names = c("LASSO", "SVM", "RF"),
  filename = NULL,
  output = TRUE,
  fill = c("#e41a1c", "#377eb8", "#4daf4a"),
  alpha = 0.5,
  cex = 1.2,
  fontfamily = "sans",
  cat.cex = 1.2,
  cat.fontfamily = "sans",
  cat.default.pos = "outer"
)
grid.draw(venn.plot)
dev.off()

tiff(file.path(fig_dir, "Figure_6B_Feature_gene_venn_traditional.tiff"),
     width = 7, height = 7, units = "in", res = 600, compression = "lzw")
grid.draw(venn.plot)
dev.off()

# ======================================
# 3. Intersection of all three methods
# ======================================
cat("\nComputing gene intersections...\n")
intersection_three <- Reduce(intersect, list(lasso_genes, svm_genes, rf_genes))
intersection_lasso_svm <- intersect(lasso_genes, svm_genes)
intersection_lasso_rf <- intersect(lasso_genes, rf_genes)
intersection_svm_rf <- intersect(svm_genes, rf_genes)

cat("Intersection of all 3 methods:", length(intersection_three), "\n")
if (length(intersection_three) > 0) {
  cat("Intersected genes:", paste(intersection_three, collapse = ", "), "\n")
}

cat("LASSO ∩ SVM:", length(intersection_lasso_svm), "\n")
cat("LASSO ∩ RF:", length(intersection_lasso_rf), "\n")
cat("SVM ∩ RF:", length(intersection_svm_rf), "\n")

# Save final feature gene list
final_genes <- intersection_three

# If no triple intersection, use union or pairwise intersections
if (length(final_genes) == 0) {
  cat("No genes selected by all three methods.\n")
  cat("Using intersection of LASSO and SVM as final features.\n")
  final_genes <- intersection_lasso_svm
  
  if (length(final_genes) == 0) {
    cat("LASSO ∩ SVM is empty. Using RF top 20 genes.\n")
    rf_importance <- ml_results$rf_importance
    final_genes <- as.character(rf_importance$Gene[1:min(20, nrow(rf_importance))])
  }
}

cat("Final feature genes:", length(final_genes), "\n")

writeLines(final_genes, file.path(tables_dir, "final_feature_genes.txt"))
write.csv(data.frame(FeatureGene = final_genes), 
          file.path(tables_dir, "final_feature_genes.csv"), row.names = FALSE)

# Complete intersection table
intersection_summary <- data.frame(
  Comparison = c("LASSO ∩ SVM ∩ RF", "LASSO ∩ SVM", "LASSO ∩ RF", "SVM ∩ RF",
                  "LASSO only", "SVM only", "RF only"),
  Count = c(
    length(intersection_three),
    length(intersection_lasso_svm),
    length(intersection_lasso_rf),
    length(intersection_svm_rf),
    length(setdiff(lasso_genes, union(svm_genes, rf_genes))),
    length(setdiff(svm_genes, union(lasso_genes, rf_genes))),
    length(setdiff(rf_genes, union(lasso_genes, svm_genes)))
  ),
  stringsAsFactors = FALSE
)
write.csv(intersection_summary, 
          file.path(tables_dir, "feature_gene_intersection_summary.csv"), row.names = FALSE)

# ======================================
# 4. Box plots: expression of final feature genes
# ======================================
cat("Generating feature gene box plots...\n")

if (length(final_genes) > 0 && all(final_genes %in% rownames(exprs_after))) {
  # Create a combined data frame for box plots
  plot_data <- data.frame()
  for (gene in final_genes) {
    gene_df <- data.frame(
      Expression = exprs_after[gene, ],
      Gene = gene,
      Group = pheno$group,
      stringsAsFactors = FALSE
    )
    plot_data <- rbind(plot_data, gene_df)
  }
  
  # Box plot for each gene
  n_plots <- length(final_genes)
  n_col <- min(4, n_plots)
  n_row <- ceiling(n_plots / n_col)
  
  p_box <- ggplot(plot_data, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(outlier.size = 0.8, alpha = 0.7, notch = FALSE) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.4) +
    facet_wrap(~ Gene, scales = "free_y", ncol = n_col) +
    scale_fill_manual(values = c("Normal" = "#4daf4a", "Tumor" = "#e41a1c")) +
    labs(title = "Expression of Final Feature Genes: Tumor vs Normal",
         x = "", y = "Expression Level") +
    theme_publication +
    theme(
      strip.text = element_text(size = 9, face = "italic"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
    )
  
  save_publication_figure(p_box, "Figure_6C_Feature_gene_boxplot", 
                          width = min(12, 3 * n_col), 
                          height = min(10, 3 * n_row))
}

# ======================================
# 5. Volcano plot highlighting final feature genes
# ======================================
cat("Generating volcano plot with highlighted feature genes...\n")

deg_table$is_feature <- ifelse(rownames(deg_table) %in% final_genes, "Feature Gene", "Other")
deg_table$is_feature <- factor(deg_table$is_feature, levels = c("Other", "Feature Gene"))

p_volcano <- ggplot(deg_table, aes(x = logFC, y = -log10(adj.P.Val), 
                                    color = is_feature)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Other" = "grey60", "Feature Gene" = "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", linewidth = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue", linewidth = 0.5) +
  labs(title = "Volcano Plot Highlighting Final Feature Genes",
       x = expression(Log[2]~Fold~Change),
       y = expression(-Log[10]~(Adjusted~P~value)),
       color = "") +
  theme_publication +
  theme(legend.position = "bottom")

if (length(final_genes) > 0 && length(final_genes) <= 20) {
  p_volcano <- p_volcano +
    ggrepel::geom_text_repel(
      data = deg_table[deg_table$is_feature == "Feature Gene", ],
      aes(label = rownames(deg_table[deg_table$is_feature == "Feature Gene", ])),
      size = 3, max.overlaps = 20, box.padding = 0.5
    )
}

save_publication_figure(p_volcano, "Figure_6D_Feature_gene_volcano", 
                        width = 9, height = 7)

# ======================================
# 6. Circos plot (if final genes are available)
# ======================================
cat("Generating circos plot...\n")

if (length(final_genes) >= 3 && length(final_genes) <= 50) {
  tryCatch({
    # Prepare data for circos
    final_gene_expr <- exprs_after[final_genes, ]
    
    # Calculate logFC for each gene
    tumor_mean <- rowMeans(final_gene_expr[, pheno$group == "Tumor", drop = FALSE])
    normal_mean <- rowMeans(final_gene_expr[, pheno$group == "Normal", drop = FALSE])
    gene_logFC <- tumor_mean - normal_mean
    
    # Calculate -log10(pval)
    gene_pval <- -log10(deg_table[final_genes, "adj.P.Val"])
    
    # Create circos data
    circos_data <- data.frame(
      gene = final_genes,
      logFC = round(gene_logFC, 2),
      pval = round(gene_pval, 2),
      stringsAsFactors = FALSE
    )
    
    write.csv(circos_data, file.path(tables_dir, "feature_gene_circos_data.csv"), 
              row.names = FALSE)
    
    # Generate circos plot
    pdf(file.path(fig_dir, "Figure_6E_Feature_gene_circos.pdf"), 
        width = 10, height = 10)
    
    circos.clear()
    circos.par(start.degree = 90, gap.degree = 4)
    
    # Colors
    fc_colors <- ifelse(circos_data$logFC > 0, "#e41a1c", "#377eb8")
    names(fc_colors) <- circos_data$gene
    
    # Initialize
    chordDiagram(
      x = circos_data[, c("gene", "gene")],
      grid.col = fc_colors,
      transparency = 0.3,
      annotationTrack = "grid",
      preDraw = function() {
        par(family = "sans")
      },
      direction.type = c("arrows", "diffHeight"),
      directional = 1,
      diffHeight = 0.04
    )
    
    # Add labels
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.1,
                  CELL_META$sector.index,
                  facing = "clockwise",
                  niceFacing = TRUE,
                  adj = c(0, 0.5),
                  cex = 0.6)
    }, bg.border = NA)
    
    title("Feature Gene Circos Plot\nRed = Upregulated, Blue = Downregulated")
    
    dev.off()
    
    tiff(file.path(fig_dir, "Figure_6E_Feature_gene_circos.tiff"), 
         width = 10, height = 10, units = "in", res = 600, compression = "lzw")
    
    circos.clear()
    circos.par(start.degree = 90, gap.degree = 4)
    
    chordDiagram(
      x = circos_data[, c("gene", "gene")],
      grid.col = fc_colors,
      transparency = 0.3,
      annotationTrack = "grid",
      direction.type = c("arrows", "diffHeight"),
      directional = 1,
      diffHeight = 0.04
    )
    
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.1,
                  CELL_META$sector.index,
                  facing = "clockwise",
                  niceFacing = TRUE,
                  adj = c(0, 0.5),
                  cex = 0.6)
    }, bg.border = NA)
    
    title("Feature Gene Circos Plot\nRed = Upregulated, Blue = Downregulated")
    
    dev.off()
    
    circos.clear()
    cat("  Circos plot saved\n")
    
  }, error = function(e) {
    cat("  Circos plot error:", e$message, "\n")
    cat("  Skipping circos plot.\n")
  })
}

cat("\n========================================================\n")
cat("FEATURE GENE INTERSECTION COMPLETE\n")
cat("========================================================\n")
cat("Final feature genes:", length(final_genes), "\n")
cat("Results saved to results/figures/Feature_Gene_Intersection/ and results/tables/\n")

# Save final feature genes for downstream
saveRDS(final_genes, file.path(data_processed_dir, "final_feature_genes.rds"))
