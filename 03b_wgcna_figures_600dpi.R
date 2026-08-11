###############################################################################
# WGCNA Figure Regeneration: All figures at 600 dpi TIFF
###############################################################################
library(WGCNA)
library(ggplot2)
library(pheatmap)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

OUT_DIR <- "F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)"
FIG_DIR <- file.path(OUT_DIR, "figures/WGCNA")
RDA_DIR <- file.path(OUT_DIR, "RData")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

plot_wgcna_figures <- function(wgcna_obj, dataset_name) {
  datExpr <- wgcna_obj$datExpr
  moduleColors <- wgcna_obj$moduleColors
  MEs <- wgcna_obj$MEs
  geneTree <- wgcna_obj$geneTree
  moduleTraitCor <- wgcna_obj$moduleTraitCor
  moduleTraitPval <- wgcna_obj$moduleTraitPval
  GS <- wgcna_obj$GS
  MM <- wgcna_obj$MM
  sft <- wgcna_obj$sft
  power <- wgcna_obj$power
  intraModDF <- wgcna_obj$intraModDF
  module_sizes <- wgcna_obj$module_sizes
  
  cat(sprintf("Plotting %s (%d genes, %d modules, power=%d)\n",
              dataset_name, ncol(datExpr), length(unique(moduleColors)), power))
  
  # 1. Sample dendrogram
  tiff(file.path(FIG_DIR, paste0(dataset_name, "_sample_dendrogram.tiff")),
       width = 12, height = 6, units = "in", res = 600, compression = "lzw")
  sampleTree <- hclust(dist(datExpr), method = "average")
  par(mar = c(4, 4, 2, 1))
  plot(sampleTree, main = paste(dataset_name, "- Sample Dendrogram"),
       xlab = "", sub = "", cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.5)
  abline(h = mean(sampleTree$height) + 2 * sd(sampleTree$height), col = "red", lty = 2)
  dev.off()
  
  # 2. Soft threshold
  tiff(file.path(FIG_DIR, paste0(dataset_name, "_soft_threshold.tiff")),
       width = 12, height = 6, units = "in", res = 600, compression = "lzw")
  par(mfrow = c(1, 2))
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
       type = "n", main = paste("Scale independence"))
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = sft$fitIndices[, 1], col = "red")
  abline(h = 0.85, col = "red", lty = 2)
  
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
       type = "n", main = paste("Mean connectivity"))
  text(sft$fitIndices[, 1], sft$fitIndices[, 5],
       labels = sft$fitIndices[, 1], col = "red")
  dev.off()
  
  # 3. Gene dendrogram with module colors
  tiff(file.path(FIG_DIR, paste0(dataset_name, "_gene_dendrogram.tiff")),
       width = 12, height = 6, units = "in", res = 600, compression = "lzw")
  plotDendroAndColors(geneTree, moduleColors[geneTree$order],
                      "Module Colors", dendroLabels = FALSE,
                      main = paste(dataset_name, "- Gene Dendrogram"),
                      marAll = c(1, 6, 3, 1), cex.colorLabels = 0.8, cex.main = 1.5)
  dev.off()
  
  # 4. Module-trait heatmap
  tiff(file.path(FIG_DIR, paste0(dataset_name, "_module_trait_heatmap.tiff")),
       width = 8, height = 10, units = "in", res = 600, compression = "lzw")
  par(mar = c(8, 10, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor, xLabels = "TumorVsNormal",
                 yLabels = names(MEs), ySymbols = names(MEs),
                 colorLabels = FALSE, colors = blueWhiteRed(50),
                 textMatrix = paste0(sprintf("%.2f", moduleTraitCor), "\n(",
                                     sprintf("%.2e", moduleTraitPval), ")"),
                 cex.text = 0.6, cex.lab = 0.8,
                 main = paste("Module-Trait Relationship -", dataset_name))
  dev.off()
  
  # 5. GS vs MM plots (top 8 modules)
  top_modules <- names(sort(as.numeric(module_sizes), decreasing = TRUE))[1:min(10, length(module_sizes))]
  for (mod in top_modules) {
    if (mod == "grey") next
    mod_genes <- names(moduleColors[moduleColors == mod])
    if (length(mod_genes) < 5) next
    
    mm_col <- paste0("MM.", mod)
    gs_col <- "GS.TumorVsNormal"
    
    if (!(mm_col %in% colnames(MM))) next
    
    plot_df <- data.frame(
      MM = MM[mod_genes, mm_col],
      GS = GS[mod_genes, gs_col]
    )
    cor_val <- cor(plot_df$MM, plot_df$GS, method = "pearson", use = "pairwise")
    
    p <- ggplot(plot_df, aes(x = MM, y = GS)) +
      geom_point(alpha = 0.6, size = 1.5, color = mod) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
      annotate("text", x = min(plot_df$MM, na.rm = TRUE), y = max(plot_df$GS, na.rm = TRUE),
               label = paste("cor =", round(cor_val, 3)), hjust = 0, size = 4) +
      labs(title = paste(dataset_name, "- Module", mod),
           x = paste("Module Membership in", mod),
           y = "Gene Significance (Tumor vs Normal)") +
      theme_bw()
    ggsave(file.path(FIG_DIR, paste0(dataset_name, "_GS_vs_MM_", mod, ".tiff")),
           p, width = 6, height = 6, dpi = 600, compression = "lzw")
  }
  
  # 6. Module sizes
  mod_df <- data.frame(Module = names(module_sizes), Size = as.numeric(module_sizes))
  mod_df$Module <- factor(mod_df$Module, levels = mod_df$Module[order(mod_df$Size)])
  p <- ggplot(mod_df, aes(x = Module, y = Size, fill = Module)) +
    geom_bar(stat = "identity") + coord_flip() +
    labs(title = paste(dataset_name, "- Module Sizes"), x = "Module", y = "Number of Genes") +
    theme_bw() + theme(legend.position = "none")
  ggsave(file.path(FIG_DIR, paste0(dataset_name, "_module_sizes.tiff")),
         p, width = 10, height = 6, dpi = 600, compression = "lzw")
  
  # 7. ME Adjacency heatmap
  tiff(file.path(FIG_DIR, paste0(dataset_name, "_ME_adjacency.tiff")),
       width = 8, height = 8, units = "in", res = 600, compression = "lzw")
  ME_cor <- cor(MEs, use = "pairwise")
  pheatmap(ME_cor, main = paste(dataset_name, "- Module Eigengene Adjacency"),
           fontsize = 8, display_numbers = TRUE, number_format = "%.2f",
           color = colorRampPalette(c("blue", "white", "red"))(50))
  dev.off()
  
  # 8. TOM heatmap (sample of 500 genes)
  tom_file <- file.path(OUT_DIR, "networks", paste0("TOM_", dataset_name, ".RData"))
  if (file.exists(tom_file)) {
    load(tom_file)
    n_genes <- ncol(datExpr)
    if (n_genes > 500) {
      set.seed(123)
      sample_idx <- sample(n_genes, 500)
    } else {
      sample_idx <- 1:n_genes
    }
    tom_subset <- as.matrix(TOM[sample_idx, sample_idx])
    rownames(tom_subset) <- colnames(datExpr)[sample_idx]
    colnames(tom_subset) <- colnames(datExpr)[sample_idx]
    geneTree_sub <- hclust(as.dist(1 - tom_subset), method = "average")
    tiff(file.path(FIG_DIR, paste0(dataset_name, "_TOM_heatmap.tiff")),
         width = 10, height = 10, units = "in", res = 600, compression = "lzw")
    TOMplot(tom_subset, geneTree_sub, moduleColors[sample_idx],
            main = paste(dataset_name, "- TOM Heatmap"))
    dev.off()
    rm(TOM)
  } else {
    cat(sprintf("  TOM file not found: %s\n", tom_file))
  }
}

# Load and plot each dataset
datasets <- c("GSE42568", "GSE15852", "GSE29431")
for (ds in datasets) {
  cat(sprintf("\n========== %s ==========\n", ds))
  load(file.path(RDA_DIR, paste0("03_WGCNA_", ds, ".RData")))
  plot_wgcna_figures(wgcna_result, ds)
  rm(wgcna_result)
}

cat("\nWGCNA figure regeneration complete.\n")
