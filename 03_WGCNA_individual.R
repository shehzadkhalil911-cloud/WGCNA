###############################################################################
# SCRIPT 3: Individual WGCNA Networks
# Multi-Cohort WGCNA: Breast Cancer Hub Gene Discovery
###############################################################################
load("F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)/RData/01_processed_data_full.RData")
load("F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)/RData/02_DE_analysis.RData")

OUT_DIR <- "F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)"
FIG_DIR <- file.path(OUT_DIR, "figures/WGCNA")
TAB_DIR <- file.path(OUT_DIR, "tables")
RDA_DIR <- file.path(OUT_DIR, "RData")
NET_DIR <- file.path(OUT_DIR, "networks")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

library(WGCNA)
library(flashClust)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

allowWGCNAThreads()
options(stringsAsFactors = FALSE)

###############################################################################
# WGCNA Pipeline for One Dataset
###############################################################################
run_WGCNA <- function(expr, meta, dataset_name, 
                       powers = c(1:10, seq(12, 20, by = 2)),
                       minModuleSize = 30,
                       mergeCutHeight = 0.25,
                       deepSplit = 2,
                       maxBlockSize = 5000) {
  
  cat(sprintf("\n%s========== WGCNA: %s ==========%s\n", 
              paste(rep("=", 20), collapse=""), dataset_name, paste(rep("=", 20), collapse="")))
  cat(sprintf("Genes: %d, Samples: %d\n", nrow(expr), ncol(expr)))
  
  # Transpose: WGCNA expects samples in rows, genes in columns
  datExpr <- t(expr)
  datExpr[] <- apply(datExpr, 2, as.numeric)
  
  # --- Sample clustering ---
  cat("Sample clustering...\n")
  sampleTree <- flashClust(dist(datExpr), method = "average")
  
  png(file.path(FIG_DIR, paste0(dataset_name, "_sample_dendrogram.png")), 
      width = 1200, height = 600)
  par(cex = 0.6)
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = paste("Sample clustering -", dataset_name),
       sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  traitColors <- ifelse(meta$class == "Tumor", "red", "blue")
  plotOrder <- sampleTree$order
  traitColors_ordered <- traitColors[plotOrder]
  par(mar = c(1,4,0,0))
  barplot(rep(1, length(traitColors_ordered)), col = traitColors_ordered, 
          border = FALSE, axes = FALSE, space = 0)
  mtext("Tumor/Normal", side = 2, line = 1)
  dev.off()
  
  # --- Soft-threshold power selection ---
  cat("Selecting soft-threshold power...\n")
  sft <- pickSoftThreshold(datExpr, powerVector = powers, 
                           networkType = "signed", verbose = 2)
  
  scale_indep <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]
  power_selected <- powers[which.max(scale_indep >= 0.80)]
  if (length(power_selected) == 0 || power_selected < 4) {
    power_selected <- 6
  }
  cat(sprintf("Selected power: %d\n", power_selected))
  
  # Scale independence plot
  png(file.path(FIG_DIR, paste0(dataset_name, "_soft_threshold.png")), 
      width = 1200, height = 600)
  par(mfrow = c(1,2))
  plot(sft$fitIndices[,1], scale_indep,
       xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
       type = "n", main = paste("Scale independence -", dataset_name))
  text(sft$fitIndices[,1], scale_indep, labels = powers, col = "red", cex = 1.2)
  abline(h = 0.80, col = "red", lty = 2)
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
       type = "n", main = paste("Mean connectivity -", dataset_name))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = "red", cex = 1.2)
  dev.off()
  
  # --- Network construction using blockwise approach with TOM files ---
  cat("Constructing network via blockwiseModules...\n")
  bwnet <- blockwiseModules(datExpr, power = power_selected,
                            networkType = "signed",
                            TOMType = "signed",
                            minModuleSize = minModuleSize,
                            reassignThreshold = 0, 
                            mergeCutHeight = mergeCutHeight,
                            numericLabels = TRUE,
                            pamRespectsDendro = FALSE,
                            saveTOMs = TRUE,
                            saveTOMFileBase = file.path(NET_DIR, paste0(dataset_name, "_TOM")),
                            maxBlockSize = maxBlockSize,
                            deepSplit = deepSplit,
                            verbose = 3)
  geneTree <- bwnet$dendrograms[[1]]
  dynamicColors <- labels2colors(bwnet$colors)
  cat(sprintf("Blockwise complete. %d initial modules\n", length(unique(dynamicColors))))
  
  # --- Merge similar modules ---
  MERGe <- mergeCloseModules(datExpr, dynamicColors, 
                              cutHeight = mergeCutHeight, verbose = 3)
  mergedColors <- MERGe$colors
  moduleColors <- mergedColors
  mergedMEs <- MERGe$newMEs
  
  cat(sprintf("Number of modules detected: %d\n", length(unique(moduleColors))))
  cat("Modules:", unique(moduleColors), "\n")
  
  # --- Gene dendrogram with module colors ---
  png(file.path(FIG_DIR, paste0(dataset_name, "_gene_dendrogram.png")), 
      width = 1200, height = 600)
  plotDendroAndColors(geneTree, cbind(dynamicColors, moduleColors),
                      c("Dynamic Tree Cut", "Merged Modules"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = paste("Gene Dendrogram -", dataset_name))
  dev.off()
  
  # --- Calculate module eigengenes ---
  MEs <- moduleEigengenes(datExpr, moduleColors, excludeGrey = FALSE)$eigengenes
  MEs_order <- orderMEs(MEs)
  
  # --- Module-trait relationships ---
  trait <- as.data.frame(meta$class_binary)
  colnames(trait) <- "Tumor_vs_Normal"
  rownames(trait) <- rownames(datExpr)
  
  moduleTraitCor <- cor(MEs_order, trait, use = "p")
  moduleTraitPval <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
  
  # Module-trait heatmap
  png(file.path(FIG_DIR, paste0(dataset_name, "_module_trait_heatmap.png")), 
      width = 800, height = 1000)
  textMatrix <- paste(sprintf("%.2f", moduleTraitCor), "\n(",
                      sprintf("(%.2e)", moduleTraitPval), ")", sep = "")
  dim(textMatrix) <- dim(moduleTraitCor)
  par(mar = c(6, 10, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = colnames(trait),
                 yLabels = names(MEs_order),
                 ySymbols = names(MEs_order),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-Trait Relationships -", dataset_name))
  dev.off()
  
  # --- Gene significance and module membership ---
  GS <- as.data.frame(cor(datExpr, trait, use = "p"))
  colnames(GS) <- "GS_Tumor"
  GS_Pval <- as.data.frame(corPvalueStudent(as.matrix(GS), nrow(datExpr)))
  colnames(GS_Pval) <- "GS_Pval"
  
  # Module Membership (MM)
  MM <- as.data.frame(cor(datExpr, MEs_order, use = "p"))
  MM_Pval <- as.data.frame(corPvalueStudent(as.matrix(MM), nrow(datExpr)))
  
  # --- GS vs MM plots for top tumor-associated modules ---
  modules_by_cor <- moduleTraitCor[,1]
  names(modules_by_cor) <- rownames(moduleTraitCor)
  modules_to_plot <- names(sort(abs(modules_by_cor), decreasing = TRUE))
  modules_to_plot <- setdiff(modules_to_plot, "MEgrey")[1:min(10, length(modules_to_plot))]
  
  for (mod_name in modules_to_plot) {
    mod <- gsub("^ME", "", mod_name)
    if (mod == "grey") next
    module_genes <- moduleColors == mod
    if (sum(module_genes) < 5) next
    mod_ME <- MEs_order[[mod_name]]
    
    mm_values <- abs(cor(datExpr[, module_genes, drop=FALSE], mod_ME, use = "p"))[,1]
    gs_values <- abs(GS[module_genes, 1])
    
    png(file.path(FIG_DIR, paste0(dataset_name, "_GS_vs_MM_", mod, ".png")), 
        width = 600, height = 600)
    verboseScatterplot(mm_values, gs_values,
                       xlab = "Module Membership", 
                       ylab = "Gene Significance for Tumor",
                       main = paste("Module", mod, "-", dataset_name),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,
                       col = mod, abline = TRUE)
    dev.off()
  }
  
  # --- Module size bar plot ---
  module_sizes <- table(moduleColors)
  module_sizes <- sort(module_sizes, decreasing = TRUE)
  
  png(file.path(FIG_DIR, paste0(dataset_name, "_module_sizes.png")), 
      width = 1000, height = 600)
  par(mar = c(8, 5, 4, 2))
  barplot(module_sizes, las = 2, col = names(module_sizes),
          main = paste("Module Sizes -", dataset_name),
          ylab = "Number of Genes", cex.names = 0.7)
  dev.off()
  
  # --- Module eigengene adjacency heatmap ---
  png(file.path(FIG_DIR, paste0(dataset_name, "_ME_adjacency.png")), 
      width = 800, height = 800)
  ME_cor <- cor(MEs_order, use = "p")
  par(mar = c(8, 10, 3, 3))
  labeledHeatmap(Matrix = ME_cor,
                 xLabels = names(MEs_order),
                 yLabels = names(MEs_order),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module Eigengene Adjacency -", dataset_name))
  dev.off()
  
  # --- Intramodular connectivity ---
  intraModConnectivity <- list()
  for (mod in setdiff(unique(moduleColors), "grey")) {
    in_module <- moduleColors == mod
    if (sum(in_module) < 3) next
    mod_adj <- adjacency(t(datExpr[, in_module, drop=FALSE]), 
                          power = power_selected, type = "signed")
    kWithin <- apply(mod_adj, 2, sum) - 1
    intraModConnectivity[[mod]] <- data.frame(
      Gene = colnames(datExpr)[in_module],
      kWithin = kWithin,
      Module = mod
    )
  }
  intraModDF <- do.call(rbind, intraModConnectivity)
  write.csv(intraModDF, file.path(TAB_DIR, paste0(dataset_name, "_intramodular_connectivity.csv")), 
            row.names = FALSE)
  
  # --- Topological Overlap Heatmap (sample of 300 genes) ---
  cat("Generating TOM heatmap sample...\n")
  set.seed(123)
  n_genes <- ncol(datExpr)
  sample_idx <- sample(1:n_genes, min(300, n_genes))
  sample_dissTOM <- 1 - TOMsimilarityFromExpr(datExpr[, sample_idx, drop=FALSE],
                                                power = power_selected,
                                                networkType = "signed",
                                                verbose = 1)
  sample_colors <- moduleColors[sample_idx]
  png(file.path(FIG_DIR, paste0(dataset_name, "_TOM_heatmap.png")), 
      width = 1000, height = 1000)
  TOMplot(sample_dissTOM, geneTree, as.character(sample_colors),
          main = paste("TOM Heatmap -", dataset_name))
  dev.off()
  
  # Return results
  list(
    dataset = dataset_name,
    datExpr = datExpr,
    geneTree = geneTree,
    moduleColors = moduleColors,
    dynamicColors = dynamicColors,
    MEs = MEs_order,
    moduleTraitCor = moduleTraitCor,
    moduleTraitPval = moduleTraitPval,
    GS = GS,
    GS_Pval = GS_Pval,
    MM = MM,
    MM_Pval = MM_Pval,
    power = power_selected,
    sft = sft,
    intraModDF = intraModDF,
    n_modules = length(unique(setdiff(moduleColors, "grey"))),
    module_sizes = module_sizes
  )
}

###############################################################################
# Run WGCNA on each discovery dataset
###############################################################################
cat("\nRunning WGCNA on GSE42568...\n")
wgcna_42568 <- run_WGCNA(processed_42568$expr, processed_42568$meta, "GSE42568",
                          powers = c(1:10, seq(12, 20, by = 2)),
                          minModuleSize = 50, mergeCutHeight = 0.20, deepSplit = 2,
                          maxBlockSize = 3000)

cat("\nRunning WGCNA on GSE15852...\n")
wgcna_15852 <- run_WGCNA(processed_15852$expr, processed_15852$meta, "GSE15852",
                          powers = c(1:10, seq(12, 20, by = 2)),
                          minModuleSize = 30, mergeCutHeight = 0.25, deepSplit = 2,
                          maxBlockSize = 3000)

cat("\nRunning WGCNA on GSE29431...\n")
wgcna_29431 <- run_WGCNA(processed_29431$expr, processed_29431$meta, "GSE29431",
                          powers = c(1:10, seq(12, 20, by = 2)),
                          minModuleSize = 50, mergeCutHeight = 0.20, deepSplit = 2,
                          maxBlockSize = 3000)

###############################################################################
# Summary table
###############################################################################
cat("\n\n========== WGCNA Summary ==========\n")
for (ds_name in c("wgcna_42568", "wgcna_15852", "wgcna_29431")) {
  w <- get(ds_name)
  cat(sprintf("%s: %d modules, power=%d, genes=%d\n", 
              w$dataset, w$n_modules, w$power, ncol(w$datExpr)))
  cat("  Module sizes:", paste(names(w$module_sizes), "=", w$module_sizes, collapse=", "), "\n")
}

save(wgcna_42568, wgcna_15852, wgcna_29431,
     file = file.path(RDA_DIR, "03_WGCNA_individual.RData"))
cat("\nScript 3 complete.\n")
