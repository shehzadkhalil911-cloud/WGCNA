###############################################################################
# SCRIPT 3: Individual WGCNA Networks - Run Separately Per Dataset
# Multi-Cohort WGCNA: Breast Cancer Hub Gene Discovery
###############################################################################
args <- commandArgs(trailingOnly = TRUE)
dataset_arg <- ifelse(length(args) > 0, args[1], "GSE42568")

library(WGCNA)
library(flashClust)
library(RColorBrewer)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

OUT_DIR <- "F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)"
FIG_DIR <- file.path(OUT_DIR, "figures/WGCNA")
TAB_DIR <- file.path(OUT_DIR, "tables")
RDA_DIR <- file.path(OUT_DIR, "RData")
NET_DIR <- file.path(OUT_DIR, "networks")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)
dir.create(NET_DIR, recursive=TRUE, showWarnings=FALSE)

load(file.path(RDA_DIR, "01_processed_data_full.RData"))
load(file.path(RDA_DIR, "02_DE_analysis.RData"))

dataset_params <- list(
  GSE42568 = list(name = "42568", minModSize = 50, mergeCut = 0.20, deepSplit = 2),
  GSE15852 = list(name = "15852", minModSize = 30, mergeCut = 0.25, deepSplit = 2),
  GSE29431 = list(name = "29431", minModSize = 50, mergeCut = 0.20, deepSplit = 2)
)

params <- dataset_params[[dataset_arg]]
expr_orig <- get(paste0("processed_", params$name))$expr
meta <- get(paste0("processed_", params$name))$meta

cat(sprintf("\n========== WGCNA: %s ==========\n", dataset_arg))
cat(sprintf("Genes: %d, Samples: %d\n", nrow(expr_orig), ncol(expr_orig)))

# Use common discovery genes
load(file.path(RDA_DIR, "01_processed_data_common.RData"))
expr <- expr_orig[intersect(rownames(expr_orig), common_genes_discovery), , drop = FALSE]
cat(sprintf("Using %d common genes\n", nrow(expr)))

datExpr <- t(expr)
datExpr[] <- apply(datExpr, 2, as.numeric)

# Sample clustering
cat("Sample clustering...\n")
sampleTree <- flashClust(dist(datExpr), method = "average")
png(file.path(FIG_DIR, paste0(dataset_arg, "_sample_dendrogram.png")), width = 1600, height = 800)
par(cex = 0.6, mar = c(0,4,2,0))
plot(sampleTree, main = paste("Sample clustering -", dataset_arg), sub = "", xlab = "", 
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
traitColors <- ifelse(meta$class == "Tumor", "red", "blue")[sampleTree$order]
par(mar = c(1,4,0,0))
barplot(rep(1, length(traitColors)), col = traitColors, border = FALSE, axes = FALSE, space = 0)
mtext("Tumor/Normal", side = 2, line = 1)
dev.off()
cat("Sample dendrogram saved.\n")

# Soft threshold
cat("Selecting soft-threshold power...\n")
powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 1)
scale_indep <- -sign(sft$fitIndices[,3]) * sft$fitIndices[,2]
power_selected <- powers[which.max(scale_indep >= 0.80)]
if (length(power_selected) == 0 || power_selected < 4) power_selected <- 6
cat(sprintf("Selected power: %d\n", power_selected))

png(file.path(FIG_DIR, paste0(dataset_arg, "_soft_threshold.png")), width = 1600, height = 800)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], scale_indep, xlab = "Soft Threshold (power)", 
     ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", 
     main = paste("Scale independence -", dataset_arg))
text(sft$fitIndices[,1], scale_indep, labels = powers, col = "red", cex = 1.2)
abline(h = 0.80, col = "red", lty = 2)
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity -", dataset_arg))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = "red", cex = 1.2)
dev.off()
cat("Soft threshold plots saved.\n")

# Network construction
cat("Constructing network...\n")
bwnet <- blockwiseModules(datExpr, power = power_selected,
                          networkType = "signed", TOMType = "signed",
                          minModuleSize = params$minModSize,
                          reassignThreshold = 0, mergeCutHeight = params$mergeCut,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = TRUE, 
                          saveTOMFileBase = file.path(NET_DIR, paste0(dataset_arg, "_TOM")),
                          maxBlockSize = 2000, deepSplit = params$deepSplit,
                          verbose = 3)

geneTree <- bwnet$dendrograms[[1]]
moduleColors_initial <- labels2colors(bwnet$colors)

cat(sprintf("Initial modules: %d\n", length(unique(setdiff(moduleColors_initial, "grey")))))
cat("Colors:", unique(moduleColors_initial), "\n")

# Merge similar modules
cat("Merging similar modules...\n")
MERGe <- mergeCloseModules(datExpr, moduleColors_initial, cutHeight = params$mergeCut, verbose = 2)
moduleColors <- MERGe$colors
mergedMEs <- MERGe$newMEs
cat(sprintf("After merging: %d modules\n", length(unique(setdiff(moduleColors, "grey")))))

# Gene dendrogram - use the first block's dendrogram
if (length(bwnet$dendrograms) >= 1) {
  # For blockwise, use the first dendrogram with some genes
  png(file.path(FIG_DIR, paste0(dataset_arg, "_gene_dendrogram.png")), width = 1600, height = 800)
  plotDendroAndColors(geneTree, moduleColors[1:length(geneTree$order)],
                      "Module Colors",
                      dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                      main = paste("Gene Dendrogram -", dataset_arg))
  dev.off()
  cat("Gene dendrogram saved.\n")
}

# Module eigengenes
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs_order <- orderMEs(MEs)

# Module-trait relationships
trait <- as.data.frame(meta$class_binary)
colnames(trait) <- "Tumor_vs_Normal"
moduleTraitCor <- cor(MEs_order, trait, use = "p")
moduleTraitPval <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

png(file.path(FIG_DIR, paste0(dataset_arg, "_module_trait_heatmap.png")), width = 800, height = 1000)
textMatrix <- paste(sprintf("%.2f", moduleTraitCor), "\n(", sprintf("%.1e", moduleTraitPval), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = colnames(trait), yLabels = names(MEs_order),
               ySymbols = names(MEs_order), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-Trait Relationships -", dataset_arg))
dev.off()
cat("Module-trait heatmap saved.\n")

# Gene significance
GS <- as.data.frame(cor(datExpr, trait, use = "p"))
colnames(GS) <- "GS_Tumor"
GS_Pval <- as.data.frame(corPvalueStudent(as.matrix(GS), nrow(datExpr)))
MM <- as.data.frame(cor(datExpr, MEs_order, use = "p"))

# GS vs MM for top modules
modules_by_cor <- abs(moduleTraitCor[,1])
names(modules_by_cor) <- rownames(moduleTraitCor)
top_mods <- names(sort(modules_by_cor, decreasing = TRUE))
top_mods <- setdiff(top_mods, "MEgrey")[1:min(10, length(top_mods))]

for (mod_name in top_mods) {
  mod <- gsub("^ME", "", mod_name)
  module_genes <- moduleColors == mod
  if (sum(module_genes) < 5) next
  mm_vals <- abs(cor(datExpr[, module_genes, drop=FALSE], MEs_order[[mod_name]], use = "p"))[,1]
  gs_vals <- abs(GS[module_genes, 1])
  png(file.path(FIG_DIR, paste0(dataset_arg, "_GS_vs_MM_", mod, ".png")), width = 600, height = 600)
  verboseScatterplot(mm_vals, gs_vals, xlab = "Module Membership", ylab = "Gene Significance",
                     main = paste("Module", mod, "-", dataset_arg), col = mod, abline = TRUE)
  dev.off()
}

# Module sizes
module_sizes <- sort(table(moduleColors), decreasing = TRUE)
png(file.path(FIG_DIR, paste0(dataset_arg, "_module_sizes.png")), width = 1000, height = 600)
par(mar = c(8, 5, 4, 2))
barplot(module_sizes, las = 2, col = names(module_sizes),
        main = paste("Module Sizes -", dataset_arg), ylab = "Number of Genes", cex.names = 0.7)
dev.off()

# Module eigengene adjacency
ME_cor <- cor(MEs_order, use = "p")
png(file.path(FIG_DIR, paste0(dataset_arg, "_ME_adjacency.png")), width = 800, height = 800)
par(mar = c(8, 10, 3, 3))
labeledHeatmap(Matrix = ME_cor, xLabels = names(MEs_order), yLabels = names(MEs_order),
               colorLabels = FALSE, colors = blueWhiteRed(50), setStdMargins = FALSE,
               cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module Eigengene Adjacency -", dataset_arg))
dev.off()

# TOM heatmap sample
cat("Generating TOM heatmap sample...\n")
set.seed(123)
sample_idx <- sample(1:ncol(datExpr), min(300, ncol(datExpr)))
sample_TOM <- TOMsimilarityFromExpr(datExpr[, sample_idx], power = power_selected, 
                                     networkType = "signed", verbose = 0)
sample_dissTOM <- 1 - sample_TOM
# Create a simple hclust tree for the sample
sample_tree <- flashClust(as.dist(sample_dissTOM), method = "average")
png(file.path(FIG_DIR, paste0(dataset_arg, "_TOM_heatmap.png")), width = 1000, height = 1000)
TOMplot(sample_dissTOM, sample_tree, as.character(moduleColors[sample_idx]),
        main = paste("TOM Heatmap -", dataset_arg))
dev.off()

# Intramodular connectivity
intraModList <- list()
for (mod in setdiff(unique(moduleColors), "grey")) {
  in_mod <- moduleColors == mod
  if (sum(in_mod) < 3) next
  mod_adj <- adjacency(datExpr[, in_mod, drop=FALSE], power = power_selected, type = "signed")
  kWithin <- apply(mod_adj, 2, sum) - 1
  intraModList[[mod]] <- data.frame(Gene = colnames(datExpr)[in_mod], kWithin = kWithin, Module = mod)
}
intraModDF <- do.call(rbind, intraModList)
write.csv(intraModDF, file.path(TAB_DIR, paste0(dataset_arg, "_intramodular_connectivity.csv")), row.names = FALSE)

# Save results
wgcna_result <- list(
  dataset = dataset_arg, datExpr = datExpr, geneTree = geneTree,
  moduleColors = moduleColors, MEs = MEs_order,
  moduleTraitCor = moduleTraitCor, moduleTraitPval = moduleTraitPval,
  GS = GS, GS_Pval = GS_Pval, MM = MM,
  power = power_selected, sft = sft, intraModDF = intraModDF,
  n_modules = length(unique(setdiff(moduleColors, "grey"))),
  module_sizes = module_sizes
)
save(wgcna_result, file = file.path(RDA_DIR, paste0("03_WGCNA_", dataset_arg, ".RData")))
cat(sprintf("\n%s WGCNA complete. %d modules, power=%d\n", dataset_arg, wgcna_result$n_modules, power_selected))
