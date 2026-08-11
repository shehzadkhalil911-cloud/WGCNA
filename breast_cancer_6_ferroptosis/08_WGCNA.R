source("R/utils.R")
load_libraries()
set.seed(123)

BASE_DIR <- "F:/independent projects/cancer/breast cancer"
OUT_DIR <- file.path(BASE_DIR, "Breast cancer 6 Ferroptosis")

log_msg("=== Step 8: WGCNA and Co-expression Analysis ===")

combat_expr <- readRDS(file.path(OUT_DIR, "data_processed/discovery/discovery_combat.rds"))
metadata <- read.csv(file.path(OUT_DIR, "data_processed/discovery/discovery_metadata.csv"))
ferr_score <- readRDS(file.path(OUT_DIR, "results/ferroptosis_score/ferroptosis_score_vector.rds"))

library(WGCNA)

# Filter to most variable genes (top 3000 for computational feasibility)
gene_vars <- apply(combat_expr, 1, var, na.rm = TRUE)
gene_vars_sorted <- sort(gene_vars, decreasing = TRUE)
n_top <- min(3000, length(gene_vars_sorted))
expr_filt <- combat_expr[names(gene_vars_sorted)[1:n_top], ]
log_msg(sprintf("Genes retained after variance filter: %d/%d (top %d most variable)", nrow(expr_filt), nrow(combat_expr), n_top))

# Check for missing values
gsg <- goodSamplesGenes(t(expr_filt), verbose = 0)
if (!gsg$allOK) {
  log_msg("Removing problematic genes/samples...")
  expr_filt <- expr_filt[gsg$goodGenes, gsg$goodSamples]
}

# Sample clustering and outlier detection
sample_tree <- hclust(dist(t(expr_filt)), method = "average")
open_tiff(file.path(OUT_DIR, "results/figures/Figure5A_sample_dendrogram.tiff"), width = 10, height = 5)
par(mar = c(4, 4, 2, 1))
plot(sample_tree, labels = FALSE, main = "Sample Clustering Dendrogram",
     xlab = "", sub = "")
# No outlier removal unless clearly justified
abline(h = 60, col = "red", lty = 2)
dev.off()

# Choose soft-threshold power
powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(t(expr_filt), powerVector = powers,
                          networkType = "signed", verbose = 0)

open_tiff(file.path(OUT_DIR, "results/figures/Figure5B_soft_threshold.tiff"), width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale Independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red", cex = 1)
abline(h = 0.85, col = "red", lty = 2)

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red", cex = 1)
dev.off()

softPower <- sft$powerEstimate
if (is.na(softPower) || softPower < 1) {
  softPower <- 6
  log_msg(sprintf("Power estimate ambiguous, using default: %d", softPower))
}
log_msg(sprintf("Selected soft threshold power: %d", softPower))

# Manual network construction
cor_st <- function(x, y = NULL, use = "everything", method = "pearson") {
  stats::cor(x, y = y, use = use, method = method)
}
adj_mat <- adjacency(t(expr_filt), power = softPower,
                      type = "signed", corFnc = "cor_st",
                      corOptions = "use = 'p'")
tom <- TOMsimilarity(adj_mat, TOMType = "signed")
diss_tom <- 1 - tom
gene_tree <- hclust(as.dist(diss_tom), method = "average")
min_module_size <- 30
module_labels <- cutreeDynamic(dendro = gene_tree, distM = diss_tom,
                                deepSplit = 2, pamRespectsDendro = FALSE,
                                minClusterSize = min_module_size)
module_labels <- labels2colors(module_labels)

# Merge similar modules
me_list <- moduleEigengenes(t(expr_filt), colors = module_labels, excludeGrey = TRUE)
mes <- me_list$eigengenes
merge <- mergeCloseModules(t(expr_filt), module_labels, cutHeight = 0.25, verbose = 0)
module_labels <- merge$colors
module_colors <- module_labels
n_modules <- length(unique(module_colors))

log_msg(sprintf("Network constructed: %d modules", n_modules))

# Gene dendrogram with module colors
open_tiff(file.path(OUT_DIR, "results/figures/Figure5C_gene_dendrogram.tiff"), width = 10, height = 5)
plotDendroAndColors(gene_tree, module_colors,
                    "Module Colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main = "Gene Dendrogram and Module Colors")
dev.off()

# Module eigengenes
MEs <- moduleEigengenes(t(expr_filt), colors = module_colors, excludeGrey = TRUE)$eigengenes
rownames(MEs) <- colnames(expr_filt)

# Trait data
trait_data <- data.frame(
  Tumor = ifelse(metadata$group == "Tumor", 1, 0),
  FerroptosisScore = ferr_score[metadata$sample_id],
  row.names = metadata$sample_id
)
trait_data <- trait_data[colnames(expr_filt), , drop = FALSE]

# Module-trait relationships
module_trait_cor <- cor(MEs, trait_data, method = "spearman", use = "pairwise")
module_trait_pval <- matrix(NA, nrow = ncol(MEs), ncol = ncol(trait_data))
rownames(module_trait_pval) <- colnames(MEs)
colnames(module_trait_pval) <- colnames(trait_data)

for (i in 1:ncol(MEs)) {
  for (j in 1:ncol(trait_data)) {
    if (length(unique(trait_data[[j]])) > 2) {
      module_trait_pval[i, j] <- cor.test(MEs[[i]], trait_data[[j]], method = "spearman")$p.value
    } else {
      module_trait_pval[i, j] <- wilcox.test(MEs[[i]] ~ trait_data[[j]])$p.value
    }
  }
}

# Module-trait heatmap
open_tiff(file.path(OUT_DIR, "results/figures/Figure5D_module_trait_heatmap.tiff"), width = 8, height = 6)
labeledHeatmap(Matrix = module_trait_cor, xLabels = colnames(trait_data),
               yLabels = gsub("^ME", "", names(MEs)), ySymbols = gsub("^ME", "", names(MEs)),
               colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = sprintf("%.2f\n(p=%.3f)", module_trait_cor, module_trait_pval),
               cex.text = 0.5, cex.lab = 0.7, setStdMargins = FALSE,
               main = "Module-Trait Relationships")
dev.off()

# Module eigengene adjacency heatmap
ME_cor <- cor(MEs, method = "spearman")
open_tiff(file.path(OUT_DIR, "results/figures/Figure5E_ME_network.tiff"), width = 7, height = 6)
pheatmap::pheatmap(ME_cor, display_numbers = TRUE, fontsize = 7,
                   main = "Module Eigengene Adjacency",
                   color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

# Gene significance and module membership
tumor_status <- trait_data$Tumor
GS_tumor <- as.numeric(apply(expr_filt, 1, function(x) {
  if (var(x) == 0) return(0)
  -log10(wilcox.test(x ~ tumor_status)$p.value)
}))
GS_ferroptosis <- as.numeric(apply(expr_filt, 1, function(x) {
  if (var(x) == 0) return(0)
  cor(x, trait_data$FerroptosisScore, method = "spearman", use = "pairwise")
}))

# Module membership
MM <- as.data.frame(cor(t(expr_filt), MEs, method = "spearman", use = "pairwise"))

# Find modules with highest ferroptosis association
ferr_module_assoc <- module_trait_cor[, "FerroptosisScore", drop = FALSE]
ferr_module_assoc <- cbind(ferr_module_assoc, Pval = module_trait_pval[, "FerroptosisScore"])
ferr_module_assoc <- ferr_module_assoc[order(abs(ferr_module_assoc[, 1]), decreasing = TRUE), ]
write.csv(ferr_module_assoc, file.path(OUT_DIR, "results/WGCNA/module_ferroptosis_association.csv"))

# Top modules associated with both tumor and ferroptosis
top_modules <- intersect(
  rownames(module_trait_cor)[order(abs(module_trait_cor[, "Tumor"]), decreasing = TRUE)][1:5],
  rownames(module_trait_cor)[order(abs(module_trait_cor[, "FerroptosisScore"]), decreasing = TRUE)][1:5]
)
log_msg(sprintf("Top modules associated with both Tumor and Ferroptosis: %s",
                paste(top_modules, collapse = ", ")))

# GS vs MM plots for top modules
for (mod in top_modules[1:min(4, length(top_modules))]) {
  mod_idx <- module_colors == mod
  mod_genes <- rownames(expr_filt)[mod_idx]
  if (length(mod_genes) < 3) next

  plot_data <- data.frame(
    MM = MM[mod_genes, mod, drop = TRUE],
    GS_tumor = GS_tumor[mod_genes],
    GS_ferr = GS_ferroptosis[mod_genes]
  )

  p_mm_gs <- ggplot(plot_data, aes(MM, GS_tumor)) +
    geom_point(alpha = 0.5, color = mod) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    labs(title = paste("Module", mod, "- Module Membership vs Gene Significance"),
         x = "Module Membership", y = "Gene Significance (Tumor)") +
    theme_minimal()
  save_tiff_600dpi(p_mm_gs, file.path(OUT_DIR, sprintf("results/figures/Figure5F_MM_vs_GS_%s.tiff", mod)),
                   width = 6, height = 5)
}

# Intramodular connectivity
adj <- adjacency(t(expr_filt), power = softPower, type = "signed")
intramodular_conn <- data.frame(Gene = rownames(expr_filt), stringsAsFactors = FALSE)
for (mod in unique(module_colors)) {
  mod_idx <- which(module_colors == mod)
  if (length(mod_idx) > 1) {
    conn <- apply(adj[mod_idx, mod_idx, drop = FALSE], 2, sum, na.rm = TRUE) / (length(mod_idx) - 1)
    conn_vec <- rep(NA, nrow(expr_filt))
    conn_vec[mod_idx] <- conn
    intramodular_conn[[mod]] <- conn_vec
  }
}
write.csv(intramodular_conn, file.path(OUT_DIR, "results/WGCNA/intramodular_connectivity.csv"), row.names = FALSE)

# Save all WGCNA results
module_info <- data.frame(
  Gene = rownames(expr_filt),
  Module = module_colors,
  ModuleLabel = module_labels,
  GS_Tumor = GS_tumor,
  GS_Ferroptosis = GS_ferroptosis,
  stringsAsFactors = FALSE
)
module_info <- cbind(module_info, MM)

write.csv(module_info, file.path(OUT_DIR, "results/WGCNA/module_assignments.csv"), row.names = FALSE)
saveRDS(list(tree = gene_tree, diss_tom = diss_tom, adj = adj_mat),
        file.path(OUT_DIR, "results/WGCNA/WGCNA_network.rds"))
saveRDS(MEs, file.path(OUT_DIR, "results/WGCNA/module_eigengenes.rds"))
saveRDS(module_trait_cor, file.path(OUT_DIR, "results/WGCNA/module_trait_correlations.rds"))
write.csv(as.data.frame(module_trait_cor), file.path(OUT_DIR, "results/WGCNA/module_trait_correlations.csv"))
write.csv(as.data.frame(module_trait_pval), file.path(OUT_DIR, "results/WGCNA/module_trait_pvalues.csv"))

log_msg("Step 8 complete: WGCNA analysis finished")
