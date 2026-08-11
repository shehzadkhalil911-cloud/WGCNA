# =============================================================================
# Script 06: WGCNA Network Analysis
# Breast Cancer 7 Cell Death Crosstalk
# =============================================================================

source(file.path("..", "config", "analysis_config.R"))

library(WGCNA)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

log_message("Starting WGCNA analysis.")

# --- 1. LOAD DATA ---
combat_matrix <- readRDS(file.path(OUTPUT$batch_correction, "combat_corrected_matrix.rds"))
discovery_metadata <- readRDS(file.path(OUTPUT$batch_correction, "discovery_metadata.rds"))
scores <- readRDS(file.path(OUTPUT$activity, "program_scores.rds"))

# --- 2. FILTER GENES FOR WGCNA ---
# Use most variable genes
gene_var <- apply(combat_matrix, 1, var, na.rm = TRUE)
gene_mean <- apply(combat_matrix, 1, mean, na.rm = TRUE)
gene_cv <- gene_var / (gene_mean + 0.01)

# Select top 5000 most variable genes
n_top_genes <- 5000
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(n_top_genes, length(gene_var))]
expr_data <- t(combat_matrix[top_genes, ])

log_message(sprintf("WGCNA input: %d genes x %d samples", ncol(expr_data), nrow(expr_data)))

# --- 3. SAMPLE DENDOGRAM ---
sample_tree <- hclust(dist(expr_data), method = "average")

open_tiff(file.path(OUTPUT$wgcna, "sample_dendrogram.tiff"), width = 10, height = 5)
par(mar = c(4, 4, 2, 2))
plot(sample_tree, main = "Sample Dendrogram (Discovery)", 
     xlab = "", sub = "", cex = 0.6, labels = FALSE)
dev.off()

# --- 4. SOFT-THRESHOLD SELECTION ---
powers <- c(1:10, seq(12, 20, by = 2))
sft <- pickSoftThreshold(expr_data, powerVector = powers, verbose = 0, networkType = "signed")

# Scale independence plot
sft_df <- data.frame(
  Power = sft$fitIndices[, 1],
  SFT_Rsq = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  MeanK = sft$fitIndices[, 5]
)

p_scale <- ggplot(sft_df, aes(x = Power, y = SFT_Rsq)) +
  geom_point(size = 3, color = "#0072B2") +
  geom_text(aes(label = Power), vjust = -0.8, size = 3) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_hline(yintercept = 0.9, linetype = "dashed", color = "orange", linewidth = 0.5) +
  labs(title = "Scale Independence", x = "Soft Threshold (Power)", y = "SFT R-squared") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_tiff_600dpi(p_scale, file.path(OUTPUT$wgcna, "soft_threshold_scale.tiff"))

p_connect <- ggplot(sft_df, aes(x = Power, y = MeanK)) +
  geom_point(size = 3, color = "#D55E00") +
  geom_text(aes(label = Power), vjust = -0.8, size = 3) +
  labs(title = "Mean Connectivity", x = "Soft Threshold (Power)", y = "Mean Connectivity") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_tiff_600dpi(p_connect, file.path(OUTPUT$wgcna, "soft_threshold_connectivity.tiff"))

# Select soft power
softPower <- sft$powerEstimate
if (is.na(softPower) || softPower < 1) softPower <- 6
log_message(sprintf("Selected soft power: %d", softPower))

# --- 5. NETWORK CONSTRUCTION ---
net <- blockwiseModules(expr_data,
                        power = softPower,
                        TOMType = "signed",
                        minModuleSize = 30,
                        maxBlockSize = 5000,
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = file.path(OUTPUT$wgcna, "TOM"),
                        verbose = 0)

module_labels <- net$colors
module_colors <- labels2colors(module_labels)
module_names <- paste0("ME", module_labels)

log_message(sprintf("Number of modules: %d", length(unique(module_labels))))
log_message(sprintf("Module sizes: %s", 
                    paste(sort(table(module_labels), decreasing = TRUE), collapse = ", ")))

# Save module assignments
module_df <- data.frame(
  Gene = top_genes,
  Module = module_labels,
  ModuleColor = module_colors,
  stringsAsFactors = FALSE
)
write.csv(module_df, file.path(OUTPUT$wgcna, "module_assignments.csv"), row.names = FALSE)
saveRDS(net, file.path(OUTPUT$wgcna, "WGCNA_network.rds"))

# --- 6. GENE DENDOGRAM WITH MODULE COLORS ---
open_tiff(file.path(OUTPUT$wgcna, "gene_dendrogram.tiff"), width = 12, height = 6)
plotDendroAndColors(net$dendrograms[[1]], 
                    colors = cbind(module_colors[net$blockGenes[[1]]]),
                    groupLabels = "Modules",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene Dendrogram with Module Colors")
dev.off()

# --- 7. MODULE EIGENGENES ---
MEs <- moduleEigengenes(expr_data, colors = module_colors)$eigengenes
colnames(MEs) <- gsub("ME", "", colnames(MEs))
# Reorder modules (handle non-numeric names like "grey")
mod_order <- order(suppressWarnings(as.numeric(gsub("ME", "", colnames(MEs)))), na.last = TRUE)
MEs <- MEs[, mod_order]

saveRDS(MEs, file.path(OUTPUT$wgcna, "module_eigengenes.rds"))
write.csv(MEs, file.path(OUTPUT$wgcna, "module_eigengenes.csv"))

# --- 8. MODULE-TRAIT CORRELATIONS ---
traits <- data.frame(
  Group = as.numeric(discovery_metadata$group) - 1,  # 0=Normal, 1=Tumor
  Ferroptosis = scores["Ferroptosis", ],
  Necroptosis = scores["Necroptosis", ],
  Cuproptosis = scores["Cuproptosis", ],
  stringsAsFactors = FALSE
)
rownames(traits) <- colnames(scores)

# Match sample order
traits <- traits[rownames(MEs), ]

module_trait_cor <- cor(MEs, traits, method = "spearman", use = "p")
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(MEs))

write.csv(module_trait_cor, file.path(OUTPUT$wgcna, "module_trait_correlations.csv"))
write.csv(module_trait_pval, file.path(OUTPUT$wgcna, "module_trait_pvalues.csv"))
saveRDS(list(cor = module_trait_cor, pval = module_trait_pval), 
        file.path(OUTPUT$wgcna, "module_trait_correlations.rds"))

# Module-trait heatmap
p_cor_text <- ifelse(module_trait_pval < 0.001, "***",
              ifelse(module_trait_pval < 0.01, "**",
              ifelse(module_trait_pval < 0.05, "*", "")))

ht_mt <- Heatmap(module_trait_cor,
                 name = "Correlation",
                 col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(sprintf("%.2f%s", module_trait_cor[i, j], p_cor_text[i, j]), 
                             x, y, gp = gpar(fontsize = 8))
                 },
                 rect_gp = gpar(col = "grey60", lwd = 0.5),
                 row_title = "Module Eigengene",
                 column_title = "Module-Trait Correlations",
                 row_names_gp = gpar(fontsize = 8),
                 column_names_gp = gpar(fontsize = 10))

open_tiff(file.path(OUTPUT$wgcna, "module_trait_heatmap.tiff"), width = 8, height = 10)
draw(ht_mt)
dev.off()

# --- 9. MODULE-MODULE ADJACENCY ---
ME_adj <- cor(MEs, method = "spearman")

ht_me <- Heatmap(ME_adj,
                 name = "Correlation",
                 col = colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B")),
                 rect_gp = gpar(col = "grey60", lwd = 0.5),
                 column_title = "Module Eigengene Adjacency",
                 row_names_gp = gpar(fontsize = 8),
                 column_names_gp = gpar(fontsize = 8))

open_tiff(file.path(OUTPUT$wgcna, "module_eigengene_adjacency.tiff"), width = 8, height = 7)
draw(ht_me)
dev.off()

# --- 10. MODULE SIZE BAR PLOT ---
module_size <- data.frame(
  Module = names(table(module_labels)),
  Size = as.numeric(table(module_labels)),
  stringsAsFactors = FALSE
)

p_mod_size <- ggplot(module_size, aes(x = reorder(Module, -Size), y = Size, fill = as.factor(Size))) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_viridis_d() +
  labs(title = "Module Sizes", x = "Module", y = "Number of Genes") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
save_tiff_600dpi(p_mod_size, file.path(OUTPUT$wgcna, "module_sizes.tiff"))

# --- 11. GENE SIGNIFICANCE vs MODULE MEMBERSHIP ---
gs_group <- as.numeric(discovery_metadata$group) - 1
gene_significance <- cor(expr_data, gs_group, method = "spearman", use = "p")
gene_significance <- as.numeric(gene_significance)

module_membership <- cor(expr_data, MEs, method = "spearman", use = "p")

# Per module: GS vs MM plots
for (mod in unique(module_labels)) {
  mod_name <- paste0("ME", mod)
  if (mod_name %in% colnames(module_membership)) {
    in_module <- module_labels == mod
    mm_values <- module_membership[, mod_name]
    
    plot_df <- data.frame(
      MM = mm_values,
      GS = gene_significance,
      InModule = in_module
    )
    
    cor_val <- cor(mm_values[in_module], gene_significance[in_module], method = "spearman")
    
    p_gs_mm <- ggplot(plot_df, aes(x = MM, y = GS, color = InModule)) +
      geom_point(size = 0.5, alpha = 0.5) +
      scale_color_manual(values = c("grey70", "red3")) +
      annotate("text", x = min(mm_values), y = max(gene_significance), 
               label = sprintf("r = %.3f", cor_val), hjust = 0, size = 4) +
      labs(title = sprintf("Module %s: GS vs MM", mod_name),
           x = "Module Membership", y = "Gene Significance") +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "none")
    save_tiff_600dpi(p_gs_mm, file.path(OUTPUT$wgcna, sprintf("GS_vs_MM_module_%s.tiff", mod_name)))
  }
}

# --- 12. INTRAMODULAR CONNECTIVITY ---
intramodular_connectivity <- intramodularConnectivity(adjacency(expr_data, power = softPower, type = "signed"), 
                                                      module_colors)
intramodular_connectivity$Gene <- top_genes
intramodular_connectivity$Module <- module_labels
intramodular_connectivity$ModuleColor <- module_colors
write.csv(intramodular_connectivity, file.path(OUTPUT$wgcna, "intramodular_connectivity.csv"))

# --- 13. TOM HEATMAP (for a subset of genes) ---
if (file.exists(file.path(OUTPUT$wgcna, "TOM-block.1.RData"))) {
  tryCatch({
    tom_env <- new.env()
    load(file.path(OUTPUT$wgcna, "TOM-block.1.RData"), envir = tom_env)
    if (exists("TOM", where = tom_env)) {
      TOM <- get("TOM", envir = tom_env)
      # Subset for visualization
      top_modules <- sort(unique(module_labels))[1:5]
      top_mod_genes <- which(module_labels %in% top_modules)
      if (length(top_mod_genes) > 10) {
        set.seed(123)
        n_sel <- min(200, length(top_mod_genes))
        selected <- sample(top_mod_genes, n_sel)
        selected <- sort(selected)
        
        tom_subset <- TOM[selected, selected]
        if (nrow(tom_subset) > 2 && ncol(tom_subset) > 2) {
          rownames(tom_subset) <- top_genes[selected]
          colnames(tom_subset) <- top_genes[selected]
          
          ht_tom <- Heatmap(tom_subset,
                            name = "TOM",
                            col = colorRamp2(c(0, 0.5, 1), c("white", "orange", "red")),
                            show_row_names = FALSE,
                            show_column_names = FALSE,
                            column_title = "Topological Overlap Matrix")
          
          open_tiff(file.path(OUTPUT$wgcna, "TOM_heatmap.tiff"), width = 8, height = 7)
          draw(ht_tom)
          dev.off()
        }
      }
    }
  }, error = function(e) log_message(sprintf("TOM heatmap error: %s", e$message), "WARNING"))
}

log_message("WGCNA analysis complete.")
