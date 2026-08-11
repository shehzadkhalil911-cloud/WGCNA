# 22_figures.R
source("config/analysis_config.R")
source("R/load_functions.R")
suppressMessages({ library(ggplot2); library(ggpubr); library(patchwork); library(ComplexHeatmap); library(circlize) })

## ---- Figure 1: Discovery workflow summary composites ----
## Panel A: ComBat batch correction PCA (use existing pca figure if present)
pca_before <- list.files(paths$res_pca, pattern = "pca_before", full.names = TRUE)[1]
pca_after <- list.files(paths$res_pca, pattern = "pca_after", full.names = TRUE)[1]

## Panel B: DEG volcano
deg <- readRDS(file.path(paths$res_de, "discovery_DEG_full.rds"))
if (is.list(deg) && !is.null(deg$table)) degt <- deg$table else degt <- deg
if (!"logFC" %in% colnames(degt) && "log2FC" %in% colnames(degt)) degt$logFC <- degt$log2FC
degt$sig <- with(degt, adj.P.Val < 0.05 & abs(logFC) >= 1)
vol <- ggplot(degt, aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
  geom_point(size = 0.6, alpha = 0.6) + scale_color_manual(values = c("grey", "#E15759")) +
  labs(title = "Discovery DEGs (842)", x = "log2 FC", y = "-log10 adj.P") + pub_theme() +
  theme(legend.position = "none")
save_publication_plot(vol, file.path(paths$res_fig, "fig1B_volcano"), width = 6, height = 5)

## Panel C: consensus evidence heatmap
cons <- readRDS(file.path(paths$res_ml, "consensus", "consensus_genes.rds"))
ev <- cons$evidence
evm <- as.matrix(ev[, c("DEG_status", "CRG_status", "LASSO", "SVM", "RF")])
evm[, "DEG_status"] <- as.numeric(evm[, "DEG_status"])
col_fun <- colorRamp2(c(0, 1), c("white", "steelblue"))
grDevices::png(file.path(paths$res_fig, "fig1C_consensus.png"), width = 900, height = 700, res = 130)
ht <- ComplexHeatmap::Heatmap(evm, name = "flag", col = col_fun,
  show_row_names = TRUE, show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 9), column_names_gp = grid::gpar(fontsize = 10))
ComplexHeatmap::draw(ht)
grDevices::dev.off()
grDevices::tiff(file.path(paths$res_fig, "fig1C_consensus.tiff"), width = 7, height = 5,
                units = "in", res = 600, compression = "lzw")
ComplexHeatmap::draw(ht)
grDevices::dev.off()

## Panel D: model comparison ROC (discovery)
if (file.exists(file.path(paths$res_model, "combined_roc.pdf"))) {
  file.copy(file.path(paths$res_model, "combined_roc.pdf"),
            file.path(paths$res_fig, "fig1D_model_roc.pdf"), overwrite = TRUE)
}

## ---- Figure 2: Cuproptosis crosstalk ----
cr <- readRDS(file.path(paths$res_cross, "crosstalk_correlations.rds"))
edges <- cr$decrg_edges
if (!is.null(edges) && nrow(edges) > 0) {
  top_e <- head(edges[order(abs(edges$r), decreasing = TRUE), ], 30)
  top_e$lab <- paste(top_e$from, top_e$to, sep = " - ")
  pcr <- ggplot(top_e, aes(x = reorder(lab, abs(r)), y = r, fill = r > 0)) +
    geom_bar(stat = "identity") + coord_flip() + scale_fill_manual(values = c("blue", "red")) +
    labs(title = "Top CRG-DEG correlations", x = "", y = "Spearman rho") + pub_theme() +
    theme(legend.position = "none", axis.text.y = element_text(size = 6))
  save_publication_plot(pcr, file.path(paths$res_fig, "fig2_crosstalk"), width = 7, height = 8)
}

## ---- Figure 3: Immune infiltration ----
if (file.exists(file.path(paths$res_immune, "immune_boxplots.pdf")))
  file.copy(file.path(paths$res_immune, "immune_boxplots.pdf"),
            file.path(paths$res_fig, "fig3A_immune_boxplots.pdf"), overwrite = TRUE)
if (file.exists(file.path(paths$res_immune, "immune_heatmap.png")))
  file.copy(file.path(paths$res_immune, "immune_heatmap.png"),
            file.path(paths$res_fig, "fig3B_immune_heatmap.png"), overwrite = TRUE)
ct <- read.csv(file.path(paths$res_immune, "gene_immune_crosstalk.csv"))
if (nrow(ct) > 0) {
  sig_ct <- ct[ct$adj_p < 0.05, ]
  pc <- ggplot(sig_ct, aes(x = gene, y = immune_feature, fill = rho)) +
    geom_tile() + scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    labs(title = "Gene-immune crosstalk (significant)") + pub_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 7))
  save_publication_plot(pc, file.path(paths$res_fig, "fig3C_gene_immune"), width = 10, height = 8)
}

## ---- Figure 4: Validation & SHAP ----
if (file.exists(file.path(paths$res_ext, "validation_roc.pdf")))
  file.copy(file.path(paths$res_ext, "validation_roc.pdf"),
            file.path(paths$res_fig, "fig4A_validation_roc.pdf"), overwrite = TRUE)
if (file.exists(file.path(paths$res_ext, "validation_shap_importance.pdf")))
  file.copy(file.path(paths$res_ext, "validation_shap_importance.pdf"),
            file.path(paths$res_fig, "fig4B_shap_importance.pdf"), overwrite = TRUE)

cat("Composite figures assembled in", paths$res_fig, "\n")
quit(save = "no", status = 0, runLast = FALSE)
