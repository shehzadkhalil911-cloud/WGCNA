# 06_batch_correction_and_PCA.R
source("config/analysis_config.R")
source("R/load_functions.R")
ensure_paths()
suppressMessages({ library(ggplot2); library(ggpubr); library(pheatmap); library(cowplot); library(uwot) })

disc <- readRDS(file.path(paths$discovery, "discovery_cohort_raw.rds"))
expr <- disc$expression; meta <- disc$metadata

## ComBat preserving tumor/normal
combat <- combat_correct(expr, meta, batch_var = "dataset", group_var = "class",
                          par.prior = TRUE)
combat <- as.matrix(combat)
saveRDS(combat, file.path(paths$discovery, "discovery_combat.rds"))
write.csv(as.data.frame(combat), file.path(paths$discovery, "discovery_combat.csv"))

## batch assessment
bass <- batch_assessment(expr, combat, meta, "dataset", "class")
write.csv(as.data.frame(bass), file.path(paths$res_pca, "batch_assessment.csv"), row.names = FALSE)

## PCA before/after coloured by dataset and group
pca_fun <- function(mat, colour_by, shape_by = NULL, title) {
  pca <- prcomp(t(mat), scale. = FALSE)
  pc <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
  pc$col <- meta[[colour_by]]
  if (!is.null(shape_by)) pc$shp <- meta[[shape_by]]
  p <- ggplot(pc, aes(x = PC1, y = PC2, colour = col)) +
    geom_point(size = 2) + labs(title = title, colour = colour_by) + pub_theme()
  if (!is.null(shape_by)) p <- p + aes(shape = shp)
  p
}
p_ds_before <- pca_fun(expr, "dataset", NULL, "PCA before correction (by dataset)")
p_ds_after  <- pca_fun(combat, "dataset", NULL, "PCA after correction (by dataset)")
p_gr_before <- pca_fun(expr, "class", NULL, "PCA before correction (by group)")
p_gr_after  <- pca_fun(combat, "class", NULL, "PCA after correction (by group)")

## before/after correlation heatmaps
cor_before <- cor(expr); cor_before[!is.finite(cor_before)] <- 0
cor_after  <- cor(combat); cor_after[!is.finite(cor_after)] <- 0
ord <- order(meta$dataset)
pheatmap::pheatmap(cor_before[ord, ord], main = "Sample correlation before correction",
                   filename = file.path(paths$res_pca, "corr_before.png"),
                   show_rownames = FALSE, show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE)
pheatmap::pheatmap(cor_after[ord, ord], main = "Sample correlation after correction",
                   filename = file.path(paths$res_pca, "corr_after.png"),
                   show_rownames = FALSE, show_colnames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE)

## UMAP (supplementary)
set.seed(GLOBAL_SEED)
um <- uwot::umap(t(combat), n_neighbors = 15, seed = GLOBAL_SEED)
ump <- data.frame(UMAP1 = um[,1], UMAP2 = um[,2], group = meta$class)
p_umap <- ggplot(ump, aes(x = UMAP1, y = UMAP2, colour = group)) + geom_point() +
  labs(title = "UMAP after correction") + pub_theme()

## combined PCA comparison
p_combined <- cowplot::plot_grid(p_ds_before, p_ds_after, p_gr_before, p_gr_after,
                                 labels = c("A","B","C","D"), nrow = 2)
save_publication_plot(p_combined, file.path(paths$res_pca, "PCA_combined"), width = 10, height = 8)
save_publication_plot(p_ds_before, file.path(paths$res_pca, "PCA_before_dataset"))
save_publication_plot(p_ds_after, file.path(paths$res_pca, "PCA_after_dataset"))
save_publication_plot(p_gr_before, file.path(paths$res_pca, "PCA_before_group"))
save_publication_plot(p_gr_after, file.path(paths$res_pca, "PCA_after_group"))
save_publication_plot(p_umap, file.path(paths$res_pca, "UMAP_after"))

cat("Batch correction + PCA complete.\n")
cat("Batch var explained before:", round(bass$batch_var_before,3), "after:", round(bass$batch_var_after,3), "\n")
cat("Group var explained before:", round(bass$group_var_before,3), "after:", round(bass$group_var_after,3), "\n")
quit(save = "no", status = 0, runLast = FALSE)
