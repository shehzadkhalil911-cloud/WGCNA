# ============================================================
# 12 - SHAP analysis
# SHAP values for the final model on holdout samples (shapviz).
# Representative-sample waterfalls via a median-|SHAP|-distance
# rule; beeswarm, top-feature dependence, and SHAP heatmap.
# ============================================================
PROJECT_DIR <- "F:/independent projects/cancer/breast cancer/Breast cancer 8"
source(file.path(PROJECT_DIR, "config/config.R"))
source(file.path(PROJECT_DIR, "R/utils.R"))
suppressMessages({library(dplyr); library(ggplot2); library(shapviz)})

logmsg <- function(...) cat(sprintf(...), "\n")

final <- readRDS(file.path(DIR_MODELS, "final_model.rds"))
split <- readRDS(file.path(DIR_MODELS, "primary_split.rds"))
disc  <- load_discovery()

expr_ho <- disc$expr[final$genes, split$hold_samples, drop = FALSE]
meta_ho <- disc$meta[match(split$hold_samples, disc$meta$sample_id), ]
truth_ho <- meta_ho$class

X <- as.matrix(t(expr_ho))
stopifnot(identical(colnames(X), final$genes))

set.seed(SEED)
shp <- shapviz(final$model, X_pred = X)
logmsg("SHAP computed for %d samples x %d features", nrow(X), ncol(X))

# ---------------- aggregate SHAP ----------------
S <- get_shap_values(shp)
baseline <- get_baseline(shp)
mean_shap <- colMeans(abs(S))

shap_summary <- data.frame(
  gene = names(sort(mean_shap, decreasing = TRUE)),
  mean_abs_shap = sort(mean_shap, decreasing = TRUE),
  rank = seq_along(mean_shap)
)
write.csv(shap_summary, file.path(DIR_TABLES, "Table_19_SHAP_Feature_Contributions.csv"),
          row.names = FALSE)

# per-class summary (top 40)
top40 <- shap_summary$gene[1:40]
per_class <- do.call(rbind, lapply(c("Tumor", "Normal"), function(cl) {
  idx <- which(truth_ho == cl)
  data.frame(class = cl, gene = top40,
             mean_shap = colMeans(S[idx, top40, drop = FALSE]),
             mean_abs_shap = colMeans(abs(S[idx, top40, drop = FALSE])))
}))
write.csv(per_class, file.path(DIR_TABLES, "Table_20_SHAP_Per_Class.csv"), row.names = FALSE)

# ---------------- representative samples ----------------
med_profile <- colMeans(abs(S))
dist_to_med <- apply(S, 1, function(s) mean(abs(abs(s) - med_profile)))
rep_idx <- order(dist_to_med)
# pick 2 representative samples per class closest to median profile
rep_tumor <- which(truth_ho == "Tumor")[order(dist_to_med[truth_ho == "Tumor"])[1:2]]
rep_normal <- which(truth_ho == "Normal")[order(dist_to_med[truth_ho == "Normal"])[1:2]]
rep_idx <- c(rep_tumor, rep_normal)
rep_names <- rownames(X)[rep_idx]
logmsg("Representative samples: %s", paste(rep_names, collapse = ", "))

# ---------------- figures (Figure_11) ----------------
# beeswarm (top 30)
p_swarm <- sv_importance(shp, kind = "beeswarm", max_display = 30) +
  theme_bc(base_size = 9) +
  labs(title = "SHAP beeswarm (holdout)") +
  theme(legend.position = "right")
save_tiff_600dpi(p_swarm, "Figure_11_SHAP_Beeswarm.tiff", width = 8, height = 7)

# mean |SHAP| bar (top 30)
p_bar <- sv_importance(shp, kind = "bar", max_display = 30) +
  theme_bc(base_size = 9) +
  labs(title = "Mean |SHAP| feature importance (holdout)")
save_tiff_600dpi(p_bar, "Figure_11_SHAP_Importance.tiff", width = 6.5, height = 7)

# dependence plots for top 3
for (i in 1:3) {
  g <- top40[i]
  p_dep <- sv_dependence(shp, v = g, color_var = "auto") +
    theme_bc(base_size = 9) +
    labs(title = sprintf("SHAP dependence: %s", g))
  save_tiff_600dpi(p_dep, sprintf("Figure_11_SHAP_Dependence_%02d_%s.tiff", i, g),
                   width = 6.5, height = 5)
}

# SHAP heatmap (top 40, holdout samples)
open_tiff_base("Figure_11_SHAP_Heatmap.tiff", width = 8, height = 8)
shp40 <- S[, top40]
heat_data <- t(shp40)
heat_data[!is.finite(heat_data)] <- 0
ann_col <- data.frame(Class = truth_ho, row.names = rownames(X))
ComplexHeatmap::Heatmap(
  heat_data, name = "SHAP value",
  col = circlize::colorRamp2(c(-1, 0, 1) * max(abs(heat_data)) / 3,
                             c("#2166AC", "#F7F7F7", "#B2182B")),
  top_annotation = ComplexHeatmap::HeatmapAnnotation(
    df = ann_col, col = list(Class = c(Tumor = COL_TUMOR, Normal = COL_NORMAL))),
  show_column_names = FALSE, row_names_gp = grid::gpar(fontsize = 6),
  column_split = factor(truth_ho, levels = CLASS_LEVELS),
  use_raster = TRUE, raster_quality = 5,
  heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 9), labels_gp = grid::gpar(fontsize = 8)))
dev.off()

# waterfall for representative samples
for (i in seq_along(rep_idx)) {
  s <- rep_idx[i]
  p_w <- sv_waterfall(shp, row_id = s, max_display = 15) +
    theme_bc(base_size = 9) +
    labs(title = sprintf("SHAP waterfall: %s (%s)", rownames(X)[s], truth_ho[s]))
  save_tiff_600dpi(p_w, sprintf("Figure_11_SHAP_Waterfall_%02d.tiff", i), width = 7.5, height = 6)
}

saveRDS(list(shap = S, summary = shap_summary, per_class = per_class,
             representative = rep_names),
        file.path(DIR_MODELS, "shap_results.rds"))
logmsg("SHAP analysis complete.")
