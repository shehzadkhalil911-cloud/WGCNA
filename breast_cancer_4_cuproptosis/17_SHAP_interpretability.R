# 17_SHAP_interpretability.R
source("config/analysis_config.R")
source("R/load_functions.R")
ensure_paths()
suppressMessages({ library(shapviz); library(kernelshap); library(ranger); library(ggplot2); library(pheatmap) })

pm <- readRDS(file.path(paths$res_model, "predictive_models.rds"))
combat <- readRDS(file.path(paths$discovery, "discovery_combat.rds"))
meta <- readRDS(file.path(paths$discovery, "discovery_cohort_raw.rds"))$metadata
final_genes <- pm$final_genes
best <- pm$best_model
y_test <- pm$y_test; te <- pm$test_idx
x_test <- t(combat[final_genes, te, drop = FALSE])

## best model object -> compute SHAP via kernelshap then wrap with shapviz
if (best == "RF") {
  fit <- pm$fits$rf
  ks <- kernelshap(fit, X = as.data.frame(x_test))
  sv <- shapviz(ks)
} else if (best == "LASSO" || best == "Logistic") {
  fit <- pm$fits$cvfit
  ks <- kernelshap(fit, X = as.matrix(x_test), bg_X = as.matrix(x_test))
  sv <- shapviz(ks)
} else {
  fit <- pm$fits$svm
  ks <- kernelshap(fit, X = as.data.frame(x_test))
  sv <- shapviz(ks)
}

shap_vals <- get_shap_values(sv)
## for probability models, shap_vals is a list per class; take the 'Tumor' class
if (is.list(shap_vals)) shap_mat <- as.matrix(shap_vals[["Tumor"]]) else shap_mat <- as.matrix(shap_vals)
write.csv(shap_mat, file.path(paths$res_shap, "SHAP_values.csv"), row.names = FALSE)
saveRDS(sv, file.path(paths$res_shap, "SHAP_object.rds"))

## bar plot (mean abs importance)
p_bar <- sv_importance(sv, kind = "bar", max_vars = length(final_genes))
save_publication_plot(p_bar, file.path(paths$res_shap, "SHAP_bar"), width = 7, height = 5)
## beeswarm
p_bee <- sv_importance(sv, kind = "beeswarm", max_vars = length(final_genes), seed = GLOBAL_SEED)
save_publication_plot(p_bee, file.path(paths$res_shap, "SHAP_beeswarm"), width = 7, height = 5)
## dependence plots for each gene
for (g in final_genes) {
  p <- sv_dependence(sv, v = g)
  save_publication_plot(p, file.path(paths$res_shap, paste0("SHAP_dependence_", g)), width = 6, height = 5)
}
## waterfall for representative tumor (highest prob) and normal (lowest prob)
probs <- pm$preds[[best]]
ord <- order(probs, decreasing = TRUE)
p_wt <- sv_waterfall(sv, row_id = ord[1])
save_publication_plot(p_wt, file.path(paths$res_shap, "SHAP_waterfall_tumor"), width = 7, height = 5)
p_wn <- sv_waterfall(sv, row_id = ord[length(ord)])
save_publication_plot(p_wn, file.path(paths$res_shap, "SHAP_waterfall_normal"), width = 7, height = 5)

## sample x feature SHAP heatmap
pheatmap::pheatmap(shap_mat, main = "SHAP values (test set)",
                   filename = file.path(paths$res_shap, "SHAP_heatmap.png"), show_rownames = FALSE)

## explanation table
exp_tab <- data.frame(
  sample_id = colnames(x_test),
  predicted_prob = probs,
  observed = as.character(y_test),
  stringsAsFactors = FALSE
)
write.csv(exp_tab, file.path(paths$res_shap, "SHAP_prediction_explanation.csv"), row.names = FALSE)

cat("SHAP analysis complete (best model:", best, ").\n")
cat(" Representative tumor sample idx:", ord[1], " normal sample idx:", ord[length(ord)], "\n")
quit(save = "no", status = 0, runLast = FALSE)
