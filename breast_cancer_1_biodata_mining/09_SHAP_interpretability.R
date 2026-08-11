###############################################################################
# Script 09: SHAP Model Interpretability Analysis
###############################################################################

set.seed(123)

library(xgboost)
library(shapviz)
library(ggplot2)
library(tidyverse)
library(caret)
library(pROC)
library(viridis)
library(patchwork)

base_dir <- normalizePath("breast_cancer_integrated_analysis")
data_processed_dir <- file.path(base_dir, "data_processed")
fig_dir <- file.path(base_dir, "results", "figures", "SHAP")
tables_dir <- file.path(base_dir, "results", "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helper: Save 600 dpi ---
save_publication_figure <- function(plot, filename, width = 8, height = 6) {
  ggsave(file.path(fig_dir, paste0(filename, ".tiff")), plot = plot,
         width = width, height = height, dpi = 600, compression = "lzw")
  ggsave(file.path(fig_dir, paste0(filename, ".pdf")), plot = plot,
         width = width, height = height, dpi = 600)
  cat("  Saved:", filename, "\n")
}

theme_publication <- theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

# --- Load data ---
cat("Loading data...\n")
exprs_after <- readRDS(file.path(data_processed_dir, "merged_exprs_after_batch.rds"))
pheno <- readRDS(file.path(data_processed_dir, "merged_pheno_with_batch.rds"))
final_genes <- readRDS(file.path(data_processed_dir, "final_feature_genes.rds"))

cat("Final feature genes:", length(final_genes), "\n")
cat("Genes:", paste(final_genes, collapse = ", "), "\n")

# Filter to genes that exist in expression matrix
final_genes <- intersect(final_genes, rownames(exprs_after))
cat("Genes available in expression matrix:", length(final_genes), "\n")

if (length(final_genes) < 2) {
  stop("Need at least 2 final feature genes for model building.")
}

# --- Prepare data for XGBoost ---
x_data <- t(exprs_after[final_genes, ])
y_data <- factor(pheno$group, levels = c("Normal", "Tumor"))
y_binary <- ifelse(y_data == "Tumor", 1, 0)

cat("\nClass distribution:\n")
print(table(y_data))

# --- Train-test split ---
set.seed(123)
train_idx <- createDataPartition(y_data, p = 0.7, list = FALSE)
x_train <- x_data[train_idx, ]
x_test <- x_data[-train_idx, ]
y_train <- y_binary[train_idx]
y_test <- y_binary[-train_idx]
y_train_factor <- y_data[train_idx]
y_test_factor <- y_data[-train_idx]

cat("\nTraining samples:", nrow(x_train), "\n")
cat("Testing samples:", nrow(x_test), "\n")

# --- XGBoost model ---
dtrain <- xgb.DMatrix(data = as.matrix(x_train), label = y_train)
dtest <- xgb.DMatrix(data = as.matrix(x_test), label = y_test)

params <- list(
  objective = "binary:logistic",
  eval_metric = "auc",
  max_depth = 4,
  eta = 0.05,
  subsample = 0.8,
  colsample_bytree = 0.8,
  min_child_weight = 3,
  gamma = 0.1
)

set.seed(123)
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 500,
  early_stopping_rounds = 20,
  watchlist = list(train = dtrain, test = dtest),
  verbose = 0
)

cat("\nXGBoost best iteration:", xgb_model$best_iteration, "\n")
cat("Best score:", xgb_model$best_score, "\n")

# --- Predictions ---
pred_train <- predict(xgb_model, dtrain)
pred_test <- predict(xgb_model, dtest)
pred_class <- ifelse(pred_test > 0.5, "Tumor", "Normal")
pred_class <- factor(pred_class, levels = c("Normal", "Tumor"))

# --- ROC curve ---
cat("\nGenerating ROC curve...\n")
roc_obj <- roc(y_test, pred_test)
cat("AUC:", round(auc(roc_obj), 4), "\n")

roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(color = "#377eb8", linewidth = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.75, y = 0.25, 
           label = paste0("AUC = ", round(auc(roc_obj), 3)),
           size = 5, fontface = "bold") +
  labs(title = "ROC Curve: XGBoost Classifier",
       x = "1 - Specificity (False Positive Rate)",
       y = "Sensitivity (True Positive Rate)") +
  coord_fixed() +
  theme_publication

save_publication_figure(p_roc, "Figure_7A_ROC_curve", width = 7, height = 6)

# --- Model evaluation metrics ---
conf_mat <- confusionMatrix(pred_class, y_test_factor)
cat("\nConfusion Matrix:\n")
print(conf_mat$table)

metrics <- data.frame(
  Metric = c("Accuracy", "Sensitivity", "Specificity", "Precision", 
             "Recall", "F1 Score", "AUC"),
  Value = c(
    conf_mat$overall["Accuracy"],
    conf_mat$byClass["Sensitivity"],
    conf_mat$byClass["Specificity"],
    conf_mat$byClass["Precision"],
    conf_mat$byClass["Recall"],
    conf_mat$byClass["F1"],
    auc(roc_obj)
  ),
  stringsAsFactors = FALSE
)
cat("\nModel Performance Metrics:\n")
print(metrics)
write.csv(metrics, file.path(tables_dir, "XGBoost_model_metrics.csv"), row.names = FALSE)

# --- SHAP Analysis ---
cat("\nComputing SHAP values...\n")

# Prepare data for shapviz
x_train_df <- as.data.frame(x_train)
x_test_df <- as.data.frame(x_test)

# Create shapviz object
shp <- shapviz(xgb_model, X_pred = as.matrix(x_test_df), X = x_test_df)

# A. SHAP bar plot (global importance)
cat("Generating SHAP bar plot...\n")
p_bar <- sv_importance(shp, kind = "bar") +
  labs(title = "SHAP Feature Importance (Bar)") +
  theme_publication +
  scale_fill_viridis_d()
save_publication_figure(p_bar, "Figure_7B_SHAP_barplot", width = 8, height = 6)

# B. SHAP beeswarm plot
cat("Generating SHAP beeswarm plot...\n")
p_beeswarm <- sv_importance(shp, kind = "bee") +
  labs(title = "SHAP Feature Importance (Beeswarm)") +
  theme_publication +
  scale_color_viridis(option = "D")
save_publication_figure(p_beeswarm, "Figure_7C_SHAP_beeswarm", width = 8, height = 6)

# C. SHAP dependence plots for top genes
cat("Generating SHAP dependence plots...\n")
top_genes <- colnames(x_test_df)[1:min(6, ncol(x_test_df))]

p_dep_list <- list()
for (i in seq_along(top_genes)) {
  p_dep <- sv_dependence(shp, v = top_genes[i], color = "auto") +
    labs(title = paste0("SHAP Dependence: ", top_genes[i])) +
    theme_publication
  p_dep_list[[i]] <- p_dep
}

# Combine dependence plots
if (length(p_dep_list) >= 2) {
  n_row <- ceiling(length(p_dep_list) / 2)
  combined_dep <- wrap_plots(p_dep_list, ncol = 2) +
    plot_annotation(title = "SHAP Dependence Plots for Top Feature Genes",
                    theme = theme(plot.title = element_text(size = 14, face = "bold",
                                                            hjust = 0.5)))
  save_publication_figure(combined_dep, "Figure_7D_SHAP_scatter", 
                          width = 12, height = 4 * n_row)
}

# D. SHAP waterfall for one representative sample
cat("Generating SHAP waterfall plot...\n")
tryCatch({
  p_waterfall <- sv_waterfall(shp, row_id = 1) +
    labs(title = "SHAP Waterfall Plot: Representative Sample") +
    theme_publication
  save_publication_figure(p_waterfall, "Figure_7E_SHAP_waterfall", width = 8, height = 6)
}, error = function(e) {
  cat("  Waterfall plot error:", e$message, "\n")
})

# E. SHAP force plot
cat("Generating SHAP force plot...\n")
tryCatch({
  # Save static force plot as image
  p_force <- sv_force(shp, row_id = 1) +
    labs(title = "SHAP Force Plot: Representative Sample") +
    theme_publication
  save_publication_figure(p_force, "Figure_7F_SHAP_force_plot", width = 10, height = 4)
  
  # Also create a combined force plot for multiple samples
  if (nrow(x_test_df) >= 4) {
    p_force_multi <- sv_force(shp, row_id = 1:4) +
      labs(title = "SHAP Force Plot: Multiple Samples") +
      theme_publication
    save_publication_figure(p_force_multi, "Figure_7F_SHAP_force_multi", 
                            width = 10, height = 6)
  }
}, error = function(e) {
  cat("  Force plot error:", e$message, "\n")
})

# --- Save model and SHAP objects ---
saveRDS(xgb_model, file.path(data_processed_dir, "xgb_model.rds"))
saveRDS(shp, file.path(data_processed_dir, "shap_object.rds"))
saveRDS(metrics, file.path(data_processed_dir, "model_metrics.rds"))

cat("\n========================================================\n")
cat("SHAP INTERPRETABILITY ANALYSIS COMPLETE\n")
cat("========================================================\n")
cat("Model AUC:", round(auc(roc_obj), 4), "\n")
cat("Model Accuracy:", round(conf_mat$overall["Accuracy"], 4), "\n")
cat("Results saved to results/figures/SHAP/ and results/tables/\n")
