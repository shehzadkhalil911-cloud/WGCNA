setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 17_internal_holdout_evaluation.R
# Evaluate all models on untouched internal hold-out set
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")
source("R/stacking_functions.R")
source("R/evaluation_functions.R")
source("R/plotting_functions.R")
load_required_packages()

base_models <- readRDS(file.path(BASE_MODELS_DIR, "base_models.rds"))
holdout_df <- readRDS(file.path(BASE_MODELS_DIR, "holdout_data.rds"))
stacked_model <- readRDS(file.path(STACKED_DIR, "stacked_elastic_net.rds"))

cat("Internal hold-out evaluation\n")
cat("Hold-out samples:", nrow(holdout_df), "\n")
print(table(holdout_df$Class))

algo_names <- c("rf", "gbm", "svm", "knn", "nnet")
display_names <- c("RF", "GB", "SVM", "KNN", "NN")

# Predict on hold-out
all_probs <- data.frame(row.names = rownames(holdout_df))
all_metrics <- list()

cat("\nIndividual model predictions:\n")
for (i in seq_along(base_models)) {
    prob <- predict(base_models[[i]], newdata = holdout_df, type = "prob")[, POSITIVE_CLASS]
    all_probs[[display_names[i]]] <- prob
    metrics <- calc_metrics(prob, holdout_df$Class)
    all_metrics[[display_names[i]]] <- metrics
    cat(sprintf("  %s AUC: %.4f, Acc: %.4f\n", 
                display_names[i], metrics$roc_auc, metrics$accuracy))
}

# Stacked model predictions
meta_features <- generate_holdout_meta_features(base_models, holdout_df)
stacked_probs <- predict_stacked(stacked_model, meta_features)
all_probs$Stacked <- stacked_probs
all_metrics[["Stacked"]] <- calc_metrics(stacked_probs, holdout_df$Class)
cat(sprintf("  Stacked AUC: %.4f, Acc: %.4f\n", 
            all_metrics[["Stacked"]]$roc_auc, all_metrics[["Stacked"]]$accuracy))

# Combined metrics table
metrics_table <- do.call(rbind, lapply(names(all_metrics), function(m) {
    data.frame(
        Model = m,
        Accuracy = all_metrics[[m]]$accuracy,
        Balanced_Accuracy = all_metrics[[m]]$balanced_accuracy,
        Sensitivity = all_metrics[[m]]$sensitivity,
        Specificity = all_metrics[[m]]$specificity,
        PPV = all_metrics[[m]]$ppv,
        NPV = all_metrics[[m]]$npv,
        F1 = all_metrics[[m]]$f1,
        MCC = all_metrics[[m]]$mcc,
        Kappa = all_metrics[[m]]$kappa,
        ROC_AUC = all_metrics[[m]]$roc_auc,
        PR_AUC = all_metrics[[m]]$pr_auc,
        Brier = all_metrics[[m]]$brier_score,
        LogLoss = all_metrics[[m]]$log_loss,
        stringsAsFactors = FALSE
    )
}))
cat("\nHold-out Performance Summary:\n")
print(metrics_table, row.names = FALSE)

write.csv(metrics_table, file.path(INTVAL_DIR, "holdout_metrics.csv"), row.names = FALSE)
saveRDS(all_probs, file.path(INTVAL_DIR, "holdout_predictions.rds"))
saveRDS(all_metrics, file.path(INTVAL_DIR, "holdout_metrics.rds"))
saveRDS(meta_features, file.path(INTVAL_DIR, "holdout_meta_features.rds"))

cat("\nInternal hold-out evaluation complete.\n")

