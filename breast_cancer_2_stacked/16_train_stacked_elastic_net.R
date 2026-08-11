setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 16_train_stacked_elastic_net.R
# Train stacked Elastic Net meta-learner
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")
source("R/stacking_functions.R")
source("R/evaluation_functions.R")
source("R/plotting_functions.R")
load_required_packages()

oof_matrix <- readRDS(file.path(STACKED_DIR, "oof_predictions.rds"))

cat("Training stacked Elastic Net meta-learner...\n")

# Train stacked Elastic Net
stacked_model <- train_stacked_elastic_net(oof_matrix, alpha_seq = ELASTIC_NET_ALPHA)

cat(sprintf("Best alpha: %.1f\n", stacked_model$best_alpha))
cat(sprintf("Best lambda: %.6f\n", stacked_model$best_lambda))
cat(sprintf("CV AUC: %.4f\n", stacked_model$auc))

# Coefficients
coefs <- as.matrix(coef(stacked_model$final_model, s = stacked_model$best_lambda))
cat("\nMeta-learner coefficients:\n")
for (i in 1:nrow(coefs)) {
    cat(sprintf("  %s: %.4f\n", rownames(coefs)[i], coefs[i, 1]))
}

# Evaluate on OOF data
prob_cols <- grep("_prob$", colnames(oof_matrix), value = TRUE)
oof_meta <- oof_matrix[, prob_cols]
stacked_probs <- predict_stacked(stacked_model, oof_meta)

# OOF evaluation
oof_metrics <- calc_metrics(stacked_probs, oof_matrix$Observed)
cat(sprintf("\nOOF Performance:\n"))
cat(sprintf("  ROC AUC: %.4f\n", oof_metrics$roc_auc))
cat(sprintf("  Accuracy: %.4f\n", oof_metrics$accuracy))
cat(sprintf("  Sensitivity: %.4f\n", oof_metrics$sensitivity))
cat(sprintf("  Specificity: %.4f\n", oof_metrics$specificity))

# Add OOF predictions to stacked model object
stacked_model$oof_probs <- stacked_probs
stacked_model$oof_metrics <- oof_metrics

# Save
saveRDS(stacked_model, file.path(STACKED_DIR, "stacked_elastic_net.rds"))
write.csv(as.data.frame(as.matrix(coefs)), 
          file.path(STACKED_DIR, "meta_learner_coefficients.csv"))

cat("\nStacked Elastic Net training complete.\n")

