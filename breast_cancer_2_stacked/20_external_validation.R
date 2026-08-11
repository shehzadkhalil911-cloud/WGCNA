setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 20_external_validation.R
# Independent external validation using locked models (direct prediction)
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")
source("R/evaluation_functions.R")
load_required_packages()

cat("External Validation\n")

base_models <- readRDS(file.path(BASE_MODELS_DIR, "base_models.rds"))
stacked_model <- readRDS(file.path(STACKED_DIR, "stacked_elastic_net.rds"))
consensus <- readRDS(file.path(CONSENSUS_DIR, "consensus_genes.rds"))
final_genes <- consensus$final_consensus
train_df <- readRDS(file.path(BASE_MODELS_DIR, "training_data.rds"))

val_data <- readRDS(file.path(EXTVAL_DIR, "validation_cohort.rds"))
val_expr <- val_data$expression
val_meta <- val_data$metadata
cat("Samples:", ncol(val_expr), "\n")
cat("Available genes in validation:", nrow(val_expr), "\n")

sig_genes_available <- intersect(final_genes, rownames(val_expr))
missing_genes <- setdiff(final_genes, rownames(val_expr))
cat(sprintf("Signature genes available: %d/%d\n", length(sig_genes_available), length(final_genes)))
if (length(missing_genes) > 0) {
    cat("Missing signature genes:", paste(missing_genes, collapse = ", "), "\n")
    write.csv(data.frame(Missing_Gene = missing_genes),
              file.path(EXTVAL_RES_DIR, "missing_signature_genes.csv"), row.names = FALSE)
}

# Build validation predictor matrix
gene_means <- sapply(train_df[, final_genes, drop = FALSE], mean, na.rm = TRUE)
val_matrix <- matrix(NA_real_, nrow = ncol(val_expr), ncol = length(final_genes))
colnames(val_matrix) <- final_genes
rownames(val_matrix) <- colnames(val_expr)
for (g in final_genes) {
    if (g %in% rownames(val_expr)) {
        val_matrix[, g] <- as.numeric(val_expr[g, ])
    } else {
        val_matrix[, g] <- gene_means[g]
    }
}
val_pred_df <- as.data.frame(val_matrix)

# Impute remaining NAs with training gene means (external validation data has NAs)
na_count <- sum(is.na(val_pred_df))
if (na_count > 0) {
    cat("Imputing", na_count, "NAs with training gene means...\n")
    for (g in colnames(val_pred_df)) {
        na_idx <- is.na(val_pred_df[[g]])
        if (any(na_idx)) {
            val_pred_df[na_idx, g] <- gene_means[g]
        }
    }
}
cat("Validation predictor matrix:", nrow(val_pred_df), "x", ncol(val_pred_df), "\n")
cat("Remaining NAs:", sum(is.na(val_pred_df)), "\n")

# Predict using caret predict on all models
cat("Generating predictions...\n")
display_names <- c("RF", "GB", "SVM", "KNN", "NN")
val_probs <- data.frame(row.names = rownames(val_pred_df))

for (i in seq_along(base_models)) {
    model_name <- names(base_models)[i]
    cat(sprintf("  %s...", display_names[i]))
    pred_prob <- predict(base_models[[model_name]], newdata = val_pred_df, type = "prob")
    cat(" [returned", nrow(pred_prob), "rows]")
    val_probs[[display_names[i]]] <- pred_prob[, POSITIVE_CLASS]
    cat(" OK\n")
}

# Stacked - average of base models
val_probs$Stacked <- rowMeans(val_probs[, display_names], na.rm = TRUE)
cat("  Stacked (average) OK\n")

# Evaluate
val_metrics <- list()
for (model_name in colnames(val_probs)) {
    val_metrics[[model_name]] <- calc_metrics(val_probs[[model_name]],
        factor(val_meta$class, levels = c(NEGATIVE_CLASS, POSITIVE_CLASS)))
}

val_metrics_table <- do.call(rbind, lapply(names(val_metrics), function(m) {
    data.frame(
        Model = m,
        Accuracy = val_metrics[[m]]$accuracy,
        Balanced_Accuracy = val_metrics[[m]]$balanced_accuracy,
        Sensitivity = val_metrics[[m]]$sensitivity,
        Specificity = val_metrics[[m]]$specificity,
        PPV = val_metrics[[m]]$ppv,
        NPV = val_metrics[[m]]$npv,
        F1 = val_metrics[[m]]$f1,
        MCC = val_metrics[[m]]$mcc,
        ROC_AUC = val_metrics[[m]]$roc_auc,
        PR_AUC = val_metrics[[m]]$pr_auc,
        Brier = val_metrics[[m]]$brier_score,
        stringsAsFactors = FALSE
    )
}))

cat("\nExternal Validation Performance:\n")
print(val_metrics_table, row.names = FALSE)

write.csv(val_metrics_table, file.path(EXTVAL_RES_DIR, "external_validation_metrics.csv"), row.names = FALSE)
saveRDS(val_probs, file.path(EXTVAL_RES_DIR, "external_predictions.rds"))
saveRDS(val_metrics, file.path(EXTVAL_RES_DIR, "external_metrics.rds"))
cat("\nExternal validation complete.\n")
