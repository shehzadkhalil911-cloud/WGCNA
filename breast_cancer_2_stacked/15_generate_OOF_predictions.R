setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 15_generate_OOF_predictions.R
# Generate out-of-fold predictions from base classifiers
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")
source("R/stacking_functions.R")
load_required_packages()

base_models <- readRDS(file.path(BASE_MODELS_DIR, "base_models.rds"))
train_df <- readRDS(file.path(BASE_MODELS_DIR, "training_data.rds"))

cat("Extracting out-of-fold predictions...\n")

algo_names <- c("rf", "gbm", "svm", "knn", "nnet")
display_names <- c("RF", "GB", "SVM", "KNN", "NN")

# Extract OOF predictions from each model
oof_list <- list()
for (i in seq_along(base_models)) {
    cat(sprintf("  Extracting %s OOF...\n", display_names[i]))
    oof <- extract_oof_predictions(base_models[[i]], display_names[i])
    oof_list[[i]] <- oof
}

# Merge into single OOF matrix
oof_matrix <- Reduce(function(x, y) merge(x, y, by = "rowIndex", all = TRUE), oof_list)
oof_matrix$rowIndex <- NULL
rownames(oof_matrix) <- rownames(train_df)
oof_matrix$Observed <- as.character(train_df$Class)
oof_matrix$Sample_ID <- rownames(train_df)

cat("\nOOF matrix dimensions:", nrow(oof_matrix), "x", ncol(oof_matrix), "\n")
cat("OOF columns:", paste(colnames(oof_matrix), collapse = ", "), "\n")

# Basic diagnostics
prob_cols <- grep("_prob$", colnames(oof_matrix), value = TRUE)
cat("\nOOF probability ranges:\n")
for (col in prob_cols) {
    cat(sprintf("  %s: %.3f - %.3f (mean: %.3f)\n", 
                col, min(oof_matrix[[col]], na.rm = TRUE), 
                max(oof_matrix[[col]], na.rm = TRUE),
                mean(oof_matrix[[col]], na.rm = TRUE)))
}

saveRDS(oof_matrix, file.path(STACKED_DIR, "oof_predictions.rds"))
write.csv(oof_matrix, file.path(STACKED_DIR, "oof_predictions.csv"), row.names = TRUE)

cat("\nOut-of-fold predictions saved.\n")

