setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 14_train_base_classifiers.R
# Train five tuned base classifiers using 80/20 split
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")
source("R/resampling_functions.R")
source("R/modeling_functions.R")
source("R/plotting_functions.R")
load_required_packages()

combat_expr <- readRDS(file.path(DISCOVERY_DIR, "discovery_combat.rds"))
meta <- readRDS(file.path(DISCOVERY_DIR, "discovery_cohort.rds"))$metadata
consensus <- readRDS(file.path(CONSENSUS_DIR, "consensus_genes.rds"))
final_genes <- consensus$final_consensus

cat("Training base classifiers using 80/20 split\n")
cat("Consensus genes:", length(final_genes), "\n")

# Use 80/20 split for primary analysis
set.seed(GLOBAL_SEED)
train_idx <- createDataPartition(meta$class, p = 0.8, list = TRUE)[[1]]
test_idx <- setdiff(seq_len(nrow(meta)), train_idx)

# Prepare training and hold-out data
train_expr <- t(combat_expr[final_genes, train_idx])
test_expr <- t(combat_expr[final_genes, test_idx])

train_df <- as.data.frame(train_expr)
train_df$Class <- factor(meta$class[train_idx], levels = c(NEGATIVE_CLASS, POSITIVE_CLASS))

holdout_df <- as.data.frame(test_expr)
holdout_df$Class <- factor(meta$class[test_idx], levels = c(NEGATIVE_CLASS, POSITIVE_CLASS))

cat(sprintf("Training samples: %d, Hold-out: %d\n", nrow(train_df), nrow(holdout_df)))

# Generate shared CV indices for training data
set.seed(GLOBAL_SEED)
shared_cv <- generate_shared_cv_indices(train_df$Class, CV_FOLDS, CV_REPEATS, GLOBAL_SEED)
tc <- create_train_control(shared_cv, twoClassSummary)

# Train all five models
cat("\nTraining Random Forest...\n")
rf_model <- train_rf_model(train_df, tc)
cat("  Best tune:", paste(names(rf_model$bestTune), rf_model$bestTune, sep = "=", collapse = ", "), "\n")

cat("Training Gradient Boosting...\n")
gb_model <- train_gbm_model(train_df, tc)
cat("  Best tune:", paste(names(gb_model$bestTune), gb_model$bestTune, sep = "=", collapse = ", "), "\n")

cat("Training SVM...\n")
svm_model <- train_svm_model(train_df, tc)
cat("  Best tune:", paste(names(svm_model$bestTune), svm_model$bestTune, sep = "=", collapse = ", "), "\n")

cat("Training KNN...\n")
knn_model <- train_knn_model(train_df, tc)
cat("  Best tune:", paste(names(knn_model$bestTune), knn_model$bestTune, sep = "=", collapse = ", "), "\n")

cat("Training Neural Network...\n")
nn_model <- train_nn_model(train_df, tc)
cat("  Best tune:", paste(names(nn_model$bestTune), nn_model$bestTune, sep = "=", collapse = ", "), "\n")

base_models <- list(
  rf = rf_model, gbm = gb_model, svm = svm_model,
  knn = knn_model, nnet = nn_model
)

# Save models and data
saveRDS(base_models, file.path(BASE_MODELS_DIR, "base_models.rds"))
saveRDS(train_df, file.path(BASE_MODELS_DIR, "training_data.rds"))
saveRDS(holdout_df, file.path(BASE_MODELS_DIR, "holdout_data.rds"))
saveRDS(list(train = train_idx, test = test_idx), 
        file.path(BASE_MODELS_DIR, "primary_split_indices.rds"))

cat("\nBase classifier training complete.\n")

