setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 08_RF_feature_selection.R
# Random Forest feature selection across allocations
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")
source("R/resampling_functions.R")
source("R/modeling_functions.R")
source("R/feature_selection_functions.R")
load_required_packages()

combat_expr <- readRDS(file.path(DISCOVERY_DIR, "discovery_combat.rds"))
meta <- readRDS(file.path(DISCOVERY_DIR, "discovery_cohort.rds"))$metadata
deg_genes <- readRDS(file.path(DE_DIR, "deg_gene_list.rds"))
split_indices <- readRDS(file.path(SPLITS_DIR, "train_test_indices.rds"))
cv_indices_list <- readRDS(file.path(SPLITS_DIR, "cv_indices_list.rds"))

cat("RF Feature Selection\n")

rf_selected <- list()
rf_importance_list <- list()

for (pct in names(split_indices)) {
    cat(sprintf("\n--- Allocation: %s ---\n", pct))
    train_idx <- split_indices[[pct]]$train
    
    # Prepare training data using DEGs only
    train_expr <- t(combat_expr[deg_genes, train_idx])
    train_df <- as.data.frame(train_expr)
    train_df$Class <- factor(meta$class[train_idx], 
                             levels = c(NEGATIVE_CLASS, POSITIVE_CLASS))
    
    cat(sprintf("Training samples: %d, Features: %d\n", nrow(train_df), ncol(train_df) - 1))
    
    # Train RF with CV
    cl <- parallel::makeCluster(N_CORES)
    doParallel::registerDoParallel(cl)
    
    tc <- create_train_control(cv_indices_list[[pct]], twoClassSummary)
    
    model <- caret::train(
        Class ~ .,
        data = train_df,
        method = "ranger",
        trControl = tc,
        tuneGrid = RF_TUNE_GRID,
        metric = "ROC",
        maximize = TRUE,
        importance = "permutation",
        verbose = FALSE
    )
    
    parallel::stopCluster(cl)
    
    # Extract importance
    imp <- ranger::importance(model$finalModel)
    if (is.null(imp)) {
        imp <- varImp(model)$importance
        imp <- imp[, 1]
    }
    imp_df <- data.frame(Gene = names(imp), Importance = as.numeric(imp))
    imp_df <- imp_df[order(-imp_df$Importance), ]
    
    # Select top features (importance > 0)
    selected <- imp_df$Gene[imp_df$Importance > 0]
    if (length(selected) > TOP_RANKED_FEATURES) {
        selected <- imp_df$Gene[1:TOP_RANKED_FEATURES]
    }
    
    rf_selected[[pct]] <- selected
    rf_importance_list[[pct]] <- imp_df
    
    cat(sprintf("  Selected %d genes\n", length(selected)))
    
    # Save results
    write.csv(imp_df, file.path(FS_DIR, "RF", paste0("rf_importance_", pct, ".csv")), 
              row.names = FALSE)
    write.csv(data.frame(Gene = selected), 
              file.path(FS_DIR, "RF", paste0("rf_selected_", pct, ".csv")), 
              row.names = FALSE)
    saveRDS(model, file.path(FS_DIR, "RF", paste0("rf_model_", pct, ".rds")))
}

# Compute selection frequency
all_rf_selected <- unlist(rf_selected)
rf_freq <- compute_selection_frequency(rf_selected, length(split_indices))
rf_freq <- rf_freq[order(-rf_freq$Frequency), ]
write.csv(rf_freq, file.path(FS_DIR, "RF", "rf_selection_frequency.csv"), row.names = FALSE)

# Stable genes
rf_stable <- identify_stable_genes(rf_freq, MIN_SELECTION_FREQUENCY)
write.csv(data.frame(Gene = rf_stable), 
          file.path(FS_DIR, "RF", "rf_stable_genes.csv"), row.names = FALSE)

cat(sprintf("\nRF stable genes (freq >= %.2f): %d\n", 
            MIN_SELECTION_FREQUENCY, length(rf_stable)))
cat("RF feature selection complete.\n")

