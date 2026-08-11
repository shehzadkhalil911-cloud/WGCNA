setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 12_NN_feature_selection.R
# Neural Network feature selection across allocations
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")
source("R/resampling_functions.R")
source("R/feature_selection_functions.R")
load_required_packages()

combat_expr <- readRDS(file.path(DISCOVERY_DIR, "discovery_combat.rds"))
meta <- readRDS(file.path(DISCOVERY_DIR, "discovery_cohort.rds"))$metadata
deg_genes <- readRDS(file.path(DE_DIR, "deg_gene_list.rds"))
split_indices <- readRDS(file.path(SPLITS_DIR, "train_test_indices.rds"))
cv_indices_list <- readRDS(file.path(SPLITS_DIR, "cv_indices_list.rds"))

cat("NN Feature Selection\n")

nn_selected <- list()
nn_importance_list <- list()

for (pct in names(split_indices)) {
    cat(sprintf("\n--- Allocation: %s ---\n", pct))
    train_idx <- split_indices[[pct]]$train
    train_x <- t(combat_expr[deg_genes, train_idx])
    train_y <- factor(meta$class[train_idx], levels = c(NEGATIVE_CLASS, POSITIVE_CLASS))
    train_df <- as.data.frame(train_x)
    train_df$Class <- train_y
    
    # Pre-filter to top features by variance for NN efficiency
    max_nn_features <- 200
    if (ncol(train_x) > max_nn_features) {
        vars <- apply(train_x, 2, var, na.rm = TRUE)
        top_var <- names(sort(vars, decreasing = TRUE))[1:max_nn_features]
        train_df <- train_df[, c(top_var, "Class")]
    }
    
    cat(sprintf("Training samples: %d, Features: %d\n", nrow(train_df), ncol(train_df) - 1))
    
    cl <- parallel::makeCluster(N_CORES)
    doParallel::registerDoParallel(cl)
    
    tc <- create_train_control(cv_indices_list[[pct]], twoClassSummary)
    
    model <- caret::train(
        Class ~ .,
        data = train_df,
        method = "nnet",
        trControl = tc,
        tuneGrid = NN_TUNE_GRID,
        metric = "ROC",
        maximize = TRUE,
        preProcess = c("center", "scale"),
        trace = FALSE,
        verbose = FALSE,
        MaxNWts = 10000
    )
    
    parallel::stopCluster(cl)
    
    # Feature importance
    imp <- varImp(model)
    if (!is.null(imp) && nrow(imp$importance) > 0) {
        imp_df <- data.frame(
            Gene = rownames(imp$importance),
            Importance = abs(imp$importance[, 1]),
            stringsAsFactors = FALSE
        )
        imp_df <- imp_df[order(-imp_df$Importance), ]
        selected <- imp_df$Gene[1:min(TOP_RANKED_FEATURES, nrow(imp_df))]
    } else {
        imp_df <- data.frame(Gene = character(0), Importance = numeric(0))
        selected <- colnames(train_x)
    }
    
    nn_selected[[pct]] <- selected
    nn_importance_list[[pct]] <- imp_df
    
    cat(sprintf("  Selected %d genes\n", length(selected)))
    
    write.csv(imp_df, file.path(FS_DIR, "NN", paste0("nn_importance_", pct, ".csv")), 
              row.names = FALSE)
    write.csv(data.frame(Gene = selected), 
              file.path(FS_DIR, "NN", paste0("nn_selected_", pct, ".csv")), 
              row.names = FALSE)
    saveRDS(model, file.path(FS_DIR, "NN", paste0("nn_model_", pct, ".rds")))
}

nn_freq <- compute_selection_frequency(nn_selected, length(split_indices))
nn_freq <- nn_freq[order(-nn_freq$Frequency), ]
write.csv(nn_freq, file.path(FS_DIR, "NN", "nn_selection_frequency.csv"), row.names = FALSE)

nn_stable <- identify_stable_genes(nn_freq, MIN_SELECTION_FREQUENCY)
write.csv(data.frame(Gene = nn_stable), 
          file.path(FS_DIR, "NN", "nn_stable_genes.csv"), row.names = FALSE)

cat(sprintf("\nNN stable genes (freq >= %.2f): %d\n", 
            MIN_SELECTION_FREQUENCY, length(nn_stable)))
cat("NN feature selection complete.\n")

