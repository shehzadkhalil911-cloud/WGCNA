setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 10_SVM_feature_selection.R
# SVM feature selection across allocations
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

cat("SVM Feature Selection\n")

svm_selected <- list()
svm_importance_list <- list()

for (pct in names(split_indices)) {
    cat(sprintf("\n--- Allocation: %s ---\n", pct))
    train_idx <- split_indices[[pct]]$train
    train_x <- t(combat_expr[deg_genes, train_idx])
    train_y <- factor(meta$class[train_idx], levels = c(NEGATIVE_CLASS, POSITIVE_CLASS))
    train_df <- as.data.frame(train_x)
    train_df$Class <- train_y
    
    # Pre-filter by variance if needed
    if (ncol(train_x) > MAX_FEATURE_CANDIDATES) {
        vars <- apply(train_x, 2, var)
        top_var <- names(sort(vars, decreasing = TRUE))[1:MAX_FEATURE_CANDIDATES]
        train_df <- train_df[, c(top_var, "Class")]
        train_x <- train_x[, top_var]
    }
    
    cat(sprintf("Training samples: %d, Features: %d\n", nrow(train_df), ncol(train_df) - 1))
    
    cl <- parallel::makeCluster(N_CORES)
    doParallel::registerDoParallel(cl)
    
    tc <- create_train_control(cv_indices_list[[pct]], twoClassSummary)
    
    model <- caret::train(
        Class ~ .,
        data = train_df,
        method = "svmLinear",
        trControl = tc,
        tuneGrid = SVM_TUNE_GRID,
        metric = "ROC",
        maximize = TRUE,
        verbose = FALSE
    )
    
    parallel::stopCluster(cl)
    
    # Feature importance from SVM weights
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
        selected <- colnames(train_x)  # fallback: use all
    }
    
    svm_selected[[pct]] <- selected
    svm_importance_list[[pct]] <- imp_df
    
    cat(sprintf("  Selected %d genes\n", length(selected)))
    
    write.csv(imp_df, file.path(FS_DIR, "SVM", paste0("svm_importance_", pct, ".csv")), 
              row.names = FALSE)
    write.csv(data.frame(Gene = selected), 
              file.path(FS_DIR, "SVM", paste0("svm_selected_", pct, ".csv")), 
              row.names = FALSE)
    saveRDS(model, file.path(FS_DIR, "SVM", paste0("svm_model_", pct, ".rds")))
}

svm_freq <- compute_selection_frequency(svm_selected, length(split_indices))
svm_freq <- svm_freq[order(-svm_freq$Frequency), ]
write.csv(svm_freq, file.path(FS_DIR, "SVM", "svm_selection_frequency.csv"), row.names = FALSE)

svm_stable <- identify_stable_genes(svm_freq, MIN_SELECTION_FREQUENCY)
write.csv(data.frame(Gene = svm_stable), 
          file.path(FS_DIR, "SVM", "svm_stable_genes.csv"), row.names = FALSE)

cat(sprintf("\nSVM stable genes (freq >= %.2f): %d\n", 
            MIN_SELECTION_FREQUENCY, length(svm_stable)))
cat("SVM feature selection complete.\n")

