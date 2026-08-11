setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 09_GB_feature_selection.R
# Gradient Boosting feature selection across allocations
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

cat("GB Feature Selection\n")

gb_selected <- list()
gb_importance_list <- list()

for (pct in names(split_indices)) {
    cat(sprintf("\n--- Allocation: %s ---\n", pct))
    train_idx <- split_indices[[pct]]$train
    train_expr <- t(combat_expr[deg_genes, train_idx])
    train_df <- as.data.frame(train_expr)
    train_df$Class <- factor(meta$class[train_idx], 
                             levels = c(NEGATIVE_CLASS, POSITIVE_CLASS))
    
    # Limit features for GB to avoid excessive runtime
    if (ncol(train_df) - 1 > MAX_FEATURE_CANDIDATES) {
        # Pre-filter by variance
        vars <- apply(train_expr, 2, var)
        top_var <- names(sort(vars, decreasing = TRUE))[1:MAX_FEATURE_CANDIDATES]
        train_df <- train_df[, c(top_var, "Class")]
    }
    
    cat(sprintf("Training samples: %d, Features: %d\n", nrow(train_df), ncol(train_df) - 1))
    
    cl <- parallel::makeCluster(N_CORES)
    doParallel::registerDoParallel(cl)
    
    tc <- create_train_control(cv_indices_list[[pct]], twoClassSummary)
    
    model <- caret::train(
        Class ~ .,
        data = train_df,
        method = "gbm",
        trControl = tc,
        tuneGrid = GB_TUNE_GRID,
        metric = "ROC",
        maximize = TRUE,
        verbose = FALSE
    )
    
    parallel::stopCluster(cl)
    
    # Extract importance
    imp <- summary(model$finalModel, plotit = FALSE)
    imp_df <- data.frame(Gene = as.character(imp$var), 
                         Importance = imp$rel.inf,
                         stringsAsFactors = FALSE)
    imp_df <- imp_df[order(-imp_df$Importance), ]
    
    selected <- imp_df$Gene[1:min(TOP_RANKED_FEATURES, nrow(imp_df))]
    gb_selected[[pct]] <- selected
    gb_importance_list[[pct]] <- imp_df
    
    cat(sprintf("  Selected %d genes\n", length(selected)))
    
    write.csv(imp_df, file.path(FS_DIR, "GB", paste0("gb_importance_", pct, ".csv")), 
              row.names = FALSE)
    write.csv(data.frame(Gene = selected), 
              file.path(FS_DIR, "GB", paste0("gb_selected_", pct, ".csv")), 
              row.names = FALSE)
    saveRDS(model, file.path(FS_DIR, "GB", paste0("gb_model_", pct, ".rds")))
}

gb_freq <- compute_selection_frequency(gb_selected, length(split_indices))
gb_freq <- gb_freq[order(-gb_freq$Frequency), ]
write.csv(gb_freq, file.path(FS_DIR, "GB", "gb_selection_frequency.csv"), row.names = FALSE)

gb_stable <- identify_stable_genes(gb_freq, MIN_SELECTION_FREQUENCY)
write.csv(data.frame(Gene = gb_stable), 
          file.path(FS_DIR, "GB", "gb_stable_genes.csv"), row.names = FALSE)

cat(sprintf("\nGB stable genes (freq >= %.2f): %d\n", 
            MIN_SELECTION_FREQUENCY, length(gb_stable)))
cat("GB feature selection complete.\n")

