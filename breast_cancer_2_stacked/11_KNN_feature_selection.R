setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 11_KNN_feature_selection.R
# KNN feature selection across allocations
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

cat("KNN Feature Selection\n")

knn_selected <- list()
knn_importance_list <- list()

for (pct in names(split_indices)) {
    cat(sprintf("\n--- Allocation: %s ---\n", pct))
    train_idx <- split_indices[[pct]]$train
    train_x <- t(combat_expr[deg_genes, train_idx])
    train_y <- factor(meta$class[train_idx], levels = c(NEGATIVE_CLASS, POSITIVE_CLASS))
    
    # KNN feature selection using univariate ranking (training data only)
    cat("  Performing univariate feature ranking...\n")
    n_rank <- min(TOP_RANKED_FEATURES, ncol(train_x))
    selected <- univariate_ranking(train_x, train_y, n_top = n_rank)
    
    cat(sprintf("  Selected %d genes\n", length(selected)))
    
    # Create importance-like ranking
    imp_df <- data.frame(
        Gene = selected,
        Importance = rev(seq_len(length(selected))),
        stringsAsFactors = FALSE
    )
    
    knn_selected[[pct]] <- selected
    knn_importance_list[[pct]] <- imp_df
    
    write.csv(imp_df, file.path(FS_DIR, "KNN", paste0("knn_ranking_", pct, ".csv")), 
              row.names = FALSE)
    write.csv(data.frame(Gene = selected), 
              file.path(FS_DIR, "KNN", paste0("knn_selected_", pct, ".csv")), 
              row.names = FALSE)
}

knn_freq <- compute_selection_frequency(knn_selected, length(split_indices))
knn_freq <- knn_freq[order(-knn_freq$Frequency), ]
write.csv(knn_freq, file.path(FS_DIR, "KNN", "knn_selection_frequency.csv"), row.names = FALSE)

knn_stable <- identify_stable_genes(knn_freq, MIN_SELECTION_FREQUENCY)
write.csv(data.frame(Gene = knn_stable), 
          file.path(FS_DIR, "KNN", "knn_stable_genes.csv"), row.names = FALSE)

cat(sprintf("\nKNN stable genes (freq >= %.2f): %d\n", 
            MIN_SELECTION_FREQUENCY, length(knn_stable)))
cat("KNN feature selection complete.\n")

