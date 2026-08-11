setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 07_generate_data_splits.R
# Generate stratified train/test allocations and CV indices
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")
source("R/resampling_functions.R")
load_required_packages()

meta <- readRDS(file.path(DISCOVERY_DIR, "discovery_cohort.rds"))$metadata
labels <- meta$class

cat("Generating train/test splits for allocations:", 
    paste(ALLOCATIONS * 100, "%", sep = ""), "\n")

# Generate split indices
split_indices <- generate_split_indices(labels, allocations = ALLOCATIONS, seed = GLOBAL_SEED)

# Summary table
split_summary <- data.frame()
for (pct in names(split_indices)) {
    train_idx <- split_indices[[pct]]$train
    test_idx <- split_indices[[pct]]$test
    train_tab <- table(labels[train_idx])
    test_tab <- table(labels[test_idx])
    split_summary <- rbind(split_summary, data.frame(
        Allocation = paste0(as.numeric(pct) * 100, "%"),
        Train_Total = length(train_idx),
        Train_Tumor = train_tab["Tumor"],
        Train_Normal = train_tab["Normal"],
        Test_Total = length(test_idx),
        Test_Tumor = test_tab["Tumor"],
        Test_Normal = test_tab["Normal"],
        stringsAsFactors = FALSE
    ))
}

cat("\nTrain/Test Split Summary:\n")
print(split_summary)
write.csv(split_summary, file.path(SPLITS_DIR, "split_summary.csv"), row.names = FALSE)

# Generate shared CV indices for EACH allocation
cv_indices_list <- list()
for (pct in names(split_indices)) {
    train_idx <- split_indices[[pct]]$train
    train_labels <- labels[train_idx]
    cv_indices <- generate_shared_cv_indices(train_labels, 
                                              n_folds = CV_FOLDS, 
                                              n_repeats = CV_REPEATS,
                                              seed = GLOBAL_SEED)
    cv_indices_list[[pct]] <- cv_indices
    cat(sprintf("\nAllocation %s: %d CV resamples generated\n", 
                pct, length(cv_indices)))
}

# Save all
saveRDS(split_indices, file.path(SPLITS_DIR, "train_test_indices.rds"))
saveRDS(cv_indices_list, file.path(SPLITS_DIR, "cv_indices_list.rds"))
saveRDS(split_summary, file.path(SPLITS_DIR, "split_summary.rds"))

cat("\nData splits and CV indices saved.\n")

