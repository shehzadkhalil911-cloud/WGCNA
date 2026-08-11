# =============================================================================
# resampling_functions.R
# Functions for train/test splitting and cross-validation indices
# =============================================================================

#' Generate stratified train/test indices for multiple allocations
generate_split_indices <- function(labels, allocations = ALLOCATIONS, seed = GLOBAL_SEED) {
    set.seed(seed)
    indices <- list()
    for (p in allocations) {
        train_idx <- caret::createDataPartition(labels, p = p, list = TRUE)[[1]]
        test_idx <- setdiff(seq_along(labels), train_idx)
        attr(train_idx, "allocation") <- p
        attr(test_idx, "allocation") <- p
        indices[[as.character(p)]] <- list(
            train = train_idx,
            test = test_idx
        )
    }
    return(indices)
}

#' Generate shared repeated CV fold indices
generate_shared_cv_indices <- function(train_labels, n_folds = CV_FOLDS,
                                        n_repeats = CV_REPEATS, seed = GLOBAL_SEED) {
    set.seed(seed)
    all_indices <- list()
    n_samples <- length(train_labels)
    fold_num <- 1

    for (r in 1:n_repeats) {
        set.seed(seed + r)
        folds <- caret::createFolds(train_labels, k = n_folds, list = TRUE)
        for (f in 1:n_folds) {
            # For caret: each element is the held-out samples (test indices within CV)
            all_indices[[fold_num]] <- folds[[f]]
            names(all_indices)[fold_num] <- paste0("Rep", r, ".Fold", f)
            fold_num <- fold_num + 1
        }
    }
    return(all_indices)
}

#' Create trainControl object for repeated CV
create_train_control <- function(cv_indices, summary_fn = twoClassSummary) {
    trainControl(
        method = "repeatedcv",
        number = CV_FOLDS,
        repeats = CV_REPEATS,
        classProbs = TRUE,
        summaryFunction = summary_fn,
        savePredictions = "final",
        sampling = "up",
        allowParallel = TRUE,
        index = cv_indices,
        returnResamp = "final",
        verboseIter = FALSE
    )
}

#' Stratified bootstrap resampling
stratified_bootstrap <- function(data, n_boot = N_BOOTSTRAP, seed = GLOBAL_SEED) {
    set.seed(seed)
    boot_indices <- list()
    classes <- data$Observed
    for (i in 1:n_boot) {
        idx_tumor <- which(classes == POSITIVE_CLASS)
        idx_normal <- which(classes == NEGATIVE_CLASS)
        boot_tumor <- sample(idx_tumor, size = length(idx_tumor), replace = TRUE)
        boot_normal <- sample(idx_normal, size = length(idx_normal), replace = TRUE)
        boot_indices[[i]] <- c(boot_tumor, boot_normal)
    }
    return(boot_indices)
}
