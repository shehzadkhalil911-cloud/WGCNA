# =============================================================================
# stacking_functions.R
# Functions for out-of-fold prediction generation and stacked ensemble
# =============================================================================

#' Extract out-of-fold predictions from a caret model
extract_oof_predictions <- function(model, algorithm_name, positive_class = POSITIVE_CLASS) {
    if (is.null(model$pred)) {
        stop("Model must have savePredictions = 'final' or 'all'")
    }
    preds <- model$pred
    # Determine the column name for the positive class probability
    prob_col <- if (positive_class %in% colnames(preds)) {
        positive_class
    } else {
        grep(positive_class, colnames(preds), value = TRUE)[1]
    }
    if (is.na(prob_col)) {
        prob_col <- grep("Tumor|Cancer|Case", colnames(preds), value = TRUE)[1]
    }
    if (is.na(prob_col)) {
        prob_col <- colnames(preds)[grep("^prob", colnames(preds))[1]]
    }

    # Aggregate by row index (take mean if multiple repeats/folds)
    oof <- aggregate(preds[, prob_col],
                      by = list(rowIndex = preds$rowIndex),
                      FUN = mean, na.rm = TRUE)
    colnames(oof)[2] <- paste0(algorithm_name, "_prob")
    return(oof)
}

#' Build out-of-fold meta-feature matrix
build_oof_matrix <- function(models, algorithm_names, original_data) {
    oof_list <- list()
    for (i in seq_along(models)) {
        oof <- extract_oof_predictions(models[[i]], algorithm_names[i])
        oof_list[[i]] <- oof
    }
    # Merge all OOF predictions
    oof_matrix <- Reduce(function(x, y) merge(x, y, by = "rowIndex", all = TRUE), oof_list)
    oof_matrix$Observed <- original_data$Class[oof_matrix$rowIndex]
    oof_matrix$Sample_ID <- rownames(original_data)[oof_matrix$rowIndex]
    return(oof_matrix)
}

#' Train stacked Elastic Net meta-learner
train_stacked_elastic_net <- function(oof_matrix, alpha_seq = ELASTIC_NET_ALPHA,
                                       seed = GLOBAL_SEED) {
    set.seed(seed)
    x <- as.matrix(oof_matrix[, grep("_prob$", colnames(oof_matrix))])
    y <- oof_matrix$Observed
    y_bin <- ifelse(y == POSITIVE_CLASS, 1, 0)

    # Cross-validated glmnet
    cv_fit <- cv.glmnet(
        x = x,
        y = y_bin,
        family = "binomial",
        alpha = 0.5,  # Elastic Net (midpoint)
        type.measure = "auc",
        nfolds = 5,
        parallel = USE_PARALLEL
    )

    # Tune over alpha
    results <- list()
    for (a in alpha_seq) {
        set.seed(seed)
        cv <- cv.glmnet(
            x = x,
            y = y_bin,
            family = "binomial",
            alpha = a,
            type.measure = "auc",
            nfolds = 5,
            parallel = USE_PARALLEL
        )
        results[[as.character(a)]] <- list(
            alpha = a,
            lambda_min = cv$lambda.min,
            lambda_1se = cv$lambda.1se,
            auc_max = max(cv$cvm),
            cv_fit = cv
        )
    }

    # Select best alpha
    best_alpha_idx <- which.max(sapply(results, function(r) r$auc_max))
    best_alpha <- as.numeric(names(results)[best_alpha_idx])
    best_result <- results[[best_alpha_idx]]

    # Final model
    final_model <- glmnet(
        x = x,
        y = y_bin,
        family = "binomial",
        alpha = best_alpha,
        lambda = best_result$lambda_min
    )

    list(
        final_model = final_model,
        best_alpha = best_alpha,
        best_lambda = best_result$lambda_min,
        cv_results = results,
        cv_fit = best_result$cv_fit,
        auc = best_result$auc_max,
        features = colnames(x)
    )
}

#' Predict using stacked ensemble
predict_stacked <- function(stacked_model, meta_features) {
    x <- as.matrix(meta_features[, stacked_model$features])
    prob <- predict(stacked_model$final_model, newx = x,
                    type = "response", s = stacked_model$best_lambda)
    as.numeric(prob[, 1])
}

#' Generate hold-out meta-features
generate_holdout_meta_features <- function(models, holdout_data) {
    display_names <- c("RF", "GB", "SVM", "KNN", "NN")
    meta <- data.frame(row.names = rownames(holdout_data))
    for (i in seq_along(models)) {
        prob <- predict(models[[i]], newdata = holdout_data, type = "prob")[, POSITIVE_CLASS]
        meta[[paste0(display_names[i], "_prob")]] <- prob
    }
    return(meta)
}
