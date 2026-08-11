# =============================================================================
# feature_selection_functions.R
# Feature selection functions for each ML algorithm
# =============================================================================

#' Extract feature importance from caret models
extract_feature_importance <- function(model, algorithm) {
    tryCatch({
        if (algorithm == "rf") {
            imp <- ranger::importance(model$finalModel)
            if (is.null(imp)) {
                imp <- varImp(model)$importance
            }
            imp <- data.frame(
                Gene = rownames(imp),
                Importance = imp[, 1],
                row.names = NULL
            )
        } else if (algorithm == "gbm") {
            imp <- summary(model$finalModel, plotit = FALSE)
            imp <- data.frame(
                Gene = as.character(imp$var),
                Importance = imp$rel.inf,
                row.names = NULL
            )
        } else if (algorithm == "svmLinear") {
            imp <- varImp(model)
            if (!is.null(imp)) {
                imp <- data.frame(
                    Gene = rownames(imp$importance),
                    Importance = abs(imp$importance[, 1]),
                    row.names = NULL
                )
            } else {
                imp <- data.frame(Gene = character(0), Importance = numeric(0))
            }
        } else if (algorithm == "knn") {
            imp <- varImp(model)
            if (!is.null(imp) && nrow(imp$importance) > 0) {
                imp <- data.frame(
                    Gene = rownames(imp$importance),
                    Importance = abs(imp$importance[, 1]),
                    row.names = NULL
                )
            } else {
                imp <- data.frame(Gene = character(0), Importance = numeric(0))
            }
        } else if (algorithm == "nnet") {
            imp <- varImp(model)
            if (!is.null(imp) && nrow(imp$importance) > 0) {
                imp <- data.frame(
                    Gene = rownames(imp$importance),
                    Importance = abs(imp$importance[, 1]),
                    row.names = NULL
                )
            } else {
                imp <- data.frame(Gene = character(0), Importance = numeric(0))
            }
        } else {
            imp <- data.frame(Gene = character(0), Importance = numeric(0))
        }
        imp <- imp[order(-imp$Importance), ]
        return(imp)
    }, error = function(e) {
        return(data.frame(Gene = character(0), Importance = numeric(0)))
    })
}

#' Compute selection frequency across allocations
compute_selection_frequency <- function(selected_lists, total_splits) {
    all_genes <- unique(unlist(selected_lists))
    freq <- sapply(all_genes, function(g) {
        sum(sapply(selected_lists, function(lst) g %in% lst)) / total_splits
    })
    data.frame(Gene = all_genes, Frequency = freq, row.names = NULL)
}

#' Identify stable genes based on frequency threshold
identify_stable_genes <- function(freq_df, threshold = MIN_SELECTION_FREQUENCY) {
    freq_df[freq_df$Frequency >= threshold, "Gene", drop = TRUE]
}

#' Recursive feature elimination wrapper for SVM and KNN
rfe_feature_selection <- function(train_x, train_y, algorithm = "svmLinear",
                                   sizes = c(5, 10, 15, 20, 30, 50, 100)) {
    if (algorithm == "svmLinear") {
        model_func <- caret::svmLinear
    } else if (algorithm == "knn") {
        model_func <- caret::knn3
    } else {
        return(colnames(train_x))
    }

    control <- rfeControl(
        functions = caretFuncs,
        method = "cv",
        number = 5,
        verbose = FALSE
    )

    result <- rfe(
        x = train_x,
        y = train_y,
        sizes = sizes,
        rfeControl = control,
        method = algorithm,
        metric = "ROC",
        maximize = TRUE
    )
    return(predictors(result))
}

#' Univariate feature ranking (for KNN when RFE is too expensive)
univariate_ranking <- function(train_x, train_y, n_top = 200) {
    p_vals <- apply(train_x, 2, function(x) {
        tryCatch({
            wilcox.test(x ~ train_y)$p.value
        }, error = function(e) 1)
    })
    ranked <- names(sort(p_vals))[1:min(n_top, length(p_vals))]
    return(ranked)
}
