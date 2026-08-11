# =============================================================================
# evaluation_functions.R
# Functions for model evaluation, metrics, bootstrap, and DeLong tests
# =============================================================================

#' Calculate comprehensive classification metrics
calc_metrics <- function(probabilities, observed, positive_class = POSITIVE_CLASS,
                         threshold = 0.5) {
    observed_bin <- ifelse(observed == positive_class, 1, 0)
    predicted_class <- ifelse(probabilities >= threshold, positive_class,
                              setdiff(unique(observed), positive_class)[1])
    predicted_bin <- ifelse(predicted_class == positive_class, 1, 0)

    cm <- caret::confusionMatrix(
        factor(predicted_class, levels = unique(observed)),
        factor(observed, levels = unique(observed)),
        positive = positive_class
    )

    roc_obj <- tryCatch(pROC::roc(observed_bin, probabilities, quiet = TRUE),
                        error = function(e) NULL)
    pr_obj <- tryCatch(PRROC::pr.curve(
        scores.class0 = probabilities[observed_bin == 1],
        scores.class1 = probabilities[observed_bin == 0]
    ), error = function(e) NULL)

    auc_val <- if (!is.null(roc_obj)) as.numeric(roc_obj$auc) else NA
    prauc_val <- if (!is.null(pr_obj)) pr_obj$auc.integral else NA

    brier <- mean((probabilities - observed_bin)^2, na.rm = TRUE)
    logloss <- -mean(
        observed_bin * log(pmax(probabilities, 1e-15)) +
            (1 - observed_bin) * log(pmax(1 - probabilities, 1e-15)),
        na.rm = TRUE
    )

    tp <- sum(predicted_bin == 1 & observed_bin == 1)
    tn <- sum(predicted_bin == 0 & observed_bin == 0)
    fp <- sum(predicted_bin == 1 & observed_bin == 0)
    fn <- sum(predicted_bin == 0 & observed_bin == 1)

    sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA
    ppv <- if ((tp + fp) > 0) tp / (tp + fp) else NA
    npv <- if ((tn + fn) > 0) tn / (tn + fn) else NA
    precision <- ppv
    recall <- sensitivity
    f1 <- if ((precision + recall) > 0) 2 * precision * recall / (precision + recall) else NA

    n_total <- length(observed_bin)
    accuracy <- (tp + tn) / n_total
    balanced_acc <- (sensitivity + specificity) / 2

    mcc_num <- (tp * tn) - (fp * fn)
    mcc_den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc <- if (mcc_den > 0) mcc_num / mcc_den else NA

    # Cohen's kappa
    p_observed <- accuracy
    p_expected <- ((tp + fp) * (tp + fn) + (tn + fp) * (tn + fn)) / (n_total^2)
    kappa <- (p_observed - p_expected) / (1 - p_expected)

    list(
        accuracy = accuracy,
        balanced_accuracy = balanced_acc,
        sensitivity = sensitivity,
        specificity = specificity,
        ppv = ppv,
        npv = npv,
        precision = precision,
        recall = recall,
        f1 = f1,
        mcc = mcc,
        kappa = kappa,
        roc_auc = auc_val,
        pr_auc = prauc_val,
        brier_score = brier,
        log_loss = logloss,
        confusion_matrix = cm$table,
        roc_obj = roc_obj,
        pr_obj = pr_obj,
        tp = tp, tn = tn, fp = fp, fn = fn
    )
}

#' Bootstrap confidence intervals for metrics
bootstrap_metrics <- function(probabilities, observed, n_boot = N_BOOTSTRAP,
                               seed = GLOBAL_SEED, threshold = 0.5) {
    set.seed(seed)
    observed_bin <- ifelse(observed == POSITIVE_CLASS, 1, 0)
    n <- length(probabilities)
    boot_stats <- list()

    metrics_to_calc <- c("accuracy", "balanced_accuracy", "sensitivity",
                         "specificity", "ppv", "npv", "f1", "mcc",
                         "roc_auc", "pr_auc", "brier_score")

    boot_results <- matrix(NA, nrow = n_boot, ncol = length(metrics_to_calc))
    colnames(boot_results) <- metrics_to_calc
    boot_failures <- 0

    for (b in 1:n_boot) {
        idx_tumor <- which(observed_bin == 1)
        idx_normal <- which(observed_bin == 0)
        if (length(idx_tumor) > 0 && length(idx_normal) > 0) {
            boot_tumor <- sample(idx_tumor, size = length(idx_tumor), replace = TRUE)
            boot_normal <- sample(idx_normal, size = length(idx_normal), replace = TRUE)
            boot_idx <- c(boot_tumor, boot_normal)
        } else {
            boot_idx <- sample(1:n, size = n, replace = TRUE)
        }

        tryCatch({
            m <- calc_metrics(probabilities[boot_idx], observed[boot_idx],
                              threshold = threshold)
            for (j in seq_along(metrics_to_calc)) {
                boot_results[b, j] <- m[[metrics_to_calc[j]]]
            }
        }, error = function(e) {
            boot_failures <<- boot_failures + 1
        })
    }

    # Calculate CIs
    ci_list <- list()
    for (j in seq_along(metrics_to_calc)) {
        metric_name <- metrics_to_calc[j]
        vals <- boot_results[!is.na(boot_results[, j]), j]
        ci_list[[metric_name]] <- list(
            mean = mean(vals, na.rm = TRUE),
            median = median(vals, na.rm = TRUE),
            sd = sd(vals, na.rm = TRUE),
            ci_lower = quantile(vals, 0.025, na.rm = TRUE),
            ci_upper = quantile(vals, 0.975, na.rm = TRUE),
            all_values = vals
        )
    }

    list(
        ci = ci_list,
        boot_matrix = boot_results,
        n_bootstrap = n_boot,
        failures = boot_failures
    )
}

#' Paired DeLong test between two models
delong_test <- function(prob1, prob2, observed, paired = TRUE) {
    observed_bin <- ifelse(observed == POSITIVE_CLASS, 1, 0)
    roc1 <- pROC::roc(observed_bin, prob1, quiet = TRUE)
    roc2 <- pROC::roc(observed_bin, prob2, quiet = TRUE)
    test <- pROC::roc.test(roc1, roc2, method = "delong", paired = paired)
    list(
        roc1 = roc1,
        roc2 = roc2,
        auc1 = as.numeric(roc1$auc),
        auc2 = as.numeric(roc2$auc),
        auc_diff = abs(as.numeric(roc1$auc) - as.numeric(roc2$auc)),
        statistic = test$statistic,
        p_value = test$p.value
    )
}

#' Run all pairwise DeLong tests
run_delong_comparisons <- function(prob_matrix, observed) {
    model_names <- colnames(prob_matrix)
    comparisons <- list()
    k <- 1
    for (i in 1:(length(model_names) - 1)) {
        for (j in (i + 1):length(model_names)) {
            test <- delong_test(
                prob_matrix[, i],
                prob_matrix[, j],
                observed
            )
            comparisons[[k]] <- data.frame(
                Model_A = model_names[i],
                Model_B = model_names[j],
                AUC_A = test$auc1,
                AUC_B = test$auc2,
                AUC_Diff = test$auc_diff,
                Statistic = test$statistic,
                P_Value = test$p_value,
                stringsAsFactors = FALSE
            )
            k <- k + 1
        }
    }
    result <- do.call(rbind, comparisons)
    result$Adj_P_Value <- p.adjust(result$P_Value, method = "BH")
    return(result)
}

#' Determine optimal probability threshold
optimize_threshold <- function(probabilities, observed, metric = "youden") {
    observed_bin <- ifelse(observed == POSITIVE_CLASS, 1, 0)
    roc_obj <- pROC::roc(observed_bin, probabilities, quiet = TRUE)
    if (metric == "youden") {
        best <- pROC::coords(roc_obj, "best", ret = c("threshold", "sensitivity",
                                                        "specificity"))
        return(best$threshold)
    } else {
        return(0.5)
    }
}
