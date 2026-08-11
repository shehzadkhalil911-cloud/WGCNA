# =============================================================================
# preprocessing_functions.R
# Functions for data preprocessing, normalization, and quality control
# =============================================================================

#' Assess whether expression data needs log2 transformation
assess_log2_transform <- function(expr_matrix) {
    vals <- as.numeric(as.matrix(expr_matrix))
    vals <- vals[!is.na(vals) & is.finite(vals)]
    summary_stats <- list(
        min = min(vals),
        max = max(vals),
        median = median(vals),
        mean = mean(vals),
        q25 = quantile(vals, 0.25),
        q75 = quantile(vals, 0.75),
        has_negative = any(vals < 0),
        range = diff(range(vals))
    )
    # Heuristic: if max > 100 and min >= 0, likely raw intensities -> log2
    # If range < 30 and values between 0-20, likely already log2
    needs_log2 <- FALSE
    if (summary_stats$max > 100 && summary_stats$min >= 0) {
        needs_log2 <- TRUE
    }
    if (summary_stats$range < 30 && summary_stats$max < 25) {
        needs_log2 <- FALSE
    }
    summary_stats$needs_log2 <- needs_log2
    return(summary_stats)
}

#' Apply log2 transformation safely
safe_log2 <- function(x) {
    x <- as.matrix(x)
    x[x < 0] <- NA
    log2(x + 1)
}

#' Quantile normalize expression matrix
quantile_normalize <- function(expr_matrix) {
    expr_matrix <- as.matrix(expr_matrix)
    normalized <- normalize.quantiles(expr_matrix, copy = TRUE)
    rownames(normalized) <- rownames(expr_matrix)
    colnames(normalized) <- colnames(expr_matrix)
    return(normalized)
}

#' Remove low-expression probes
filter_low_expression <- function(expr_matrix, min_threshold = NULL) {
    if (is.null(min_threshold)) {
        min_threshold <- quantile(as.numeric(expr_matrix), 0.05, na.rm = TRUE)
    }
    keep <- apply(expr_matrix, 1, function(x) {
        sum(x > min_threshold, na.rm = TRUE) > (0.5 * ncol(expr_matrix))
    })
    return(expr_matrix[keep, , drop = FALSE])
}

#' Detect and optionally remove sample outliers
detect_outliers <- function(expr_matrix, method = "correlation", threshold = 0.8) {
    cor_mat <- cor(expr_matrix, method = "spearman", use = "pairwise.complete.obs")
    mean_cor <- colMeans(cor_mat, na.rm = TRUE)
    outliers <- names(mean_cor[mean_cor < threshold])
    return(list(
        outlier_samples = outliers,
        mean_correlations = mean_cor,
        cor_matrix = cor_mat
    ))
}

#' Calculate RLE (Relative Log Expression)
calc_rle <- function(expr_matrix) {
    expr_matrix <- as.matrix(expr_matrix)
    pseudo_ref <- apply(expr_matrix, 1, median, na.rm = TRUE)
    rle_matrix <- sweep(expr_matrix, 1, pseudo_ref, "-")
    return(rle_matrix)
}

#' Center and scale using training parameters
center_scale_from_train <- function(x, center_vec = NULL, scale_vec = NULL) {
    if (!is.null(center_vec) && !is.null(scale_vec)) {
        x <- sweep(x, 1, center_vec, "-")
        x <- sweep(x, 1, scale_vec, "/")
    }
    return(x)
}

#' Get centering/scaling params from training data
get_center_scale_params <- function(train_matrix) {
    list(
        center = rowMeans(train_matrix, na.rm = TRUE),
        scale = apply(train_matrix, 1, sd, na.rm = TRUE)
    )
}
