# ============================================================
# Shared utility functions for the Breast cancer 8 study
# ============================================================

#' Save a ggplot as a 600 dpi LZW-compressed TIFF in the figures folder
save_tiff_600dpi <- function(
    plot,
    filename,
    width = 7,
    height = 6
) {
  ggplot2::ggsave(
    filename = file.path(PROJECT_DIR, "figures", filename),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = FIG_DPI,
    compression = "lzw",
    bg = "white",
    device = "tiff"
  )
}

#' Open a base-graphics TIFF device in the figures folder
open_tiff_base <- function(filename, width = 7, height = 6) {
  grDevices::tiff(
    filename = file.path(PROJECT_DIR, "figures", filename),
    width = width,
    height = height,
    units = "in",
    res = FIG_DPI,
    compression = "lzw",
    bg = "white"
  )
}

#' Consistent publication theme
theme_bc <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.border     = ggplot2::element_rect(colour = "grey40", fill = NA),
      legend.position  = "bottom",
      legend.title     = ggplot2::element_text(size = ggplot2::rel(0.9)),
      legend.text      = ggplot2::element_text(size = ggplot2::rel(0.85)),
      axis.text        = ggplot2::element_text(colour = "black"),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA)
    )
}

#' Compute full classification performance metrics from a probability vector
compute_metrics <- function(prob, truth, threshold = 0.5, positive = "Tumor") {
  truth_f <- factor(truth, levels = CLASS_LEVELS)
  truth01 <- as.numeric(truth_f == positive)
  pred01  <- as.numeric(prob >= threshold)
  pred_f  <- factor(ifelse(pred01 == 1, positive, setdiff(CLASS_LEVELS, positive)),
                    levels = CLASS_LEVELS)

  cm <- table(Predicted = pred_f, Actual = truth_f)
  TP <- sum(pred01 == 1 & truth01 == 1)
  FP <- sum(pred01 == 1 & truth01 == 0)
  FN <- sum(pred01 == 0 & truth01 == 1)
  TN <- sum(pred01 == 0 & truth01 == 0)

  eps <- 1e-9
  acc      <- (TP + TN) / max(TP + TN + FP + FN, eps)
  bal_acc  <- 0.5 * (TP / max(TP + FN, eps) + TN / max(TN + FP, eps))
  sens     <- TP / max(TP + FN, eps)
  spec     <- TN / max(TN + FP, eps)
  ppv      <- TP / max(TP + FP, eps)
  npv      <- TN / max(TN + FN, eps)
  prec     <- ppv
  rec      <- sens
  f1       <- 2 * prec * rec / max(prec + rec, eps)
  mcc_num  <- (TP * TN - FP * FN)
  mcc_den  <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc      <- ifelse(mcc_den == 0, 0, mcc_num / mcc_den)

  # Cohen's kappa
  po <- acc
  pe <- ((TP + FN) * (TP + FP) + (TN + FP) * (TN + FN)) / max((TP + TN + FP + FN)^2, eps)
  kappa <- (po - pe) / max(1 - pe, eps)

  # ROC / PR AUC
  if (length(unique(truth01)) > 1 && all(is.finite(prob)) && stats::var(prob) > 0) {
    roc_obj <- pROC::roc(truth01, prob, quiet = TRUE, direction = "<")
    roc_auc <- as.numeric(pROC::auc(roc_obj))
    pr_obj  <- PRROC::pr.curve(scores.class0 = prob[truth01 == 1],
                               scores.class1 = prob[truth01 == 0], curve = FALSE)
    pr_auc  <- pr_obj$auc.integral
    roc_ci  <- tryCatch(as.numeric(pROC::ci.auc(roc_obj, method = "delong")), error = function(e) c(NA, NA, NA))
  } else {
    roc_auc <- NA_real_; pr_auc <- NA_real_; roc_ci <- c(NA, NA, NA)
  }

  # log loss and Brier
  prob_c <- pmin(pmax(prob, eps), 1 - eps)
  logloss <- -mean(truth01 * log(prob_c) + (1 - truth01) * log(1 - prob_c))
  brier   <- mean((truth01 - prob)^2)

  # calibration intercept/slope (logistic recalibration on logit(p))
  cal <- tryCatch({
    m <- stats::glm(truth01 ~ log(prob_c / (1 - prob_c)), family = stats::binomial())
    c(intercept = unname(stats::coef(m)[1]),
      slope     = unname(stats::coef(m)[2]))
  }, error = function(e) c(intercept = NA_real_, slope = NA_real_))

  data.frame(
    metric = c("ROC_AUC", "PR_AUC", "accuracy", "balanced_accuracy", "sensitivity",
               "specificity", "PPV", "NPV", "precision", "recall", "F1", "MCC",
               "cohens_kappa", "brier_score", "log_loss",
               "calibration_intercept", "calibration_slope"),
    value = c(roc_auc, pr_auc, acc, bal_acc, sens, spec, ppv, npv, prec, rec, f1, mcc,
              kappa, brier, logloss, cal["intercept"], cal["slope"]),
    stringsAsFactors = FALSE
  )
}

#' Build ROC + PR curve data frames for plotting
curve_data <- function(prob, truth, positive = "Tumor") {
  truth01 <- as.numeric(factor(truth, levels = CLASS_LEVELS) == positive)
  roc_obj <- pROC::roc(truth01, prob, quiet = TRUE, direction = "<")
  roc_df  <- data.frame(specificity = roc_obj$specificities,
                        sensitivity = roc_obj$sensitivities)
  roc_df <- roc_df[!duplicated(roc_df), ]
  pr <- PRROC::pr.curve(scores.class0 = prob[truth01 == 1],
                        scores.class1 = prob[truth01 == 0], curve = TRUE)
  pr_df <- data.frame(recall = pr$curve[, 1], precision = pr$curve[, 2])
  list(roc = roc_df, pr = pr_df,
       auc_roc = as.numeric(pROC::auc(roc_obj)), auc_pr = pr$auc.integral)
}

#' Youden optimal threshold
youden_threshold <- function(prob, truth, positive = "Tumor") {
  truth01 <- as.numeric(factor(truth, levels = CLASS_LEVELS) == positive)
  roc_obj <- pROC::roc(truth01, prob, quiet = TRUE, direction = "<")
  th  <- pROC::coords(roc_obj, x = "best", best.method = "youden",
                      ret = c("threshold", "specificity", "sensitivity"))
  thr <- unname(th[["threshold"]])
  if (length(thr) > 1) thr <- thr[1]
  thr
}

#' Probe-to-gene collapse using the probe with the highest IQR
collapse_probes_to_genes <- function(expr_mat, probe2gene, method = "maxIQR") {
  stopifnot(nrow(expr_mat) == length(probe2gene))
  sym  <- probe2gene
  emat <- expr_mat
  uniq <- unique(sym)
  iqr  <- apply(emat, 1, function(x) unname(stats::IQR(x, na.rm = TRUE)))
  ord  <- order(match(sym, uniq), -iqr)
  emat <- emat[ord, , drop = FALSE]
  dup  <- duplicated(sym[ord])
  emat <- emat[!dup, , drop = FALSE]
  rownames(emat) <- sym[ord][!dup]
  emat
}

#' Near-zero variance filtering (caret C implementation)
near_zero_variance <- function(mat, freq_cut = 95/5, unique_cut = 10) {
  nz <- caret::nearZeroVar(as.matrix(t(mat)), freqCut = freq_cut,
                           uniqueCut = unique_cut, names = TRUE, allowParallel = FALSE)
  alln <- rownames(mat)
  keep <- rep(TRUE, length(alln))
  names(keep) <- alln
  if (length(nz)) keep[nz] <- FALSE
  keep
}

#' Correlation filter: greedily remove highly correlated rows
correlation_filter <- function(mat, cutoff = 0.90) {
  cm <- stats::cor(t(mat), method = "pearson")
  cm[!is.finite(cm)] <- 0
  keep <- rep(TRUE, nrow(cm))
  for (i in seq_len(nrow(cm))) {
    if (!keep[i]) next
    hc <- which(keep & (abs(cm[i, ]) > cutoff))
    hc <- hc[hc != i]
    if (length(hc)) keep[hc] <- FALSE
  }
  rownames(cm)[keep]
}

#' Training-fold-only limma differential ranking
limma_rank_genes <- function(expr_mat, group, patient_id = NULL, top = NULL) {
  if (is.numeric(group)) group <- ifelse(group == 1, "Tumor", "Normal")
  group <- factor(group, levels = CLASS_LEVELS)
  if (length(unique(group)) < 2) stop("Need >= 2 groups for limma ranking")
  if (!is.null(patient_id) && length(unique(patient_id)) > 1) {
    block <- factor(patient_id)
    design <- stats::model.matrix(~ 0 + group)
    colnames(design) <- gsub("group", "", colnames(design))
    dupcor <- limma::duplicateCorrelation(expr_mat, design = design, block = block)
    fit <- limma::lmFit(expr_mat, design = design, block = block,
                        correlation = dupcor$consensus)
    cont <- makeContrasts("Tumor-Normal", levels = design)
    fit2 <- limma::contrasts.fit(fit, cont)
  } else {
    design <- stats::model.matrix(~ 0 + group)
    colnames(design) <- gsub("group", "", colnames(design))
    fit2 <- limma::lmFit(expr_mat, design = design)
  }
  fit2 <- limma::eBayes(fit2)
  tt <- limma::topTable(fit2, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt <- tt[order(tt$adj.P.Val, -abs(tt$logFC)), ]
  if (!is.null(top)) tt <- head(tt, top)
  tt
}

#' Staged leakage-safe feature selection on a training matrix
select_stable_features <- function(expr_tr, group_tr, patient_tr,
                                   top_deg = FOLD_TOP_DEG,
                                   corr_cutoff = CORR_CUTOFF,
                                   nzv = TRUE) {
  if (nzv) {
    nzv_keep <- near_zero_variance(expr_tr)
    expr_tr <- expr_tr[nzv_keep, , drop = FALSE]
  }
  tt <- limma_rank_genes(expr_tr, group_tr, patient_tr, top = top_deg)
  feats <- tt$gene
  if (length(feats) > 2) {
    sub <- expr_tr[feats, , drop = FALSE]
    sub <- sub[stats::complete.cases(sub), , drop = FALSE]
    if (nrow(sub) > 2) {
      keep <- correlation_filter(sub, cutoff = corr_cutoff)
      feats <- intersect(feats, keep)
    }
  }
  feats
}

#' Fit an xgboost binary classifier returning the model and predictions
fit_xgboost <- function(expr_tr, y_tr, expr_te = NULL, params = NULL,
                        nrounds = NULL, early_stopping = TUNE_EARLY_STOP,
                        seed = SEED) {
  if (is.null(params)) {
    params <- list(objective = "binary:logistic", eval_metric = "auc",
                   eta = 0.05, max_depth = 4, min_child_weight = 2,
                   subsample = 0.8, colsample_bytree = 0.7, gamma = 0,
                   alpha = 0, lambda = 1, nthread = 1)
  }
  if (is.null(nrounds)) nrounds <- TUNE_NROUNDS_MAX
  set.seed(seed)
  dtrain <- xgboost::xgb.DMatrix(as.matrix(t(expr_tr)), label = y_tr)
  model <- xgboost::xgb.train(
    params = params, data = dtrain,
    nrounds = nrounds, verbose = 0,
    early_stopping_rounds = early_stopping
  )
  if (!is.null(expr_te)) {
    dtest <- xgboost::xgb.DMatrix(as.matrix(t(expr_te)))
    pred  <- stats::predict(model, dtest)
  } else pred <- NULL
  list(model = model, pred = pred)
}

#' Predict from an xgboost model using a gene-by-sample matrix
predict_xgb <- function(model, expr_mat) {
  dm <- xgboost::xgb.DMatrix(as.matrix(t(expr_mat)))
  stats::predict(model, dm)
}

#' Fit probability calibration
fit_calibration <- function(prob, truth, method = "isotonic") {
  truth01 <- as.numeric(factor(truth, levels = CLASS_LEVELS) == "Tumor")
  if (method == "isotonic") {
    cal <- stats::isoreg(prob, truth01)
  } else if (method == "platt") {
    prob_c <- pmin(pmax(prob, 1e-6), 1 - 1e-6)
    df <- data.frame(y = truth01, z = log(prob_c / (1 - prob_c)))
    cal <- stats::glm(y ~ z, data = df, family = stats::binomial())
  } else stop("Unknown calibration method")
  list(model = cal, method = method)
}

#' Apply fitted calibration to new probabilities
apply_calibration <- function(cal, prob) {
  if (cal$method == "isotonic") {
    p <- stats::approx(cal$model$x, cal$model$yf, xout = prob, rule = 2)$y
    pmin(pmax(p, 0), 1)
  } else if (cal$method == "platt") {
    prob_c <- pmin(pmax(prob, 1e-6), 1 - 1e-6)
    z <- log(prob_c / (1 - prob_c))
    eta <- stats::predict(cal$model, newdata = data.frame(z = z), type = "response")
    as.numeric(eta)
  }
}

#' Calibration curve data (binned reliability)
calibration_curve_df <- function(prob, truth, n_bins = 10) {
  truth01 <- as.numeric(factor(truth, levels = CLASS_LEVELS) == "Tumor")
  bins <- cut(prob, breaks = seq(0, 1, length.out = n_bins + 1), include.lowest = TRUE)
  df <- data.frame(prob = prob, truth = truth01, bin = bins)
  agg <- stats::aggregate(truth ~ bin, df, mean)
  pred <- stats::aggregate(prob ~ bin, df, mean)
  out <- data.frame(bin = agg$bin, mean_pred = pred$prob, mean_obs = agg$truth)
  out[stats::complete.cases(out), ]
}

#' ECE (expected calibration error)
expected_calibration_error <- function(prob, truth, n_bins = 10) {
  cd <- calibration_curve_df(prob, truth, n_bins)
  truth01 <- as.numeric(factor(truth, levels = CLASS_LEVELS) == "Tumor")
  bins <- cut(prob, breaks = seq(0, 1, length.out = n_bins + 1), include.lowest = TRUE)
  tab <- as.numeric(table(bins)) / length(truth01)
  w <- tab[match(as.character(cd$bin), levels(bins))]
  sum(w * abs(cd$mean_pred - cd$mean_obs), na.rm = TRUE)
}

#' Decision curve analysis
decision_curve_df <- function(prob, truth) {
  truth01 <- as.numeric(factor(truth, levels = CLASS_LEVELS) == "Tumor")
  n <- length(truth01)
  prev <- mean(truth01)
  th <- seq(0.01, 0.99, 0.01)
  res <- do.call(rbind, lapply(th, function(t) {
    pred <- as.numeric(prob >= t)
    tp <- sum(pred == 1 & truth01 == 1)
    fp <- sum(pred == 1 & truth01 == 0)
    nb  <- tp / n - fp / n * (t / (1 - t))
    data.frame(threshold = t, net_benefit = nb)
  }))
  res <- rbind(data.frame(threshold = 0.01, net_benefit = 0), res)
  list(model = res, all = data.frame(threshold = th, net_benefit = prev), none = 0)
}

#' XGBoost CV tuning for a list of parameter configs using custom folds
xgb_cv_configs <- function(expr_tr, y_tr, configs, folds, nrounds_max = TUNE_NROUNDS_MAX,
                           early_stop = TUNE_EARLY_STOP, seed = SEED) {
  requireNamespace("xgboost", quietly = TRUE)
  dtrain <- xgboost::xgb.DMatrix(as.matrix(t(expr_tr)), label = y_tr)
  res <- lapply(seq_along(configs), function(i) {
    prm <- configs[[i]]
    prm$objective <- "binary:logistic"
    prm$eval_metric <- "auc"
    prm$nthread <- 1
    set.seed(seed + i)
    cv <- tryCatch(
      xgboost::xgb.cv(
        params = prm, data = dtrain, nfold = length(folds),
        folds = folds, nrounds = nrounds_max, verbose = 0,
        early_stopping_rounds = early_stop, prediction = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(cv)) return(NULL)
    tr <- cv$evaluation_log
    if (nrow(tr) == 0) return(NULL)
    best_iter <- tr$iter[which.max(tr$test_auc_mean)]
    best_auc  <- tr$test_auc_mean[which.max(tr$test_auc_mean)]
    data.frame(
      config_id = i,
      best_iter = best_iter,
      best_auc  = best_auc,
      eta = prm$eta, max_depth = prm$max_depth,
      min_child_weight = prm$min_child_weight, subsample = prm$subsample,
      colsample_bytree = prm$colsample_bytree, gamma = prm$gamma,
      alpha = prm$alpha, lambda = prm$lambda
    )
  })
  do.call(rbind, res)
}

#' Generate grouped CV folds (patients within a fold together)
grouped_folds <- function(n, patient_id, k, seed = SEED) {
  set.seed(seed)
  patients <- unique(patient_id)
  pidx <- sample(seq_along(patients))
  grp <- rep(seq_len(k), length.out = length(patients))[rank(pidx, ties.method = "first")]
  grp <- grp[order(pidx)]
  folds <- lapply(seq_len(k), function(i) {
    pats <- patients[grp == i]
    which(patient_id %in% pats)
  })
  folds
}
