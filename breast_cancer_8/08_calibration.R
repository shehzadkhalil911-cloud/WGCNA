# ============================================================
# 08 - Probability calibration
# Fits isotonic vs Platt scaling on nested-CV OOF development
# predictions only; selects method by internal CV of the
# calibration fit; locks method + Youden threshold on OOF.
# ============================================================
PROJECT_DIR <- "F:/independent projects/cancer/breast cancer/Breast cancer 8"
source(file.path(PROJECT_DIR, "config/config.R"))
source(file.path(PROJECT_DIR, "R/utils.R"))
suppressMessages({library(dplyr); library(ggplot2)})

logmsg <- function(...) cat(sprintf(...), "\n")

nc  <- readRDS(file.path(DIR_MODELS, "nested_cv_results.rds"))
disc <- load_discovery()
split <- readRDS(file.path(DIR_MODELS, "primary_split.rds"))
meta_tr <- disc$meta[match(split$train_samples, disc$meta$sample_id), ]
pat_map <- setNames(meta_tr$patient_id, meta_tr$sample_id)

oof <- nc$oof
oof$patient_id <- pat_map[oof$sample_id]
stopifnot(all(!is.na(oof$patient_id)))

# ---------------- compare calibration methods by internal CV ----------------
set.seed(SEED)
cvv <- 5
folds <- grouped_folds(nrow(oof), oof$patient_id, cvv, seed = SEED)

methods <- c("isotonic", "platt")
cv_out <- do.call(rbind, lapply(methods, function(mth) {
  per_fold <- do.call(rbind, lapply(seq_len(cvv), function(k) {
    te <- folds[[k]]
    tr <- setdiff(seq_len(nrow(oof)), te)
    cal <- fit_calibration(oof$oof_prob[tr], oof$truth[tr], method = mth)
    p <- apply_calibration(cal, oof$oof_prob[te])
    data.frame(
      fold = k, method = mth,
      ece   = expected_calibration_error(p, oof$truth[te]),
      brier = mean((as.numeric(oof$truth[te] == "Tumor") - p)^2),
      logloss = -mean(as.numeric(oof$truth[te] == "Tumor") * log(pmin(pmax(p, 1e-6), 1 - 1e-6)) +
                      as.numeric(oof$truth[te] != "Tumor") * log(1 - pmin(pmax(p, 1e-6), 1 - 1e-6)))
    )
  }))
  data.frame(method = mth,
             ece_mean = mean(per_fold$ece), ece_sd = sd(per_fold$ece),
             brier_mean = mean(per_fold$brier), brier_sd = sd(per_fold$brier),
             logloss_mean = mean(per_fold$logloss), logloss_sd = sd(per_fold$logloss))
}))

chosen <- cv_out$method[which.min(cv_out$ece_mean)]
logmsg("Calibration method chosen by CV-ECE: %s", chosen)

# ---------------- refit on full OOF, lock threshold ----------------
cal <- fit_calibration(oof$oof_prob, oof$truth, method = chosen)
oof$cal_prob <- apply_calibration(cal, oof$oof_prob)
thr <- youden_threshold(oof$cal_prob, oof$truth)
logmsg("Youden threshold on calibrated OOF: %.4f", thr)

# metrics at locked threshold on OOF (descriptive)
m_raw <- compute_metrics(oof$oof_prob, oof$truth, threshold = 0.5)
m_cal <- compute_metrics(oof$cal_prob, oof$truth, threshold = thr)

# ---------------- calibration metrics table ----------------
full_metrics <- do.call(rbind, lapply(c("raw", "isotonic", "platt"), function(mth) {
  p <- if (mth == "raw") oof$oof_prob else apply_calibration(cal, oof$oof_prob)
  data.frame(
    method = mth,
    ECE = expected_calibration_error(p, oof$truth),
    brier = mean((as.numeric(oof$truth == "Tumor") - p)^2),
    logloss = -mean(as.numeric(oof$truth == "Tumor") * log(pmin(pmax(p, 1e-6), 1 - 1e-6)) +
                    as.numeric(oof$truth != "Tumor") * log(1 - pmin(pmax(p, 1e-6), 1 - 1e-6))),
    mcm = tryCatch({
      m <- glm(as.numeric(oof$truth == "Tumor") ~ log(pmin(pmax(p, 1e-6), 1 - 1e-6) /
                                                          (1 - pmin(pmax(p, 1e-6), 1 - 1e-6))),
               family = binomial())
      c(intercept = unname(coef(m)[1]), slope = unname(coef(m)[2]))
    }, error = function(e) c(intercept = NA, slope = NA))["slope"]
  )
}))
full_metrics$method <- factor(full_metrics$method, levels = c("raw", "isotonic", "platt"))

# ---------------- save ----------------
saveRDS(list(cal = cal, method = chosen, threshold = thr,
             cv_comparison = cv_out, full_metrics = full_metrics,
             oof_cal = oof[, c("sample_id", "truth", "oof_prob", "cal_prob")]),
        file.path(DIR_MODELS, "calibration.rds"))

write.csv(cv_out, file.path(DIR_TABLES, "Table_15_Calibration_Comparison.csv"), row.names = FALSE)

final_metrics <- data.frame(
  stage = c("pre_calibration_0.5", "post_calibration_youden"),
  threshold = c(0.5, thr),
  ROC_AUC = c(m_raw$value[m_raw$metric == "ROC_AUC"], m_cal$value[m_cal$metric == "ROC_AUC"]),
  PR_AUC = c(m_raw$value[m_raw$metric == "PR_AUC"], m_cal$value[m_cal$metric == "PR_AUC"]),
  accuracy = c(m_raw$value[m_raw$metric == "accuracy"], m_cal$value[m_cal$metric == "accuracy"]),
  balanced_accuracy = c(m_raw$value[m_raw$metric == "balanced_accuracy"], m_cal$value[m_cal$metric == "balanced_accuracy"]),
  sensitivity = c(m_raw$value[m_raw$metric == "sensitivity"], m_cal$value[m_cal$metric == "sensitivity"]),
  specificity = c(m_raw$value[m_raw$metric == "specificity"], m_cal$value[m_cal$metric == "specificity"]),
  MCC = c(m_raw$value[m_raw$metric == "MCC"], m_cal$value[m_cal$metric == "MCC"])
)
write.csv(final_metrics, file.path(DIR_TABLES, "Table_15b_Calibrated_OOF_Performance.csv"), row.names = FALSE)
write.csv(data.frame(calibration_method = chosen, threshold = thr),
          file.path(DIR_TABLES, "Table_15c_Calibration_Settings.csv"), row.names = FALSE)

# ---------------- figures ----------------
cal_df <- rbind(
  cbind(calibration_curve_df(oof$oof_prob, oof$truth), method = "raw (uncalibrated)"),
  cbind(calibration_curve_df(oof$cal_prob, oof$truth), method = sprintf("%s-calibrated", chosen))
)
p_cal <- ggplot(cal_df, aes(x = mean_pred, y = mean_obs, colour = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
  geom_point(size = 2) + geom_line() +
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
  theme_bc(base_size = 9) +
  labs(x = "Mean predicted probability", y = "Observed proportion (Tumor)",
       colour = NULL, title = "Reliability curves (OOF development)") +
  theme(legend.position = "bottom")
save_tiff_600dpi(p_cal, "Figure_11_Calibration_Reliability.tiff", width = 6.5, height = 6)

hist_df <- data.frame(
  prob = c(oof$oof_prob, oof$cal_prob),
  stage = rep(c("raw", sprintf("%s-calibrated", chosen)), each = nrow(oof))
)
p_hist <- ggplot(hist_df, aes(x = prob, fill = stage)) +
  geom_histogram(bins = 40, colour = "black", alpha = 0.7, position = "identity") +
  geom_vline(xintercept = if (chosen == "isotonic") oof$cal_prob else thr, linetype = 2) +
  theme_bc(base_size = 9) +
  labs(x = "Predicted probability", y = "Count", fill = NULL,
       title = "Distribution of predicted probabilities") +
  theme(legend.position = "bottom")
save_tiff_600dpi(p_hist, "Figure_11_Probability_Distribution.tiff", width = 6.5, height = 5)

logmsg("Calibration done: method=%s threshold=%.4f", chosen, thr)
logmsg("ECE: raw=%.4f -> %s=%.4f", full_metrics$ECE[full_metrics$method == "raw"],
       chosen, full_metrics$ECE[full_metrics$method == chosen])
