# ============================================================
# 14 - External validation (multi-cohort)
# Independent cohorts: GSE15852 (GPL96), GSE29431 (GPL570),
# GSE42568 (GPL570). Sample-level and patient-level performance
# of the locked final model with bootstrap CIs (N_BOOT_EXT).
# ============================================================
PROJECT_DIR <- "F:/independent projects/cancer/breast cancer/Breast cancer 8"
source(file.path(PROJECT_DIR, "config/config.R"))
source(file.path(PROJECT_DIR, "R/utils.R"))
suppressMessages({library(dplyr); library(ggplot2)})

logmsg <- function(...) cat(sprintf(...), "\n")

final <- readRDS(file.path(DIR_MODELS, "final_model.rds"))

cohort_paths <- list(
  GSE15852 = file.path(DIR_VAL_PROC, "validation_processed.rds"),
  GSE29431 = file.path(DIR_EXT_PROC, "GSE29431_processed.rds"),
  GSE42568 = file.path(DIR_EXT_PROC, "GSE42568_processed.rds")
)

metrics_of <- function(pr, tr, t) {
  mm <- compute_metrics(pr, tr, threshold = t)
  setNames(mm$value, mm$metric)
}

all_pred <- list()
sum_tab  <- list()

for (cohort in names(cohort_paths)) {
  logmsg("===== %s =====", cohort)
  obj  <- readRDS(cohort_paths[[cohort]])
  expr_v <- obj$expr[final$genes, , drop = FALSE]
  meta_v <- obj$meta
  truth_v <- meta_v$class
  logmsg("  %d samples (%d normal / %d tumor), %d genes",
         ncol(expr_v), sum(truth_v == "Normal"), sum(truth_v == "Tumor"), nrow(expr_v))

  pred <- final$predict(expr_v)
  prob <- pred$calibrated
  thr  <- final$threshold

  # patient-level: paired cohorts -> within-patient consistency;
  # unpaired cohorts -> sample-level is patient-level
  pt <- meta_v$patient_id
  paired <- any(table(pt) > 1)
  if (paired) {
    wid <- as.data.frame(cbind(patient_id = pt, truth = truth_v, prob = prob))
    pw <- tidyr::pivot_wider(wid, id_cols = patient_id, names_from = truth,
                             values_from = prob)
    pw <- pw[stats::complete.cases(pw), ]
    paired_consistency <- mean(pw$Tumor > pw$Normal)
    pt_cons <- sum(pw$Tumor > pw$Normal) / nrow(pw)
    m_p <- list(paired_consistency = pt_cons, n_patients = nrow(pw))
    logmsg("  Paired consistency (tumor prob > normal prob): %.3f (%d/%d patients)",
           pt_cons, sum(pw$Tumor > pw$Normal), nrow(pw))
  } else {
    pt_cons <- NA
    m_p <- NULL
  }

  m_s <- metrics_of(prob, truth_v, thr)

  # per-cohort Youden-optimal operating point (threshold-dependent view)
  thr_y <- youden_threshold(prob, truth_v)
  m_y <- metrics_of(prob, truth_v, thr_y)

  # bootstrap CIs (patient-level resampling)
  patients <- unique(meta_v$patient_id)
  n_pat <- length(patients)
  boot <- do.call(rbind, lapply(seq_len(N_BOOT_EXT), function(b) {
    set.seed(SEED + b)
    pats <- sample(patients, n_pat, replace = TRUE)
    s <- unlist(lapply(pats, function(p) which(meta_v$patient_id == p)))
    mm <- metrics_of(prob[s], truth_v[s], thr)
    mm[c("ROC_AUC", "PR_AUC", "accuracy", "balanced_accuracy", "sensitivity",
         "specificity", "PPV", "NPV", "F1", "MCC", "brier_score")]
  }))
  boot <- as.data.frame(boot)

  ci_tab <- do.call(rbind, lapply(names(m_s), function(mname) {
    if (!(mname %in% colnames(boot))) return(NULL)
    data.frame(metric = mname, estimate = unname(m_s[mname]),
               ci_lower = unname(quantile(boot[[mname]], 0.025, na.rm = TRUE)),
               ci_upper = unname(quantile(boot[[mname]], 0.975, na.rm = TRUE)),
               bootstrap_sd = sd(boot[[mname]], na.rm = TRUE))
  }))
  ci_tab$cohort <- cohort
  sum_tab[[cohort]] <- ci_tab

  out_s <- data.frame(cohort = cohort, sample_id = colnames(expr_v),
                      patient_id = meta_v$patient_id, truth = truth_v,
                      prob_raw = round(pred$raw, 6), prob_calibrated = round(prob, 6),
                      predicted_class = ifelse(prob >= thr, "Tumor", "Normal"),
                      stringsAsFactors = FALSE)
  all_pred[[cohort]] <- list(prob = prob, truth = truth_v, sample = out_s,
                             paired = paired, paired_consistency = pt_cons,
                             youden_metrics = m_y, youden_threshold = thr_y,
                             locked_metrics = m_s)

  if (paired) {
    logmsg("  Sample ROC AUC = %.4f | accuracy = %.3f | MCC = %.3f | paired consistency = %.3f",
           m_s["ROC_AUC"], m_s["accuracy"], m_s["MCC"], pt_cons)
  } else {
    logmsg("  Sample ROC AUC = %.4f | accuracy = %.3f | MCC = %.3f (locked threshold)",
           m_s["ROC_AUC"], m_s["accuracy"], m_s["MCC"])
  }
  write.csv(out_s, file.path(DIR_TABLES, sprintf("Table_23_%s_Sample_Predictions.csv", cohort)),
            row.names = FALSE)
}

# combined summary table
comb <- do.call(rbind, sum_tab)
comb <- comb[, c("cohort", "metric", "estimate", "ci_lower", "ci_upper", "bootstrap_sd")]
comb_wide <- tidyr::pivot_wider(comb, id_cols = metric, names_from = cohort,
                                values_from = c(estimate, ci_lower, ci_upper))
write.csv(comb, file.path(DIR_TABLES, "Table_22_External_Validation_Performance.csv"),
          row.names = FALSE)
write.csv(comb_wide, file.path(DIR_TABLES, "Table_22b_External_Validation_Summary_Wide.csv"),
          row.names = FALSE)
print(comb[comb$metric %in% c("ROC_AUC", "PR_AUC", "accuracy", "balanced_accuracy",
                              "sensitivity", "specificity", "MCC"), ])

# cohort-specific Youden-optimal operating points
youden_tab <- do.call(rbind, lapply(names(all_pred), function(cn) {
  y <- all_pred[[cn]]$youden_metrics
  data.frame(cohort = cn, youden_threshold = all_pred[[cn]]$youden_threshold,
             metric = names(y), value = unname(y))
}))
write.csv(youden_tab,
          file.path(DIR_TABLES, "Table_22c_External_Youden_Optimal_Points.csv"), row.names = FALSE)

# paired consistency summary
pair_tab <- do.call(rbind, lapply(names(all_pred), function(cn) {
  a <- all_pred[[cn]]
  data.frame(cohort = cn, paired = a$paired,
             paired_consistency = a$paired_consistency)
}))
write.csv(pair_tab, file.path(DIR_TABLES, "Table_24_External_Paired_Consistency.csv"),
          row.names = FALSE)

# ---------------- figures (Figure_13) ----------------
cd_roc <- lapply(names(all_pred), function(cn) {
  cd <- curve_data(all_pred[[cn]]$prob, all_pred[[cn]]$truth)
  data.frame(cohort = cn, specificity = cd$roc$specificity,
             sensitivity = cd$roc$sensitivity, auc = cd$auc_roc)
})
roc_all <- do.call(rbind, cd_roc)
auc_lab <- unique(roc_all[, c("cohort", "auc")])
p_roc <- ggplot(roc_all, aes(x = 1 - specificity, y = sensitivity, colour = cohort)) +
  geom_line(linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  geom_label(data = auc_lab, aes(x = 0.72, y = rev(seq(0.30, 0.55, length.out = nrow(auc_lab))),
                                 label = sprintf("%s (AUC=%.3f)", cohort, auc), colour = cohort),
             size = 3, show.legend = FALSE) +
  theme_bc(base_size = 9) +
  labs(x = "1 - Specificity", y = "Sensitivity", colour = NULL,
       title = "External validation ROC (sample-level)") +
  theme(legend.position = "bottom")
save_tiff_600dpi(p_roc, "Figure_13_External_ROC_All.tiff", width = 7, height = 6)

for (cn in names(all_pred)) {
  d <- all_pred[[cn]]
  cd <- curve_data(d$prob, d$truth)
  p1 <- ggplot(cd$roc, aes(x = 1 - specificity, y = sensitivity)) +
    geom_line(colour = "#0072B2", linewidth = 1.1) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
    theme_bc(base_size = 9) +
    labs(x = "1 - Specificity", y = "Sensitivity",
         title = sprintf("%s ROC (AUC = %.3f)", cn, cd$auc_roc))
  save_tiff_600dpi(p1, sprintf("Figure_13_%s_ROC.tiff", cn), width = 6, height = 5.5)

  p2 <- ggplot(cd$pr, aes(x = recall, y = precision)) +
    geom_line(colour = "#D55E00", linewidth = 1.1) +
    geom_hline(yintercept = mean(d$truth == "Tumor"), linetype = 2, colour = "grey50") +
    theme_bc(base_size = 9) +
    labs(x = "Recall", y = "Precision", title = sprintf("%s PR curve (AUC = %.3f)", cn, cd$auc_pr))
  save_tiff_600dpi(p2, sprintf("Figure_13_%s_PR_Curve.tiff", cn), width = 6, height = 5.5)

  cal_df <- calibration_curve_df(d$prob, d$truth)
  p3 <- ggplot(cal_df, aes(x = mean_pred, y = mean_obs)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
    geom_point(size = 2.5, colour = "#0072B2") + geom_line(colour = "#0072B2") +
    scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
    theme_bc(base_size = 9) +
    labs(x = "Mean predicted probability", y = "Observed proportion (Tumor)",
         title = sprintf("%s calibration (ECE = %.3f)", cn,
                         expected_calibration_error(d$prob, d$truth)))
  save_tiff_600dpi(p3, sprintf("Figure_13_%s_Calibration.tiff", cn), width = 6, height = 5.5)

  prob_df <- data.frame(prob = d$prob, truth = d$truth)
  p4 <- ggplot(prob_df, aes(x = prob, fill = truth)) +
    geom_histogram(bins = 30, colour = "black", alpha = 0.8, position = "dodge") +
    geom_vline(xintercept = final$threshold, linetype = 2, colour = "black") +
    scale_fill_manual(values = c(Tumor = COL_TUMOR, Normal = COL_NORMAL)) +
    theme_bc(base_size = 9) +
    labs(x = "Calibrated probability", y = "Count", fill = NULL,
         title = sprintf("%s predicted probabilities", cn)) +
    theme(legend.position = "bottom")
  save_tiff_600dpi(p4, sprintf("Figure_13_%s_Probability_Distribution.tiff", cn),
                   width = 6.5, height = 5)
}

saveRDS(list(predictions = lapply(all_pred, function(x) x[c("prob", "truth")]),
             sample_tables = lapply(all_pred, `[[`, "sample"),
             patient_tables = lapply(all_pred, `[[`, "patient"),
             summary = comb),
        file.path(DIR_MODELS, "external_validation_results.rds"))
logmsg("External validation complete (all cohorts).")
