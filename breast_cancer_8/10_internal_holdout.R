# ============================================================
# 10 - Internal holdout evaluation
# Evaluates the locked final model on the patient-level held-out
# samples (never used for feature selection, tuning, fitting, or
# calibration).
# ============================================================
PROJECT_DIR <- "F:/independent projects/cancer/breast cancer/Breast cancer 8"
source(file.path(PROJECT_DIR, "config/config.R"))
source(file.path(PROJECT_DIR, "R/utils.R"))
suppressMessages({library(dplyr); library(ggplot2)})

logmsg <- function(...) cat(sprintf(...), "\n")

final <- readRDS(file.path(DIR_MODELS, "final_model.rds"))
split <- readRDS(file.path(DIR_MODELS, "primary_split.rds"))
disc  <- load_discovery()
nc    <- readRDS(file.path(DIR_MODELS, "nested_cv_results.rds"))

expr_ho <- disc$expr[final$genes, split$hold_samples, drop = FALSE]
meta_ho <- disc$meta[match(split$hold_samples, disc$meta$sample_id), ]
truth_ho <- meta_ho$class
logmsg("Holdout: %d samples, %d genes", ncol(expr_ho), nrow(expr_ho))
logmsg("Holdout class balance: Tumor=%d Normal=%d",
       sum(truth_ho == "Tumor"), sum(truth_ho == "Normal"))

pred <- final$predict(expr_ho)
prob <- pred$calibrated
thr  <- final$threshold

# ---------------- metrics ----------------
m <- compute_metrics(prob, truth_ho, threshold = thr)
write.csv(m, file.path(DIR_TABLES, "Table_17_Holdout_Performance.csv"), row.names = FALSE)
print(m)

out_df <- data.frame(sample_id = colnames(expr_ho),
                     patient_id = meta_ho$patient_id,
                     truth = truth_ho,
                     prob_raw = round(pred$raw, 6),
                     prob_calibrated = round(prob, 6),
                     predicted_class = ifelse(prob >= thr, "Tumor", "Normal"),
                     stringsAsFactors = FALSE)
write.csv(out_df, file.path(DIR_TABLES, "Table_17b_Holdout_Predictions.csv"), row.names = FALSE)

# OOF comparison reference
logmsg("OOF (nested CV) ROC AUC = %.4f", pROC::roc(nc$oof$truth, nc$oof$oof_prob, quiet = TRUE)$auc)

# ---------------- figures (Figure_09) ----------------
cd <- curve_data(prob, truth_ho)

# ROC
roc_df <- cd$roc
p_roc <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(colour = "#0072B2", linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  theme_bc(base_size = 9) +
  labs(x = "1 - Specificity", y = "Sensitivity",
       title = sprintf("Holdout ROC (AUC = %.3f)", cd$auc_roc))
save_tiff_600dpi(p_roc, "Figure_09_Holdout_ROC.tiff", width = 6, height = 5.5)

# PR
pr_df <- cd$pr
p_pr <- ggplot(pr_df, aes(x = recall, y = precision)) +
  geom_line(colour = "#D55E00", linewidth = 1.1) +
  geom_hline(yintercept = mean(truth_ho == "Tumor"), linetype = 2, colour = "grey50") +
  theme_bc(base_size = 9) +
  labs(x = "Recall", y = "Precision",
       title = sprintf("Holdout PR curve (AUC = %.3f)", cd$auc_pr))
save_tiff_600dpi(p_pr, "Figure_09_Holdout_PR_Curve.tiff", width = 6, height = 5.5)

# calibration curve on holdout
cal_df <- calibration_curve_df(prob, truth_ho)
p_cal <- ggplot(cal_df, aes(x = mean_pred, y = mean_obs)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey40") +
  geom_point(size = 2.5, colour = "#0072B2") + geom_line(colour = "#0072B2") +
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) +
  theme_bc(base_size = 9) +
  labs(x = "Mean predicted probability", y = "Observed proportion (Tumor)",
       title = sprintf("Holdout calibration (ECE = %.3f)",
                       expected_calibration_error(prob, truth_ho)))
save_tiff_600dpi(p_cal, "Figure_09_Holdout_Calibration.tiff", width = 6, height = 5.5)

# confusion matrix
cm <- table(Predicted = ifelse(prob >= thr, "Tumor", "Normal"), Actual = truth_ho)
cm_df <- as.data.frame(cm)
p_cm <- ggplot(cm_df, aes(x = Actual, y = Predicted, fill = Freq)) +
  geom_tile(colour = "black") +
  geom_text(aes(label = Freq), size = 7, colour = "white") +
  scale_fill_gradient(low = "#9ECAE1", high = "#08519C") +
  theme_bc(base_size = 10) +
  labs(title = sprintf("Holdout confusion matrix (threshold = %.3f)", thr),
       x = "Actual", y = "Predicted") +
  theme(legend.position = "none")
save_tiff_600dpi(p_cm, "Figure_09_Holdout_Confusion.tiff", width = 5.5, height = 5)

# decision curve
dca <- decision_curve_df(prob, truth_ho)
dca_df <- data.frame(threshold = dca$model$threshold, net_benefit = dca$model$net_benefit,
                     strategy = "Model")
dca_all <- data.frame(threshold = dca$all$threshold, net_benefit = dca$all$net_benefit,
                      strategy = "Treat all")
dca_none <- data.frame(threshold = seq(0.01, 0.99, 0.01), net_benefit = 0, strategy = "Treat none")
dca_df <- rbind(dca_df, dca_all, dca_none)
p_dca <- ggplot(dca_df, aes(x = threshold, y = net_benefit, colour = strategy)) +
  geom_line(linewidth = 1) +
  theme_bc(base_size = 9) +
  labs(x = "Threshold probability", y = "Net benefit", colour = NULL,
       title = "Decision curve analysis (holdout)") +
  theme(legend.position = "bottom")
save_tiff_600dpi(p_dca, "Figure_09_Holdout_DCA.tiff", width = 6.5, height = 5.5)

# probability distribution by class
prob_df <- data.frame(prob = prob, truth = truth_ho)
p_hist <- ggplot(prob_df, aes(x = prob, fill = truth)) +
  geom_histogram(bins = 30, colour = "black", alpha = 0.8, position = "dodge") +
  geom_vline(xintercept = thr, linetype = 2, colour = "black") +
  scale_fill_manual(values = c(Tumor = COL_TUMOR, Normal = COL_NORMAL)) +
  theme_bc(base_size = 9) +
  labs(x = "Calibrated probability", y = "Count", fill = NULL,
       title = "Holdout predicted probabilities by class") +
  theme(legend.position = "bottom")
save_tiff_600dpi(p_hist, "Figure_09_Holdout_Probability_Distribution.tiff", width = 6.5, height = 5)

# scatter raw vs calibrated
sc_df <- data.frame(raw = pred$raw, calibrated = prob, truth = truth_ho)
p_sc <- ggplot(sc_df, aes(x = raw, y = calibrated, colour = truth)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey50") +
  scale_colour_manual(values = c(Tumor = COL_TUMOR, Normal = COL_NORMAL)) +
  theme_bc(base_size = 9) +
  labs(x = "Raw XGBoost probability", y = "Platt-calibrated probability",
       colour = NULL, title = "Raw vs calibrated probabilities (holdout)") +
  theme(legend.position = "bottom")
save_tiff_600dpi(p_sc, "Figure_09_Holdout_Raw_vs_Calibrated.tiff", width = 6, height = 5.5)

logmsg("Holdout evaluation complete.")
