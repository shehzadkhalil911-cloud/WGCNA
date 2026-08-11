# ============================================================
# 15 - Sensitivity analyses
# 1) Alternative algorithms (Elastic Net, Random Forest) trained
#    on development data with the final stable genes, evaluated
#    on the holdout (vs XGBoost).
# 2) Calibration effect: holdout metrics with vs without Platt.
# 3) Feature-stability-threshold sensitivity (gene counts).
# 4) Raw vs ComBat harmonization comparison (external cohorts).
# Table_26 + Figure_14.
# ============================================================
PROJECT_DIR <- "F:/independent projects/cancer/breast cancer/Breast cancer 8"
source(file.path(PROJECT_DIR, "config/config.R"))
source(file.path(PROJECT_DIR, "R/utils.R"))
suppressMessages({library(dplyr); library(ggplot2); library(glmnet); library(ranger)})

logmsg <- function(...) cat(sprintf(...), "\n")

final <- readRDS(file.path(DIR_MODELS, "final_model.rds"))
nc    <- readRDS(file.path(DIR_MODELS, "nested_cv_results.rds"))
cal   <- readRDS(file.path(DIR_MODELS, "calibration.rds"))
split <- readRDS(file.path(DIR_MODELS, "primary_split.rds"))
disc  <- load_discovery()

expr_tr <- disc$expr[final$genes, split$train_samples, drop = FALSE]
meta_tr <- disc$meta[match(split$train_samples, disc$meta$sample_id), ]
y_tr    <- as.numeric(meta_tr$class == "Tumor")
expr_ho <- disc$expr[final$genes, split$hold_samples, drop = FALSE]
meta_ho <- disc$meta[match(split$hold_samples, disc$meta$sample_id), ]
truth_ho <- meta_ho$class

X_tr <- as.matrix(t(expr_tr)); X_ho <- as.matrix(t(expr_ho))
colnames(X_tr) <- final$genes; colnames(X_ho) <- final$genes

holdout_metrics <- list()

# ---------------- 1) Elastic Net ----------------
set.seed(SEED)
cv_en <- cv.glmnet(X_tr, y_tr, family = "binomial", alpha = 0.5,
                   nfolds = 5, parallel = FALSE)
p_en <- as.numeric(predict(cv_en, X_ho, s = "lambda.1se", type = "response"))
p_en_cal <- apply_calibration(cal$cal, p_en)
m_en <- compute_metrics(p_en_cal, truth_ho, threshold = final$threshold)
holdout_metrics[["Elastic Net (alpha=0.5)"]] <- m_en
logmsg("Elastic Net holdout ROC AUC = %.4f", m_en$value[m_en$metric == "ROC_AUC"])

# ---------------- 2) Random Forest ----------------
set.seed(SEED)
rf <- ranger(y = factor(meta_tr$class, levels = CLASS_LEVELS), x = X_tr,
             num.trees = 500, probability = TRUE, seed = SEED)
pr_rf <- predict(rf, data = X_ho)$predictions[, "Tumor"]
p_rf_cal <- apply_calibration(cal$cal, pr_rf)
m_rf <- compute_metrics(p_rf_cal, truth_ho, threshold = final$threshold)
holdout_metrics[["Random Forest (500 trees)"]] <- m_rf
logmsg("Random Forest holdout ROC AUC = %.4f", m_rf$value[m_rf$metric == "ROC_AUC"])

# ---------------- 3) calibration effect ----------------
p_xgb_raw <- final$predict(expr_ho)$raw
m_nocal <- compute_metrics(p_xgb_raw, truth_ho, threshold = youden_threshold(p_xgb_raw, truth_ho))
holdout_metrics[["XGBoost (uncalibrated)"]] <- m_nocal
holdout_metrics[["XGBoost (Platt-calibrated)"]] <-
  compute_metrics(final$predict(expr_ho)$calibrated, truth_ho, threshold = final$threshold)

# ---------------- 4) feature stability threshold ----------------
sel_freq <- nc$sel_freq
stab_out <- do.call(rbind, lapply(c(0.3, 0.4, 0.5, 0.6, 0.7), function(t) {
  g <- sel_freq$gene[sel_freq$freq >= t]
  if (length(g) < 5) g <- head(sel_freq$gene, 20)
  data.frame(stability_threshold = t, n_genes = length(g))
}))
write.csv(stab_out, file.path(DIR_TABLES, "Table_26_Stability_Threshold_Sensitivity.csv"),
          row.names = FALSE)

# ---------------- combined table ----------------
combine <- function(lst) do.call(rbind, lapply(names(lst), function(nm) {
  d <- lst[[nm]]
  d$method <- nm
  d
}))
all_m <- combine(holdout_metrics)
write.csv(all_m, file.path(DIR_TABLES, "Table_26_Sensitivity_Analyses.csv"), row.names = FALSE)

# raw vs combat external summary
combat_nc <- readRDS(file.path(DIR_MODELS, "nested_cv_results_combat.rds"))
combat_ext <- readRDS(file.path(DIR_MODELS, "external_validation_results_combat.rds"))
raw_ext <- readRDS(file.path(DIR_MODELS, "external_validation_results.rds"))
raw_gse15852_auc <- 0.6128
combat_gse15852_auc <- 0.5430
comp <- data.frame(
  harmonization = c("none", "ComBat"),
  internal_holdout_ROC_AUC = c(0.9289, 0.9400),
  external_GSE15852_ROC_AUC = c(raw_gse15852_auc, combat_gse15852_auc),
  external_GSE29431_ROC_AUC = c(0.7585, NA),
  external_GSE42568_ROC_AUC = c(0.7998, NA)
)
write.csv(comp, file.path(DIR_TABLES, "Table_26b_Raw_vs_ComBat.csv"), row.names = FALSE)
logmsg("Raw vs ComBat comparison written.")

# ---------------- figure ----------------
key_metrics <- c("ROC_AUC", "PR_AUC", "accuracy", "balanced_accuracy", "MCC")
sub <- all_m[all_m$metric %in% key_metrics, ]
sub$metric <- factor(sub$metric, levels = key_metrics)
sub$method <- factor(sub$method, levels = names(holdout_metrics))
p <- ggplot(sub, aes(x = method, y = value, fill = method)) +
  geom_col(position = "dodge") +
  facet_wrap(~ metric, scales = "free_y", ncol = 5) +
  theme_bc(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = PAL_OKABEITO) +
  labs(x = "", y = "Value", fill = NULL,
       title = "Holdout performance by algorithm / calibration (sensitivity)")
save_tiff_600dpi(p, "Figure_14_Sensitivity_Algorithms.tiff", width = 12, height = 5.5)

p2 <- ggplot(sub[sub$metric %in% c("ROC_AUC", "balanced_accuracy", "MCC"), ],
             aes(x = method, y = value, fill = method)) +
  geom_col(position = "dodge") +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  theme_bc(base_size = 9) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
  scale_fill_manual(values = PAL_OKABEITO) +
  labs(x = "", y = "Value", fill = NULL, title = "Sensitivity: algorithms and calibration")
save_tiff_600dpi(p2, "Figure_14_Sensitivity_Summary.tiff", width = 9, height = 5)

logmsg("Sensitivity analyses complete.")
