# ============================================================
# 11 - Bootstrap confidence intervals
# Patient-level stratified bootstrap (2000 replicates) of the
# holdout sample to quantify uncertainty around performance
# metrics. Table_18 + bootstrap distribution figures.
# ============================================================
PROJECT_DIR <- "F:/independent projects/cancer/breast cancer/Breast cancer 8"
source(file.path(PROJECT_DIR, "config/config.R"))
source(file.path(PROJECT_DIR, "R/utils.R"))
suppressMessages({library(dplyr); library(ggplot2)})

logmsg <- function(...) cat(sprintf(...), "\n")

final <- readRDS(file.path(DIR_MODELS, "final_model.rds"))
split <- readRDS(file.path(DIR_MODELS, "primary_split.rds"))
disc  <- load_discovery()

expr_ho <- disc$expr[final$genes, split$hold_samples, drop = FALSE]
meta_ho <- disc$meta[match(split$hold_samples, disc$meta$sample_id), ]
truth_ho <- meta_ho$class
pred <- final$predict(expr_ho)
prob <- pred$calibrated
thr  <- final$threshold

patients <- unique(meta_ho$patient_id)
pat_class <- sapply(patients, function(p) unique(meta_ho$class[meta_ho$patient_id == p]))
n_pat <- length(patients)
logmsg("Bootstrap: %d patients, %d replicates", n_pat, N_BOOT)

metrics_of <- function(pr, tr, t) {
  m <- compute_metrics(pr, tr, threshold = t)
  setNames(m$value, m$metric)
}

set.seed(SEED)
obs <- metrics_of(prob, truth_ho, thr)

boot_list <- lapply(seq_len(N_BOOT), function(b) {
  set.seed(SEED + b)
  pats <- sample(patients, n_pat, replace = TRUE)
  s <- unlist(lapply(pats, function(p) which(meta_ho$patient_id == p)))
  mm <- metrics_of(prob[s], truth_ho[s], thr)
  mm[c("ROC_AUC", "PR_AUC", "accuracy", "balanced_accuracy", "sensitivity",
       "specificity", "PPV", "NPV", "F1", "MCC", "cohens_kappa", "brier_score", "log_loss")]
})
boot <- as.data.frame(do.call(rbind, boot_list))

ci_tab <- do.call(rbind, lapply(names(obs), function(metric) {
  v <- obs[metric]
  if (!(metric %in% colnames(boot))) return(NULL)
  bv <- boot[[metric]]
  data.frame(metric = metric, estimate = unname(v),
             ci_lower = unname(quantile(bv, 0.025, na.rm = TRUE)),
             ci_upper = unname(quantile(bv, 0.975, na.rm = TRUE)),
             bootstrap_sd = sd(bv, na.rm = TRUE),
             n_boot = N_BOOT)
}))
write.csv(ci_tab, file.path(DIR_TABLES, "Table_18_Bootstrap_Confidence_Intervals.csv"),
          row.names = FALSE)
print(ci_tab)

# ---------------- figures ----------------
boot_long <- tidyr::pivot_longer(boot, everything(), names_to = "metric", values_to = "value")
boot_long$metric <- factor(boot_long$metric,
                           levels = c("ROC_AUC", "PR_AUC", "accuracy", "balanced_accuracy",
                                      "sensitivity", "specificity", "PPV", "NPV", "F1",
                                      "MCC", "cohens_kappa", "brier_score", "log_loss"))
obs_df <- data.frame(metric = names(obs), value = unname(obs))
obs_df <- obs_df[obs_df$metric %in% levels(boot_long$metric), ]

p_boot <- ggplot(boot_long, aes(x = value)) +
  geom_histogram(bins = 40, fill = "#0072B2", colour = "black", alpha = 0.85) +
  geom_vline(data = obs_df, aes(xintercept = value), colour = "#D55E00", linewidth = 0.8, linetype = 2) +
  facet_wrap(~ metric, scales = "free", ncol = 4) +
  theme_bc(base_size = 8) +
  labs(x = "Bootstrap value", y = "Count",
       title = sprintf("Bootstrap distributions (n = %d, patient-level)", N_BOOT))
save_tiff_600dpi(p_boot, "Figure_10_Bootstrap_Distributions.tiff", width = 10, height = 7)

# AUC bootstrap focus
auc_boot <- data.frame(roc_auc = boot$ROC_AUC, pr_auc = boot$PR_AUC)
p_auc <- ggplot(auc_boot, aes(x = roc_auc)) +
  geom_histogram(bins = 40, fill = "#0072B2", colour = "black", alpha = 0.85) +
  geom_vline(xintercept = obs[["ROC_AUC"]], colour = "#D55E00", linewidth = 0.8, linetype = 2) +
  theme_bc(base_size = 9) +
  labs(x = "Bootstrap ROC AUC", y = "Count",
       title = sprintf("ROC AUC bootstrap (95%% CI: %.3f - %.3f)",
                       ci_tab$ci_lower[ci_tab$metric == "ROC_AUC"],
                       ci_tab$ci_upper[ci_tab$metric == "ROC_AUC"]))
save_tiff_600dpi(p_auc, "Figure_10_Bootstrap_AUC.tiff", width = 6.5, height = 5)

logmsg("Bootstrap CI complete.")
