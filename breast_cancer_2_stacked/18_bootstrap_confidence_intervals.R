setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 18_bootstrap_confidence_intervals.R
# Bootstrap confidence intervals for all models on hold-out set
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")
source("R/evaluation_functions.R")
load_required_packages()

all_probs <- readRDS(file.path(INTVAL_DIR, "holdout_predictions.rds"))
holdout_data <- readRDS(file.path(BASE_MODELS_DIR, "holdout_data.rds"))
observed <- holdout_data$Class

cat("Bootstrap confidence interval estimation (", N_BOOTSTRAP, " replicates)\n")

all_bootstrap <- list()

for (model_name in colnames(all_probs)) {
    cat(sprintf("  Bootstrapping %s...\n", model_name))
    boot_result <- bootstrap_metrics(all_probs[[model_name]], observed, 
                                      n_boot = N_BOOTSTRAP, seed = GLOBAL_SEED)
    all_bootstrap[[model_name]] <- boot_result
}

# Compile CI table
ci_table <- data.frame()
for (model_name in names(all_bootstrap)) {
    for (metric_name in names(all_bootstrap[[model_name]]$ci)) {
        ci_info <- all_bootstrap[[model_name]]$ci[[metric_name]]
        ci_table <- rbind(ci_table, data.frame(
            Model = model_name,
            Metric = metric_name,
            Mean = ci_info$mean,
            Median = ci_info$median,
            SD = ci_info$sd,
            CI_Lower = ci_info$ci_lower,
            CI_Upper = ci_info$ci_upper,
            stringsAsFactors = FALSE
        ))
    }
}

cat("\nBootstrap CI Summary (ROC AUC):\n")
auc_ci <- ci_table[ci_table$Metric == "roc_auc", ]
print(auc_ci, row.names = FALSE)

write.csv(ci_table, file.path(BOOT_DIR, "bootstrap_ci_table.csv"), row.names = FALSE)
saveRDS(all_bootstrap, file.path(BOOT_DIR, "bootstrap_results.rds"))

cat("\nBootstrap analysis complete.\n")

