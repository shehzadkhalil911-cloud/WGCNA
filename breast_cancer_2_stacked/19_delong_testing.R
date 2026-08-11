setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 19_delong_testing.R
# DeLong's test for ROC AUC comparisons
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")
source("R/evaluation_functions.R")
load_required_packages()

all_probs <- readRDS(file.path(INTVAL_DIR, "holdout_predictions.rds"))
holdout_data <- readRDS(file.path(BASE_MODELS_DIR, "holdout_data.rds"))
observed <- holdout_data$Class

cat("DeLong's ROC comparison tests\n")

# Primary comparisons: Stacked vs each base model
primary_comparisons <- list(
  c("Stacked", "RF"),
  c("Stacked", "GB"),
  c("Stacked", "SVM"),
  c("Stacked", "KNN"),
  c("Stacked", "NN")
)

primary_results <- list()
for (comp in primary_comparisons) {
    cat(sprintf("  Testing: %s vs %s\n", comp[1], comp[2]))
    test <- delong_test(all_probs[[comp[1]]], all_probs[[comp[2]]], observed)
    primary_results[[paste(comp[1], "vs", comp[2])]] <- data.frame(
        Model_A = comp[1],
        Model_B = comp[2],
        AUC_A = test$auc1,
        AUC_B = test$auc2,
        AUC_Diff = test$auc_diff,
        Statistic = test$statistic,
        P_Value = test$p_value,
        stringsAsFactors = FALSE
    )
}
primary_df <- do.call(rbind, primary_results)
primary_df$Adj_P_Value <- p.adjust(primary_df$P_Value, method = "BH")

cat("\nPrimary DeLong Results:\n")
print(primary_df, row.names = FALSE)

# All pairwise comparisons (supplementary)
prob_matrix <- as.matrix(all_probs)
all_pairwise <- run_delong_comparisons(prob_matrix, observed)
cat("\nAll Pairwise DeLong Results:\n")
print(all_pairwise, row.names = FALSE)

write.csv(primary_df, file.path(DELONG_DIR, "delong_primary_comparisons.csv"), row.names = FALSE)
write.csv(all_pairwise, file.path(DELONG_DIR, "delong_all_pairwise.csv"), row.names = FALSE)
saveRDS(list(primary = primary_df, all = all_pairwise), 
        file.path(DELONG_DIR, "delong_results.rds"))

cat("\nDeLong analysis complete.\n")

