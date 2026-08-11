###############################################################################
# Script 07: Machine Learning Feature Selection
# Methods: LASSO, SVM-RFE, Random Forest
###############################################################################

set.seed(123)

library(glmnet)
library(e1071)
library(randomForest)
library(caret)
library(ggplot2)
library(tidyverse)
library(viridis)
library(pROC)

base_dir <- normalizePath("breast_cancer_integrated_analysis")
data_processed_dir <- file.path(base_dir, "data_processed")
fig_dir <- file.path(base_dir, "results", "figures", "Machine_Learning")
tables_dir <- file.path(base_dir, "results", "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helper: Save 600 dpi ---
save_publication_figure <- function(plot, filename, width = 8, height = 6) {
  ggsave(file.path(fig_dir, paste0(filename, ".tiff")), plot = plot,
         width = width, height = height, dpi = 600, compression = "lzw")
  ggsave(file.path(fig_dir, paste0(filename, ".pdf")), plot = plot,
         width = width, height = height, dpi = 600)
}

theme_publication <- theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  )

# --- Load data ---
cat("Loading expression and DEG data...\n")
exprs_after <- readRDS(file.path(data_processed_dir, "merged_exprs_after_batch.rds"))
pheno <- readRDS(file.path(data_processed_dir, "merged_pheno_with_batch.rds"))
deg_genes <- readRDS(file.path(data_processed_dir, "deg_genes.rds"))

cat("DEG genes:", length(deg_genes), "\n")

if (length(deg_genes) == 0) {
  # If no DEGs found at strict threshold, use relaxed
  cat("No DEGs at strict threshold. Checking for relaxed DEGs...\n")
  deg_table <- readRDS(file.path(data_processed_dir, "limma_deg_table.rds"))
  deg_genes <- rownames(deg_table[deg_table$adj.P.Val < 0.05 & abs(deg_table$logFC) > 0.5, ])
  cat("Relaxed DEG genes:", length(deg_genes), "\n")
  
  if (length(deg_genes) == 0) {
    deg_genes <- rownames(head(deg_table[deg_table$adj.P.Val < 0.05, ], 100))
    cat("Using top DEGs by p-value:", length(deg_genes), "\n")
  }
  
  if (length(deg_genes) == 0) {
    stop("No DEGs found at any threshold. Cannot proceed with ML.")
  }
}

# Use DEG expression data
common_genes <- intersect(deg_genes, rownames(exprs_after))
deg_exprs_raw <- exprs_after[common_genes, ]

# Remove genes with missing values (NAs from ComBat)
na_count_per_gene <- apply(deg_exprs_raw, 1, function(x) sum(is.na(x)))
cat("Genes with NAs:", sum(na_count_per_gene > 0), "/", length(common_genes), "\n")
if (sum(na_count_per_gene > 0) > 0) {
  cat("  Max NAs per gene:", max(na_count_per_gene), "\n")
}
keep_genes <- na_count_per_gene == 0
deg_exprs <- t(deg_exprs_raw[keep_genes, ])

missing_genes <- setdiff(deg_genes, rownames(exprs_after))
if (length(missing_genes) > 0) {
  cat("Missing genes in expression matrix:", length(missing_genes), "\n")
}
cat("Using", ncol(deg_exprs), "DEGs for ML (after removing NAs)\n")

# Create response variable
response <- factor(pheno$group, levels = c("Normal", "Tumor"))
response_binary <- ifelse(response == "Tumor", 1, 0)

cat("Class distribution:\n")
print(table(response))

# --- Stratified train-test split ---
set.seed(123)
train_idx <- createDataPartition(response, p = 0.7, list = FALSE)
train_data <- deg_exprs[train_idx, ]
test_data <- deg_exprs[-train_idx, ]
train_response <- response[train_idx]
test_response <- response[-train_idx]
train_binary <- response_binary[train_idx]
test_binary <- response_binary[-train_idx]

cat("\nTraining samples:", nrow(train_data), "\n")
cat("Testing samples:", nrow(test_data), "\n")

# ======================================
# 1. LASSO Logistic Regression
# ======================================
cat("\n========================================================\n")
cat("LASSO Logistic Regression\n")
cat("========================================================\n")

# Prepare matrices for glmnet
x_train <- as.matrix(train_data)
y_train <- train_binary

# Cross-validation for lambda
set.seed(123)
cv_lasso <- cv.glmnet(x_train, y_train, family = "binomial", 
                       alpha = 1, nfolds = 5, type.measure = "class")

# LASSO coefficient profile
lasso_fit <- glmnet(x_train, y_train, family = "binomial", alpha = 1)

# Plot: LASSO coefficient profile
cat("Generating LASSO coefficient profile...\n")
pdf(file.path(fig_dir, "Figure_4A_LASSO_coefficient_profile.pdf"), 
    width = 10, height = 6)
plot(lasso_fit, xvar = "lambda", label = TRUE)
title("LASSO Coefficient Profile", line = 2.5)
dev.off()

tiff(file.path(fig_dir, "Figure_4A_LASSO_coefficient_profile.tiff"), 
     width = 10, height = 6, units = "in", res = 600, compression = "lzw")
plot(lasso_fit, xvar = "lambda", label = TRUE)
title("LASSO Coefficient Profile", line = 2.5)
dev.off()

# Plot: CV plot
cat("Generating LASSO cross-validation plot...\n")
pdf(file.path(fig_dir, "Figure_4B_LASSO_cross_validation.pdf"), 
    width = 8, height = 6)
plot(cv_lasso)
title("LASSO Cross-Validation (5-fold)", line = 2.5)
dev.off()

tiff(file.path(fig_dir, "Figure_4B_LASSO_cross_validation.tiff"), 
     width = 8, height = 6, units = "in", res = 600, compression = "lzw")
plot(cv_lasso)
title("LASSO Cross-Validation (5-fold)", line = 2.5)
dev.off()

# Extract selected features (lambda.min or lambda.1se)
lasso_coef <- coef(cv_lasso, s = "lambda.min")
lasso_genes <- rownames(lasso_coef)[which(lasso_coef[, 1] != 0)]
lasso_genes <- setdiff(lasso_genes, "(Intercept)")
cat("LASSO selected genes (lambda.min):", length(lasso_genes), "\n")
cat("LASSO genes:", lasso_genes, "\n")

# Also try lambda.1se
lasso_coef_1se <- coef(cv_lasso, s = "lambda.1se")
lasso_genes_1se <- rownames(lasso_coef_1se)[which(lasso_coef_1se[, 1] != 0)]
lasso_genes_1se <- setdiff(lasso_genes_1se, "(Intercept)")

# Save selection results
lasso_results <- list(
  genes = lasso_genes,
  genes_1se = lasso_genes_1se,
  cv_fit = cv_lasso,
  fit = lasso_fit
)
saveRDS(lasso_results, file.path(data_processed_dir, "lasso_results.rds"))
writeLines(lasso_genes, file.path(tables_dir, "LASSO_selected_genes.txt"))

svm_genes <- character(0)

# ======================================
# 2. SVM Feature Selection
# ======================================
cat("\n========================================================\n")
cat("SVM Feature Selection\n")
cat("========================================================\n")

# SVM with recursive feature elimination would be complex; use SVM weights instead
# Train SVM with linear kernel
svm_fit <- tryCatch({
  svm(x = train_data, y = train_response, kernel = "linear", 
      scale = TRUE, cost = 1, probability = TRUE, cross = 5)
}, error = function(e) {
  cat("SVM training error:", e$message, "\n")
  # Try with smaller subset
  cat("Retrying with top 50 features...\n")
  svm(x = train_data[, 1:min(50, ncol(train_data))], 
      y = train_response, kernel = "linear", 
      scale = TRUE, cost = 1, probability = TRUE, cross = 5)
})

# Extract feature weights
if (!is.null(svm_fit)) {
  svm_weights <- t(svm_fit$coefs) %*% svm_fit$SV
  svm_importance <- abs(as.numeric(svm_weights))
  names(svm_importance) <- colnames(svm_fit$SV)
  svm_importance <- sort(svm_importance, decreasing = TRUE)
  
  # Select top features (heuristic: weight > median or top N)
  svm_genes <- names(svm_importance[svm_importance > median(svm_importance)])
  cat("SVM selected genes:", length(svm_genes), "\n")
  
  # Save
  writeLines(svm_genes, file.path(tables_dir, "SVM_selected_genes.txt"))
  
  # CV error and accuracy from the SVM (cross = 5)
  cat("SVM CV accuracy:", svm_fit$accuracies, "\n")
  cat("SVM CV mean accuracy:", mean(svm_fit$accuracies), "\n")
  cat("SVM total accuracy:", svm_fit$tot.accuracy, "\n")
  
  # Create CV results data frame
  cv_results <- data.frame(
    Fold = 1:length(svm_fit$accuracies),
    Accuracy = svm_fit$accuracies,
    Error = 100 - svm_fit$accuracies
  )
  
  # SVM CV error plot
  p_svm_cv <- ggplot(cv_results, aes(x = factor(Fold), y = Error)) +
    geom_bar(stat = "identity", fill = "#e41a1c", alpha = 0.8) +
    geom_hline(yintercept = mean(cv_results$Error), linetype = "dashed", 
               color = "darkred", linewidth = 1) +
    labs(title = "SVM 5-Fold Cross-Validation Error",
         x = "Fold", y = "Error Rate (%)") +
    theme_publication +
    ylim(0, max(cv_results$Error) * 1.2)
  save_publication_figure(p_svm_cv, "Figure_5A_SVM_CV_error")
  
  # SVM CV accuracy plot
  p_svm_acc <- ggplot(cv_results, aes(x = factor(Fold), y = Accuracy)) +
    geom_bar(stat = "identity", fill = "#377eb8", alpha = 0.8) +
    geom_hline(yintercept = mean(cv_results$Accuracy), linetype = "dashed", 
               color = "darkblue", linewidth = 1) +
    labs(title = "SVM 5-Fold Cross-Validation Accuracy",
         x = "Fold", y = "Accuracy (%)") +
    theme_publication +
    ylim(0, 100)
  save_publication_figure(p_svm_acc, "Figure_5B_SVM_CV_accuracy")
}

# ======================================
# 3. Random Forest Feature Selection
# ======================================
cat("\n========================================================\n")
cat("Random Forest Feature Selection\n")
cat("========================================================\n")

# Train Random Forest
rf_fit <- randomForest(
  x = train_data, y = train_response,
  ntree = 500, mtry = max(1, floor(sqrt(ncol(train_data)))),
  importance = TRUE, nodesize = 5
)

cat("RF OOB error rate:", rf_fit$err.rate[nrow(rf_fit$err.rate), 1], "\n")

# Feature importance
rf_importance <- importance(rf_fit)
rf_importance_df <- data.frame(
  Gene = rownames(rf_importance),
  MeanDecreaseGini = rf_importance[, "MeanDecreaseGini"],
  MeanDecreaseAccuracy = rf_importance[, "MeanDecreaseAccuracy"]
)
rf_importance_df <- rf_importance_df[order(rf_importance_df$MeanDecreaseGini, decreasing = TRUE), ]

# Select top features (MeanDecreaseGini > median)
rf_genes <- as.character(rf_importance_df$Gene[
  rf_importance_df$MeanDecreaseGini > median(rf_importance_df$MeanDecreaseGini)
])
cat("RF selected genes:", length(rf_genes), "\n")

# Save
writeLines(rf_genes, file.path(tables_dir, "RF_selected_genes.txt"))
saveRDS(rf_fit, file.path(data_processed_dir, "rf_fit.rds"))

# Mean Decrease Gini plot
rf_importance_top <- rf_importance_df[1:min(20, nrow(rf_importance_df)), ]
rf_importance_top$Gene <- factor(rf_importance_top$Gene, 
                                  levels = rev(rf_importance_top$Gene))

p_rf_gini <- ggplot(rf_importance_top, aes(x = Gene, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "#4daf4a", alpha = 0.85) +
  coord_flip() +
  labs(title = "Random Forest: Mean Decrease Gini (Top 20)",
       x = "", y = "Mean Decrease Gini") +
  theme_publication
save_publication_figure(p_rf_gini, "Figure_6A_RF_MeanDecreaseGini")

# Mean Decrease Accuracy plot
p_rf_acc <- ggplot(rf_importance_top, aes(x = Gene, y = MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "#984ea3", alpha = 0.85) +
  coord_flip() +
  labs(title = "Random Forest: Mean Decrease Accuracy (Top 20)",
       x = "", y = "Mean Decrease Accuracy") +
  theme_publication
save_publication_figure(p_rf_acc, "Figure_6B_RF_MeanDecreaseAccuracy")

# Top feature importance bar plot
p_rf_top <- ggplot(rf_importance_top, aes(x = Gene, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = viridis(n = 20, option = "D"), alpha = 0.85) +
  coord_flip() +
  labs(title = "Top 20 Important Features (Random Forest)",
       x = "", y = "Mean Decrease Gini") +
  theme_publication
save_publication_figure(p_rf_top, "Figure_6C_RF_Top_Features", width = 8, height = 7)

# Save all ML results
ml_results <- list(
  lasso_genes = lasso_genes,
  svm_genes = svm_genes,
  rf_genes = rf_genes,
  rf_importance = rf_importance_df
)
saveRDS(ml_results, file.path(data_processed_dir, "ml_selection_results.rds"))

cat("\n========================================================\n")
cat("MACHINE LEARNING FEATURE SELECTION COMPLETE\n")
cat("========================================================\n")
cat("LASSO selected:", length(lasso_genes), "genes\n")
cat("SVM selected:", length(svm_genes), "genes\n")
cat("RF selected:", length(rf_genes), "genes\n")
cat("Results saved to:\n")
cat("  - results/figures/Machine_Learning/\n")
cat("  - results/tables/\n")
