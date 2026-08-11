# ============================================================
# Step 9: Machine Learning Feature Selection
# ============================================================
set.seed(123)
library(glmnet)
library(e1071)
library(randomForest)
library(caret)
library(ROCR)
library(pROC)
library(ggplot2)

base_dir <- normalizePath(file.path(".", "Breast Cancer 5 Necroptosis"))
output_dir <- base_dir

sink(file.path(output_dir, "logs", "09_machine_learning.log"), split=TRUE)
cat("MACHINE LEARNING FEATURE SELECTION\n")
cat("==================================\n\n")

save_tiff_600dpi <- function(plot, filename, width = 7, height = 6) {
  ggsave(filename=filename, plot=plot, width=width, height=height,
         units="in", dpi=600, compression="lzw", bg="white", device="tiff")
}

combat_expr <- readRDS(file.path(output_dir, "data_processed", "discovery_combat_corrected.rds"))
merged_meta <- readRDS(file.path(output_dir, "data_processed", "discovery_metadata.rds"))
de_nrgs <- readLines(file.path(output_dir, "results", "necroptosis_genes", "DE_NRGs.txt"))

# Use DE-NRGs as feature pool
features <- intersect(de_nrgs, rownames(combat_expr))
cat("DE-NRG features:", length(features), "\n")

if(length(features) < 5) {
  cat("WARNING: Too few DE-NRGs. Using top DEGs instead.\n")
  deg_sig <- read.csv(file.path(output_dir, "results", "tables", "DEG_significant.csv"), row.names=1)
  features <- intersect(rownames(deg_sig)[1:min(100, nrow(deg_sig))], rownames(combat_expr))
}

# Prepare data
X <- t(combat_expr[features, ])
y <- factor(merged_meta$class, levels=c("Normal", "Tumor"))
y_binary <- ifelse(y == "Tumor", 1, 0)

cat("Training samples:", nrow(X), "Features:", ncol(X), "\n")
cat("Class distribution:", table(y), "\n")

# Stratified 80/20 split
set.seed(123)
train_index <- createDataPartition(y, p=0.8, list=FALSE, times=1)
X_train <- X[train_index, ]
y_train <- y[train_index]
y_train_bin <- y_binary[train_index]
X_test <- X[-train_index, ]
y_test <- y[-train_index]
y_test_bin <- y_binary[-train_index]

cat("Train:", nrow(X_train), "Test:", nrow(X_test), "\n")
cat("Train class:", table(y_train), "\n")
cat("Test class:", table(y_test), "\n")

# Save indices
saveRDS(train_index, file.path(output_dir, "model_objects", "train_indices.rds"))

# Repeated CV indices (10-fold, 10 repeats)
set.seed(123)
cv_folds <- createMultiFolds(y_train, k=10, times=10)
saveRDS(cv_folds, file.path(output_dir, "model_objects", "cv_folds.rds"))

# Helper: plot ROC
plot_roc <- function(pred, truth, title) {
  roc_obj <- roc(truth, pred)
  df <- data.frame(FPR=1-roc_obj$specificities, TPR=roc_obj$sensitivities)
  auc_val <- auc(roc_obj)
  p <- ggplot(df, aes(x=FPR, y=TPR)) +
    geom_line(color="red", size=1) +
    geom_abline(intercept=0, slope=1, linetype="dashed", color="grey") +
    theme_bw() +
    labs(title=paste(title, "(AUC =", round(auc_val, 3), ")"),
         x="False Positive Rate", y="True Positive Rate") +
    annotate("text", x=0.75, y=0.25, label=paste("AUC =", round(auc_val, 3)), size=5)
  list(plot=p, auc=auc_val, roc=roc_obj)
}

# ============ LASSO ============
cat("\n============ LASSO LOGISTIC REGRESSION ============\n")
set.seed(123)

# Cross-validation for lambda
cv_lasso <- cv.glmnet(as.matrix(X_train), y_train_bin, family="binomial",
                       alpha=1, nfolds=10, type.measure="deviance",
                       parallel=FALSE)

saveRDS(cv_lasso, file.path(output_dir, "model_objects", "lasso_cv.rds"))

# Lambda plot
tiff(file.path(output_dir, "results", "machine_learning", "lasso_cv_error.tiff"),
     width=8, height=6, units="in", res=600, compression="lzw", bg="white")
plot(cv_lasso)
title("LASSO Cross-Validation Error", line=3)
dev.off()

# Selected features at lambda.1se
lambda_1se <- cv_lasso$lambda.1se
lasso_coef <- coef(cv_lasso, s=lambda_1se)
lasso_selected <- rownames(lasso_coef)[which(lasso_coef[,1] != 0)][-1]  # Remove intercept
cat("LASSO selected features at lambda.1se:", length(lasso_selected), "\n")
cat("  ", paste(lasso_selected, collapse=", "), "\n")

writeLines(lasso_selected, file.path(output_dir, "results", "machine_learning", "lasso_selected_genes.txt"))

# Coefficient plot
lasso_coef_full <- coef(cv_lasso, s=lambda_1se)
lasso_coef_df <- data.frame(
  Gene = rownames(lasso_coef_full),
  Coefficient = as.numeric(lasso_coef_full[,1])
)
lasso_coef_df <- lasso_coef_df[lasso_coef_df$Coefficient != 0 & lasso_coef_df$Gene != "(Intercept)", ]
lasso_coef_df <- lasso_coef_df[order(abs(lasso_coef_df$Coefficient), decreasing=TRUE), ]

p_lasso_coef <- ggplot(lasso_coef_df, aes(x=reorder(Gene, Coefficient), y=Coefficient, fill=Coefficient>0)) +
  geom_bar(stat="identity") + theme_bw() +
  coord_flip() + labs(title="LASSO Selected Coefficients", x="", y="Coefficient") +
  scale_fill_manual(values=c("TRUE"="red", "FALSE"="blue"), guide="none")
save_tiff_600dpi(p_lasso_coef, file.path(output_dir, "results", "machine_learning", "lasso_coefficients.tiff"), width=8, height=max(5, length(lasso_selected)*0.3))

# LASSO predictions
lasso_pred_prob <- predict(cv_lasso, as.matrix(X_test), s=lambda_1se, type="response")
lasso_pred_class <- ifelse(lasso_pred_prob > 0.5, "Tumor", "Normal")

# LASSO ROC
lasso_roc <- plot_roc(as.numeric(lasso_pred_prob), y_test, "LASSO")
save_tiff_600dpi(lasso_roc$plot, file.path(output_dir, "results", "machine_learning", "lasso_ROC.tiff"))

# LASSO confusion matrix
cm_lasso <- confusionMatrix(factor(lasso_pred_class, levels=c("Normal", "Tumor")), y_test)
cat("LASSO Confusion Matrix:\n")
print(cm_lasso$table)
write.csv(as.matrix(cm_lasso$table), file.path(output_dir, "results", "tables", "lasso_confusion_matrix.csv"))

# ============ SVM ============
cat("\n============ SVM ============\n")
set.seed(123)

# Use all DE-NRGs as features for SVM
svm_selected <- features
cat("SVM features:", length(svm_selected), "\n")
writeLines(svm_selected, file.path(output_dir, "results", "machine_learning", "svm_selected_genes.txt"))

# Cost tuning
set.seed(123)
svm_tune <- tune(svm, train.x=X_train[, svm_selected, drop=FALSE], train.y=y_train,
                 kernel="linear", ranges=list(cost=10^seq(-3, 2, length=10)),
                 tunecontrol=tune.control(sampling="cross", cross=10))
saveRDS(svm_tune, file.path(output_dir, "model_objects", "svm_tune.rds"))

tiff(file.path(output_dir, "results", "machine_learning", "svm_cost_tuning.tiff"),
     width=8, height=6, units="in", res=600, compression="lzw", bg="white")
plot(svm_tune)
title("SVM Cost Tuning", line=3)
dev.off()

# Best SVM model
svm_best_cost <- svm_tune$best.parameters$cost
svm_model <- svm(X_train[, svm_selected, drop=FALSE], y_train, kernel="linear",
                  cost=svm_best_cost, probability=TRUE)
saveRDS(svm_model, file.path(output_dir, "model_objects", "svm_final.rds"))

# SVM predictions
svm_pred_prob <- attr(predict(svm_model, X_test[, svm_selected, drop=FALSE], probability=TRUE), "probabilities")[, "Tumor"]
svm_pred_class <- predict(svm_model, X_test[, svm_selected, drop=FALSE])

# SVM ROC
svm_roc <- plot_roc(svm_pred_prob, y_test, "SVM")
save_tiff_600dpi(svm_roc$plot, file.path(output_dir, "results", "machine_learning", "svm_ROC.tiff"))

# SVM confusion matrix
cm_svm <- confusionMatrix(svm_pred_class, y_test)
cat("SVM Confusion Matrix:\n")
print(cm_svm$table)
write.csv(as.matrix(cm_svm$table), file.path(output_dir, "results", "tables", "svm_confusion_matrix.csv"))

# ============ Random Forest ============
cat("\n============ RANDOM FOREST ============\n")

# Tune mtry
set.seed(123)
tune_grid <- expand.grid(mtry = seq(2, max(5, floor(sqrt(ncol(X_train)))), by=2))
rf_cv <- train(X_train, y_train, method="rf", metric="Accuracy",
               trControl=trainControl(method="cv", number=10, verboseIter=FALSE),
               tuneGrid=tune_grid, ntree=1000)
saveRDS(rf_cv, file.path(output_dir, "model_objects", "rf_tune.rds"))

p_rf_tune <- ggplot(rf_cv$results, aes(x=mtry, y=Accuracy)) +
  geom_line(size=1) + geom_point(size=2) + theme_bw() +
  labs(title="Random Forest: mtry Tuning", x="mtry", y="Accuracy")
save_tiff_600dpi(p_rf_tune, file.path(output_dir, "results", "machine_learning", "rf_tuning.tiff"), width=7, height=5)

# Best RF model
rf_model <- randomForest(X_train, y_train, ntree=1000,
                          mtry=rf_cv$bestTune$mtry, importance=TRUE)
saveRDS(rf_model, file.path(output_dir, "model_objects", "rf_final.rds"))

# Variable importance
rf_imp <- importance(rf_model)
rf_imp_df <- data.frame(Gene=rownames(rf_imp), MeanDecreaseGini=rf_imp[, "MeanDecreaseGini"])
rf_imp_df <- rf_imp_df[order(rf_imp_df$MeanDecreaseGini, decreasing=TRUE), ]
write.csv(rf_imp_df, file.path(output_dir, "results", "machine_learning", "rf_importance.csv"), row.names=FALSE)

rf_selected <- rf_imp_df$Gene[1:min(50, nrow(rf_imp_df))]
writeLines(rf_selected, file.path(output_dir, "results", "machine_learning", "rf_selected_genes.txt"))

# Top importance plot
p_rf_imp <- ggplot(head(rf_imp_df, 20), aes(x=reorder(Gene, MeanDecreaseGini), y=MeanDecreaseGini)) +
  geom_bar(stat="identity", fill="forestgreen") + theme_bw() +
  coord_flip() + labs(title="Random Forest: Top 20 Important Features", x="", y="Mean Decrease Gini")
save_tiff_600dpi(p_rf_imp, file.path(output_dir, "results", "machine_learning", "rf_importance_plot.tiff"), width=8, height=6)

# RF OOB error plot
tiff(file.path(output_dir, "results", "machine_learning", "rf_oob_error.tiff"),
     width=8, height=6, units="in", res=600, compression="lzw", bg="white")
plot(rf_model, main="Random Forest OOB Error")
legend("topright", legend=colnames(rf_model$err.rate), col=1:3, lty=1, cex=0.8)
dev.off()

# RF predictions
rf_pred_prob <- predict(rf_model, X_test, type="prob")[, "Tumor"]
rf_pred_class <- predict(rf_model, X_test)

# RF ROC
rf_roc <- plot_roc(rf_pred_prob, y_test, "Random Forest")
save_tiff_600dpi(rf_roc$plot, file.path(output_dir, "results", "machine_learning", "rf_ROC.tiff"))

# RF confusion matrix
cm_rf <- confusionMatrix(rf_pred_class, y_test)
cat("RF Confusion Matrix:\n")
print(cm_rf$table)
write.csv(as.matrix(cm_rf$table), file.path(output_dir, "results", "tables", "rf_confusion_matrix.csv"))

# Save all predictions
predictions_df <- data.frame(
  Sample = rownames(X_test),
  True = y_test,
  LASSO_prob = as.numeric(lasso_pred_prob),
  LASSO_class = lasso_pred_class,
  SVM_prob = svm_pred_prob,
  SVM_class = svm_pred_class,
  RF_prob = rf_pred_prob,
  RF_class = rf_pred_class
)
write.csv(predictions_df, file.path(output_dir, "results", "machine_learning", "all_predictions.csv"), row.names=FALSE)

# Performance comparison
perf_df <- data.frame(
  Method = c("LASSO", "SVM", "RF"),
  AUC = c(lasso_roc$auc, svm_roc$auc, rf_roc$auc),
  Accuracy = c(cm_lasso$overall["Accuracy"], cm_svm$overall["Accuracy"], cm_rf$overall["Accuracy"]),
  Sensitivity = c(cm_lasso$byClass["Sensitivity"], cm_svm$byClass["Sensitivity"], cm_rf$byClass["Sensitivity"]),
  Specificity = c(cm_lasso$byClass["Specificity"], cm_svm$byClass["Specificity"], cm_rf$byClass["Specificity"])
)
write.csv(perf_df, file.path(output_dir, "results", "tables", "ML_performance_comparison.csv"), row.names=FALSE)
cat("\nML Performance Summary:\n")
print(perf_df)

cat("\n\nMACHINE LEARNING COMPLETE\n")
sink()
