# ============================================================
# Step 11: Predictive Model and SHAP Interpretation
# ============================================================
set.seed(123)
library(glmnet)
library(e1071)
library(randomForest)
library(caret)
library(pROC)
library(PRROC)
library(ggplot2)

base_dir <- normalizePath(file.path(".", "Breast Cancer 5 Necroptosis"))
output_dir <- base_dir

sink(file.path(output_dir, "logs", "11_predictive_model.log"), split=TRUE)
cat("PREDICTIVE MODEL AND SHAP\n")
cat("=========================\n\n")

save_tiff_600dpi <- function(plot, filename, width = 7, height = 6) {
  ggsave(filename=filename, plot=plot, width=width, height=height,
         units="in", dpi=600, compression="lzw", bg="white", device="tiff")
}

combat_expr <- readRDS(file.path(output_dir, "data_processed", "discovery_combat_corrected.rds"))
merged_meta <- readRDS(file.path(output_dir, "data_processed", "discovery_metadata.rds"))
final_genes <- readLines(file.path(output_dir, "results", "consensus_genes", "final_consensus_genes.txt"))
train_index <- readRDS(file.path(output_dir, "model_objects", "train_indices.rds"))

final_genes <- intersect(final_genes, rownames(combat_expr))
cat("Final consensus genes available:", length(final_genes), "\n")

if(length(final_genes) < 3) {
  cat("WARNING: Too few consensus genes. Using DE-NRGs instead.\n")
  de_nrgs <- readLines(file.path(output_dir, "results", "necroptosis_genes", "DE_NRGs.txt"))
  final_genes <- intersect(de_nrgs, rownames(combat_expr))
}

X <- t(combat_expr[final_genes, ])
y <- factor(merged_meta$class, levels=c("Normal", "Tumor"))

X_train <- X[train_index, ]
y_train <- y[train_index]
X_test <- X[-train_index, ]
y_test <- y[-train_index]

cat("Train:", nrow(X_train), "Test:", nrow(X_test), "\n")

# ============ Penalized Logistic Regression ============
cat("\n--- Penalized Logistic Regression ---\n")
set.seed(123)
cv_plr <- cv.glmnet(as.matrix(X_train), y_train, family="binomial", alpha=0.5,
                     nfolds=10, type.measure="deviance")
plr_model <- glmnet(as.matrix(X_train), y_train, family="binomial", alpha=0.5,
                     lambda=cv_plr$lambda.1se)
saveRDS(plr_model, file.path(output_dir, "model_objects", "plr_final.rds"))

plr_prob <- predict(plr_model, as.matrix(X_test), type="response")[,1]
plr_class <- ifelse(plr_prob > 0.5, "Tumor", "Normal")

# ============ SVM ============
cat("\n--- SVM ---\n")
set.seed(123)
svm_tune <- tune(svm, X_train, y_train, kernel="linear",
                 ranges=list(cost=10^seq(-3, 1, length=5)),
                 tunecontrol=tune.control(sampling="cross", cross=10))
svm_final <- svm(X_train, y_train, kernel="linear", 
                  cost=svm_tune$best.parameters$cost, probability=TRUE)
saveRDS(svm_final, file.path(output_dir, "model_objects", "svm_final_model.rds"))

svm_prob <- attr(predict(svm_final, X_test, probability=TRUE), "probabilities")[, "Tumor"]
svm_class <- predict(svm_final, X_test)

# ============ Random Forest ============
cat("\n--- Random Forest ---\n")
set.seed(123)
rf_final <- randomForest(X_train, y_train, ntree=1000,
                          mtry=max(1, floor(sqrt(ncol(X_train)))), importance=TRUE)
saveRDS(rf_final, file.path(output_dir, "model_objects", "rf_final_model.rds"))

rf_prob <- predict(rf_final, X_test, type="prob")[, "Tumor"]
rf_class <- predict(rf_final, X_test)

# ============ Evaluation ============
eval_model <- function(prob, pred_class, true_class, model_name) {
  roc_obj <- roc(true_class, prob)
  pr_obj <- pr.curve(scores.class0=prob[true_class=="Tumor"], 
                     scores.class1=prob[true_class=="Normal"], curve=TRUE)
  
  cm <- confusionMatrix(factor(pred_class, levels=c("Normal", "Tumor")), true_class)
  
  # Brier score
  y_bin <- ifelse(true_class == "Tumor", 1, 0)
  brier <- mean((prob - y_bin)^2)
  
  # Calibration
  cal_df <- data.frame(pred=prob, obs=y_bin)
  cal_loess <- loess(obs ~ pred, data=cal_df, degree=1)
  
  list(
    Model = model_name,
    AUC = as.numeric(auc(roc_obj)),
    AUPRC = pr_obj$auc.integral,
    Accuracy = cm$overall["Accuracy"],
    BalancedAccuracy = cm$byClass["Balanced Accuracy"],
    Sensitivity = cm$byClass["Sensitivity"],
    Specificity = cm$byClass["Specificity"],
    PPV = cm$byClass["Pos Pred Value"],
    NPV = cm$byClass["Neg Pred Value"],
    F1 = cm$byClass["F1"],
    MCC = { tp <- sum(pred_class=="Tumor" & true_class=="Tumor")
            tn <- sum(pred_class=="Normal" & true_class=="Normal")
            fp <- sum(pred_class=="Tumor" & true_class=="Normal")
            fn <- sum(pred_class=="Normal" & true_class=="Tumor")
            denom <- sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
            if(denom == 0) 0 else (tp*tn - fp*fn) / denom },
    Brier = brier,
    ROC = roc_obj,
    PR = pr_obj,
    Confusion = cm$table,
    Calibration = cal_df
  )
}

results <- list(
  PLR = eval_model(plr_prob, plr_class, y_test, "Penalized Logistic Regression"),
  SVM = eval_model(svm_prob, svm_class, y_test, "SVM"),
  RF = eval_model(rf_prob, rf_class, y_test, "Random Forest")
)

metrics_df <- do.call(rbind, lapply(results, function(r) {
  data.frame(Model=r$Model, AUC=r$AUC, AUPRC=r$AUPRC, Accuracy=r$Accuracy,
             BalancedAccuracy=r$BalancedAccuracy, Sensitivity=r$Sensitivity,
             Specificity=r$Specificity, PPV=r$PPV, NPV=r$NPV, F1=r$F1,
             MCC=ifelse(is.na(r$MCC), 0, r$MCC), Brier=r$Brier)
}))
rownames(metrics_df) <- NULL
cat("\nModel Performance:\n")
print(metrics_df)
write.csv(metrics_df, file.path(output_dir, "results", "predictive_model", "model_comparison.csv"), row.names=FALSE)

# Combined ROC curves
roc_df <- data.frame()
for(nm in names(results)) {
  r <- results[[nm]]
  df <- data.frame(FPR=1-r$ROC$specificities, TPR=r$ROC$sensitivities, Model=nm)
  roc_df <- rbind(roc_df, df)
}
p_roc <- ggplot(roc_df, aes(x=FPR, y=TPR, color=Model)) +
  geom_line(size=1) + theme_bw() +
  geom_abline(intercept=0, slope=1, linetype="dashed", color="grey") +
  labs(title="ROC Curves - Hold-out Test Set", x="False Positive Rate", y="True Positive Rate") +
  annotate("text", x=0.7, y=0.3, hjust=0,
           label=paste0("PLR AUC=", round(results$PLR$AUC,3), "\n",
                        "SVM AUC=", round(results$SVM$AUC,3), "\n",
                        "RF AUC=", round(results$RF$AUC,3)))
save_tiff_600dpi(p_roc, file.path(output_dir, "results", "predictive_model", "combined_ROC.tiff"), width=8, height=7)

# PR curves
pr_df <- data.frame()
for(nm in names(results)) {
  r <- results[[nm]]
  df <- data.frame(Recall=r$PR$curve[,1], Precision=r$PR$curve[,2], Model=nm)
  pr_df <- rbind(pr_df, df)
}
p_pr <- ggplot(pr_df, aes(x=Recall, y=Precision, color=Model)) +
  geom_line(size=1) + theme_bw() +
  labs(title="Precision-Recall Curves", x="Recall", y="Precision") +
  annotate("text", x=0.7, y=0.9, hjust=0,
           label=paste0("PLR AUPRC=", round(results$PLR$AUPRC,3), "\n",
                        "SVM AUPRC=", round(results$SVM$AUPRC,3), "\n",
                        "RF AUPRC=", round(results$RF$AUPRC,3)))
save_tiff_600dpi(p_pr, file.path(output_dir, "results", "predictive_model", "PR_curves.tiff"), width=8, height=7)

# Metric comparison plot
metrics_long <- reshape2::melt(metrics_df[, c("Model", "AUC", "Accuracy", "F1", "MCC", "Brier")], id="Model")
p_metrics <- ggplot(metrics_long, aes(x=Model, y=value, fill=Model)) +
  geom_bar(stat="identity", position="dodge") + theme_bw() +
  facet_wrap(~variable, scales="free_y", nrow=2) +
  labs(title="Model Performance Metrics", x="", y="Value") +
  theme(axis.text.x=element_text(angle=45, hjust=1))
save_tiff_600dpi(p_metrics, file.path(output_dir, "results", "predictive_model", "metrics_comparison.tiff"), width=10, height=8)

# Confusion matrices
for(nm in names(results)) {
  cm_df <- as.data.frame(results[[nm]]$Confusion)
  write.csv(cm_df, file.path(output_dir, "results", "tables", paste0("cm_", nm, ".csv")))
}

# Calibration plots
cal_df <- data.frame()
for(nm in names(results)) {
  r <- results[[nm]]
  cal_df <- rbind(cal_df, data.frame(Predicted=r$Calibration$pred, Observed=r$Calibration$obs, Model=nm))
}
p_cal <- ggplot(cal_df, aes(x=Predicted, y=Observed, color=Model)) +
  geom_point(alpha=0.3) + geom_smooth(method="loess", se=FALSE) +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  theme_bw() + labs(title="Calibration Plot", x="Predicted Probability", y="Observed Frequency")
save_tiff_600dpi(p_cal, file.path(output_dir, "results", "predictive_model", "calibration_plot.tiff"), width=8, height=7)

# Decision curve analysis
dca_df <- data.frame()
for(nm in names(results)) {
  prob <- if(nm=="PLR") plr_prob else if(nm=="SVM") svm_prob else rf_prob
  y_bin <- ifelse(y_test == "Tumor", 1, 0)
  
  thresholds <- seq(0, 1, by=0.01)
  net_benefit <- sapply(thresholds, function(t) {
    tp <- sum(prob >= t & y_bin == 1)
    fp <- sum(prob >= t & y_bin == 0)
    n <- length(y_bin)
    tp/n - fp/n * (t/(1-t))
  })
  dca_df <- rbind(dca_df, data.frame(Threshold=thresholds, NetBenefit=net_benefit, Model=nm))
}
# Add treat all / treat none
y_bin <- ifelse(y_test == "Tumor", 1, 0)
prevalence <- mean(y_bin)
dca_df <- rbind(dca_df, data.frame(Threshold=thresholds, NetBenefit=prevalence - (1-prevalence)*thresholds/(1-thresholds), Model="Treat All"))
dca_df <- rbind(dca_df, data.frame(Threshold=thresholds, NetBenefit=0, Model="Treat None"))

p_dca <- ggplot(dca_df, aes(x=Threshold, y=NetBenefit, color=Model)) +
  geom_line(size=1) + theme_bw() +
  labs(title="Decision Curve Analysis", x="Threshold Probability", y="Net Benefit") +
  ylim(c(min(dca_df$NetBenefit), max(dca_df$NetBenefit)))
save_tiff_600dpi(p_dca, file.path(output_dir, "results", "predictive_model", "decision_curve.tiff"), width=8, height=7)

# ============ SHAP Analysis ============
cat("\n--- SHAP Analysis ---\n")

# Use the best model (by AUC)
best_model_name <- names(which.max(sapply(results, function(r) r$AUC)))
cat("Best model for SHAP:", best_model_name, "\n")

# Model-agnostic SHAP with fastshap (fallback to coefficient proxy)
if(require(fastshap, quietly=TRUE)) {
  shap_func <- function(object, newdata) {
    if(inherits(object, "randomForest")) {
      predict(object, newdata, type="prob")[, "Tumor"]
    } else if(inherits(object, "svm")) {
      attr(predict(object, newdata, probability=TRUE), "probabilities")[, "Tumor"]
    } else {
      predict(object, newdata, type="response")
    }
  }
  
  best_model <- switch(best_model_name,
    "RF" = rf_final,
    "SVM" = svm_final,
    "PLR" = plr_model
  )
  
  cat("Computing SHAP values with fastshap (nsim=50)...\n")
  shap_values <- fastshap::explain(best_model, X=X_train, pred_wrapper=shap_func, nsim=50)
  shap_values <- matrix(as.numeric(shap_values), nrow=nrow(shap_values))
  colnames(shap_values) <- colnames(X_train)
  rownames(shap_values) <- rownames(X_train)
  saveRDS(shap_values, file.path(output_dir, "results", "SHAP", "shap_values.rds"))
  
  # SHAP importance
  shap_imp <- data.frame(Gene=colnames(shap_values), Importance=colMeans(abs(shap_values)))
  shap_imp <- shap_imp[order(shap_imp$Importance, decreasing=TRUE), ]
  write.csv(shap_imp, file.path(output_dir, "results", "SHAP", "shap_importance.csv"), row.names=FALSE)
  
  p_shap_bar <- ggplot(head(shap_imp, 20), aes(x=reorder(Gene, Importance), y=Importance)) +
    geom_bar(stat="identity", fill="steelblue") + theme_bw() + coord_flip() +
    labs(title="SHAP Feature Importance", x="", y="Mean |SHAP Value|")
  save_tiff_600dpi(p_shap_bar, file.path(output_dir, "results", "SHAP", "shap_importance_bar.tiff"), width=8, height=6)
  
  # SHAP beeswarm
  tiff(file.path(output_dir, "results", "SHAP", "shap_beeswarm.tiff"),
       width=10, height=8, units="in", res=600, compression="lzw", bg="white")
  top_n <- min(20, ncol(shap_values))
  top_feats <- shap_imp$Gene[1:top_n]
  plot_df <- do.call(rbind, lapply(top_feats, function(g) {
    data.frame(Gene=g, SHAP=shap_values[, g], FeatureValue=X_train[, g])
  }))
  plot_df$Gene <- factor(plot_df$Gene, levels=rev(top_feats))
  print(ggplot(plot_df, aes(x=SHAP, y=Gene, color=FeatureValue)) +
    geom_jitter(width=0, height=0.15, alpha=0.6, size=0.8) +
    scale_color_gradient2(low="blue", mid="purple", high="red", midpoint=median(X_train)) +
    theme_bw() + labs(title="SHAP Beeswarm Plot", x="SHAP Value", y="", color="Feature\nValue"))
  dev.off()
  
  # SHAP dependence plots (top 5)
  for(gene in shap_imp$Gene[1:min(5, nrow(shap_imp))]) {
    dp_df <- data.frame(Feature=X_train[, gene], SHAP=shap_values[, gene])
    p_dp <- ggplot(dp_df, aes(x=Feature, y=SHAP)) +
      geom_point(alpha=0.5, size=0.8) + geom_smooth(method="loess", se=TRUE, color="red") +
      theme_bw() + labs(title=paste("SHAP Dependence:", gene), x=gene, y="SHAP Value")
    save_tiff_600dpi(p_dp, file.path(output_dir, "results", "SHAP", paste0("shap_dependence_", gene, ".tiff")), width=7, height=6)
  }
  
  # SHAP waterfall
  sample_idx <- which(y_train == "Tumor")[1]
  shap_waterfall <- data.frame(Gene=colnames(shap_values), SHAP=as.numeric(shap_values[sample_idx, ]))
  shap_waterfall <- shap_waterfall[order(abs(shap_waterfall$SHAP), decreasing=TRUE), ]
  shap_waterfall$Gene <- factor(shap_waterfall$Gene, levels=rev(shap_waterfall$Gene))
  p_waterfall <- ggplot(head(shap_waterfall, 15), aes(x=SHAP, y=Gene)) +
    geom_bar(stat="identity", aes(fill=SHAP>0)) + theme_bw() +
    scale_fill_manual(values=c("TRUE"="red", "FALSE"="blue"), guide="none") +
    labs(title=paste("SHAP Waterfall - Sample", sample_idx), x="SHAP Value", y="")
  save_tiff_600dpi(p_waterfall, file.path(output_dir, "results", "SHAP", "shap_waterfall.tiff"), width=8, height=6)
  
  # SHAP heatmap
  tiff(file.path(output_dir, "results", "SHAP", "shap_heatmap.tiff"),
       width=10, height=8, units="in", res=600, compression="lzw", bg="white")
  pheatmap::pheatmap(t(shap_values[, head(shap_imp$Gene, 30)]),
                     main="SHAP Values Heatmap",
                     show_rownames=TRUE, show_colnames=FALSE,
                     color=colorRampPalette(c("navy", "white", "firebrick3"))(100),
                     fontsize_row=8)
  dev.off()
  
  cat("SHAP analysis complete.\n")
} else {
  cat("fastshap not installed. Using coefficient-based importance proxy.\n")
  if(best_model_name == "RF") {
    imp_df <- data.frame(Gene=rownames(rf_final$importance), 
                         Importance=rf_final$importance[, "MeanDecreaseGini"])
  } else if(best_model_name == "SVM") {
    if(require(kernlab, quietly=TRUE)) {
      coef_weights <- t(svm_final$coefs) %*% svm_final$SV
      imp_df <- data.frame(Gene=colnames(coef_weights), Importance=abs(as.numeric(coef_weights)))
    } else {
      imp_df <- data.frame(Gene=colnames(X_train), Importance=rep(1, ncol(X_train)))
    }
  } else {
    coef_mat <- as.matrix(coef(plr_model))
    imp_df <- data.frame(Gene=rownames(coef_mat)[-1], Importance=abs(coef_mat[-1,1]))
  }
  imp_df <- imp_df[order(imp_df$Importance, decreasing=TRUE), ]
  write.csv(imp_df, file.path(output_dir, "results", "SHAP", "shap_importance.csv"), row.names=FALSE)
}

cat("\n\nPREDICTIVE MODEL AND SHAP COMPLETE\n")
sink()
