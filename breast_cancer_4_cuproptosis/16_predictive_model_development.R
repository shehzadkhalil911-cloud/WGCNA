# 16_predictive_model_development.R
source("config/analysis_config.R")
source("R/load_functions.R")
ensure_paths()
suppressMessages({ library(caret); library(glmnet); library(e1071); library(ranger); library(pROC); library(PRROC); library(ggplot2); library(cowplot) })

combat <- readRDS(file.path(paths$discovery, "discovery_combat.rds"))
meta <- readRDS(file.path(paths$discovery, "discovery_cohort_raw.rds"))$metadata
cons <- readRDS(file.path(paths$res_ml, "consensus", "consensus_genes.rds"))
final_genes <- cons$final_genes
stopifnot(length(final_genes) >= 2)

x <- t(combat[final_genes, ])
y <- factor(meta$class, levels = c("Normal", "Tumor"))
split <- readRDS(file.path(paths$res_ml, "LASSO", "train_test_split.rds"))
tr <- split$train; te <- split$test

metrics_fun <- function(prob, obs, pos = "Tumor") {
  pred <- factor(ifelse(prob >= 0.5, pos, "Normal"), levels = c("Normal", pos))
  cm <- confusionMatrix(pred, obs, positive = pos)
  roc_obj <- roc(obs, prob, levels = c("Normal", pos))
  auc_v <- auc(roc_obj)
  tp <- cm$table[pos, pos]; fp <- cm$table["Normal", pos]
  fn <- cm$table[pos, "Normal"]; tn <- cm$table["Normal", "Normal"]
  precision_v <- tp / (tp + fp); recall_v <- tp / (tp + fn)
  f1_v <- 2 * precision_v * recall_v / (precision_v + recall_v)
  ## MCC
  denom <- sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  mcc_v <- ifelse(denom == 0, 0, (tp*tn - fp*fn)/denom)
  kappa_v <- cm$overall["Kappa"]
  list(accuracy = cm$overall["Accuracy"], bal_acc = cm$byClass["Balanced Accuracy"],
       sensitivity = cm$byClass["Sensitivity"], specificity = cm$byClass["Specificity"],
       ppv = cm$byClass["Pos Pred Value"], npv = cm$byClass["Neg Pred Value"],
       precision = precision_v, recall = recall_v, f1 = f1_v, mcc = mcc_v, kappa = kappa_v,
       roc_auc = auc_v,
       precision_recall_auc = pr.curve(scores.class0 = prob[obs == pos],
                                       scores.class1 = prob[obs != pos])$auc.integral,
       brier = mean((prob - (as.numeric(obs == pos)))^2))
}

## fit 4 models on training
set.seed(GLOBAL_SEED)
## logistic / LASSO (glmnet binomial)
fit_glmnet <- glmnet(as.matrix(x[tr,]), y[tr], family = "binomial", alpha = 1)
cvfit <- cv.glmnet(as.matrix(x[tr,]), y[tr], family = "binomial", alpha = 1, type.measure = "auc", nfolds = CV_FOLDS)
pred_lasso <- predict(cvfit, newx = as.matrix(x[te,]), s = "lambda.min", type = "response")[,1]
## logistic (no penalty via glmnet alpha=0)
pred_logit <- predict(glmnet(as.matrix(x[tr,]), y[tr], family = "binomial", alpha = 0),
                      newx = as.matrix(x[te,]), s = 0, type = "response")[,1]
## SVM
svm_fit <- svm(x[tr,], y[tr], kernel = "linear", probability = TRUE)
pred_svm <- attr(predict(svm_fit, x[te,], probability = TRUE), "probabilities")[, "Tumor"]
## RF
tr_df <- data.frame(x[tr,,drop=FALSE], class = y[tr], stringsAsFactors = FALSE)
rf_fit <- ranger(class ~ ., data = tr_df, num.trees = 1000, probability = TRUE, seed = GLOBAL_SEED)
pred_rf <- predict(rf_fit, data.frame(x[te,,drop=FALSE]))$predictions[, "Tumor"]

preds <- list(LASSO = pred_lasso, Logistic = pred_logit, SVM = pred_svm, RF = pred_rf)
mtab <- data.frame()
metric_names <- c("accuracy","bal_acc","sensitivity","specificity","ppv","npv","precision","recall",
                 "f1","mcc","kappa","roc_auc","pr_auc","brier")
for (nm in names(preds)) {
  m <- metrics_fun(preds[[nm]], y[te])
  mtab <- rbind(mtab, data.frame(model = nm, t(as.numeric(m))))
}
names(mtab) <- c("model", metric_names)
write.csv(mtab, file.path(paths$res_model, "model_performance.csv"), row.names = FALSE)

## best model by test ROC AUC
best <- names(preds)[which.max(sapply(preds, function(p) metrics_fun(p, y[te])$roc_auc))]
cat("Best model by hold-out ROC AUC:", best, "\n")
saveRDS(list(preds = preds, metrics = mtab, best_model = best,
             fits = list(glmnet = fit_glmnet, cvfit = cvfit, svm = svm_fit, rf = rf_fit),
             final_genes = final_genes, test_idx = te, y_test = y[te]),
        file.path(paths$res_model, "predictive_models.rds"))

## combined ROC
roc_df <- data.frame()
for (nm in names(preds)) {
  roc_df <- rbind(roc_df, data.frame(model = nm,
    fpr = 1 - rev(roc(y[te], preds[[nm]], levels = c("Normal","Tumor"))$spec),
    tpr = rev(roc(y[te], preds[[nm]], levels = c("Normal","Tumor"))$sen)))
}
p_roc <- ggplot(roc_df, aes(x = fpr, y = tpr, colour = model)) + geom_line(size = 1) +
  geom_abline(linetype = 2, colour = "grey") + labs(x = "1-Specificity", y = "Sensitivity", title = "Combined ROC curves") + pub_theme()
save_publication_plot(p_roc, file.path(paths$res_model, "combined_ROC"), width = 7, height = 6)

cat("Predictive model development complete.\n")
print(mtab)
quit(save = "no", status = 0, runLast = FALSE)
