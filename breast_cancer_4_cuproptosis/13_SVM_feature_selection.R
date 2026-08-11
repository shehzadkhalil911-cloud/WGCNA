# 13_SVM_feature_selection.R
source("config/analysis_config.R")
source("R/load_functions.R")
ensure_paths()
suppressMessages({ library(caret); library(e1071); library(pROC); library(ggplot2) })

combat <- readRDS(file.path(paths$discovery, "discovery_combat.rds"))
meta <- readRDS(file.path(paths$discovery, "discovery_cohort_raw.rds"))$metadata
crg_sum <- readRDS(file.path(paths$res_cross, "DE_CRG_summary.rds"))
pool <- crg_sum$candidate_pool
x <- t(combat[pool, ])
y <- factor(meta$class, levels = c("Normal", "Tumor"))

split <- readRDS(file.path(paths$res_ml, "LASSO", "train_test_split.rds"))
train_idx <- split$train; test_idx <- split$test

ctrl <- rfeControl(functions = lrFuncs, method = "repeatedcv", repeats = 5, number = CV_FOLDS,
                   seeds = as.list(rep(list(1:5), 56)), returnResamp = "final")
set.seed(GLOBAL_SEED)
rfe_res <- rfe(x[train_idx, ], y[train_idx], sizes = c(2,4,6,8,10),
               rfeControl = ctrl)
best <- predictors(rfe_res)
svm_fit <- svm(x[train_idx, best, drop = FALSE], y[train_idx], kernel = "linear", probability = TRUE)
pred_test <- attr(predict(svm_fit, x[test_idx, best, drop = FALSE], probability = TRUE), "probabilities")[, "Tumor"]
roc_obj <- roc(y[test_idx], pred_test)
auc_val <- auc(roc_obj)

res <- list(rfe = rfe_res, selected = best, svm_fit = svm_fit,
            pred_test = pred_test, y_test = y[test_idx], auc_test = auc_val)
saveRDS(res, file.path(paths$res_ml, "SVM", "svm_result.rds"))
write.csv(data.frame(gene = best, rank = seq_along(best)),
          file.path(paths$res_ml, "SVM", "selected_genes.csv"), row.names = FALSE)

## plots
tiff(file.path(paths$res_ml, "SVM", "error_vs_features.tiff"), width = 7, height = 5,
     units = "in", res = 600, compression = "lzw")
print(plot(rfe_res, type = c("g", "o")))
dev.off()
pdf(file.path(paths$res_ml, "SVM", "error_vs_features.pdf"), width = 7, height = 5)
print(plot(rfe_res, type = c("g", "o")))
dev.off()
cat("SVM-RFE complete. Selected:", paste(best, collapse = ", "), " Test AUC:", round(auc_val,3), "\n")
quit(save = "no", status = 0, runLast = FALSE)
