# 12_LASSO_feature_selection.R
source("config/analysis_config.R")
source("R/load_functions.R")
ensure_paths()
suppressMessages({ library(ggplot2); library(glmnet); library(pROC); library(caret); library(doParallel) })

combat <- readRDS(file.path(paths$discovery, "discovery_combat.rds"))
meta <- readRDS(file.path(paths$discovery, "discovery_cohort_raw.rds"))$metadata
crg_sum <- readRDS(file.path(paths$res_cross, "DE_CRG_summary.rds"))
pool <- crg_sum$candidate_pool

## feature matrix = candidate CRG pool
x <- t(combat[pool, ])
y <- factor(meta$class, levels = c("Normal", "Tumor"))
cat("LASSO pool:", length(pool), "samples:", nrow(x), "\n")

## shared train/test split (saved for all algorithms)
split <- shared_split(y, train_prop = TRAIN_PROPORTION, seed = GLOBAL_SEED)
saveRDS(split, file.path(paths$res_ml, "LASSO", "train_test_split.rds"))
train_idx <- split$train; test_idx <- split$test

set.seed(GLOBAL_SEED)
cvfit <- cv.glmnet(as.matrix(x[train_idx, ]), y[train_idx], family = "binomial",
                   alpha = 1, type.measure = "auc", nfolds = CV_FOLDS)
fit_min <- glmnet(as.matrix(x[train_idx, ]), y[train_idx], family = "binomial", alpha = 1, lambda = cvfit$lambda.min)
fit_1se <- glmnet(as.matrix(x[train_idx, ]), y[train_idx], family = "binomial", alpha = 1, lambda = cvfit$lambda.1se)
coef_min <- as.matrix(coef(fit_min)); coef_1se <- as.matrix(coef(fit_1se))
sel_min <- rownames(coef_min)[coef_min[,1] != 0][-1]
sel_1se <- rownames(coef_1se)[coef_1se[,1] != 0][-1]

pred_test <- predict(cvfit, newx = as.matrix(x[test_idx, ]), s = "lambda.min", type = "response")[,1]
roc_obj <- roc(y[test_idx], pred_test)
auc_val <- auc(roc_obj)
pred_train <- predict(cvfit, newx = as.matrix(x[train_idx, ]), s = "lambda.min", type = "response")[,1]
roc_tr <- roc(y[train_idx], pred_train)
auc_tr <- auc(roc_tr)

res <- list(cvfit = cvfit, fit_min = fit_min, fit_1se = fit_1se,
            selected_min = sel_min, selected_1se = sel_1se,
            coef_min = coef_min, coef_1se = coef_1se,
            pred_test = pred_test, y_test = y[test_idx],
            lambda_min = cvfit$lambda.min, lambda_1se = cvfit$lambda.1se,
            auc_train = auc_tr, auc_test = auc_val)
saveRDS(res, file.path(paths$res_ml, "LASSO", "lasso_result.rds"))
write.csv(data.frame(gene = sel_min, coefficient = coef_min[sel_min,1]),
          file.path(paths$res_ml, "LASSO", "selected_genes_lambda_min.csv"), row.names = FALSE)
write.csv(data.frame(gene = sel_1se, coefficient = coef_1se[sel_1se,1]),
          file.path(paths$res_ml, "LASSO", "selected_genes_lambda_1se.csv"), row.names = FALSE)

## plots
pdf(file.path(paths$res_ml, "LASSO", "coefficient_profile.pdf"), width = 7, height = 5)
plot(cvfit$glmnet.fit, xvar = "lambda", label = TRUE); title("LASSO coefficient profile")
dev.off()
tiff(file.path(paths$res_ml, "LASSO", "coefficient_profile.tiff"), width = 7, height = 5,
     units = "in", res = 600, compression = "lzw")
plot(cvfit$glmnet.fit, xvar = "lambda", label = TRUE); title("LASSO coefficient profile")
dev.off()
pdf(file.path(paths$res_ml, "LASSO", "cv_plot.pdf"), width = 7, height = 5)
plot(cvfit); abline(v = log(cvfit$lambda.min), col = "red", lty = 2)
abline(v = log(cvfit$lambda.1se), col = "blue", lty = 2); title("LASSO CV error")
dev.off()
tiff(file.path(paths$res_ml, "LASSO", "cv_plot.tiff"), width = 7, height = 5,
     units = "in", res = 600, compression = "lzw")
plot(cvfit); abline(v = log(cvfit$lambda.min), col = "red", lty = 2)
abline(v = log(cvfit$lambda.1se), col = "blue", lty = 2); title("LASSO CV error")
dev.off()
cat("LASSO complete. Selected (lambda.min):", length(sel_min), " Selected (lambda.1se):", length(sel_1se),
    " Test AUC:", round(auc_val,3), "\n")
quit(save = "no", status = 0, runLast = FALSE)
