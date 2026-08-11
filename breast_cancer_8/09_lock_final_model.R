# ============================================================
# 09 - Lock final model
# Aggregates nested-CV tuning to choose final hyperparameters,
# refits XGBoost on all development data with the final stable
# features, attaches the locked Platt calibration + Youden
# threshold, and stores a self-contained prediction object.
# ============================================================
PROJECT_DIR <- "F:/independent projects/cancer/breast cancer/Breast cancer 8"
source(file.path(PROJECT_DIR, "config/config.R"))
source(file.path(PROJECT_DIR, "R/utils.R"))
suppressMessages(library(dplyr))

logmsg <- function(...) cat(sprintf(...), "\n")

nc  <- readRDS(file.path(DIR_MODELS, "nested_cv_results.rds"))
cal <- readRDS(file.path(DIR_MODELS, "calibration.rds"))
disc <- load_discovery()
split <- readRDS(file.path(DIR_MODELS, "primary_split.rds"))

genes <- nc$final_stable
expr_tr <- disc$expr[genes, split$train_samples, drop = FALSE]
meta_tr <- disc$meta[match(split$train_samples, disc$meta$sample_id), ]
y_tr <- as.numeric(meta_tr$class == "Tumor")
logmsg("Development: %d samples, %d genes", ncol(expr_tr), nrow(expr_tr))

# ---------------- choose final hyperparameters ----------------
ta <- nc$tune_all
agg <- aggregate(best_auc ~ eta + max_depth + min_child_weight + subsample +
                   colsample_bytree + gamma + alpha + lambda, ta, median)
agg <- agg[order(-agg$best_auc), ]
best_cfg <- agg[1, ]
prm <- list(
  eta = best_cfg$eta, max_depth = best_cfg$max_depth,
  min_child_weight = best_cfg$min_child_weight,
  subsample = best_cfg$subsample, colsample_bytree = best_cfg$colsample_bytree,
  gamma = best_cfg$gamma, alpha = best_cfg$alpha, lambda = best_cfg$lambda
)

# nrounds: median best_iter across outer folds for the chosen config
sel_rows <- which(ta$eta == best_cfg$eta & ta$max_depth == best_cfg$max_depth &
                    ta$min_child_weight == best_cfg$min_child_weight &
                    ta$subsample == best_cfg$subsample &
                    ta$colsample_bytree == best_cfg$colsample_bytree &
                    ta$gamma == best_cfg$gamma & ta$alpha == best_cfg$alpha &
                    ta$lambda == best_cfg$lambda)
nrounds <- round(median(ta$best_iter[sel_rows], na.rm = TRUE))
if (!is.finite(nrounds) || nrounds < 10) nrounds <- 200
logmsg("Chosen config median inner AUC: %.4f, nrounds=%d", best_cfg$best_auc, nrounds)

# ---------------- refit on all development data ----------------
set.seed(SEED)
fit_prm <- prm; fit_prm$objective <- "binary:logistic"; fit_prm$eval_metric <- "auc"
fit_prm$nthread <- N_CORES
dtrain <- xgboost::xgb.DMatrix(as.matrix(t(expr_tr)), label = y_tr)
final_model <- xgboost::xgb.train(params = fit_prm, data = dtrain, nrounds = nrounds,
                                  verbose = 0)

# ---------------- prediction object ----------------
predict_final <- function(expr_mat) {
  stopifnot(all(final$genes %in% rownames(expr_mat)))
  X <- as.matrix(t(expr_mat[final$genes, , drop = FALSE]))
  p_raw <- stats::predict(final$model, xgboost::xgb.DMatrix(X))
  p_cal <- apply_calibration(final$calibration, p_raw)
  list(raw = p_raw, calibrated = p_cal, class = ifelse(p_cal >= final$threshold, "Tumor", "Normal"))
}

final <- list(
  model = final_model,
  genes = genes,
  params = prm,
  nrounds = nrounds,
  calibration = cal$cal,
  calibration_method = cal$method,
  threshold = cal$threshold,
  predict = predict_final,
  class_levels = CLASS_LEVELS,
  created = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
saveRDS(final, file.path(DIR_MODELS, "final_model.rds"))
logmsg("Final model saved: %d genes, nrounds=%d, threshold=%.4f",
       length(genes), nrounds, final$threshold)

# ---------------- Table_15_Final_Hyperparameters.csv ----------------
hp <- data.frame(
  parameter = c("algorithm", "n_genes", "eta", "max_depth", "min_child_weight",
                "subsample", "colsample_bytree", "gamma", "alpha", "lambda",
                "nrounds", "objective", "eval_metric",
                "calibration_method", "threshold"),
  value = c("XGBoost (gbtree, binary:logistic)", length(genes),
            as.character(prm$eta), as.character(prm$max_depth),
            as.character(prm$min_child_weight), as.character(prm$subsample),
            as.character(prm$colsample_bytree), as.character(prm$gamma),
            as.character(prm$alpha), as.character(prm$lambda),
            as.character(nrounds), "binary:logistic", "auc",
            cal$method, sprintf("%.4f", final$threshold))
)
write.csv(hp, file.path(DIR_TABLES, "Table_15_Final_Hyperparameters.csv"), row.names = FALSE)
logmsg("Table_15_Final_Hyperparameters.csv written.")
