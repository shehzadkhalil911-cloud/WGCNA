# 14_RF_feature_selection.R
source("config/analysis_config.R")
source("R/load_functions.R")
ensure_paths()
suppressMessages({ library(ranger); library(pROC); library(ggplot2); library(dplyr); library(caret) })

combat <- readRDS(file.path(paths$discovery, "discovery_combat.rds"))
meta <- readRDS(file.path(paths$discovery, "discovery_cohort_raw.rds"))$metadata
crg_sum <- readRDS(file.path(paths$res_cross, "DE_CRG_summary.rds"))
pool <- crg_sum$candidate_pool
x <- t(combat[pool, ])
y <- factor(meta$class, levels = c("Normal", "Tumor"))

split <- readRDS(file.path(paths$res_ml, "LASSO", "train_test_split.rds"))
train_idx <- split$train; test_idx <- split$test

set.seed(GLOBAL_SEED)
tr_df <- data.frame(x[train_idx, , drop = FALSE], class = y[train_idx], stringsAsFactors = FALSE)
rf <- ranger(class ~ ., data = tr_df, num.trees = 1000, importance = "permutation",
             probability = TRUE, seed = GLOBAL_SEED, mtry = floor(sqrt(ncol(x[train_idx,]))),
             min.node.size = 5)
imp <- importance(rf)
ord <- order(imp, decreasing = TRUE)
top <- colnames(x)[ord[1:min(TOP_RF_GENES, length(ord))]]
pred_test <- predict(rf, x[test_idx, ])$predictions[, "Tumor"]
roc_obj <- roc(y[test_idx], pred_test)
auc_val <- auc(roc_obj)

res <- list(rf = rf, importance = imp, top_genes = top,
            pred_test = pred_test, y_test = y[test_idx], auc_test = auc_val)
saveRDS(res, file.path(paths$res_ml, "Random_Forest", "rf_result.rds"))
write.csv(data.frame(gene = names(imp), importance = imp)[order(imp, decreasing = TRUE), ],
          file.path(paths$res_ml, "Random_Forest", "variable_importance.csv"), row.names = FALSE)
write.csv(data.frame(gene = top, rank = seq_along(top)),
          file.path(paths$res_ml, "Random_Forest", "top_genes.csv"), row.names = FALSE)

## plots
imp_df <- data.frame(gene = names(imp), importance = imp) %>% dplyr::top_n(TOP_RF_GENES, importance)
p_imp <- ggplot(imp_df, aes(x = reorder(gene, importance), y = importance)) +
  geom_bar(stat = "identity", fill = "#B07AA1") + coord_flip() + labs(title = "Random Forest variable importance") + pub_theme()
save_publication_plot(p_imp, file.path(paths$res_ml, "Random_Forest", "variable_importance"), width = 7, height = 7)
cat("RF complete. Top genes:", paste(head(top, 10), collapse = ", "), " Test AUC:", round(auc_val,3), "\n")
quit(save = "no", status = 0, runLast = FALSE)
