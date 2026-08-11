# ============================================================
# 07 - Nested cross-validation with leakage-safe feature selection
# Outer: 5-fold x 5 repeats (patient-grouped)
# Inner: 5-fold (patient-grouped); staged feature selection +
#        XGBoost random-search tuning inside each outer training fold.
# Produces out-of-fold predictions for calibration, feature-selection
# frequencies, tuning results, and nested-CV performance.
# ============================================================
PROJECT_DIR <- "F:/independent projects/cancer/breast cancer/Breast cancer 8"
source(file.path(PROJECT_DIR, "config/config.R"))
source(file.path(PROJECT_DIR, "R/utils.R"))
suppressMessages({
  library(dplyr); library(ggplot2); library(doParallel); library(foreach)
  library(xgboost); library(ComplexHeatmap); library(circlize)
})

logmsg <- function(...) cat(sprintf(...), "\n")

CHECKPOINT <- file.path(DIR_MODELS, "nested_cv_results.rds")
if (file.exists(CHECKPOINT)) {
  logmsg("Nested CV results already present; loading checkpoint.")
  nc <- readRDS(CHECKPOINT)
} else {

  disc <- load_discovery()
  common <- readRDS(file.path(DIR_MODELS, "common_genes.rds"))
  split  <- readRDS(file.path(DIR_MODELS, "primary_split.rds"))
  deg    <- readRDS(file.path(DIR_MODELS, "deg_results.rds"))

  expr  <- disc$expr[common, , drop = FALSE]
  meta  <- disc$meta
  train_samples <- split$train_samples

  expr_tr <- expr[, train_samples, drop = FALSE]
  meta_tr <- meta[match(train_samples, meta$sample_id), ]
  y_tr    <- as.numeric(meta_tr$class == "Tumor")

  logmsg("Development set: %d samples, %d genes", ncol(expr_tr), nrow(expr_tr))
  logmsg("Class balance in development: Tumor=%d Normal=%d", sum(y_tr == 1), sum(y_tr == 0))

  # ---------------- outer folds ----------------
  set.seed(SEED)
  outer_folds <- list()
  for (r in seq_len(OUTER_REPEATS)) {
    f <- grouped_folds(ncol(expr_tr), meta_tr$patient_id, OUTER_FOLDS,
                       seed = SEED + r * 1000)
    for (k in seq_len(OUTER_FOLDS)) {
      outer_folds[[length(outer_folds) + 1]] <- list(repeat_id = r, fold = k, idx = f[[k]])
    }
  }
  n_outer <- length(outer_folds)
  logmsg("Total outer folds: %d", n_outer)

  # ---------------- random-search configs ----------------
  set.seed(SEED)
  cfg_grid <- expand.grid(
    eta = c(0.01, 0.05, 0.1),
    max_depth = c(3, 4, 6),
    min_child_weight = c(1, 5),
    subsample = c(0.7, 0.9),
    colsample_bytree = c(0.6, 0.8),
    gamma = c(0, 0.1),
    alpha = c(0, 0.5),
    lambda = c(0.5, 1)
  )
  cfg_draws <- cfg_grid[sample(nrow(cfg_grid), TUNE_RANDOM_SEARCH_N), , drop = FALSE]
  configs <- lapply(seq_len(nrow(cfg_draws)), function(i) as.list(cfg_draws[i, ]))

  # ---------------- parallel outer loop ----------------
  cl <- parallel::makeCluster(N_CORES)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))

  t0 <- Sys.time()
  results <- foreach(o = seq_len(n_outer), .packages = c("xgboost", "limma", "caret"),
                     .export = c("grouped_folds", "select_stable_features", "limma_rank_genes",
                                 "near_zero_variance", "correlation_filter", "xgb_cv_configs",
                                 "fit_xgboost", "predict_xgb",
                                 "STABILITY_FREQ", "INNER_FOLDS", "FOLD_TOP_DEG", "CORR_CUTOFF",
                                 "TUNE_RANDOM_SEARCH_N", "TUNE_NROUNDS_MAX", "TUNE_EARLY_STOP",
                                 "SEED"),
                     .combine = list, .multicombine = TRUE,
                     .maxcombine = n_outer) %dopar% {
    of <- outer_folds[[o]]
    test_idx  <- of$idx
    train_idx <- setdiff(seq_len(ncol(expr_tr)), test_idx)

    expr_out_tr <- expr_tr[, train_idx, drop = FALSE]
    expr_out_te <- expr_tr[, test_idx, drop = FALSE]
    y_out_tr <- y_tr[train_idx]
    pat_tr   <- meta_tr$patient_id[train_idx]

    # inner folds
    inner_folds <- grouped_folds(length(train_idx), pat_tr, INNER_FOLDS,
                                 seed = SEED + o)

    # ---- staged feature selection across inner training sets ----
    feat_freq <- table(c())
    inner_feature_sets <- lapply(seq_len(INNER_FOLDS), function(inn) {
      itrain_idx <- train_idx[inner_folds[[inn]]]
      in_tr <- setdiff(train_idx, itrain_idx)
      feats <- tryCatch(
        select_stable_features(expr_tr[, in_tr, drop = FALSE],
                               y_tr[in_tr],
                               meta_tr$patient_id[in_tr],
                               top_deg = FOLD_TOP_DEG, corr_cutoff = CORR_CUTOFF),
        error = function(e) character(0)
      )
      feats
    })
    feat_freq <- sort(table(unlist(inner_feature_sets)), decreasing = TRUE)
    stable <- names(feat_freq)[feat_freq / INNER_FOLDS >= STABILITY_FREQ]
    if (length(stable) < 5) stable <- names(head(sort(feat_freq, decreasing = TRUE), 20))
    stable <- intersect(stable, rownames(expr_tr))

    # ---- tuning with inner CV ----
    inner_y <- y_tr[train_idx]
    tune_res <- tryCatch(
      xgb_cv_configs(expr_tr[stable, train_idx, drop = FALSE], inner_y,
                     configs, folds = inner_folds),
      error = function(e) NULL
    )
    if (is.null(tune_res) || nrow(tune_res) == 0) {
      best_prm <- list(eta = 0.05, max_depth = 4, min_child_weight = 2,
                       subsample = 0.8, colsample_bytree = 0.7,
                       gamma = 0, alpha = 0, lambda = 1)
      best_nr <- 300
      best_auc <- NA
    } else {
      best <- tune_res[which.max(tune_res$best_auc), ]
      best_prm <- list(eta = best$eta, max_depth = best$max_depth,
                       min_child_weight = best$min_child_weight,
                       subsample = best$subsample,
                       colsample_bytree = best$colsample_bytree,
                       gamma = best$gamma, alpha = best$alpha, lambda = best$lambda)
      best_nr <- round(best$best_iter)
      if (!is.finite(best_nr) || best_nr < 10) best_nr <- 200
      best_auc <- best$best_auc
    }

    # ---- fit and predict outer test ----
    set.seed(SEED + o)
    dtrain <- xgboost::xgb.DMatrix(as.matrix(t(expr_tr[stable, train_idx, drop = FALSE])),
                                   label = y_tr[train_idx])
    prm <- best_prm; prm$objective <- "binary:logistic"; prm$eval_metric <- "auc"; prm$nthread <- 1
    m <- xgboost::xgb.train(params = prm, data = dtrain, nrounds = best_nr, verbose = 0)
    dtest <- xgboost::xgb.DMatrix(as.matrix(t(expr_tr[stable, test_idx, drop = FALSE])))
    pred <- stats::predict(m, dtest)

    imp <- xgboost::xgb.importance(feature_names = stable, model = m)

    list(
      repeat_id = of$repeat_id, fold = of$fold,
      test_idx = test_idx, pred = pred, truth = y_tr[test_idx],
      stable_features = stable, feat_freq = feat_freq,
      tune_res = tune_res, best_auc = best_auc, best_nr = best_nr,
      importance = imp
    )
  }
  parallel::stopCluster(cl)
  on.exit(NULL)

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  logmsg("Nested CV completed in %.1f minutes", elapsed)

  # ---------------- aggregate ----------------
  oof <- do.call(rbind, lapply(results, function(r) {
    data.frame(sample_id = colnames(expr_tr)[r$test_idx],
               truth = ifelse(r$truth == 1, "Tumor", "Normal"),
               oof_prob = r$pred,
               repeat_id = r$repeat_id, outer_fold = r$fold,
               stringsAsFactors = FALSE)
  }))
  oof$truth <- factor(oof$truth, levels = CLASS_LEVELS)
  oof <- oof[!duplicated(oof$sample_id), ]
  logmsg("OOF predictions: %d unique development samples", nrow(oof))

  # feature selection frequency across outer folds
  all_freq <- do.call(rbind, lapply(results, function(r) {
    if (length(r$feat_freq) == 0) return(NULL)
    data.frame(gene = names(r$feat_freq), freq = as.numeric(r$feat_freq),
               outer_fold = r$fold, repeat_id = r$repeat_id)
  }))
  sel_freq <- aggregate(freq ~ gene, all_freq, function(x) sum(x > 0) / n_outer)
  sel_freq <- sel_freq[order(-sel_freq$freq), ]
  final_stable <- sel_freq$gene[sel_freq$freq >= FINAL_STABILITY_FREQ]
  if (length(final_stable) < 5) final_stable <- head(sel_freq$gene, 20)

  # outer CV performance
  perf <- do.call(rbind, lapply(results, function(r) {
    m <- compute_metrics(r$pred, ifelse(r$truth == 1, "Tumor", "Normal"))
    data.frame(repeat_id = r$repeat_id, outer_fold = r$fold,
               n_features = length(r$stable_features),
               metric = m$metric, value = m$value)
  }))

  # tuning results combined
  tune_all <- do.call(rbind, lapply(results, function(r) {
    if (is.null(r$tune_res)) return(NULL)
    r$tune_res$repeat_id <- r$repeat_id
    r$tune_res$outer_fold <- r$fold
    r$tune_res
  }))

  nc <- list(
    oof = oof, sel_freq = sel_freq, final_stable = final_stable,
    perf = perf, tune_all = tune_all, outer_folds = outer_folds,
    results = results, elapsed = elapsed,
    n_features_per_fold = sapply(results, function(r) length(r$stable_features))
  )
  saveRDS(nc, CHECKPOINT)
  logmsg("Saved nested CV checkpoint.")
}

# ensure objects needed downstream are available (checkpoint path)
if (!exists("expr_tr", inherits = FALSE)) {
  disc <- load_discovery()
  common <- readRDS(file.path(DIR_MODELS, "common_genes.rds"))
  split  <- readRDS(file.path(DIR_MODELS, "primary_split.rds"))
  expr_tr <- disc$expr[common, split$train_samples, drop = FALSE]
}

# ---------------- tables ----------------
logmsg("Final stable features: %d", length(nc$final_stable))
write.csv(data.frame(gene = nc$sel_freq$gene, selection_frequency = nc$sel_freq$freq),
          file.path(DIR_TABLES, "Table_12_Feature_Selection_Frequency.csv"), row.names = FALSE)
write.csv(data.frame(gene = nc$final_stable, stringsAsFactors = FALSE),
          file.path(DIR_TABLES, "Table_13_Final_Stable_Features.csv"), row.names = FALSE)

# nested CV indices
idx_df <- do.call(rbind, lapply(seq_along(nc$outer_folds), function(i) {
  of <- nc$outer_folds[[i]]
  data.frame(repeat_id = of$repeat_id, outer_fold = of$fold,
             sample_id = colnames(expr_tr)[of$idx], role = "outer_test",
             stringsAsFactors = FALSE)
}))
write.csv(idx_df, file.path(DIR_TABLES, "Table_11_Nested_CV_Indices.csv"), row.names = FALSE)

# tuning results
if (!is.null(nc$tune_all)) {
  write.csv(nc$tune_all, file.path(DIR_TABLES, "Table_14_XGBoost_Tuning_Results.csv"), row.names = FALSE)
}

# outer CV performance summary
perf_sum <- nc$perf
write.csv(perf_sum, file.path(DIR_TABLES, "Table_16_Outer_CV_Performance.csv"), row.names = FALSE)

# ---------------- figures ----------------
# feature selection frequency
freq_df <- data.frame(gene = nc$sel_freq$gene, freq = nc$sel_freq$freq)
freq_df <- freq_df[1:min(50, nrow(freq_df)), ]
freq_df$gene <- factor(freq_df$gene, levels = rev(freq_df$gene))
p_freq <- ggplot(freq_df, aes(x = gene, y = freq)) +
  geom_col(fill = "#0072B2", colour = "black", width = 0.75) +
  coord_flip() + theme_bc(base_size = 8) +
  labs(x = "", y = "Selection frequency across outer folds",
       title = "Top genes by selection frequency (stability selection)")
save_tiff_600dpi(p_freq, "Figure_09_Feature_Selection_Frequency.tiff", width = 8, height = 9)

# feature-rank stability: rank of each gene across folds (by frequency within fold)
all_freq2 <- do.call(rbind, lapply(nc$results, function(r) {
  if (length(r$feat_freq) == 0) return(NULL)
  data.frame(gene = names(r$feat_freq), freq = as.numeric(r$feat_freq),
             fold = paste0(r$repeat_id, ".", r$fold))
}))
top_genes_rank <- names(sort(table(all_freq2$gene), decreasing = TRUE))[1:40]
rank_df <- all_freq2[all_freq2$gene %in% top_genes_rank, ]
rank_df$gene <- factor(rank_df$gene, levels = rev(top_genes_rank))
p_rank <- ggplot(rank_df, aes(x = fold, y = gene, fill = freq)) +
  geom_tile() +
  scale_fill_gradient(low = "#F7F7F7", high = "#2166AC", name = "inner-fold\nselection\nfrequency") +
  theme_bc(base_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        legend.position = "bottom") +
  labs(x = "Outer fold (repeat.fold)", y = "",
       title = "Feature rank stability across outer folds")
save_tiff_600dpi(p_rank, "Figure_10_Feature_Rank_Stability.tiff", width = 10, height = 7)

# co-selection heatmap (top genes by co-occurrence)
top_cosel <- names(head(sort(table(all_freq2$gene), decreasing = TRUE), 50))
folds_genes <- lapply(unique(all_freq2$fold), function(f) all_freq2$gene[all_freq2$fold == f])
names(folds_genes) <- unique(all_freq2$fold)
co <- sapply(top_cosel, function(g) sapply(top_cosel, function(h) {
  mean(sapply(folds_genes, function(f) (g %in% f) & (h %in% f)))
}))
open_tiff_base("Figure_10_Feature_CoSelection_Heatmap.tiff", width = 8, height = 7)
ComplexHeatmap::Heatmap(co, name = "co-selection\nprobability",
                        col = colorRamp2(c(0, 0.5, 1), c("#2166AC", "#F7F7F7", "#B2182B")),
                        show_row_names = TRUE, show_column_names = TRUE,
                        row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5),
                        heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8)))
dev.off()

# final stable feature expression heatmap (development set)
disc <- load_discovery()
split <- readRDS(file.path(DIR_MODELS, "primary_split.rds"))
zmat <- t(scale(t(disc$expr[nc$final_stable, split$train_samples, drop = FALSE])))
zmat[!is.finite(zmat)] <- 0
meta_tr <- disc$meta[match(split$train_samples, disc$meta$sample_id), ]
ann_col <- data.frame(Class = meta_tr$class, row.names = colnames(zmat))
open_tiff_base("Figure_10_Final_Stable_Features_Heatmap.tiff", width = 8, height = 8)
ComplexHeatmap::Heatmap(
  zmat, name = "z-score",
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B")),
  top_annotation = HeatmapAnnotation(df = ann_col,
                                     col = list(Class = c(Tumor = COL_TUMOR, Normal = COL_NORMAL))),
  show_column_names = FALSE, row_names_gp = gpar(fontsize = 5),
  column_split = factor(meta_tr$class, levels = CLASS_LEVELS),
  use_raster = TRUE, raster_quality = 5,
  heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8)))
dev.off()

logmsg("Saved nested CV tables and figures.")
