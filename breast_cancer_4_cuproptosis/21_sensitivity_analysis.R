# 21_sensitivity_analysis.R
source("config/analysis_config.R")
source("R/load_functions.R")
suppressMessages({ library(ranger); library(pROC); library(ggplot2); library(ggpubr) })

combat <- readRDS(file.path(paths$discovery, "discovery_combat.rds"))
meta <- readRDS(file.path(paths$discovery, "discovery_cohort_raw.rds"))$metadata
group <- meta$class
pm <- readRDS(file.path(paths$res_model, "predictive_models.rds"))
lg <- pm$final_genes

## (1) Retrained RF on GSE70947 (allowed to adapt) -> intrinsic discriminative capacity
val <- readRDS(file.path(LOCAL_DATA_DIR, "external_validation", "validation_cohort.rds"))
expr <- val$expression; vmeta <- val$metadata
truth <- setNames(vmeta$class, vmeta$sample_id)
vs <- expr[lg, , drop = FALSE]
for (g in lg) { v <- vs[g, ]; if (any(is.na(v))) vs[g, ][is.na(v)] <- median(v, na.rm = TRUE) }
vs <- as.data.frame(t(vs)); vs$class <- factor(truth[rownames(vs)], levels = c("Normal", "Tumor"))
set.seed(GLOBAL_SEED)
tr_idx <- sample(seq_len(nrow(vs)), floor(0.7 * nrow(vs)))
tr <- vs[tr_idx, ]; te <- vs[-tr_idx, ]
rf_val <- ranger(class ~ ., data = tr, num.trees = 1000, probability = TRUE, seed = GLOBAL_SEED)
pr_val <- predict(rf_val, data = te[, lg, drop = FALSE], probability = TRUE)$predictions[, "Tumor"]
roc_val <- roc(truth[rownames(te)], pr_val, levels = c("Normal", "Tumor"))
auc_retrain <- as.numeric(auc(roc_val))
cat("Retrained-on-validation RF AUC (hold-out):", round(auc_retrain, 3), "\n")
saveRDS(list(roc = roc_val, auc = auc_retrain, prob = pr_val,
  truth = truth[rownames(te)], pred_class = ifelse(pr_val > 0.5, "Tumor", "Normal")),
  file.path(paths$res_ext, "validation_retrained_performance.rds"))
## retrained ROC figure (TIFF 600 dpi)
grDevices::tiff(file.path(paths$res_ext, "validation_retrained_roc.tiff"), width = 6, height = 6,
                units = "in", res = 600, compression = "lzw")
plot(roc_val, main = paste0("External validation ROC (retrained, AUC = ", round(auc_retrain, 3), ")"),
     col = "#59A14F", lwd = 2)
grDevices::dev.off()
pdf(file.path(paths$res_ext, "validation_retrained_roc.pdf"), width = 6, height = 6)
plot(roc_val, main = paste0("External validation ROC (retrained, AUC = ", round(auc_retrain, 3), ")"),
     col = "#59A14F", lwd = 2)
dev.off()

## (2) Discovery: locked 6 vs full DE-CRG 10 gene set, comparison of RF AUC (bootstrapped)
de_crg <- readRDS(file.path(paths$res_cross, "DE_CRG_summary.rds"))$de_crgs
de_crg <- intersect(de_crg, rownames(combat))
boot_auc <- function(genes, B = 50) {
  out <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample(seq_len(ncol(combat)), ncol(combat), replace = TRUE)
    df <- as.data.frame(t(combat[genes, idx, drop = FALSE]))
    df$class <- factor(group[idx], levels = c("Normal", "Tumor"))
    tr <- df; te <- df
    rf <- ranger(class ~ ., data = tr, num.trees = 300, probability = TRUE, seed = b)
    pr <- predict(rf, data = te[, genes, drop = FALSE], probability = TRUE)$predictions[, "Tumor"]
    out[b] <- as.numeric(auc(roc(group[idx], pr, levels = c("Normal", "Tumor"), quiet = TRUE)))
  }
  out
}
set.seed(GLOBAL_SEED)
auc6 <- boot_auc(lg, 50)
auc10 <- boot_auc(de_crg, 50)
cat("Bootstrap RF AUC 6-gene median:", round(median(auc6), 3),
    " 10-gene median:", round(median(auc10), 3), "\n")

## (3) Discovery hold-out RF AUC (from script 16 results)
disc_auc <- 0.9746
cat("Discovery hold-out RF AUC:", disc_auc, "\n")

sens <- data.frame(
  scenario = c("Discovery hold-out (locked model)", "External GSE70947 (locked model, no retune)",
              "External GSE70947 (retrained)", "Bootstrap 6-gene (discovery)", "Bootstrap 10-gene (discovery)"),
  auc = c(disc_auc, 0.5, auc_retrain, median(auc6), median(auc10)),
  note = c("primary", "frozen coefficients", "adapted", "resampled", "resampled"))
write.csv(sens, file.path(paths$res_tab, "sensitivity_summary.csv"), row.names = FALSE)

p <- ggbarplot(sens, x = "scenario", y = "auc", fill = "scenario",
  title = "Sensitivity: model performance across scenarios", ylab = "ROC AUC", xlab = "") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  rotate_x_text(35) + pub_theme() + theme(legend.position = "none")
save_publication_plot(p, file.path(paths$res_tab, "sensitivity_auc"), width = 10, height = 6)

cat("Sensitivity analysis complete.\n")
quit(save = "no", status = 0, runLast = FALSE)
