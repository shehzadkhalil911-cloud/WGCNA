# 20_external_validation.R
source("config/analysis_config.R")
source("R/load_functions.R")
suppressMessages({ library(ranger); library(ggplot2); library(pROC); library(kernelshap); library(shapviz) })

val <- readRDS(file.path(LOCAL_DATA_DIR, "external_validation", "validation_cohort.rds"))
expr <- val$expression
meta <- val$metadata
truth <- setNames(meta$class, meta$sample_id)
colnames(expr) <- meta$sample_id[match(colnames(expr), meta$sample_id)]
group <- truth

combat <- readRDS(file.path(paths$discovery, "discovery_combat.rds"))
pm <- readRDS(file.path(paths$res_model, "predictive_models.rds"))
rf <- pm$fits$rf
lock_genes <- pm$final_genes
cat("Locked genes:", paste(lock_genes, collapse = ", "), "\n")

## map locked genes into validation (symbols)
avail <- intersect(lock_genes, rownames(expr))
if (length(avail) < length(lock_genes)) stop("Missing locked genes in validation: ",
  paste(setdiff(lock_genes, rownames(expr)), collapse = ", "))
val_sub <- expr[lock_genes, , drop = FALSE]

## NA handling: median impute per gene (robust to negatives)
na_count <- rowSums(is.na(val_sub))
cat("NA counts per locked gene:", paste(lock_genes, na_count, sep = "=", collapse = ", "), "\n")
for (g in lock_genes) {
  v <- val_sub[g, ]
  if (any(is.na(v))) val_sub[g, ][is.na(v)] <- median(v, na.rm = TRUE)
}

## standardize to discovery distribution (z-score using discovery mean/sd)
disc_ms <- list(
  mean = rowMeans(combat[lock_genes, , drop = FALSE], na.rm = TRUE),
  sd = apply(combat[lock_genes, , drop = FALSE], 1, sd, na.rm = TRUE))
val_z <- (val_sub - disc_ms$mean) / disc_ms$sd
val_z[!is.finite(val_z)] <- 0
val_df <- as.data.frame(t(val_z))
val_df$class <- factor(truth[rownames(val_df)], levels = c("Normal", "Tumor"))

## predict with locked model (no retuning)
pred_rf <- predict(rf, data = val_df[, lock_genes, drop = FALSE], probability = TRUE)
prob_tumor <- pred_rf$predictions[, "Tumor"]
roc_obj <- roc(truth[rownames(val_df)], prob_tumor, levels = c("Normal", "Tumor"), direction = "<")
auc_val <- as.numeric(auc(roc_obj))
cat("Validation RF AUC:", round(auc_val, 4), "\n")
saveRDS(list(roc = roc_obj, auc = auc_val, prob = prob_tumor,
  truth = truth[rownames(val_df)], pred_class = ifelse(prob_tumor > 0.5, "Tumor", "Normal")),
  file.path(paths$res_ext, "validation_performance.rds"))

pdf(file.path(paths$res_ext, "validation_roc.pdf"), width = 6, height = 6)
plot(roc_obj, main = paste0("External validation ROC (AUC = ", round(auc_val, 3), ")"),
  col = "#E15759", lwd = 2)
dev.off()
grDevices::tiff(file.path(paths$res_ext, "validation_roc.tiff"), width = 6, height = 6,
                units = "in", res = 600, compression = "lzw")
plot(roc_obj, main = paste0("External validation ROC (AUC = ", round(auc_val, 3), ")"),
  col = "#E15759", lwd = 2)
grDevices::dev.off()

## SHAP on validation predictions
set.seed(GLOBAL_SEED)
bg <- val_df[sample(seq_len(nrow(val_df)), min(50, nrow(val_df))), lock_genes, drop = FALSE]
ks <- kernelshap(rf, X = as.data.frame(val_df[, lock_genes, drop = FALSE]),
  bg_X = as.data.frame(bg), species = "Tumor", parallel = FALSE)
sv <- shapviz(ks)
shap_imp <- sv_importance(sv, kind = "bee")
ggsave(file.path(paths$res_ext, "validation_shap_importance.pdf"), shap_imp, width = 8, height = 5)
save_publication_plot(shap_imp, file.path(paths$res_ext, "validation_shap_importance"), width = 8, height = 5)
pdf(file.path(paths$res_ext, "validation_shap_dependence.pdf"), width = 7, height = 5)
for (g in lock_genes) print(sv_dependence(sv, g))
dev.off()
grDevices::tiff(file.path(paths$res_ext, "validation_shap_dependence.tiff"), width = 7, height = 5,
                units = "in", res = 600, compression = "lzw")
for (g in lock_genes) print(sv_dependence(sv, g))
grDevices::dev.off()
saveRDS(sv, file.path(paths$res_ext, "validation_shap.rds"))

## immune infiltration on validation (xCell primary)
cat("Running xCell on validation...\n"); flush.console()
expr_imp <- expr
for (i in seq_len(nrow(expr_imp))) {
  r <- expr_imp[i, ]
  if (any(is.na(r))) expr_imp[i, ][is.na(r)] <- median(r, na.rm = TRUE)
}
val_immune <- run_xcell(expr_imp)
saveRDS(val_immune, file.path(paths$res_ext, "validation_xcell.rds"))
val_imm_df <- as.data.frame(t(val_immune))
write.csv(val_imm_df, file.path(paths$res_ext, "validation_immune_scores.csv"))

## correlation of locked genes with immune cells in validation
imm_cols <- colnames(val_imm_df)
ct_rows <- list()
for (g in lock_genes) {
  gv <- as.numeric(val_z[g, colnames(val_z)])
  for (ic in imm_cols) {
    iv <- as.numeric(val_imm_df[[ic]])
    if (any(!is.finite(gv)) || any(!is.finite(iv))) next
    rs <- cor.test(gv, iv, method = "spearman")
    ct_rows[[paste(g, ic, sep = "|")]] <- data.frame(gene = g, immune_feature = ic,
      rho = rs$estimate, p_value = rs$p.value)
  }
}
ct <- do.call(rbind, ct_rows)
ct$adj_p <- p.adjust(ct$p_value, method = "BH")
write.csv(ct, file.path(paths$res_ext, "validation_gene_immune_crosstalk.csv"), row.names = FALSE)
cat("Validation crosstalk significant (adj.P<0.05):", sum(ct$adj_p < 0.05), "\n")

cat("External validation complete.\n")
quit(save = "no", status = 0, runLast = FALSE)
