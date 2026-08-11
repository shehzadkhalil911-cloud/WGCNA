# 03_preprocess_individual_datasets.R
# The local preprocessed RDS files already contain normalized gene-level matrices.
# We load them, record scale assessment + annotation coverage, and re-derive QC
# figures WITHOUT overwriting the original Datasets/ files.
source("config/analysis_config.R")
source("R/load_functions.R")
ensure_paths()
suppressMessages({ library(dplyr); library(openxlsx) })

LOCAL_PROC <- file.path(LOCAL_DATA_DIR, "processed")
proc_paths <- setNames(file.path(LOCAL_PROC, paste0(c(DISCOVERY_GSE, VALIDATION_GSE), "_processed.rds")),
                       c(DISCOVERY_GSE, VALIDATION_GSE))

proc_list <- list(); scale_df <- list(); ann_df <- list()
for (g in c(DISCOVERY_GSE, VALIDATION_GSE)) {
  obj <- readRDS(proc_paths[[g]])
  proc_list[[g]] <- obj
  scale_df[[g]] <- cbind(dataset = g, as.data.frame(as.list(assess_expression_scale(obj$expression_raw))))
  ann_df[[g]] <- annotation_summary(obj)
  ## QC figures (saved under results/preprocessing, not in Datasets/)
  make_qc_plots(obj, g, file.path(paths$res_qc, g))
  ## save cleaned gene-level matrix + phenotype copy
  gene_mat <- obj$expression
  colnames(gene_mat) <- obj$phenotype$geo_accession[match(colnames(gene_mat), rownames(obj$phenotype))]
  saveRDS(gene_mat, file.path(paths$processed_ind, paste0(g, "_gene_expression.rds")))
  write.csv(as.data.frame(gene_mat), file.path(paths$processed_ind, paste0(g, "_gene_expression.csv")))
}
scale_all <- do.call(rbind, scale_df)
ann_all <- do.call(rbind, ann_df)
write.csv(scale_all, file.path(paths$res_preproc, "expression_scale_assessment.csv"), row.names = FALSE)
write.csv(ann_all, file.path(paths$res_preproc, "annotation_coverage.csv"), row.names = FALSE)
saveRDS(proc_list, file.path(paths$processed_ind, "processed_objects.rds"))
cat("Per-dataset preprocessing/QC complete.\n")
print(ann_all)
print(scale_all)
quit(save = "no", status = 0, runLast = FALSE)
