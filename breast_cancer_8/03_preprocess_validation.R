# ============================================================
# 03 - Independent preprocessing: all three validation cohorts
# GSE15852 (GPL96, 86 samples, 43 paired, MAS5)
# GSE29431 (GPL570, 66 samples, 12N/54T)
# GSE42568 (GPL570, 121 samples, 17N/104T)
# Same pipeline for each: probe filter, symbol validation,
# log2-if-needed, quantile normalization, max-IQR collapse.
# ============================================================
PROJECT_DIR <- "F:/independent projects/cancer/breast cancer/Breast cancer 8"
source(file.path(PROJECT_DIR, "config/config.R"))
source(file.path(PROJECT_DIR, "R/utils.R"))
source(file.path(PROJECT_DIR, "R/preprocess.R"))
suppressMessages({
  library(Biobase); library(ggplot2); library(dplyr); library(tidyr)
  library(ComplexHeatmap); library(org.Hs.eg.db); library(AnnotationDbi)
  library(preprocessCore)
})

logmsg <- function(...) cat(sprintf(...), "\n")
logmsg("=== 03 Preprocessing validation cohorts ===")

preprocess_validation <- function(gse, eset_path, soft_path, meta_csv = NULL,
                                  platform, title2class = NULL, fig_prefix) {
  logmsg("===== %s (%s) =====", gse, platform)
  es <- readRDS(eset_path)
  pd <- pData(es)
  exprs_raw <- exprs(es)

  ann <- parse_soft_platform(soft_path)
  stopifnot(nrow(ann) == nrow(exprs_raw))
  rowmatch <- match(rownames(exprs_raw), as.character(ann$ID))
  stopifnot(!any(is.na(rowmatch)))
  ann <- ann[rowmatch, ]
  stopifnot(identical(as.character(ann$ID), rownames(exprs_raw)))

  # ---- class / patient assignment ----
  if (!is.null(title2class)) {
    cls <- ifelse(grepl("^Normal", pd$title), "Normal", ifelse(grepl("^Cancer", pd$title), "Tumor", NA))
    stopifnot(all(!is.na(cls)))
    patient <- sub("^(Normal|Cancer) (BC[0-9]+)[NT]$", "\\2", pd$title)
  } else {
    cc <- read.csv(meta_csv, stringsAsFactors = FALSE)
    cc <- cc[match(sampleNames(es), cc$sample_id), ]
    stopifnot(all(sampleNames(es) == cc$sample_id))
    cls <- cc$class
    patient <- sampleNames(es)
  }
  groups <- cls
  names(groups) <- sampleNames(es)
  names(patient) <- sampleNames(es)

  meta <- data.frame(
    sample_id = sampleNames(es),
    dataset = gse,
    class = groups,
    patient_id = patient,
    title = pd$title,
    stringsAsFactors = FALSE
  )
  rownames(meta) <- meta$sample_id

  # ---- probe filtering ----
  gsym_raw <- first_symbol(ann$`Gene Symbol`)
  is_bad <- grepl("^AFFX", rownames(exprs_raw)) | grepl("^AFFX", ann$ID) | is.na(gsym_raw)
  expr_qc <- exprs_raw[!is_bad, , drop = FALSE]
  gsym <- gsym_raw[!is_bad]
  valid <- validate_symbols(gsym)
  expr_qc <- expr_qc[valid, , drop = FALSE]
  gsym <- gsym[valid]
  logmsg("  Control/no-symbol removed: %d; remaining probes: %d; official symbols: %d",
         sum(is_bad), nrow(expr_qc), sum(valid))

  # ---- log2 if needed + quantile normalization ----
  logmsg("  Raw scale: min %.2f max %.2f median %.2f", min(expr_qc), max(expr_qc), median(expr_qc))
  if (median(expr_qc, na.rm = TRUE) > 100) {
    logmsg("  %s: linear scale detected, applying log2", gse)
    expr_qc <- log2(expr_qc)
    needs_log2 <- TRUE
  } else needs_log2 <- FALSE
  expr_norm <- quantile_normalize(expr_qc)
  expr_gene <- collapse_probes_to_genes(expr_norm, gsym, method = COLLAPSE_METHOD)
  logmsg("  Unique genes after collapse: %d", nrow(expr_gene))

  # ---- outliers ----
  out_df <- assess_outliers(expr_gene)
  logmsg("  Samples flagged as outliers: %d", sum(out_df$is_outlier))
  meta$is_outlier <- out_df$is_outlier[match(meta$sample_id, out_df$sample_id)]

  # ---- preprocessing summary ----
  tab <- data.frame(
    dataset = gse, platform = platform,
    n_probes_raw = nrow(exprs_raw),
    n_control_no_symbol_removed = sum(is_bad),
    n_probes_no_official_symbol = sum(valid == FALSE),
    n_probes_after_qc = nrow(expr_qc),
    n_imputed_cells = 0,
    n_genes_after_collapse = nrow(expr_gene),
    n_normal = sum(meta$class == "Normal"), n_tumor = sum(meta$class == "Tumor"),
    scale = "log2",
    normalization = "log2 + quantile (preprocessCore)",
    needs_log2 = needs_log2,
    stringsAsFactors = FALSE
  )
  write.csv(tab, file.path(DIR_TABLES, sprintf("Table_05b_Preprocessing_Summary_%s.csv", gse)),
            row.names = FALSE)

  # ---- QC figures ----
  step <- function(nm, expr) {
    tryCatch({ expr; logmsg("  figure step OK: %s", nm) },
             error = function(e) { logmsg("  figure step FAILED: %s -- %s", nm, conditionMessage(e)); NULL })
  }
  qc <- qc_box_density(expr_qc[1:min(5000, nrow(expr_qc)), ], expr_gene, groups, sprintf("%s (%s)", gse, platform))
  step("box", save_tiff_600dpi(qc$box, sprintf("%s_Normalization_QC.tiff", fig_prefix), width = 12, height = 6))
  step("dens", save_tiff_600dpi(qc$dens, sprintf("%s_Density_QC.tiff", fig_prefix), width = 10, height = 5))
  step("pca", {
    p_pca <- pca_plot(expr_gene, groups, gse)
    save_tiff_600dpi(p_pca, sprintf("%s_PCA.tiff", fig_prefix), width = 7, height = 6)
  })
  step("corr", correlation_heatmap(expr_gene, groups, sprintf("%s_Correlation_Heatmap.tiff", fig_prefix), gse))
  step("dendro", dendrogram_plot(expr_gene, groups, sprintf("%s_Dendrogram.tiff", fig_prefix), gse))

  counts <- c(probes_raw = nrow(exprs_raw), control_no_symbol = sum(is_bad),
              no_official_symbol = sum(valid == FALSE), probes_kept = nrow(expr_qc),
              genes_final = nrow(expr_gene))
  step("retention", retention_plot(counts, sprintf("%s_Probe_Retention.tiff", fig_prefix), gse))
  step("outlier", outlier_plot(out_df, sprintf("%s_Outliers.tiff", fig_prefix), gse))

  list(
    expr = expr_gene,
    expr_log2 = expr_qc,
    probe_expr = expr_norm,
    meta = meta,
    annotation = ann[!is_bad, ][valid, ],
    probe2gene = gsym,
    scale_info = list(min = min(expr_gene), max = max(expr_gene),
                      needs_log2 = FALSE, median = median(expr_gene)),
    preprocessing = tab,
    outlier_assessment = out_df,
    dataset = gse
  )
}

# ---------------- GSE15852 (paired validation) ----------------
val <- preprocess_validation("GSE15852", RAW_GSE15852_ESET, GPL96_SOFT,
                             platform = "GPL96", title2class = TRUE,
                             fig_prefix = "Figure_05_Validation")
saveRDS(val, file.path(DIR_VAL_PROC, "validation_processed.rds"))
write.csv(val$meta, file.path(DIR_VAL_PROC, "validation_metadata.csv"), row.names = FALSE)
logmsg("Saved GSE15852 processed object: %d genes x %d samples", nrow(val$expr), ncol(val$expr))

# ---------------- GSE29431 / GSE42568 (unpaired validation) ----------------
v1 <- preprocess_validation("GSE29431", RAW_GSE29431_ESET, GPL570_SOFT,
                            meta_csv = file.path(META_DIR, "GSE29431_sample_classes.csv"),
                            platform = "GPL570", fig_prefix = "Figure_05_GSE29431")
saveRDS(v1, file.path(DIR_EXT_PROC, "GSE29431_processed.rds"))
write.csv(v1$meta, file.path(DIR_EXT_PROC, "GSE29431_metadata.csv"), row.names = FALSE)
logmsg("Saved GSE29431 processed object: %d genes x %d samples", nrow(v1$expr), ncol(v1$expr))

v2 <- preprocess_validation("GSE42568", RAW_GSE42568_ESET, GPL570_SOFT,
                            meta_csv = file.path(META_DIR, "GSE42568_sample_classes.csv"),
                            platform = "GPL570", fig_prefix = "Figure_05_GSE42568")
saveRDS(v2, file.path(DIR_EXT_PROC, "GSE42568_processed.rds"))
write.csv(v2$meta, file.path(DIR_EXT_PROC, "GSE42568_metadata.csv"), row.names = FALSE)
logmsg("Saved GSE42568 processed object: %d genes x %d samples", nrow(v2$expr), ncol(v2$expr))

logmsg("All three validation cohorts preprocessed and saved.")
