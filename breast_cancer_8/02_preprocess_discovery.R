# ============================================================
# 02 - Independent preprocessing: GSE70947 (discovery)
# Agilent GPL13607, 296 samples (148 paired), log2 scale.
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
logmsg("=== 02 Preprocessing GSE70947 (discovery) ===")

es <- readRDS(RAW_GSE70947_ESET)
pd <- pData(es)
exprs_raw <- exprs(es)   # 62976 x 296, log2, some NAs

# ---------- platform annotation from GEO SOFT file ----------
ann <- parse_soft_platform(GPL13607_SOFT)
stopifnot(nrow(ann) == nrow(exprs_raw))
stopifnot(as.character(ann$ID) == rownames(exprs_raw))

# ---------- sample metadata ----------
tissue <- as.character(pd$`tissue:ch1`)
groups <- ifelse(tissue == "tumor", "Tumor", "Normal")
names(groups) <- pd$geo_accession
patient <- sub("-(normal|tumor)$", "", pd$title, ignore.case = TRUE)
names(patient) <- pd$geo_accession

meta <- data.frame(
  sample_id = pd$geo_accession,
  dataset = "GSE70947",
  class = groups,
  patient_id = patient,
  title = pd$title,
  er = pd$`er:ch1`,
  her2 = pd$`her2:ch1`,
  age_in_years = as.numeric(pd$`age_in_years:ch1`),
  stringsAsFactors = FALSE
)
rownames(meta) <- meta$sample_id

# ---------- probe filtering ----------
is_control <- ann$ControlType %in% c(1, NA) & !is.na(ann$ControlType)
is_control[ann$ProbeName %in% c("GE_BrightCorner", "DarkCorner")] <- TRUE
logmsg("Control probes: %d", sum(is_control))
logmsg("Non-control probes: %d", sum(!is_control))

sym <- first_symbol(ann$GeneName)
valid <- validate_symbols(sym)
logmsg("Probes with official symbol: %d", sum(valid & !is_control))

keep_ann <- !is_control & valid
# missingness filter
frac_na <- apply(exprs_raw, 1, function(x) mean(is.na(x)))
keep_na <- frac_na <= MAX_NA_FRAC
logmsg("Probes passing missingness filter (<=%.0f%% NA): %d", MAX_NA_FRAC * 100, sum(keep_na))

keep <- keep_ann & keep_na
logmsg("Probes retained after all filters: %d", sum(keep))

expr_qc <- exprs_raw[keep, , drop = FALSE]
ann_keep <- ann[keep, , drop = FALSE]
sym_keep <- sym[keep]

# ---------- imputation of residual NAs ----------
logmsg("Remaining NAs before imputation: %d", sum(is.na(expr_qc)))
expr_imp <- impute_row_median(expr_qc)

# ---------- independent quantile normalization ----------
expr_norm <- quantile_normalize(expr_imp)

# ---------- probe -> gene collapse (max IQR) ----------
expr_gene <- collapse_probes_to_genes(expr_norm, sym_keep, method = COLLAPSE_METHOD)
logmsg("Unique genes after collapse: %d", nrow(expr_gene))

# ---------- outlier assessment ----------
out_df <- assess_outliers(expr_gene)
logmsg("Samples flagged as outliers: %d", sum(out_df$is_outlier))
meta$is_outlier <- out_df$is_outlier[match(meta$sample_id, out_df$sample_id)]

# ---------- preprocessing summary table ----------
tab5 <- data.frame(
  dataset = "GSE70947",
  platform = "GPL13607",
  n_probes_raw = nrow(exprs_raw),
  n_control_probes_removed = sum(is_control),
  n_probes_no_symbol_removed = sum(keep_ann == FALSE & !is_control),
  n_probes_high_missing_removed = sum(keep_na == FALSE & keep_ann),
  n_probes_after_qc = sum(keep),
  n_imputed_cells = sum(is.na(expr_qc)),
  n_genes_after_collapse = nrow(expr_gene),
  scale = "log2",
  normalization = "quantile (preprocessCore) + original GEO limma quantile",
  needs_log2 = FALSE,
  stringsAsFactors = FALSE
)
write.csv(tab5, file.path(DIR_TABLES, "Table_05a_Preprocessing_Summary_Discovery.csv"), row.names = FALSE)

# ---------- figures ----------
step <- function(nm, expr) {
  tryCatch({ expr; logmsg("figure step OK: %s", nm) },
           error = function(e) { logmsg("figure step FAILED: %s -- %s", nm, conditionMessage(e)); NULL })
}
qc <- qc_box_density(expr_qc, expr_gene, groups, "GSE70947 (discovery)")
step("fig03", save_tiff_600dpi(qc$box, "Figure_03_Discovery_Normalization_QC.tiff", width = 12, height = 6))
step("fig04", save_tiff_600dpi(qc$dens, "Figure_04_Discovery_Density_QC.tiff", width = 10, height = 5))

step("fig05", {
  p_pca <- pca_plot(expr_gene, groups, "GSE70947")
  save_tiff_600dpi(p_pca, "Figure_05_Discovery_PCA.tiff", width = 7, height = 6)
})

step("fig06", correlation_heatmap(expr_gene, groups, "Figure_06_Discovery_Correlation_Heatmap.tiff", "GSE70947"))
step("fig07", dendrogram_plot(expr_gene, groups, "Figure_07_Discovery_Dendrogram.tiff", "GSE70947"))

counts <- c(
  probes_raw = nrow(exprs_raw),
  controls_removed = sum(is_control),
  no_symbol_removed = sum(keep_ann == FALSE & !is_control),
  missing_removed = sum(keep_na == FALSE & keep_ann),
  probes_kept = sum(keep),
  genes_final = nrow(expr_gene)
)
step("fig08", retention_plot(counts, "Figure_03_Discovery_Probe_Retention.tiff", "GSE70947"))
step("fig09", outlier_plot(out_df, "Figure_03_Discovery_Outliers.tiff", "GSE70947"))

# ---------- save ----------
disc <- list(
  expr = expr_gene,                       # gene x sample, log2 quantile-normalized
  expr_before_norm = expr_imp,            # gene-collapse happened after norm; keep probe-level
  probe_expr = expr_norm,
  meta = meta,
  annotation = ann_keep,
  probe2gene = sym_keep,
  scale_info = list(min = min(expr_gene), max = max(expr_gene),
                    needs_log2 = FALSE, median = median(expr_gene)),
  preprocessing = tab5,
  outlier_assessment = out_df,
  dataset = "GSE70947"
)
saveRDS(disc, file.path(DIR_DISC_PROC, "discovery_processed.rds"))
write.csv(meta, file.path(DIR_DISC_PROC, "discovery_metadata.csv"), row.names = FALSE)

logmsg("Saved discovery processed object: %d genes x %d samples", nrow(expr_gene), ncol(expr_gene))
