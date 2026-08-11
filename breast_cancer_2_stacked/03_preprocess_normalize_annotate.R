setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 03_preprocess_normalize_annotate.R
# Preprocess, normalize, and annotate each GEO dataset separately
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")
source("R/preprocessing_functions.R")
source("R/annotation_functions.R")
source("R/plotting_functions.R")
load_required_packages()

geo_list <- readRDS(file.path(RAW_DIR, "geo_raw_list.rds"))

processed_list <- list()

for (gse in names(geo_list)) {
    cat("\n", paste(rep("=", 70), collapse = ""), "\n")
    cat("PROCESSING:", gse, "\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")

    gse_obj <- geo_list[[gse]]
    if (is.list(gse_obj)) gse_obj <- gse_obj[[1]]

    expr_orig <- exprs(gse_obj)
    pdata <- pData(gse_obj)
    platform_id <- annotation(gse_obj)

    cat("Platform:", platform_id, "\n")

    # --- Expression scale assessment ---
    scale_info <- assess_log2_transform(expr_orig)
    cat(sprintf("Range: %.1f - %.1f, Median: %.1f\n",
                scale_info$min, scale_info$max, scale_info$median))
    cat("Needs log2 transform:", scale_info$needs_log2, "\n")

    # Apply log2 if needed
    if (scale_info$needs_log2) {
        cat("Applying log2 transformation...\n")
        expr <- safe_log2(expr_orig)
    } else {
        expr <- as.matrix(expr_orig)
    }

    # --- Normalization ---
    cat("Applying quantile normalization...\n")
    expr_norm <- quantile_normalize(expr)

    # --- Annotation ---
    cat("Fetching platform annotation...\n")
    gpl_table <- get_platform_annotation(platform_id)

    cat("Mapping probes to genes...\n")
    gene_map <- map_probes_to_genes(rownames(expr_norm), gpl_table)

    ann_summary <- annotate_summary(gene_map)
    cat(sprintf("Mapped %d/%d probes (%.1f%%)\n",
                ann_summary$mapped_genes, ann_summary$total_probes,
                ann_summary$mapping_rate * 100))

    # Collapse duplicates
    cat("Collapsing duplicate probes (method:", DUPLICATE_PROBE_METHOD, ")...\n")
    expr_gene <- collapse_duplicates(expr_norm, gene_map, method = DUPLICATE_PROBE_METHOD)

    cat("Final gene-level matrix:", nrow(expr_gene), "genes x", ncol(expr_gene), "samples\n")

    processed_list[[gse]] <- list(
        expression = expr_gene,
        expression_raw = expr_orig,
        expression_norm = expr_norm,
        phenotype = pdata,
        platform = platform_id,
        annotation = gpl_table,
        gene_map = gene_map,
        scale_info = scale_info,
        gse_obj = gse_obj
    )

    # Save processed
    saveRDS(processed_list[[gse]], file.path(PROCESSED_DIR, paste0(gse, "_processed.rds")))
    fwrite(as.data.frame(expr_gene), file.path(PROCESSED_DIR, paste0(gse, "_expression.csv")),
           row.names = TRUE)
}

cat("\nAll datasets preprocessed and saved.\n")

