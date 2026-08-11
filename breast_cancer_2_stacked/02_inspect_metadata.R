setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 02_inspect_metadata.R
# Inspect all four GEO datasets for suitability
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")
load_required_packages()

geo_list <- readRDS(file.path(RAW_DIR, "geo_raw_list.rds"))

# Comprehensive inspection for each dataset
dataset_info <- list()

for (gse in names(geo_list)) {
    cat("\n", paste(rep("=", 60), collapse = ""), "\n")
    cat("INSPECTING:", gse, "\n")
    cat(paste(rep("=", 60), collapse = ""), "\n")

    gse_obj <- geo_list[[gse]]

    # Check if it's a list (multiple platforms)
    if (is.list(gse_obj) && length(gse_obj) > 1) {
        cat("Multiple platforms detected. Checking each...\n")
    }
    if (is.list(gse_obj)) {
        gse_obj <- gse_obj[[1]]  # Use first/platform
    }

    # Get expression data
    expr <- exprs(gse_obj)
    pdata <- pData(gse_obj)
    fdata <- fData(gse_obj)
    platform <- annotation(gse_obj)

    cat("Platform:", platform, "\n")
    cat("Samples:", ncol(expr), "\n")
    cat("Probes:", nrow(expr), "\n")

    # Examine phenotype data
    cat("\nPhenotype columns:\n")
    for (col in colnames(pdata)) {
        vals <- unique(as.character(pdata[[col]]))
        vals <- vals[!is.na(vals) & vals != ""]
        if (length(vals) <= 20) {
            cat("  ", col, ":", paste(head(vals, 10), collapse = ", "),
                if (length(vals) > 10) "...", "\n")
        }
    }

    # Try to identify tumor/normal columns
    cat("\nSearching for disease/group columns...\n")
    for (col in colnames(pdata)) {
        vals <- tolower(unique(as.character(pdata[[col]])))
        if (any(grepl("tumor|cancer|normal|disease|breast|case|control", vals))) {
            cat("  Possible class column:", col, "\n")
            cat("  Values:", paste(unique(as.character(pdata[[col]])), collapse = " | "), "\n")
        }
    }

    # Expression summary
    vals <- as.numeric(as.matrix(expr))
    cat("\nExpression range:", range(vals, na.rm = TRUE), "\n")
    cat("Expression median:", median(vals, na.rm = TRUE), "\n")
    cat("Has negative values:", any(vals < 0, na.rm = TRUE), "\n")

    dataset_info[[gse]] <- list(
        platform = platform,
        n_samples = ncol(expr),
        n_probes = nrow(expr),
        expression_range = range(vals, na.rm = TRUE),
        expression_median = median(vals, na.rm = TRUE),
        has_negative = any(vals < 0, na.rm = TRUE),
        phenotype_cols = colnames(pdata),
        sample_names = sampleNames(gse_obj)
    )
}

# Save inspection results
saveRDS(dataset_info, file.path(METADATA_DIR, "dataset_inspection.rds"))
cat("\nInspection complete. Results saved.\n")

# Print summary table
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat(sprintf("%-15s %-12s %-10s %-10s %-15s %-10s\n",
            "Dataset", "Platform", "Samples", "Probes", "Expression Range", "Has Neg"))
cat(paste(rep("-", 80), collapse = ""), "\n")
for (gse in names(dataset_info)) {
    info <- dataset_info[[gse]]
    cat(sprintf("%-15s %-12s %-10d %-10d %-15s %-10s\n",
                gse,
                info$platform,
                info$n_samples,
                info$n_probes,
                paste(round(info$expression_range, 1), collapse = "-"),
                info$has_negative))
}
cat(paste(rep("=", 80), collapse = ""), "\n")

