setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 01_load_GEO.R
# Load pre-downloaded GEO ExpressionSet objects
# =============================================================================

source("config/analysis_config.R")
load_required_packages()

all_gse <- c(DISCOVERY_GSE, VALIDATION_GSE)
cat("Loading pre-downloaded GEO datasets...\n")

geo_list <- list()

for (gse in all_gse) {
    rds_file <- file.path(RAW_DIR, paste0(gse, "_eset.rds"))
    if (!file.exists(rds_file)) {
        cat("File not found:", rds_file, "\n")
        cat("Attempting to download from GEO...\n")
        geo_list[[gse]] <- getGEO(gse, destdir = RAW_DIR, getGPL = TRUE)
        if (is.list(geo_list[[gse]])) {
            geo_list[[gse]] <- geo_list[[gse]][[1]]
        }
    } else {
        cat("Loading", rds_file, "\n")
        geo_list[[gse]] <- readRDS(rds_file)
    }
    
    cat(sprintf("  %s: Platform %s, %d samples, %d probes\n",
                gse, annotation(geo_list[[gse]]),
                ncol(exprs(geo_list[[gse]])),
                nrow(exprs(geo_list[[gse]]))))
}

# Save as combined list
saveRDS(geo_list, file.path(RAW_DIR, "geo_raw_list.rds"))
cat("\nAll datasets loaded and saved.\n")

