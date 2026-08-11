###############################################################################
# Script 01: Download GEO Data
# Purpose: Download breast cancer GEO datasets using GEOquery
# Datasets: GSE42568, GSE15852, GSE29431, GSE70947
###############################################################################

set.seed(123)

# --- Load required packages ---
library(GEOquery)
library(Biobase)
library(tidyverse)

# --- Define project paths ---
base_dir <- normalizePath(file.path(".", ".."))
data_raw_dir <- file.path(base_dir, "data_raw")
data_processed_dir <- file.path(base_dir, "data_processed")

dir.create(data_raw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_processed_dir, showWarnings = FALSE, recursive = TRUE)

# --- Define datasets ---
geo_ids <- c("GSE42568", "GSE15852", "GSE29431", "GSE70947")

# --- Download function with retries ---
download_geo_dataset <- function(gse_id, dest_dir) {
  cat("\n========================================================\n")
  cat("Processing:", gse_id, "\n")
  cat("========================================================\n")

  # Try to get the series matrix file
  eset_list <- NULL
  n_attempts <- 3
  for (i in seq_len(n_attempts)) {
    eset_list <- tryCatch({
      getGEO(gse_id, destdir = dest_dir, AnnotGPL = FALSE, 
             getGPL = FALSE)
    }, error = function(e) {
      cat("  Attempt", i, "failed:", e$message, "\n")
      NULL
    })
    if (!is.null(eset_list)) break
    Sys.sleep(5)
  }

  if (is.null(eset_list)) {
    stop("Failed to download ", gse_id, " after ", n_attempts, " attempts.")
  }

  cat("  Successfully retrieved", gse_id, "\n")
  cat("  Number of platforms:", length(eset_list), "\n")

  if (length(eset_list) == 0) {
    stop("No ExpressionSets found for ", gse_id)
  }

  # Usually the first element is the main dataset
  eset <- eset_list[[1]]

  cat("  Samples:", ncol(eset), "\n")
  cat("  Features:", nrow(eset), "\n")

  return(eset)
}

# --- Download all datasets ---
geo_data <- list()

for (gse_id in geo_ids) {
  geo_data[[gse_id]] <- download_geo_dataset(gse_id, data_raw_dir)
}

# --- Save raw ExpressionSets as RDS ---
for (gse_id in geo_ids) {
  saveRDS(geo_data[[gse_id]], 
          file = file.path(data_raw_dir, paste0(gse_id, "_eset.rds")))
  cat("Saved:", gse_id, "\n")
}

# --- Create dataset summary ---
summary_list <- lapply(geo_ids, function(gse_id) {
  eset <- geo_data[[gse_id]]
  pdata <- pData(eset)
  data.frame(
    Dataset = gse_id,
    Samples = ncol(eset),
    Features = nrow(eset),
    Platform = unique(eset@annotation),
    stringsAsFactors = FALSE
  )
})
dataset_summary <- do.call(rbind, summary_list)
cat("\nDataset Summary:\n")
print(dataset_summary)

write.csv(dataset_summary, 
          file = file.path(data_processed_dir, "dataset_summary.csv"), 
          row.names = FALSE)

# --- Save session info ---
cat("\nDownload complete. Session info saved.\n")
