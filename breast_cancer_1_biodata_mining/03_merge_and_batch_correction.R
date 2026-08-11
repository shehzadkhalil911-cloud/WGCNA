###############################################################################
# Script 03: Merge Datasets and Batch Correction
# Purpose:
#   - Find common genes across all datasets
#   - Merge expression matrices
#   - Create batch and group variables
#   - ComBat batch correction
#   - Save merged data
###############################################################################

set.seed(123)

library(sva)
library(limma)
library(tidyverse)

base_dir <- normalizePath("breast_cancer_integrated_analysis")
data_processed_dir <- file.path(base_dir, "data_processed")
results_dir <- file.path(base_dir, "results")
dir.create(results_dir, showWarnings = FALSE)
dir.create(file.path(results_dir, "tables"), showWarnings = FALSE, recursive = TRUE)

# --- Load processed data ---
processed_data <- readRDS(file.path(data_processed_dir, "all_processed_data.rds"))
geo_ids <- names(processed_data)

cat("========================================================\n")
cat("MERGING AND BATCH CORRECTION\n")
cat("========================================================\n")

# --- Find common genes ---
gene_lists <- lapply(processed_data, function(x) rownames(x$exprs))
common_genes <- Reduce(intersect, gene_lists)

cat("\nCommon genes across all datasets:", length(common_genes), "\n")
for (gse_id in geo_ids) {
  cat("  ", gse_id, ":", length(gene_lists[[gse_id]]), "genes\n")
}

common_genes <- sort(common_genes)

# --- Merge expression matrices ---
merged_exprs <- NULL
merged_pheno <- NULL
merged_batch <- NULL
merged_group <- NULL

for (gse_id in geo_ids) {
  mat <- processed_data[[gse_id]]$exprs[common_genes, , drop = FALSE]
  pdata <- processed_data[[gse_id]]$pheno
  group <- processed_data[[gse_id]]$group
  
  colnames(mat) <- paste0(gse_id, "_", colnames(mat))
  
  merged_exprs <- cbind(merged_exprs, mat)
  merged_batch <- c(merged_batch, rep(gse_id, ncol(mat)))
  merged_group <- c(merged_group, group)
  
  # Add batch and group info to pheno
  pdata$dataset <- gse_id
  pdata$group <- group
  merged_pheno <- bind_rows(merged_pheno, pdata)
}

cat("\nMerged expression matrix dimensions:", dim(merged_exprs), "\n")
cat("Samples per dataset:\n")
print(table(merged_batch))
cat("\nGroups:\n")
print(table(merged_group))

# --- Remove constant/zero-variance genes ---
gene_var <- apply(merged_exprs, 1, var, na.rm = TRUE)
zero_var <- gene_var == 0 | is.na(gene_var)
if (any(zero_var)) {
  cat("\nRemoving", sum(zero_var), "zero-variance genes\n")
  merged_exprs <- merged_exprs[!zero_var, ]
}

# --- Save pre-batch-correction data ---
saveRDS(merged_exprs, file.path(data_processed_dir, "merged_exprs_before_batch.rds"))
write.csv(merged_exprs, file.path(data_processed_dir, "merged_exprs_before_batch.csv"))

saveRDS(merged_pheno, file.path(data_processed_dir, "merged_pheno.rds"))
write.csv(merged_pheno, file.path(data_processed_dir, "merged_pheno.csv"))

batch_table <- data.frame(
  Sample = colnames(merged_exprs),
  Batch = merged_batch,
  Group = merged_group,
  stringsAsFactors = FALSE
)
write.csv(batch_table, file.path(results_dir, "tables", "batch_info.csv"), row.names = FALSE)

cat("\nPre-batch-correction data saved.\n")

# --- ComBat batch correction ---
cat("\n========================================================\n")
cat("Running ComBat batch correction...\n")
cat("========================================================\n")

# Check if we have multiple batches
n_batches <- length(unique(merged_batch))
if (n_batches < 2) {
  cat("WARNING: Only one batch found. Skipping batch correction.\n")
  corrected_exprs <- merged_exprs
} else {
  cat("Number of batches:", n_batches, "\n")
  cat("Batch distribution:\n")
  print(table(merged_batch))
  
  # Create model matrix with biological covariate (group)
  mod <- model.matrix(~ factor(merged_group))
  
  # Run ComBat
  corrected_exprs <- tryCatch({
    ComBat(dat = merged_exprs, batch = merged_batch, 
           mod = mod, par.prior = TRUE, prior.plots = FALSE)
  }, error = function(e) {
    cat("ComBat with covariates failed:", e$message, "\n")
    cat("Retrying without covariates...\n")
    ComBat(dat = merged_exprs, batch = merged_batch, 
           par.prior = TRUE, prior.plots = FALSE)
  })
  
  cat("Batch correction complete. Corrected matrix dimensions:", 
      dim(corrected_exprs), "\n")
}

# --- Save post-batch-correction data ---
saveRDS(corrected_exprs, file.path(data_processed_dir, "merged_exprs_after_batch.rds"))
write.csv(corrected_exprs, file.path(data_processed_dir, "merged_exprs_after_batch.csv"))

# --- Save merged phenotype with batch info ---
merged_pheno$batch <- merged_batch
saveRDS(merged_pheno, file.path(data_processed_dir, "merged_pheno_with_batch.rds"))

# --- Summary stats ---
batch_summary <- data.frame(
  Dataset = names(table(merged_batch)),
  Samples = as.numeric(table(merged_batch)),
  Tumor = as.numeric(table(merged_batch[merged_group == "Tumor"])),
  Normal = as.numeric(table(merged_batch[merged_group == "Normal"])),
  stringsAsFactors = FALSE
)
print(batch_summary)
write.csv(batch_summary, file.path(results_dir, "tables", "batch_summary.csv"), 
          row.names = FALSE)

# --- Save batch correction metadata for report ---
batch_correction_info <- list(
  method = "ComBat (Empirical Bayes)",
  parameters = "par.prior = TRUE, biological covariate = group",
  n_genes_before = nrow(merged_exprs),
  n_genes_after = nrow(corrected_exprs),
  n_samples = ncol(merged_exprs),
  n_batches = n_batches,
  batches = names(table(merged_batch)),
  batch_sizes = as.numeric(table(merged_batch))
)
saveRDS(batch_correction_info, 
        file.path(data_processed_dir, "batch_correction_info.rds"))

cat("\nMerge and batch correction complete.\n")
cat("Output files:\n")
cat("  - data_processed/merged_exprs_before_batch.rds\n")
cat("  - data_processed/merged_exprs_after_batch.rds\n")
cat("  - data_processed/merged_pheno.rds\n")
cat("  - results/tables/batch_summary.csv\n")
