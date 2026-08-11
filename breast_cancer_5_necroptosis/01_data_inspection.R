# ============================================================
# Step 1: Data Inspection and Suitability Report
# ============================================================
set.seed(123)
library(readxl)
library(Biobase)
library(limma)
library(dplyr)
library(tidyr)

base_dir <- normalizePath(file.path(".", "Breast Cancer 5 Necroptosis"))
data_dir <- normalizePath(file.path(".", "Datasets", "data"))
output_dir <- base_dir

log_file <- file.path(output_dir, "logs", "01_data_inspection.log")
sink(log_file, split=TRUE)

cat("========================================\n")
cat("DATA INSPECTION REPORT\n")
cat("========================================\n\n")

cat("Date:", date(), "\n\n")

# 1.1 Input file manifest
cat("--- INPUT FILE MANIFEST ---\n\n")

raw_files <- list.files(file.path(data_dir, "raw"), full.names=TRUE)
meta_files <- list.files(file.path(data_dir, "metadata"), full.names=TRUE)
proc_files <- list.files(file.path(data_dir, "processed"), full.names=TRUE)
disc_files <- list.files(file.path(data_dir, "discovery"), full.names=TRUE)
val_files <- list.files(file.path(data_dir, "external_validation"), full.names=TRUE)

cat("RAW DATA FILES:\n")
for(f in raw_files) {
  info <- file.info(f)
  cat("  ", basename(f), " | Size:", round(info$size/1024, 1), "KB | Modified:", info$mtime, "\n")
}
cat("\nMETADATA FILES:\n")
for(f in meta_files) {
  info <- file.info(f)
  cat("  ", basename(f), " | Size:", round(info$size/1024, 1), "KB\n")
}
cat("\nPROCESSED FILES:\n")
for(f in proc_files) {
  info <- file.info(f)
  cat("  ", basename(f), " | Size:", round(info$size/1024, 1), "KB\n")
}
cat("\nDISCOVERY FILES:\n")
for(f in disc_files) {
  info <- file.info(f)
  cat("  ", basename(f), " | Size:", round(info$size/1024, 1), "KB\n")
}
cat("\nVALIDATION FILES:\n")
for(f in val_files) {
  info <- file.info(f)
  cat("  ", basename(f), " | Size:", round(info$size/1024, 1), "KB\n")
}

# 1.2 Necroptosis gene list
cat("\n--- NECROPTOSIS GENE LIST ---\n\n")
ng <- read_excel("Necroptosis genes.xlsx", sheet=1)
ng_genes <- c(colnames(ng)[1], ng[,1,drop=TRUE])
cat("Raw necroptosis genes file:\n")
cat("  File: Necroptosis genes.xlsx\n")
cat("  Sheets:", excel_sheets("Necroptosis genes.xlsx"), "\n")
cat("  Rows:", nrow(ng), "Cols:", ncol(ng), "\n")
cat("  Total unique genes:", length(unique(ng_genes)), "\n")

# Standardize gene symbols
suppressPackageStartupMessages(library(org.Hs.eg.db))
gene_std <- mapIds(org.Hs.eg.db, keys=ng_genes, keytype="SYMBOL", column="SYMBOL", multiVals="first")
unmapped <- ng_genes[is.na(gene_std)]
cat("\n  Mapped to official symbols:", sum(!is.na(gene_std)), "\n")
cat("  Unmapped:", length(unmapped), "\n")
if(length(unmapped) > 0) cat("  Unmapped genes:", paste(unmapped, collapse=", "), "\n")
cat("  Duplicates in original:", sum(duplicated(ng_genes)), "\n")
cat("  Unique official symbols:", length(unique(gene_std[!is.na(gene_std)])), "\n")

ng_official <- unique(gene_std[!is.na(gene_std)])
saveRDS(ng_official, file.path(output_dir, "metadata", "necroptosis_genes_official.rds"))
writeLines(ng_official, file.path(output_dir, "metadata", "necroptosis_genes_official.txt"))
cat("  Saved necroptosis gene list to metadata/\n")

# 1.3 Sample metadata inspection
cat("\n\n--- SAMPLE METADATA INSPECTION ---\n\n")
datasets <- c("GSE42568", "GSE15852", "GSE29431", "GSE70947")
meta_summary <- data.frame()

for(ds in datasets) {
  m <- read.csv(file.path(data_dir, "metadata", paste0(ds, "_sample_classes.csv")))
  cat(ds, ":\n")
  cat("  Total samples:", nrow(m), "\n")
  cat("  Columns:", paste(colnames(m), collapse=", "), "\n")
  cat("  Class distribution:\n")
  tbl <- table(m$class)
  for(nm in names(tbl)) cat("    ", nm, ": ", tbl[nm], "\n", sep="")
  
  # Check duplicate sample IDs
  dup <- duplicated(m$sample_id)
  cat("  Duplicate sample IDs:", sum(dup), "\n")
  
  # Check missing metadata
  missing <- sapply(m, function(x) sum(is.na(x) | x == ""))
  cat("  Missing values:", paste(names(missing), missing, sep="=", collapse=", "), "\n")
  
  # Expression data check
  csv_file <- file.path(data_dir, "processed", paste0(ds, "_expression.csv"))
  if(file.exists(csv_file)) {
    e <- read.csv(csv_file, row.names=1, check.names=FALSE)
    cat("  Expression matrix:", nrow(e), "genes x", ncol(e), "samples\n")
    common <- intersect(colnames(e), m$sample_id)
    cat("  Matching sample IDs:", length(common), "/", ncol(e), "\n")
    if(length(common) < ncol(e)) {
      missing_s <- setdiff(colnames(e), m$sample_id)
      cat("  WARNING: Samples in expression not in metadata:", 
          paste(head(missing_s, 5), collapse=", "), 
          if(length(missing_s)>5) paste("... (+", length(missing_s)-5, "more)"), "\n")
    }
  }
  
  # RDS raw data check
  rds_file <- file.path(data_dir, "raw", paste0(ds, "_eset.rds"))
  if(file.exists(rds_file)) {
    eset <- readRDS(rds_file)
    cat("  Raw ExpressionSet:", nrow(exprs(eset)), "probes x", ncol(exprs(eset)), "samples\n")
    cat("  Platform:", notes(eset)$platform, "\n")
  }
  
  meta_summary <- rbind(meta_summary, data.frame(
    Dataset = ds,
    Total_Samples = nrow(m),
    Normal = sum(m$class == "Normal"),
    Tumor = sum(m$class == "Tumor"),
    Duplicates = sum(dup),
    Expression_Genes = if(file.exists(csv_file)) nrow(e) else NA,
    Raw_Probes = if(file.exists(rds_file)) nrow(exprs(eset)) else NA,
    stringsAsFactors = FALSE
  ))
  cat("\n")
}

cat("\n--- DATASET SUMMARY TABLE ---\n")
print(meta_summary)
write.csv(meta_summary, file.path(output_dir, "metadata", "dataset_summary.csv"), row.names=FALSE)

# 1.4 Cohort assignment
cat("\n--- COHORT ASSIGNMENT ---\n")
cat("Discovery cohort: GSE42568, GSE15852, GSE29431\n")
cat("Independent validation: GSE70947\n")
cat("Validation dataset will NOT be used in discovery batch correction, DE, feature selection, or model training.\n")

# 1.5 Discovery data loading
cat("\n--- LOADING DISCOVERY DATA ---\n")
discovery_expr <- list()
discovery_meta <- list()
for(ds in c("GSE42568", "GSE15852", "GSE29431")) {
  csv_file <- file.path(data_dir, "processed", paste0(ds, "_expression.csv"))
  meta_file <- file.path(data_dir, "metadata", paste0(ds, "_sample_classes.csv"))
  e <- read.csv(csv_file, row.names=1, check.names=FALSE)
  m <- read.csv(meta_file)
  discovery_expr[[ds]] <- e
  discovery_meta[[ds]] <- m
  cat(ds, ": loaded", nrow(e), "genes x", ncol(e), "samples\n")
}

# Load validation data
val_csv <- file.path(data_dir, "processed", "GSE70947_expression.csv")
val_meta_file <- file.path(data_dir, "metadata", "GSE70947_sample_classes.csv")
val_expr <- read.csv(val_csv, row.names=1, check.names=FALSE)
val_meta <- read.csv(val_meta_file)
cat("GSE70947 (validation):", nrow(val_expr), "genes x", ncol(val_expr), "samples\n")

# Find common genes
common_genes <- Reduce(intersect, lapply(discovery_expr, rownames))
cat("\nCommon genes across all 3 discovery datasets:", length(common_genes), "\n")

saveRDS(list(expression=discovery_expr, metadata=discovery_meta, common_genes=common_genes,
             val_expr=val_expr, val_meta=val_meta),
        file.path(output_dir, "data_processed", "inspection_data.rds"))

cat("\n--- INSPECTION COMPLETE ---\n")
sink()
