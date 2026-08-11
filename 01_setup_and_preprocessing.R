###############################################################################
# SCRIPT 1: Setup, Data Loading, and Preprocessing
# Multi-Cohort WGCNA: Breast Cancer Hub Gene Discovery
###############################################################################
library(WGCNA)
library(limma)
library(edgeR)
library(DESeq2)
library(sva)
library(biomaRt)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(hgu133plus2.db)
library(hgu133a.db)
library(impute)
library(preprocessCore)
library(flashClust)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(stringr)
library(gridExtra)
library(cowplot)
library(ggpubr)
library(RColorBrewer)
library(pROC)
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(pathview)
library(msigdbr)
library(GSVA)
library(DOSE)
library(circlize)
library(ggvenn)
library(tidyverse)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

BASE_DIR <- "F:/independent projects/cancer/breast cancer"
DATA_DIR <- file.path(BASE_DIR, "Datasets/data")
OUT_DIR <- file.path(BASE_DIR, "breast cancer 3 (WCGNA)")

dir.create(file.path(OUT_DIR, "figures/DEA"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "figures/WGCNA"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "figures/preservation"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "figures/enrichment"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "figures/validation"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "figures/immune"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "tables"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "networks"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(OUT_DIR, "RData"), recursive=TRUE, showWarnings=FALSE)

###############################################################################
# Load expression data and metadata
###############################################################################
load_expression <- function(csv_path) {
  mat <- as.matrix(read.csv(csv_path, row.names=1, check.names=FALSE))
  return(mat)
}

message("Loading GSE42568...")
expr_GSE42568 <- load_expression(file.path(DATA_DIR, "processed/GSE42568_expression.csv"))
meta_GSE42568 <- read.csv(file.path(DATA_DIR, "metadata/GSE42568_sample_classes.csv"))

message("Loading GSE15852...")
expr_GSE15852 <- load_expression(file.path(DATA_DIR, "processed/GSE15852_expression.csv"))
meta_GSE15852 <- read.csv(file.path(DATA_DIR, "metadata/GSE15852_sample_classes.csv"))

message("Loading GSE29431...")
expr_GSE29431 <- load_expression(file.path(DATA_DIR, "processed/GSE29431_expression.csv"))
meta_GSE29431 <- read.csv(file.path(DATA_DIR, "metadata/GSE29431_sample_classes.csv"))

message("Loading GSE70947...")
expr_GSE70947 <- load_expression(file.path(DATA_DIR, "processed/GSE70947_expression.csv"))
meta_GSE70947 <- read.csv(file.path(DATA_DIR, "metadata/GSE70947_sample_classes.csv"))

cat("\n=== SAMPLE COUNTS (from local phenotype files) ===\n")
for (ds in c("GSE42568", "GSE15852", "GSE29431", "GSE70947")) {
  meta <- get(paste0("meta_", ds))
  expr <- get(paste0("expr_", ds))
  n_normal <- sum(meta$class == "Normal")
  n_tumor <- sum(meta$class == "Tumor")
  cat(sprintf("%s: %d Normal, %d Tumor, %d total, %d probes\n", 
              ds, n_normal, n_tumor, ncol(expr), nrow(expr)))
}

###############################################################################
# Preprocessing function
###############################################################################
preprocess_dataset <- function(expr_mat, meta_df, dataset_name,
                                collapse_genes = TRUE,
                                min_expr_filter = TRUE,
                                min_var_filter = TRUE) {

  common_samples <- intersect(colnames(expr_mat), meta_df$sample_id)
  cat(sprintf("\n%s: %d/%d samples matched\n", dataset_name,
              length(common_samples), ncol(expr_mat)))
  expr_mat <- expr_mat[, common_samples, drop = FALSE]
  meta_df <- meta_df[match(common_samples, meta_df$sample_id), ]
  rownames(meta_df) <- meta_df$sample_id
  meta_df$class_binary <- ifelse(meta_df$class == "Tumor", 1, 0)

  # Remove control probes (for Agilent platform)
  control_pattern <- "^GE_BrightCorner|^DarkCorner|^ERCC|^NC1_|^NC2_"
  if (any(grepl(control_pattern, rownames(expr_mat), ignore.case = TRUE))) {
    expr_mat <- expr_mat[!grepl(control_pattern, rownames(expr_mat), ignore.case = TRUE), ]
    cat(sprintf("Removed control probes: %d remaining\n", nrow(expr_mat)))
  }

  # Log2 if needed (check if data is already log-scale)
  if (max(expr_mat, na.rm = TRUE) > 50) {
    cat("Data appears non-log2, applying log2 transformation...\n")
    expr_mat <- log2(expr_mat + 1)
  }

  # Quantile normalize
  cat("Quantile normalizing...\n")
  expr_norm <- normalize.quantiles(expr_mat)
  rownames(expr_norm) <- rownames(expr_mat)
  colnames(expr_norm) <- colnames(expr_mat)

  # For discovery datasets: row names are already gene symbols, just deduplicate
  if (collapse_genes) {
    cat("Checking for duplicate gene symbols...\n")
    gene_means <- rowMeans(expr_norm, na.rm = TRUE)
    dup_idx <- duplicated(rownames(expr_norm))
    if (any(dup_idx)) {
      # Keep the gene with highest mean expression
      keep_idx <- !dup_idx
      dup_genes <- unique(rownames(expr_norm)[dup_idx])
      for (g in dup_genes) {
        idx <- which(rownames(expr_norm) == g)
        best <- idx[which.max(gene_means[idx])]
        keep_idx[best] <- TRUE
        keep_idx[idx[idx != best]] <- FALSE
      }
      expr_norm <- expr_norm[keep_idx, ]
      cat(sprintf("After deduplication: %d genes\n", nrow(expr_norm)))
    }
  }

  if (min_expr_filter) {
    expr_cutoff <- quantile(rowMeans(expr_norm, na.rm = TRUE), probs = 0.1)
    expr_norm <- expr_norm[rowMeans(expr_norm, na.rm = TRUE) > expr_cutoff, ]
    cat(sprintf("After expression filter (>%0.2f): %d genes\n", expr_cutoff, nrow(expr_norm)))
  }

  if (min_var_filter) {
    gene_vars <- apply(expr_norm, 1, var, na.rm = TRUE)
    var_cutoff <- quantile(gene_vars, probs = 0.25, na.rm = TRUE)
    expr_norm <- expr_norm[gene_vars > var_cutoff, ]
    cat(sprintf("After variance filter (>%0.4f): %d genes\n", var_cutoff, nrow(expr_norm)))
  }

  na_idx <- apply(is.na(expr_norm), 1, any)
  if (any(na_idx)) {
    expr_norm <- expr_norm[!na_idx, ]
    cat(sprintf("After NA removal: %d genes\n", nrow(expr_norm)))
  }

  list(expr = expr_norm, meta = meta_df)
}

cat("\n\n========== PREPROCESSING DISCOVERY DATASETS ==========\n")
processed_42568 <- preprocess_dataset(expr_GSE42568, meta_GSE42568, "GSE42568")
processed_15852 <- preprocess_dataset(expr_GSE15852, meta_GSE15852, "GSE15852")
processed_29431 <- preprocess_dataset(expr_GSE29431, meta_GSE29431, "GSE29431")

cat("\n\n========== PREPROCESSING VALIDATION DATASET ==========\n")
processed_70947 <- preprocess_dataset(expr_GSE70947, meta_GSE70947, "GSE70947",
                                       collapse_genes = TRUE)

# For GSE70947, many probes already have gene symbols; keep those with valid symbols
expr_val <- processed_70947$expr
# Remove probes that are clearly not gene symbols (start with numbers, lincRNA, etc.)
# but keep everything - we'll intersect with discovery genes anyway
processed_70947$expr <- expr_val
cat(sprintf("GSE70947 final genes: %d\n", nrow(processed_70947$expr)))

# Find common genes across all datasets
common_genes_discovery <- Reduce(intersect, list(
  rownames(processed_42568$expr),
  rownames(processed_15852$expr),
  rownames(processed_29431$expr)
))
cat(sprintf("\nCommon genes across 3 discovery datasets: %d\n", length(common_genes_discovery)))

common_genes_all <- intersect(common_genes_discovery, rownames(processed_70947$expr))
cat(sprintf("Common genes across all 4 datasets: %d\n", length(common_genes_all)))

# Save full filtered datasets for individual WGCNA
save(processed_42568, processed_15852, processed_29431, processed_70947,
     common_genes_discovery, common_genes_all,
     file = file.path(OUT_DIR, "RData/01_processed_data_full.RData"))

# Create common-gene subset versions for cross-dataset comparisons
comm_42568 <- list(expr = processed_42568$expr[common_genes_all, ], meta = processed_42568$meta)
comm_15852 <- list(expr = processed_15852$expr[common_genes_all, ], meta = processed_15852$meta)
comm_29431 <- list(expr = processed_29431$expr[common_genes_all, ], meta = processed_29431$meta)
comm_70947 <- list(expr = processed_70947$expr[common_genes_all, ], meta = processed_70947$meta)

save(comm_42568, comm_15852, comm_29431, comm_70947,
     common_genes_discovery, common_genes_all,
     file = file.path(OUT_DIR, "RData/01_processed_data_common.RData"))

cat(sprintf("\nFinal counts per dataset (full):"))
cat(sprintf("\nGSE42568: %d genes, %d samples", nrow(processed_42568$expr), ncol(processed_42568$expr)))
cat(sprintf("\nGSE15852: %d genes, %d samples", nrow(processed_15852$expr), ncol(processed_15852$expr)))
cat(sprintf("\nGSE29431: %d genes, %d samples", nrow(processed_29431$expr), ncol(processed_29431$expr)))
cat(sprintf("\nCommon genes (discovery): %d", length(common_genes_discovery)))
cat(sprintf("\nCommon genes (all 4): %d", length(common_genes_all)))
cat("\n\nScript 1 complete. Processed data saved.\n")
