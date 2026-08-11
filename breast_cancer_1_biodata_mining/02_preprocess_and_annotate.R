###############################################################################
# Script 02: Preprocess and Annotate GEO Datasets
# Purpose: 
#   - Check platform information
#   - Extract expression matrices and phenotype metadata
#   - Identify tumor vs normal samples
#   - Probe-to-gene symbol annotation
#   - Collapse duplicate genes
#   - Log2 transformation and normalization
#   - Save individual processed data
###############################################################################

set.seed(123)

library(GEOquery)
library(Biobase)
library(limma)
library(annotate)
library(AnnotationDbi)
library(hgu133a.db)
library(hgu133plus2.db)
library(org.Hs.eg.db)
library(tidyverse)

# --- Paths ---
base_dir <- normalizePath("breast_cancer_integrated_analysis")
data_raw_dir <- file.path(base_dir, "data_raw")
data_processed_dir <- file.path(base_dir, "data_processed")
dir.create(data_processed_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helper: Check if data needs log2 transformation ---
needs_log2_transform <- function(exprs_mat) {
  vals <- as.numeric(exprs_mat)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) return(TRUE)
  q75 <- quantile(vals, 0.75, na.rm = TRUE)
  q25 <- quantile(vals, 0.25, na.rm = TRUE)
  iqr <- q75 - q25
  median_val <- median(vals, na.rm = TRUE)
  if (iqr == 0) return(TRUE)
  max_val <- max(vals, na.rm = TRUE)
  if (max_val > 100) {
    return(TRUE)
  }
  FALSE
}

# --- Helper: log2 transform ---
safe_log2 <- function(x) {
  x[x <= 0] <- NA
  log2(x)
}

# --- Helper: Probe-to-gene annotation ---
map_probes_to_genes <- function(exprs_mat, platform_pkg) {
  require(platform_pkg, character.only = TRUE)
  
  probe_ids <- rownames(exprs_mat)
  
  if (platform_pkg == "hgu133a.db") {
    gene_symbols <- mapIds(hgu133a.db, keys = probe_ids, 
                           column = "SYMBOL", keytype = "PROBEID", 
                           multiVals = "first")
  } else if (platform_pkg == "hgu133plus2.db") {
    gene_symbols <- mapIds(hgu133plus2.db, keys = probe_ids, 
                           column = "SYMBOL", keytype = "PROBEID", 
                           multiVals = "first")
  } else {
    stop("Unknown platform package: ", platform_pkg)
  }
  
  cat("  Mapped", sum(!is.na(gene_symbols)), "out of", length(probe_ids), 
      "probes to gene symbols\n")
  
  list(gene_symbols = gene_symbols, exprs_mat = exprs_mat)
}

# --- Helper: Collapse duplicate genes (keep highest mean expression) ---
collapse_duplicate_genes <- function(exprs_mat, gene_symbols, method = "max_mean") {
  valid <- !is.na(gene_symbols) & gene_symbols != ""
  exprs_mat <- exprs_mat[valid, , drop = FALSE]
  gene_symbols <- gene_symbols[valid]
  
  # Remove control probes (AFFX prefix)
  affx_idx <- grep("^AFFX", gene_symbols, ignore.case = TRUE)
  if (length(affx_idx) > 0) {
    exprs_mat <- exprs_mat[-affx_idx, , drop = FALSE]
    gene_symbols <- gene_symbols[-affx_idx]
    cat("  Removed", length(affx_idx), "control probes\n")
  }
  
  if (method == "max_mean") {
    gene_means <- rowMeans(exprs_mat, na.rm = TRUE)
    keep_idx <- !duplicated(gene_symbols)
    if (any(duplicated(gene_symbols))) {
      for (g in unique(gene_symbols[duplicated(gene_symbols)])) {
        idx <- which(gene_symbols == g)
        best <- idx[which.max(gene_means[idx])]
        keep_idx[idx] <- FALSE
        keep_idx[best] <- TRUE
      }
    }
  } else if (method == "max_var") {
    gene_vars <- apply(exprs_mat, 1, var, na.rm = TRUE)
    keep_idx <- !duplicated(gene_symbols)
    if (any(duplicated(gene_symbols))) {
      for (g in unique(gene_symbols[duplicated(gene_symbols)])) {
        idx <- which(gene_symbols == g)
        best <- idx[which.max(gene_vars[idx])]
        keep_idx[idx] <- FALSE
        keep_idx[best] <- TRUE
      }
    }
  } else {
    stop("Unknown collapse method: ", method)
  }
  
  exprs_mat <- exprs_mat[keep_idx, , drop = FALSE]
  gene_symbols <- gene_symbols[keep_idx]
  rownames(exprs_mat) <- gene_symbols
  
  cat("  Final genes after collapsing:", nrow(exprs_mat), "\n")
  
  exprs_mat
}

# --- Helper: Identify tumor vs normal samples ---
identify_groups <- function(pdata, dataset_id) {
  cat("\n  ---- Phenotype columns for", dataset_id, "----\n")
  cat("  Column names:", paste(colnames(pdata), collapse = ", "), "\n")
  
  # Search for relevant columns
  possible_cols <- c("title", "characteristics_ch1", "source_name_ch1")
  found_cols <- intersect(possible_cols, colnames(pdata))
  
  cat("  Found potential annotation columns:", paste(found_cols, collapse = ", "), "\n")
  
  for (col in found_cols) {
    cat("  Values in", col, ":\n")
    vals <- table(as.character(pdata[[col]]))
    print(head(vals, 30))
  }
  
  group <- rep(NA, nrow(pdata))
  
  # Priority 1: Title column is most reliable for clean classification
  if ("title" %in% colnames(pdata)) {
    titles <- trimws(as.character(pdata$title))
    for (i in seq_along(titles)) {
      t <- tolower(titles[i])
      if (grepl("^n[0-9]", t) || grepl("normal", t, fixed = TRUE) || grepl("control", t, fixed = TRUE)) {
        group[i] <- "Normal"
      } else if (grepl("^t[0-9]", t) || grepl("tumor", t, fixed = TRUE) || 
                 grepl("cancer", t, fixed = TRUE) || grepl("carcinoma", t, fixed = TRUE)) {
        group[i] <- "Tumor"
      }
    }
  }
  
  # Priority 2: Source name (check normal first to avoid "cancer" false matches)
  if ("source_name_ch1" %in% colnames(pdata)) {
    src <- trimws(as.character(pdata$source_name_ch1))
    for (i in seq_along(src)) {
      if (is.na(group[i])) {
        t <- tolower(src[i])
        if (grepl("adjacent to tumor", t) || grepl("adjacent non", t) || 
            grepl("normal adjacent", t)) {
          group[i] <- "Normal"
        } else if (grepl("normal", t, fixed = TRUE)) {
          group[i] <- "Normal"
        } else if (grepl("tumor", t, fixed = TRUE)) {
          group[i] <- "Tumor"
        }
      }
    }
  }
  
  # Priority 3: Characteristics
  if ("characteristics_ch1" %in% colnames(pdata)) {
    chars <- trimws(as.character(pdata$characteristics_ch1))
    for (i in seq_along(chars)) {
      if (!is.na(group[i])) next
      t <- tolower(chars[i])
      if (grepl("normal", t, fixed = TRUE) || grepl("control", t, fixed = TRUE)) {
        group[i] <- "Normal"
      } else if (grepl("tumor", t, fixed = TRUE) || grepl("cancer", t, fixed = TRUE) ||
                 grepl("carcinoma", t, fixed = TRUE)) {
        group[i] <- "Tumor"
      }
    }
  }
  
  # Dataset-specific: GSE70947 "adjacent to tumor" samples are Normal
  if (dataset_id == "GSE70947" && "source_name_ch1" %in% colnames(pdata)) {
    src <- trimws(as.character(pdata$source_name_ch1))
    for (i in seq_along(src)) {
      t <- tolower(src[i])
      if (grepl("adjacent", t)) {
        group[i] <- "Normal"
      }
    }
  }
  
  group
}

# --- Download/load GEO datasets ---
geo_ids <- c("GSE42568", "GSE15852", "GSE29431", "GSE70947")
geo_data <- list()

for (gse_id in geo_ids) {
  rds_file <- file.path(data_raw_dir, paste0(gse_id, "_eset.rds"))
  if (file.exists(rds_file)) {
    cat("\nLoading cached", gse_id, "...\n")
    geo_data[[gse_id]] <- readRDS(rds_file)
  } else {
    cat("\nDownloading", gse_id, "...\n")
    geo_data[[gse_id]] <- getGEO(gse_id, destdir = data_raw_dir, 
                                 AnnotGPL = FALSE, getGPL = FALSE)[[1]]
    saveRDS(geo_data[[gse_id]], rds_file)
  }
}

# --- Process each dataset ---
platform_map <- list(
  GSE15852 = list(pkg = "hgu133a.db", platform_id = "GPL96"),
  GSE42568 = list(pkg = "hgu133plus2.db", platform_id = "GPL570"),
  GSE29431 = list(pkg = "hgu133plus2.db", platform_id = "GPL570"),
  GSE70947 = list(pkg = NULL, platform_id = "GPL13607")
)

processed_data <- list()
sample_group_summary <- list()

for (gse_id in geo_ids) {
  cat("\n########################################################\n")
  cat("Processing:", gse_id, "\n")
  cat("########################################################\n")
  
  eset <- geo_data[[gse_id]]
  eset_platform <- eset@annotation
  
  cat("  Platform:", eset_platform, "\n")
  cat("  Samples:", ncol(eset), "\n")
  cat("  Probes:", nrow(eset), "\n")
  
  # Extract expression matrix
  exprs_mat <- exprs(eset)
  
  # Check log2 transformation
  if (needs_log2_transform(exprs_mat)) {
    cat("  Applying log2 transformation\n")
    exprs_mat <- safe_log2(exprs_mat)
  } else {
    cat("  Data appears already log2 transformed\n")
  }
  
  # Extract phenotype data
  pdata <- pData(eset)
  cat("\n  Phenotype columns:\n")
  cat("  ", paste(colnames(pdata), collapse = "\n  "), "\n")
  
  # Identify sample groups
  group <- identify_groups(pdata, gse_id)
  pdata$group <- group
  
  cat("\n  Group assignment:\n")
  group_table <- table(group, useNA = "ifany")
  print(group_table)
  
  # Handle ambiguous groups
  n_na <- sum(is.na(group))
  if (n_na > 0) {
    cat("\n  WARNING:", n_na, "samples have unassigned group labels\n")
    cat("  Creating group inspection table\n")
    
    inspection_table <- data.frame(
      Sample = sampleNames(eset),
      Title = as.character(pdata$title),
      Source = as.character(pdata$source_name_ch1 %||% NA),
      Characteristics = as.character(pdata$characteristics_ch1 %||% NA),
      Assigned_Group = group,
      stringsAsFactors = FALSE
    )
    
    write.csv(inspection_table,
              file = file.path(data_processed_dir, 
                               paste0(gse_id, "_group_inspection.csv")),
              row.names = FALSE)
    
    cat("  Inspection table saved to:", 
        file.path(data_processed_dir, paste0(gse_id, "_group_inspection.csv")), "\n")
  }
  
  # Remove samples with missing group labels  
  valid_samples <- !is.na(group)
  cat("\n  Valid samples with group labels:", sum(valid_samples), "/", length(group), "\n")
  
  exprs_mat <- exprs_mat[, valid_samples, drop = FALSE]
  pdata <- pdata[valid_samples, , drop = FALSE]
  group <- group[valid_samples]
  
  # For GSE70947 (Agilent), handle differently
  if (gse_id == "GSE70947") {
    cat("\n  GSE70947 uses Agilent platform - attempting custom annotation\n")
    
    # Try using the GPL annotation from the ExpressionSet
    cat("  Using feature annotation from ExpressionSet\n")
    
    fData <- fData(eset)
    cat("  Feature data columns:", paste(colnames(fData), collapse = ", "), "\n")
    
    # Check for gene annotation in fData
    gene_col <- NULL
    possible_gene_cols <- c("Gene Symbol", "GENE_SYMBOL", "Symbol", "GENE", 
                            "Gene", "gene_symbol", "gene_assignment", "GeneName")
    for (col in possible_gene_cols) {
      if (col %in% colnames(fData)) {
        gene_col <- col
        break
      }
    }
    
    if (!is.null(gene_col)) {
      cat("  Using gene annotation from column:", gene_col, "\n")
      gene_symbols <- fData[[gene_col]]
      
      # Clean gene symbols (Agilent often has format "GeneName // ..." or multiple)
      clean_symbols <- function(x) {
        x <- as.character(x)
        x <- gsub(" //.*", "", x)
        x <- gsub("---", NA, x)
        x[nchar(x) == 0] <- NA
        x
      }
      gene_symbols <- clean_symbols(gene_symbols)
      names(gene_symbols) <- rownames(exprs_mat)
      
      cat("  Mapped", sum(!is.na(gene_symbols)), "probes to gene symbols\n")
    } else {
      # Try using the platform annotation file
      cat("  No gene column found in feature data - checking platform annotation\n")
      
      gpl_file <- file.path(data_raw_dir, "GPL13607_old_annotations.txt.gz")
      if (file.exists(gpl_file)) {
        cat("  Reading platform annotation file:", gpl_file, "\n")
        # Custom parser for Agilent GPL format
        gpl_lines <- readLines(gzfile(gpl_file), warn = FALSE)
        # Find data start (skip header lines starting with ^, # or LINK)
        data_start <- which(!grepl("^[#^L]", gpl_lines))[1]
        # Column names from header comments
        col_names <- c("ID", "ProbeName", "GB_ACC", "ControlType", "accessions",
                       "GeneName", "Description", "chr_coord", "SEQUENCE", "SPOT_ID")
        gpl_annot <- read.delim(gzfile(gpl_file), skip = data_start - 1,
                                header = FALSE, stringsAsFactors = FALSE,
                                col.names = col_names,
                                quote = "\"", comment.char = "")
        cat("  GPL annotation rows:", nrow(gpl_annot), "\n")
        cat("  GPL annotation columns:", paste(colnames(gpl_annot), collapse = ", "), "\n")
        
        gene_symbols <- gpl_annot$GeneName
        names(gene_symbols) <- as.character(gpl_annot$ID)
        gene_symbols[nchar(gene_symbols) == 0 | gene_symbols == ""] <- NA
        
        # Map probes
        matched_symbols <- gene_symbols[rownames(exprs_mat)]
        names(matched_symbols) <- rownames(exprs_mat)
        gene_symbols <- matched_symbols
        cat("  Mapped", sum(!is.na(gene_symbols)), "probes to gene symbols\n")
      } else {
        stop("Cannot annotate GSE70947. No gene column in fData and no GPL annotation file.")
      }
    }
    
    # Collapse duplicates for Agilent
    exprs_mat_collapsed <- collapse_duplicate_genes(exprs_mat, gene_symbols, method = "max_mean")
    
  } else {
    # For Affymetrix datasets
    platform_pkg <- platform_map[[gse_id]]$pkg
    
    cat("\n  Annotating probes using", platform_pkg, "\n")
    annotation_result <- map_probes_to_genes(exprs_mat, platform_pkg)
    gene_symbols <- annotation_result$gene_symbols
    
    # Collapse duplicate genes
    exprs_mat_collapsed <- collapse_duplicate_genes(exprs_mat, gene_symbols, method = "max_mean")
  }
  
  # Store processed data
  processed_data[[gse_id]] <- list(
    exprs = exprs_mat_collapsed,
    pheno = pdata,
    group = group[valid_samples]
  )
  
  # Save individual processed files
  saveRDS(exprs_mat_collapsed, 
          file = file.path(data_processed_dir, paste0(gse_id, "_exprs.rds")))
  write.csv(exprs_mat_collapsed, 
            file = file.path(data_processed_dir, paste0(gse_id, "_exprs.csv")))
  
  saveRDS(pdata, 
          file = file.path(data_processed_dir, paste0(gse_id, "_pheno.rds")))
  write.csv(pdata, 
            file = file.path(data_processed_dir, paste0(gse_id, "_pheno.csv")))
  
  cat("  Processed", gse_id, "- Genes:", nrow(exprs_mat_collapsed), 
      ", Samples:", ncol(exprs_mat_collapsed), "\n")
  
  # Sample group summary
  sample_group_summary[[gse_id]] <- table(processed_data[[gse_id]]$group)
}

# --- Summary ---
cat("\n\n========================================================\n")
cat("PREPROCESSING COMPLETE - Sample Summary\n")
cat("========================================================\n")
for (gse_id in geo_ids) {
  cat(gse_id, ":\n")
  cat("  Platform:", platform_map[[gse_id]]$platform_id, "\n")
  cat("  Genes:", nrow(processed_data[[gse_id]]$exprs), "\n")
  cat("  Samples:", ncol(processed_data[[gse_id]]$exprs), "\n")
  cat("  Groups:\n")
  print(sample_group_summary[[gse_id]])
  cat("\n")
}

# Save summary
saveRDS(processed_data, file.path(data_processed_dir, "all_processed_data.rds"))

cat("\nAll individual datasets processed and saved.\n")
