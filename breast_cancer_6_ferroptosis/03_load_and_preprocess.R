source("R/utils.R")
set.seed(123)

BASE_DIR <- "F:/independent projects/cancer/breast cancer"
OUT_DIR <- file.path(BASE_DIR, "Breast cancer 6 Ferroptosis")

log_msg("=== Step 3: Load and Preprocess Datasets ===")

# Load ferroptosis genes
ferr_genes <- readRDS(file.path(OUT_DIR, "data_processed/ferroptosis_gene_list/ferroptosis_genes_final.rds"))
log_msg(sprintf("Loaded %d ferroptosis genes", length(ferr_genes)))

# Helper to load and process each dataset
preprocess_dataset <- function(ds_name, is_validation = FALSE) {
  csv_path <- file.path(BASE_DIR, "Datasets/data/processed", paste0(ds_name, "_expression.csv"))
  meta_path <- file.path(BASE_DIR, "Datasets/data/metadata", paste0(ds_name, "_sample_classes.csv"))

  if (!file.exists(csv_path)) stop(sprintf("Expression file not found: %s", csv_path))
  if (!file.exists(meta_path)) stop(sprintf("Metadata file not found: %s", meta_path))

  expr <- as.matrix(read.csv(csv_path, row.names = 1, check.names = FALSE))
  meta <- read.csv(meta_path)

  class_col <- intersect(c("class", "Group", "group", "Condition", "condition"), colnames(meta))[1]
  sample_col <- intersect(c("sample_id", "Sample", "sample", "ID", "id"), colnames(meta))[1]

  if (is.na(class_col) || is.na(sample_col)) {
    log_msg(sprintf("WARNING: %s - Using default column names", ds_name))
    class_col <- colnames(meta)[grep("class|group|condition", tolower(colnames(meta)))][1]
    sample_col <- colnames(meta)[grep("sample|id", tolower(colnames(meta)))][1]
  }

  # Check sample alignment
  common_samples <- intersect(colnames(expr), meta[[sample_col]])
  if (length(common_samples) == 0) {
    # Try matching after removing prefix
    short_names <- gsub("^X", "", colnames(expr))
    common_samples <- intersect(short_names, meta[[sample_col]])
    if (length(common_samples) > 0) {
      colnames(expr) <- short_names
    }
  }

  meta <- meta[meta[[sample_col]] %in% common_samples, ]
  expr <- expr[, meta[[sample_col]], drop = FALSE]

  # Classify samples
  group <- ifelse(grepl("tumor|Tumor|Cancer|cancer", meta[[class_col]]), "Tumor", "Normal")

  # Detect pairing
  is_paired <- FALSE
  pair_info <- NULL
  if (ds_name == "GSE15852") {
    # Alternating Normal/Tumor pattern suggests paired samples
    is_paired <- TRUE
    pair_info <- data.frame(
      sample_id = meta[[sample_col]],
      class = group,
      pair_id = rep(1:(nrow(meta)/2), each = 2),
      stringsAsFactors = FALSE
    )
  }
  if (ds_name == "GSE70947") {
    is_paired <- TRUE
    # First 148 are Normal, next 148 are Tumor
    n <- nrow(meta)
    pair_info <- data.frame(
      sample_id = meta[[sample_col]],
      class = group,
      pair_id = c(1:(n/2), 1:(n/2)),
      stringsAsFactors = FALSE
    )
  }

  # Check expression scale (log2 or not)
  expr_range <- range(expr, na.rm = TRUE)
  needs_log2 <- if (max(expr_range, na.rm = TRUE) > 100) TRUE else FALSE

  # Check for missing values
  na_count <- sum(is.na(expr))
  na_pct <- na_count / prod(dim(expr)) * 100

  log_msg(sprintf("%s: %d genes x %d samples, range [%.2f, %.2f], NAs: %.2f%%",
                  ds_name, nrow(expr), ncol(expr), expr_range[1], expr_range[2], na_pct))

  # Remove genes with excessive missingness (>50%)
  max_na <- apply(expr, 1, function(x) mean(is.na(x)))
  genes_to_keep <- max_na < 0.5
  expr <- expr[genes_to_keep, ]
  log_msg(sprintf("%s: %d/%d genes retained after NA filter", ds_name, sum(genes_to_keep), length(genes_to_keep)))

  list(
    expr = expr,
    metadata = meta,
    group = group,
    sample_col = sample_col,
    class_col = class_col,
    is_paired = is_paired,
    pair_info = pair_info,
    needs_log2 = needs_log2,
    na_pct = na_pct
  )
}

# Process discovery datasets
disc_names <- c("GSE42568", "GSE15852", "GSE29431")
discovery_list <- list()

for (ds in disc_names) {
  log_msg(sprintf("Processing %s...", ds))
  discovery_list[[ds]] <- preprocess_dataset(ds, is_validation = FALSE)
}

# Process validation independently
log_msg("Processing GSE70947 (validation)...")
val_data <- preprocess_dataset("GSE70947", is_validation = TRUE)

# Save individual processed data
for (ds in disc_names) {
  d <- discovery_list[[ds]]
  out_path <- file.path(OUT_DIR, "data_processed/individual_datasets")
  write.csv(d$expr, file.path(out_path, paste0(ds, "_expr_processed.csv")))
  saveRDS(d, file.path(out_path, paste0(ds, "_processed.rds")))
  # Metadata
  meta_df <- data.frame(
    sample_id = colnames(d$expr),
    group = d$group,
    stringsAsFactors = FALSE
  )
  write.csv(meta_df, file.path(out_path, paste0(ds, "_metadata.csv")), row.names = FALSE)
}

# Save validation
val_out <- file.path(OUT_DIR, "data_processed/validation")
write.csv(val_data$expr, file.path(val_out, "GSE70947_expr_processed.csv"))
saveRDS(val_data, file.path(val_out, "GSE70947_processed.rds"))
meta_val <- data.frame(
  sample_id = colnames(val_data$expr),
  group = val_data$group,
  stringsAsFactors = FALSE
)
write.csv(meta_val, file.path(val_out, "GSE70947_metadata.csv"), row.names = FALSE)

# Create sample composition summary
sample_summary <- do.call(rbind, lapply(names(discovery_list), function(ds) {
  d <- discovery_list[[ds]]
  tbl <- table(d$group)
  data.frame(
    Dataset = ds,
    Tumor = tbl["Tumor"] %||% 0,
    Normal = tbl["Normal"] %||% 0,
    Total = sum(tbl),
    Paired = d$is_paired,
    stringsAsFactors = FALSE
  )
}))
val_tbl <- table(val_data$group)
sample_summary <- rbind(sample_summary, data.frame(
  Dataset = "GSE70947",
  Tumor = val_tbl["Tumor"] %||% 0,
  Normal = val_tbl["Normal"] %||% 0,
  Total = sum(val_tbl),
  Paired = val_data$is_paired,
  stringsAsFactors = FALSE
))

write.csv(sample_summary, file.path(OUT_DIR, "results/input_inventory/sample_composition.csv"), row.names = FALSE)
log_msg(sprintf("Sample composition:\n%s", paste(capture.output(print(sample_summary)), collapse = "\n")))

saveRDS(discovery_list, file.path(OUT_DIR, "data_processed/discovery/discovery_list.rds"))
saveRDS(val_data, file.path(OUT_DIR, "data_processed/validation/validation_data.rds"))

log_msg("Step 3 complete: All datasets loaded and preprocessed")
