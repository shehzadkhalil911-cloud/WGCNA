source("R/utils.R")
set.seed(123)

BASE_DIR <- "F:/independent projects/cancer/breast cancer"
OUT_DIR <- file.path(BASE_DIR, "Breast cancer 6 Ferroptosis")

log_msg("=== Step 1: Input Inventory ===")

inventory <- data.frame(
  file_path = character(),
  file_type = character(),
  dataset = character(),
  purpose = character(),
  n_rows = integer(),
  n_cols = integer(),
  selected = character(),
  reason = character(),
  stringsAsFactors = FALSE
)

add_entry <- function(path, ftype, dataset, purpose, nrow, ncol, sel, reason) {
  rbind(inventory, data.frame(
    file_path = path, file_type = ftype, dataset = dataset,
    purpose = purpose, n_rows = nrow, n_cols = ncol,
    selected = sel, reason = reason, stringsAsFactors = FALSE
  ))
}

# --- Expression CSVs ---
for (ds in c("GSE42568", "GSE15852", "GSE29431", "GSE70947")) {
  csv_path <- file.path(BASE_DIR, "Datasets/data/processed", paste0(ds, "_expression.csv"))
  if (file.exists(csv_path)) {
    dat <- read.csv(csv_path, row.names = 1, check.names = FALSE)
    is_discovery <- ds %in% c("GSE42568", "GSE15852", "GSE29431")
    inventory <- add_entry(csv_path, "expression CSV (gene x sample)", ds,
                           ifelse(is_discovery, "Discovery cohort expression", "Validation cohort expression"),
                           nrow(dat), ncol(dat),
                           ifelse(is_discovery, "Discovery", "Validation"),
                           ifelse(is_discovery,
                                  sprintf("Discovery cohort: %d genes, %d samples", nrow(dat), ncol(dat)),
                                  sprintf("Independent validation: %d genes, %d samples", nrow(dat), ncol(dat))))
    log_msg(sprintf("%s: %d genes x %d samples", ds, nrow(dat), ncol(dat)))
  }
}

# --- Processed RDS ---
for (ds in c("GSE42568", "GSE15852", "GSE29431", "GSE70947")) {
  rds_path <- file.path(BASE_DIR, "Datasets/data/processed", paste0(ds, "_processed.rds"))
  if (file.exists(rds_path)) {
    inventory <- add_entry(rds_path, "RDS (processed ExpressionSet)", ds,
                           "Processed ExpressionSet object", NA, NA, "Yes",
                           "Pre-processed ExpressionSet for expression data")
  }
}

# --- Metadata CSVs ---
for (ds in c("GSE42568", "GSE15852", "GSE29431", "GSE70947")) {
  meta_path <- file.path(BASE_DIR, "Datasets/data/metadata", paste0(ds, "_sample_classes.csv"))
  if (file.exists(meta_path)) {
    meta <- read.csv(meta_path)
    tbl <- table(meta$class)
    n_tumor <- sum(grepl("tumor|Tumor", meta$class))
    n_normal <- sum(grepl("normal|Normal", meta$class))
    inventory <- add_entry(meta_path, "CSV (sample metadata)", ds,
                           "Sample class labels", nrow(meta), ncol(meta), "Yes",
                           sprintf("%d Tumor, %d Normal samples", n_tumor, n_normal))
    log_msg(sprintf("%s metadata: %d Tumor, %d Normal", ds, n_tumor, n_normal))
  }
}

# --- Raw eset RDS ---
for (ds in c("GSE42568", "GSE15852", "GSE29431", "GSE70947")) {
  eset_path <- file.path(BASE_DIR, "Datasets/data/raw", paste0(ds, "_eset.rds"))
  if (file.exists(eset_path)) {
    inventory <- add_entry(eset_path, "RDS (raw ExpressionSet)", ds,
                           "Raw GEO ExpressionSet", NA, NA, "Available as fallback",
                           "Raw expression set from GEO")
  }
}

# --- Platform files ---
for (plat in c("GPL96", "GPL570", "GPL13607")) {
  plat_path <- file.path(BASE_DIR, "Datasets/data/raw", paste0(plat, ".soft.gz"))
  if (file.exists(plat_path)) {
    inventory <- add_entry(plat_path, "SOFT.gz (platform annotation)", plat,
                           "Platform annotation", NA, NA, "Yes",
                           sprintf("Platform annotation file for %s", plat))
  }
}

# --- Discovery batch-corrected objects ---
disc_path <- file.path(BASE_DIR, "Datasets/data/discovery/discovery_cohort.rds")
if (file.exists(disc_path)) {
  inventory <- add_entry(disc_path, "RDS (merged discovery)", "Discovery",
                         "Merged discovery cohort expression", NA, NA, "Reference",
                         "Pre-merged discovery cohort expression matrix")
}

combat_path <- file.path(BASE_DIR, "Datasets/data/discovery/discovery_combat.rds")
if (file.exists(combat_path)) {
  inventory <- add_entry(combat_path, "RDS (combat-corrected)", "Discovery",
                         "Batch-corrected discovery cohort", NA, NA, "Reference",
                         "Previously ComBat-corrected discovery data")
}

val_path <- file.path(BASE_DIR, "Datasets/data/external_validation/validation_cohort.rds")
if (file.exists(val_path)) {
  inventory <- add_entry(val_path, "RDS (validation cohort)", "Validation",
                         "Validation cohort expression", NA, NA, "Reference",
                         "Validation cohort expression matrix")
}

# --- Ferroptosis gene list ---
ferr_path <- file.path(BASE_DIR, "Ferroptosis Genes.xlsx")
if (file.exists(ferr_path)) {
  ferr_sheets <- openxlsx::getSheetNames(ferr_path)
  ferr_raw <- openxlsx::read.xlsx(ferr_path, sheet = 1)
  inventory <- add_entry(ferr_path, "XLSX (ferroptosis gene list)", "Ferroptosis",
                         "Curated ferroptosis gene list", nrow(ferr_raw), ncol(ferr_raw), "Yes",
                         sprintf("%d rows, %d cols, sheets: %s", nrow(ferr_raw), ncol(ferr_raw),
                                 paste(ferr_sheets, collapse = ", ")))
  log_msg(sprintf("Ferroptosis gene list: %d rows, %d columns", nrow(ferr_raw), ncol(ferr_raw)))
  log_msg(sprintf("Column names: %s", paste(colnames(ferr_raw), collapse = ", ")))
}

# --- Raw list ---
raw_list_path <- file.path(BASE_DIR, "Datasets/data/raw/geo_raw_list.rds")
if (file.exists(raw_list_path)) {
  inventory <- add_entry(raw_list_path, "RDS (raw GEO list)", "All",
                         "List of raw GEO objects", NA, NA, "Available as fallback",
                         "Combined list of raw GEO expression sets")
}

write.csv(inventory, file.path(OUT_DIR, "results/input_inventory/input_manifest.csv"), row.names = FALSE)
writexl::write_xlsx(inventory, file.path(OUT_DIR, "results/input_inventory/input_manifest.xlsx"))

log_msg(sprintf("Input manifest saved: %d entries", nrow(inventory)))
log_msg("Step 1 complete: Input inventory created")
