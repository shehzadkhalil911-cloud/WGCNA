suppressMessages({
  library(Biobase)
  library(data.table)
})

out_dir   <- "F:/independent projects/cancer/breast cancer/Breast cancer 8"
data_root <- "F:/independent projects/cancer/breast cancer/Datasets"

manifest <- data.frame(
  file          = character(0),
  path          = character(0),
  size_bytes    = numeric(0),
  dataset       = character(0),
  role          = character(0),
  format        = character(0),
  stringsAsFactors = FALSE
)

add_file <- function(path, dataset, role) {
  if (!file.exists(path)) return()
  manifest <<- rbind(manifest, data.frame(
    file       = basename(path),
    path       = path,
    size_bytes = file.info(path)$size,
    dataset    = dataset,
    role       = role,
    format     = tools::file_ext(path),
    stringsAsFactors = FALSE
  ))
}

add_file(file.path(data_root, "data/metadata/dataset_inspection.rds"), "All", "metadata summary")
add_file(file.path(data_root, "data/metadata/GSE70947_sample_classes.csv"), "GSE70947", "sample classes")
add_file(file.path(data_root, "data/metadata/GSE15852_sample_classes.csv"), "GSE15852", "sample classes")
add_file(file.path(data_root, "data/raw/GSE70947_eset.rds"), "GSE70947", "raw eset")
add_file(file.path(data_root, "data/raw/GSE15852_eset.rds"), "GSE15852", "raw eset")
add_file(file.path(data_root, "data/raw/GPL13607.soft.gz"), "GSE70947", "platform annotation")
add_file(file.path(data_root, "data/raw/GPL96.soft.gz"), "GSE15852", "platform annotation")
add_file(file.path(data_root, "data/raw/GPL570.soft.gz"), "Other (GPL570)", "platform annotation")
add_file(file.path(data_root, "data/processed/GSE70947_expression.csv"), "GSE70947", "processed expression csv")
add_file(file.path(data_root, "data/processed/GSE70947_processed.rds"), "GSE70947", "processed rds")
add_file(file.path(data_root, "data/processed/GSE15852_expression.csv"), "GSE15852", "processed expression csv")
add_file(file.path(data_root, "data/processed/GSE15852_processed.rds"), "GSE15852", "processed rds")
add_file(file.path(data_root, "data/discovery/discovery_cohort.rds"), "GSE70947", "discovery cohort rds")
add_file(file.path(data_root, "data/discovery/discovery_combat.rds"), "GSE70947", "combat batch corrected")
add_file(file.path(data_root, "data/external_validation/validation_cohort.rds"), "GSE15852", "validation cohort rds")

write.csv(manifest, file.path(out_dir, "tables/Table_01_Input_File_Manifest.csv"), row.names = FALSE)

cat("MANIFEST ROWS:", nrow(manifest), "\n")

d_info <- readRDS(file.path(data_root, "data/metadata/dataset_inspection.rds"))

cat("\n===== DATASET INSPECTION SUMMARY =====\n")
for (nm in names(d_info)) {
  x <- d_info[[nm]]
  cat(sprintf("\n%s: platform=%s, n_samples=%d, n_probes=%d\n", nm, x$platform, x$n_samples, x$n_probes))
  cat(sprintf("  expression range: %.2f .. %.2f, median %.2f, has_negative=%s, needs_log2=%s\n",
      x$expression_range[1], x$expression_range[2], x$expression_median,
      x$has_negative, x$needs_log2))
}

cat("\n===== GSE70947 (DISCOVERY) RAW ESET =====\n")
es <- readRDS(file.path(data_root, "data/raw/GSE70947_eset.rds"))
m  <- exprs(es)
pd <- pData(es)
cat("exprs dim:", nrow(m), "x", ncol(m), "\n")
cat("exprs rownames 1:10:", head(rownames(m), 10), "\n")
cat("exprs colnames 1:6:", head(colnames(m), 6), "\n")
cat("exprs NA count:", sum(is.na(m)), "\n")
cat("exprs range:", range(m, na.rm = TRUE), "\n")

tab_tissue <- table(pd$`tissue:ch1`)
cat("tissue table:\n"); print(tab_tissue)
tab_src <- table(pd$source_name_ch1)
cat("source_name table:\n"); print(tab_src)

title_tab <- table(grepl("normal$", pd$title), useNA = "ifany")
cat("title ends in 'normal':\n"); print(title_tab)

cat("er:", table(pd$`er:ch1`), "\n")
cat("her2:", table(pd$`her2:ch1`), "\n")
cat("cohort:", table(pd$`cohort:ch1`), "\n")

cat("\n===== GSE15852 (VALIDATION) RAW ESET =====\n")
es2 <- readRDS(file.path(data_root, "data/raw/GSE15852_eset.rds"))
m2  <- exprs(es2)
pd2 <- pData(es2)
cat("exprs dim:", nrow(m2), "x", ncol(m2), "\n")
cat("exprs rownames 1:10:", head(rownames(m2), 10), "\n")
cat("exprs NA count:", sum(is.na(m2)), "\n")
cat("exprs range:", range(m2, na.rm = TRUE), "\n")

cat("title table:\n"); print(table(pd2$title[1:10]))
cat("source_name table:\n"); print(table(pd2$source_name_ch1))
cat("histopath table:\n"); print(table(pd2$`histopathological exam:ch1`))
cat("grade table:\n"); print(table(pd2$`grade:ch1`))
