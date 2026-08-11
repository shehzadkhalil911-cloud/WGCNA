# 04_prepare_cuproptosis_gene_list.R
source("config/analysis_config.R")
source("R/load_functions.R")
ensure_paths()
suppressMessages({ library(dplyr); library(openxlsx); library(officer) })

## locate the cuproptosis gene list (project root, one level up from output dir)
crg_file <- "../Cuproptosis genes.csv"
if (!file.exists(crg_file)) {
  ## fallback search in Datasets
  cand <- list.files(DATASETS_DIR, pattern = "Cuproptosis", recursive = TRUE, full.names = TRUE)
  crg_file <- cand[1]
}
raw <- read.csv(crg_file, stringsAsFactors = FALSE)[, 1]
names(raw) <- NULL
std <- standardize_cuproptosis_list(raw)

## save all components
write.csv(data.frame(original_symbol = std$original), file.path(paths$crg, "original_gene_list.csv"), row.names = FALSE)
write.csv(std$standardized, file.path(paths$crg, "standardized_cuproptosis_genes.csv"), row.names = FALSE)
write.csv(data.frame(duplicated_entry = std$duplicates), file.path(paths$crg, "duplicated_entries.csv"), row.names = FALSE)
write.csv(data.frame(alias_original = std$alias_corrected$original, official_symbol = std$alias_corrected$official_symbol),
          file.path(paths$crg, "alias_corrected.csv"), row.names = FALSE)
write.csv(data.frame(unmapped_symbol = std$unmapped), file.path(paths$crg, "unmapped_genes.csv"), row.names = FALSE)
saveRDS(std, file.path(paths$crg, "cuproptosis_standardized.rds"))

cat("Cuproptosis list prepared.\n")
cat(" Original entries:", length(std$original), "\n")
cat(" Standardized (final) CRGs:", nrow(std$standardized), "\n")
cat(" Duplicates:", length(std$duplicates), "\n")
cat(" Alias-corrected:", nrow(std$alias_corrected), "\n")
cat(" Unmapped:", length(std$unmapped), "\n")
if (length(std$unmapped) > 0) cat("  unmapped:", paste(std$unmapped, collapse = ", "), "\n")
quit(save = "no", status = 0, runLast = FALSE)
