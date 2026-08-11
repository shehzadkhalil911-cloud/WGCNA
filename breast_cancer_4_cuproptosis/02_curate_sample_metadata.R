# 02_curate_sample_metadata.R
source("config/analysis_config.R")
source("R/load_functions.R")
ensure_paths()
suppressMessages({ library(dplyr); library(openxlsx); library(officer) })

LOCAL_PROC <- file.path(LOCAL_DATA_DIR, "processed")
META_CSV  <- file.path(LOCAL_DATA_DIR, "metadata")

proc_paths <- setNames(
  file.path(LOCAL_PROC, paste0(DISCOVERY_GSE, "_processed.rds")), DISCOVERY_GSE)
proc_paths[[VALIDATION_GSE]] <- file.path(LOCAL_PROC, paste0(VALIDATION_GSE, "_processed.rds"))

meta_tables <- list()
for (g in c(DISCOVERY_GSE, VALIDATION_GSE)) {
  obj <- readRDS(proc_paths[[g]])
  sc_csv <- file.path(META_CSV, paste0(g, "_sample_classes.csv"))
  meta_tables[[g]] <- build_metadata_table(obj, g, sc_csv)
}

all_meta <- compile_metadata(meta_tables)
## confidence should be high for all (curated csv present)
low_conf <- all_meta[all_meta$confidence == "low", ]
if (nrow(low_conf) > 0) {
  cat("AMBIGUOUS GROUP ASSIGNMENT for", nrow(low_conf), "samples. STOPPING for manual review.\n")
  print(low_conf)
  stop("Manual verification required before proceeding.")
}

write.csv(all_meta, file.path(paths$metadata, "sample_group_assignment.csv"), row.names = FALSE)
wb <- createWorkbook(); addWorksheet(wb, "sample_metadata"); writeData(wb, "sample_metadata", all_meta)
saveWorkbook(wb, file.path(paths$metadata, "sample_group_assignment.xlsx"), overwrite = TRUE)
saveRDS(all_meta, file.path(paths$metadata, "sample_group_assignment.rds"))

## dataset composition
comp <- all_meta %>% group_by(dataset, assigned_group) %>% summarise(n = n(), .groups = "drop")
write.csv(comp, file.path(paths$metadata, "dataset_composition.csv"), row.names = FALSE)
cat("Metadata curation complete.\n")
print(comp)
quit(save = "no", status = 0, runLast = FALSE)
