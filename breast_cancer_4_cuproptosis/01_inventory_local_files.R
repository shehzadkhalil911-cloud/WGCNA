# 01_inventory_local_files.R
source("config/analysis_config.R")
source("R/load_functions.R")
ensure_paths()
suppressMessages({ library(openxlsx); library(dplyr); library(officer) })

inv <- inventory_local_files(DATASETS_DIR, root_dir = ".")
crg_files <- detect_cuproptosis_lists(inv)

## selection logic
inv$selected_for_analysis <- FALSE
inv$reason <- ""
## select processed rds, discovery/validation merged, metadata csvs, soft files, cuproptosis list
inv$selected_for_analysis <- grepl("processed.*\\.rds$|discovery_cohort\\.rds$|discovery_combat\\.rds$|validation_cohort\\.rds$|sample_classes\\.csv$|GPL[0-9]+\\.soft|eset\\.rds$|Cuproptosis genes\\.csv$", inv$file_path)
inv$reason[inv$selected_for_analysis] <- "usable preprocessed/annotation/cupoprtosis file"
inv$reason[!inv$selected_for_analysis & inv$readable] <- "not required for pipeline or raw auxiliary"

## save manifest CSV
write.csv(inv, file.path(paths$res_inv, "input_file_manifest.csv"), row.names = FALSE)

## save Word manifest
wb <- createWorkbook()
addWorksheet(wb, "input_manifest")
writeData(wb, "input_manifest", inv)
addWorksheet(wb, "cuproptosis_lists")
writeData(wb, "cuproptosis_lists", crg_files)
saveWorkbook(wb, file.path(paths$res_inv, "input_file_manifest.xlsx"), overwrite = TRUE)

## if multiple cuproptosis lists, compare (here only one expected)
if (nrow(crg_files) > 1) {
  cat("WARNING: multiple cuproptosis gene-list candidates detected. Review report.\n")
}
saveRDS(inv, file.path(paths$res_inv, "input_file_manifest.rds"))
cat("Inventory complete:", nrow(inv), "files,", sum(inv$selected_for_analysis), "selected.\n")
cat("Cuproptosis list file(s):\n"); print(crg_files$file_path)
quit(save = "no", status = 0, runLast = FALSE)
