# ============================================================
# Master Runner Script: Breast Cancer 5 Necroptosis Analysis
# ============================================================
set.seed(123)
options(warn=1)
options(stringsAsFactors=FALSE)

start_time <- Sys.time()
cat("========================================\n")
cat("BREAST CANCER NECROPTOSIS ANALYSIS\n")
cat("Start time:", format(start_time), "\n")
cat("R version:", R.version.string, "\n")
cat("Working directory:", getwd(), "\n")
cat("========================================\n\n")

base_dir <- normalizePath(file.path(".", "Breast Cancer 5 Necroptosis"))
scripts_dir <- file.path(base_dir, "scripts")
log_dir <- file.path(base_dir, "logs")

dir.create(log_dir, showWarnings=FALSE, recursive=TRUE)

setup_log <- file.path(log_dir, "setup.log")
sink(setup_log, split=TRUE)
cat("Running setup...\n")
source(file.path(scripts_dir, "00_setup.R"), echo=TRUE)
sink()

scripts <- c(
  "01_data_inspection.R",
  "02_preprocessing.R",
  "03_batch_correction.R",
  "04_differential_expression.R",
  "05_necroptosis_analysis.R",
  "06_necroptosis_score.R",
  "07_enrichment.R",
  "08_networks.R",
  "09_machine_learning.R",
  "10_consensus_genes.R",
  "11_predictive_model.R",
  "12_immune_infiltration.R",
  "13_immune_crosstalk.R",
  "14_external_validation.R"
)

for(script in scripts) {
  cat("\n\n========================================\n")
  cat("Running:", script, "\n")
  cat("Time:", format(Sys.time()), "\n")
  cat("========================================\n\n")
  
  script_path <- file.path(scripts_dir, script)
  log_path <- file.path(log_dir, gsub("\\.R$", ".log", script))
  
  tryCatch({
    sink(log_path, split=TRUE)
    source(script_path, echo=TRUE, local=TRUE)
    sink()
    cat("[SUCCESS]", script, "completed.\n")
  }, error=function(e) {
    sink()
    cat("[ERROR]", script, "failed:", e$message, "\n")
    error_log <- file.path(log_dir, "errors.log")
    cat(format(Sys.time()), script, ":", e$message, "\n", file=error_log, append=TRUE)
  })
}

# Save session info
session_info <- sessionInfo()
saveRDS(session_info, file.path(base_dir, "sessionInfo.rds"))
capture.output(session_info, file=file.path(base_dir, "sessionInfo.txt"))

end_time <- Sys.time()
cat("\n\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("End time:", format(end_time), "\n")
cat("Total runtime:", format(difftime(end_time, start_time)), "\n")
cat("========================================\n")
