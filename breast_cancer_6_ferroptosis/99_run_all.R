# Master script to run the complete breast cancer-ferroptosis crosstalk analysis
# Set base directory
BASE_DIR <- "F:/independent projects/cancer/breast cancer"
OUT_DIR <- file.path(BASE_DIR, "Breast cancer 6 Ferroptosis")
set.seed(123)

# Initialize logs
for (logf in c("logs/analysis.log", "logs/error.log")) {
  write(sprintf("=== Breast Cancer-Ferroptosis Crosstalk Analysis ===\nStarted: %s",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        file = file.path(OUT_DIR, logf))
}

start_time <- Sys.time()

source("R/utils.R")

log_msg("Starting full analysis pipeline", file.path(OUT_DIR, "logs/analysis.log"))

scripts <- c(
  "00_check_environment.R",
  "01_input_inventory.R",
  "02_ferroptosis_gene_list.R",
  "03_load_and_preprocess.R",
  "04_batch_correction.R",
  "05_differential_expression.R",
  "06_ferroptosis_DE.R",
  "07_ferroptosis_score.R",
  "08_WGCNA.R",
  "09_hub_genes.R",
  "10_enrichment.R",
  "11_pathway_crosstalk.R",
  "12_PPI_network.R",
  "13_immune_infiltration.R",
  "14_validation.R",
  "15_module_preservation.R"
)

for (script in scripts) {
  script_path <- file.path(OUT_DIR, "scripts", script)
  if (file.exists(script_path)) {
    log_msg(sprintf("Running script: %s", script), file.path(OUT_DIR, "logs/analysis.log"))
    tryCatch({
      source(script_path)
      log_msg(sprintf("Completed: %s", script), file.path(OUT_DIR, "logs/analysis.log"))
    }, error = function(e) {
      log_error(sprintf("Script %s failed: %s", script, e$message),
                file.path(OUT_DIR, "logs/error.log"))
    })
  } else {
    log_msg(sprintf("Script not found: %s", script), file.path(OUT_DIR, "logs/analysis.log"))
  }
}

elapsed <- difftime(Sys.time(), start_time, units = "mins")
log_msg(sprintf("Pipeline completed in %.2f minutes", elapsed),
        file.path(OUT_DIR, "logs/analysis.log"))

# Save session info
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo.txt"))

log_msg("Analysis pipeline completed")
