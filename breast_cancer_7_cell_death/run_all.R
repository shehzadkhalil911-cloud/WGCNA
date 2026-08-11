#!/usr/bin/env Rscript
# =============================================================================
# Master Runner: Breast Cancer 7 Cell Death Crosstalk
# Run all analysis scripts sequentially
# =============================================================================

# Set working directory to scripts folder
scripts_dir <- normalizePath(file.path(
  "F:", "independent projects", "cancer", "breast cancer",
  "Breast cancer 7 Cell Death Crosstalk", "scripts"
))

setwd(scripts_dir)

# Timestamp
cat("=====================================\n")
cat("Starting analysis pipeline\n")
cat(sprintf("Time: %s\n", Sys.time()))
cat("=====================================\n\n")

# Source config
source(file.path("..", "config", "analysis_config.R"))

# Log start
log_message("=== ANALYSIS PIPELINE START ===")

# Record runtime
start_time <- Sys.time()

# Script list
scripts <- c(
  "01_input_inventory_and_gene_lists.R",
  "02_metadata_and_preprocessing.R",
  "03_batch_correction_and_DE.R",
  "04_regulated_cell_death_gene_analysis.R",
  "05_activity_scores.R",
  "06_WGCNA.R",
  "07_enrichment.R",
  "08_pathway_crosstalk.R",
  "09_PPI_and_immune.R",
  "10_multi_evidence_prioritisation.R",
  "11_render_report.R"
)

# Run each script
for (script in scripts) {
  cat(sprintf("\n>>> Running: %s\n", script))
  cat(sprintf("    Start: %s\n", Sys.time()))
  
  script_start <- Sys.time()
  
  tryCatch({
    source(script, local = TRUE, echo = FALSE)
    cat(sprintf("    Completed: %s\n", Sys.time()))
    cat(sprintf("    Duration: %.1f mins\n", 
                difftime(Sys.time(), script_start, units = "mins")))
  }, error = function(e) {
    cat(sprintf("    ERROR in %s: %s\n", script, e$message))
    log_message(sprintf("Pipeline error in %s: %s", script, e$message), "ERROR")
  })
  
  cat("-------------------------------------\n")
}

# Total runtime
total_time <- difftime(Sys.time(), start_time, units = "mins")
cat(sprintf("\nTotal pipeline runtime: %.1f minutes\n", total_time))
log_message(sprintf("=== ANALYSIS PIPELINE COMPLETE: %.1f minutes ===\n", total_time))

# Print summary
cat("\n=====================================\n")
cat("ANALYSIS COMPLETE\n")
cat(sprintf("Total time: %.1f mins\n", total_time))
cat("=====================================\n")
