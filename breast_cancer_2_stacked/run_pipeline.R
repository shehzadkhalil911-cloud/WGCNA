setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# run_pipeline.R
# Master pipeline script - runs all analysis scripts sequentially
# =============================================================================

cat("========================================================\n")
cat("BREAST CANCER STACKED ENSEMBLE PIPELINE\n")
cat("========================================================\n\n")

# Record start time
start_time <- Sys.time()
cat("Start time:", format(start_time), "\n\n")

# Step 0: Install packages (if needed)
cat("\n--- Step 00: Install Packages ---\n")
source("scripts/00_install_packages.R")

# Step 1: Download GEO datasets
cat("\n--- Step 01: Download GEO Datasets ---\n")
tryCatch(source("scripts/01_download_GEO.R"), error = function(e) {
    cat("ERROR in step 01:", e$message, "\n")
    # Check for cached data
    if (!file.exists("data/raw/geo_raw_list.rds")) {
        stop("GEO download failed and no cached data available.")
    }
    cat("Using cached GEO data.\n")
})

# Step 2: Inspect metadata and cohort suitability
cat("\n--- Step 02: Inspect Metadata ---\n")
tryCatch(source("scripts/02_inspect_metadata.R"), error = function(e) {
    cat("ERROR in step 02:", e$message, "\n")
})

# Step 3: Preprocess, normalize, annotate
cat("\n--- Step 03: Preprocess, Normalize, Annotate ---\n")
tryCatch(source("scripts/03_preprocess_normalize_annotate.R"), error = function(e) {
    cat("ERROR in step 03:", e$message, "\n")
    stop("Preprocessing failed.")
})

# Step 4: Construct discovery cohort and check suitability
cat("\n--- Step 04: Construct Discovery Cohort ---\n")
tryCatch(source("scripts/04_construct_discovery_cohort.R"), error = function(e) {
    cat("ERROR in step 04:", e$message, "\n")
    stop("Cohort construction failed.")
})

# Read suitability assessment
suitability <- read.csv("results/cohort_suitability/cohort_suitability.csv")
val_info <- suitability[suitability$Dataset == VALIDATION_GSE, ]

if (!val_info$Has_Both_Classes || val_info$Common_Genes < 100) {
    cat("\n*** SUITABILITY WARNING ***\n")
    cat("Selected validation cohort", VALIDATION_GSE, "may not be suitable.\n")
    cat("Has both classes:", val_info$Has_Both_Classes, "\n")
    cat("Common genes:", val_info$Common_Genes, "\n")
    cat("Consider changing VALIDATION_GSE in config/analysis_config.R\n")
    # Ask to continue
    cat("Continuing with pipeline...\n")
}

# Step 5: Batch correction
cat("\n--- Step 05: Batch Correction & PCA ---\n")
tryCatch(source("scripts/05_batch_correction_PCA.R"), error = function(e) {
    cat("ERROR in step 05:", e$message, "\n")
    stop("Batch correction failed.")
})

# Step 6: Differential expression
cat("\n--- Step 06: Differential Expression ---\n")
tryCatch(source("scripts/06_differential_expression.R"), error = function(e) {
    cat("ERROR in step 06:", e$message, "\n")
    stop("DE analysis failed.")
})

# Step 7: Data splits
cat("\n--- Step 07: Generate Data Splits ---\n")
tryCatch(source("scripts/07_generate_data_splits.R"), error = function(e) {
    cat("ERROR in step 07:", e$message, "\n")
    stop("Data split generation failed.")
})

# Steps 8-12: Feature selection for each algorithm
cat("\n--- Step 08: RF Feature Selection ---\n")
tryCatch(source("scripts/08_RF_feature_selection.R"), error = function(e) {
    cat("ERROR in step 08:", e$message, "\n")
})

cat("\n--- Step 09: GB Feature Selection ---\n")
tryCatch(source("scripts/09_GB_feature_selection.R"), error = function(e) {
    cat("ERROR in step 09:", e$message, "\n")
})

cat("\n--- Step 10: SVM Feature Selection ---\n")
tryCatch(source("scripts/10_SVM_feature_selection.R"), error = function(e) {
    cat("ERROR in step 10:", e$message, "\n")
})

cat("\n--- Step 11: KNN Feature Selection ---\n")
tryCatch(source("scripts/11_KNN_feature_selection.R"), error = function(e) {
    cat("ERROR in step 11:", e$message, "\n")
})

cat("\n--- Step 12: NN Feature Selection ---\n")
tryCatch(source("scripts/12_NN_feature_selection.R"), error = function(e) {
    cat("ERROR in step 12:", e$message, "\n")
})

# Step 13: Consensus genes
cat("\n--- Step 13: Consensus Feature Genes ---\n")
tryCatch(source("scripts/13_consensus_feature_genes.R"), error = function(e) {
    cat("ERROR in step 13:", e$message, "\n")
    stop("Consensus gene analysis failed.")
})

# Step 14: Train base classifiers
cat("\n--- Step 14: Train Base Classifiers ---\n")
tryCatch(source("scripts/14_train_base_classifiers.R"), error = function(e) {
    cat("ERROR in step 14:", e$message, "\n")
    stop("Base classifier training failed.")
})

# Step 15: Generate OOF predictions
cat("\n--- Step 15: Generate OOF Predictions ---\n")
tryCatch(source("scripts/15_generate_OOF_predictions.R"), error = function(e) {
    cat("ERROR in step 15:", e$message, "\n")
    stop("OOF generation failed.")
})

# Step 16: Train stacked Elastic Net
cat("\n--- Step 16: Train Stacked Elastic Net ---\n")
tryCatch(source("scripts/16_train_stacked_elastic_net.R"), error = function(e) {
    cat("ERROR in step 16:", e$message, "\n")
    stop("Stacked model training failed.")
})

# Step 17: Internal hold-out evaluation
cat("\n--- Step 17: Internal Hold-Out Evaluation ---\n")
tryCatch(source("scripts/17_internal_holdout_evaluation.R"), error = function(e) {
    cat("ERROR in step 17:", e$message, "\n")
})

# Step 18: Bootstrap confidence intervals
cat("\n--- Step 18: Bootstrap CIs ---\n")
tryCatch(source("scripts/18_bootstrap_confidence_intervals.R"), error = function(e) {
    cat("ERROR in step 18:", e$message, "\n")
})

# Step 19: DeLong testing
cat("\n--- Step 19: DeLong Testing ---\n")
tryCatch(source("scripts/19_delong_testing.R"), error = function(e) {
    cat("ERROR in step 19:", e$message, "\n")
})

# Step 20: External validation
cat("\n--- Step 20: External Validation ---\n")
tryCatch(source("scripts/20_external_validation.R"), error = function(e) {
    cat("ERROR in step 20:", e$message, "\n")
})

# Step 21: GO/KEGG enrichment
cat("\n--- Step 21: GO/KEGG Enrichment ---\n")
tryCatch(source("scripts/21_GO_KEGG_enrichment.R"), error = function(e) {
    cat("ERROR in step 21:", e$message, "\n")
})

# Step 22: Pathview analysis
cat("\n--- Step 22: Pathview Analysis ---\n")
tryCatch(source("scripts/22_pathview_analysis.R"), error = function(e) {
    cat("ERROR in step 22:", e$message, "\n")
})

# Step 23: Generate figures
cat("\n--- Step 23: Generate Figures ---\n")
tryCatch(source("scripts/23_generate_figures.R"), error = function(e) {
    cat("ERROR in step 23:", e$message, "\n")
})

# Step 24: Render report
cat("\n--- Step 24: Render Report ---\n")
tryCatch(source("scripts/24_render_report.R"), error = function(e) {
    cat("ERROR in step 24:", e$message, "\n")
})

# Summary
end_time <- Sys.time()
duration <- difftime(end_time, start_time, units = "hours")

cat("\n========================================================\n")
cat("PIPELINE COMPLETE\n")
cat("========================================================\n")
cat("Start time:", format(start_time), "\n")
cat("End time:", format(end_time), "\n")
cat("Duration:", round(duration, 2), "hours\n")

# Check all results
cat("\n--- Output Summary ---\n")
result_dirs <- c(
    "cohort_suitability", "preprocessing", "quality_control", "PCA",
    "differential_expression", "data_splits", "feature_selection",
    "consensus_genes", "base_models", "stacked_model", "internal_validation",
    "bootstrap", "delong_tests", "external_validation", "enrichment",
    "pathview", "figures", "tables"
)
for (d in result_dirs) {
    dir_path <- file.path("results", d)
    files <- list.files(dir_path, pattern = "\\.(csv|rds|png|tiff|pdf)$")
    cat(sprintf("  %s: %d files\n", d, length(files)))
}

