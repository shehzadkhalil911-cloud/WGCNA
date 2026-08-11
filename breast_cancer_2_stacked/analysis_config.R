# =============================================================================
# analysis_config.R
# Master configuration file for Breast Cancer Stacked Ensemble Analysis
# =============================================================================

# --- Cohort Assignment -------------------------------------------------------
DISCOVERY_GSE    <- c("GSE42568", "GSE15852", "GSE29431")
VALIDATION_GSE   <- "GSE70947"

# --- Reproducibility ---------------------------------------------------------
GLOBAL_SEED      <- 123

# --- Directories -------------------------------------------------------------
BASE_DIR         <- getwd()
DATA_DIR         <- file.path(BASE_DIR, "data")
RAW_DIR          <- file.path(DATA_DIR, "raw")
METADATA_DIR     <- file.path(DATA_DIR, "metadata")
PROCESSED_DIR    <- file.path(DATA_DIR, "processed")
DISCOVERY_DIR    <- file.path(DATA_DIR, "discovery")
EXTVAL_DIR       <- file.path(DATA_DIR, "external_validation")
R_DIR            <- file.path(BASE_DIR, "R")
SCRIPTS_DIR      <- file.path(BASE_DIR, "scripts")
RESULTS_DIR      <- file.path(BASE_DIR, "results")
REPORT_DIR       <- file.path(BASE_DIR, "report")

# Sub-results directories
COHORT_SUIT_DIR  <- file.path(RESULTS_DIR, "cohort_suitability")
PREPROC_DIR      <- file.path(RESULTS_DIR, "preprocessing")
QC_DIR           <- file.path(RESULTS_DIR, "quality_control")
PCA_DIR          <- file.path(RESULTS_DIR, "PCA")
DE_DIR           <- file.path(RESULTS_DIR, "differential_expression")
SPLITS_DIR       <- file.path(RESULTS_DIR, "data_splits")
FS_DIR           <- file.path(RESULTS_DIR, "feature_selection")
CONSENSUS_DIR    <- file.path(RESULTS_DIR, "consensus_genes")
BASE_MODELS_DIR  <- file.path(RESULTS_DIR, "base_models")
STACKED_DIR      <- file.path(RESULTS_DIR, "stacked_model")
INTVAL_DIR       <- file.path(RESULTS_DIR, "internal_validation")
BOOT_DIR         <- file.path(RESULTS_DIR, "bootstrap")
DELONG_DIR       <- file.path(RESULTS_DIR, "delong_tests")
EXTVAL_RES_DIR   <- file.path(RESULTS_DIR, "external_validation")
ENRICH_DIR       <- file.path(RESULTS_DIR, "enrichment")
PATHVIEW_DIR     <- file.path(RESULTS_DIR, "pathview")
FIGURES_DIR      <- file.path(RESULTS_DIR, "figures")
TABLES_DIR       <- file.path(RESULTS_DIR, "tables")

# --- Differential Expression Thresholds --------------------------------------
ADJ_P_THRESHOLD  <- 0.05
LOG2FC_THRESHOLD <- 1.0

# --- Feature Selection -------------------------------------------------------
MIN_SELECTION_FREQUENCY <- 0.75
MAX_FEATURE_CANDIDATES  <- 2000
TOP_RANKED_FEATURES     <- 500

# --- Duplicate Probe Resolution ----------------------------------------------
DUPLICATE_PROBE_METHOD  <- "highest_IQR"

# --- Cross-Validation --------------------------------------------------------
CV_FOLDS         <- 10
CV_REPEATS       <- 10
N_BOOTSTRAP      <- 2000

# --- Train/Test Allocations --------------------------------------------------
ALLOCATIONS <- c(0.60, 0.70, 0.80, 0.90)

# --- Machine Learning Algorithms ---------------------------------------------
ALGORITHMS <- c("rf", "gbm", "svmLinear", "knn", "nnet")

# --- Parallel Computing ------------------------------------------------------
N_CORES          <- max(1, parallel::detectCores() - 1)
USE_PARALLEL     <- TRUE

# --- File paths for saved objects --------------------------------------------
TRAIN_INDICES_FILE   <- file.path(SPLITS_DIR, "train_indices.rds")
TEST_INDICES_FILE    <- file.path(SPLITS_DIR, "test_indices.rds")
CV_INDICES_FILE      <- file.path(SPLITS_DIR, "cv_indices.rds")
OOF_FILE             <- file.path(STACKED_DIR, "oof_predictions.rds")
STACKED_MODEL_FILE   <- file.path(STACKED_DIR, "stacked_elastic_net.rds")
CONSENSUS_FILE       <- file.path(CONSENSUS_DIR, "consensus_genes.rds")
FINAL_MODEL_FILE     <- file.path(BASE_MODELS_DIR, "final_models.rds")
META_LEARNER_FILE    <- file.path(STACKED_DIR, "meta_learner.rds")

# --- Plotting Defaults -------------------------------------------------------
FIGURE_WIDTH        <- 10
FIGURE_HEIGHT       <- 8
FIGURE_DPI          <- 600
TUMOR_COLOR         <- "#E41A1C"
NORMAL_COLOR        <- "#377EB8"
MODEL_COLORS        <- c("RF" = "#E41A1C", "GB" = "#377EB8",
                         "SVM" = "#4DAF4A", "KNN" = "#984EA3",
                         "NN" = "#FF7F00", "Stacked" = "#A65628")
PALETTE_DISCRETE    <- "Set1"
PALETTE_CONTINUOUS  <- "RdBu"

# --- Model Training Control Defaults -----------------------------------------
RF_TUNE_GRID <- expand.grid(
    mtry = c(2, 4, 6, 8, 10, 15, 20),
    splitrule = "gini",
    min.node.size = c(1, 3, 5)
)

GB_TUNE_GRID <- expand.grid(
    n.trees = c(100, 500, 1000),
    interaction.depth = c(1, 3, 5, 7),
    shrinkage = c(0.01, 0.05, 0.1),
    n.minobsinnode = c(5, 10, 20)
)

SVM_TUNE_GRID <- expand.grid(C = 2^(-5:5))

KNN_TUNE_GRID <- expand.grid(
    k = c(3, 5, 7, 9, 11, 15, 21, 31)
)

NN_TUNE_GRID <- expand.grid(
    size = c(1, 3, 5, 7, 9, 11),
    decay = c(0, 0.001, 0.01, 0.1, 0.5)
)

# --- Elastic Net Meta-Learner Tuning -----------------------------------------
ELASTIC_NET_ALPHA <- seq(0, 1, by = 0.1)
ELASTIC_NET_LAMBDA <- NULL  # Let glmnet determine

# --- Positive Class ----------------------------------------------------------
POSITIVE_CLASS  <- "Tumor"
NEGATIVE_CLASS  <- "Normal"
