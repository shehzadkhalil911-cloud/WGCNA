# ============================================================
# Configuration for the Breast cancer 8 study
# Stability-selected, calibrated XGBoost classifier
# distinguishing breast tumor vs normal tissue
# ============================================================

# ---------- paths ----------
PROJECT_DIR <- "F:/independent projects/cancer/breast cancer/Breast cancer 8"
DATA_DIR    <- "F:/independent projects/cancer/breast cancer/Datasets"

DIR_RAW          <- file.path(PROJECT_DIR, "data_processed", "raw")
DIR_DISC_PROC    <- file.path(PROJECT_DIR, "data_processed", "discovery")
DIR_VAL_PROC     <- file.path(PROJECT_DIR, "data_processed", "validation")
DIR_EXT_PROC     <- file.path(PROJECT_DIR, "data_processed", "external")
DIR_COMBAT_PROC  <- file.path(PROJECT_DIR, "data_processed", "combat")
DIR_COMBINED     <- file.path(PROJECT_DIR, "data_processed", "combined")
DIR_MODELS       <- file.path(PROJECT_DIR, "model_objects")
DIR_LOGS         <- file.path(PROJECT_DIR, "logs")
DIR_TABLES       <- file.path(PROJECT_DIR, "tables")
DIR_FIGURES      <- file.path(PROJECT_DIR, "figures")
DIR_REPORT       <- file.path(PROJECT_DIR, "report")
DIR_CONFIG       <- file.path(PROJECT_DIR, "config")

for (d in c(DIR_RAW, DIR_DISC_PROC, DIR_VAL_PROC, DIR_EXT_PROC, DIR_COMBAT_PROC,
            DIR_COMBINED, DIR_MODELS, DIR_LOGS, DIR_TABLES, DIR_FIGURES,
            DIR_REPORT, DIR_CONFIG)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---------- analysis data source ----------
# "raw"   = per-cohort quantile-normalized data (no cross-cohort correction)
# "combat"= combined-cohort ComBat batch-corrected data (primary analysis)
DATA_SOURCE <- "raw"

load_discovery <- function() {
  if (identical(DATA_SOURCE, "combat")) {
    readRDS(file.path(DIR_COMBAT_PROC, "discovery_combat.rds"))
  } else {
    readRDS(file.path(DIR_DISC_PROC, "discovery_processed.rds"))
  }
}
load_validation <- function() {
  if (identical(DATA_SOURCE, "combat")) {
    readRDS(file.path(DIR_COMBAT_PROC, "validation_combat.rds"))
  } else {
    readRDS(file.path(DIR_VAL_PROC, "validation_processed.rds"))
  }
}

# ---------- data sources (read-only) ----------
RAW_GSE70947_ESET <- file.path(DATA_DIR, "data/raw/GSE70947_eset.rds")
RAW_GSE15852_ESET <- file.path(DATA_DIR, "data/raw/GSE15852_eset.rds")
RAW_GSE29431_ESET <- file.path(DATA_DIR, "data/raw/GSE29431_eset.rds")
RAW_GSE42568_ESET <- file.path(DATA_DIR, "data/raw/GSE42568_eset.rds")
GPL13607_SOFT    <- file.path(DATA_DIR, "data/raw/GPL13607.soft.gz")
GPL96_SOFT       <- file.path(DATA_DIR, "data/raw/GPL96.soft.gz")
GPL570_SOFT      <- file.path(DATA_DIR, "data/raw/GPL570.soft.gz")
META_DIR         <- file.path(DATA_DIR, "data/metadata")

# ---------- global randomness ----------
SEED <- 123
set.seed(SEED)

# ---------- preprocessing ----------
MAX_NA_FRAC         <- 0.20   # probe kept if <= 20% missing values
COLLAPSE_METHOD     <- "maxIQR"  # collapse duplicate probes by max IQR
LOG2_FLOOR          <- 0.0
VALIDATE_SYMBOLS    <- TRUE  # validate symbols against org.Hs.eg.db

# ---------- differential expression ----------
DEG_PVAL    <- 0.05
DEG_LOGFC   <- 1.0
DEG_PAIRED  <- TRUE   # use paired analysis (verified 148 pairs)

# ---------- primary split ----------
HOLDOUT_FRACTION <- 0.20

# ---------- nested CV ----------
OUTER_FOLDS   <- 5
OUTER_REPEATS <- 5
INNER_FOLDS   <- 5

# ---------- leakage-safe feature selection ----------
FOLD_TOP_DEG         <- 1200  # top genes by fold-specific limma
CORR_CUTOFF          <- 0.90  # drop genes with pairwise |r| > cutoff
STABILITY_FREQ       <- 0.50  # retain genes selected in >=50% of inner folds
FINAL_STABILITY_FREQ <- 0.50  # retain genes selected in >=50% of outer folds
NZV_FREQCUT          <- 95/100  # near-zero variance frequency cutoff
NZV_UNIQUECUT        <- 10/100  # near-zero variance unique-value cutoff

# ---------- xgboost tuning ----------
TUNE_RANDOM_SEARCH_N <- 8        # random-search draws per outer fold
TUNE_NROUNDS_MAX     <- 250
TUNE_EARLY_STOP      <- 15
TUNE_METRIC          <- "auc"    # primary tuning metric

# ---------- calibration ----------
CALIBRATION_METHOD <- "isotonic" # locked after comparison (see 08)
# alternative: "platt" (logistic regression on logit(p))

# ---------- threshold ----------
THRESHOLD_METHOD <- "youden"     # Youden J on OOF development predictions

# ---------- bootstrap ----------
N_BOOT       <- 2000   # internal hold-out bootstrap
N_BOOT_EXT   <- 2000   # external validation bootstrap

# ---------- parallel ----------
N_CORES <- parallel::detectCores()
if (is.na(N_CORES) || N_CORES < 1) N_CORES <- 2
N_CORES <- min(N_CORES, 4)

# ---------- class/colour constants ----------
CLASS_LEVELS   <- c("Tumor", "Normal")
COL_TUMOR      <- "#D55E00"
COL_NORMAL     <- "#0072B2"
PAL_OKABEITO   <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                    "#0072B2", "#D55E00", "#CC79A7", "#000000")

# ---------- figure defaults ----------
FIG_DPI    <- 600
FIG_WHITE  <- TRUE
