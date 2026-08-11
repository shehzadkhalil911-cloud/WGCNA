# ============================================================
# 06 - Primary discovery split: patient-level 80/20
# Model development set vs untouched internal hold-out set.
# ============================================================
PROJECT_DIR <- "F:/independent projects/cancer/breast cancer/Breast cancer 8"
source(file.path(PROJECT_DIR, "config/config.R"))
source(file.path(PROJECT_DIR, "R/utils.R"))
suppressMessages({library(dplyr)})

logmsg <- function(...) cat(sprintf(...), "\n")
logmsg("=== 06 Primary discovery split ===")

disc <- readRDS(file.path(DIR_DISC_PROC, "discovery_processed.rds"))
meta <- disc$meta

set.seed(SEED)
patients <- unique(meta$patient_id)
n_hold <- round(length(patients) * HOLDOUT_FRACTION)

# stratified by patient-level ER subtype for balance
subtype <- sapply(split(meta$er, meta$patient_id), function(x) x[1])
set.seed(SEED)
hold_patients <- caret::createDataPartition(
  factor(subtype[patients]), p = HOLDOUT_FRACTION, list = FALSE
)
hold_patients <- patients[hold_patients]
train_patients <- setdiff(patients, hold_patients)

hold_samples <- meta$sample_id[meta$patient_id %in% hold_patients]
train_samples <- meta$sample_id[meta$patient_id %in% train_patients]

logmsg("Patients: total %d, development %d, hold-out %d",
       length(patients), length(train_patients), length(hold_patients))
logmsg("Samples: development %d, hold-out %d", length(train_samples), length(hold_samples))

stopifnot(length(intersect(train_samples, hold_samples)) == 0)
stopifnot(all(c(train_samples, hold_samples) %in% meta$sample_id))

# ---------- split table ----------
tab10 <- data.frame(
  sample_id = meta$sample_id,
  patient_id = meta$patient_id,
  class = meta$class,
  split = ifelse(meta$sample_id %in% train_samples, "development", "holdout"),
  stringsAsFactors = FALSE
)
tab10$index_in_discovery <- match(tab10$sample_id, colnames(disc$expr))
write.csv(tab10, file.path(DIR_TABLES, "Table_10_Discovery_Train_Holdout_Indices.csv"), row.names = FALSE)

saveRDS(list(train_samples = train_samples, hold_samples = hold_samples,
             train_patients = train_patients, hold_patients = hold_patients,
             seed = SEED, holdout_fraction = HOLDOUT_FRACTION),
        file.path(DIR_MODELS, "primary_split.rds"))

logmsg("Saved Table_10 and split indices.")
