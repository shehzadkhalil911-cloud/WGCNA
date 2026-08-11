setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 04_construct_discovery_cohort.R
# Build the discovery cohort from three datasets and prepare cohort suitability
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")
source("R/annotation_functions.R")
source("R/plotting_functions.R")
load_required_packages()

# --- Identify tumor/normal samples for each dataset ---
# This requires manual inspection of phenotype data
# We'll use the known sample annotations for these GEO series

identify_sample_classes <- function(pdata, gse) {
    # Dynamic class identification based on available metadata columns
    all_text <- apply(pdata, 2, function(x) paste(as.character(x), collapse = " | "))
    combined <- paste(all_text, collapse = " ")
    
    # Strategy: search systematically across all columns for tumor/normal indicators
    class <- rep(NA, nrow(pdata))
    
    # First, check specific known columns for each dataset
    if (gse == "GSE42568") {
        source <- as.character(pdata$source_name_ch1)
        class[grepl("cancer", source, ignore.case = TRUE)] <- "Tumor"
        class[grepl("normal", source, ignore.case = TRUE)] <- "Normal"
    } else if (gse == "GSE15852") {
        # Use source_name_ch1 which has "normal breast tissue" or "breast tumor tissue"
        source <- as.character(pdata$source_name_ch1)
        class[grepl("tumor", source, ignore.case = TRUE)] <- "Tumor"
        class[grepl("normal", source, ignore.case = TRUE)] <- "Normal"
    } else if (gse == "GSE29431") {
        title <- as.character(pdata$title)
        class[grepl("tumor|cancer|breast cancer", title, ignore.case = TRUE)] <- "Tumor"
        class[grepl("normal|healthy", title, ignore.case = TRUE)] <- "Normal"
    } else if (gse == "GSE70947") {
        # Use the tissue:ch1 column if available
        if ("tissue:ch1" %in% colnames(pdata)) {
            tissue <- as.character(pdata[["tissue:ch1"]])
            class[grepl("tumor|cancer", tissue, ignore.case = TRUE)] <- "Tumor"
            class[grepl("normal", tissue, ignore.case = TRUE)] <- "Normal"
        } else {
            title <- as.character(pdata$title)
            class[grepl("tumor|cancer", title, ignore.case = TRUE)] <- "Tumor"
            class[grepl("normal", title, ignore.case = TRUE)] <- "Normal"
        }
    } else {
        # Generic: search all columns
        for (col in colnames(pdata)) {
            vals <- as.character(pdata[[col]])
            tumor_idx <- grepl("tumor|cancer", vals, ignore.case = TRUE)
            normal_idx <- grepl("normal", vals, ignore.case = TRUE)
            if (any(tumor_idx) || any(normal_idx)) {
                class[tumor_idx] <- "Tumor"
                class[normal_idx] <- "Normal"
                if (sum(!is.na(class)) > nrow(pdata) * 0.5) break
            }
        }
    }
    
    # Remove conflicting assignments
    conflict_idx <- which(!is.na(class) & class == "Tumor" & 
                           grepl("normal", apply(pdata, 1, paste, collapse = " "), ignore.case = TRUE) &
                           !grepl("tumor|cancer", apply(pdata, 1, paste, collapse = " "), ignore.case = TRUE))
    class[conflict_idx] <- NA
    
    return(data.frame(sample_id = rownames(pdata), class = class,
                      stringsAsFactors = FALSE))
}

# Load processed datasets from saved files
all_gse <- c(DISCOVERY_GSE, VALIDATION_GSE)
sample_info_list <- list()
expression_list <- list()
processed_list <- list()

for (gse in all_gse) {
    rds_path <- file.path(PROCESSED_DIR, paste0(gse, "_processed.rds"))
    if (!file.exists(rds_path)) {
        stop("Processed data not found for ", gse, ". Run script 03 first.")
    }
    processed_list[[gse]] <- readRDS(rds_path)
    data <- processed_list[[gse]]
    pdata <- data$phenotype
    expr <- data$expression

    # Identify classes
    sample_classes <- identify_sample_classes(pdata, gse)
    cat("\n", gse, "sample classes:\n")
    print(table(sample_classes$class, useNA = "ifany"))

    sample_info_list[[gse]] <- sample_classes
    expression_list[[gse]] <- expr

    # Save sample mapping
    write.csv(sample_classes,
              file.path(METADATA_DIR, paste0(gse, "_sample_classes.csv")),
              row.names = FALSE)
}

# --- Check common genes ---
common_genes <- Reduce(intersect, lapply(expression_list, rownames))
cat("\nCommon genes across all datasets:", length(common_genes), "\n")

# Discovery-only common genes
disc_expr_list <- expression_list[names(expression_list) %in% DISCOVERY_GSE]
disc_common_genes <- Reduce(intersect, lapply(disc_expr_list, rownames))
cat("Common genes across discovery datasets:", length(disc_common_genes), "\n")

# --- Build discovery cohort ---
discovery_samples <- do.call(rbind, sample_info_list[DISCOVERY_GSE])
discovery_metadata <- discovery_samples
discovery_metadata$dataset <- rep(DISCOVERY_GSE,
                                  times = sapply(sample_info_list[DISCOVERY_GSE], nrow))
discovery_metadata$platform <- sapply(discovery_samples$sample_id, function(sid) {
    for (gse in names(processed_list)) {
        if (sid %in% colnames(processed_list[[gse]]$expression)) {
            return(processed_list[[gse]]$platform)
        }
    }
    return(NA)
})

# Restrict to common genes
discovery_expr <- do.call(cbind, lapply(DISCOVERY_GSE, function(gse) {
    expression_list[[gse]][disc_common_genes, ]
}))
colnames(discovery_expr) <- discovery_samples$sample_id

# Remove any samples with NA class
valid_idx <- !is.na(discovery_metadata$class)
discovery_expr <- discovery_expr[, valid_idx]
discovery_metadata <- discovery_metadata[valid_idx, ]

cat("\nFinal discovery cohort:\n")
cat("  Samples:", ncol(discovery_expr), "\n")
cat("  Genes:", nrow(discovery_expr), "\n")
print(table(discovery_metadata$class, discovery_metadata$dataset))

# --- Build validation cohort ---
val_gse <- VALIDATION_GSE
val_samples <- sample_info_list[[val_gse]]
common_val_genes <- intersect(disc_common_genes, rownames(expression_list[[val_gse]]))
cat("Genes common between discovery and validation:", length(common_val_genes), "\n")
val_expr <- expression_list[[val_gse]][common_val_genes, ]
val_valid <- !is.na(val_samples$class)
val_expr <- val_expr[, val_valid]
val_samples <- val_samples[val_valid, ]

cat("\nValidation cohort:\n")
cat("  Samples:", ncol(val_expr), "\n")
cat("  Genes:", nrow(val_expr), "\n")
print(table(val_samples$class))

# --- Cohort suitability check ---
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("COHORT SUITABILITY ASSESSMENT\n")
cat(paste(rep("=", 70), collapse = ""), "\n")

suitability <- data.frame(
    Dataset = c(DISCOVERY_GSE, VALIDATION_GSE),
    Role = c(rep("Discovery", 3), "Validation"),
    Platform = sapply(processed_list, function(x) x$platform),
    Total_Samples = sapply(sample_info_list, nrow),
    Tumor_Samples = sapply(sample_info_list, function(x) sum(x$class == "Tumor", na.rm = TRUE)),
    Normal_Samples = sapply(sample_info_list, function(x) sum(x$class == "Normal", na.rm = TRUE)),
    Excluded_Samples = sapply(sample_info_list, function(x) sum(is.na(x$class))),
    Probes = sapply(processed_list, function(x) nrow(x$expression_norm)),
    Mapped_Genes = sapply(processed_list, function(x) nrow(x$expression)),
    Common_Genes = c(rep(length(disc_common_genes), 3), length(common_val_genes)),
    Has_Both_Classes = sapply(sample_info_list, function(x) {
        "Tumor" %in% x$class && "Normal" %in% x$class
    })
)

cat("\nCohort Suitability Table:\n")
print(suitability)

write.csv(suitability, file.path(COHORT_SUIT_DIR, "cohort_suitability.csv"), row.names = FALSE)
saveRDS(suitability, file.path(COHORT_SUIT_DIR, "cohort_suitability.rds"))

# Check validation suitability
val_suitable <- suitability[suitability$Dataset == VALIDATION_GSE, ]
if (val_suitable$Has_Both_Classes && val_suitable$Common_Genes > 100) {
    cat("\nVALIDATION COHORT SUITABLE:", VALIDATION_GSE, "\n")
    cat("  Has both tumor and normal samples:", val_suitable$Has_Both_Classes, "\n")
    cat("  Common genes:", val_suitable$Common_Genes, "\n")
} else {
    cat("\n*** WARNING: Validation cohort may not be suitable ***\n")
    cat("  Has both classes:", val_suitable$Has_Both_Classes, "\n")
    cat("  Common genes:", val_suitable$Common_Genes, "\n")
}

# Save discovery and validation data
saveRDS(list(expression = discovery_expr, metadata = discovery_metadata),
        file.path(DISCOVERY_DIR, "discovery_cohort.rds"))
saveRDS(list(expression = val_expr, metadata = val_samples),
        file.path(EXTVAL_DIR, "validation_cohort.rds"))

cat("\nDiscovery and validation cohorts saved.\n")

