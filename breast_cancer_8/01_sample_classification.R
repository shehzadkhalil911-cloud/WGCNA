# ============================================================
# 01 - Sample classification and verification
# Builds Table_02/03/04 and Figure_02 for all four cohorts
# (GSE70947 development; GSE15852/GSE29431/GSE42568 validation).
# STOPS with an error if tumor/normal assignment is unreliable
# for the paired cohorts.
# ============================================================
PROJECT_DIR <- "F:/independent projects/cancer/breast cancer/Breast cancer 8"
source(file.path(PROJECT_DIR, "config/config.R"))
source(file.path(PROJECT_DIR, "R/utils.R"))
suppressMessages({library(Biobase); library(ggplot2); library(dplyr)})

logmsg <- function(...) cat(sprintf(...), "\n")

logmsg("=== 01 Sample classification ===")

# ---------------- GSE70947 (discovery) ----------------
es1 <- readRDS(RAW_GSE70947_ESET)
pd1 <- pData(es1)

tissue1 <- as.character(pd1$`tissue:ch1`)
stopifnot(all(tissue1 %in% c("normal", "tumor")))
class1 <- ifelse(tissue1 == "tumor", "Tumor", "Normal")
patient1 <- sub("-(normal|tumor)$", "", pd1$title, ignore.case = TRUE)

# ---------------- GSE15852 (validation, paired) ----------------
es2 <- readRDS(RAW_GSE15852_ESET)
pd2 <- pData(es2)
cls2 <- ifelse(grepl("^Normal", pd2$title), "Normal", ifelse(grepl("^Cancer", pd2$title), "Tumor", NA))
stopifnot(all(!is.na(cls2)))
patient2 <- sub("^(Normal|Cancer) (BC[0-9]+)[NT]$", "\\2", pd2$title)

# ---------------- GSE29431 / GSE42568 (validation, unpaired) ----------------
cls_csv <- function(path, es) {
  cc <- read.csv(path, stringsAsFactors = FALSE)
  cc <- cc[match(sampleNames(es), cc$sample_id), ]
  stopifnot(all(sampleNames(es) == cc$sample_id))
  cc$class
}

es3 <- readRDS(RAW_GSE29431_ESET)
es4 <- readRDS(RAW_GSE42568_ESET)
cls3 <- cls_csv(file.path(META_DIR, "GSE29431_sample_classes.csv"), es3)
cls4 <- cls_csv(file.path(META_DIR, "GSE42568_sample_classes.csv"), es4)
stopifnot(all(cls3 %in% CLASS_LEVELS), all(cls4 %in% CLASS_LEVELS))
patient3 <- sampleNames(es3)
patient4 <- sampleNames(es4)

# ---------------- build assignment table ----------------
mk_table3 <- function(dataset, gsm, title, source_name, characteristics, assigned_class,
                      patient_id, inclusion = "Included", excl_reason = NA, confidence = "High") {
  data.frame(
    dataset            = dataset,
    sample_id          = gsm,
    original_title     = title,
    source_name        = source_name,
    characteristics    = characteristics,
    assigned_class     = assigned_class,
    patient_id         = patient_id,
    inclusion_status   = inclusion,
    exclusion_reason   = excl_reason,
    confidence         = confidence,
    stringsAsFactors   = FALSE
  )
}

t3a <- mk_table3(
  dataset = "GSE70947", gsm = pd1$geo_accession, title = pd1$title,
  source_name = pd1$source_name_ch1,
  characteristics = paste0("tissue=", tissue1, "; sample_id=", pd1$`sample id:ch1`, "; er=", pd1$`er:ch1`, "; her2=", pd1$`her2:ch1`),
  assigned_class = class1, patient_id = patient1
)
t3b <- mk_table3(
  dataset = "GSE15852", gsm = pd2$geo_accession, title = pd2$title,
  source_name = pd2$source_name_ch1,
  characteristics = paste0("histopath=", pd2$`histopathological exam:ch1`, "; grade=", pd2$`grade:ch1`, "; age=", pd2$`age in years:ch1`),
  assigned_class = cls2, patient_id = patient2
)
t3c <- mk_table3(
  dataset = "GSE29431", gsm = sampleNames(es3), title = pData(es3)$title,
  source_name = pData(es3)$source_name_ch1,
  characteristics = "unpaired; class from GEO series metadata",
  assigned_class = cls3, patient_id = patient3
)
t3d <- mk_table3(
  dataset = "GSE42568", gsm = sampleNames(es4), title = pData(es4)$title,
  source_name = pData(es4)$source_name_ch1,
  characteristics = "unpaired; class from GEO series metadata",
  assigned_class = cls4, patient_id = patient4
)

tab3 <- rbind(t3a, t3b, t3c, t3d)

# ---------------- verify pairing / duplicates (paired cohorts only) ----------------
paired_ok <- function(ds, patient_vec, cls_vec) {
  v1 <- tapply(cls_vec, patient_vec, function(z) sort(z))
  bad <- names(which(!sapply(v1, function(z) identical(unname(z), c("Normal", "Tumor")))))
  bad
}
bad1 <- paired_ok("GSE70947", tab3$patient_id[tab3$dataset == "GSE70947"],
                  tab3$assigned_class[tab3$dataset == "GSE70947"])
bad2 <- paired_ok("GSE15852", tab3$patient_id[tab3$dataset == "GSE15852"],
                  tab3$assigned_class[tab3$dataset == "GSE15852"])

if (length(bad1) || length(bad2)) {
  write.csv(tab3, file.path(DIR_TABLES, "Table_03_Sample_Group_Assignments.csv"), row.names = FALSE)
  stop("PAIRING VIOLATION: patients without exactly one tumor + one normal. See Table_03.")
}

dup_gsm <- tab3$sample_id[duplicated(tab3$sample_id)]
if (length(dup_gsm)) {
  write.csv(tab3, file.path(DIR_TABLES, "Table_03_Sample_Group_Assignments.csv"), row.names = FALSE)
  stop("DUPLICATE SAMPLE IDS: ", paste(dup_gsm, collapse = ", "))
}

for (ds in c("GSE70947", "GSE15852", "GSE29431", "GSE42568")) {
  sub <- tab3[tab3$dataset == ds, ]
  logmsg("%s: %d samples, %d patients", ds, nrow(sub), length(unique(sub$patient_id)))
}

# ---------------- Table_04 excluded samples ----------------
tab4 <- tab3[tab3$inclusion_status != "Included", , drop = FALSE]
if (nrow(tab4) == 0) {
  tab4 <- data.frame(note = "No samples excluded; all samples verified and included.",
                     stringsAsFactors = FALSE)
}

# ---------------- Table_02 dataset characteristics ----------------
tab2 <- data.frame(
  dataset              = c("GSE70947", "GSE15852", "GSE29431", "GSE42568"),
  role                 = c("Development", "Validation", "Validation", "Validation"),
  platform             = c("GPL13607 (Agilent SurePrint G3 Human GE 8x60K)",
                           "GPL96 (Affymetrix Human Genome U133A)",
                           "GPL570 (Affymetrix Human Genome U133 Plus 2.0)",
                           "GPL570 (Affymetrix Human Genome U133 Plus 2.0)"),
  n_samples_total      = c(296, 86, 66, 121),
  n_patients           = c(148, 43, 66, 121),
  n_tumor              = c(148, 43, 54, 104),
  n_normal             = c(148, 43, 12, 17),
  paired               = c(TRUE, TRUE, FALSE, FALSE),
  tumor_label_source   = c("tissue:ch1 == tumor", "title prefix 'Cancer' + source_name",
                           "GEO series metadata", "GEO series metadata"),
  normal_label_source  = c("tissue:ch1 == normal", "title prefix 'Normal' + source_name",
                           "GEO series metadata", "GEO series metadata"),
  expression_scale     = c("log2 (background-corrected, quantile-normalized)",
                           "MAS5 linear (needs log2)",
                           "log2 (processed)", "log2 (processed)"),
  n_probes_raw         = c(62976, 22283, 54675, 54675),
  stringsAsFactors     = FALSE
)

# ---------------- Figure_02 sample composition ----------------
plot_tab <- tab3 %>%
  group_by(dataset, assigned_class) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(dataset = factor(dataset, levels = c("GSE70947", "GSE15852", "GSE29431", "GSE42568")),
         assigned_class = factor(assigned_class, levels = CLASS_LEVELS))

fig2 <- ggplot(plot_tab, aes(x = dataset, y = n, fill = assigned_class)) +
  geom_col(position = "dodge", colour = "black", width = 0.7) +
  geom_text(aes(label = n), position = position_dodge(width = 0.7), vjust = -0.4, size = 4) +
  scale_fill_manual(values = c(Tumor = COL_TUMOR, Normal = COL_NORMAL),
                    name = "Assigned class") +
  labs(x = "Dataset", y = "Number of samples",
       title = "Sample composition by dataset and class") +
  theme_bc() +
  theme(legend.position = "bottom")
save_tiff_600dpi(fig2, "Figure_02_Sample_Composition.tiff", width = 8, height = 6)

# ---------------- save ----------------
write.csv(tab2, file.path(DIR_TABLES, "Table_02_Dataset_Characteristics.csv"), row.names = FALSE)
write.csv(tab3, file.path(DIR_TABLES, "Table_03_Sample_Group_Assignments.csv"), row.names = FALSE)
write.csv(tab4, file.path(DIR_TABLES, "Table_04_Excluded_Samples.csv"), row.names = FALSE)

saveRDS(tab3, file.path(DIR_MODELS, "sample_assignments.rds"))
logmsg("Saved Table_02/03/04 and Figure_02. All checks passed.")
