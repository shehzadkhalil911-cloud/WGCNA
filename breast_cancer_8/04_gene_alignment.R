# ============================================================
# 04 - Cross-cohort gene alignment (4 cohorts)
# Table_05 (combined preprocessing), Table_06 (annotation),
# Table_07 (common genes across all four cohorts) and
# alignment figures.
# ============================================================
PROJECT_DIR <- "F:/independent projects/cancer/breast cancer/Breast cancer 8"
source(file.path(PROJECT_DIR, "config/config.R"))
source(file.path(PROJECT_DIR, "R/utils.R"))
suppressMessages({library(ggplot2); library(dplyr); library(ggVennDiagram)})

logmsg <- function(...) cat(sprintf(...), "\n")
logmsg("=== 04 Gene alignment (4 cohorts) ===")

disc <- readRDS(file.path(DIR_DISC_PROC, "discovery_processed.rds"))
val  <- readRDS(file.path(DIR_VAL_PROC, "validation_processed.rds"))
ext1 <- readRDS(file.path(DIR_EXT_PROC, "GSE29431_processed.rds"))
ext2 <- readRDS(file.path(DIR_EXT_PROC, "GSE42568_processed.rds"))

genes <- list(
  Discovery = rownames(disc[["expr"]]),
  GSE15852  = rownames(val[["expr"]]),
  GSE29431  = rownames(ext1[["expr"]]),
  GSE42568  = rownames(ext2[["expr"]])
)
common    <- Reduce(intersect, genes)
n_common  <- length(common)
absent_each <- sapply(genes, function(g) length(setdiff(genes$Discovery, g)))

for (nm in names(genes)) {
  logmsg("%-9s genes: %d", nm, length(genes[[nm]]))
}
logmsg("Common genes across all four cohorts: %d", n_common)

# ---------- Table_05 combined preprocessing summary ----------
t5a <- read.csv(file.path(DIR_TABLES, "Table_05a_Preprocessing_Summary_Discovery.csv"))
t5b1 <- read.csv(file.path(DIR_TABLES, "Table_05b_Preprocessing_Summary_GSE15852.csv"))
t5b2 <- read.csv(file.path(DIR_TABLES, "Table_05b_Preprocessing_Summary_GSE29431.csv"))
t5b3 <- read.csv(file.path(DIR_TABLES, "Table_05b_Preprocessing_Summary_GSE42568.csv"))
tab5 <- dplyr::bind_rows(t5a, t5b1, t5b2, t5b3)
write.csv(tab5, file.path(DIR_TABLES, "Table_05_Preprocessing_Summary.csv"), row.names = FALSE)
for (f in c("Table_05a_Preprocessing_Summary_Discovery.csv",
            "Table_05b_Preprocessing_Summary_GSE15852.csv",
            "Table_05b_Preprocessing_Summary_GSE29431.csv",
            "Table_05b_Preprocessing_Summary_GSE42568.csv")) {
  file.remove(file.path(DIR_TABLES, f))
}

# ---------- Table_06 annotation summary ----------
tab6 <- data.frame(
  item = c(paste0("platform_", names(genes)),
           paste0("genes_final_", names(genes)),
           "genes_common_all_four"),
  value = c("GPL13607", "GPL96", "GPL570", "GPL570",
            sapply(genes, length),
            n_common),
  stringsAsFactors = FALSE
)
write.csv(tab6, file.path(DIR_TABLES, "Table_06_Gene_Annotation_Summary.csv"), row.names = FALSE)

# ---------- Table_07 common genes ----------
tab7 <- data.frame(gene = common, stringsAsFactors = FALSE)
write.csv(tab7, file.path(DIR_TABLES, "Table_07_Common_Genes.csv"), row.names = FALSE)
saveRDS(common, file.path(DIR_MODELS, "common_genes.rds"))

# ---------- alignment figures (4-set Venn) ----------
p <- ggVennDiagram(genes, label_alpha = 0, set_size = 3.5, label = "count") +
  scale_fill_gradient2(low = "#D55E00", mid = "#EEDDCC", high = "#0072B2") +
  ggtitle("Gene overlap across all four cohorts") +
  theme_bc() + theme(legend.position = "none")
save_tiff_600dpi(p, "Figure_06_Gene_Overlap.tiff", width = 8, height = 6)

# bar: discovery genes absent from each validation cohort
absent_df <- data.frame(
  cohort = names(absent_each),
  count = unname(absent_each)
)
p2 <- ggplot(absent_df, aes(x = cohort, y = count, fill = cohort)) +
  geom_col(colour = "black", width = 0.6) +
  geom_text(aes(label = count), vjust = -0.4, size = 5) +
  scale_fill_manual(values = c(GSE15852 = "#CC79A7", GSE29431 = "#0072B2", GSE42568 = "#009E73")) +
  theme_bc() + labs(x = "", y = "Number of genes",
                    title = "Discovery genes not measured on each validation platform") +
  theme(legend.position = "none")
save_tiff_600dpi(p2, "Figure_06_Absent_Genes.tiff", width = 7, height = 6)

logmsg("Saved Table_05/06/07 and alignment figures.")
