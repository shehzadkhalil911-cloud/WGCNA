source("R/utils.R")
set.seed(123)
options(stringsAsFactors = FALSE)

required <- c(
  "limma", "Biobase", "GEOquery", "sva", "WGCNA", "ggplot2", "pheatmap",
  "reshape2", "dplyr", "tidyr", "tibble", "readxl", "writexl", "openxlsx",
  "org.Hs.eg.db", "AnnotationDbi", "clusterProfiler", "enrichplot", "ReactomePA",
  "GSVA", "ComplexHeatmap", "circlize", "RColorBrewer", "corrplot",
  "igraph", "stringr", "gridExtra", "cowplot", "ggrepel", "VennDiagram",
  "UpSetR", "scales"
)

missing <- c()
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing <- c(missing, pkg)
  }
}

if (length(missing) > 0) {
  log_msg(sprintf("Missing packages: %s", paste(missing, collapse = ", ")))
  log_msg("Installing missing packages...")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(missing, update = FALSE, ask = FALSE)
}

for (pkg in required) {
  library(pkg, character.only = TRUE)
}

log_msg("All required packages loaded successfully")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
