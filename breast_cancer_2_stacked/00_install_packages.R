setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 00_install_packages.R
# Install all required packages
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")

cran_pkgs <- c(
    "tidyverse", "data.table", "caret", "caretEnsemble",
    "ranger", "gbm", "kernlab", "e1071", "kknn", "nnet",
    "glmnet", "ROSE", "recipes", "themis", "rsample",
    "doParallel", "foreach", "doRNG", "pROC", "PRROC",
    "yardstick", "boot", "MLmetrics", "ggplot2", "ggpubr",
    "ggrepel", "patchwork", "cowplot", "pheatmap",
    "circlize", "VennDiagram", "ggvenn", "viridis",
    "RColorBrewer", "officer", "flextable", "quarto",
    "rmarkdown", "renv", "corrplot", "gridExtra", "scales"
)

bioc_pkgs <- c(
    "GEOquery", "Biobase", "limma", "sva",
    "preprocessCore", "AnnotationDbi", "annotate",
    "org.Hs.eg.db", "ComplexHeatmap", "EnhancedVolcano",
    "clusterProfiler", "enrichplot", "DOSE", "pathview",
    "KEGGREST", "ComplexUpset", "DiagrammeR"
)

install_cran_packages(cran_pkgs)
install_bioc_packages(bioc_pkgs)

# Verify installation
missing <- setdiff(c(cran_pkgs, bioc_pkgs),
                   rownames(installed.packages()))
if (length(missing) > 0) {
    warning("Missing packages: ", paste(missing, collapse = ", "))
} else {
    message("All packages installed successfully.")
}

# Save session info
sink("sessionInfo.txt")
sessionInfo()
sink()

