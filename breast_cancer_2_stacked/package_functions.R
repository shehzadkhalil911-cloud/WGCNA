# =============================================================================
# package_functions.R
# Package installation and loading functions
# =============================================================================

#' Install missing CRAN packages
install_cran_packages <- function(pkgs, repos = "https://cran.r-project.org") {
    installed <- rownames(installed.packages())
    to_install <- setdiff(pkgs, installed)
    if (length(to_install) > 0) {
        message("Installing CRAN packages: ", paste(to_install, collapse = ", "))
        install.packages(to_install, dependencies = TRUE, repos = repos)
    }
    invisible(to_install)
}

#' Install missing Bioconductor packages
install_bioc_packages <- function(pkgs) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    installed <- rownames(installed.packages())
    to_install <- setdiff(pkgs, installed)
    if (length(to_install) > 0) {
        message("Installing Bioconductor packages: ", paste(to_install, collapse = ", "))
        BiocManager::install(to_install, update = FALSE, ask = FALSE)
    }
    invisible(to_install)
}

#' Load all required packages
load_required_packages <- function() {
    cran_pkgs <- c(
        "tidyverse", "data.table", "caret", "caretEnsemble",
        "ranger", "gbm", "kernlab", "e1071", "kknn", "nnet",
        "glmnet", "ROSE", "recipes", "themis", "rsample",
        "doParallel", "foreach", "doRNG", "pROC", "PRROC",
        "yardstick", "boot", "MLmetrics", "ggplot2", "ggpubr",
        "ggrepel", "patchwork", "cowplot", "pheatmap",
        "circlize", "VennDiagram", "ggvenn", "viridis",
        "RColorBrewer", "officer", "flextable", "quarto",
        "rmarkdown", "renv"
    )
    bioc_pkgs <- c(
        "GEOquery", "Biobase", "limma", "sva",
        "preprocessCore", "AnnotationDbi", "annotate",
        "org.Hs.eg.db", "ComplexHeatmap", "EnhancedVolcano",
        "clusterProfiler", "enrichplot", "DOSE", "pathview",
        "KEGGREST", "ComplexUpset", "DiagrammeR"
    )
    # Check for complex packages
    extra_cran <- c("corrplot", "gridExtra", "scales", "graphics")
    cran_pkgs <- unique(c(cran_pkgs, extra_cran))

    install_cran_packages(cran_pkgs)
    install_bioc_packages(bioc_pkgs)

    library(tidyverse)
    library(data.table)
    library(caret)
    library(ranger)
    library(gbm)
    library(kernlab)
    library(kknn)
    library(nnet)
    library(glmnet)
    library(pROC)
    library(preprocessCore)
    library(limma)
    library(sva)
    library(GEOquery)
    library(Biobase)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    library(clusterProfiler)
    library(enrichplot)
    library(pathview)
    library(ComplexHeatmap)
    library(patchwork)
    library(cowplot)
    library(viridis)
    library(flextable)
    library(officer)

    message("All packages loaded successfully.")
}
