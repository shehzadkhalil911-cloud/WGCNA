# ============================================================
# Setup: Install packages, define paths, load libraries
# Breast Cancer 5 Necroptosis Analysis
# ============================================================
set.seed(123)

base_dir <- normalizePath(file.path(".", "Breast Cancer 5 Necroptosis"))
data_raw_dir <- normalizePath(file.path(".", "Datasets", "data"))
output_dir <- base_dir

cat("Base directory:", base_dir, "\n")
cat("Data directory:", data_raw_dir, "\n")

dir.create(file.path(output_dir, "config"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "scripts"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "data_processed"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "metadata"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "results", "preprocessing"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "results", "batch_correction"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "results", "differential_expression"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "results", "necroptosis_genes"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "results", "necroptosis_score"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "results", "enrichment"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "results", "networks"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "results", "machine_learning"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "results", "consensus_genes"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "results", "predictive_model"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "results", "SHAP"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "results", "immune_infiltration"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "results", "external_validation"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "results", "tables"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "results", "figures"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "model_objects"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "logs"), showWarnings=FALSE, recursive=TRUE)
dir.create(file.path(output_dir, "report"), showWarnings=FALSE, recursive=TRUE)

pkg_required <- c(
  "readxl", "Biobase", "limma", "sva", "GSVA", "glmnet", "e1071",
  "randomForest", "pheatmap", "ggplot2", "ComplexHeatmap", "circlize",
  "igraph", "enrichR", "clusterProfiler", "org.Hs.eg.db", "pathview",
  "xCell", "caret", "ROCR", "pROC", "VennDiagram", "UpSetR",
  "corrplot", "reshape2", "dplyr", "tidyr", "RColorBrewer",
  "gridExtra", "ggrepel", "plotly"
)

for(pkg in pkg_required) {
  if(!require(pkg, quietly=TRUE, character.only=TRUE)) {
    cat("Installing:", pkg, "\n")
    if(pkg %in% c("Biobase", "limma", "sva", "GSVA", "ComplexHeatmap",
                   "clusterProfiler", "org.Hs.eg.db", "pathview", "xCell")) {
      BiocManager::install(pkg, ask=FALSE, update=FALSE)
    } else {
      install.packages(pkg, repos="https://cloud.r-project.org")
    }
    library(pkg, character.only=TRUE)
  }
}

if(!require("MCPcounter", quietly=TRUE)) {
  try(devtools::install_github("ebecht/MCPcounter", ref="master", subdir="MCPcounter"), silent=TRUE)
}

if(!require("SHAPforxgboost", quietly=TRUE)) {
  try(devtools::install_github("liuyanguu/SHAPforxgboost"), silent=TRUE)
}

if(!require("fastshap", quietly=TRUE)) {
  try(install.packages("fastshap", repos=c("https://bgreenwell.r-universe.dev", "https://cloud.r-project.org")), silent=TRUE)
}

# Helper function for saving TIFF figures
save_tiff_600dpi <- function(plot, filename, width = 7, height = 6) {
  ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white",
    device = "tiff"
  )
}

cat("Setup complete.\n")
