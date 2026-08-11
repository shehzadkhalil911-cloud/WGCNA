# 00_check_environment.R
suppressMessages({
  library(sessioninfo)
})
source("config/analysis_config.R")
source("R/load_functions.R")
ensure_paths()
cat("=== Environment check ===\n")
cat("R version:", R.version.string, "\n")
cat("OUTPUT_DIR:", OUTPUT_DIR, "\n")
cat("Datasets exist:", dir.exists(DATASETS_DIR), "\n")
pkgs_expected <- c("GEOquery","Biobase","limma","sva","preprocessCore","AnnotationDbi","annotate",
  "org.Hs.eg.db","tidyverse","data.table","janitor","stringr","readxl","openxlsx","caret","glmnet",
  "ranger","e1071","kernlab","pROC","PRROC","yardstick","ROSE","recipes","doParallel","doRNG","foreach",
  "clusterProfiler","enrichplot","DOSE","ReactomePA","msigdbr","pathview","KEGGREST","STRINGdb","igraph",
  "ggraph","ComplexHeatmap","circlize","EnhancedVolcano","pheatmap","ggplot2","ggpubr","ggrepel",
  "patchwork","cowplot","ggvenn","ComplexUpset","shapviz","DALEX","ingredients","immunedeconv","GSVA",
  "GSEABase","xCell","corrplot","Hmisc","psych","rstatix","boot","officer","flextable","quarto","rmarkdown","renv")
inst <- rownames(installed.packages())
missing <- setdiff(pkgs_expected, inst)
cat("Missing expected packages:", if (length(missing) == 0) "NONE" else paste(missing, collapse = ", "), "\n")
cat("\nDocumented package substitutions (R 4.6.0 / Bioc 3.23):\n")
for (k in names(PKG_SUBSTITUTIONS)) cat(" ", k, "->", PKG_SUBSTITUTIONS[[k]], "\n")
saveRDS(list(r_version = R.version.string, missing = missing,
             substitutions = PKG_SUBSTITUTIONS, timestamp = Sys.time()),
        file.path(paths$res_logs, "environment_check.rds"))
cat("\nEnvironment check complete.\n")
quit(save = "no", status = 0, runLast = FALSE)
