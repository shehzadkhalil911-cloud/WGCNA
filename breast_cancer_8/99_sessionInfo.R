suppressMessages({
  library(xgboost); library(limma); library(glmnet); library(ranger)
  library(clusterProfiler); library(shapviz); library(ComplexHeatmap)
  library(ggplot2); library(pROC); library(PRROC); library(sva)
  library(enrichplot); library(rmarkdown); library(knitr); library(magick)
})
con <- file("F:/independent projects/cancer/breast cancer/Breast cancer 8/report/sessionInfo.txt", open = "w")
writeLines(c(
  "Breast cancer 8 study - R session info",
  paste0("Generated: ", Sys.time()),
  "=======================================",
  ""
), con)
writeLines(capture.output(utils::sessionInfo()), con)
close(con)
cat("sessionInfo written.\n")
