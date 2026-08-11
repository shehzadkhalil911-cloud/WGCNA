load_libraries <- function() {
  library(ggplot2)
  library(limma)
  library(Biobase)
  library(sva)
  library(WGCNA)
  library(pheatmap)
  library(reshape2)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(writexl)
  library(openxlsx)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(enrichplot)
  library(ReactomePA)
  library(GSVA)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(corrplot)
  library(igraph)
  library(stringr)
  library(gridExtra)
  library(cowplot)
  library(ggrepel)
  library(VennDiagram)
  library(scales)
}

save_tiff_600dpi <- function(plot, filename, width = 7, height = 6) {
  ggplot2::ggsave(
    filename = filename, plot = plot, width = width, height = height,
    units = "in", dpi = 600, compression = "lzw", bg = "white", device = "tiff"
  )
}

open_tiff <- function(filename, width = 7, height = 6) {
  grDevices::tiff(filename, width = width, height = height,
                  units = "in", res = 600, compression = "lzw", bg = "white")
}

log_msg <- function(msg, logfile = "logs/analysis.log") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  write(sprintf("[%s] %s", timestamp, msg), file = logfile, append = TRUE)
  cat(sprintf("[%s] %s\n", timestamp, msg))
}

log_error <- function(msg, logfile = "logs/error.log") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  write(sprintf("[%s] ERROR: %s", timestamp, msg), file = logfile, append = TRUE)
  cat(sprintf("[%s] ERROR: %s\n", timestamp, msg))
}

.get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
