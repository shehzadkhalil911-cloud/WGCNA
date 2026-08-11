###############################################################################
# SCRIPT: Combined enrichment bar plot (2x2: GO BP, GO CC, GO MF, KEGG)
###############################################################################
OUT_DIR <- "F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)"
FIG_DIR <- file.path(OUT_DIR, "figures/enrichment")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(gridExtra)

load(file.path(OUT_DIR, "RData/06_enrichment.RData"))

enrich_list <- list(
  BP = enrich_hubs$go_bp,
  CC = enrich_hubs$go_cc,
  MF = enrich_hubs$go_mf,
  KEGG = enrich_hubs$kegg
)

plot_list <- list()

for (nm in names(enrich_list)) {
  obj <- enrich_list[[nm]]
  if (is.null(obj) || nrow(obj@result) == 0) next
  
  df <- obj@result %>% head(5) %>%
    mutate(
      Description = factor(Description, levels = rev(Description)),
      neg_log10 = -log10(p.adjust)
    )
  
  p <- ggplot(df, aes(x = Count, y = Description, fill = neg_log10)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_gradient(low = "darkblue", high = "firebrick3", name = "-log10(FDR)") +
    labs(
      x = "Gene Count", y = NULL,
      title = nm
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.key.height = unit(0.4, "cm"),
      plot.title = element_text(face = "bold", size = 12),
      axis.text.y = element_text(size = 8)
    ) +
    geom_text(aes(label = Count), hjust = -0.3, size = 3)
  
  plot_list[[nm]] <- p
}

combined <- grid.arrange(
  plot_list[["BP"]], plot_list[["CC"]],
  plot_list[["MF"]], plot_list[["KEGG"]],
  nrow = 2, ncol = 2,
  top = "Functional Enrichment of Hub Genes"
)

ggsave(file.path(FIG_DIR, "enrichment_combined_2x2_barplot.tiff"),
       combined, width = 16, height = 10, dpi = 600, compression = "lzw")
ggsave(file.path(FIG_DIR, "enrichment_combined_2x2_barplot.png"),
       combined, width = 16, height = 10, dpi = 300)

# Also save to main figure directory
ggsave(file.path(OUT_DIR, "figures/Figure4_enrichment_2x2.tiff"),
       combined, width = 16, height = 10, dpi = 600, compression = "lzw")

cat("Combined enrichment bar plot saved.\n")
