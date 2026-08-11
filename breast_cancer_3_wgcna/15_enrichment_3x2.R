###############################################################################
# SCRIPT: Combined enrichment bar plot (3x2: BP, CC, MF, KEGG, Reactome, Hallmark)
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
  "GO Biological Process" = enrich_hubs$go_bp,
  "GO Cellular Component" = enrich_hubs$go_cc,
  "GO Molecular Function" = enrich_hubs$go_mf,
  "KEGG Pathway" = enrich_hubs$kegg,
  "Reactome" = enrich_hubs$reactome,
  "Hallmark" = enrich_hubs$hallmark
)

plot_list <- list()

for (nm in names(enrich_list)) {
  obj <- enrich_list[[nm]]
  if (is.null(obj) || nrow(obj@result) == 0) next
  
  df <- obj@result %>% head(5) %>%
    mutate(
      neg_log10 = -log10(p.adjust),
      Description = gsub("^HALLMARK_", "", Description),
      Description = gsub("_", " ", Description)
    ) %>%
    arrange(p.adjust) %>%
    mutate(Description = factor(Description, levels = rev(Description)))
  
  p <- ggplot(df, aes(x = Count, y = Description, fill = neg_log10)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_gradient(low = "darkblue", high = "firebrick3", name = "-log10(FDR)") +
    labs(x = "Gene Count", y = NULL, title = nm) +
    theme_bw(base_size = 6) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.key.height = unit(0.25, "cm"),
      legend.text = element_text(size = 5),
      legend.title = element_text(size = 5, face = "bold"),
      plot.title = element_text(face = "bold", size = 7),
      axis.title.x = element_text(size = 6, face = "bold"),
      axis.text.y = element_text(face = "bold", size = 5),
      axis.text.x = element_text(size = 5)
    ) +
    geom_text(aes(label = Count), hjust = -0.3, size = 1.8, fontface = "bold")
  
  plot_list[[nm]] <- p
}

combined <- grid.arrange(
  plot_list[["GO Biological Process"]],
  plot_list[["GO Cellular Component"]],
  plot_list[["GO Molecular Function"]],
  plot_list[["KEGG Pathway"]],
  plot_list[["Reactome"]],
  plot_list[["Hallmark"]],
  nrow = 3, ncol = 2,
  top = grid::textGrob("Functional Enrichment of Hub Genes",
                        gp = grid::gpar(fontsize = 9, fontface = "bold"))
)

ggsave(file.path(FIG_DIR, "enrichment_combined_3x2_barplot.tiff"),
       combined, width = 7, height = 5, dpi = 1200, compression = "lzw")
ggsave(file.path(FIG_DIR, "enrichment_combined_3x2_barplot.png"),
       combined, width = 7, height = 5, dpi = 600)

ggsave(file.path(OUT_DIR, "figures/Figure4_enrichment_3x2.tiff"),
       combined, width = 7, height = 5, dpi = 1200, compression = "lzw")

cat("Combined 3x2 enrichment bar plot saved.\n")
