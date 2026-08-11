###############################################################################
# SCRIPT: Combined immune infiltration box plots (Figure 6A-C)
###############################################################################
OUT_DIR <- "F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)"
FIG_DIR <- file.path(OUT_DIR, "figures/immune")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

load(file.path(OUT_DIR, "RData/08_immune_infiltration.RData"))
load(file.path(OUT_DIR, "RData/01_processed_data_full.RData"))

# Select top biologically meaningful cell types per dataset
select_cells <- list(
  GSE42568 = c("Adipocytes", "Astrocytes", "Epithelial cells", 
               "Pericytes", "Th1 cells", "StromaScore"),
  GSE15852 = c("Adipocytes", "Epithelial cells", "Hepatocytes",
               "Neutrophils", "mv Endothelial cells", "StromaScore"),
  GSE29431 = c("Adipocytes", "Th2 cells", "Macrophages M1",
               "mv Endothelial cells", "ly Endothelial cells", "CD4+ T-cells")
)

full_expr <- list(
  GSE42568 = processed_42568,
  GSE15852 = processed_15852,
  GSE29431 = processed_29431
)

plot_list <- list()

for (ds in names(select_cells)) {
  scores <- immune_results[[ds]][["scores"]]
  meta <- full_expr[[ds]][["meta"]]
  cells_present <- intersect(select_cells[[ds]], rownames(scores))
  
  sub <- scores[cells_present, , drop = FALSE]
  
  # Build data frame: rows=samples, cols=cell types
  df <- as.data.frame(t(sub))
  df$sample_id <- rownames(df)
  df <- merge(df, meta[, c("sample_id", "class")], by = "sample_id", all.x = TRUE)
  
  long <- df %>%
    pivot_longer(cols = all_of(cells_present), names_to = "CellType", values_to = "Score") %>%
    mutate(CellType = factor(CellType, levels = cells_present))
  
  # Get FDR labels from diff table
  diff_df <- immune_results[[ds]][["diff"]]
  
  p <- ggplot(long, aes(x = class, y = Score, fill = class)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.7, width = 0.5) +
    facet_wrap(~ CellType, nrow = 2, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = c(Normal = "steelblue", Tumor = "tomato")) +
    labs(
      x = NULL, y = "xCell Score",
      title = ds
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      panel.grid.major.x = element_blank(),
      plot.title = element_text(face = "bold", size = 12)
    )
  
  plot_list[[ds]] <- p
}

# Combine into one figure
combined <- grid.arrange(
  plot_list[["GSE42568"]], plot_list[["GSE15852"]], plot_list[["GSE29431"]],
  nrow = 3,
  top = "Immune Infiltration Landscape in Breast Cancer"
)

ggsave(file.path(FIG_DIR, "Figure6_ABC_immune_boxplots.tiff"),
       combined, width = 14, height = 14, dpi = 600, compression = "lzw")
ggsave(file.path(FIG_DIR, "Figure6_ABC_immune_boxplots.png"),
       combined, width = 14, height = 14, dpi = 300)

cat("Figure 6A-C saved to:", FIG_DIR, "\n")
