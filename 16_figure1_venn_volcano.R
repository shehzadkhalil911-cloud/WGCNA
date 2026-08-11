###############################################################################
# SCRIPT: 2x2 figure - DEG Venn + Enhanced Volcano plots (3 datasets)
###############################################################################
OUT_DIR <- "F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)"
FIG_DIR <- file.path(OUT_DIR, "figures/DEA")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

library(ggplot2)
library(ggrepel)
library(ggvenn)
library(dplyr)
library(gridExtra)
library(grid)

load(file.path(OUT_DIR, "RData/02_DE_analysis.RData"))

# 50 hub genes for highlighting
load(file.path(OUT_DIR, "RData/05_hub_genes.RData"))
hubs <- final_hubs

# ---- DEG lists (p < 0.05) ----
deg_42568 <- de_42568$results %>% filter(adj.P.Val < 0.05) %>% pull(Gene)
deg_15852 <- de_15852$results %>% filter(adj.P.Val < 0.05) %>% pull(Gene)
deg_29431 <- de_29431$results %>% filter(adj.P.Val < 0.05) %>% pull(Gene)

# ---- Panel 1: Venn diagram ----
venn_data <- list(GSE42568 = deg_42568, GSE15852 = deg_15852, GSE29431 = deg_29431)
p_venn <- ggvenn(venn_data,
  fill_color = c("#E41A1C", "#377EB8", "#4DAF4A"),
  stroke_size = 0.3, set_name_size = 3, text_size = 2.2
) +
  labs(title = "DEG Overlap Across Datasets") +
  theme(plot.title = element_text(face = "bold", size = 8, hjust = 0.5))

# ---- Panels 2-4: Enhanced volcano plots ----
de_list <- list(
  GSE42568 = de_42568$results,
  GSE15852 = de_15852$results,
  GSE29431 = de_29431$results
)

names_pretty <- c(GSE42568 = "GSE42568", GSE15852 = "GSE15852", GSE29431 = "GSE29431")

volcano_list <- list()

for (nm in names(de_list)) {
  df <- de_list[[nm]]
  df$sig <- "NS"
  df$sig[df$adj.P.Val < 0.05 & df$logFC > 1] <- "Up"
  df$sig[df$adj.P.Val < 0.05 & df$logFC < -1] <- "Down"
  df$sig <- factor(df$sig, levels = c("Up", "Down", "NS"))
  df$neg_log10 <- -log10(df$adj.P.Val)
  
  p <- ggplot(df, aes(x = logFC, y = neg_log10)) +
    geom_point(aes(color = sig, alpha = sig, size = sig)) +
    scale_color_manual(values = c(Up = "red", Down = "blue", NS = "grey75"),
                       labels = c("Upregulated", "Downregulated", "Not significant")) +
    scale_alpha_manual(values = c(Up = 0.6, Down = 0.6, NS = 0.25),
                       labels = c("Upregulated", "Downregulated", "Not significant")) +
    scale_size_manual(values = c(Up = 0.6, Down = 0.6, NS = 0.4),
                      labels = c("Upregulated", "Downregulated", "Not significant")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.3, linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.3, linewidth = 0.3) +
    labs(
      x = "log2(Fold Change)", y = "-log10(Adj. P)",
      title = names_pretty[nm],
      color = NULL, alpha = NULL, size = NULL
    ) +
    theme_bw(base_size = 6) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 5),
      legend.key.size = unit(0.25, "cm"),
      plot.title = element_text(face = "bold", size = 7, hjust = 0.5),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 5.5),
      panel.grid.minor = element_blank()
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 2, alpha = 1), nrow = 1),
      alpha = "none", size = "none"
    ) +
    xlim(-8, 8)
  
  volcano_list[[nm]] <- p
}

# ---- Combine 2x2 ----
combined <- grid.arrange(
  p_venn, volcano_list[["GSE42568"]],
  volcano_list[["GSE15852"]], volcano_list[["GSE29431"]],
  nrow = 2, ncol = 2
)

ggsave(file.path(FIG_DIR, "Figure1_venn_volcano_2x2.tiff"),
       combined, width = 8, height = 6, dpi = 1200, compression = "lzw")
ggsave(file.path(FIG_DIR, "Figure1_venn_volcano_2x2.png"),
       combined, width = 8, height = 6, dpi = 600)

cat("Figure 1 (Venn + Volcanoes) saved.\n")
