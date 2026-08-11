###############################################################################
# SCRIPT: Generate Figure 3 - Hub Gene Expression, Volcano, and Meta-logFC
###############################################################################
OUT_DIR <- "F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)"
FIG_DIR <- file.path(OUT_DIR, "figures/hubs")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

library(grid)
library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(circlize)
library(RColorBrewer)

# Load data
load(file.path(OUT_DIR, "RData/01_processed_data_full.RData"))
load(file.path(OUT_DIR, "RData/05_hub_genes.RData"))
load(file.path(OUT_DIR, "RData/02_DE_analysis.RData"))

# 50 final hub genes
hubs <- final_hubs  # character vector

# Discovery common genes for subsetting
disc_genes <- common_genes_discovery

# ---------------------------------------------------------------------------
# Figure 3A: Heatmap of 50 hub gene expression across 3 datasets
# ---------------------------------------------------------------------------
build_hub_heatmap <- function() {
  datasets <- list(GSE42568 = processed_42568, GSE15852 = processed_15852, GSE29431 = processed_29431)
  
  # Extract expression for hub genes, z-score per dataset
  combined_expr <- NULL
  combined_meta <- NULL
  
  for (nm in names(datasets)) {
    full_expr <- datasets[[nm]]$expr
    # Subset to discovery common genes first, then hub genes
    common_idx <- intersect(disc_genes, rownames(full_expr))
    expr_common <- full_expr[common_idx, , drop = FALSE]
    meta <- datasets[[nm]]$meta
    present_hubs <- intersect(hubs, rownames(expr_common))
    sub <- expr_common[present_hubs, , drop = FALSE]
    sub <- t(scale(t(sub)))
    colnames(sub) <- paste0(nm, "_", colnames(sub))
    combined_expr <- cbind(combined_expr, sub)
    combined_meta <- rbind(combined_meta, data.frame(
      sample = colnames(sub),
      dataset = nm,
      condition = meta$class,
      stringsAsFactors = FALSE
    ))
  }
  
  combined_expr[is.na(combined_expr)] <- 0
  
  # Column annotations
  col_ha <- HeatmapAnnotation(
    Dataset = combined_meta$dataset,
    Condition = combined_meta$condition,
    col = list(
      Dataset = c(GSE42568 = "#E41A1C", GSE15852 = "#377EB8", GSE29431 = "#4DAF4A"),
      Condition = c(Normal = "blue", Tumor = "red")
    ),
    annotation_name_gp = gpar(fontsize = 10)
  )
  
  # Row annotation: mean logFC direction
  hub_mean <- rowMeans(cbind(
    hub_summary$logFC_42568[match(rownames(combined_expr), hub_summary$Gene)],
    hub_summary$logFC_15852[match(rownames(combined_expr), hub_summary$Gene)],
    hub_summary$logFC_29431[match(rownames(combined_expr), hub_summary$Gene)]
  ), na.rm = TRUE)
  
  row_ha <- rowAnnotation(
    Direction = ifelse(hub_mean > 0, "Up", "Down"),
    col = list(Direction = c(Up = "red", Down = "blue")),
    annotation_name_gp = gpar(fontsize = 10)
  )
  
  # Color scale
  col_fun <- colorRamp2(c(-2, 0, 2), c("navy", "white", "firebrick3"))
  
  ht <- Heatmap(
    combined_expr,
    name = "Z-score",
    col = col_fun,
    top_annotation = col_ha,
    right_annotation = row_ha,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 6),
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    column_split = combined_meta$dataset,
    column_title_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(direction = "vertical"),
    row_dend_reorder = TRUE,
    column_dend_reorder = TRUE
  )
  
  tiff(file.path(FIG_DIR, "Figure3A_hub_heatmap.tiff"),
       width = 14, height = 10, units = "in", res = 600, compression = "lzw")
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  
  # Also save PNG for quick view
  png(file.path(FIG_DIR, "Figure3A_hub_heatmap.png"),
      width = 14, height = 10, units = "in", res = 300)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

# ---------------------------------------------------------------------------
# Figure 3B: Combined volcano plot highlighting hub genes
# ---------------------------------------------------------------------------
build_combined_volcano <- function() {
  de_list <- list(GSE42568 = de_42568$results, 
                  GSE15852 = de_15852$results, 
                  GSE29431 = de_29431$results)
  
  volcano_data <- bind_rows(lapply(names(de_list), function(nm) {
    df <- de_list[[nm]]
    df$Dataset <- nm
    df$is_hub <- df$Gene %in% hubs
    df
  }))
  
  volcano_data$Significance <- "NS"
  volcano_data$Significance[volcano_data$adj.P.Val < 0.05 & abs(volcano_data$logFC) > 1] <- "DEG"
  volcano_data$Significance[volcano_data$is_hub] <- "Hub"
  volcano_data$Significance <- factor(volcano_data$Significance, levels = c("NS", "DEG", "Hub"))
  
  p <- ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = Significance, alpha = Significance, size = Significance)) +
    scale_color_manual(values = c(NS = "grey75", DEG = "grey40", Hub = "red")) +
    scale_alpha_manual(values = c(NS = 0.3, DEG = 0.3, Hub = 1)) +
    scale_size_manual(values = c(NS = 1, DEG = 1, Hub = 2.5)) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.3) +
    facet_wrap(~ Dataset, nrow = 1) +
    labs(x = "log2(Fold Change)", y = "-log10(Adjusted P-value)",
         title = "Differential Expression with Hub Genes Highlighted") +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey95"),
      panel.grid.minor = element_blank()
    ) +
    guides(
      color = guide_legend(override.aes = list(size = 3, alpha = 1))
    )
  
  # Label hub genes
  hub_labels <- volcano_data %>% filter(is_hub) %>%
    group_by(Gene) %>% slice_max(abs(logFC), n = 1) %>% ungroup()
  
  p <- p + geom_text_repel(
    data = hub_labels,
    aes(label = Gene),
    size = 2.5,
    max.overlaps = 30,
    box.padding = 0.3,
    segment.color = "grey30",
    segment.size = 0.3,
    min.segment.length = 0
  )
  
  ggsave(file.path(FIG_DIR, "Figure3B_combined_volcano.tiff"),
         p, width = 16, height = 8, dpi = 600, compression = "lzw")
  ggsave(file.path(FIG_DIR, "Figure3B_combined_volcano.png"),
         p, width = 16, height = 8, dpi = 300)
}

# ---------------------------------------------------------------------------
# Figure 3C: Bar plot of hub gene meta-logFC
# ---------------------------------------------------------------------------
build_meta_logFC_barplot <- function() {
  # Compute mean logFC from the 3 datasets
  plot_data <- hub_summary %>%
    mutate(
      mean_logFC = (logFC_42568 + logFC_15852 + logFC_29431) / 3,
      Direction = ifelse(mean_logFC > 0, "Upregulated", "Downregulated")
    ) %>%
    arrange(mean_logFC) %>%
    mutate(Gene = factor(Gene, levels = Gene))
  
  p <- ggplot(plot_data, aes(x = Gene, y = mean_logFC, fill = Direction)) +
    geom_bar(stat = "identity", width = 0.85) +
    scale_fill_manual(values = c(Downregulated = "steelblue", Upregulated = "tomato")) +
    labs(
      x = "Hub Gene",
      y = "Mean log2(Fold Change)",
      title = "Meta-Analysis of Hub Gene Expression",
      subtitle = paste0(length(hubs), " multi-cohort consensus hub genes (", 
                        sum(plot_data$Direction == "Downregulated"), " down, ",
                        sum(plot_data$Direction == "Upregulated"), " up)")
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
      axis.title.x = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "top"
    ) +
    geom_hline(yintercept = 0, linewidth = 0.5)
  
  ggsave(file.path(FIG_DIR, "Figure3C_meta_logFC_barplot.tiff"),
         p, width = 14, height = 7, dpi = 600, compression = "lzw")
  ggsave(file.path(FIG_DIR, "Figure3C_meta_logFC_barplot.png"),
         p, width = 14, height = 7, dpi = 300)
}

# Run all
cat("Building Figure 3A: Hub gene expression heatmap...\n")
build_hub_heatmap()

cat("Building Figure 3B: Combined volcano plot...\n")
build_combined_volcano()

cat("Building Figure 3C: Meta-logFC bar plot...\n")
build_meta_logFC_barplot()

cat("All Figure 3 panels generated in:", FIG_DIR, "\n")
