###############################################################################
# SCRIPT: Generate all main figures in EMF format
###############################################################################
OUT_DIR <- "F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)"
FIG_DIR <- file.path(OUT_DIR, "figures")
library(devEMF)

# Load all required data
load(file.path(OUT_DIR, "RData/01_processed_data_full.RData"))
load(file.path(OUT_DIR, "RData/05_hub_genes.RData"))
load(file.path(OUT_DIR, "RData/02_DE_analysis.RData"))
load(file.path(OUT_DIR, "RData/07_validation.RData"))
load(file.path(OUT_DIR, "RData/08_immune_infiltration.RData"))

hubs <- final_hubs
disc_genes <- common_genes_discovery
validated_genes <- c("ACADL", "ACSM5", "ADH1B", "C10orf10",
                     "CACNA2D1", "CEBPA", "CHRDL1", "DLC1")

suppressPackageStartupMessages({
  library(grid)
  library(ComplexHeatmap)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(circlize)
  library(gridExtra)
})

# ---------------------------------------------------------------------------
# Figure 3B: Combined volcano plot (EMF)
# ---------------------------------------------------------------------------
cat("Figure 3B: Combined volcano...\n")
de_list <- list(GSE42568 = de_42568$results, 
                GSE15852 = de_15852$results, 
                GSE29431 = de_29431$results)

volcano_data <- bind_rows(lapply(names(de_list), function(nm) {
  df <- de_list[[nm]]; df$Dataset <- nm; df$is_hub <- df$Gene %in% hubs; df
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
  labs(x = "log2(Fold Change)", y = "-log10(Adjusted P-value)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "grey95"),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

hub_labels <- volcano_data %>% filter(is_hub) %>%
  group_by(Gene) %>% slice_max(abs(logFC), n = 1) %>% ungroup()
p <- p + geom_text_repel(data = hub_labels, aes(label = Gene), size = 2.5,
                         max.overlaps = 30, box.padding = 0.3,
                         segment.color = "grey30", segment.size = 0.3)

emf(file.path(FIG_DIR, "hubs/Figure3B_combined_volcano.emf"), width = 16, height = 8,
    emfPlus = FALSE, coordDPI = 600)
print(p); dev.off()

# ---------------------------------------------------------------------------
# Figure 3C: Meta-logFC bar plot (EMF)
# ---------------------------------------------------------------------------
cat("Figure 3C: Meta-logFC bar plot...\n")
plot_data <- hub_summary %>%
  mutate(mean_logFC = (logFC_42568 + logFC_15852 + logFC_29431) / 3,
         Direction = ifelse(mean_logFC > 0, "Upregulated", "Downregulated")) %>%
  arrange(mean_logFC) %>% mutate(Gene = factor(Gene, levels = Gene))

p <- ggplot(plot_data, aes(x = Gene, y = mean_logFC, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.85) +
  scale_fill_manual(values = c(Downregulated = "steelblue", Upregulated = "tomato")) +
  labs(x = "Hub Gene", y = "Mean log2(Fold Change)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.title.x = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "top") +
  geom_hline(yintercept = 0, linewidth = 0.5)

emf(file.path(FIG_DIR, "hubs/Figure3C_meta_logFC_barplot.emf"), width = 14, height = 7,
    emfPlus = FALSE, coordDPI = 600)
print(p); dev.off()

# ---------------------------------------------------------------------------
# Figure 5A: Validation boxplots (EMF)
# ---------------------------------------------------------------------------
cat("Figure 5A: Validation boxplots...\n")
expr <- processed_70947[["expr"]]
meta <- processed_70947[["meta"]]
sub <- expr[validated_genes, , drop = FALSE]
df <- as.data.frame(t(sub))
df$sample_id <- rownames(df)
df <- merge(df, meta[, c("sample_id", "class")], by = "sample_id")

long <- df %>% pivot_longer(cols = all_of(validated_genes), names_to = "Gene", values_to = "Expression") %>%
  mutate(Gene = factor(Gene, levels = validated_genes))

p <- ggplot(long, aes(x = class, y = Expression, fill = class)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7, width = 0.5) +
  facet_wrap(~ Gene, nrow = 2, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = c(Normal = "steelblue", Tumor = "tomato")) +
  labs(x = NULL, y = "Expression Level") +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank())

emf(file.path(FIG_DIR, "validation/Figure5A_validation_boxplots.emf"), width = 14, height = 10,
    emfPlus = FALSE, coordDPI = 600)
print(p); dev.off()

# ---------------------------------------------------------------------------
# Figure 6A-C: Immune boxplots (EMF)
# ---------------------------------------------------------------------------
cat("Figure 6A-C: Immune boxplots...\n")
select_cells <- list(
  GSE42568 = c("Adipocytes", "Astrocytes", "Epithelial cells", "Pericytes", "Th1 cells", "StromaScore"),
  GSE15852 = c("Adipocytes", "Epithelial cells", "Hepatocytes", "Neutrophils", "mv Endothelial cells", "StromaScore"),
  GSE29431 = c("Adipocytes", "Th2 cells", "Macrophages M1", "mv Endothelial cells", "ly Endothelial cells", "CD4+ T-cells")
)
full_expr <- list(GSE42568 = processed_42568, GSE15852 = processed_15852, GSE29431 = processed_29431)
plot_list <- list()

for (ds in names(select_cells)) {
  scores <- immune_results[[ds]][["scores"]]
  meta <- full_expr[[ds]][["meta"]]
  cells_present <- intersect(select_cells[[ds]], rownames(scores))
  sub_mat <- scores[cells_present, , drop = FALSE]
  df <- as.data.frame(t(sub_mat))
  df$sample_id <- rownames(df)
  df <- merge(df, meta[, c("sample_id", "class")], by = "sample_id")
  long <- df %>% pivot_longer(cols = all_of(cells_present), names_to = "CellType", values_to = "Score") %>%
    mutate(CellType = factor(CellType, levels = cells_present))
  
  p <- ggplot(long, aes(x = class, y = Score, fill = class)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.7, width = 0.5) +
    facet_wrap(~ CellType, nrow = 2, ncol = 3, scales = "free_y") +
    scale_fill_manual(values = c(Normal = "steelblue", Tumor = "tomato")) +
    labs(x = NULL, y = "xCell Score", title = ds) +
    theme_bw(base_size = 10) +
    theme(strip.background = element_rect(fill = "grey95"),
          strip.text = element_text(face = "bold", size = 9),
          legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          panel.grid.major.x = element_blank(), plot.title = element_text(face = "bold", size = 12))
  plot_list[[ds]] <- p
}

combined <- grid.arrange(plot_list[["GSE42568"]], plot_list[["GSE15852"]], plot_list[["GSE29431"]],
                          nrow = 3, top = "Immune Infiltration Landscape in Breast Cancer")

emf(file.path(FIG_DIR, "immune/Figure6_ABC_immune_boxplots.emf"), width = 14, height = 14,
    emfPlus = FALSE, coordDPI = 600)
grid.draw(combined); dev.off()

cat("All EMF figures generated.\n")
