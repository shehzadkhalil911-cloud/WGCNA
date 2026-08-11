###############################################################################
# SCRIPT: Combined boxplot of 8 validated hub genes in GSE70947
###############################################################################
OUT_DIR <- "F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)"
FIG_DIR <- file.path(OUT_DIR, "figures/validation")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

library(ggplot2)
library(dplyr)
library(tidyr)

load(file.path(OUT_DIR, "RData/01_processed_data_full.RData"))
load(file.path(OUT_DIR, "RData/07_validation.RData"))

validated_genes <- c("ACADL", "ACSM5", "ADH1B", "C10orf10",
                     "CACNA2D1", "CEBPA", "CHRDL1", "DLC1")

expr <- processed_70947[["expr"]]
meta <- processed_70947[["meta"]]

# Extract expression for validated genes and reshape
sub <- expr[validated_genes, , drop = FALSE]
df <- as.data.frame(t(sub))
df$sample <- rownames(df)
df <- merge(df, meta[, c("sample_id", "class")], by.x = "sample", by.y = "sample_id")

long <- df %>%
  pivot_longer(cols = all_of(validated_genes), names_to = "Gene", values_to = "Expression") %>%
  mutate(Gene = factor(Gene, levels = validated_genes))

# Compute stats for significance annotation
stats <- long %>%
  group_by(Gene) %>%
  summarise(
    p_label = paste0("FDR = ", format(
      val_results$FDR[match(unique(Gene), val_results$Gene)], digits = 3, scientific = TRUE
    )),
    y_max = max(Expression)
  )

p <- ggplot(long, aes(x = class, y = Expression, fill = class)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7, width = 0.5) +
  facet_wrap(~ Gene, nrow = 2, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = c(Normal = "steelblue", Tumor = "tomato")) +
  labs(
    x = NULL, y = "Expression Level",
    title = "Validation of Hub Genes in GSE70947 (148 Normal, 148 Tumor)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )

# Add significance annotations
p <- p + geom_text(
  data = stats,
  aes(x = 1.5, y = y_max, label = p_label),
  size = 2.8, vjust = -0.5, inherit.aes = FALSE
)

ggsave(file.path(FIG_DIR, "Figure5A_validation_boxplots.tiff"),
       p, width = 14, height = 10, dpi = 600, compression = "lzw")
ggsave(file.path(FIG_DIR, "Figure5A_validation_boxplots.png"),
       p, width = 14, height = 10, dpi = 300)

cat("Validation boxplots saved to:", FIG_DIR, "\n")
