# ============================================================
# 05 - Discovery-cohort differential expression (GSE70947)
# Paired tumor-vs-normal limma analysis (148 verified pairs).
# The validation cohort is NOT used here.
# Feature universe for modeling = genes measured on both platforms.
# ============================================================
PROJECT_DIR <- "F:/independent projects/cancer/breast cancer/Breast cancer 8"
source(file.path(PROJECT_DIR, "config/config.R"))
source(file.path(PROJECT_DIR, "R/utils.R"))
suppressMessages({
  library(limma); library(ggplot2); library(dplyr); library(tidyr)
  library(ComplexHeatmap); library(circlize)
})

logmsg <- function(...) cat(sprintf(...), "\n")
logmsg("=== 05 Discovery differential expression (paired) ===")

disc  <- readRDS(file.path(DIR_DISC_PROC, "discovery_processed.rds"))
common <- readRDS(file.path(DIR_MODELS, "common_genes.rds"))

expr  <- disc$expr[common, , drop = FALSE]
meta  <- disc$meta
expr  <- expr[, meta$sample_id]
meta  <- meta[match(colnames(expr), meta$sample_id), ]

group   <- factor(meta$class, levels = CLASS_LEVELS)
patient <- factor(meta$patient_id)
stopifnot(all(sort(unique(patient)) == unique(patient)))

# ---------- paired design via duplicateCorrelation ----------
design <- model.matrix(~ 0 + group)
colnames(design) <- gsub("group", "", colnames(design))
dupcor <- duplicateCorrelation(expr, design = design, block = patient)
logmsg("Consensus within-patient correlation: %.4f", dupcor$consensus)

fit <- lmFit(expr, design = design, block = patient, correlation = dupcor$consensus)
cont <- makeContrasts("Tumor-Normal", levels = design)
fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2)

tt <- topTable(fit2, number = Inf, sort.by = "P")
tt$gene <- rownames(tt)
tt <- tt[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
names(tt) <- c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Value", "B")
tt <- tt[order(tt$P.Value), ]

# ---------- unpaired comparison (for sensitivity) ----------
design_u <- model.matrix(~ group)
fit_u <- lmFit(expr, design = design_u)
fit_u <- eBayes(fit_u)
tt_u <- topTable(fit_u, number = Inf, sort.by = "P", coef = 2)
tt_u$gene <- rownames(tt_u)
tt_u <- tt_u[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
names(tt_u) <- c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Value", "B")
tt_u <- tt_u[order(tt_u$P.Value), ]

# ---------- significant DEGs ----------
sig <- tt[tt$adj.P.Value < DEG_PVAL & abs(tt$logFC) >= DEG_LOGFC, ]
up   <- sig[sig$logFC > 0, ]
down <- sig[sig$logFC < 0, ]
logmsg("Significant DEGs (adj p<%.3f, |log2FC|>=%.1f): %d (up %d / down %d)",
       DEG_PVAL, DEG_LOGFC, nrow(sig), nrow(up), nrow(down))

# ---------- tables ----------
write.csv(tt, file.path(DIR_TABLES, "Table_08_Complete_DEG_Results.csv"), row.names = FALSE)
write.csv(tt_u, file.path(DIR_TABLES, "Table_08b_Complete_DEG_Results_Unpaired.csv"), row.names = FALSE)
write.csv(sig, file.path(DIR_TABLES, "Table_09_Significant_DEGs.csv"), row.names = FALSE)
write.csv(data.frame(rank = seq_len(nrow(tt)), tt), file.path(DIR_TABLES, "Table_09b_Ranked_Gene_List.csv"), row.names = FALSE)

deg <- list(paired = tt, unpaired = tt_u, significant = sig, up = up, down = down,
            dupcor = dupcor$consensus)
saveRDS(deg, file.path(DIR_MODELS, "deg_results.rds"))

# ---------- volcano ----------
vol <- tt
vol$sig <- ifelse(vol$adj.P.Value < DEG_PVAL & abs(vol$logFC) >= DEG_LOGFC,
                  ifelse(vol$logFC > 0, "Up", "Down"), "NS")
vol$sig <- factor(vol$sig, levels = c("Up", "Down", "NS"))
vol$neglogp <- -log10(vol$P.Value)

top_lab <- vol[vol$sig %in% c("Up", "Down"), ][order(vol[vol$sig %in% c("Up", "Down"), "P.Value"]), ][1:15, ]
vol_plot <- ggplot(vol, aes(x = logFC, y = neglogp, colour = sig)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey40") +
  geom_vline(xintercept = c(-DEG_LOGFC, DEG_LOGFC), linetype = 2, colour = "grey40") +
  ggrepel::geom_text_repel(data = top_lab, aes(label = gene), size = 2.2,
                           max.overlaps = 20, segment.size = 0.2,
                           show.legend = FALSE) +
  scale_colour_manual(values = c(Up = "#D55E00", Down = "#0072B2", NS = "grey75")) +
  theme_bc() +
  labs(x = "log2 fold-change (Tumor vs Normal)", y = "-log10 P-value",
       title = "Discovery differential expression (GSE70947)",
       subtitle = sprintf("Paired limma; %d significant DEGs (up %d / down %d)",
                          nrow(sig), nrow(up), nrow(down))) +
  theme(legend.position = "bottom")
save_tiff_600dpi(vol_plot, "Figure_04_Discovery_Volcano.tiff", width = 7, height = 6)

# ---------- MA plot ----------
ma <- tt
ma$sig <- ifelse(ma$adj.P.Value < DEG_PVAL & abs(ma$logFC) >= DEG_LOGFC,
                 ifelse(ma$logFC > 0, "Up", "Down"), "NS")
ma$sig <- factor(ma$sig, levels = c("Up", "Down", "NS"))
ma_plot <- ggplot(ma, aes(x = AveExpr, y = logFC, colour = sig)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "black") +
  scale_colour_manual(values = c(Up = "#D55E00", Down = "#0072B2", NS = "grey75")) +
  theme_bc() +
  labs(x = "Average expression", y = "log2 fold-change",
       title = "MA plot (discovery)", colour = "DEG") +
  theme(legend.position = "bottom")
save_tiff_600dpi(ma_plot, "Figure_04_Discovery_MA.tiff", width = 7, height = 6)

# ---------- top-50 DEG heatmap ----------
top50 <- sig[order(sig$P.Value), ][1:50, "gene"]
zmat <- t(scale(t(expr[top50, ])))
zmat[!is.finite(zmat)] <- 0
ann_col <- data.frame(Class = meta$class, row.names = colnames(zmat))
col_ann <- list(Class = c(Tumor = COL_TUMOR, Normal = COL_NORMAL))
open_tiff_base("Figure_04_Discovery_DEG_Heatmap.tiff", width = 8, height = 9)
h <- ComplexHeatmap::Heatmap(
  zmat, name = "z-score",
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B")),
  top_annotation = HeatmapAnnotation(df = ann_col, col = col_ann,
                                     annotation_name_side = "left"),
  show_column_names = FALSE, show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_split = meta$class,
  cluster_column_slices = FALSE,
  use_raster = TRUE, raster_quality = 5,
  heatmap_legend_param = list(title_gp = gpar(fontsize = 9), labels_gp = gpar(fontsize = 8)))
draw(h, annotation_legend_side = "bottom")
dev.off()

# ---------- top-DEG expression boxplots ----------
topg <- c(head(up$gene, 4), head(down$gene, 4))
box_long <- as.data.frame(t(expr[topg, , drop = FALSE]))
box_long$sample_id <- rownames(box_long)
box_dat <- tidyr::pivot_longer(box_long, -sample_id, names_to = "gene", values_to = "value") |>
  dplyr::mutate(class = meta$class[match(sample_id, meta$sample_id)])
box_dat$gene <- factor(box_dat$gene, levels = topg)
bp <- ggplot(box_dat, aes(x = gene, y = value, fill = class)) +
  geom_boxplot(outlier.size = 0.4, lwd = 0.3) +
  scale_fill_manual(values = c(Tumor = COL_TUMOR, Normal = COL_NORMAL)) +
  theme_bc(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "log2 expression",
       title = "Top differential genes: expression by class") +
  theme(legend.position = "bottom")
save_tiff_600dpi(bp, "Figure_04_Discovery_TopDEG_Boxplots.tiff", width = 10, height = 6)

# ---------- up/down counts ----------
count_df <- data.frame(direction = c("Upregulated", "Downregulated"),
                       n = c(nrow(up), nrow(down)))
cp <- ggplot(count_df, aes(x = direction, y = n, fill = direction)) +
  geom_col(colour = "black", width = 0.6) +
  geom_text(aes(label = n), vjust = -0.4, size = 5) +
  scale_fill_manual(values = c(Upregulated = COL_TUMOR, Downregulated = COL_NORMAL)) +
  theme_bc() +
  labs(x = "", y = "Number of genes",
       title = "Significant DEG counts (discovery)") +
  theme(legend.position = "none")
save_tiff_600dpi(cp, "Figure_04_Discovery_DEG_Counts.tiff", width = 6, height = 6)

logmsg("Saved DEG tables and figures.")
