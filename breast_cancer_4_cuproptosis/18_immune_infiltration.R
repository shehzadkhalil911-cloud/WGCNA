# 18_immune_infiltration.R
source("config/analysis_config.R")
source("R/load_functions.R")
ensure_paths()
suppressMessages({ library(ggplot2); library(ggpubr); library(pheatmap); library(reshape2); library(rstatix) })

combat <- readRDS(file.path(paths$discovery, "discovery_combat.rds"))
meta <- readRDS(file.path(paths$discovery, "discovery_cohort_raw.rds"))$metadata
group <- meta$class

## Primary: xCell
cat("Running xCell...\n"); flush.console()
xcell <- run_xcell(combat)
saveRDS(xcell, file.path(paths$res_immune, "xcell_scores.rds"))
## Sensitivity: ssGSEA immune
cat("Running ssGSEA immune...\n"); flush.console()
ssg <- run_ssgsea_immune(combat)
saveRDS(ssg, file.path(paths$res_immune, "ssgsea_immune.rds"))
## Manual ESTIMATE
est <- manual_estimate_score(combat)
saveRDS(est, file.path(paths$res_immune, "estimate_scores.rds"))

## combine primary (xCell) with ESTIMATE for differential testing
immune_primary <- as.data.frame(t(xcell))   # samples x cells
immune_primary$ESTIMATE_immune_score <- est$ESTIMATE_immune_score[rownames(immune_primary)]
immune_primary$ESTIMATE_stromal_score <- est$ESTIMATE_stromal_score[rownames(immune_primary)]
write.csv(immune_primary, file.path(paths$res_immune, "immune_infiltration_primary.csv"))

## differential immune cell abundance (tumor vs normal)
cells <- setdiff(colnames(immune_primary), c("ESTIMATE_immune_score","ESTIMATE_stromal_score"))
res_rows <- list()
for (c in cells) {
  vals <- immune_primary[[c]]
  tt <- wilcox.test(vals[group == "Tumor"], vals[group == "Normal"])
  res_rows[[c]] <- data.frame(cell = c, mean_tumor = mean(vals[group=="Tumor"]),
    mean_normal = mean(vals[group=="Normal"]), p_value = tt$p.value,
    log2FC = log2(mean(vals[group=="Tumor"]+1e-6)/mean(vals[group=="Normal"]+1e-6)))
}
diff_immune <- do.call(rbind, res_rows)
diff_immune$adj_p <- p.adjust(diff_immune$p_value, method = "BH")
write.csv(diff_immune, file.path(paths$res_immune, "differential_immune_cells.csv"), row.names = FALSE)

## box plots
samp_names <- rownames(immune_primary)
group_ord <- group[match(samp_names, colnames(combat))]
idf <- cbind(sample = samp_names, immune_primary[samp_names, cells, drop = FALSE])
idf <- reshape2::melt(idf, id.vars = "sample")
colnames(idf) <- c("sample","cell","score")
idf$group <- group_ord[match(idf$sample, samp_names)]
idf$group <- factor(idf$group, levels = c("Normal","Tumor"))
p_box <- ggplot(idf, aes(x = cell, y = score, fill = group)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title = "Immune cell abundance (xCell)") + pub_theme()
save_publication_plot(p_box, file.path(paths$res_immune, "immune_boxplots"), width = 10, height = 6)
## heatmap (ComplexHeatmap; pheatmap annotation broken on R4.6)
suppressMessages(library(ComplexHeatmap))
mat <- t(scale(t(immune_primary[samp_names, cells, drop = FALSE])))
colnames(mat) <- cells; rownames(mat) <- samp_names
ha <- ComplexHeatmap::rowAnnotation(group = group_ord,
  col = list(group = c(Normal = "#4E79A7", Tumor = "#E15759")))
ht <- ComplexHeatmap::Heatmap(mat, name = "z-score", left_annotation = ha,
  show_row_names = FALSE, show_column_names = TRUE, column_names_gp = grid::gpar(fontsize = 8))
grDevices::png(file.path(paths$res_immune, "immune_heatmap.png"), width = 1200, height = 1400, res = 150)
ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
grDevices::dev.off()
grDevices::tiff(file.path(paths$res_immune, "immune_heatmap.tiff"), width = 8, height = 9,
                units = "in", res = 600, compression = "lzw")
ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
grDevices::dev.off()
## correlation heatmap
corr <- cor(immune_primary[samp_names, cells, drop = FALSE])
corr[!is.finite(corr)] <- 0
hc <- ComplexHeatmap::Heatmap(corr, name = "rho", show_row_names = FALSE, show_column_names = FALSE,
  col = colorRampPalette(c("blue","white","red"))(50))
grDevices::png(file.path(paths$res_immune, "immune_correlation_heatmap.png"), width = 800, height = 800, res = 120)
ComplexHeatmap::draw(hc)
grDevices::dev.off()
grDevices::tiff(file.path(paths$res_immune, "immune_correlation_heatmap.tiff"), width = 6, height = 6,
                units = "in", res = 600, compression = "lzw")
ComplexHeatmap::draw(hc)
grDevices::dev.off()

cat("Immune infiltration complete. Significant (adj.P<0.05) cells:", sum(diff_immune$adj_p < 0.05), "\n")
quit(save = "no", status = 0, runLast = FALSE)
