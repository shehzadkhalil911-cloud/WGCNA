# 19_gene_immune_crosstalk.R
source("config/analysis_config.R")
source("R/load_functions.R")
suppressMessages({ library(ggplot2); library(ggpubr); library(pheatmap); library(reshape2) })

combat <- readRDS(file.path(paths$discovery, "discovery_combat.rds"))
meta <- readRDS(file.path(paths$discovery, "discovery_cohort_raw.rds"))$metadata
group <- meta$class

xcell <- readRDS(file.path(paths$res_immune, "xcell_scores.rds"))
est <- readRDS(file.path(paths$res_immune, "estimate_scores.rds"))
consensus <- readRDS(file.path(paths$res_ml, "consensus", "consensus_genes.rds"))
locked <- consensus$final_genes
deg_crg <- readRDS(file.path(paths$res_cross, "DE_CRG_summary.rds"))

genes_of_interest <- unique(c(locked, deg_crg$de_crgs[deg_crg$de_crgs %in% rownames(combat)]))
genes_of_interest <- intersect(genes_of_interest, rownames(combat))
immune_mat <- as.data.frame(t(xcell))   # samples x cells

## correlation of each gene with each immune feature (spearman, by group)
samp <- colnames(combat)
immune_mat <- immune_mat[match(samp, rownames(immune_mat)), , drop = FALSE]
immune_mat$ESTIMATE_immune_score <- est$ESTIMATE_immune_score[samp]
immune_mat$ESTIMATE_stromal_score <- est$ESTIMATE_stromal_score[samp]
imm_cols <- colnames(immune_mat)
rho_mat <- matrix(NA, nrow = length(genes_of_interest), ncol = length(imm_cols),
  dimnames = list(genes_of_interest, imm_cols))
cor_rows <- list()
for (g in genes_of_interest) {
  gv <- as.numeric(combat[g, samp])
  for (ic in imm_cols) {
    iv <- as.numeric(immune_mat[[ic]])
    if (any(!is.finite(gv)) || any(!is.finite(iv))) next
    rs <- cor.test(gv, iv, method = "spearman")
    rho_mat[g, ic] <- rs$estimate
    cor_rows[[paste(g, ic, sep = "|")]] <- data.frame(gene = g, immune_feature = ic,
      rho = rs$estimate, p_value = rs$p.value)
  }
}
crosstalk <- do.call(rbind, cor_rows)
crosstalk$adj_p <- p.adjust(crosstalk$p_value, method = "BH")
crosstalk$group_effect <- ifelse(crosstalk$rho > 0, "positive", "negative")
write.csv(crosstalk, file.path(paths$res_immune, "gene_immune_crosstalk.csv"), row.names = FALSE)
suppressMessages(library(ComplexHeatmap)); suppressMessages(library(circlize))
grDevices::png(file.path(paths$res_immune, "gene_immune_crosstalk_heatmap.png"), width = 1200, height = 1000, res = 130)
col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
hc <- ComplexHeatmap::Heatmap(rho_mat, name = "rho", col = col_fun, na_col = "grey",
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 10), column_names_gp = grid::gpar(fontsize = 9))
ComplexHeatmap::draw(hc)
grDevices::dev.off()
grDevices::tiff(file.path(paths$res_immune, "gene_immune_crosstalk_heatmap.tiff"), width = 10, height = 8,
                units = "in", res = 600, compression = "lzw")
ComplexHeatmap::draw(hc)
grDevices::dev.off()

## scatter for top significant associations
sig <- crosstalk[crosstalk$adj_p < 0.05, ]
if (nrow(sig) > 0) {
  top <- head(sig[order(abs(sig$rho), decreasing = TRUE), ], 6)
  pl <- list()
  for (i in seq_len(nrow(top))) {
    g <- top$gene[i]; ic <- top$immune_feature[i]
    df <- data.frame(gene = gv, immune = immune_mat[[ic]], group = group)
    pl[[i]] <- ggplot(df, aes(x = gene, y = immune, color = group)) + geom_point(alpha = 0.5) +
      geom_smooth(method = "lm", se = FALSE) + labs(title = paste(g, "~", ic),
      x = g, y = ic) + pub_theme() + theme(legend.position = "none")
  }
  pcomb <- ggarrange(plotlist = pl, ncol = 3, nrow = 2)
  ggsave(file.path(paths$res_immune, "top_gene_immune_scatter.pdf"), pcomb, width = 12, height = 8)
}

cat("Gene-immune crosstalk complete. Genes:", length(genes_of_interest),
    " Significant (adj.P<0.05) associations:", nrow(sig), "\n")
quit(save = "no", status = 0, runLast = FALSE)
