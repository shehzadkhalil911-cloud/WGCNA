# 07_differential_expression.R
source("config/analysis_config.R")
source("R/load_functions.R")
ensure_paths()
suppressMessages({ library(ggplot2); library(ggpubr); library(EnhancedVolcano); library(pheatmap); library(limma) })

disc <- readRDS(file.path(paths$discovery, "discovery_cohort_raw.rds"))
combat <- readRDS(file.path(paths$discovery, "discovery_combat.rds"))
meta <- disc$metadata
group <- meta$class

## discovery DE on ComBat matrix
deg <- limma_deg(combat, group, adj_p = ADJ_P_THRESHOLD, lfc = LOG2FC_THRESHOLD)
deg$symbol <- rownames(deg)
sig <- deg[deg$sig, ]
up <- deg[deg$sig & deg$log2FC > 0, ]
down <- deg[deg$sig & deg$log2FC < 0, ]

saveRDS(deg, file.path(paths$res_de, "discovery_DEG_full.rds"))
write.csv(deg, file.path(paths$res_de, "discovery_DEG_full.csv"), row.names = FALSE)
write.csv(sig, file.path(paths$res_de, "discovery_DEG_significant.csv"), row.names = FALSE)
write.csv(up, file.path(paths$res_de, "discovery_DEG_upregulated.csv"), row.names = FALSE)
write.csv(down, file.path(paths$res_de, "discovery_DEG_downregulated.csv"), row.names = FALSE)

cat("Discovery DEGs (adj.P<", ADJ_P_THRESHOLD, ", |log2FC|>=", LOG2FC_THRESHOLD, "):", nrow(sig),
    " (up:", nrow(up), ", down:", nrow(down), ")\n")

## volcano
p_volcano <- EnhancedVolcano::EnhancedVolcano(deg, lab = deg$symbol,
  x = "log2FC", y = "adj.P.Val", pCutoff = ADJ_P_THRESHOLD, FCcutoff = LOG2FC_THRESHOLD,
  title = "Discovery differential expression", subtitle = "ComBat-corrected")
save_publication_plot(p_volcano, file.path(paths$res_de, "volcano"), width = 8, height = 7)

## MA plot
deg$deg_status <- ifelse(deg$sig, ifelse(deg$log2FC>0,"up","down"),"ns")
p_ma <- ggplot(deg, aes(x = AveExpr, y = log2FC, colour = deg_status)) + geom_point(size = 0.5) +
  scale_colour_manual(values = c(ns = "grey", up = group_cols["Tumor"], down = group_cols["Normal"])) +
  labs(title = "MA plot") + pub_theme()
save_publication_plot(p_ma, file.path(paths$res_de, "MAplot"))

## top DEG heatmaps
for (n in c(25, 50, 100)) {
  topg <- head(deg$symbol[order(deg$adj.P.Val)], n)
  mat <- combat[topg, ]
  mat <- scale(mat)
  pheatmap::pheatmap(mat, annotation_col = data.frame(group = group, row.names = colnames(combat)),
    show_rownames = FALSE, main = paste0("Top ", n, " DEGs"),
    filename = file.path(paths$res_de, paste0("heatmap_top", n, ".png")))
}

## per-dataset DE for consistency
LOCAL_PROC <- file.path(LOCAL_DATA_DIR, "processed")
de_list <- list()
for (g in DISCOVERY_GSE) {
  obj <- readRDS(file.path(LOCAL_PROC, paste0(g, "_processed.rds")))
  gm <- obj$expression
  sc <- read.csv(file.path(LOCAL_DATA_DIR, "metadata", paste0(g, "_sample_classes.csv")), stringsAsFactors = FALSE)
  cls <- sc$class[match(colnames(gm), sc$sample_id)]
  tt <- limma_deg(gm, cls, adj_p = ADJ_P_THRESHOLD, lfc = LOG2FC_THRESHOLD)
  de_list[[g]] <- tt
  write.csv(tt, file.path(paths$res_de, paste0("DEG_", g, ".csv")), row.names = FALSE)
}
union_genes <- Reduce(union, lapply(de_list, function(x) x$gene))
cons <- dataset_consistency(de_list, union_genes)
write.csv(as.data.frame(cons$fc_matrix), file.path(paths$res_de, "per_dataset_foldchange_matrix.csv"))
saveRDS(cons, file.path(paths$res_de, "per_dataset_consistency.rds"))

## DEG count barplot
dc <- data.frame(dataset = DISCOVERY_GSE,
                 nsig = sapply(de_list, function(x) sum(x$sig)))
p_bar <- ggplot(dc, aes(x = dataset, y = nsig)) + geom_bar(stat = "identity", fill = "#59A14F") +
  labs(title = "Significant DEGs per dataset", y = "n DEGs") + pub_theme()
save_publication_plot(p_bar, file.path(paths$res_de, "DEG_count_barplot"))

cat("Differential expression complete.\n")
quit(save = "no", status = 0, runLast = FALSE)
