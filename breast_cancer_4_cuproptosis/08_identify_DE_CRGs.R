# 08_identify_DE_CRGs.R
source("config/analysis_config.R")
source("R/load_functions.R")
ensure_paths()
suppressMessages({ library(ggplot2); library(ggpubr); library(pheatmap); library(ggvenn); library(UpSetR) })

deg <- readRDS(file.path(paths$res_de, "discovery_DEG_full.rds"))
crg <- readRDS(file.path(paths$crg, "cuproptosis_standardized.rds"))
crg_sym <- crg$standardized$official_symbol
disc <- readRDS(file.path(paths$discovery, "discovery_cohort_raw.rds"))
combat <- readRDS(file.path(paths$discovery, "discovery_combat.rds"))
meta <- disc$metadata

## CRGs detected on platform
crg_detected <- crg_sym[crg_sym %in% rownames(combat)]
crg_not_detected <- setdiff(crg_sym, rownames(combat))
cat("CRGs total:", length(crg_sym), " detected:", length(crg_detected), " not detected:", length(crg_not_detected), "\n")

## DE-CRGs = significant DEG ∩ CRG
deg_sig_symbols <- deg$symbol[deg$sig]
de_crgs <- intersect(crg_detected, deg_sig_symbols)
de_crgs_up <- de_crgs[de_crgs %in% deg$symbol[deg$sig & deg$log2FC > 0]]
de_crgs_down <- setdiff(de_crgs, de_crgs_up)
non_sig_crg <- setdiff(crg_detected, deg_sig_symbols)

cat("DE-CRGs:", length(de_crgs), " (up:", length(de_crgs_up), ", down:", length(de_crgs_down), ")\n")

## sensitivity expansion if too few
if (length(de_crgs) < MIN_CRG_OVERLAP && RELAX_FC_FOR_CRG_POOL) {
  de_crgs_relax <- intersect(crg_detected, deg$symbol[deg$adj.P.Val < ADJ_P_THRESHOLD])
  cat("SENSITIVITY: relaxing FC threshold ->", length(de_crgs_relax), "DE-CRGs (adj.P<", ADJ_P_THRESHOLD, "only)\n")
  candidate_crg_pool <- de_crgs_relax
  attr(candidate_crg_pool, "relaxed") <- TRUE
} else {
  candidate_crg_pool <- de_crgs
  attr(candidate_crg_pool, "relaxed") <- FALSE
}

## save tables
crg_expr <- combat[crg_detected, ]
saveRDS(crg_expr, file.path(paths$res_cross, "CRG_expression.rds"))
write.csv(as.data.frame(crg_expr), file.path(paths$res_cross, "CRG_expression.csv"))
de_crg_tab <- deg[deg$symbol %in% de_crgs, c("symbol","log2FC","AveExpr","t","P.Value","adj.P.Val","B")]
write.csv(de_crg_tab, file.path(paths$res_cross, "DE_CRGs.csv"), row.names = FALSE)
write.csv(data.frame(gene = de_crgs_up), file.path(paths$res_cross, "DE_CRGs_upregulated.csv"), row.names = FALSE)
write.csv(data.frame(gene = de_crgs_down), file.path(paths$res_cross, "DE_CRGs_downregulated.csv"), row.names = FALSE)
write.csv(data.frame(gene = crg_not_detected), file.path(paths$res_cross, "CRGs_not_detected.csv"), row.names = FALSE)
write.csv(data.frame(gene = non_sig_crg), file.path(paths$res_cross, "CRGs_nonsignificant.csv"), row.names = FALSE)
saveRDS(list(de_crgs = de_crgs, de_crgs_up = de_crgs_up, de_crgs_down = de_crgs_down,
             crg_detected = crg_detected, crg_not_detected = crg_not_detected,
             candidate_pool = candidate_crg_pool),
        file.path(paths$res_cross, "DE_CRG_summary.rds"))

## Venn diagram
png(file.path(paths$res_cross, "DEG_CRG_venn.png"), width = 6, height = 6, units = "in", res = FIG_DPI)
ggvenn::ggvenn(list(DEGs = deg_sig_symbols, CRGs = crg_detected))
dev.off()
tiff(file.path(paths$res_cross, "DEG_CRG_venn.tiff"), width = 6, height = 6, units = "in", res = 600, compression = "lzw")
ggvenn::ggvenn(list(DEGs = deg_sig_symbols, CRGs = crg_detected))
dev.off()
pdf(file.path(paths$res_cross, "DEG_CRG_venn.pdf"), width = 6, height = 6)
ggvenn::ggvenn(list(DEGs = deg_sig_symbols, CRGs = crg_detected))
dev.off()

## DE-CRG heatmap
mat <- scale(combat[de_crgs, ])
pheatmap::pheatmap(mat, annotation_col = data.frame(group = meta$class, row.names = colnames(combat)),
  main = "DE-CRG heatmap", filename = file.path(paths$res_cross, "DE_CRG_heatmap.png"), show_rownames = TRUE)

## DE-CRG box plots
boxdf <- data.frame(gene = rep(de_crgs, each = ncol(combat)),
                    expr = as.vector(t(combat[de_crgs, ])),
                    group = rep(meta$class, times = length(de_crgs)))
p_box <- ggplot(boxdf, aes(x = gene, y = expr, fill = group)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(title = "DE-CRG expression") + pub_theme()
save_publication_plot(p_box, file.path(paths$res_cross, "DE_CRG_boxplots"), width = 8, height = 6)

## DE-CRG fold-change bar plot
fcdf <- deg[deg$symbol %in% de_crgs, c("symbol","log2FC")]
p_fc <- ggplot(fcdf, aes(x = reorder(symbol, log2FC), y = log2FC, fill = log2FC > 0)) +
  geom_bar(stat = "identity") + coord_flip() + labs(title = "DE-CRG log2FC") + pub_theme()
save_publication_plot(p_fc, file.path(paths$res_cross, "DE_CRG_foldchange"), width = 6, height = 6)

cat("DE-CRG identification complete. Candidate pool size:", length(candidate_crg_pool),
    if (attr(candidate_crg_pool,"relaxed")) " (RELAXED FC)" else "", "\n")
quit(save = "no", status = 0, runLast = FALSE)
