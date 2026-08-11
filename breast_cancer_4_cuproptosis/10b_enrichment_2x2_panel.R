# 10b_enrichment_2x2_panel.R
# Build the 2x2 DEGs/DE-CRGs x BP/CC/MF/KEGG panel from saved enrichment,
# computing only the (previously-missing) KEGG terms for the two gene sets.
source("config/analysis_config.R")
source("R/load_functions.R")
suppressMessages({ library(ggplot2); library(enrichplot); library(patchwork); library(clusterProfiler); library(org.Hs.eg.db) })

base <- paths$res_enrich
e <- readRDS(file.path(base, "enrichment_all.rds"))

deg <- readRDS(file.path(paths$res_de, "discovery_DEG_full.rds"))
crg_sum <- readRDS(file.path(paths$res_cross, "DE_CRG_summary.rds"))
all_deg_genes <- unique(deg$symbol[deg$sig])
decrg_genes <- crg_sum$de_crgs

run_kegg <- function(genes) {
  eg <- clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  if (nrow(eg) == 0) return(NULL)
  clusterProfiler::enrichKEGG(gene = eg$ENTREZID, organism = "hsa", keyType = "ncbi-geneid",
                              pvalueCutoff = 0.05, pAdjustMethod = "BH")
}
cat("Computing KEGG for all_DEG...\n"); flush.console()
e$all_DEG$KEGG <- run_kegg(all_deg_genes)
cat("Computing KEGG for DE_CRG...\n"); flush.console()
e$DE_CRG$KEGG <- run_kegg(decrg_genes)
saveRDS(e, file.path(base, "enrichment_all.rds"))

panel_bar <- function(res, title) {
  d <- if (is.null(res)) data.frame() else as.data.frame(res)
  if (is.null(d) || nrow(d) == 0)
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No significant terms", size = 4) +
             labs(title = title) + theme_void())
  barplot(res, showCategory = 5, title = title) + pub_theme()
}
row_sets <- list(DEGs = "all_DEG", `DE-CRGs` = "DE_CRG")
col_onts <- list(BP = "GO_BP", CC = "GO_CC", MF = "GO_MF", KEGG = "KEGG")
pl <- list()
for (rn in names(row_sets)) for (cn in names(col_onts)) {
  pl[[paste0(rn, "__", cn)]] <- panel_bar(e[[row_sets[[rn]]]][[col_onts[[cn]]]], cn)
}
## one figure per gene set: 1x4 panel (BP | CC | MF | KEGG), barplots top 5 per pathway
p_deg <- pl[["DEGs__BP"]] + pl[["DEGs__CC"]] + pl[["DEGs__MF"]] + pl[["DEGs__KEGG"]] +
  plot_annotation(title = "Functional enrichment of differentially expressed genes (DEGs)")
p_decrg <- pl[["DE-CRGs__BP"]] + pl[["DE-CRGs__CC"]] + pl[["DE-CRGs__MF"]] + pl[["DE-CRGs__KEGG"]] +
  plot_annotation(title = "Functional enrichment of downregulated DE-CRGs")
save_publication_plot(p_deg, file.path(base, "enrichment_DEG_bar"), width = 16, height = 4.5)
save_publication_plot(p_decrg, file.path(base, "enrichment_DECRG_bar"), width = 16, height = 4.5)

## also rewrite the per-set KEGG csvs now that KEGG exists
for (nm in c("all_DEG","DE_CRG")) {
  kk <- e[[nm]]$KEGG
  d <- if (is.null(kk)) data.frame() else as.data.frame(kk)
  if (nrow(d) > 0) write.csv(d, file.path(base, paste0(nm, "_KEGG.csv")), row.names = FALSE)
}
cat("Separate enrichment panels written: DEGs (enrichment_DEG_1x4) and DE-CRGs (enrichment_DECRG_1x4).\n")
quit(save = "no", status = 0, runLast = FALSE)
