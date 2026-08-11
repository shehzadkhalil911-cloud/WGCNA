# 10_functional_enrichment.R
source("config/analysis_config.R")
source("R/load_functions.R")
ensure_paths()
suppressMessages({ library(ggplot2); library(enrichplot); library(clusterProfiler); library(DOSE) })

deg <- readRDS(file.path(paths$res_de, "discovery_DEG_full.rds"))
crg_sum <- readRDS(file.path(paths$res_cross, "DE_CRG_summary.rds"))
universe <- rownames(readRDS(file.path(paths$discovery, "discovery_combat.rds")))

gene_sets <- list(
  all_DEG = deg$symbol[deg$sig],
  up_DEG = deg$symbol[deg$sig & deg$log2FC > 0],
  down_DEG = deg$symbol[deg$sig & deg$log2FC < 0],
  DE_CRG = crg_sum$de_crgs
)

## add consensus if available
cons_file <- file.path(paths$res_ml, "consensus", "consensus_genes.rds")
if (file.exists(cons_file)) {
  cg <- readRDS(cons_file)$final_genes
  if (length(cg) >= 5) gene_sets$consensus <- cg
}

enr_all <- list()
for (nm in names(gene_sets)) {
  cat("Enriching:", nm, " (", length(gene_sets[[nm]]), " genes)\n")
  res <- run_enrichment(gene_sets[[nm]], universe = universe)
  enr_all[[nm]] <- res
  for (ont in c("GO_BP","GO_CC","GO_MF","KEGG","Reactome","Hallmark")) {
    if (!is.null(res[[ont]]) && nrow(as_df(res[[ont]])) > 0) {
      write.csv(as_df(res[[ont]]),
                file.path(paths$res_enrich, paste0(nm, "_", ont, ".csv")), row.names = FALSE)
    }
  }
}
saveRDS(enr_all, file.path(paths$res_enrich, "enrichment_all.rds"))

suppressMessages(library(patchwork))
## key plots for all_DEG GO_BP and KEGG
plot_enr <- function(res, ont, fname, show = 10) {
  d <- as_df(res[[ont]])
  if (nrow(d) == 0) return(invisible(NULL))
  p <- enrichplot::dotplot(res[[ont]], showCategory = show, title = paste(ont, "enrichment"))
  save_publication_plot(p, file.path(paths$res_enrich, fname))
  invisible(NULL)
}
if (!is.null(enr_all$all_DEG$GO_BP)) plot_enr(enr_all$all_DEG, "GO_BP", "allDEG_GO_BP_dotplot")
if (!is.null(enr_all$all_DEG$KEGG)) plot_enr(enr_all$all_DEG, "KEGG", "allDEG_KEGG_barplot")
if (!is.null(enr_all$DE_CRG$GO_BP)) plot_enr(enr_all$DE_CRG, "GO_BP", "DECRG_GO_BP_dotplot")
if (!is.null(enr_all$DE_CRG$KEGG)) plot_enr(enr_all$DE_CRG, "KEGG", "DECRG_KEGG_dotplot")

## 2x2 combined enrichment panel: rows = DEGs / DE-CRGs, cols = BP / CC / MF / KEGG
panel_dot <- function(res, ont, title) {
  d <- as_df(res[[ont]])
  if (is.null(d) || nrow(d) == 0) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No significant terms",
                               size = 4) + labs(title = title) + theme_void())
  }
  enrichplot::dotplot(res[[ont]], showCategory = 8, title = title) + pub_theme(base_size = 9)
}
row_sets <- list(DEGs = "all_DEG", `DE-CRGs` = "DE_CRG")
col_onts <- list(BP = "GO_BP", CC = "GO_CC", MF = "GO_MF", KEGG = "KEGG")
pl_list <- list()
for (rn in names(row_sets)) {
  for (cn in names(col_onts)) {
    key <- paste0(rn, "__", cn)
    pl_list[[key]] <- panel_dot(enr_all[[row_sets[[rn]]]][[col_onts[[cn]]],
                                paste0(rn, " - ", cn))
  }
}
p_2x2 <- (pl_list[["DEGs__BP"]] + pl_list[["DEGs__CC"]] + pl_list[["DEGs__MF"]] + pl_list[["DEGs__KEGG"]]) /
         (pl_list[["DE-CRGs__BP"]] + pl_list[["DE-CRGs__CC"]] + pl_list[["DE-CRGs__MF"]] + pl_list[["DE-CRGs__KEGG"]])
save_publication_plot(p_2x2, file.path(paths$res_enrich, "enrichment_2x2_DEG_DECRG"), width = 16, height = 9)
cat("2x2 enrichment panel (DEGs vs DE-CRGs x BP/CC/MF/KEGG) written.\n")

## combined bar plot: number of significant terms per category for all DEGs (5 categories)
cat_categories <- c(GO_BP = "GO_BP", GO_CC = "GO_CC", GO_MF = "GO_MF", KEGG = "KEGG", Reactome = "Reactome")
bar_counts <- sapply(cat_categories, function(ont) {
  d <- as_df(enr_all$all_DEG[[ont]])
  if (is.null(d) || nrow(d) == 0) 0L else nrow(d)
})
bar_df <- data.frame(category = names(cat_categories), n_terms = as.integer(bar_counts),
                     stringsAsFactors = FALSE)
write.csv(bar_df, file.path(paths$res_enrich, "allDEG_category_counts.csv"), row.names = FALSE)
p_bar <- ggplot(bar_df, aes(x = reorder(category, n_terms), y = n_terms, fill = category)) +
  geom_bar(stat = "identity") + coord_flip() +
  labs(title = "Enriched terms per category (all DEGs)", x = "Category", y = "Number of significant terms") +
  scale_fill_brewer(palette = "Set2") + pub_theme() + theme(legend.position = "none")
save_publication_plot(p_bar, file.path(paths$res_enrich, "allDEG_combined_barplot"), width = 7, height = 5)

cat("Functional enrichment complete for", length(gene_sets), "gene sets.\n")
quit(save = "no", status = 0, runLast = FALSE)
