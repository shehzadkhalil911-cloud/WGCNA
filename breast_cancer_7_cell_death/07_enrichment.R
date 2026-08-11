# =============================================================================
# Script 07: Functional Enrichment Analysis
# Breast Cancer 7 Cell Death Crosstalk
# =============================================================================

source(file.path("..", "config", "analysis_config.R"))

library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(msigdbr)
library(org.Hs.eg.db)
library(tidyverse)

log_message("Starting enrichment analysis.")

# --- 1. LOAD DATA ---
de_results <- readRDS(file.path(OUTPUT$de, "discovery_DEG_full.rds"))
deg_sig <- read.csv(file.path(OUTPUT$de, "discovery_DEG_significant.csv"), row.names = 1)
shared_info <- readRDS(file.path(OUTPUT$shared, "shared_gene_classification.rds"))

# --- 2. GENE LISTS ---
gene_sets <- list(
  All_DEGs = rownames(deg_sig),
  Up_DEGs = rownames(deg_sig[deg_sig$logFC > 0, ]),
  Down_DEGs = rownames(deg_sig[deg_sig$logFC < 0, ]),
  DE_FRGs = shared_info$Gene[shared_info$DE_FRG],
  DE_NRGs = shared_info$Gene[shared_info$DE_NRG],
  DE_CRGs = shared_info$Gene[shared_info$DE_CRG]
)
gene_sets <- gene_sets[sapply(gene_sets, length) >= 5]

# --- 3. ENTREZ MAPPING ---
map_to_entrez <- function(genes) {
  entrez <- mapIds(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
  as.character(entrez[!is.na(entrez)])
}

# --- 4. RUN ENRICHMENT ---
all_enrichment <- list()

msig_h <- msigdbr(species = "Homo sapiens", category = "H")

for (gs_name in names(gene_sets)) {
  genes <- gene_sets[[gs_name]]
  entrez <- map_to_entrez(genes)
  log_message(sprintf("%s: %d symbols -> %d Entrez IDs", gs_name, length(genes), length(entrez)))
  if (length(entrez) < 5) next

  res <- list()
  types <- if (gs_name == "All_DEGs") c("GO_BP", "Reactome", "Hallmark") else "GO_BP"

  for (t in types) {
    r <- switch(t,
      GO_BP = tryCatch(enrichGO(entrez, org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                                 pvalueCutoff = 0.05, readable = TRUE), error = function(e) NULL),
      Reactome = tryCatch(enrichPathway(entrez, organism = "human", pAdjustMethod = "BH",
                                        pvalueCutoff = 0.05, readable = TRUE), error = function(e) NULL),
      Hallmark = tryCatch(enricher(genes,
                                    TERM2GENE = data.frame(gs_name = msig_h$gs_name, gene_symbol = msig_h$gene_symbol),
                                    pAdjustMethod = "BH", pvalueCutoff = 0.05), error = function(e) NULL)
    )
    res[[t]] <- r
    if (!is.null(r) && nrow(r@result) > 0) {
      write.csv(r@result, file.path(OUTPUT$enrichment, sprintf("%s_%s.csv", gs_name, t)), row.names = FALSE)
    }
  }

  all_enrichment[[gs_name]] <- res
  log_message(sprintf("%s: GO_BP=%d, Reactome=%d, Hallmark=%d", gs_name,
                      ifelse(!is.null(res$GO_BP), nrow(res$GO_BP@result), 0),
                      ifelse(!is.null(res$Reactome), nrow(res$Reactome@result), 0),
                      ifelse(!is.null(res$Hallmark), nrow(res$Hallmark@result), 0)))
}
saveRDS(all_enrichment, file.path(OUTPUT$enrichment, "all_enrichment_results.rds"))

# --- 5. GSEA ---
ranked_list <- de_results$logFC
names(ranked_list) <- rownames(de_results)
ranked_list <- sort(ranked_list, decreasing = TRUE)
rank_entrez <- mapIds(org.Hs.eg.db, keys = names(ranked_list), keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
names(ranked_list) <- as.character(rank_entrez[names(ranked_list)])
ranked_list <- ranked_list[!is.na(names(ranked_list)) & !duplicated(names(ranked_list))]

if (length(ranked_list) < 10) {
  log_message("Too few genes for GSEA, skipping.", "WARNING")
  gsea_go <- NULL
} else {
  gsea_go <- tryCatch(gseGO(ranked_list, org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                           pvalueCutoff = 0.05, verbose = FALSE), error = function(e) NULL)
  if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
    write.csv(gsea_go@result, file.path(OUTPUT$enrichment, "GSEA_GO_BP.csv"), row.names = FALSE)
    saveRDS(gsea_go, file.path(OUTPUT$enrichment, "GSEA_GO_BP.rds"))
  }
  gsea_kegg <- tryCatch(gseKEGG(ranked_list, organism = "hsa", pAdjustMethod = "BH",
                                pvalueCutoff = 0.05, verbose = FALSE, use_internal_data = TRUE),
                        error = function(e) NULL)
  if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
    write.csv(gsea_kegg@result, file.path(OUTPUT$enrichment, "GSEA_KEGG.csv"), row.names = FALSE)
    saveRDS(gsea_kegg, file.path(OUTPUT$enrichment, "GSEA_KEGG.rds"))
  }
}

# --- 6. FIGURES ---
if (!is.null(all_enrichment$All_DEGs$GO_BP) && nrow(all_enrichment$All_DEGs$GO_BP@result) > 0) {
  p <- dotplot(all_enrichment$All_DEGs$GO_BP, showCategory = 20, font.size = 10) +
    ggtitle("GO Biological Process (All DEGs)")
  save_tiff_600dpi(p, file.path(OUTPUT$enrichment, "GO_BP_dotplot_all_DEGs.tiff"))
}
if (!is.null(all_enrichment$All_DEGs$Reactome) && nrow(all_enrichment$All_DEGs$Reactome@result) > 0) {
  p <- dotplot(all_enrichment$All_DEGs$Reactome, showCategory = 20, font.size = 10) +
    ggtitle("Reactome Pathways (All DEGs)")
  save_tiff_600dpi(p, file.path(OUTPUT$enrichment, "Reactome_dotplot_all_DEGs.tiff"))
}
if (!is.null(all_enrichment$All_DEGs$Hallmark) && nrow(all_enrichment$All_DEGs$Hallmark@result) > 0) {
  p <- dotplot(all_enrichment$All_DEGs$Hallmark, showCategory = 20, font.size = 10) +
    ggtitle("Hallmark Gene Sets (All DEGs)")
  save_tiff_600dpi(p, file.path(OUTPUT$enrichment, "Hallmark_dotplot_all_DEGs.tiff"))
}

# GSEA curves
if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
  top_gsea <- head(gsea_go@result$ID, 10)
  p <- gseaplot2(gsea_go, geneSetID = top_gsea, title = "Top 10 GSEA GO Terms", pvalue_table = FALSE)
  ggsave(file.path(OUTPUT$enrichment, "GSEA_GO_curves.tiff"),
         plot = p, width = 10, height = 8, dpi = 600, compression = "lzw")
}

log_message("Enrichment analysis complete.")
