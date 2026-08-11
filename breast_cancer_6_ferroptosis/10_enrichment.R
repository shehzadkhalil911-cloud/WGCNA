source("R/utils.R")
load_libraries()
set.seed(123)

BASE_DIR <- "F:/independent projects/cancer/breast cancer"
OUT_DIR <- file.path(BASE_DIR, "Breast cancer 6 Ferroptosis")

log_msg("=== Step 10: Functional Enrichment Analysis ===")

deg_table <- readRDS(file.path(OUT_DIR, "results/differential_expression/combined_DEG_results.rds"))
de_frg <- readRDS(file.path(OUT_DIR, "results/ferroptosis_genes/DE_FRG_results.rds"))
module_info <- read.csv(file.path(OUT_DIR, "results/WGCNA/module_assignments.csv"), row.names = 1)
top_hubs <- readRDS(file.path(OUT_DIR, "results/WGCNA/top_hub_genes.rds"))

library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(org.Hs.eg.db)

# Get Entrez IDs for genes
get_entrez <- function(genes) {
  suppressMessages(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db))
}

run_enrichment <- function(gene_set, name, filename_prefix) {
  if (length(gene_set) < 5) {
    log_msg(sprintf("  %s: too few genes (%d), skipping", name, length(gene_set)))
    return(NULL)
  }

  entrez <- get_entrez(gene_set)
  if (is.null(entrez) || nrow(entrez) < 3) {
    log_msg(sprintf("  %s: too few mapped genes (%d), skipping", name, nrow(entrez)))
    return(NULL)
  }

  results <- list()

  # GO BP
  go_bp <- tryCatch(enrichGO(entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                              ont = "BP", pAdjustMethod = "BH",
                              pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                              readable = TRUE),
                    error = function(e) NULL)
  # GO CC
  go_cc <- tryCatch(enrichGO(entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                              ont = "CC", pAdjustMethod = "BH",
                              pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                              readable = TRUE),
                    error = function(e) NULL)
  # GO MF
  go_mf <- tryCatch(enrichGO(entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                              ont = "MF", pAdjustMethod = "BH",
                              pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                              readable = TRUE),
                    error = function(e) NULL)
  # KEGG
  kegg <- tryCatch(enrichKEGG(entrez$ENTREZID, organism = "hsa",
                               pAdjustMethod = "BH",
                               pvalueCutoff = 0.05, qvalueCutoff = 0.2),
                   error = function(e) NULL)
  # Reactome
  reactome <- tryCatch(enrichPathway(entrez$ENTREZID, organism = "human",
                                      pAdjustMethod = "BH",
                                      pvalueCutoff = 0.05, qvalueCutoff = 0.2),
                       error = function(e) NULL)

  # Hallmark (via msigdbr - if available)
  hallmark <- NULL
  if (requireNamespace("msigdbr", quietly = TRUE)) {
    h_all <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    h_gene_sets <- split(h_all$gene_symbol, h_all$gs_name)
    hallmark <- tryCatch(enricher(entrez$SYMBOL, TERM2GENE = stack(h_gene_sets)[, c(2, 1)],
                                   pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
                          error = function(e) NULL)
  }

  # Save results
  enrich_dir <- file.path(OUT_DIR, "results/enrichment")
  if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
    write.csv(go_bp@result, file.path(enrich_dir, paste0(filename_prefix, "_GO_BP.csv")), row.names = FALSE)
    results$go_bp <- go_bp
  }
  if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
    write.csv(go_cc@result, file.path(enrich_dir, paste0(filename_prefix, "_GO_CC.csv")), row.names = FALSE)
    results$go_cc <- go_cc
  }
  if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
    write.csv(go_mf@result, file.path(enrich_dir, paste0(filename_prefix, "_GO_MF.csv")), row.names = FALSE)
    results$go_mf <- go_mf
  }
  if (!is.null(kegg) && nrow(kegg@result) > 0) {
    write.csv(kegg@result, file.path(enrich_dir, paste0(filename_prefix, "_KEGG.csv")), row.names = FALSE)
    results$kegg <- kegg
  }
  if (!is.null(reactome) && nrow(reactome@result) > 0) {
    write.csv(reactome@result, file.path(enrich_dir, paste0(filename_prefix, "_REACTOME.csv")), row.names = FALSE)
    results$reactome <- reactome
  }
  if (!is.null(hallmark) && nrow(hallmark@result) > 0) {
    write.csv(hallmark@result, file.path(enrich_dir, paste0(filename_prefix, "_HALLMARK.csv")), row.names = FALSE)
    results$hallmark <- hallmark
  }

  log_msg(sprintf("  %s: GO_BP=%d, GO_CC=%d, GO_MF=%d, KEGG=%d, Reactome=%d, Hallmark=%d",
                  name,
                  ifelse(is.null(go_bp), 0, sum(go_bp@result$p.adjust < 0.05)),
                  ifelse(is.null(go_cc), 0, sum(go_cc@result$p.adjust < 0.05)),
                  ifelse(is.null(go_mf), 0, sum(go_mf@result$p.adjust < 0.05)),
                  ifelse(is.null(kegg), 0, sum(kegg@result$p.adjust < 0.05)),
                  ifelse(is.null(reactome), 0, sum(reactome@result$p.adjust < 0.05)),
                  ifelse(is.null(hallmark), 0, sum(hallmark@result$p.adjust < 0.05))))

  return(results)
}

# Run enrichment for various gene sets
sig_deg_genes <- rownames(deg_table)[deg_table$Significant]
up_deg_genes <- rownames(deg_table)[deg_table$logFC > 0 & deg_table$adj.P.Val < 0.05]
down_deg_genes <- rownames(deg_table)[deg_table$logFC < 0 & deg_table$adj.P.Val < 0.05]
de_frg_genes <- rownames(de_frg)
hub_genes <- top_hubs$Gene

log_msg("Running enrichment for all DEGs...")
e_all <- run_enrichment(sig_deg_genes, "All DEGs", "all_DEGs")

log_msg("Running enrichment for upregulated DEGs...")
e_up <- run_enrichment(up_deg_genes, "Upregulated DEGs", "up_DEGs")

log_msg("Running enrichment for downregulated DEGs...")
e_down <- run_enrichment(down_deg_genes, "Downregulated DEGs", "down_DEGs")

log_msg("Running enrichment for DE-FRGs...")
e_defrg <- run_enrichment(de_frg_genes, "DE-FRGs", "DE_FRGs")

log_msg("Running enrichment for hub genes...")
e_hubs <- run_enrichment(hub_genes, "Hub genes", "hub_genes")

# GSEA
ranked_genes <- deg_table$t
names(ranked_genes) <- rownames(deg_table)
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

gsea_res <- tryCatch(gseGO(ranked_genes, OrgDb = org.Hs.eg.db,
                            ont = "BP", keyType = "SYMBOL",
                            pAdjustMethod = "BH", pvalueCutoff = 0.05),
                     error = function(e) NULL)
if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
  write.csv(gsea_res@result, file.path(OUT_DIR, "results/enrichment/GSEA_GO_BP.csv"), row.names = FALSE)
  saveRDS(gsea_res, file.path(OUT_DIR, "results/enrichment/GSEA_results.rds"))

  # GSEA ridge plot
  if (nrow(gsea_res@result) > 1) {
    p_ridge <- ridgeplot(gsea_res, showCategory = 15) +
      labs(title = "GSEA Enrichment Results")
    save_tiff_600dpi(p_ridge, file.path(OUT_DIR, "results/figures/Figure8E_GSEA_ridge.tiff"),
                     width = 9, height = 7)

    # Top GSEA enrichment curves
    top_gsea <- head(gsea_res@result$ID, 6)
    if (length(top_gsea) > 0) {
      p_gsea_curves <- gseaplot2(gsea_res, geneSetID = top_gsea,
                                  title = "Top GSEA Enrichment Curves")
      save_tiff_600dpi(p_gsea_curves, file.path(OUT_DIR, "results/figures/Figure8F_GSEA_curves.tiff"),
                       width = 9, height = 7)
    }
  }
}

# GO dot plot (using first available GO result)
for (res_name in c("go_bp")) {
  res <- NULL
  if (!is.null(e_all[[res_name]])) res <- e_all[[res_name]]
  else if (!is.null(e_up[[res_name]])) res <- e_up[[res_name]]

  if (!is.null(res) && nrow(res@result) > 0) {
    p_dot <- dotplot(res, showCategory = 15) + labs(title = "GO Biological Process Enrichment")
    save_tiff_600dpi(p_dot, file.path(OUT_DIR, "results/figures/Figure8A_GO_dot.tiff"))

    p_bar <- barplot(res, showCategory = 15) + labs(title = "GO BP Bar Plot")
    save_tiff_600dpi(p_bar, file.path(OUT_DIR, "results/figures/Figure8B_GO_bar.tiff"))
  }
}

# KEGG dot plot
if (!is.null(e_all$kegg) && nrow(e_all$kegg@result) > 0) {
  p_kegg_dot <- dotplot(e_all$kegg, showCategory = 15) + labs(title = "KEGG Pathway Enrichment")
  save_tiff_600dpi(p_kegg_dot, file.path(OUT_DIR, "results/figures/Figure8C_KEGG_dot.tiff"))

  p_kegg_bar <- barplot(e_all$kegg, showCategory = 15) + labs(title = "KEGG Bar Plot")
  save_tiff_600dpi(p_kegg_bar, file.path(OUT_DIR, "results/figures/Figure8D_KEGG_bar.tiff"))
}

# Reactome dot plot
if (!is.null(e_all$reactome) && nrow(e_all$reactome@result) > 0) {
  p_reactome <- dotplot(e_all$reactome, showCategory = 15) + labs(title = "Reactome Enrichment")
  save_tiff_600dpi(p_reactome, file.path(OUT_DIR, "results/figures/Figure8E_Reactome_dot.tiff"))
}

# Module-wise enrichment heatmap
modules_of_interest <- unique(module_info$Module)
module_enrichment_list <- list()
for (mod in modules_of_interest) {
  mod_genes <- module_info$Gene[module_info$Module == mod]
  if (length(mod_genes) >= 5) {
    entrez <- get_entrez(mod_genes)
    if (!is.null(entrez) && nrow(entrez) >= 3) {
      kegg_mod <- tryCatch(enrichKEGG(entrez$ENTREZID, organism = "hsa",
                                       pAdjustMethod = "BH", pvalueCutoff = 0.1),
                            error = function(e) NULL)
      if (!is.null(kegg_mod) && nrow(kegg_mod@result) > 0) {
        module_enrichment_list[[mod]] <- kegg_mod@result$Description[1:min(5, nrow(kegg_mod@result))]
      }
    }
  }
}

if (length(module_enrichment_list) > 0) {
  saveRDS(module_enrichment_list, file.path(OUT_DIR, "results/enrichment/module_enrichment_list.rds"))
}

log_msg("Step 10 complete: Functional enrichment analysis finished")
