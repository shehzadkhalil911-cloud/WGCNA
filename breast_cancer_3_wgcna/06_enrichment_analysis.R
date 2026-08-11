###############################################################################
# SCRIPT 6: Functional Enrichment Analysis
# Only dotplot + circular cnetplot, TIFF 600 dpi
###############################################################################
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(msigdbr)
library(DOSE)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

options(stringsAsFactors = FALSE)

OUT_DIR <- "F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)"
FIG_DIR <- file.path(OUT_DIR, "figures/enrichment")
TAB_DIR <- file.path(OUT_DIR, "tables")
RDA_DIR <- file.path(OUT_DIR, "RData")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

load(file.path(RDA_DIR, "05_hub_genes.RData"))

run_enrichment <- function(gene_list, name, output_prefix) {
  if (length(gene_list) < 5) {
    cat(sprintf("Skipping %s: too few genes (%d)\n", name, length(gene_list)))
    return(NULL)
  }

  gene_df <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db, drop = TRUE)
  if (is.null(gene_df) || nrow(gene_df) < 3) {
    cat(sprintf("Skipping %s: too few mapped genes\n", name))
    return(NULL)
  }

  entrez_ids <- gene_df$ENTREZID
  results <- list()

  results$go_bp <- enrichGO(entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP",
                             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  results$go_cc <- enrichGO(entrez_ids, OrgDb = org.Hs.eg.db, ont = "CC",
                             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  results$go_mf <- enrichGO(entrez_ids, OrgDb = org.Hs.eg.db, ont = "MF",
                             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  results$kegg <- enrichKEGG(entrez_ids, organism = "hsa",
                              pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  results$reactome <- enrichPathway(entrez_ids, organism = "human",
                                     pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  msig_h <- msigdbr(species = "Homo sapiens", category = "H")
  results$hallmark <- enricher(entrez_ids, TERM2GENE = msig_h[, c("gs_name", "entrez_gene")],
                                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)

  # Save tables
  for (onto in names(results)) {
    if (!is.null(results[[onto]]) && nrow(results[[onto]]) > 0) {
      write.csv(as.data.frame(results[[onto]]),
                file.path(TAB_DIR, paste0(output_prefix, "_", onto, ".csv")), row.names = FALSE)
    }
  }

  # Dotplot + circular Cnetplot for all 6 ontologies
  all_ontos <- c("go_bp", "go_cc", "go_mf", "kegg", "reactome", "hallmark")
  for (onto in all_ontos) {
    if (!is.null(results[[onto]]) && nrow(results[[onto]]) > 0) {
      # Dotplot
      p_dot <- dotplot(results[[onto]], showCategory = 15, title = paste(name, "-", onto))
      ggsave(file.path(FIG_DIR, paste0(output_prefix, "_", onto, "_dotplot.tiff")),
             p_dot, width = 10, height = 8, dpi = 300, compression = "lzw", limitsize = FALSE)

      # Cnetplot
      tryCatch({
        p_cnet <- cnetplot(results[[onto]], showCategory = 3)
        tiff(file.path(FIG_DIR, paste0(output_prefix, "_", onto, "_cnet.tiff")),
             width = 8, height = 6, units = "in", res = 300, compression = "lzw")
        print(p_cnet)
        dev.off()
      }, error = function(e) {
        try(dev.off(), silent = TRUE)
        cat(sprintf("  cnetplot failed for %s: %s\n", onto, e$message))
      })
    } else {
      cat(sprintf("  Skipping %s: no significant terms\n", onto))
    }
  }

  results
}

cat("Running enrichment on final hub genes...\n")
enrich_hubs <- run_enrichment(final_hubs, "Final Hub Genes", "hub_genes")

cat("\nEnrichment complete. TIFF plots saved at 300 dpi.\n")
save(enrich_hubs, file = file.path(RDA_DIR, "06_enrichment.RData"))
cat("Script 6 complete.\n")
