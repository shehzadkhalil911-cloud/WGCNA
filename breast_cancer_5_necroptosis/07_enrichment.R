# ============================================================
# Step 7: Functional Enrichment Analysis
# ============================================================
set.seed(123)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichR)
library(DOSE)

base_dir <- normalizePath(file.path(".", "Breast Cancer 5 Necroptosis"))
output_dir <- base_dir

sink(file.path(output_dir, "logs", "07_enrichment.log"), split=TRUE)
cat("FUNCTIONAL ENRICHMENT ANALYSIS\n")
cat("==============================\n\n")

save_tiff_600dpi <- function(plot, filename, width = 7, height = 6) {
  ggsave(filename=filename, plot=plot, width=width, height=height,
         units="in", dpi=600, compression="lzw", bg="white", device="tiff")
}

deg_sig <- read.csv(file.path(output_dir, "results", "tables", "DEG_significant.csv"), row.names=1)
de_nrgs <- readLines(file.path(output_dir, "results", "necroptosis_genes", "DE_NRGs.txt"))
deg_up <- read.csv(file.path(output_dir, "results", "tables", "DEG_upregulated.csv"), row.names=1)
deg_down <- read.csv(file.path(output_dir, "results", "tables", "DEG_downregulated.csv"), row.names=1)
ranked_list <- readRDS(file.path(output_dir, "results", "differential_expression", "ranked_gene_list.rds"))

# Helper: run enrichment for a gene list
run_enrichment <- function(gene_list, label, output_subdir) {
  cat("\n--- Enrichment for:", label, "---\n")
  cat("Genes:", length(gene_list), "\n")
  
  if(length(gene_list) < 10) {
    cat("  Too few genes for meaningful enrichment. Skipping.\n")
    return()
  }
  
  # Convert symbols to ENTREZ IDs
  gene_entrez <- bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
  entrez_ids <- unique(gene_entrez$ENTREZID)
  cat("  Mapped to ENTREZ:", length(entrez_ids), "\n")
  
  if(length(entrez_ids) < 5) {
    cat("  Too few ENTREZ IDs. Skipping clusterProfiler.\n")
    return()
  }
  
  # GO BP
  ego_bp <- tryCatch(enrichGO(gene=entrez_ids, OrgDb=org.Hs.eg.db,
                               ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05,
                               qvalueCutoff=0.2, readable=TRUE), error=function(e) NULL)
  if(!is.null(ego_bp) && nrow(ego_bp) > 0) {
    write.csv(as.data.frame(ego_bp), file.path(output_dir, "results", "enrichment", output_subdir,
              paste0(label, "_GO_BP.csv")), row.names=FALSE)
    p <- dotplot(ego_bp, showCategory=15) + ggtitle(paste("GO BP -", label))
    save_tiff_600dpi(p, file.path(output_dir, "results", "enrichment", output_subdir,
                    paste0(label, "_GO_BP_dotplot.tiff")), width=9, height=7)
    p2 <- barplot(ego_bp, showCategory=15) + ggtitle(paste("GO BP -", label))
    save_tiff_600dpi(p2, file.path(output_dir, "results", "enrichment", output_subdir,
                     paste0(label, "_GO_BP_bar.tiff")), width=9, height=7)
  }
  
  # GO CC
  ego_cc <- tryCatch(enrichGO(gene=entrez_ids, OrgDb=org.Hs.eg.db,
                               ont="CC", pAdjustMethod="BH", pvalueCutoff=0.05,
                               qvalueCutoff=0.2, readable=TRUE), error=function(e) NULL)
  if(!is.null(ego_cc) && nrow(ego_cc) > 0) {
    write.csv(as.data.frame(ego_cc), file.path(output_dir, "results", "enrichment", output_subdir,
              paste0(label, "_GO_CC.csv")), row.names=FALSE)
    p <- dotplot(ego_cc, showCategory=15) + ggtitle(paste("GO CC -", label))
    save_tiff_600dpi(p, file.path(output_dir, "results", "enrichment", output_subdir,
                    paste0(label, "_GO_CC_dotplot.tiff")), width=9, height=7)
  }
  
  # GO MF
  ego_mf <- tryCatch(enrichGO(gene=entrez_ids, OrgDb=org.Hs.eg.db,
                               ont="MF", pAdjustMethod="BH", pvalueCutoff=0.05,
                               qvalueCutoff=0.2, readable=TRUE), error=function(e) NULL)
  if(!is.null(ego_mf) && nrow(ego_mf) > 0) {
    write.csv(as.data.frame(ego_mf), file.path(output_dir, "results", "enrichment", output_subdir,
              paste0(label, "_GO_MF.csv")), row.names=FALSE)
  }
  
  # KEGG
  kk <- tryCatch(enrichKEGG(gene=entrez_ids, organism="hsa", pvalueCutoff=0.05,
                             pAdjustMethod="BH"), error=function(e) NULL)
  if(!is.null(kk) && nrow(kk) > 0) {
    write.csv(as.data.frame(kk), file.path(output_dir, "results", "enrichment", output_subdir,
              paste0(label, "_KEGG.csv")), row.names=FALSE)
    p <- dotplot(kk, showCategory=15) + ggtitle(paste("KEGG -", label))
    save_tiff_600dpi(p, file.path(output_dir, "results", "enrichment", output_subdir,
                    paste0(label, "_KEGG_dotplot.tiff")), width=9, height=7)
    # Pathview for top KEGG pathways (only if KEGG permission exists)
    if(nrow(kk) > 0) {
      top_pathways <- head(kk$ID, 3)
      for(pid in top_pathways) {
        tryCatch({
          pathview(gene.data=ranked_list, pathway.id=pid, species="hsa", 
                   out.suffix=paste0(label, "_", pid), kegg.native=TRUE,
                   key.pos="topright")
          # Move the PNG output
          png_files <- list.files(pattern=paste0(pid, ".*\\.png$"))
          for(pf in png_files) {
            file.rename(pf, file.path(output_dir, "results", "enrichment", output_subdir, pf))
          }
        }, error=function(e) cat("Pathview error for", pid, ":", e$message, "\n"))
      }
    }
  }
  
  # Reactome via enrichR
  dbs <- c("Reactome_2022", "KEGG_2021_Human")
  tryCatch({
    er <- enrichr(gene_list, dbs)
    for(db in names(er)) {
      if(nrow(er[[db]]) > 0) {
        write.csv(er[[db]], file.path(output_dir, "results", "enrichment", output_subdir,
                  paste0(label, "_", gsub(" ", "_", db), ".csv")), row.names=FALSE)
      }
    }
  }, error=function(e) cat("enrichR error:", e$message, "\n"))
  
  # Cnetplot
  if(!is.null(ego_bp) && nrow(ego_bp) > 0) {
    tiff(file.path(output_dir, "results", "enrichment", output_subdir,
         paste0(label, "_cnetplot.tiff")),
         width=12, height=10, units="in", res=600, compression="lzw", bg="white")
    print(cnetplot(ego_bp, showCategory=5))
    dev.off()
  }
}

# Run enrichment for all DEGs
dir.create(file.path(output_dir, "results", "enrichment", "all_DEGs"), showWarnings=FALSE, recursive=TRUE)
run_enrichment(rownames(deg_sig), "all_DEGs", "all_DEGs")

# Run enrichment for DE-NRGs
dir.create(file.path(output_dir, "results", "enrichment", "DE_NRGs"), showWarnings=FALSE, recursive=TRUE)
if(length(de_nrgs) >= 5) {
  run_enrichment(de_nrgs, "DE_NRGs", "DE_NRGs")
}

# GSEA using ranked gene list
cat("\n--- GSEA using ranked DEG list ---\n")
gsea_rank <- ranked_list[!is.na(names(ranked_list))]
gsea_rank <- gsea_rank[!duplicated(names(gsea_rank))]

# Convert to ENTREZ for GSEA
rank_df <- data.frame(SYMBOL=names(gsea_rank), Score=gsea_rank)
rank_entrez <- bitr(names(gsea_rank), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
rank_entrez <- rank_entrez[!duplicated(rank_entrez$SYMBOL), ]
ranked_entrez <- gsea_rank[rank_entrez$SYMBOL]
names(ranked_entrez) <- rank_entrez$ENTREZID
ranked_entrez <- sort(ranked_entrez, decreasing=TRUE)

gsea_go <- tryCatch(gseGO(geneList=ranked_entrez, OrgDb=org.Hs.eg.db,
                           ont="BP", minGSSize=10, maxGSSize=500,
                           pvalueCutoff=0.05, verbose=FALSE), error=function(e) NULL)
if(!is.null(gsea_go) && nrow(gsea_go) > 0) {
  write.csv(as.data.frame(gsea_go), file.path(output_dir, "results", "enrichment", "GSEA_GO_BP.csv"), row.names=FALSE)
  p <- dotplot(gsea_go, showCategory=15) + ggtitle("GSEA - GO BP")
  save_tiff_600dpi(p, file.path(output_dir, "results", "enrichment", "GSEA_GO_BP_dotplot.tiff"), width=10, height=8)
}

gsea_kegg <- tryCatch(gseKEGG(geneList=ranked_entrez, organism="hsa",
                               minGSSize=10, maxGSSize=500,
                               pvalueCutoff=0.05, verbose=FALSE), error=function(e) NULL)
if(!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
  write.csv(as.data.frame(gsea_kegg), file.path(output_dir, "results", "enrichment", "GSEA_KEGG.csv"), row.names=FALSE)
}

cat("\n\nENRICHMENT ANALYSIS COMPLETE\n")
sink()
