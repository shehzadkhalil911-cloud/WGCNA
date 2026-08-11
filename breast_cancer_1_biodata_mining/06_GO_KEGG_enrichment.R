###############################################################################
# Script 06: GO and KEGG Enrichment Analysis
###############################################################################

set.seed(123)

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(viridis)
library(pathview)
library(circlize)
library(patchwork)

base_dir <- normalizePath("breast_cancer_integrated_analysis")
data_processed_dir <- file.path(base_dir, "data_processed")
fig_dir <- file.path(base_dir, "results", "figures", "Enrichment")
tables_dir <- file.path(base_dir, "results", "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helper: Save 600 dpi ---
save_publication_figure <- function(plot, filename, width = 8, height = 7) {
  ggsave(file.path(fig_dir, paste0(filename, ".tiff")), plot = plot, 
         width = width, height = height, dpi = 600, compression = "lzw")
  ggsave(file.path(fig_dir, paste0(filename, ".pdf")), plot = plot, 
         width = width, height = height, dpi = 600)
}

theme_publication <- theme_classic() +
  theme(
    text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

# --- Load DEG data ---
cat("Loading DEG results...\n")
deg_table <- readRDS(file.path(data_processed_dir, "limma_deg_table.rds"))
deg_genes <- readRDS(file.path(data_processed_dir, "deg_genes.rds"))

cat("Total DEGs:", length(deg_genes), "\n")

if (length(deg_genes) == 0) {
  stop("No DEG genes found. Cannot perform enrichment analysis.")
}

# --- Convert gene symbols to Entrez IDs ---
entrez_ids <- bitr(deg_genes, fromType = "SYMBOL", toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)
cat("Converted", nrow(entrez_ids), "genes to Entrez IDs\n")

if (nrow(entrez_ids) == 0) {
  stop("No Entrez IDs could be mapped.")
}

entrez_list <- entrez_ids$ENTREZID

# --- Background genes (all genes tested) ---
all_entrez <- bitr(rownames(deg_table), fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)
all_entrez_list <- all_entrez$ENTREZID

# --- GO Enrichment Analysis ---
cat("\nPerforming GO enrichment analysis...\n")

go_bp <- enrichGO(gene = entrez_list, universe = all_entrez_list,
                  OrgDb = org.Hs.eg.db, ont = "BP", 
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                  readable = TRUE)
go_cc <- enrichGO(gene = entrez_list, universe = all_entrez_list,
                  OrgDb = org.Hs.eg.db, ont = "CC",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                  readable = TRUE)
go_mf <- enrichGO(gene = entrez_list, universe = all_entrez_list,
                  OrgDb = org.Hs.eg.db, ont = "MF",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                  readable = TRUE)

# --- KEGG Enrichment ---
cat("Performing KEGG enrichment analysis...\n")
kegg_enrich <- enrichKEGG(gene = entrez_list, universe = all_entrez_list,
                          organism = "hsa", pAdjustMethod = "BH",
                          pvalueCutoff = 0.05, qvalueCutoff = 0.2)

# --- Save enrichment results ---
enrich_results <- list(GO_BP = go_bp, GO_CC = go_cc, GO_MF = go_mf, KEGG = kegg_enrich)
saveRDS(enrich_results, file.path(data_processed_dir, "enrichment_results.rds"))

# Write tables
if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
  write.csv(go_bp@result, file.path(tables_dir, "GO_BP_enrichment.csv"), row.names = FALSE)
}
if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
  write.csv(go_cc@result, file.path(tables_dir, "GO_CC_enrichment.csv"), row.names = FALSE)
}
if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
  write.csv(go_mf@result, file.path(tables_dir, "GO_MF_enrichment.csv"), row.names = FALSE)
}
if (!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
  write.csv(kegg_enrich@result, file.path(tables_dir, "KEGG_enrichment.csv"), row.names = FALSE)
}

cat("Enrichment results saved to tables/\n")

# --- Generate Figures ---

# A. GO/KEGG dot plot
cat("\nGenerating enrichment figures...\n")

if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
  p_go_dot <- dotplot(go_bp, showCategory = 20, title = "GO Biological Process") +
    theme_publication
  save_publication_figure(p_go_dot, "Figure_3A_GO_dotplot")
}

if (!is.null(go_cc) && nrow(go_cc@result) > 0) {
  p_cc_dot <- dotplot(go_cc, showCategory = 15, title = "GO Cellular Component") +
    theme_publication
  save_publication_figure(p_cc_dot, "Figure_3A_GO_CC_dotplot", width = 8, height = 6)
}

if (!is.null(go_mf) && nrow(go_mf@result) > 0) {
  p_mf_dot <- dotplot(go_mf, showCategory = 15, title = "GO Molecular Function") +
    theme_publication
  save_publication_figure(p_mf_dot, "Figure_3A_GO_MF_dotplot", width = 8, height = 6)
}

# B. GO circle plot (barplot)
if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
  p_go_bar <- barplot(go_bp, showCategory = 20, title = "GO Biological Process Enrichment") +
    theme_publication
  save_publication_figure(p_go_bar, "Figure_3B_GO_barplot")
}

# C. KEGG bar plot
if (!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
  p_kegg_bar <- barplot(kegg_enrich, showCategory = 20, 
                         title = "KEGG Pathway Enrichment") +
    theme_publication
  save_publication_figure(p_kegg_bar, "Figure_3C_KEGG_barplot")
  
  # KEGG dot plot
  p_kegg_dot <- dotplot(kegg_enrich, showCategory = 20, 
                         title = "KEGG Pathway Enrichment") +
    theme_publication
  save_publication_figure(p_kegg_dot, "Figure_3C_KEGG_dotplot")
}

# D. cnetplot (pathway-network)
cat("Generating network plots...\n")
if (!is.null(go_bp) && sum(go_bp@result$p.adjust < 0.05) > 0) {
  tryCatch({
    p_cnet <- cnetplot(go_bp, showCategory = 5) +
      theme_publication
    save_publication_figure(p_cnet, "Figure_3D_GO_cnetplot", width = 12, height = 10)
  }, error = function(e) cat("cnetplot error:", e$message, "\n"))
}

# E. emapplot
if (!is.null(go_bp) && sum(go_bp@result$p.adjust < 0.05) > 0) {
  tryCatch({
    go_bp_pairwise <- pairwise_termsim(go_bp)
    p_emap <- emapplot(go_bp_pairwise, showCategory = 30) +
      theme_publication
    save_publication_figure(p_emap, "Figure_3D_GO_emapplot", width = 10, height = 8)
  }, error = function(e) cat("emapplot error:", e$message, "\n"))
}

# F. ridgeplot
if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
  tryCatch({
    p_ridge <- ridgeplot(go_bp, showCategory = 20) +
      theme_publication
    save_publication_figure(p_ridge, "Figure_3D_GO_ridgeplot", width = 10, height = 8)
  }, error = function(e) cat("ridgeplot error:", e$message, "\n"))
}

# G. KEGG pathway maps using pathview
cat("\nGenerating KEGG pathway maps...\n")
if (!is.null(kegg_enrich) && nrow(kegg_enrich@result) > 0) {
  sig_kegg <- kegg_enrich@result[kegg_enrich@result$p.adjust < 0.05, ]
  
  if (nrow(sig_kegg) > 0) {
    top_kegg <- head(sig_kegg$ID, 5)
    
    # Prepare fold change data
    fc_data <- deg_table$logFC
    names(fc_data) <- rownames(deg_table)
    
    # Map to Entrez
    gene_map <- bitr(names(fc_data), fromType = "SYMBOL", toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
    fc_entrez <- fc_data[gene_map$SYMBOL]
    names(fc_entrez) <- gene_map$ENTREZID
    
    for (kegg_id in top_kegg) {
      cat("  Processing pathway:", kegg_id, "\n")
      tryCatch({
        pathview(gene.data = fc_entrez, pathway.id = kegg_id,
                 species = "hsa", out.suffix = gsub(":", "_", kegg_id),
                 kegg.native = TRUE, same.layer = TRUE,
                 plot.col.key = TRUE)
        # Move files to figures directory
        pv_files <- list.files(pattern = paste0(kegg_id, ".*\\.png"), path = ".", full.names = FALSE)
        for (f in pv_files) {
          file.rename(f, file.path(fig_dir, f))
        }
      }, error = function(e) {
        cat("    pathview error for", kegg_id, ":", e$message, "\n")
      })
    }
  }
}

# --- Create combined GO plot ---
cat("\nCreating combined enrichment summary...\n")
if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
  # Combine multiple GO plots
  p_go_dot <- dotplot(go_bp, showCategory = 15, title = "GO Biological Process") +
    theme_publication
  p_go_bar <- barplot(go_bp, showCategory = 15, title = "GO BP Enrichment") +
    theme_publication
  
  combined_go <- p_go_dot + p_go_bar + 
    plot_layout(ncol = 2) +
    plot_annotation(title = "Gene Ontology Enrichment Analysis",
                    theme = theme(plot.title = element_text(size = 14, face = "bold",
                                                            hjust = 0.5)))
  save_publication_figure(combined_go, "Figure_3AB_GO_combined", width = 14, height = 8)
}

cat("\n========================================================\n")
cat("GO/KEGG ENRICHMENT ANALYSIS COMPLETE\n")
cat("========================================================\n")
cat("GO BP terms enriched:", ifelse(!is.null(go_bp), sum(go_bp@result$p.adjust < 0.05), 0), "\n")
cat("GO CC terms enriched:", ifelse(!is.null(go_cc), sum(go_cc@result$p.adjust < 0.05), 0), "\n")
cat("GO MF terms enriched:", ifelse(!is.null(go_mf), sum(go_mf@result$p.adjust < 0.05), 0), "\n")
cat("KEGG pathways enriched:", ifelse(!is.null(kegg_enrich), sum(kegg_enrich@result$p.adjust < 0.05), 0), "\n")
cat("Results saved to results/figures/Enrichment/ and results/tables/\n")
