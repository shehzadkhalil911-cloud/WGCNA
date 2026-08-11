###############################################################################
# Script 10: GSEA and GSVA Analysis
###############################################################################

set.seed(123)

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(msigdbr)
library(GSVA)
library(GSEABase)
library(ggplot2)
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(pheatmap)

base_dir <- normalizePath("breast_cancer_integrated_analysis")
data_processed_dir <- file.path(base_dir, "data_processed")
fig_dir <- file.path(base_dir, "results", "figures", "GSEA_GSVA")
tables_dir <- file.path(base_dir, "results", "tables")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helper: Save 600 dpi ---
save_publication_figure <- function(plot, filename, width = 9, height = 7) {
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
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  )

# --- Load data ---
cat("Loading data...\n")
exprs_after <- readRDS(file.path(data_processed_dir, "merged_exprs_after_batch.rds"))
pheno <- readRDS(file.path(data_processed_dir, "merged_pheno_with_batch.rds"))
deg_table <- readRDS(file.path(data_processed_dir, "limma_deg_table.rds"))

# ================================================================
# PART 1: GSEA
# ================================================================
cat("\n========================================================\n")
cat("GENE SET ENRICHMENT ANALYSIS (GSEA)\n")
cat("========================================================\n")

# --- Prepare ranked gene list ---
gene_list <- deg_table[["logFC"]]
names(gene_list) <- rownames(deg_table)
# Remove NA values and sort
gene_list <- gene_list[!is.na(gene_list)]
gene_list <- sort(gene_list, decreasing = TRUE)

cat("Ranked gene list length:", length(gene_list), "\n")

# Convert gene symbols to Entrez for GSEA
gene_map <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)
gene_map <- gene_map[!duplicated(gene_map$SYMBOL), ]
rownames(gene_map) <- gene_map$SYMBOL

# For GSEA we need gene list names to match gene set gene symbols
# Since hallmark/go use gene_symbol and KEGG uses entrez, keep both
# Use gene symbols as primary (match with hallmark/go)
gene_list_symbols <- gene_list
names(gene_list_symbols) <- names(gene_list)
gene_list_symbols <- gene_list_symbols[!is.na(names(gene_list_symbols))]

# Also create entrez version for KEGG
gene_list_entrez <- gene_list[gene_map$SYMBOL]
names(gene_list_entrez) <- gene_map$ENTREZID
gene_list_entrez <- gene_list_entrez[!is.na(names(gene_list_entrez))]

cat("Ranked gene list with gene symbols:", length(gene_list_symbols), "\n")
cat("Ranked gene list with Entrez IDs:", length(gene_list_entrez), "\n")

# --- Load MSigDB gene sets ---
cat("Loading MSigDB gene sets...\n")

# Hallmarks
hallmark_set <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_genes <- split(hallmark_set$gene_symbol, hallmark_set$gs_name)

# GO BP
go_set_all <- msigdbr(species = "Homo sapiens", collection = "C5")
go_set <- go_set_all[go_set_all$gs_subcollection == "GO:BP", ]
go_genes <- split(go_set$gene_symbol, go_set$gs_name)

# KEGG 
kegg_set_all <- msigdbr(species = "Homo sapiens", collection = "C2")
kegg_set <- kegg_set_all[kegg_set_all$gs_subcollection == "CP:KEGG", ]
kegg_genes <- split(kegg_set$gene_symbol, kegg_set$gs_name)

# --- Run GSEA ---
run_gsea <- function(gene_list, gene_sets, min_size = 15, max_size = 500) {
  set.seed(123)
  # Clean gene sets: remove NA names, empty sets, and duplicates
  gene_sets <- gene_sets[!is.na(names(gene_sets))]
  gene_sets <- gene_sets[sapply(gene_sets, length) > 0]
  gene_sets <- gene_sets[!duplicated(names(gene_sets))]
  
  term2gene <- data.frame(
    term = rep(names(gene_sets), sapply(gene_sets, length)),
    gene = unlist(unname(gene_sets)),
    stringsAsFactors = FALSE
  )
  # Remove any rows with NA
  term2gene <- term2gene[complete.cases(term2gene), ]
  
  res <- GSEA(gene_list, 
              TERM2GENE = term2gene,
              minGSSize = min_size,
              maxGSSize = max_size,
              pvalueCutoff = 0.05,
              pAdjustMethod = "BH",
              verbose = FALSE)
  res
}

gsea_hallmark <- run_gsea(gene_list_symbols, hallmark_genes)
gsea_go <- run_gsea(gene_list_symbols, go_genes)
gsea_kegg <- run_gsea(gene_list_symbols, kegg_genes)

# --- Save GSEA results ---
gsea_results <- list(Hallmark = gsea_hallmark, GO_BP = gsea_go, KEGG = gsea_kegg)
saveRDS(gsea_results, file.path(data_processed_dir, "gsea_results.rds"))

if (!is.null(gsea_hallmark) && nrow(gsea_hallmark@result) > 0) {
  write.csv(gsea_hallmark@result, file.path(tables_dir, "GSEA_hallmark.csv"), row.names = FALSE)
}
if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
  write.csv(gsea_go@result, file.path(tables_dir, "GSEA_GO_BP.csv"), row.names = FALSE)
}
if (!is.null(gsea_kegg) && nrow(gsea_kegg@result) > 0) {
  write.csv(gsea_kegg@result, file.path(tables_dir, "GSEA_KEGG.csv"), row.names = FALSE)
}

# --- GSEA dot plot ---
cat("\nGenerating GSEA figures...\n")

if (!is.null(gsea_hallmark) && nrow(gsea_hallmark@result) > 0) {
  p_gsea_dot <- dotplot(gsea_hallmark, showCategory = 20, 
                         split = ".sign") + 
    facet_grid(. ~ .sign) +
    labs(title = "GSEA: Hallmark Gene Sets") +
    theme_publication
  save_publication_figure(p_gsea_dot, "Figure_8A_GSEA_dotplot_hallmark")
}

if (!is.null(gsea_go) && nrow(gsea_go@result) > 0) {
  p_go_dot <- dotplot(gsea_go, showCategory = 15, 
                       split = ".sign") + 
    facet_grid(. ~ .sign) +
    labs(title = "GSEA: GO Biological Process") +
    theme_publication
  save_publication_figure(p_go_dot, "Figure_8A_GSEA_dotplot_GO", width = 10, height = 8)
}

# --- GSEA ridge plot ---
if (!is.null(gsea_hallmark) && nrow(gsea_hallmark@result) > 0) {
  tryCatch({
    p_ridge <- ridgeplot(gsea_hallmark, showCategory = 20) +
      labs(title = "GSEA: Hallmark Gene Sets (Ridge Plot)") +
      theme_publication
    save_publication_figure(p_ridge, "Figure_8A_GSEA_ridgeplot", width = 10, height = 8)
  }, error = function(e) cat("  ridgeplot error:", e$message, "\n"))
}

# --- GSEA enrichment plots for top pathways ---
if (!is.null(gsea_hallmark) && nrow(gsea_hallmark@result) > 0) {
  sig_hallmark <- gsea_hallmark@result[gsea_hallmark@result$p.adjust < 0.05, ]
  if (nrow(sig_hallmark) > 0) {
    top_paths <- head(sig_hallmark$ID, 5)
    for (i in seq_along(top_paths)) {
      tryCatch({
        p_enrich <- enrichplot::gseaplot2(gsea_hallmark, geneSetID = top_paths[i],
                                           title = top_paths[i])
        save_publication_figure(p_enrich, 
                                paste0("Figure_8A_GSEA_enrichment_", gsub("/", "_", top_paths[i])),
                                width = 8, height = 6)
      }, error = function(e) cat("  gseaplot error:", e$message, "\n"))
    }
  }
}

# ================================================================
# PART 2: GSVA
# ================================================================
cat("\n========================================================\n")
cat("GENE SET VARIATION ANALYSIS (GSVA)\n")
cat("========================================================\n")

# --- Prepare GSVA input ---
cat("Preparing gene sets for GSVA...\n")

# Convert gene expression matrix to matrix format
gsva_input <- as.matrix(exprs_after)

# Use Hallmark gene sets for GSVA (convert to list of gene symbols)
hallmark_genes_list <- split(hallmark_set$gene_symbol, hallmark_set$gs_name)

# Filter to only include genes present in expression matrix
hallmark_gsva <- lapply(hallmark_genes_list, function(genes) {
  intersect(genes, rownames(gsva_input))
})
hallmark_gsva <- hallmark_gsva[sapply(hallmark_gsva, length) >= 10]
cat("Hallmark gene sets for GSVA:", length(hallmark_gsva), "\n")

# --- Run GSVA ---
cat("Running GSVA...\n")
gsva_results <- gsva(gsva_input, hallmark_gsva, 
                      method = "gsva",
                      kcdf = "Gaussian",
                      min.sz = 10,
                      max.sz = 500,
                      verbose = TRUE)

cat("GSVA result dimensions:", dim(gsva_results), "\n")

# --- Compare GSVA scores: Tumor vs Normal ---
gsva_scores <- as.data.frame(t(gsva_results))
gsva_scores$Group <- pheno$group

# Differential GSVA analysis
gsva_diff <- data.frame()
for (pathway in rownames(gsva_results)) {
  tumor_scores <- gsva_results[pathway, pheno$group == "Tumor"]
  normal_scores <- gsva_results[pathway, pheno$group == "Normal"]
  
  t_test <- tryCatch(t.test(tumor_scores, normal_scores), error = function(e) NULL)
  
  if (!is.null(t_test)) {
    gsva_diff <- rbind(gsva_diff, data.frame(
      Pathway = pathway,
      Mean_Tumor = mean(tumor_scores),
      Mean_Normal = mean(normal_scores),
      Mean_Diff = mean(tumor_scores) - mean(normal_scores),
      P_value = t_test$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

if (nrow(gsva_diff) > 0) {
  gsva_diff$Adj_P <- p.adjust(gsva_diff$P_value, method = "BH")
  gsva_diff <- gsva_diff[order(gsva_diff$P_value), ]
  
  write.csv(gsva_diff, file.path(tables_dir, "GSVA_differential_pathways.csv"), 
            row.names = FALSE)
  cat("Differential GSVA analysis complete.\n")
}

# --- GSVA heatmap ---
cat("Generating GSVA heatmap...\n")

# Select top variable pathways
pathway_var <- apply(gsva_results, 1, var)
top_pathways <- names(sort(pathway_var, decreasing = TRUE)[1:min(25, length(pathway_var))])

annotation_col <- data.frame(
  Group = pheno$group,
  row.names = colnames(exprs_after)
)
annotation_colors <- list(
  Group = c(Normal = "#4daf4a", Tumor = "#e41a1c")
)

pheatmap(
  gsva_results[top_pathways, ],
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  main = "GSVA: Hallmark Pathway Scores",
  filename = file.path(fig_dir, "Figure_8B_GSVA_heatmap.tiff"),
  width = 10,
  height = 8,
  resolution = 600
)

pheatmap(
  gsva_results[top_pathways, ],
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  main = "GSVA: Hallmark Pathway Scores",
  filename = file.path(fig_dir, "Figure_8B_GSVA_heatmap.pdf"),
  width = 10,
  height = 8
)

# --- GSVA box plots for top pathways ---
cat("Generating GSVA box plots...\n")

if (nrow(gsva_diff) > 0) {
  top_diff_paths <- head(gsva_diff$Pathway, 9)
  
  gsva_plot_data <- data.frame()
  for (path in top_diff_paths) {
    path_df <- data.frame(
      Score = gsva_results[path, ],
      Pathway = gsub("HALLMARK_", "", path),
      Group = pheno$group,
      stringsAsFactors = FALSE
    )
    gsva_plot_data <- rbind(gsva_plot_data, path_df)
  }
  
  p_gsva_box <- ggplot(gsva_plot_data, aes(x = Group, y = Score, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.4) +
    facet_wrap(~ Pathway, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c("Normal" = "#4daf4a", "Tumor" = "#e41a1c")) +
    labs(title = "GSVA Pathway Scores: Tumor vs Normal",
         x = "", y = "GSVA Score") +
    theme_publication +
    theme(strip.text = element_text(size = 8, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  
  save_publication_figure(p_gsva_box, "Figure_8C_GSVA_boxplot", 
                          width = 12, height = 10)
}

# --- Save GSVA results ---
saveRDS(gsva_results, file.path(data_processed_dir, "gsva_results.rds"))
saveRDS(gsva_diff, file.path(data_processed_dir, "gsva_differential.rds"))

cat("\n========================================================\n")
cat("GSEA/GSVA ANALYSIS COMPLETE\n")
cat("========================================================\n")
cat("GSEA significant hallmark pathways:", 
    ifelse(!is.null(gsea_hallmark), sum(gsea_hallmark@result$p.adjust < 0.05), 0), "\n")
cat("GSVA pathways tested:", nrow(gsva_results), "\n")
cat("GSVA significantly different:", sum(gsva_diff$Adj_P < 0.05, na.rm = TRUE), "\n")
cat("Results saved to results/figures/GSEA_GSVA/ and results/tables/\n")
