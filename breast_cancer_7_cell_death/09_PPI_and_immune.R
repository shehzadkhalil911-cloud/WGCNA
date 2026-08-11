# =============================================================================
# Script 09: PPI Network and Immune Infiltration
# Breast Cancer 7 Cell Death Crosstalk
# =============================================================================

source(file.path("..", "config", "analysis_config.R"))

library(tidyverse)
library(igraph)
library(ggraph)
library(STRINGdb)
library(xCell)
library(ComplexHeatmap)
library(circlize)

log_message("Starting PPI and immune infiltration analysis.")

# --- 1. LOAD DATA ---
shared_info <- readRDS(file.path(OUTPUT$shared, "shared_gene_classification.rds"))
deg_sig <- read.csv(file.path(OUTPUT$de, "discovery_DEG_significant.csv"), row.names = 1)
combat_matrix <- readRDS(file.path(OUTPUT$batch_correction, "combat_corrected_matrix.rds"))
discovery_metadata <- readRDS(file.path(OUTPUT$batch_correction, "discovery_metadata.rds"))
scores <- readRDS(file.path(OUTPUT$activity, "program_scores.rds"))

# --- 2. PPI NETWORK USING STRINGdb ---
build_ppi <- function(genes, name, output_dir) {
  if (length(genes) < 3) {
    log_message(sprintf("PPI for %s: too few genes (%d)", name, length(genes)), "WARNING")
    return(NULL)
  }
  
  # Map to Entrez
  entrez_ids <- tryCatch({
    mapIds(org.Hs.eg.db, keys = genes, keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
  }, error = function(e) NULL)
  entrez_ids <- as.character(entrez_ids[!is.na(entrez_ids)])
  
  if (length(entrez_ids) < 3) {
    log_message(sprintf("PPI for %s: too few Entrez IDs (%d)", name, length(entrez_ids)), "WARNING")
    return(NULL)
  }
  
  # Use STRINGdb
  string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 400, 
                             input_directory = file.path(OUTPUT$ppi, "string_db"))
  
  mapped <- string_db$map(data.frame(gene = entrez_ids), "gene", removeUnmappedRows = TRUE)
  
  if (nrow(mapped) < 3) {
    log_message(sprintf("PPI for %s: too few STRING mappings (%d)", name, nrow(mapped)), "WARNING")
    return(NULL)
  }
  
  # Get interactions
  interactions <- string_db$get_interactions(mapped$STRING_id)
  
  if (is.null(interactions) || nrow(interactions) == 0) {
    log_message(sprintf("PPI for %s: no interactions found", name), "WARNING")
    return(NULL)
  }
  
  # Create graph
  nodes <- data.frame(name = unique(c(interactions$from, interactions$to)), stringsAsFactors = FALSE)
  # Map back to gene symbols
  node_map <- setNames(mapped$gene, mapped$STRING_id)
  nodes$gene_symbol <- node_map[nodes$name]
  nodes$gene_symbol[is.na(nodes$gene_symbol)] <- nodes$name[is.na(nodes$gene_symbol)]
  
  edges <- data.frame(
    from = interactions$from,
    to = interactions$to,
    weight = interactions$score / 1000,
    stringsAsFactors = FALSE
  )
  
  g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  
  # Calculate centrality
  V(g)$degree <- degree(g, normalized = TRUE)
  V(g)$betweenness <- betweenness(g, normalized = TRUE, weights = NA)
  V(g)$closeness <- closeness(g, normalized = TRUE, weights = NA)
  V(g)$eigen <- eigen_centrality(g, weights = NA)$vector
  
  # Communities
  communities <- cluster_walktrap(g)
  V(g)$community <- membership(communities)
  
  # Save
  node_metrics <- data.frame(
    STRING_id = V(g)$name,
    gene_symbol = V(g)$gene_symbol,
    degree = V(g)$degree,
    betweenness = V(g)$betweenness,
    closeness = V(g)$closeness,
    eigenvector = V(g)$eigen,
    community = V(g)$community,
    stringsAsFactors = FALSE
  )
  write.csv(node_metrics, file.path(output_dir, sprintf("%s_node_metrics.csv", name)), row.names = FALSE)
  write.csv(edges, file.path(output_dir, sprintf("%s_edge_list.csv", name)), row.names = FALSE)
  saveRDS(g, file.path(output_dir, sprintf("%s_igraph.rds", name)))
  
  # Network plot
  p_net <- ggraph(g, layout = "kk") +
    geom_edge_link(aes(alpha = weight), color = "grey60", show.legend = FALSE) +
    geom_node_point(aes(size = degree, color = as.factor(community)), alpha = 0.8) +
    geom_node_text(aes(label = gene_symbol), repel = TRUE, size = 2.5) +
    scale_size_continuous(range = c(1, 6)) +
    labs(title = sprintf("%s PPI Network", name)) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
  
  save_tiff_600dpi(p_net, file.path(output_dir, sprintf("%s_network.tiff", name)), width = 10, height = 8)
  
  log_message(sprintf("PPI for %s: %d nodes, %d edges", name, nrow(node_metrics), nrow(edges)))
  
  return(list(graph = g, nodes = node_metrics, edges = edges))
}

# Build PPI for each programme
ppi_results <- list()

# DE-FRGs
frg_genes <- shared_info$Gene[shared_info$DE_FRG]
if (length(frg_genes) >= 3) {
  ppi_results$DE_FRG <- build_ppi(frg_genes, "DE_FRG", OUTPUT$ppi)
}

# DE-NRGs
nrg_genes <- shared_info$Gene[shared_info$DE_NRG]
if (length(nrg_genes) >= 3) {
  ppi_results$DE_NRG <- build_ppi(nrg_genes, "DE_NRG", OUTPUT$ppi)
}

# DE-CRGs
crg_genes <- shared_info$Gene[shared_info$DE_CRG]
if (length(crg_genes) >= 3) {
  ppi_results$DE_CRG <- build_ppi(crg_genes, "DE_CRG", OUTPUT$ppi)
}

# Shared genes
shared_genes <- shared_info$Gene[shared_info$Category %in% c("All Three", "Ferroptosis + Necroptosis", 
                                                              "Ferroptosis + Cuproptosis", "Necroptosis + Cuproptosis")]
if (length(shared_genes) >= 3) {
  ppi_results$Shared <- build_ppi(shared_genes, "Shared", OUTPUT$ppi)
}

# All DE-regulated cell death genes
all_de_rcd <- unique(c(frg_genes, nrg_genes, crg_genes))
if (length(all_de_rcd) >= 3) {
  ppi_results$All_RCD <- build_ppi(all_de_rcd, "All_RCD", OUTPUT$ppi)
}

# --- 3. IMMUNE INFILTRATION (xCell) ---
log_message("Running xCell immune deconvolution...")

# xCell needs expression matrix with gene symbols as rownames
xcell_matrix <- combat_matrix

# Run xCell
xcell_scores <- tryCatch({
  xCell::xCellAnalysis(xcell_matrix, 
                        rnaseq = FALSE, 
                        file.name = file.path(OUTPUT$immune, "xcell_results.txt"),
                        parallel.sz = 1)
}, error = function(e) {
  log_message(sprintf("xCell failed: %s", e$message), "WARNING")
  return(NULL)
})

if (!is.null(xcell_scores)) {
  saveRDS(xcell_scores, file.path(OUTPUT$immune, "xcell_scores.rds"))
  write.csv(xcell_scores, file.path(OUTPUT$immune, "xcell_scores.csv"))
  
  log_message(sprintf("xCell scores: %d cell types x %d samples", nrow(xcell_scores), ncol(xcell_scores)))
  
  # --- 4. IMMUNE CELL COMPARISON: Tumour vs Normal ---
  immune_df <- as.data.frame(t(xcell_scores))
  immune_df$Sample <- rownames(immune_df)
  immune_df$Group <- discovery_metadata$group[match(immune_df$Sample, discovery_metadata$merged_id)]
  immune_df$Dataset <- discovery_metadata$dataset[match(immune_df$Sample, discovery_metadata$merged_id)]
  
  immune_long <- immune_df %>%
    pivot_longer(cols = -c(Sample, Group, Dataset), names_to = "CellType", values_to = "Score")
  
  # Boxplots
  immune_cell_types <- colnames(immune_df)[!colnames(immune_df) %in% c("Sample", "Group", "Dataset")]
  
  p_immune_boxes <- ggplot(immune_long, aes(x = CellType, y = Score, fill = Group)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.7, width = 0.6) +
    scale_fill_manual(values = c("Normal" = "#0072B2", "Tumor" = "#D55E00")) +
    labs(title = "Immune Cell Infiltration (xCell)", 
         subtitle = "Tumour vs Normal", y = "Enrichment Score") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))
  save_tiff_600dpi(p_immune_boxes, file.path(OUTPUT$immune, "immune_boxplots.tiff"), width = 14, height = 8)
  
  # Differential immune infiltration
  immune_diff <- data.frame()
  for (ct in immune_cell_types) {
    tumor_scores <- immune_df[[ct]][immune_df$Group == "Tumor"]
    normal_scores <- immune_df[[ct]][immune_df$Group == "Normal"]
    
    if (length(tumor_scores) > 3 && length(normal_scores) > 3) {
      wt <- wilcox.test(tumor_scores, normal_scores)
      immune_diff <- rbind(immune_diff, data.frame(
        CellType = ct,
        Tumor_Mean = mean(tumor_scores, na.rm = TRUE),
        Normal_Mean = mean(normal_scores, na.rm = TRUE),
        log2FC = log2(mean(tumor_scores, na.rm = TRUE) / (mean(normal_scores, na.rm = TRUE) + 0.001)),
        P_value = wt$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
  immune_diff$FDR <- p.adjust(immune_diff$P_value, method = "BH")
  write.csv(immune_diff, file.path(OUTPUT$immune, "differential_immune.csv"), row.names = FALSE)
  
  # Immune infiltration heatmap
  ht_immune <- Heatmap(xcell_scores,
                       name = "Score",
                       col = colorRamp2(c(0, 0.5, 1), c("#2166AC", "white", "#B2182B")),
                       top_annotation = HeatmapAnnotation(
                         Group = discovery_metadata$group,
                         col = list(Group = c(Normal = "#0072B2", Tumor = "#D55E00"))
                       ),
                       show_column_names = FALSE,
                       show_row_names = TRUE,
                       row_names_gp = gpar(fontsize = 8),
                       column_title = "Immune Infiltration")
  
  open_tiff(file.path(OUTPUT$immune, "immune_heatmap.tiff"), width = 10, height = 6)
  draw(ht_immune)
  dev.off()
  
  # --- 5. SCORE-IMMUNE CORRELATION ---
  score_immune_cor <- matrix(0, nrow = 3, ncol = length(immune_cell_types))
  rownames(score_immune_cor) <- c("Ferroptosis", "Necroptosis", "Cuproptosis")
  colnames(score_immune_cor) <- immune_cell_types
  score_immune_pval <- score_immune_cor
  
  for (prog in c("Ferroptosis", "Necroptosis", "Cuproptosis")) {
    for (ct in immune_cell_types) {
      ct_scores <- immune_df[[ct]]
      prog_scores <- scores[prog, immune_df$Sample]
      valid <- !is.na(ct_scores) & !is.na(prog_scores)
      if (sum(valid) > 5) {
        ct_res <- cor.test(ct_scores[valid], prog_scores[valid], method = "spearman")
        score_immune_cor[prog, ct] <- ct_res$estimate
        score_immune_pval[prog, ct] <- ct_res$p.value
      }
    }
  }
  
  write.csv(score_immune_cor, file.path(OUTPUT$immune, "score_immune_correlation.csv"))
  
  # Score-immune heatmap
  ht_si <- Heatmap(score_immune_cor,
                   name = "Spearman r",
                   col = colorRamp2(c(-0.5, 0, 0.5), c("#2166AC", "white", "#B2182B")),
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     pval <- score_immune_pval[i, j]
                     sig <- ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", 
                                   ifelse(pval < 0.05, "*", "")))
                     grid.text(sprintf("%.2f%s", score_immune_cor[i, j], sig), 
                               x, y, gp = gpar(fontsize = 8))
                   },
                   rect_gp = gpar(col = "grey60", lwd = 0.5),
                   row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                   column_names_gp = gpar(fontsize = 8, angle = 45, hjust = 1),
                   column_title = "Score-Immune Correlation",
                   row_title = "Program")
  
  open_tiff(file.path(OUTPUT$immune, "score_immune_heatmap.tiff"), width = 10, height = 4)
  draw(ht_si)
  dev.off()
  
  # --- 6. GENE-IMMUNE CORRELATION ---
  all_de_rcd_genes <- unique(c(frg_genes, nrg_genes, crg_genes))
  all_de_rcd_genes <- intersect(all_de_rcd_genes, rownames(combat_matrix))
  
  if (length(all_de_rcd_genes) > 3 && length(immune_cell_types) > 3) {
    gene_immune_cor <- cor(t(combat_matrix[all_de_rcd_genes, ]), 
                           as.matrix(immune_df[, immune_cell_types]),
                           method = "spearman")
    
    write.csv(gene_immune_cor, file.path(OUTPUT$immune, "gene_immune_correlation.csv"))
    
    if (nrow(gene_immune_cor) > 2) {
      ht_gi <- Heatmap(gene_immune_cor,
                       name = "Spearman r",
                       col = colorRamp2(c(-0.5, 0, 0.5), c("#2166AC", "white", "#B2182B")),
                       show_row_names = FALSE,
                       show_column_names = TRUE,
                       column_names_gp = gpar(fontsize = 8, angle = 45, hjust = 1),
                       column_title = "Gene-Immune Correlation",
                       row_title = "DE RCD Genes",
                       row_km = 3)
      
      open_tiff(file.path(OUTPUT$immune, "gene_immune_heatmap.tiff"), width = 10, height = 8)
      draw(ht_gi)
      dev.off()
    }
  }
  
  # --- 7. MODULE-IMMUNE CORRELATION ---
  MEs <- readRDS(file.path(OUTPUT$wgcna, "module_eigengenes.rds"))
  MEs <- MEs[immune_df$Sample, ]
  
  module_immune_cor <- cor(MEs, as.matrix(immune_df[, immune_cell_types]), method = "spearman")
  
  write.csv(module_immune_cor, file.path(OUTPUT$immune, "module_immune_correlation.csv"))
  
  ht_mi <- Heatmap(module_immune_cor,
                   name = "Spearman r",
                   col = colorRamp2(c(-0.5, 0, 0.5), c("#2166AC", "white", "#B2182B")),
                   rect_gp = gpar(col = "grey60", lwd = 0.5),
                   row_names_gp = gpar(fontsize = 8),
                   column_names_gp = gpar(fontsize = 8, angle = 45, hjust = 1),
                   column_title = "Module-Immune Correlation",
                   row_title = "Module")
  
  open_tiff(file.path(OUTPUT$immune, "module_immune_heatmap.tiff"), width = 10, height = 6)
  draw(ht_mi)
  dev.off()
  
  log_message("Immune infiltration analysis complete.")
} else {
  log_message("xCell analysis failed. Using alternative approach with ssGSEA immune signatures.", "WARNING")
  
  # Alternative: ssGSEA with immune signatures
  cibersort_genes <- msigdbr(species = "Homo sapiens", category = "C7")
  # Focus on immune-related gene sets
  immune_sigs <- cibersort_genes[grep("IMMUNE|TCELL|BCELL|MACROPHAGE|DENDRITIC|NK_CELL|NEUTROPHIL", 
                                       cibersort_genes$gs_name, ignore.case = TRUE), ]
  
  if (nrow(immune_sigs) > 0) {
    immune_gs_list <- split(immune_sigs$gene_symbol, immune_sigs$gs_name)
    immune_gs_list <- lapply(immune_gs_list, function(x) intersect(x, rownames(combat_matrix)))
    immune_gs_list <- immune_gs_list[sapply(immune_gs_list, length) >= 10]
    
    if (length(immune_gs_list) > 0) {
      gsva_par <- gsvaParam(combat_matrix, immune_gs_list, minSize = 10, maxSize = 500)
      ssgsea_immune <- gsva(gsva_par, method = "ssgsea")
      
      saveRDS(ssgsea_immune, file.path(OUTPUT$immune, "ssgsea_immune_scores.rds"))
      
      ht_imm <- Heatmap(ssgsea_immune,
                       name = "Score",
                       col = colorRamp2(c(-0.3, 0, 0.3), c("#2166AC", "white", "#B2182B")),
                       top_annotation = HeatmapAnnotation(
                         Group = discovery_metadata$group,
                         col = list(Group = c(Normal = "#0072B2", Tumor = "#D55E00"))
                       ),
                       show_column_names = FALSE,
                       show_row_names = FALSE,
                       column_title = "ssGSEA Immune Signatures")
      
      open_tiff(file.path(OUTPUT$immune, "immune_heatmap_ssgsea.tiff"), width = 10, height = 6)
      draw(ht_imm)
      dev.off()
      
      log_message("ssGSEA immune analysis complete.")
    }
  }
}

log_message("PPI and immune infiltration analysis complete.")
