# =============================================================================
# Script 08: Pathway Crosstalk Analysis
# Breast Cancer 7 Cell Death Crosstalk
# =============================================================================

source(file.path("..", "config", "analysis_config.R"))

library(tidyverse)
library(clusterProfiler)
library(ComplexHeatmap)
library(circlize)
library(igraph)
library(ggraph)
library(org.Hs.eg.db)

log_message("Starting pathway crosstalk analysis.")

# --- 1. LOAD DATA ---
shared_info <- readRDS(file.path(OUTPUT$shared, "shared_gene_classification.rds"))
deg_sig <- read.csv(file.path(OUTPUT$de, "discovery_DEG_significant.csv"), row.names = 1)

# --- 2. CROSSTALK PATHWAYS ---
crosstalk_pathways <- list(
  `Iron Metabolism` = c("TFRC", "SLC40A1", "HFE", "FTH1", "FTL", "HAMP", "HMOX1", "BLVRA", "STEAP3", "CP", "TF"),
  `Lipid Peroxidation` = c("GPX4", "ACSL4", "LPCAT3", "ALOX12", "ALOX15", "PTGS2", "CYP1B1", "CYP2E1", "NOX1", "NOX4"),
  `Glutathione Metabolism` = c("GCLC", "GCLM", "GSR", "GPX4", "GPX1", "GPX2", "GPX3", "SLC7A11", "SLC3A2", "GSTP1", "MGST1"),
  `Copper Metabolism` = c("ATP7A", "ATP7B", "SLC31A1", "COX17", "SOD1", "SOD3", "LOX", "CCS", "ATOX1", "COMMD1"),
  `TCA Cycle` = c("CS", "ACO2", "IDH3A", "IDH3B", "IDH3G", "OGDH", "DLST", "SUCLA2", "SUCLG1", "SDHA", "SDHB", "FH", "MDH2"),
  `Mitochondrial Metabolism` = c("COX4I1", "COX5A", "COX6A1", "COX7A2", "COX8A", "COX10", "COX15", "SCO1", "SCO2", "NDUFS1", "NDUFS8", "UQCRC1", "UQCRC2"),
  `Oxidative Stress` = c("NFE2L2", "KEAP1", "HMOX1", "NQO1", "SOD1", "SOD2", "CAT", "PRDX1", "PRDX6", "TXN", "TXNRD1", "GSR"),
  `TNF Signalling` = c("TNF", "TNFRSF1A", "TNFRSF1B", "TRADD", "FADD", "RIPK1", "TRAF2", "TRAF5", "MAP3K7", "IKBKB", "CHUK", "NFKB1"),
  `RIPK1-RIPK3-MLKL` = c("RIPK1", "RIPK3", "MLKL", "RIPK2", "CYLD", "SPATA2", "BIRC2", "BIRC3", "TRAF2", "TRIM11"),
  `p53 Signalling` = c("TP53", "CDKN1A", "MDM2", "MDM4", "BAX", "BBC3", "PMAIP1", "SESN1", "SESN2", "TIGAR", "GADD45A"),
  `NF-kB Signalling` = c("NFKB1", "RELA", "IKBKB", "IKBKG", "NFKBIA", "NFKBIE", "RELB", "REL", "MAP3K7", "BIRC3"),
  `Hypoxia` = c("HIF1A", "EPAS1", "VEGFA", "SLC2A1", "LDHA", "PDK1", "BNIP3", "CA9", "EGLN1", "EGLN3", "ADM"),
  `Autophagy` = c("BECN1", "ATG5", "ATG7", "ATG12", "SQSTM1", "ULK1", "ULK2", "LAMP1", "LAMP2", "CTSD", "MAP1LC3A", "GABARAP"),
  `Apoptosis` = c("BAX", "BAK1", "BCL2", "BCL2L1", "MCL1", "BID", "BAD", "BIM", "CASP3", "CASP8", "CASP9", "CYCS", "APAF1"),
  `Inflammation` = c("IL6", "IL1B", "TNF", "CXCL8", "CCL2", "CCL5", "CXCL10", "IL10", "TGFB1", "IL18", "NLRP3"),
  `Immune Signalling` = c("CD8A", "CD4", "CD274", "PDCD1", "CTLA4", "IFNG", "GZMB", "PRF1", "CD86", "CD80", "ICOS", "FOXP3")
)

# --- 3. FILTER TO GENES IN OUR DATA ---
all_de_genes <- rownames(deg_sig)

crosstalk_filtered <- list()
for (pw in names(crosstalk_pathways)) {
  genes_in_de <- intersect(crosstalk_pathways[[pw]], all_de_genes)
  if (length(genes_in_de) > 0) crosstalk_filtered[[pw]] <- genes_in_de
}

log_message(sprintf("Crosstalk pathways with DE genes: %d", length(crosstalk_filtered)))
for (pw in names(crosstalk_filtered)) {
  log_message(sprintf("  %s: %d DE genes", pw, length(crosstalk_filtered[[pw]])))
}

# --- 4. PATHWAY OVERLAP HEATMAP ---
pathway_overlap <- matrix(0, nrow = length(crosstalk_filtered), ncol = length(crosstalk_filtered))
rownames(pathway_overlap) <- names(crosstalk_filtered)
colnames(pathway_overlap) <- names(crosstalk_filtered)

for (i in seq_along(crosstalk_filtered)) {
  for (j in seq_along(crosstalk_filtered)) {
    pathway_overlap[i, j] <- length(intersect(crosstalk_filtered[[i]], crosstalk_filtered[[j]]))
  }
}

ht_pw <- Heatmap(pathway_overlap,
                 name = "Gene Overlap",
                 col = colorRamp2(c(0, max(pathway_overlap)), c("white", "#D55E00")),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   if (pathway_overlap[i, j] > 0)
                     grid.text(pathway_overlap[i, j], x, y, gp = gpar(fontsize = 6))
                 },
                 rect_gp = gpar(col = "grey60", lwd = 0.5),
                 row_names_gp = gpar(fontsize = 8),
                 column_names_gp = gpar(fontsize = 8),
                 column_title = "Pathway-Pathway Overlap",
                 row_title = "Pathway",
                 cluster_rows = TRUE,
                 cluster_columns = TRUE)

open_tiff(file.path(OUTPUT$crosstalk, "pathway_overlap_heatmap.tiff"), width = 8, height = 7)
draw(ht_pw)
dev.off()

# --- 5. PATHWAY-PROGRAMME OVERLAP ---
# Which pathways contain DE-FRGs, DE-NRGs, DE-CRGs
program_pathway <- data.frame()

for (pw in names(crosstalk_filtered)) {
  for (prog_name in c("DE_FRGs", "DE_NRGs", "DE_CRGs")) {
    prog_genes <- shared_info$Gene[shared_info[[prog_name]]]
    overlap <- intersect(crosstalk_filtered[[pw]], prog_genes)
    if (length(overlap) > 0) {
      program_pathway <- rbind(program_pathway, data.frame(
        Pathway = pw,
        Program = gsub("DE_", "", gsub("Gs$", "", prog_name)),
        Overlap_Genes = length(overlap),
        Genes = paste(overlap, collapse = ","),
        stringsAsFactors = FALSE
      ))
    }
  }
}

write.csv(program_pathway, file.path(OUTPUT$crosstalk, "program_pathway_overlap.csv"), row.names = FALSE)

# --- 6. THREE-PROGRAMME PATHWAY NETWORK ---
if (nrow(program_pathway) > 0) {
  nodes <- data.frame(
    name = unique(c(program_pathway$Pathway, program_pathway$Program)),
    type = c(rep("Pathway", length(unique(program_pathway$Pathway))),
             rep("Program", length(unique(program_pathway$Program)))),
    stringsAsFactors = FALSE
  )
  
  edges <- data.frame(
    from = program_pathway$Pathway,
    to = program_pathway$Program,
    weight = program_pathway$Overlap_Genes,
    stringsAsFactors = FALSE
  )
  
  g_pw <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)
  
  p_net <- ggraph(g_pw, layout = "kk") +
    geom_edge_link(aes(width = weight), alpha = 0.6, color = "grey60") +
    geom_node_point(aes(color = type), size = 8) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_edge_width_continuous(range = c(0.5, 3)) +
    scale_color_manual(values = c("Program" = "#CC79A7", "Pathway" = "#0072B2")) +
    labs(title = "Program-Pathway Network") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  save_tiff_600dpi(p_net, file.path(OUTPUT$crosstalk, "program_pathway_network.tiff"), width = 10, height = 8)
}

# --- 7. GENE-PATHWAY CHORD DIAGRAM ---
# For top shared genes
gene_pathway_edges <- data.frame()
for (pw in names(crosstalk_filtered)) {
  for (gene in crosstalk_filtered[[pw]]) {
    gene_pathway_edges <- rbind(gene_pathway_edges, data.frame(
      Gene = gene, Pathway = pw, stringsAsFactors = FALSE
    ))
  }
}

write.csv(gene_pathway_edges, file.path(OUTPUT$crosstalk, "gene_pathway_edges.csv"), row.names = FALSE)

# --- 8. MODULE-PATHWAY HEATMAP ---
module_df <- read.csv(file.path(OUTPUT$wgcna, "module_assignments.csv"))
module_pathway <- matrix(0, nrow = length(unique(module_df$Module)), ncol = length(crosstalk_filtered))
rownames(module_pathway) <- paste0("M", sort(unique(module_df$Module)))
colnames(module_pathway) <- names(crosstalk_filtered)

for (i in seq_len(nrow(module_df))) {
  gene <- module_df$Gene[i]
  mod <- paste0("M", module_df$Module[i])
  for (pw in names(crosstalk_filtered)) {
    if (gene %in% crosstalk_filtered[[pw]]) {
      module_pathway[mod, pw] <- module_pathway[mod, pw] + 1
    }
  }
}

# Only show modules with at least one hit
module_pathway <- module_pathway[rowSums(module_pathway) > 0, , drop = FALSE]
module_pathway <- module_pathway[, colSums(module_pathway) > 0, drop = FALSE]

if (nrow(module_pathway) > 0 && ncol(module_pathway) > 0) {
  ht_mp <- Heatmap(module_pathway,
                   name = "Genes",
                   col = colorRamp2(c(0, max(module_pathway)), c("white", "#009E73")),
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     if (module_pathway[i, j] > 0)
                       grid.text(module_pathway[i, j], x, y, gp = gpar(fontsize = 7))
                   },
                   rect_gp = gpar(col = "grey60", lwd = 0.5),
                   row_names_gp = gpar(fontsize = 8),
                   column_names_gp = gpar(fontsize = 8, angle = 45, hjust = 1),
                   column_title = "Module-Pathway Overlap",
                   row_title = "Module")
  
  open_tiff(file.path(OUTPUT$crosstalk, "module_pathway_heatmap.tiff"), width = 9, height = 7)
  draw(ht_mp)
  dev.off()
}

# --- 9. PATHWAY SIMILARITY NETWORK ---
# Jaccard similarity between pathways
pathway_jaccard <- matrix(0, nrow = length(crosstalk_filtered), ncol = length(crosstalk_filtered))
rownames(pathway_jaccard) <- names(crosstalk_filtered)
colnames(pathway_jaccard) <- names(crosstalk_filtered)

for (i in seq_along(crosstalk_filtered)) {
  for (j in seq_along(crosstalk_filtered)) {
    a <- length(crosstalk_filtered[[i]])
    b <- length(crosstalk_filtered[[j]])
    int <- length(intersect(crosstalk_filtered[[i]], crosstalk_filtered[[j]]))
    pathway_jaccard[i, j] <- ifelse(a + b - int > 0, int / (a + b - int), 0)
  }
}

# Network of pathways with Jaccard > 0.1
pw_edges <- which(pathway_jaccard > 0.1 & upper.tri(pathway_jaccard), arr.ind = TRUE)
if (nrow(pw_edges) > 0) {
  pw_net_df <- data.frame(
    from = rownames(pathway_jaccard)[pw_edges[, 1]],
    to = colnames(pathway_jaccard)[pw_edges[, 2]],
    similarity = pathway_jaccard[pw_edges],
    stringsAsFactors = FALSE
  )
  
  g_pw_sim <- graph_from_data_frame(pw_net_df, directed = FALSE)
  
  p_pw_net <- ggraph(g_pw_sim, layout = "fr") +
    geom_edge_link(aes(width = similarity, alpha = similarity), color = "grey60") +
    geom_node_point(size = 5, color = "#0072B2") +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_edge_width_continuous(range = c(0.5, 2)) +
    labs(title = "Pathway Similarity Network (Jaccard > 0.1)") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  save_tiff_600dpi(p_pw_net, file.path(OUTPUT$crosstalk, "pathway_similarity_network.tiff"), width = 10, height = 8)
}

log_message("Pathway crosstalk analysis complete.")
