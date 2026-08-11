# ============================================================
# Step 8: Crosstalk and Network Analysis
# ============================================================
set.seed(123)
library(igraph)
library(ggplot2)
library(corrplot)
library(circlize)
library(clusterProfiler)
library(org.Hs.eg.db)

base_dir <- normalizePath(file.path(".", "Breast Cancer 5 Necroptosis"))
output_dir <- base_dir

sink(file.path(output_dir, "logs", "08_networks.log"), split=TRUE)
cat("CROSSTALK AND NETWORK ANALYSIS\n")
cat("==============================\n\n")

save_tiff_600dpi <- function(plot, filename, width = 7, height = 6) {
  ggsave(filename=filename, plot=plot, width=width, height=height,
         units="in", dpi=600, compression="lzw", bg="white", device="tiff")
}

combat_expr <- readRDS(file.path(output_dir, "data_processed", "discovery_combat_corrected.rds"))
merged_meta <- readRDS(file.path(output_dir, "data_processed", "discovery_metadata.rds"))
de_nrgs <- readLines(file.path(output_dir, "results", "necroptosis_genes", "DE_NRGs.txt"))
deg_sig <- read.csv(file.path(output_dir, "results", "tables", "DEG_significant.csv"), row.names=1)
deg_all <- read.csv(file.path(output_dir, "results", "tables", "DEG_all.csv"), row.names=1)

# Find breast cancer-associated genes from DEGs (top 100 most significant)
bc_genes <- rownames(deg_sig)[order(deg_sig$adj.P.Val)][1:min(100, nrow(deg_sig))]

# Spearman correlation between DE-NRGs and breast cancer genes
if(length(de_nrgs) > 1 && length(bc_genes) > 1) {
  common <- intersect(rownames(combat_expr), c(de_nrgs, bc_genes))
  de_nrg_common <- intersect(de_nrgs, common)
  bc_common <- intersect(bc_genes, common)
  
  sub_expr <- combat_expr[common, ]
  
  cat("DE-NRGs in network:", length(de_nrg_common), "\n")
  cat("BC genes in network:", length(bc_common), "\n")
  
  # Correlation matrix
  cor_net <- cor(t(sub_expr), method="spearman")
  
  # Save full correlation
  write.csv(cor_net, file.path(output_dir, "results", "networks", "gene_correlation_matrix.csv"))
  
  # Significant correlations (|r| > 0.6, p < 0.01)
  cor_pval <- cor.test.helper <- function(x, y) {
    t <- cor.test(x, y, method="spearman")
    c(t$estimate, t$p.value)
  }
  
  n_genes <- nrow(sub_expr)
  cor_pval_mat <- matrix(1, n_genes, n_genes)
  rownames(cor_pval_mat) <- rownames(sub_expr)
  colnames(cor_pval_mat) <- rownames(sub_expr)
  
  # Fast p-value using approximated t-distribution for Spearman
  n_samp <- ncol(sub_expr)
  cor_t <- cor_net * sqrt((n_samp - 2) / (1 - cor_net^2))
  cor_pval <- 2 * pt(-abs(cor_t), df=n_samp - 2)
  cor_pval_adj <- matrix(p.adjust(as.vector(cor_pval), "BH"), n_genes, n_genes)
  rownames(cor_pval_adj) <- rownames(sub_expr)
  colnames(cor_pval_adj) <- rownames(sub_expr)
  
  # Significant correlation pairs
  sig_pairs <- which(abs(cor_net) > 0.6 & cor_pval_adj < 0.01 & upper.tri(cor_net), arr.ind=TRUE)
  cat("Significant correlation pairs:", nrow(sig_pairs), "\n")
  
  if(nrow(sig_pairs) > 0) {
    edge_list <- data.frame(
      from = rownames(cor_net)[sig_pairs[,1]],
      to = colnames(cor_net)[sig_pairs[,2]],
      weight = cor_net[sig_pairs],
      abs_weight = abs(cor_net[sig_pairs])
    )
    write.csv(edge_list, file.path(output_dir, "results", "networks", "significant_correlation_edges.csv"), row.names=FALSE)
    
    # Network graph
    net_graph <- graph_from_data_frame(edge_list, directed=FALSE)
    cat("Network vertices:", vcount(net_graph), "Edges:", ecount(net_graph), "\n")
    
    # Centrality measures (use igraph:: explicitly to avoid circlize masking)
    cent_degree <- as.numeric(igraph::degree(net_graph))
    names(cent_degree) <- names(igraph::degree(net_graph))
    cent_betweenness <- tryCatch(as.numeric(igraph::betweenness(net_graph, directed=FALSE, weights=abs(E(net_graph)$weight))), 
                                  error=function(e) rep(NA, vcount(net_graph)))
    cent_closeness <- tryCatch(as.numeric(igraph::closeness(net_graph, weights=abs(E(net_graph)$weight))),
                                error=function(e) rep(NA, vcount(net_graph)))
    cent_eigen <- tryCatch(as.numeric(igraph::eigen_centrality(net_graph)$vector),
                            error=function(e) rep(NA, vcount(net_graph)))
    names(cent_betweenness) <- names(cent_closeness) <- names(cent_eigen) <- names(cent_degree)
    
    centrality_df <- data.frame(
      Gene = names(cent_degree),
      Degree = cent_degree,
      Betweenness = cent_betweenness,
      Closeness = cent_closeness,
      Eigenvector = cent_eigen
    )
    write.csv(centrality_df, file.path(output_dir, "results", "networks", "centrality_measures.csv"), row.names=FALSE)
    
    # Centrality ranking plot
    top_degree <- head(centrality_df[order(centrality_df$Degree, decreasing=TRUE), ], 20)
    p_cent <- ggplot(top_degree, aes(x=reorder(Gene, Degree), y=Degree)) +
      geom_bar(stat="identity", fill="steelblue") + theme_bw() +
      coord_flip() + labs(title="Top 20 Genes by Degree Centrality", x="", y="Degree")
    save_tiff_600dpi(p_cent, file.path(output_dir, "results", "networks", "centrality_ranking.tiff"), width=8, height=6)
    
    # Node and edge tables for Cytoscape
    node_table <- centrality_df
    write.csv(node_table, file.path(output_dir, "results", "networks", "cytoscape_nodes.csv"), row.names=FALSE)
    cytoscape_edges <- edge_list[, c("from", "to", "weight")]
    colnames(cytoscape_edges)[3] <- "Spearman_r"
    write.csv(cytoscape_edges, file.path(output_dir, "results", "networks", "cytoscape_edges.csv"), row.names=FALSE)
    cat("Cytoscape-compatible node/edge tables saved.\n")
  }
}

# Gene-pathway network using enrichment results
if(length(de_nrgs) >= 5) {
  tryCatch({
    gene_entrez <- bitr(de_nrgs, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
    ego <- enrichGO(gene=gene_entrez$ENTREZID, OrgDb=org.Hs.eg.db,
                    ont="BP", pvalueCutoff=0.05, readable=TRUE)
    if(!is.null(ego) && nrow(ego) > 0) {
      tiff(file.path(output_dir, "results", "networks", "gene_pathway_cnet.tiff"),
           width=12, height=10, units="in", res=600, compression="lzw", bg="white")
      print(cnetplot(ego, showCategory=5, node_label="all"))
      dev.off()
    }
  }, error=function(e) cat("  Gene-pathway network error:", e$message, "\n"))
}

# Chord/Circos plot of top correlations (simplified)
if(exists("edge_list") && nrow(edge_list) > 0) {
  top_edges <- edge_list[order(edge_list$abs_weight, decreasing=TRUE), ][1:min(30, nrow(edge_list)), ]
  
  tiff(file.path(output_dir, "results", "networks", "circos_plot.tiff"),
       width=10, height=10, units="in", res=600, compression="lzw", bg="white")
  tryCatch({
    circos.clear()
    chordDiagram(top_edges[, c("from", "to", "weight")], 
                 col=ifelse(top_edges$weight > 0, rgb(1,0,0,0.5), rgb(0,0,1,0.5)),
                 transparency=0.3,
                 annotationTrack=c("grid", "axis"))
    title("Gene-Gene Crosstalk (Circos)")
  }, error=function(e) cat("  Circos error:", e$message, "\n"))
  dev.off()
  circos.clear()
}

cat("\n\nNETWORK ANALYSIS COMPLETE\n")
sink()
