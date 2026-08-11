source("R/utils.R")
load_libraries()
set.seed(123)

BASE_DIR <- "F:/independent projects/cancer/breast cancer"
OUT_DIR <- file.path(BASE_DIR, "Breast cancer 6 Ferroptosis")

log_msg("=== Step 12: PPI Network Analysis ===")

de_frg <- readRDS(file.path(OUT_DIR, "results/ferroptosis_genes/DE_FRG_results.rds"))
top_hubs <- readRDS(file.path(OUT_DIR, "results/WGCNA/top_hub_genes.rds"))
combat_expr <- readRDS(file.path(OUT_DIR, "data_processed/discovery/discovery_combat.rds"))

# Build PPI from STRINGdb if available, otherwise use correlation-based network
library(igraph)

# Correlation-based PPI
hub_genes_final <- unique(c(rownames(de_frg), head(top_hubs$Gene, 30)))
hub_genes_final <- intersect(hub_genes_final, rownames(combat_expr))
log_msg(sprintf("Building PPI for %d genes", length(hub_genes_final)))

if (length(hub_genes_final) >= 3) {
  expr_hub <- combat_expr[hub_genes_final, ]
  cor_mat <- cor(t(expr_hub), method = "spearman")
  cor_mat[is.na(cor_mat)] <- 0

  # Threshold for edges
  cor_threshold <- quantile(abs(cor_mat), probs = 0.75, na.rm = TRUE)
  adj_mat <- abs(cor_mat) >= cor_threshold & upper.tri(cor_mat)

  # Build graph
  graph <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
  if (vcount(graph) < 2) {
    log_msg("Network too small, lowering threshold")
    cor_threshold <- quantile(abs(cor_mat), probs = 0.5, na.rm = TRUE)
    adj_mat <- abs(cor_mat) >= cor_threshold & upper.tri(cor_mat)
    graph <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", diag = FALSE)
  }

  log_msg(sprintf("PPI network: %d nodes, %d edges", vcount(graph), ecount(graph)))

  if (vcount(graph) > 1 && ecount(graph) > 0) {
    # Centrality measures
    node_metrics <- data.frame(
      Node = V(graph)$name,
      Degree = degree(graph),
      Betweenness = betweenness(graph, directed = FALSE),
      Closeness = closeness(graph, mode = "all"),
      Eigenvector = eigen_centrality(graph, directed = FALSE)$vector,
      stringsAsFactors = FALSE
    )
    write.csv(node_metrics, file.path(OUT_DIR, "results/PPI_network/PPI_node_metrics.csv"), row.names = FALSE)

    # Edge list
    edge_list <- as.data.frame(as_edgelist(graph))
    colnames(edge_list) <- c("Node1", "Node2")
    edge_list$Weight <- cor_mat[as.matrix(edge_list)]
    write.csv(edge_list, file.path(OUT_DIR, "results/PPI_network/PPI_edge_list.csv"), row.names = FALSE)

    # Community detection
    community <- cluster_walktrap(graph)
    V(graph)$community <- membership(community)

    # PPI network plot
    open_tiff(file.path(OUT_DIR, "results/figures/Figure10A_PPI_network.tiff"), width = 10, height = 8)
    set.seed(123)
    layout <- layout_with_fr(graph)
    plot(community, graph, layout = layout,
         vertex.label.cex = 0.6, vertex.size = log(node_metrics$Degree + 1) * 3 + 2,
         vertex.label.color = "black", edge.width = 0.5,
         main = "PPI Network of Candidate Hub Genes")
    dev.off()

    # Centrality plot
    p_centrality <- ggplot(node_metrics, aes(Degree, Betweenness, size = Eigenvector, color = Closeness)) +
      geom_point(alpha = 0.7) +
      geom_text_repel(aes(label = Node), size = 2, max.overlaps = 15) +
      scale_color_gradient(low = "blue", high = "red") +
      labs(title = "Node Centrality Metrics",
           x = "Degree", y = "Betweenness") +
      theme_minimal()
    save_tiff_600dpi(p_centrality, file.path(OUT_DIR, "results/figures/Figure10B_centrality_plot.tiff"),
                     width = 8, height = 6)

    # Final hub-gene network (top 15 by degree)
    top_degree <- head(node_metrics[order(node_metrics$Degree, decreasing = TRUE), ], 15)
    top_graph <- induced_subgraph(graph, top_degree$Node)

    if (vcount(top_graph) > 1) {
      open_tiff(file.path(OUT_DIR, "results/figures/Figure10E_final_hub_network.tiff"), width = 8, height = 7)
      set.seed(123)
      top_layout <- layout_with_fr(top_graph)
      V(top_graph)$color <- ifelse(V(top_graph)$name %in% rownames(de_frg), "red", "steelblue")
      plot(top_graph, layout = top_layout,
           vertex.label.cex = 0.7, vertex.size = 12,
           vertex.label.color = "black", vertex.color = V(top_graph)$color,
           edge.width = 0.5, main = "Top Hub Gene Network")
      legend("topleft", legend = c("DE-FRG", "Other Hub"), fill = c("red", "steelblue"), cex = 0.8)
      dev.off()
    }

    # Save graph
    saveRDS(graph, file.path(OUT_DIR, "results/PPI_network/PPI_graph.rds"))
    saveRDS(node_metrics, file.path(OUT_DIR, "results/PPI_network/node_metrics.rds"))
    write.csv(data.frame(community = membership(community)), file.path(OUT_DIR, "results/PPI_network/community_assignments.csv"))
  }
}

log_msg("Step 12 complete: PPI network analysis finished")
