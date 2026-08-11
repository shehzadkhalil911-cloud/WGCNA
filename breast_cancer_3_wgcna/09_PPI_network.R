###############################################################################
# SCRIPT 9: PPI Network Analysis of Hub Genes
# Multi-Cohort WGCNA: Breast Cancer Hub Gene Discovery
###############################################################################
library(STRINGdb)
library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)
library(ggplot2)

options(stringsAsFactors = FALSE)

OUT_DIR <- "F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)"
FIG_DIR <- file.path(OUT_DIR, "figures/enrichment")
TAB_DIR <- file.path(OUT_DIR, "tables")
RDA_DIR <- file.path(OUT_DIR, "RData")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

load(file.path(RDA_DIR, "05_hub_genes.RData"))

cat("========== PPI Network Analysis ==========\n")
cat(sprintf("Hub genes: %d\n", length(final_hubs)))

# Initialize STRINGdb
string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 400, input_directory = "")

# Map genes to STRING identifiers
string_ids <- string_db$map(data.frame(Gene = final_hubs), "Gene", removeUnmappedRows = TRUE)
cat(sprintf("Mapped to STRING: %d/%d\n", nrow(string_ids), length(final_hubs)))

if (nrow(string_ids) >= 3) {
  # Get PPI network
  ppi <- string_db$get_subnetwork(string_ids$STRING_id)
  
  if (!is.null(ppi) && vcount(ppi) > 0) {
    cat(sprintf("PPI network: %d nodes, %d edges\n", vcount(ppi), ecount(ppi)))
    
    # Add gene symbols as node attributes
    gene_map <- setNames(string_ids$Gene, string_ids$STRING_id)
    V(ppi)$name <- gene_map[V(ppi)$name]
    V(ppi)$name[is.na(V(ppi)$name)] <- "Unknown"
    
    # Add hub score info
    hub_score_map <- setNames(hub_summary$meta_FDR, hub_summary$Gene)
    V(ppi)$meta_FDR <- hub_score_map[V(ppi)$name]
    # Add direction info
    dir_map <- setNames(rowMeans(hub_summary[, c("logFC_42568", "logFC_15852", "logFC_29431")]),
                        hub_summary$Gene)
    V(ppi)$mean_logFC <- dir_map[V(ppi)$name]
    
    # Compute network metrics
    V(ppi)$degree <- degree(ppi, mode = "all")
    V(ppi)$betweenness <- betweenness(ppi, directed = FALSE)
    V(ppi)$closeness <- closeness(ppi, mode = "all")
    
    # Save node table
    node_df <- data.frame(
      Gene = V(ppi)$name,
      Degree = V(ppi)$degree,
      Betweenness = V(ppi)$betweenness,
      Closeness = V(ppi)$closeness,
      Mean_logFC = V(ppi)$mean_logFC,
      Meta_FDR = V(ppi)$meta_FDR,
      stringsAsFactors = FALSE
    )
    node_df <- node_df[order(-node_df$Degree), ]
    write.csv(node_df, file.path(TAB_DIR, "PPI_node_metrics.csv"), row.names = FALSE)
    
    # Save edge table
    edge_df <- as.data.frame(as_edgelist(ppi))
    colnames(edge_df) <- c("Node1", "Node2")
    edge_df$CombinedScore <- E(ppi)$combined_score
    write.csv(edge_df, file.path(TAB_DIR, "PPI_edge_list.csv"), row.names = FALSE)
    
    # Convert to tidygraph for ggraph
    ppi_tidy <- as_tbl_graph(ppi)
    
    # PPI Network plot - main
    set.seed(123)
    p_main <- ggraph(ppi_tidy, layout = "fr") +
      geom_edge_link(aes(alpha = combined_score), show.legend = FALSE, width = 0.5) +
      geom_node_point(aes(size = degree, color = mean_logFC)) +
      geom_node_text(aes(label = name), size = 2.5, repel = TRUE, max.overlaps = 20) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                            name = "Mean logFC") +
      scale_size_continuous(name = "Degree") +
      labs(title = "Hub Gene PPI Network (STRING)") +
      theme_graph() + theme(plot.title = element_text(hjust = 0.5, size = 16))
    
    ggsave(file.path(FIG_DIR, "PPI_network_main.tiff"), p_main,
           width = 14, height = 12, dpi = 300, compression = "lzw")
    
    # Hub-highlighted plot (top 10 by degree)
    V(ppi_tidy)$is_top <- V(ppi_tidy)$name %in% head(node_df$Gene, 10)
    p_hubs <- ggraph(ppi_tidy, layout = "fr") +
      geom_edge_link(aes(alpha = combined_score), show.legend = FALSE, width = 0.3) +
      geom_node_point(aes(size = degree, color = is_top), alpha = 0.8) +
      geom_node_text(aes(label = ifelse(is_top, name, ""), color = is_top),
                     size = 3, repel = TRUE, max.overlaps = 15, fontface = "bold") +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey60")) +
      scale_size_continuous(name = "Degree") +
      labs(title = "PPI Network - Top 10 Hub Genes Highlighted") +
      theme_graph() + theme(plot.title = element_text(hjust = 0.5, size = 16))
    
    ggsave(file.path(FIG_DIR, "PPI_network_top10.tiff"), p_hubs,
           width = 12, height = 10, dpi = 300, compression = "lzw")
    
    # Degree distribution plot
    p_deg <- ggplot(node_df, aes(x = Degree)) +
      geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
      labs(title = "Degree Distribution", x = "Degree", y = "Count") +
      theme_bw()
    ggsave(file.path(FIG_DIR, "PPI_degree_distribution.tiff"), p_deg,
           width = 8, height = 6, dpi = 300, compression = "lzw")
    
    cat(sprintf("\nTop 10 hub genes by degree:\n"))
    print(head(node_df, 10))
    
    # Save PPI object
    save(ppi, ppi_tidy, node_df, edge_df, file = file.path(RDA_DIR, "09_PPI_network.RData"))
    
  } else {
    cat("PPI network is empty.\n")
  }
} else {
  cat("Too few genes mapped to STRING.\n")
}

cat("\nScript 9 complete.\n")
