# 11_PPI_and_network_analysis.R
source("config/analysis_config.R")
source("R/load_functions.R")
ensure_paths()
suppressMessages({ library(ggplot2); library(igraph); library(ggraph); library(dplyr) })

crg_sum <- readRDS(file.path(paths$res_cross, "DE_CRG_summary.rds"))
deg <- readRDS(file.path(paths$res_de, "discovery_DEG_full.rds"))
cons_file <- file.path(paths$res_ml, "consensus", "consensus_genes.rds")
genes_ppi <- unique(c(crg_sum$de_crgs, if (file.exists(cons_file)) readRDS(cons_file)$final_genes else character()))

## PPI via STRINGdb
ppi <- get_ppi_network(genes_ppi)
if (!ppi$available) {
  cat("STRINGdb unavailable; writing report and using correlation network only.\n")
  writeLines("STRINGdb API unavailable for this environment; PPI analysis substituted with correlation-network centrality.",
             file.path(paths$res_ppi, "STRINGdb_unavailable_report.txt"))
  ## fallback: use DE-CRG correlation edges for centrality
  edges <- read.csv(file.path(paths$res_cross, "DECRG_DEG_significant_edges.csv"), stringsAsFactors = FALSE)
  if (nrow(edges) > 0) {
    g <- graph_from_data_frame(edges[, c("from","to")], directed = FALSE)
    cent <- data.frame(gene = V(g)$name, degree = degree(g), stringsAsFactors = FALSE)
    write.csv(cent, file.path(paths$res_ppi, "PPI_centrality.csv"), row.names = FALSE)
  }
} else {
  ## remap edges and nodes from STRING/Ensembl ids to gene symbols
  sym_by_string <- setNames(ppi$nodes$gene, ppi$nodes$STRING_id)
  edges_sym <- ppi$edges
  edges_sym$from <- sym_by_string[edges_sym$from]
  edges_sym$to <- sym_by_string[edges_sym$to]
  edges_sym <- edges_sym[!is.na(edges_sym$from) & !is.na(edges_sym$to), ]
  nodes_sym <- ppi$nodes; nodes_sym$name <- nodes_sym$gene
  cent <- ppi_centrality(edges_sym, nodes_sym)
  cent$gene <- as.character(cent$gene)
  write.csv(cent, file.path(paths$res_ppi, "PPI_centrality.csv"), row.names = FALSE)
  ## PPI network plot (gene symbols)
  g <- graph_from_data_frame(edges_sym[, c("from","to")], directed = FALSE, vertices = nodes_sym)
  V(g)$deg <- degree(g)
  p_ppi <- ggraph(g, layout = "fr") + geom_edge_link(alpha = 0.3) +
    geom_node_point(aes(size = deg), colour = "#E15759") +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    labs(title = "PPI network (DE-CRGs + consensus)") + pub_theme()
  save_publication_plot(p_ppi, file.path(paths$res_ppi, "PPI_network"), width = 8, height = 7)
  ## node-degree plot (gene symbols)
  if (nrow(cent) > 0) {
    p_deg <- ggplot(cent %>% dplyr::top_n(20, degree), aes(x = reorder(gene, degree), y = degree)) +
      geom_bar(stat = "identity", fill = "#59A14F") + coord_flip() + labs(title = "PPI node degree") + pub_theme()
    save_publication_plot(p_deg, file.path(paths$res_ppi, "PPI_node_degree"))
  }
  ## export cytoscape tables (gene symbols)
  write.csv(nodes_sym[, c("gene","STRING_id")], file.path(paths$res_ppi, "PPI_nodes_cytoscape.csv"), row.names = FALSE)
  write.csv(edges_sym[, c("from","to","combined_score","from_name","to_name")],
            file.path(paths$res_ppi, "PPI_edges_cytoscape.csv"), row.names = FALSE)
  ppi_sym <- list(available = TRUE, nodes = nodes_sym, edges = edges_sym, string_db = ppi$string_db)
  saveRDS(ppi_sym, file.path(paths$res_ppi, "PPI_object.rds"))
}
cat("PPI & network analysis complete.\n")
quit(save = "no", status = 0, runLast = FALSE)
