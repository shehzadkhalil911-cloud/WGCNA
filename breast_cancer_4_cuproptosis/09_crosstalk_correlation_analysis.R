# 09_crosstalk_correlation_analysis.R
source("config/analysis_config.R")
source("R/load_functions.R")
ensure_paths()
suppressMessages({ library(ggplot2); library(ggpubr); library(pheatmap); library(Hmisc); library(igraph); library(ggraph); library(corrplot) })

combat <- readRDS(file.path(paths$discovery, "discovery_combat.rds"))
deg <- readRDS(file.path(paths$res_de, "discovery_DEG_full.rds"))
crg_sum <- readRDS(file.path(paths$res_cross, "DE_CRG_summary.rds"))
crg_det <- crg_sum$crg_detected
de_crgs <- crg_sum$de_crgs
meta <- readRDS(file.path(paths$discovery, "discovery_cohort_raw.rds"))$metadata

## top non-CRG DEGs
top_deg <- deg$symbol[deg$sig & !(deg$symbol %in% crg_det)]
top_deg <- head(top_deg[order(deg$adj.P.Val[deg$sig & !(deg$symbol %in% crg_det)])], 100)

## 1) CRG correlation (all detected)
crg_mat <- combat[crg_det, ]
c_crg <- rcorr(t(crg_mat), type = COR_METHOD)
crg_corr_edges <- correlation_network(crg_mat, method = COR_METHOD, p_thr = COR_P_THRESHOLD)
write.csv(crg_corr_edges, file.path(paths$res_cross, "CRG_correlation_edges.csv"), row.names = FALSE)

## CRG correlation heatmap
pheatmap::pheatmap(c_crg$r, main = "CRG correlation",
                   filename = file.path(paths$res_cross, "CRG_correlation_heatmap.png"),
                   show_rownames = FALSE)

## 2) DE-CRG correlation network
decrg_mat <- combat[de_crgs, ]
decrg_edges <- correlation_network(decrg_mat, method = COR_METHOD, p_thr = COR_P_THRESHOLD)
write.csv(decrg_edges, file.path(paths$res_cross, "DECRG_correlation_edges.csv"), row.names = FALSE)
if (nrow(decrg_edges) > 0) {
  g <- graph_from_data_frame(decrg_edges[, c("from","to")], directed = FALSE)
  p_net <- ggraph(g, layout = "fr") + geom_edge_link() + geom_node_point(size = 4) +
    geom_node_text(aes(label = name), repel = TRUE) + labs(title = "DE-CRG correlation network") + pub_theme()
  save_publication_plot(p_net, file.path(paths$res_cross, "DECRG_network"), width = 7, height = 6)
}

## 3) top DEG-CRG correlation heatmap (CRGs vs top non-CRG DEGs)
sub <- combat[c(de_crgs, top_deg), ]
c_sub <- rcorr(t(sub), type = COR_METHOD)
pheatmap::pheatmap(c_sub$r, main = "DE-CRG vs top DEG correlation",
                   filename = file.path(paths$res_cross, "CRG_DEG_correlation_heatmap.png"),
                   show_rownames = FALSE)

## 4) significant positive/negative correlation network (DE-CRGs with top DEGs)
edges <- data.frame(
  from = character(), to = character(), r = numeric(), p = numeric(), stringsAsFactors = FALSE)
rmat <- c_sub$r; pmat <- c_sub$P
for (i in seq_len(nrow(rmat))) {
  for (j in seq_len(ncol(rmat))) {
    if (rownames(rmat)[i] %in% de_crgs && colnames(rmat)[j] %in% top_deg) {
      if (abs(rmat[i,j]) > 0.3 && pmat[i,j] < COR_P_THRESHOLD) {
        edges <- rbind(edges, data.frame(from = rownames(rmat)[i], to = colnames(rmat)[j],
                                         r = rmat[i,j], p = pmat[i,j], stringsAsFactors = FALSE))
      }
    }
  }
}
write.csv(edges, file.path(paths$res_cross, "DECRG_DEG_significant_edges.csv"), row.names = FALSE)

## 5) pairwise scatter for strongest CRG associations
if (nrow(crg_corr_edges) > 0) {
  top_edge <- crg_corr_edges[which.max(abs(crg_corr_edges$r)), ]
  d1 <- data.frame(x = as.vector(crg_mat[top_edge$from, ]), y = as.vector(crg_mat[top_edge$to, ]))
  colnames(d1) <- c(top_edge$from, top_edge$to)
  p_sc <- ggplot(d1, aes_string(x = top_edge$from, y = top_edge$to)) + geom_point(alpha = 0.5) +
    labs(title = paste("Strongest CRG correlation:", top_edge$from, "-", top_edge$to)) + pub_theme()
  save_publication_plot(p_sc, file.path(paths$res_cross, "strongest_CRG_scatter"))
}

saveRDS(list(crg_corr = c_crg, decrg_edges = decrg_edges, deg_assoc_edges = edges),
        file.path(paths$res_cross, "crosstalk_correlations.rds"))
cat("Crosstalk correlation analysis complete.\n")
cat(" CRG corr edges:", nrow(crg_corr_edges), " DE-CRG corr edges:", nrow(decrg_edges),
    " DE-CRG-DEG edges:", nrow(edges), "\n")
quit(save = "no", status = 0, runLast = FALSE)
