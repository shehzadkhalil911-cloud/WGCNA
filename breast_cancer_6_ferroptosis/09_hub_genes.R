source("R/utils.R")
load_libraries()
set.seed(123)

BASE_DIR <- "F:/independent projects/cancer/breast cancer"
OUT_DIR <- file.path(BASE_DIR, "Breast cancer 6 Ferroptosis")

log_msg("=== Step 9: Hub Gene Prioritization ===")

module_info <- read.csv(file.path(OUT_DIR, "results/WGCNA/module_assignments.csv"))
rownames(module_info) <- module_info$Gene
deg_table <- readRDS(file.path(OUT_DIR, "results/differential_expression/combined_DEG_results.rds"))
ferr_genes <- readRDS(file.path(OUT_DIR, "data_processed/ferroptosis_gene_list/ferroptosis_genes_final.rds"))
de_frg <- readRDS(file.path(OUT_DIR, "results/ferroptosis_genes/DE_FRG_results.rds"))
combat_expr <- readRDS(file.path(OUT_DIR, "data_processed/discovery/discovery_combat.rds"))
metadata <- read.csv(file.path(OUT_DIR, "data_processed/discovery/discovery_metadata.csv"))
ferr_score <- readRDS(file.path(OUT_DIR, "results/ferroptosis_score/ferroptosis_score_vector.rds"))

# Build evidence table
hub_genes <- data.frame(
  Gene = module_info$Gene,
  Module = module_info$Module,
  stringsAsFactors = FALSE
)

# 1. Differential expression evidence
hub_genes$DEG_significant <- hub_genes$Gene %in% rownames(deg_table)[deg_table$Significant]
hub_genes$DEG_logFC <- deg_table[hub_genes$Gene, "logFC"]
hub_genes$DEG_adjP <- deg_table[hub_genes$Gene, "adj.P.Val"]

# 2. Ferroptosis gene membership
hub_genes$Is_Ferroptosis_Gene <- hub_genes$Gene %in% ferr_genes
hub_genes$Is_DE_FRG <- hub_genes$Gene %in% rownames(de_frg)

# 3. Module membership (highest MM per gene across modules)
mm_cols <- grep("^ME", colnames(module_info), value = TRUE)
if (length(mm_cols) > 0) {
  mm_data <- as.matrix(module_info[, mm_cols])
  hub_genes$Max_MM <- apply(abs(mm_data), 1, max, na.rm = TRUE)
  # Which module gives the highest MM
  hub_genes$Best_Module <- mm_cols[apply(abs(mm_data), 1, which.max)]
}

# 4. Gene significance
hub_genes$GS_Tumor <- module_info$GS_Tumor
hub_genes$GS_Ferroptosis <- module_info$GS_Ferroptosis

# 5. Intramodular connectivity
conn_cols <- setdiff(colnames(module_info), c(colnames(hub_genes), mm_cols))
conn_cols <- intersect(conn_cols, unique(module_info$Module))
if (length(conn_cols) > 0) {
  conn_data <- as.matrix(module_info[, conn_cols, drop = FALSE])
  hub_genes$Intramodular_Connectivity <- apply(conn_data, 1, function(x) max(x, na.rm = TRUE))
}

# 6. Correlation with ferroptosis score
hub_genes$Ferroptosis_Correlation <- apply(combat_expr[hub_genes$Gene, , drop = FALSE], 1, function(x) {
  cor(x, ferr_score[colnames(combat_expr)], method = "spearman", use = "pairwise")
})

# 7. Consistent FC direction
fc_cols <- intersect(c("GSE42568", "GSE15852", "GSE29431"),
                     colnames(read.csv(file.path(OUT_DIR, "results/differential_expression/per_dataset_fold_changes.csv"),
                                       nrows = 1)))
if (length(fc_cols) >= 2) {
  fc_data <- read.csv(file.path(OUT_DIR, "results/differential_expression/per_dataset_fold_changes.csv"),
                       row.names = 1)
  fc_data <- fc_data[intersect(hub_genes$Gene, rownames(fc_data)), fc_cols, drop = FALSE]
  hub_genes$FC_Consistency <- ifelse(hub_genes$Gene %in% rownames(fc_data), 
    apply(sign(fc_data[hub_genes$Gene, , drop = FALSE]), 1, function(x) all(x == x[1], na.rm = TRUE)), NA)
} else {
  hub_genes$FC_Consistency <- NA
}

# Composite score
score_cols <- c("DEG_significant", "Is_Ferroptosis_Gene", "Is_DE_FRG",
                "Max_MM", "GS_Tumor", "GS_Ferroptosis", "Intramodular_Connectivity",
                "Ferroptosis_Correlation", "FC_Consistency")

for (col in score_cols) {
  if (!col %in% colnames(hub_genes)) next
  if (is.numeric(hub_genes[[col]])) {
    hub_genes[[paste0(col, "_scaled")]] <- scale(hub_genes[[col]])[, 1]
  }
}

# Simple composite (sum of standardized scores)
score_cols_scaled <- grep("_scaled$", colnames(hub_genes), value = TRUE)
if (length(score_cols_scaled) > 0) {
  hub_genes$Evidence_Score <- rowSums(abs(as.matrix(hub_genes[, score_cols_scaled])), na.rm = TRUE)
} else {
  hub_genes$Evidence_Score <- rowSums(hub_genes[, c("DEG_significant", "Is_Ferroptosis_Gene")] * 1, na.rm = TRUE)
}

hub_genes <- hub_genes[order(hub_genes$Evidence_Score, decreasing = TRUE), ]

# Top hub genes (top 30 by evidence score)
top_hubs <- head(hub_genes, 30)

write.csv(hub_genes, file.path(OUT_DIR, "results/WGCNA/hub_gene_evidence_table.csv"), row.names = FALSE)
write.csv(top_hubs, file.path(OUT_DIR, "results/WGCNA/top_hub_genes.csv"), row.names = FALSE)
writeLines(top_hubs$Gene, file.path(OUT_DIR, "results/WGCNA/final_hub_genes.txt"))

log_msg(sprintf("Top hub genes identified: %s",
                paste(head(top_hubs$Gene, 10), collapse = ", ")))

# Evidence heatmap for top hubs
if (nrow(top_hubs) > 1) {
  evidence_cols <- intersect(c("DEG_significant", "Is_Ferroptosis_Gene", "Is_DE_FRG",
                                "FC_Consistency", "Max_MM", "GS_Tumor_adj",
                                "Intramodular_Connectivity", "Ferroptosis_Correlation"),
                              colnames(top_hubs))

  heat_data <- top_hubs[, evidence_cols, drop = FALSE]
  heat_data_numeric <- as.data.frame(lapply(heat_data, function(x) {
    if (is.logical(x)) as.numeric(x) else x
  }))
  rownames(heat_data_numeric) <- top_hubs$Gene

  # Handle NA values
  heat_data_numeric[is.na(heat_data_numeric)] <- 0

  if (ncol(heat_data_numeric) > 1 && nrow(heat_data_numeric) > 1) {
    open_tiff(file.path(OUT_DIR, "results/figures/Figure7A_hub_evidence_heatmap.tiff"), width = 8, height = 7)
    pheatmap::pheatmap(as.matrix(heat_data_numeric),
                       scale = "column", display_numbers = TRUE,
                       fontsize = 7, main = "Hub Gene Evidence Matrix",
                       color = colorRampPalette(c("blue", "white", "red"))(100))
    dev.off()
  }
}

# Candidate gene expression box plots
if (length(top_hubs$Gene) > 0) {
  top_genes <- top_hubs$Gene[1:min(9, nrow(top_hubs))]
  expr_long <- reshape2::melt(as.matrix(combat_expr[top_genes, ]))
  colnames(expr_long) <- c("Gene", "Sample", "Expression")
  expr_long$Group <- metadata$group[match(expr_long$Sample, metadata$sample_id)]
  expr_long$Gene <- factor(expr_long$Gene, levels = top_genes)

  p_hub_box <- ggplot(expr_long, aes(Group, Expression, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    facet_wrap(~ Gene, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c(Tumor = "#E41A1C", Normal = "#377EB8")) +
    labs(title = "Hub Gene Expression (Tumor vs Normal)") +
    theme_minimal() + theme(strip.text = element_text(size = 7))
  save_tiff_600dpi(p_hub_box, file.path(OUT_DIR, "results/figures/Figure7B_hub_boxplots.tiff"),
                   width = 9, height = 7)
}

# Ranking plot
p_rank <- ggplot(head(hub_genes, 50), aes(reorder(Gene, Evidence_Score), Evidence_Score)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  coord_flip() +
  labs(title = "Top 50 Hub Gene Rankings",
       x = "Gene", y = "Evidence Score") +
  theme_minimal()
save_tiff_600dpi(p_rank, file.path(OUT_DIR, "results/figures/Figure7C_hub_ranking.tiff"),
                 width = 7, height = 8)

# Consensus gene heatmap
if (length(top_genes) > 1) {
  mat_top <- combat_expr[top_genes, ]
  mat_scaled <- t(scale(t(mat_top)))
  anno_col <- data.frame(Group = metadata$group, row.names = colnames(combat_expr))

  open_tiff(file.path(OUT_DIR, "results/figures/Figure7D_consensus_heatmap.tiff"), width = 9, height = 6)
  pheatmap::pheatmap(mat_scaled, annotation_col = anno_col,
                     show_rownames = TRUE, show_colnames = FALSE,
                     fontsize_row = 7, cutree_cols = 2,
                     main = "Top Hub Gene Expression")
  dev.off()
}

saveRDS(top_hubs, file.path(OUT_DIR, "results/WGCNA/top_hub_genes.rds"))

log_msg("Step 9 complete: Hub gene prioritization finished")
