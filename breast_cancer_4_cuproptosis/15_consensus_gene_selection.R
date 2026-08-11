# 15_consensus_gene_selection.R
source("config/analysis_config.R")
source("R/load_functions.R")
ensure_paths()
suppressMessages({ library(ggplot2); library(UpSetR); library(ggpubr); library(dplyr); library(circlize) })

lasso <- readRDS(file.path(paths$res_ml, "LASSO", "lasso_result.rds"))
svm <- readRDS(file.path(paths$res_ml, "SVM", "svm_result.rds"))
rf <- readRDS(file.path(paths$res_ml, "Random_Forest", "rf_result.rds"))
deg <- readRDS(file.path(paths$res_de, "discovery_DEG_full.rds"))
crg_sum <- readRDS(file.path(paths$res_cross, "DE_CRG_summary.rds"))
pool <- crg_sum$candidate_pool

lasso_genes <- lasso$selected_min
svm_genes <- svm$selected
rf_genes <- rf$top_genes

## method-by-gene matrix
m <- matrix(0, nrow = length(pool), ncol = 3, dimnames = list(pool, c("LASSO","SVM","RF")))
m[rownames(m) %in% lasso_genes, "LASSO"] <- 1
m[rownames(m) %in% svm_genes, "SVM"] <- 1
m[rownames(m) %in% rf_genes, "RF"] <- 1
ml_count <- rowSums(m)

## PPI degree
ppi_cent <- tryCatch(read.csv(file.path(paths$res_ppi, "PPI_centrality.csv"), stringsAsFactors = FALSE), error = function(e) NULL)
ppi_deg <- if (!is.null(ppi_cent)) setNames(ppi_cent$degree, ppi_cent$gene) else rep(0, length(pool))
names(ppi_deg) <- sub(".*\\|", "", names(ppi_deg))  # string ids may contain |

## DE status
deg_status <- pool %in% deg$symbol[deg$sig]
names(deg_status) <- pool
crg_status <- rep(1, length(pool)); names(crg_status) <- pool
lfc <- setNames(deg$log2FC[match(pool, deg$symbol)], pool)
adjp <- setNames(deg$adj.P.Val[match(pool, deg$symbol)], pool)

## evidence table
ev <- data.frame(
  gene = pool,
  DEG_status = deg_status,
  log2FC = round(lfc, 3),
  adj_P = signif(adjp, 3),
  CRG_status = crg_status,
  LASSO = m[, "LASSO"],
  SVM = m[, "SVM"],
  RF = m[, "RF"],
  ML_count = ml_count,
  PPI_degree = round(ppi_deg[match(pool, names(ppi_deg))], 2),
  stringsAsFactors = FALSE
)
ev$consensus_strict <- ml_count >= 3
ev$consensus_2ml <- ml_count >= CONSENSUS_ML_MIN
write.csv(ev, file.path(paths$res_ml, "consensus", "evidence_table.csv"), row.names = FALSE)

## final locked signature: >= CONSENSUS_ML_MIN (2) ML methods
final_genes <- pool[ev$consensus_2ml]
if (length(final_genes) == 0) final_genes <- pool[ev$consensus_strict]
if (length(final_genes) == 0) final_genes <- pool   # fallback to full candidate pool
saveRDS(list(final_genes = final_genes, evidence = ev, method_matrix = m,
             lasso_genes = lasso_genes, svm_genes = svm_genes, rf_genes = rf_genes),
        file.path(paths$res_ml, "consensus", "consensus_genes.rds"))
write.csv(data.frame(gene = final_genes), file.path(paths$res_ml, "consensus", "final_consensus_genes.csv"), row.names = FALSE)

## plots
sel_count <- data.frame(method = c("LASSO","SVM","RF"), n = c(length(lasso_genes), length(svm_genes), length(rf_genes)))
p_bar <- ggplot(sel_count, aes(x = method, y = n)) + geom_bar(stat = "identity", fill = "#4E79A7") +
  labs(title = "Selected feature count by method", y = "n genes") + pub_theme()
save_publication_plot(p_bar, file.path(paths$res_ml, "consensus", "selected_count"))

png(file.path(paths$res_ml, "consensus", "venn.png"), width = 6, height = 6, units = "in", res = FIG_DPI)
ggvenn::ggvenn(list(LASSO = lasso_genes, SVM = svm_genes, RF = rf_genes))
dev.off()
tiff(file.path(paths$res_ml, "consensus", "venn.tiff"), width = 6, height = 6, units = "in", res = 600, compression = "lzw")
ggvenn::ggvenn(list(LASSO = lasso_genes, SVM = svm_genes, RF = rf_genes))
dev.off()
pdf(file.path(paths$res_ml, "consensus", "venn.pdf"), width = 6, height = 6)
ggvenn::ggvenn(list(LASSO = lasso_genes, SVM = svm_genes, RF = rf_genes))
dev.off()

melt_df <- data.frame(gene = rep(rownames(m), ncol(m)),
                        method = rep(colnames(m), each = nrow(m)),
                        value = as.vector(m))
p_heat <- ggplot(melt_df, aes(x = method, y = gene, fill = as.factor(value))) + geom_tile() +
  scale_fill_manual(values = c("white","#59A14F")) + labs(title = "Method-by-gene selection", fill = "selected") + pub_theme()
save_publication_plot(p_heat, file.path(paths$res_ml, "consensus", "method_gene_heatmap"), width = 5, height = 6)

cat("Consensus selection complete.\n")
cat(" LASSO:", length(lasso_genes), " SVM:", length(svm_genes), " RF(top):", length(rf_genes), "\n")
cat(" Final consensus genes (>=", CONSENSUS_ML_MIN, "ML methods):", length(final_genes), "->", paste(final_genes, collapse = ", "), "\n")
quit(save = "no", status = 0, runLast = FALSE)
