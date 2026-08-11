setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 06_differential_expression.R
# Differential expression analysis using limma
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")
source("R/plotting_functions.R")
load_required_packages()

combat_expr <- readRDS(file.path(DISCOVERY_DIR, "discovery_combat.rds"))
meta <- readRDS(file.path(DISCOVERY_DIR, "discovery_cohort.rds"))$metadata

cat("Running differential expression analysis...\n")
cat("Comparison: Tumor vs Normal\n")

# Create design matrix
group <- factor(meta$class, levels = c("Normal", "Tumor"))
design <- model.matrix(~ group)
colnames(design) <- levels(group)

# Fit linear model
fit <- lmFit(combat_expr, design)
fit <- eBayes(fit)

# Extract results
de_results <- topTable(fit, coef = "Tumor", number = Inf, sort.by = "P")
de_results$Gene <- rownames(de_results)

# Identify DEGs
de_results$DEG <- "NS"
de_results$DEG[de_results$adj.P.Val < ADJ_P_THRESHOLD & 
                de_results$logFC > LOG2FC_THRESHOLD] <- "Up"
de_results$DEG[de_results$adj.P.Val < ADJ_P_THRESHOLD & 
                de_results$logFC < -LOG2FC_THRESHOLD] <- "Down"

n_up <- sum(de_results$DEG == "Up")
n_down <- sum(de_results$DEG == "Down")
n_total <- n_up + n_down
cat(sprintf("\nDEGs found: %d up, %d down, %d total\n", n_up, n_down, n_total))

# Save results
write.csv(de_results, file.path(DE_DIR, "differential_expression_full.csv"), row.names = FALSE)
write.csv(de_results[de_results$DEG != "NS", ],
          file.path(DE_DIR, "significant_DEGs.csv"), row.names = FALSE)
write.csv(de_results[de_results$DEG == "Up", ],
          file.path(DE_DIR, "upregulated_genes.csv"), row.names = FALSE)
write.csv(de_results[de_results$DEG == "Down", ],
          file.path(DE_DIR, "downregulated_genes.csv"), row.names = FALSE)

# Ranked gene list for GSEA
ranked_list <- de_results$logFC
names(ranked_list) <- de_results$Gene
ranked_list <- sort(ranked_list, decreasing = TRUE)
saveRDS(ranked_list, file.path(DE_DIR, "ranked_gene_list.rds"))

# Volcano plot
p_volcano <- plot_volcano(de_results, title = "Discovery Cohort DEGs")
save_publication_plot(p_volcano, "volcano_plot", width = 9, height = 7)

# MA plot
de_results$AveExpr <- rowMeans(combat_expr[, meta$class == "Normal"], na.rm = TRUE)
p_ma <- ggplot(de_results, aes(x = AveExpr, y = logFC, color = DEG)) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_manual(values = c("Down" = "#377EB8", "NS" = "grey70", "Up" = "#E41A1C")) +
    geom_hline(yintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD), 
               linetype = "dashed", color = "grey40") +
    labs(title = "MA Plot", x = "Average Expression", y = "log2 Fold Change") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_publication_plot(p_ma, "MA_plot", width = 9, height = 7)

# DEG count bar plot
deg_counts <- data.frame(
    Direction = c("Upregulated", "Downregulated"),
    Count = c(n_up, n_down)
)
p_bar <- ggplot(deg_counts, aes(x = Direction, y = Count, fill = Direction)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8")) +
    labs(title = "Differentially Expressed Genes", x = "", y = "Count") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
save_publication_plot(p_bar, "DEG_counts", width = 7, height = 6)

# Save DEG list
deg_genes <- de_results$Gene[de_results$DEG != "NS"]
saveRDS(deg_genes, file.path(DE_DIR, "deg_gene_list.rds"))

cat("\nDifferential expression analysis complete.\n")
cat("Significant DEGs:", n_total, "\n")

