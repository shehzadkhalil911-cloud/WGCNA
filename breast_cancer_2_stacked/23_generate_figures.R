setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 23_generate_figures.R
# Generate all publication-quality figures
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")
source("R/plotting_functions.R")
source("R/evaluation_functions.R")
load_required_packages()

cat("Generating all publication figures...\n")

# --- Figure 1: Study design ---
cat("  Figure 1: Study design...\n")
suitability <- read.csv(file.path(COHORT_SUIT_DIR, "cohort_suitability.csv"))
suitability$Role <- factor(suitability$Role, levels = c("Discovery", "Validation"))

p1b <- ggplot(suitability, aes(x = Dataset, y = Tumor_Samples, fill = Role)) +
    geom_col(aes(y = Tumor_Samples), position = "dodge") +
    geom_col(aes(y = -Normal_Samples), position = "dodge", alpha = 0.7) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    annotate("text", x = 0.5, y = max(suitability$Tumor_Samples) * 0.9, 
             label = "Tumor", hjust = 0, size = 4) +
    annotate("text", x = 0.5, y = -max(suitability$Normal_Samples) * 0.5, 
             label = "Normal", hjust = 0, size = 4) +
    labs(title = "B. Sample Composition", x = "", y = "Sample Count") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
save_publication_plot(p1b, "Figure1B_sample_composition", width = 8, height = 6)

# --- Figure 2: Preprocessing ---
cat("  Figure 2: Preprocessing...\n")
# We'll create density plots from each dataset
all_gse <- c(DISCOVERY_GSE, VALIDATION_GSE)
density_data <- data.frame()
for (gse in all_gse) {
    proc <- readRDS(file.path(PROCESSED_DIR, paste0(gse, "_processed.rds")))
    expr <- proc$expression_norm
    for (s in sample(colnames(expr), min(10, ncol(expr)))) {
        density_data <- rbind(density_data, data.frame(
            Sample = s, Dataset = gse,
            Expression = as.numeric(expr[, s]),
            stringsAsFactors = FALSE
        ))
    }
}
p2c <- ggplot(density_data, aes(x = Expression, color = Sample)) +
    geom_density(linewidth = 0.5, show.legend = FALSE) +
    facet_wrap(~ Dataset, scales = "free_y") +
    labs(title = "C. Expression Density After Normalization",
         x = "Log2 Expression", y = "Density") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_publication_plot(p2c, "Figure2C_density_after_norm", width = 10, height = 8)

# --- Figure 3: Batch correction ---
cat("  Figure 3: Batch correction PCA...\n")
pca_before <- readRDS(file.path(PCA_DIR, "pca_before.rds"))
pca_after <- readRDS(file.path(PCA_DIR, "pca_after.rds"))
meta <- readRDS(file.path(DISCOVERY_DIR, "discovery_cohort.rds"))$metadata

make_pca_df <- function(pca_obj, meta_df) {
    data.frame(
        PC1 = pca_obj$x[, 1], PC2 = pca_obj$x[, 2],
        Dataset = meta_df$dataset, Class = meta_df$class,
        stringsAsFactors = FALSE
    )
}
df_before <- make_pca_df(pca_before, meta)
df_after <- make_pca_df(pca_after, meta)

p3a <- ggplot(df_before, aes(x = PC1, y = PC2, color = Dataset)) +
    geom_point(size = 2.5, alpha = 0.8) +
    theme_bw(base_size = 12) +
    labs(title = "A. Before Correction by Dataset") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p3b <- ggplot(df_after, aes(x = PC1, y = PC2, color = Dataset)) +
    geom_point(size = 2.5, alpha = 0.8) +
    theme_bw(base_size = 12) +
    labs(title = "B. After Correction by Dataset") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p3 <- p3a + p3b + plot_layout(ncol = 2)
save_publication_plot(p3, "Figure3_batch_correction_PCA", width = 14, height = 6)

# --- Figure 4: Differential expression ---
cat("  Figure 4: Differential expression...\n")
de_results <- read.csv(file.path(DE_DIR, "differential_expression_full.csv"))
de_results$DEG <- "NS"
de_results$DEG[de_results$adj.P.Val < ADJ_P_THRESHOLD & 
                de_results$logFC > LOG2FC_THRESHOLD] <- "Up"
de_results$DEG[de_results$adj.P.Val < ADJ_P_THRESHOLD & 
                de_results$logFC < -LOG2FC_THRESHOLD] <- "Down"

p4a <- ggplot(de_results, aes(x = logFC, y = -log10(adj.P.Val), color = DEG)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("Down" = "#377EB8", "NS" = "grey70", "Up" = "#E41A1C")) +
    geom_hline(yintercept = -log10(ADJ_P_THRESHOLD), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD), 
               linetype = "dashed", color = "grey40") +
    labs(title = "A. Volcano Plot", x = "log2 Fold Change", y = "-log10 Adjusted P-value") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_publication_plot(p4a, "Figure4A_volcano", width = 9, height = 7)

# --- Figure 5-9: Model-specific performance ---
cat("  Figures 5-9: Model performance plots...\n")
all_metrics <- readRDS(file.path(INTVAL_DIR, "holdout_metrics.rds"))

# Combined ROC curves figure
all_probs <- readRDS(file.path(INTVAL_DIR, "holdout_predictions.rds"))
observed <- readRDS(file.path(BASE_MODELS_DIR, "holdout_data.rds"))$Class
roc_list <- list()
for (model_name in colnames(all_probs)) {
    obs_bin <- ifelse(observed == POSITIVE_CLASS, 1, 0)
    roc_list[[model_name]] <- pROC::roc(obs_bin, all_probs[[model_name]], quiet = TRUE)
}

p_roc_combined <- plot_combined_roc(roc_list, names(roc_list), 
                                     "Hold-out ROC Comparison")
save_publication_plot(p_roc_combined, "Figure11A_combined_ROC", width = 9, height = 8)

# --- Figure 10: Consensus genes ---
cat("  Figure 10: Consensus genes...\n")
consensus <- readRDS(file.path(CONSENSUS_DIR, "consensus_genes.rds"))
algo_genes <- consensus$algorithm_genes
gene_counts <- data.frame(
    Algorithm = names(algo_genes),
    Count = sapply(algo_genes, length)
)
p10a <- ggplot(gene_counts, aes(x = reorder(Algorithm, -Count), y = Count, fill = Algorithm)) +
    geom_bar(stat = "identity") +
    labs(title = "A. Genes Selected by Each Algorithm", x = "", y = "Number of Genes") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
save_publication_plot(p10a, "Figure10A_gene_counts", width = 8, height = 6)

# --- Figure 12: Stacked model ---
cat("  Figure 12: Stacked model...\n")
stacked_model <- readRDS(file.path(STACKED_DIR, "stacked_elastic_net.rds"))
coefs <- as.matrix(coef(stacked_model$final_model, s = stacked_model$best_lambda))
coef_df <- data.frame(
    Base_Learner = rownames(coefs),
    Coefficient = as.numeric(coefs[, 1])
)
coef_df <- coef_df[coef_df$Coefficient != 0, ]
p12d <- ggplot(coef_df, aes(x = reorder(Base_Learner, Coefficient), 
                              y = Coefficient, fill = Coefficient > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = "D. Meta-Learner Coefficients", x = "", y = "Coefficient") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
save_publication_plot(p12d, "Figure12D_meta_learner_coefs", width = 8, height = 5)

# --- Figure 13: Bootstrap and DeLong ---
cat("  Figure 13: Bootstrap/DeLong...\n")
ci_table <- read.csv(file.path(BOOT_DIR, "bootstrap_ci_table.csv"))
auc_ci <- ci_table[ci_table$Metric == "roc_auc", ]
p13a <- ggplot(auc_ci, aes(x = reorder(Model, -Mean), y = Mean, color = Model)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
    labs(title = "A. Bootstrap AUC with 95% CI", x = "", y = "AUC") +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_manual(values = MODEL_COLORS)
save_publication_plot(p13a, "Figure13A_bootstrap_AUC", width = 8, height = 6)

# --- Figure 14: External validation ---
cat("  Figure 14: External validation...\n")
val_probs <- readRDS(file.path(EXTVAL_RES_DIR, "external_predictions.rds"))
val_data <- readRDS(file.path(EXTVAL_DIR, "validation_cohort.rds"))
val_meta <- val_data$metadata

# Validation ROC
val_roc_list <- list()
for (model_name in colnames(val_probs)) {
    obs_bin <- ifelse(val_meta$class == POSITIVE_CLASS, 1, 0)
    val_roc_list[[model_name]] <- pROC::roc(obs_bin, val_probs[[model_name]], quiet = TRUE)
}
p14d <- plot_combined_roc(val_roc_list, names(val_roc_list), 
                           "D. External Validation ROC Curves")
save_publication_plot(p14d, "Figure14D_external_ROC", width = 9, height = 8)

# --- Figure 15: Enrichment ---
cat("  Figure 15: Enrichment...\n")
go_bp_file <- file.path(ENRICH_DIR, "Consensus_GO_BP.csv")
if (file.exists(go_bp_file)) {
    go_bp <- read.csv(go_bp_file)
    if (nrow(go_bp) > 0) {
        top_go <- head(go_bp, 15)
        p15b <- ggplot(top_go, aes(x = reorder(Description, -log10(pvalue)), 
                                    y = -log10(pvalue), fill = Count)) +
            geom_bar(stat = "identity") +
            coord_flip() +
            labs(title = "B. GO Biological Process Enrichment", 
                 x = "", y = "-log10 P-value") +
            theme_bw(base_size = 12) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        save_publication_plot(p15b, "Figure15B_GO_enrichment", width = 10, height = 8)
    }
}

kegg_file <- file.path(ENRICH_DIR, "Consensus_KEGG.csv")
if (file.exists(kegg_file)) {
    kegg <- read.csv(kegg_file)
    if (nrow(kegg) > 0) {
        top_kegg <- head(kegg, 15)
        p15c <- ggplot(top_kegg, aes(x = reorder(Description, -log10(pvalue)), 
                                      y = -log10(pvalue), fill = Count)) +
            geom_bar(stat = "identity") +
            coord_flip() +
            labs(title = "C. KEGG Pathway Enrichment", 
                 x = "", y = "-log10 P-value") +
            theme_bw(base_size = 12) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        save_publication_plot(p15c, "Figure15C_KEGG_enrichment", width = 10, height = 8)
    }
}

cat("\nAll figures generated.\n")

