setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
source("config/analysis_config.R")
source("R/package_functions.R")
source("R/plotting_functions.R")
source("R/annotation_functions.R")
load_required_packages()

save_fig <- function(name, w, h, expr) {
    pdf(file.path(FIGURES_DIR, paste0(name, ".pdf")), width = w, height = h)
    eval(substitute(expr))
    dev.off()
    tiff(file.path(FIGURES_DIR, paste0(name, ".tiff")), width = w, height = h,
         units = "in", res = FIGURE_DPI, compression = "lzw")
    eval(substitute(expr))
    dev.off()
}

cat("Supplementary figures...\n")

# ----------------------------------------------------------------------
# 1. Batch correction barplot — variance explained
# ----------------------------------------------------------------------
cat("  1. Batch correction barplot\n")
pca_before <- readRDS(file.path(PCA_DIR, "pca_before.rds"))
pca_after  <- readRDS(file.path(PCA_DIR, "pca_after.rds"))
meta <- readRDS(file.path(DISCOVERY_DIR, "discovery_cohort.rds"))$metadata

make_var_df <- function(pca_obj, meta_df) {
    ve <- pca_obj$sdev^2 / sum(pca_obj$sdev^2)
    pc_df <- data.frame(PC1 = pca_obj$x[, 1], PC2 = pca_obj$x[, 2],
                        dataset = meta_df$dataset, class = meta_df$class)
    r2_batch <- summary(lm(PC1 ~ dataset, pc_df))$r.squared * ve[1] +
                summary(lm(PC2 ~ dataset, pc_df))$r.squared * ve[2]
    r2_class <- summary(lm(PC1 ~ class, pc_df))$r.squared * ve[1] +
                summary(lm(PC2 ~ class, pc_df))$r.squared * ve[2]
    data.frame(Source = c("Batch (Dataset)", "Class (Tumor vs Normal)"),
               Variance = c(r2_batch, r2_class))
}
df_before <- make_var_df(pca_before, meta)
df_after  <- make_var_df(pca_after, meta)
df_before$Condition <- "Before Correction"
df_after$Condition  <- "After Correction"
var_df <- rbind(df_before, df_after)

p_batch_bar <- ggplot(var_df, aes(x = Source, y = Variance, fill = Condition)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.6) +
    scale_fill_manual(values = c("Before Correction" = "#E41A1C",
                                  "After Correction" = "#377EB8")) +
    labs(title = "Variance Explained by Batch vs Class",
         x = "", y = "Weighted R-squared (PC1+PC2)") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "top")
save_fig("Supplementary_batch_correction_barplot", 8, 5, print(p_batch_bar))

# ---- 1b. Batch correction boxplot — expression distributions ----
cat("  1b. Batch correction boxplot\n")
disc <- readRDS(file.path(DISCOVERY_DIR, "discovery_cohort.rds"))
expr_before <- disc$expression
expr_after  <- readRDS(file.path(DISCOVERY_DIR, "discovery_combat.rds"))

make_box_df <- function(expr_mat, meta_df, condition) {
    medians <- apply(expr_mat, 2, median, na.rm = TRUE)
    data.frame(Sample = names(medians), MedianExpr = medians,
               Dataset = meta_df$dataset, Class = meta_df$class,
               Condition = condition, stringsAsFactors = FALSE)
}
box_before <- make_box_df(expr_before, meta, "Before Correction")
box_after  <- make_box_df(expr_after,  meta, "After Correction")
box_df <- rbind(box_before, box_after)
box_df$Condition <- factor(box_df$Condition, c("Before Correction", "After Correction"))

p_box <- ggplot(box_df, aes(x = Dataset, y = MedianExpr, fill = Dataset)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~ Condition, ncol = 2) +
    labs(title = "Sample Median Expression Before vs After ComBat",
         x = "", y = "Median log2 Expression") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1))
save_fig("Supplementary_batch_correction_boxplot", 10, 6, print(p_box))

# ---- 1c. PCA with color = dataset, shape = class ----
cat("  1c. PCA color/shape\n")
make_pca_plot <- function(pca_obj, meta_df, title_prefix, var_exp) {
    df <- data.frame(
        PC1 = pca_obj$x[, 1], PC2 = pca_obj$x[, 2],
        Dataset = meta_df$dataset, Class = meta_df$class)
    ggplot(df, aes(x = PC1, y = PC2, color = Dataset, shape = Class)) +
        geom_point(size = 3, alpha = 0.85) +
        stat_ellipse(aes(group = Dataset), level = 0.68, linewidth = 0.6, show.legend = FALSE) +
        scale_shape_manual(values = c("Tumor" = 16, "Normal" = 17)) +
        labs(title = paste0(title_prefix, " (PC1: ", var_exp[1], "%, PC2: ", var_exp[2], "%)"),
             x = paste0("PC1 (", var_exp[1], "%)"),
             y = paste0("PC2 (", var_exp[2], "%)")) +
        theme_bw(base_size = 13) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
}
ve_before <- round(100 * pca_before$sdev^2 / sum(pca_before$sdev^2), 2)
ve_after  <- round(100 * pca_after$sdev^2  / sum(pca_after$sdev^2),  2)
p_before <- make_pca_plot(pca_before, meta, "Before Batch Correction", ve_before)
p_after  <- make_pca_plot(pca_after,  meta, "After Batch Correction",  ve_after)
p_pca <- p_before + p_after + patchwork::plot_layout(ncol = 2) +
    patchwork::plot_annotation(title = "PCA: Color by Dataset, Shape by Class",
                               theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
save_fig("Supplementary_PCA_color_by_dataset_shape_by_class", 14, 6, print(p_pca))

# ----------------------------------------------------------------------
# 2. Top DEGs heatmap
# ----------------------------------------------------------------------
cat("  2. Top DEGs heatmap\n")
de <- read.csv(file.path(DE_DIR, "differential_expression_full.csv"))
de <- de[order(de$adj.P.Val), ]
top50 <- head(de$Gene, 50)
combat_expr <- readRDS(file.path(DISCOVERY_DIR, "discovery_combat.rds"))
common <- intersect(top50, rownames(combat_expr))
heat_data <- combat_expr[common, , drop = FALSE]
zmat <- t(scale(t(as.matrix(heat_data))))

anno_col <- data.frame(
    Class   = meta$class,
    Dataset = meta$dataset,
    row.names = colnames(heat_data)
)
anno_colors <- list(
    Class   = c("Tumor" = TUMOR_COLOR, "Normal" = NORMAL_COLOR),
    Dataset = c(stats::setNames(c("#66C2A5", "#FC8D62", "#8DA0CB"), DISCOVERY_GSE))
)
ht <- pheatmap(zmat, annotation_col = anno_col, annotation_colors = anno_colors,
               show_rownames = TRUE, show_colnames = FALSE,
               cluster_rows = TRUE, cluster_cols = TRUE,
               color = colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(100),
               border_color = NA, fontsize_row = 7,
               main = "Top 50 DEGs — Expression Heatmap")
save_fig("Supplementary_DEG_heatmap_top50", 12, 10, print(ht))

# ----------------------------------------------------------------------
# 3. Variable importance — each model separately
# ----------------------------------------------------------------------
cat("  3. Variable importance plots\n")
base_models <- readRDS(file.path(BASE_MODELS_DIR, "base_models.rds"))
consensus <- readRDS(file.path(CONSENSUS_DIR, "consensus_genes.rds"))
final_genes <- consensus$final_consensus
train_df <- readRDS(file.path(BASE_MODELS_DIR, "training_data.rds"))

safe_varimp <- function(model, name, n_top = 15) {
    imp <- tryCatch({
        if (name == "GB") {
            s <- summary(model$finalModel, plotit = FALSE)
            stats::setNames(s$rel.inf, s$var)
        } else {
            imp_df <- caret::varImp(model)$importance
            stats::setNames(imp_df[, 1], rownames(imp_df))
        }
    }, error = function(e) NULL)
    if (is.null(imp) || length(imp) == 0) return(NULL)
    imp <- sort(imp, decreasing = TRUE)
    imp_df <- data.frame(
        Gene = factor(names(imp), levels = rev(names(imp))),
        Importance = as.numeric(imp))
    head(imp_df, n_top)
}

for (mn in names(base_models)) {
    disp <- switch(mn, rf = "RF", gbm = "GB", svm = "SVM", knn = "KNN", nnet = "NN")
    cat("    ", disp, "\n")
    vi <- safe_varimp(base_models[[mn]], disp)
    if (is.null(vi)) next
    p <- ggplot(vi, aes(x = Importance, y = Gene)) +
        geom_bar(stat = "identity", fill = MODEL_COLORS[disp], width = 0.7) +
        labs(title = paste(disp, "— Top Variable Importance"), x = "Importance", y = "") +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    save_fig(paste0("Supplementary_importance_", disp), 7, max(4, nrow(vi) * 0.3 + 1), print(p))
}

# KNN permutation importance
cat("    KNN (permutation)\n")
knn_pred_fun <- function(x, y) {
    pred <- suppressWarnings(kknn::kknn(y ~ ., data.frame(x, y = y), data.frame(x), k = 7))
    pred$prob[, 2]
}
knn_auc_orig <- pROC::roc(ifelse(train_df$Class == POSITIVE_CLASS, 1, 0),
                          knn_pred_fun(train_df[, final_genes], train_df$Class), quiet = TRUE)$auc
knn_imp <- sapply(final_genes, function(g) {
    perm <- train_df
    perm[[g]] <- sample(perm[[g]])
    pROC::roc(ifelse(train_df$Class == POSITIVE_CLASS, 1, 0),
              knn_pred_fun(perm[, final_genes], perm$Class), quiet = TRUE)$auc
})
knn_imp <- knn_auc_orig - knn_imp
knn_imp <- sort(knn_imp, decreasing = TRUE)
knn_top <- data.frame(
    Gene = factor(names(head(knn_imp, 15)), levels = rev(names(head(knn_imp, 15)))),
    Importance = as.numeric(head(knn_imp, 15)))
p <- ggplot(knn_top, aes(x = Importance, y = Gene)) +
    geom_bar(stat = "identity", fill = MODEL_COLORS["KNN"], width = 0.7) +
    labs(title = "KNN — Permutation Importance (AUC drop)", x = "AUC Drop", y = "") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_fig("Supplementary_importance_KNN", 7, 5, print(p))

# ----------------------------------------------------------------------
# 4. Combined importance heatmap
# ----------------------------------------------------------------------
cat("  4. Combined importance heatmap\n")
all_imp <- data.frame(row.names = final_genes)
for (mn in names(base_models)) {
    disp <- switch(mn, rf = "RF", gbm = "GB", svm = "SVM", knn = "KNN", nnet = "NN")
    if (disp == "KNN") {
        imp_vec <- knn_imp[final_genes]
        imp_vec[is.na(imp_vec)] <- 0
    } else {
        vi <- safe_varimp(base_models[[mn]], disp)
        if (is.null(vi)) { all_imp[[disp]] <- 0; next }
        imp_vec <- stats::setNames(vi$Importance, as.character(vi$Gene))
    }
    all_imp[[disp]] <- imp_vec[final_genes]
    all_imp[[disp]][is.na(all_imp[[disp]])] <- 0
}
imp_mat <- as.matrix(all_imp)
# Remove rows with zero variance (all same value) — these become NaN after scaling
row_sd <- apply(imp_mat, 1, sd, na.rm = TRUE)
imp_mat <- imp_mat[row_sd > 0, , drop = FALSE]
imp_mat <- t(scale(t(imp_mat)))
imp_mat[is.nan(imp_mat)] <- 0

ht <- ComplexHeatmap::Heatmap(imp_mat, name = "Importance\n(Z-score)",
    cluster_rows = TRUE, cluster_columns = TRUE,
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    column_title = "Variable Importance Across Models (Z-scored per Gene)")

save_fig("Supplementary_importance_combined_heatmap", 8, 14, {
    ComplexHeatmap::draw(ht)
})

# ----------------------------------------------------------------------
# 5. Venn intersection — 5 models
# ----------------------------------------------------------------------
cat("  5. Venn intersection\n")
algo_genes <- consensus$algorithm_genes

if (requireNamespace("ggvenn", quietly = TRUE)) {
    p <- ggvenn::ggvenn(algo_genes, stroke_color = "white",
                        fill_color = MODEL_COLORS[names(algo_genes)],
                        set_name_size = 4, text_size = 3) +
        labs(title = "Gene Selection Overlap Across Models") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    save_fig("Supplementary_ggvenn_5way", 8, 7, print(p))
} else if (requireNamespace("ggVennDiagram", quietly = TRUE)) {
    p <- ggVennDiagram::ggVennDiagram(algo_genes, color = "white") +
        scale_fill_gradient(low = "grey90", high = "#E41A1C") +
        labs(title = "Gene Selection Overlap Across Models") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    save_fig("Supplementary_ggvenn_5way", 8, 7, print(p))
} else {
    up_df <- data.frame(Gene = unique(unlist(algo_genes)),
        RF  = unique(unlist(algo_genes)) %in% algo_genes[["RF"]],
        GB  = unique(unlist(algo_genes)) %in% algo_genes[["GB"]],
        SVM = unique(unlist(algo_genes)) %in% algo_genes[["SVM"]],
        KNN = unique(unlist(algo_genes)) %in% algo_genes[["KNN"]],
        NN  = unique(unlist(algo_genes)) %in% algo_genes[["NN"]])
    write.csv(up_df, file.path(FIGURES_DIR, "gene_selection_overlap_matrix.csv"), row.names = FALSE)
    cat("    Venn packages unavailable — saved overlap matrix CSV\n")
}

# ----------------------------------------------------------------------
# 6. Dotplot enrichment — BP, CC, MF, KEGG combined
# ----------------------------------------------------------------------
cat("  6. Enrichment dotplot\n")
enrich_obj <- readRDS(file.path(ENRICH_DIR, "Consensus_enrichment.rds"))

cat_names <- c("go_BP" = "GO Biological Process", "go_CC" = "GO Cellular Component",
               "go_MF" = "GO Molecular Function", "KEGG" = "KEGG Pathways")
dot_plots <- list()
for (ont in c("go_BP", "go_CC", "go_MF", "KEGG")) {
    if (ont == "KEGG") {
        obj <- enrich_obj$kegg
    } else {
        obj <- enrich_obj[[ont]]
    }
    if (is.null(obj) || nrow(obj@result) == 0) next
    df <- obj@result[obj@result$p.adjust < 0.05, , drop = FALSE]
    if (nrow(df) == 0) df <- head(obj@result, 10)
    df <- head(df, 10)
    df$logP <- -log10(df$p.adjust)
    df$Description <- factor(df$Description, levels = rev(df$Description))
    p <- ggplot(df, aes(x = Count, y = Description, size = logP)) +
        geom_point(color = switch(ont, go_BP = "#E41A1C", go_CC = "#377EB8",
                                  go_MF = "#4DAF4A", "KEGG" = "#984EA3")) +
        scale_size_continuous(range = c(2, 7)) +
        labs(title = cat_names[ont], x = "Count", y = "", size = "-log10(p.adj)") +
        theme_bw(base_size = 11) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
              axis.text.y = element_text(size = 8))
    dot_plots[[ont]] <- p
}
if (length(dot_plots) > 0) {
    ncol <- min(length(dot_plots), 2)
    nrow <- ceiling(length(dot_plots) / ncol)
    p_grid <- patchwork::wrap_plots(dot_plots, ncol = ncol) +
        patchwork::plot_annotation(title = "Enrichment Analysis — Consensus Genes",
                                   theme = theme(plot.title = element_text(hjust = 0.5, face = "bold")))
    save_fig("Supplementary_dotplot_enrichment_side_by_side", 14, 5 * nrow, print(p_grid))
} else {
    cat("    No enrichment results for dotplot\n")
}

# ----------------------------------------------------------------------
# 7. Regulation barplot — consensus genes
# ----------------------------------------------------------------------
cat("  7. Regulation barplot\n")
de_results <- read.csv(file.path(DE_DIR, "differential_expression_full.csv"))
de_map <- stats::setNames(de_results$logFC, de_results$Gene)
fc_vec <- de_map[final_genes]
fc_vec[is.na(fc_vec)] <- 0
fc_df <- data.frame(
    Gene = factor(final_genes, levels = rev(final_genes[order(fc_vec)])),
    logFC = fc_vec,
    Regulation = ifelse(fc_vec > 0, "Up", "Down"))

p_fc <- ggplot(fc_df, aes(x = logFC, y = Gene, fill = Regulation)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8")) +
    labs(title = "Consensus Genes — log2 Fold Change",
         x = "log2 Fold Change (Tumor vs Normal)", y = "") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "top", axis.text.y = element_text(size = 7))
save_fig("Supplementary_consensus_genes_regulation_barplot", 9, 14, print(p_fc))

# ----------------------------------------------------------------------
cat("\nAll supplementary figures saved to results/figures/\n")
