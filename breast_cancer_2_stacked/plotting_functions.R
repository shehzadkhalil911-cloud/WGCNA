# =============================================================================
# plotting_functions.R
# Publication-quality plotting functions
# =============================================================================

#' Save publication-quality plots
save_publication_plot <- function(
    plot,
    filename,
    width = 7,
    height = 6,
    results_dir = FIGURES_DIR
) {
    tiff_file <- file.path(results_dir, paste0(filename, ".tiff"))
    pdf_file  <- file.path(results_dir, paste0(filename, ".pdf"))

    ggplot2::ggsave(
        filename = tiff_file,
        plot = plot,
        width = width,
        height = height,
        units = "in",
        dpi = FIGURE_DPI,
        compression = "lzw",
        bg = "white"
    )

    ggplot2::ggsave(
        filename = pdf_file,
        plot = plot,
        width = width,
        height = height,
        units = "in",
        device = "pdf",
        bg = "white"
    )
}

#' Create PCA plot with customizable coloring
plot_pca <- function(expr_matrix, metadata, color_var, shape_var = NULL,
                     title = "PCA Plot", ellipse = FALSE) {
    pca <- prcomp(t(expr_matrix), center = TRUE, scale. = TRUE)
    pca_df <- as.data.frame(pca$x[, 1:3])
    pca_df$Sample <- colnames(expr_matrix)
    pca_df <- merge(pca_df, metadata, by.x = "Sample", by.y = "sample_id", all.x = TRUE)
    var_exp <- round(100 * pca$sdev[1:2]^2 / sum(pca$sdev^2), 1)

    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[color_var]])) +
        geom_point(size = 3, alpha = 0.8) +
        labs(
            title = title,
            x = paste0("PC1 (", var_exp[1], "%)"),
            y = paste0("PC2 (", var_exp[2], "%)")
        ) +
        theme_bw(base_size = 14) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "right"
        ) +
        scale_color_manual(values = c(Tumor = TUMOR_COLOR, Normal = NORMAL_COLOR))

    if (!is.null(shape_var)) {
        p <- p + aes(shape = .data[[shape_var]])
    }
    if (ellipse) {
        p <- p + stat_ellipse(aes(color = .data[[color_var]]), level = 0.95)
    }
    return(p)
}

#' Create volcano plot
plot_volcano <- function(de_results, highlight_genes = NULL, title = "Volcano Plot") {
    de_results$Significance <- "NS"
    de_results$Significance[de_results$adj.P.Val < ADJ_P_THRESHOLD &
                            de_results$logFC > LOG2FC_THRESHOLD] <- "Up"
    de_results$Significance[de_results$adj.P.Val < ADJ_P_THRESHOLD &
                            de_results$logFC < -LOG2FC_THRESHOLD] <- "Down"
    de_results$Significance <- factor(de_results$Significance,
                                      levels = c("Down", "NS", "Up"))

    p <- ggplot(de_results, aes(x = logFC, y = -log10(adj.P.Val))) +
        geom_point(aes(color = Significance), alpha = 0.6, size = 1.5) +
        scale_color_manual(values = c("Down" = "#377EB8", "NS" = "grey70", "Up" = "#E41A1C")) +
        geom_hline(yintercept = -log10(ADJ_P_THRESHOLD), linetype = "dashed", color = "grey40") +
        geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD),
                   linetype = "dashed", color = "grey40") +
        labs(title = title, x = "log2 Fold Change", y = "-log10 Adjusted P-value") +
        theme_bw(base_size = 14) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    if (!is.null(highlight_genes)) {
        hl <- de_results[de_results$Gene %in% highlight_genes, ]
        if (nrow(hl) > 0) {
            p <- p + geom_point(data = hl, color = "black", size = 3, shape = 1) +
                ggrepel::geom_text_repel(
                    data = hl,
                    aes(label = Gene),
                    size = 3,
                    max.overlaps = 20,
                    box.padding = 0.5
                )
        }
    }
    return(p)
}

#' Create correlation heatmap
plot_correlation_heatmap <- function(expr_matrix, metadata = NULL,
                                      title = "Sample Correlation") {
    cor_mat <- cor(expr_matrix, method = "spearman", use = "pairwise.complete.obs")
    p <- pheatmap(
        cor_mat,
        annotation_col = metadata,
        main = title,
        show_colnames = FALSE,
        show_rownames = FALSE,
        color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100),
        border_color = NA,
        silent = TRUE
    )
    return(p)
}

#' Create ROC plot for multiple models
plot_combined_roc <- function(roc_list, model_names, title = "ROC Curves") {
    df_list <- list()
    for (i in seq_along(roc_list)) {
        df <- data.frame(
            Specificity = rev(roc_list[[i]]$specificities),
            Sensitivity = rev(roc_list[[i]]$sensitivities),
            Model = model_names[i],
            AUC = round(roc_list[[i]]$auc, 3)
        )
        df_list[[i]] <- df
    }
    roc_df <- do.call(rbind, df_list)
    roc_df$label <- paste0(roc_df$Model, " (AUC = ", roc_df$AUC, ")")

    p <- ggplot(roc_df, aes(x = 1 - Specificity, y = Sensitivity,
                            color = Model)) +
        geom_line(linewidth = 1) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
        labs(title = title, x = "1 - Specificity", y = "Sensitivity") +
        theme_bw(base_size = 14) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "right"
        ) +
        scale_color_manual(values = MODEL_COLORS) +
        coord_fixed()
    return(p)
}

#' Create calibration plot
plot_calibration <- function(pred_df, true_class, model_name = "Model") {
    cal_data <- data.frame(
        pred = pred_df,
        obs = ifelse(true_class == POSITIVE_CLASS, 1, 0)
    )
    cal_curve <- loess(obs ~ pred, data = cal_data, span = 0.5)
    cal_x <- seq(0, 1, length.out = 100)
    cal_y <- predict(cal_curve, newdata = data.frame(pred = cal_x))

    p <- ggplot() +
        geom_line(aes(x = cal_x, y = cal_y), color = MODEL_COLORS[model_name],
                  linewidth = 1.2) +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
        xlim(0, 1) + ylim(0, 1) +
        labs(title = paste(model_name, "Calibration"),
             x = "Predicted Probability", y = "Observed Frequency") +
        theme_bw(base_size = 14) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    return(p)
}

#' Create confusion matrix heatmap
plot_confusion_matrix <- function(cm, title = "Confusion Matrix") {
    cm_df <- as.data.frame(cm$table)
    cm_df$Prediction <- factor(cm_df$Prediction, levels = rev(levels(cm_df$Prediction)))
    cm_df$Percent <- round(100 * cm_df$Freq / sum(cm_df$Freq), 1)

    p <- ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
        geom_tile() +
        geom_text(aes(label = paste0(Freq, "\n(", Percent, "%)")), size = 5) +
        scale_fill_gradient(low = "white", high = "#377EB8") +
        labs(title = title, x = "Actual", y = "Predicted") +
        theme_bw(base_size = 14) +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "none"
        )
    return(p)
}
