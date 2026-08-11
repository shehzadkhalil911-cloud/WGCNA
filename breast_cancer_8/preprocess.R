# ============================================================
# Preprocessing helpers
# ============================================================

#' Parse a GEO SOFT platform table
parse_soft_platform <- function(soft_file, max_lines = 800000) {
  con <- gzfile(soft_file, "rt")
  L <- readLines(con, n = max_lines)
  close(con)
  beg <- grep("^!platform_table_begin", L)
  end <- grep("^!platform_table_end", L)
  if (!length(beg) || !length(end)) stop("Could not find platform table in ", soft_file)
  hdr <- strsplit(L[beg + 1], "\t")[[1]]
  tab <- read.table(text = L[(beg + 2):(end - 1)], sep = "\t", header = FALSE,
                    quote = "\"", stringsAsFactors = FALSE, fill = TRUE,
                    na.strings = c("NA", ""))
  colnames(tab) <- hdr[1:ncol(tab)]
  tab
}

#' Validate gene symbols against org.Hs.eg.db
validate_symbols <- function(syms) {
  map <- AnnotationDbi::mappedkeys(org.Hs.egSYMBOL2EG)
  syms %in% map
}

#' Extract a single official symbol from an annotation cell
#' Handles "A /// B" multi-symbol strings by keeping the first.
first_symbol <- function(x) {
  out <- sub("\\s*///.*$", "", x)
  out <- trimws(out)
  out[out == ""] <- NA_character_
  out
}

#' Quantile normalize a gene x sample matrix (samples in columns)
quantile_normalize <- function(mat) {
  requireNamespace("preprocessCore", quietly = TRUE)
  qn <- preprocessCore::normalize.quantiles(as.matrix(mat))
  dimnames(qn) <- dimnames(mat)
  qn
}

#' Row-median imputation of residual missing values
impute_row_median <- function(mat) {
  for (i in seq_len(nrow(mat))) {
    r <- mat[i, ]
    nas <- is.na(r)
    if (any(nas)) {
      md <- stats::median(r[!nas])
      if (!is.finite(md)) md <- 0
      mat[i, nas] <- md
    }
  }
  mat
}

#' Assess sample outliers based on mean pairwise correlation
assess_outliers <- function(mat) {
  cm <- stats::cor(mat)
  diag(cm) <- NA
  mean_cor <- apply(cm, 2, function(x) mean(x, na.rm = TRUE))
  m <- stats::median(mean_cor, na.rm = TRUE)
  s <- stats::sd(mean_cor, na.rm = TRUE)
  flag <- mean_cor < (m - 3 * s)
  data.frame(sample_id = colnames(mat), mean_pairwise_correlation = mean_cor,
             is_outlier = flag)
}

#' Plot QC panels for a dataset (boxplots + densities)
qc_box_density <- function(mat_before, mat_after, sample_groups, tag) {
  melt1 <- function(m, lab) {
    df <- as.data.frame(m)
    df$gene <- rownames(m)
    tidyr::pivot_longer(df, -gene, names_to = "sample_id", values_to = "value") |>
      dplyr::mutate(stage = lab)
  }
  b <- melt1(mat_before, "Before")
  a <- melt1(mat_after, "After")
  dat <- rbind(b, a)
  dat$stage <- factor(dat$stage, levels = c("Before", "After"))
  dat$class <- sample_groups[dat$sample_id]

  p_box <- ggplot(dat, aes(x = sample_id, y = value, fill = class)) +
    geom_boxplot(outlier.size = 0.2, lwd = 0.15) +
    facet_wrap(~ stage, scales = "free_y") +
    scale_fill_manual(values = c(Tumor = COL_TUMOR, Normal = COL_NORMAL)) +
    theme_bc(base_size = 8) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "bottom") +
    labs(x = "Sample", y = "Expression",
         title = sprintf("%s expression box plots", tag))

  p_dens <- ggplot(dat, aes(x = value, colour = class)) +
    geom_density(lwd = 0.4) +
    facet_wrap(~ stage, scales = "free") +
    scale_colour_manual(values = c(Tumor = COL_TUMOR, Normal = COL_NORMAL)) +
    theme_bc(base_size = 8) +
    labs(x = "Expression", y = "Density",
         title = sprintf("%s expression density plots", tag))
  list(box = p_box, dens = p_dens)
}

#' PCA plot for a gene x sample matrix
pca_plot <- function(mat, sample_groups, tag, pcs = c(1, 2)) {
  rv <- apply(mat, 1, stats::var, na.rm = TRUE)
  mat <- mat[is.finite(rv) & rv > 0, , drop = FALSE]
  p <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  var_expl <- summary(p)$importance[2, ] * 100
  df <- data.frame(PC1 = p$x[, pcs[1]], PC2 = p$x[, pcs[2]],
                   class = sample_groups[colnames(mat)],
                   sample_id = colnames(mat))
  ggplot(df, aes(x = PC1, y = PC2, colour = class)) +
    geom_point(size = 1.6, alpha = 0.75) +
    scale_colour_manual(values = c(Tumor = COL_TUMOR, Normal = COL_NORMAL)) +
    theme_bc() +
    labs(x = sprintf("PC%d (%.1f%%)", pcs[1], var_expl[pcs[1]]),
         y = sprintf("PC%d (%.1f%%)", pcs[2], var_expl[pcs[2]]),
         title = sprintf("%s PCA", tag), colour = "Class")
}

#' Correlation heatmap of samples
correlation_heatmap <- function(mat, sample_groups, filename, tag) {
  cm <- stats::cor(mat)
  grp_col <- ifelse(sample_groups[colnames(mat)] == "Tumor", COL_TUMOR, COL_NORMAL)
  h <- ComplexHeatmap::Heatmap(
    cm, name = "Pearson r",
    show_row_names = FALSE, show_column_names = FALSE,
    column_split = factor(sample_groups[colnames(mat)], levels = CLASS_LEVELS),
    row_split = factor(sample_groups[colnames(mat)], levels = CLASS_LEVELS),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      class = sample_groups[colnames(mat)],
      col = list(class = c(Tumor = COL_TUMOR, Normal = COL_NORMAL))),
    col = circlize::colorRamp2(c(0, 0.5, 1), c("#2166AC", "#F7F7F7", "#B2182B")),
    use_raster = TRUE, raster_quality = 3,
    heatmap_legend_param = list(title_gp = grid::gpar(fontsize = 9),
                                labels_gp = grid::gpar(fontsize = 8)))
  open_tiff_base(filename, width = 8, height = 7)
  draw(h, heatmap_legend_side = "right", annotation_legend_side = "bottom")
  dev.off()
  invisible(h)
}

#' Hierarchical clustering dendrogram of samples
dendrogram_plot <- function(mat, sample_groups, filename, tag) {
  d <- stats::dist(t(mat))
  hc <- stats::hclust(d, method = "average")
  lab_col <- ifelse(sample_groups[colnames(mat)] == "Tumor", COL_TUMOR, COL_NORMAL)
  open_tiff_base(filename, width = 9, height = 6)
  graphics::plot(hc, main = sprintf("%s sample clustering (average linkage)", tag),
                 xlab = "", sub = "", labels = FALSE)
  stats::rect.hclust(hc, k = 2)
  graphics::legend("topright", legend = c("Tumor", "Normal"),
                   col = c(COL_TUMOR, COL_NORMAL), pch = 15, cex = 0.8, bty = "n")
  dev.off()
  invisible(hc)
}

#' Probe-to-gene retention bar plot
retention_plot <- function(counts, filename, tag) {
  df <- data.frame(stage = names(counts), n = as.numeric(counts),
                   stringsAsFactors = FALSE)
  df$stage <- factor(df$stage, levels = names(counts))
  p <- ggplot(df, aes(x = stage, y = n)) +
    geom_col(fill = "#0072B2", colour = "black", width = 0.7) +
    geom_text(aes(label = n), vjust = -0.3, size = 3.2) +
    theme_bc(base_size = 10) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
    labs(x = "Processing stage", y = "Features",
         title = sprintf("%s probe-to-gene retention", tag))
  save_tiff_600dpi(p, filename, width = 8, height = 6)
  p
}

#' Outlier assessment plot
outlier_plot <- function(out_df, filename, tag) {
  p <- ggplot(out_df, aes(x = sample_id, y = mean_pairwise_correlation,
                          colour = is_outlier)) +
    geom_point(size = 1) +
    geom_hline(yintercept = median(out_df$mean_pairwise_correlation, na.rm = TRUE),
               linetype = 2) +
    scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2")) +
    theme_bc(base_size = 9) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(x = "Sample", y = "Mean pairwise Pearson correlation",
         title = sprintf("%s sample outlier assessment", tag), colour = "Outlier")
  save_tiff_600dpi(p, filename, width = 9, height = 5)
  p
}
