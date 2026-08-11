# =============================================================================
# Script 02: Metadata Curation and Independent Preprocessing
# Breast Cancer 7 Cell Death Crosstalk
# =============================================================================

source(file.path("..", "config", "analysis_config.R"))

library(tidyverse)
library(data.table)
library(limma)
library(Biobase)
library(ComplexHeatmap)
library(circlize)

log_message("Starting metadata curation and preprocessing.")

# --- 1. LOAD SAMPLE CLASS FILES ---
load_sample_classes <- function(dataset) {
  f <- SAMPLE_CLASS_FILES[[dataset]]
  df <- read.csv(f, stringsAsFactors = FALSE)
  df$dataset <- dataset
  return(df)
}

all_classes <- lapply(names(SAMPLE_CLASS_FILES), load_sample_classes)
names(all_classes) <- names(SAMPLE_CLASS_FILES)

# Save sample composition
sample_composition <- bind_rows(all_classes)
write.csv(sample_composition, file.path(OUTPUT$preprocessing, "all_sample_classes.csv"), row.names = FALSE)

# Summarise
for (ds in names(all_classes)) {
  tab <- table(all_classes[[ds]]$class)
  log_message(sprintf("%s: %s", ds, paste(names(tab), tab, sep = "=", collapse = ", ")))
}

# --- 2. LOAD EXPRESSION DATA ---
load_expression <- function(dataset) {
  f <- EXPRESSION_FILES[[dataset]]
  log_message(sprintf("Loading %s expression: %s", dataset, f))
  mat <- fread(f, data.table = FALSE)
  # First column is gene/probe identifier
  rownames(mat) <- mat[[1]]
  mat <- mat[, -1]
  return(as.matrix(mat))
}

expr_list <- lapply(names(EXPRESSION_FILES), load_expression)
names(expr_list) <- names(EXPRESSION_FILES)

# --- 3. ASSESS EXPRESSION SCALE ---
for (ds in names(expr_list)) {
  mat <- expr_list[[ds]]
  rng <- range(mat, na.rm = TRUE)
  med <- median(mat, na.rm = TRUE)
  log_message(sprintf("%s expression: range=[%.2f, %.2f], median=%.2f, log2=%s",
                      ds, rng[1], rng[2], med,
                      ifelse(med > 20, "likely not log2", ifelse(med < 5, "likely log2", "ambiguous"))))
}

# GSE15852 expression seems log2-scale (median ~8)
# GSE42568 seems log2-scale as well
# GSE29431 and GSE70947 likely similar
# All appear to be already log2 transformed

# --- 4. PROCESS EACH DATASET INDIVIDUALLY ---
process_dataset <- function(dataset) {
  log_message(sprintf("Processing %s individually...", dataset))

  mat <- expr_list[[dataset]]
  classes <- all_classes[[dataset]]
  
  # Match samples
  common_samples <- intersect(colnames(mat), classes$sample_id)
  if (length(common_samples) == 0) {
    log_message(sprintf("WARNING: %s - No matching samples between expression and metadata!", dataset), "ERROR")
    return(NULL)
  }
  
  mat <- mat[, common_samples, drop = FALSE]
  classes <- classes[match(common_samples, classes$sample_id), ]
  
  # Ensure proper ordering
  group <- factor(classes$class, levels = c("Normal", "Tumor"))
  
  # Remove constant/zero-variance columns
  sample_var <- apply(mat, 2, var, na.rm = TRUE)
  valid_samples <- which(sample_var > 1e-10)
  if (length(valid_samples) < 4) {
    log_message(sprintf("WARNING: %s - too few variable samples (%d)", dataset, length(valid_samples)), "WARNING")
    valid_samples <- seq_len(ncol(mat))
  }
  mat_pca <- mat[, valid_samples, drop = FALSE]
  
  # PCA before normalisation (no scaling to avoid constant-column issues)
  pca <- tryCatch(prcomp(t(mat_pca), scale = FALSE, center = TRUE), error = function(e) {
    log_message(sprintf("PCA failed for %s: %s", dataset, e$message), "WARNING")
    return(NULL)
  })
  if (is.null(pca)) return(NULL)
  pca_var <- summary(pca)$importance[2, 1:2] * 100
  
  pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Sample = colnames(mat_pca),
    Group = group[valid_samples],
    Dataset = dataset,
    stringsAsFactors = FALSE
  )
  
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 2.5, alpha = 0.8) +
    stat_ellipse(level = 0.95, show.legend = FALSE) +
    scale_color_manual(values = c("Normal" = "#0072B2", "Tumor" = "#D55E00")) +
    labs(title = sprintf("%s PCA (Before Normalisation)", dataset),
         x = sprintf("PC1 (%.1f%%)", pca_var[1]),
         y = sprintf("PC2 (%.1f%%)", pca_var[2])) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  save_tiff_600dpi(p_pca, file.path(OUTPUT$preprocessing, sprintf("%s_PCA_before.tiff", dataset)))
  
  # Expression density plot
  expr_long <- as.data.frame(mat) %>% 
    pivot_longer(everything(), names_to = "Sample", values_to = "Expression")
  expr_long$Group <- group[expr_long$Sample]
  
  p_density <- ggplot(expr_long, aes(x = Expression, color = Group)) +
    geom_density(alpha = 0.5) +
    scale_color_manual(values = c("Normal" = "#0072B2", "Tumor" = "#D55E00")) +
    labs(title = sprintf("%s Expression Density", dataset), x = "Expression Level", y = "Density") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  save_tiff_600dpi(p_density, file.path(OUTPUT$preprocessing, sprintf("%s_density.tiff", dataset)))
  
  # Boxplot
  p_box <- ggplot(expr_long, aes(x = Sample, y = Expression, fill = Group)) +
    geom_boxplot(outlier.size = 0.3) +
    scale_fill_manual(values = c("Normal" = "#0072B2", "Tumor" = "#D55E00")) +
    labs(title = sprintf("%s Expression Boxplot", dataset), x = "Sample", y = "Expression") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
  save_tiff_600dpi(p_box, file.path(OUTPUT$preprocessing, sprintf("%s_boxplot.tiff", dataset)))
  
  # Sample correlation heatmap
  cor_mat <- cor(mat, method = "spearman")
  
  ann_df <- data.frame(Group = group, row.names = colnames(mat))
  ann_colors <- list(Group = c(Normal = "#0072B2", Tumor = "#D55E00"))
  
  ht <- Heatmap(cor_mat,
                name = "Spearman r",
                col = colorRamp2(c(0.5, 0.75, 1), c("#2166AC", "#F7F7F7", "#B2182B")),
                top_annotation = HeatmapAnnotation(
                  Group = group,
                  col = list(Group = c(Normal = "#0072B2", Tumor = "#D55E00"))
                ),
                show_row_names = FALSE,
                show_column_names = FALSE,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                column_title = sprintf("%s Sample Correlation", dataset))
  
  open_tiff(file.path(OUTPUT$preprocessing, sprintf("%s_correlation_heatmap.tiff", dataset)))
  draw(ht)
  dev.off()
  
  # Clustering dendrogram
  hc <- hclust(dist(t(mat)), method = "average")
  hc_dend <- as.dendrogram(hc)
  
  open_tiff(file.path(OUTPUT$preprocessing, sprintf("%s_dendrogram.tiff", dataset)))
  par(mar = c(4, 4, 2, 8))
  plot(hc, main = sprintf("%s Sample Dendrogram", dataset), 
       xlab = "Sample", sub = "", cex = 0.5)
  dev.off()
  
  # Missing value assessment
  na_count <- sum(is.na(mat))
  na_per_gene <- rowSums(is.na(mat))
  na_per_sample <- colSums(is.na(mat))
  
  log_message(sprintf("%s: NA values = %d (%.4f%% of matrix)", 
                      dataset, na_count, 100 * na_count / prod(dim(mat))))
  
  if (na_count > 0) {
    # Missing value heatmap
    na_mat <- is.na(mat)
    ht_na <- Heatmap(na_mat * 1,
                     name = "Missing",
                     col = c("0" = "white", "1" = "red"),
                     show_row_names = FALSE,
                     show_column_names = FALSE,
                     column_title = sprintf("%s Missing Values", dataset))
    
    open_tiff(file.path(OUTPUT$preprocessing, sprintf("%s_missing_values.tiff", dataset)))
    draw(ht_na)
    dev.off()
  }
  
  # Store results
  result <- list(
    expression = mat,
    metadata = classes,
    group = group,
    pca = pca,
    pca_df = pca_df,
    correlation = cor_mat,
    hclust = hc,
    missing = list(total_na = na_count, per_gene = na_per_gene, per_sample = na_per_sample)
  )
  
  saveRDS(result, file.path(OUTPUT$preprocessing, sprintf("%s_processed.rds", dataset)))
  
  log_message(sprintf("%s processing complete: %d genes x %d samples", 
                      dataset, nrow(mat), ncol(mat)))
  
  return(result)
}

processed_list <- lapply(names(expr_list), process_dataset)
names(processed_list) <- names(expr_list)

# Save all processed objects
saveRDS(processed_list, file.path(OUTPUT$preprocessing, "all_processed_datasets.rds"))

# --- 5. CREATE DISCOVERY COHORT METADATA SUMMARY ---
discovery_meta <- bind_rows(all_classes[names(all_classes) %in% DISCOVERY_DATASETS])
write.csv(discovery_meta, file.path(OUTPUT$preprocessing, "discovery_metadata_summary.csv"), row.names = FALSE)

validation_meta <- all_classes[[VALIDATION_DATASET]]
write.csv(validation_meta, file.path(OUTPUT$preprocessing, "validation_metadata_summary.csv"), row.names = FALSE)

# --- 6. PLATFORM ANNOTATION CHECK ---
log_message("Checking platform compatibility...")
# GSE42568 and GSE15852: GPL570 (Affymetrix Human Genome U133 Plus 2.0)
# GSE29431: GPL96 (Affymetrix Human Genome U133A)
# GSE70947: GPL13607 (Agilent-028004 sureprint G3 Human GE 8x60k)
# All are different platforms, so probe-level merging is not possible
# All processed to gene-level already, so we can merge on gene symbols

# Check gene overlap among discovery datasets
gse42568_genes <- rownames(processed_list$GSE42568$expression)
gse15852_genes <- rownames(processed_list$GSE15852$expression)
gse29431_genes <- rownames(processed_list$GSE29431$expression)

common_discovery <- Reduce(intersect, list(gse42568_genes, gse15852_genes, gse29431_genes))
log_message(sprintf("Common genes among discovery datasets: %d", length(common_discovery)))
write.csv(data.frame(Gene = common_discovery), 
          file.path(OUTPUT$preprocessing, "discovery_common_genes.csv"), row.names = FALSE)

log_message("Metadata curation and preprocessing complete.")
