# ============================================================
# Step 14: Independent Validation in GSE70947
# ============================================================
set.seed(123)
library(limma)
library(ggplot2)
library(pheatmap)
library(pROC)
library(caret)
library(PRROC)
library(randomForest)

base_dir <- normalizePath(file.path(".", "Breast Cancer 5 Necroptosis"))
data_dir <- normalizePath(file.path(".", "Datasets", "data"))
output_dir <- base_dir

sink(file.path(output_dir, "logs", "14_external_validation.log"), split=TRUE)
cat("INDEPENDENT VALIDATION (GSE70947)\n")
cat("==================================\n\n")

save_tiff_600dpi <- function(plot, filename, width = 7, height = 6) {
  ggsave(filename=filename, plot=plot, width=width, height=height,
         units="in", dpi=600, compression="lzw", bg="white", device="tiff")
}

# Load validation data
val_csv <- file.path(data_dir, "processed", "GSE70947_expression.csv")
val_meta_file <- file.path(data_dir, "metadata", "GSE70947_sample_classes.csv")
val_expr <- read.csv(val_csv, row.names=1, check.names=FALSE)
val_meta <- read.csv(val_meta_file, stringsAsFactors=FALSE, check.names=FALSE)

# Check for paired samples
if("pair" %in% colnames(val_meta) || "patient_id" %in% colnames(val_meta)) {
  cat("Paired samples detected\n")
  paired_available <- TRUE
} else {
  paired_available <- FALSE
  cat("No paired sample information - using unpaired analysis\n")
}

# Common samples between expression and metadata
common_val <- intersect(colnames(val_expr), val_meta[["sample_id"]])
val_expr <- val_expr[, common_val, drop=FALSE]
val_meta <- val_meta[match(common_val, val_meta[["sample_id"]]), , drop=FALSE]
rownames(val_meta) <- val_meta[["sample_id"]]

# Normalize validation data independently
val_mat <- as.matrix(val_expr)
max_val <- max(val_mat, na.rm=TRUE)
if(max_val > 30) {
  val_mat <- log2(val_mat + 1)
  cat("Applied log2(x+1) transformation to validation data\n")
}
val_norm <- normalizeBetweenArrays(val_mat, method="quantile")

# Load discovery data references
combat_expr <- readRDS(file.path(output_dir, "data_processed", "discovery_combat_corrected.rds"))
merged_meta <- readRDS(file.path(output_dir, "data_processed", "discovery_metadata.rds"))
de_nrgs <- readLines(file.path(output_dir, "results", "necroptosis_genes", "DE_NRGs.txt"))
final_genes <- readLines(file.path(output_dir, "results", "consensus_genes", "final_consensus_genes.txt"))
ng_official <- readRDS(file.path(output_dir, "metadata", "necroptosis_genes_official.rds"))

# Common genes between validation and reference
common_genes_all <- intersect(rownames(val_norm), rownames(combat_expr))
cat("Common genes (validation vs discovery):", length(common_genes_all), "\n")

# ============ Validate DE-NRG Expression ============
val_de_nrgs <- intersect(de_nrgs, rownames(val_norm))
cat("DE-NRGs found in validation:", length(val_de_nrgs), "\n")

if(length(val_de_nrgs) > 0) {
  val_de_nrg_expr <- val_norm[val_de_nrgs, , drop=FALSE]
  
  # Load discovery logFC for comparison
  deg_all <- read.csv(file.path(output_dir, "results", "tables", "DEG_all.csv"), row.names=1)
  
  # Compute validation fold changes
  val_tumor <- rowMeans(val_de_nrg_expr[, val_meta[["class"]] == "Tumor", drop=FALSE])
  val_normal <- rowMeans(val_de_nrg_expr[, val_meta[["class"]] == "Normal", drop=FALSE])
  val_lfc <- val_tumor - val_normal
  
  # Discovery LFC
  disc_lfc <- deg_all[val_de_nrgs, "logFC"]
  
  # Comparison
  lfc_comp <- data.frame(
    Gene = val_de_nrgs,
    Discovery_logFC = disc_lfc,
    Validation_logFC = val_lfc,
    Direction_agreement = sign(disc_lfc) == sign(val_lfc)
  )
  lfc_comp <- lfc_comp[order(abs(lfc_comp$Validation_logFC), decreasing=TRUE), ]
  write.csv(lfc_comp, file.path(output_dir, "results", "external_validation", "DE_NRG_foldchange_comparison.csv"), row.names=FALSE)
  
  cat("Fold change direction agreement:", mean(lfc_comp$Direction_agreement)*100, "%\n")
  
  # Discovery vs Validation FC plot
  p_fc_comp <- ggplot(lfc_comp, aes(x=Discovery_logFC, y=Validation_logFC, color=Direction_agreement)) +
    geom_point(size=2, alpha=0.7) + theme_bw() +
    geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed") +
    geom_smooth(method="lm", se=TRUE, color="black", linetype="dotted") +
    scale_color_manual(values=c("TRUE"="forestgreen", "FALSE"="red")) +
    labs(title="Discovery vs Validation: Fold Change Comparison", x="Discovery log2FC", y="Validation log2FC")
  save_tiff_600dpi(p_fc_comp, file.path(output_dir, "results", "external_validation", "FC_comparison.tiff"), width=8, height=7)
  
  # Validation DE-NRG boxplots
  if(length(val_de_nrgs) <= 9) {
    plot_list <- list()
    for(gene in val_de_nrgs) {
      df <- data.frame(
        Expression = val_de_nrg_expr[gene, ],
        Group = val_meta[["class"]]
      )
      p <- ggplot(df, aes(x=Group, y=Expression, fill=Group)) +
        geom_boxplot() + theme_bw() +
        scale_fill_manual(values=c("Normal"="blue", "Tumor"="red")) +
        labs(title=gene, x="", y="Expression") +
        theme(legend.position="none", plot.title=element_text(hjust=0.5))
      plot_list[[gene]] <- p
    }
    if(length(plot_list) > 0) {
      library(gridExtra)
      p_grid <- grid.arrange(grobs=plot_list, ncol=3, top="Validation: DE-NRG Expression")
      ggsave(file.path(output_dir, "results", "external_validation", "validation_DE_NRG_boxplots.tiff"),
             plot=p_grid, width=10, height=3*ceiling(length(plot_list)/3), dpi=600, compression="lzw", bg="white")
    }
  }
  
  # Validation heatmap
  val_annot <- data.frame(Group=val_meta[["class"]], row.names=colnames(val_norm))
  val_colors <- list(Group=c(Normal="blue", Tumor="red"))
  
  tiff(file.path(output_dir, "results", "external_validation", "validation_DE_NRG_heatmap.tiff"),
       width=10, height=max(6, length(val_de_nrgs)*0.3), units="in", res=600, compression="lzw", bg="white")
  pheatmap(val_de_nrg_expr, annotation_col=val_annot, annotation_colors=val_colors,
           scale="row", main="Validation: DE-NRG Expression",
           show_rownames=TRUE, show_colnames=FALSE,
           color=colorRampPalette(c("navy", "white", "firebrick3"))(100),
           fontsize_row=8, cutree_cols=2)
  dev.off()
}

# ============ Validate Necroptosis Score ============
val_ng <- intersect(ng_official, rownames(val_norm))
cat("Necroptosis genes in validation:", length(val_ng), "\n")

if(length(val_ng) > 0) {
  library(GSVA)
  gsva_params <- gsvaParam(as.matrix(val_norm), list(necroptosis=val_ng), maxDiff=TRUE, kcdf="Gaussian")
  val_score <- gsva(gsva_params)
  val_gsva_score <- as.numeric(val_score[1, ])
  names(val_gsva_score) <- colnames(val_norm)
  
  val_score_df <- data.frame(
    Score = val_gsva_score,
    Group = val_meta[["class"]]
  )
  
  p_score_val <- ggplot(val_score_df, aes(x=Group, y=Score, fill=Group)) +
    geom_boxplot() + theme_bw() +
    scale_fill_manual(values=c("Normal"="blue", "Tumor"="red")) +
    labs(title="Validation: Necroptosis Score", x="", y="GSVA Score")
  save_tiff_600dpi(p_score_val, file.path(output_dir, "results", "external_validation", "validation_score_boxplot.tiff"), width=6, height=5)
  
  # Compare score distributions
  grp <- val_meta[["class"]]
  cat("Validation groups:", length(grp), "samples:", sum(grp=="Normal"), "Normal,", sum(grp=="Tumor"), "Tumor\n")
  wt <- tryCatch(wilcox.test(val_gsva_score ~ factor(grp)), error=function(e) NULL)
  if(!is.null(wt)) {
    cat("Validation score p-value:", wt$p.value, "\n")
  } else {
    cat("Validation score comparison skipped\n")
  }
  
  # Discovery vs Validation score comparison
  disc_score <- readRDS(file.path(output_dir, "results", "necroptosis_score", "necroptosis_score.rds"))
  disc_score_df <- data.frame(Score=disc_score, Group="Discovery", Class=merged_meta$class)
  val_score_df2 <- data.frame(Score=val_gsva_score, Group="Validation", Class=val_meta[["class"]])
  comp_score <- rbind(disc_score_df, val_score_df2)
  
  p_score_comp <- ggplot(comp_score, aes(x=Class, y=Score, fill=interaction(Group, Class))) +
    geom_boxplot() + theme_bw() +
    labs(title="Necroptosis Score: Discovery vs Validation", x="", y="GSVA Score") +
    scale_fill_manual(values=c("Discovery.Normal"="steelblue", "Discovery.Tumor"="salmon",
                                "Validation.Normal"="lightblue", "Validation.Tumor"="pink"))
  save_tiff_600dpi(p_score_comp, file.path(output_dir, "results", "external_validation", "score_discovery_vs_validation.tiff"), width=8, height=6)
}

# ============ Validate Diagnostic Model ============
val_final_genes <- intersect(final_genes, rownames(val_norm))
cat("Final consensus genes in validation:", length(val_final_genes), "\n")

if(length(val_final_genes) >= 3) {
  # Apply locked model (from RF which was best in discovery)
  rf_model_path <- file.path(output_dir, "model_objects", "rf_final_model.rds")
  if(file.exists(rf_model_path)) {
    rf_final <- readRDS(rf_model_path)
    
    # Check feature overlap
    rf_features <- rownames(rf_final$importance)
    common_features <- intersect(rf_features, rownames(val_norm))
    cat("Common RF features in validation:", length(common_features), "/", length(rf_features), "\n")
    
    if(length(common_features) >= 3) {
      val_X <- t(val_norm[common_features, , drop=FALSE])
      
      # Impute NA values with row medians
      for(j in seq_len(ncol(val_X))) {
        if(any(is.na(val_X[, j]))) {
          val_X[is.na(val_X[, j]), j] <- median(val_X[, j], na.rm=TRUE)
        }
      }
      
      # Need to ensure same columns as training
      missing_cols <- setdiff(rf_features, common_features)
      if(length(missing_cols) > 0) {
        cat("Missing features:", length(missing_cols), "- filling with 0\n")
        for(mc in missing_cols) {
          val_X <- cbind(val_X, rep(0, nrow(val_X)))
          colnames(val_X)[ncol(val_X)] <- mc
        }
      }
      val_X <- val_X[, rf_features, drop=FALSE]
      
      val_pred_prob <- predict(rf_final, val_X, type="prob")[, "Tumor"]
      val_pred_class <- predict(rf_final, val_X)
      
      y_val <- factor(val_meta[["class"]], levels=c("Normal", "Tumor"))
      
      # Validation ROC
      roc_val <- roc(y_val, val_pred_prob)
      cat("Validation AUC:", auc(roc_val), "\n")
      
      # ROC plot
      roc_val_df <- data.frame(FPR=1-roc_val$specificities, TPR=roc_val$sensitivities)
      p_roc_val <- ggplot(roc_val_df, aes(x=FPR, y=TPR)) +
        geom_line(color="darkgreen", size=1) + theme_bw() +
        geom_abline(intercept=0, slope=1, linetype="dashed") +
        labs(title=paste("Validation ROC (AUC =", round(auc(roc_val), 3), ")"),
             x="False Positive Rate", y="True Positive Rate")
      save_tiff_600dpi(p_roc_val, file.path(output_dir, "results", "external_validation", "validation_ROC.tiff"), width=7, height=6)
      
      # PR curve
      pr_val <- pr.curve(scores.class0=val_pred_prob[y_val=="Tumor"],
                         scores.class1=val_pred_prob[y_val=="Normal"], curve=TRUE)
      
      # Confusion matrix
      cm_val <- confusionMatrix(val_pred_class, y_val)
      cat("Validation Confusion Matrix:\n")
      print(cm_val$table)
      cat("Accuracy:", cm_val$overall["Accuracy"], "\n")
      write.csv(as.matrix(cm_val$table), file.path(output_dir, "results", "external_validation", "validation_confusion_matrix.csv"))
      
      # Calibration
      y_val_bin <- ifelse(y_val == "Tumor", 1, 0)
      cal_val_df <- data.frame(pred=val_pred_prob, obs=y_val_bin)
      
      p_cal_val <- ggplot(cal_val_df, aes(x=pred, y=obs)) +
        geom_point(alpha=0.3) + geom_smooth(method="loess", se=TRUE, color="darkgreen") +
        geom_abline(intercept=0, slope=1, linetype="dashed") +
        theme_bw() + labs(title="Validation Calibration", x="Predicted Probability", y="Observed")
      save_tiff_600dpi(p_cal_val, file.path(output_dir, "results", "external_validation", "validation_calibration.tiff"), width=7, height=6)
      
      # Discovery vs Validation performance
      perf_comp <- data.frame(
        Metric = c("AUC", "Accuracy", "Sensitivity", "Specificity"),
        Discovery = c(NA, NA, NA, NA),
        Validation = c(auc(roc_val), cm_val$overall["Accuracy"], 
                       cm_val$byClass["Sensitivity"], cm_val$byClass["Specificity"])
      )
      disc_metrics <- read.csv(file.path(output_dir, "results", "predictive_model", "model_comparison.csv"))
      if("RF" %in% disc_metrics$Model) {
        rf_metrics <- disc_metrics[disc_metrics$Model == "RF", ]
        perf_comp$Discovery <- c(rf_metrics$AUC, rf_metrics$Accuracy,
                                  rf_metrics$Sensitivity, rf_metrics$Specificity)
      }
      write.csv(perf_comp, file.path(output_dir, "results", "external_validation", "discovery_vs_validation_performance.csv"), row.names=FALSE)
      
      p_perf_comp <- ggplot(perf_comp, aes(x=Metric)) +
        geom_point(aes(y=Discovery, color="Discovery"), size=3) +
        geom_point(aes(y=Validation, color="Validation"), size=3) +
        theme_bw() + labs(title="Discovery vs Validation Performance", x="", y="Value") +
        scale_color_manual(values=c("Discovery"="blue", "Validation"="red"))
      save_tiff_600dpi(p_perf_comp, file.path(output_dir, "results", "external_validation", "performance_comparison.tiff"), width=8, height=6)
    }
  }
}

# ============ Validate Immune Infiltration ============
if(require(xCell, quietly=TRUE)) {
  cat("Running xCell on validation data...\n")
  val_xcell <- xCellAnalysis(as.matrix(val_norm), rnaseq=FALSE, parallel.sz=1)
  val_xcell_filt <- val_xcell[rowSums(val_xcell) > 0, ]
  
  if(nrow(val_xcell_filt) > 0) {
    # Tumor vs normal comparison
    val_imm_diff <- data.frame()
    for(ct in rownames(val_xcell_filt)) {
      t_scores <- as.numeric(val_xcell_filt[ct, val_meta[["class"]] == "Tumor"])
      n_scores <- as.numeric(val_xcell_filt[ct, val_meta[["class"]] == "Normal"])
      wt <- tryCatch(wilcox.test(t_scores, n_scores), error=function(e) NULL)
      if(!is.null(wt)) {
        val_imm_diff <- rbind(val_imm_diff, data.frame(
          CellType = ct,
          Tumor_Mean = mean(t_scores),
          Normal_Mean = mean(n_scores),
          p_value = wt$p.value
        ))
      }
    }
    val_imm_diff$FDR <- p.adjust(val_imm_diff$p_value, "BH")
    write.csv(val_imm_diff, file.path(output_dir, "results", "external_validation", "validation_immune_differential.csv"), row.names=FALSE)
    
    # Validation immune heatmap
    tiff(file.path(output_dir, "results", "external_validation", "validation_immune_heatmap.tiff"),
         width=12, height=6, units="in", res=600, compression="lzw", bg="white")
    val_annot <- data.frame(Group=val_meta[["class"]], row.names=colnames(val_norm))
    pheatmap(val_xcell_filt, annotation_col=val_annot,
             scale="row", main="Validation: Immune Infiltration",
             show_rownames=TRUE, show_colnames=FALSE, fontsize_row=8,
             color=colorRampPalette(c("navy", "white", "firebrick3"))(100))
    dev.off()
    
    # Gene-immune correlation consistency
    disc_immune_file <- file.path(output_dir, "results", "immune_infiltration", "DE_NRG_immune_correlation.csv")
    if(file.exists(disc_immune_file) && length(val_de_nrgs) > 0) {
      disc_gene_imm <- as.matrix(read.csv(disc_immune_file, row.names=1))
      val_common_genes <- intersect(rownames(disc_gene_imm), rownames(val_norm))
      val_common_imm <- intersect(colnames(disc_gene_imm), rownames(val_xcell_filt))
      
      if(length(val_common_genes) > 0 && length(val_common_imm) > 0) {
        val_gene_imm_cor <- matrix(NA, length(val_common_genes), length(val_common_imm))
        rownames(val_gene_imm_cor) <- val_common_genes
        colnames(val_gene_imm_cor) <- val_common_imm
        
        for(i in seq_along(val_common_genes)) {
          for(j in seq_along(val_common_imm)) {
            ct <- tryCatch(cor.test(val_norm[val_common_genes[i], ],
                                     as.numeric(val_xcell_filt[val_common_imm[j], ]),
                                     method="spearman"),
                           error=function(e) NULL)
            if(!is.null(ct)) val_gene_imm_cor[i,j] <- ct$estimate
          }
        }
        
        write.csv(val_gene_imm_cor, file.path(output_dir, "results", "external_validation", "validation_gene_immune_correlation.csv"))
        
        # Consistency with discovery
        disc_sub <- disc_gene_imm[val_common_genes, val_common_imm, drop=FALSE]
        val_sub <- val_gene_imm_cor
        
        if(all(dim(disc_sub) == dim(val_sub))) {
          cor_consistency <- cor(as.vector(disc_sub), as.vector(val_sub), use="pairwise.complete.obs", method="spearman")
          cat("Gene-immune correlation consistency (Spearman r):", round(cor_consistency, 3), "\n")
        }
      }
    }
  }
}

cat("\n\nEXTERNAL VALIDATION COMPLETE\n")
sink()
