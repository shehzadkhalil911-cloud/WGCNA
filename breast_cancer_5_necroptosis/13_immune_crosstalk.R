# ============================================================
# Step 13: Necroptosis–Immune Crosstalk
# ============================================================
set.seed(123)
library(ggplot2)
library(pheatmap)
library(circlize)

base_dir <- normalizePath(file.path(".", "Breast Cancer 5 Necroptosis"))
output_dir <- base_dir

sink(file.path(output_dir, "logs", "13_immune_crosstalk.log"), split=TRUE)
cat("NECROPTOSIS-IMMUNE CROSSTALK\n")
cat("============================\n\n")

save_tiff_600dpi <- function(plot, filename, width = 7, height = 6) {
  ggsave(filename=filename, plot=plot, width=width, height=height,
         units="in", dpi=600, compression="lzw", bg="white", device="tiff")
}

combat_expr <- readRDS(file.path(output_dir, "data_processed", "discovery_combat_corrected.rds"))
merged_meta <- readRDS(file.path(output_dir, "data_processed", "discovery_metadata.rds"))
de_nrgs <- readLines(file.path(output_dir, "results", "necroptosis_genes", "DE_NRGs.txt"))
final_genes <- readLines(file.path(output_dir, "results", "consensus_genes", "final_consensus_genes.txt"))
necroptosis_score <- readRDS(file.path(output_dir, "results", "necroptosis_score", "necroptosis_score.rds"))
xcell_results <- readRDS(file.path(output_dir, "results", "immune_infiltration", "xcell_results.rds"))

# Filter immune cells
xcell_filt <- xcell_results[rowSums(xcell_results) > 0, ]
immune_cells <- rownames(xcell_filt)

# Align samples
common_samples <- intersect(colnames(combat_expr), names(necroptosis_score))
common_samples <- intersect(common_samples, colnames(xcell_filt))
cat("Common samples:", length(common_samples), "\n")

necroptosis_score <- necroptosis_score[common_samples]
xcell_filt <- xcell_filt[, common_samples, drop=FALSE]
combat_expr <- combat_expr[, common_samples, drop=FALSE]
merged_meta <- merged_meta[common_samples, , drop=FALSE]

# ============ DE-NRG vs Immune Cells ============
de_nrg_common <- intersect(de_nrgs, rownames(combat_expr))
cat("DE-NRGs for immune correlation:", length(de_nrg_common), "\n")

if(length(de_nrg_common) > 0 && length(immune_cells) > 0) {
  gene_immune_cor <- matrix(NA, length(de_nrg_common), length(immune_cells))
  rownames(gene_immune_cor) <- de_nrg_common
  colnames(gene_immune_cor) <- immune_cells
  gene_immune_pval <- gene_immune_cor
  
  for(i in seq_along(de_nrg_common)) {
    for(j in seq_along(immune_cells)) {
      ct <- tryCatch(cor.test(combat_expr[de_nrg_common[i], ], 
                               as.numeric(xcell_filt[j, ]), 
                               method="spearman"),
                     error=function(e) NULL)
      if(!is.null(ct)) {
        gene_immune_cor[i,j] <- ct$estimate
        gene_immune_pval[i,j] <- ct$p.value
      }
    }
  }
  
  # Multiple testing correction
  gene_immune_fdr <- matrix(p.adjust(as.vector(gene_immune_pval), "BH"),
                            nrow=nrow(gene_immune_pval))
  rownames(gene_immune_fdr) <- rownames(gene_immune_pval)
  colnames(gene_immune_fdr) <- colnames(gene_immune_pval)
  
  # Save significant associations
  sig_pairs <- which(abs(gene_immune_cor) > 0.3 & gene_immune_fdr < 0.05, arr.ind=TRUE)
  cat("DE-NRG-Immune significant associations:", nrow(sig_pairs), "\n")
  
  # Gene-immune heatmap
  tiff(file.path(output_dir, "results", "immune_infiltration", "DE_NRG_immune_heatmap.tiff"),
       width=10, height=max(6, nrow(gene_immune_cor)*0.3), units="in", res=600, compression="lzw", bg="white")
  pheatmap(gene_immune_cor, 
           main="DE-NRG vs Immune Cell Correlation",
           color=colorRampPalette(c("navy", "white", "firebrick3"))(100),
           display_numbers=ifelse(gene_immune_fdr < 0.05, "*", ""),
           fontsize_number=6, fontsize_row=8, fontsize_col=8,
           cluster_rows=TRUE, cluster_cols=TRUE)
  dev.off()
  
  write.csv(gene_immune_cor, file.path(output_dir, "results", "immune_infiltration", "DE_NRG_immune_correlation.csv"))
  write.csv(gene_immune_fdr, file.path(output_dir, "results", "immune_infiltration", "DE_NRG_immune_FDR.csv"))
}

# ============ Consensus Genes vs Immune Cells ============
consensus_common <- intersect(final_genes, rownames(combat_expr))
if(length(consensus_common) > 0 && length(immune_cells) > 0) {
  cons_immune_cor <- matrix(NA, length(consensus_common), length(immune_cells))
  rownames(cons_immune_cor) <- consensus_common
  colnames(cons_immune_cor) <- immune_cells
  cons_immune_pval <- cons_immune_cor
  
  for(i in seq_along(consensus_common)) {
    for(j in seq_along(immune_cells)) {
      ct <- tryCatch(cor.test(combat_expr[consensus_common[i], ], 
                               as.numeric(xcell_filt[j, ]), 
                               method="spearman"),
                     error=function(e) NULL)
      if(!is.null(ct)) {
        cons_immune_cor[i,j] <- ct$estimate
        cons_immune_pval[i,j] <- ct$p.value
      }
    }
  }
  
  cons_immune_fdr <- matrix(p.adjust(as.vector(cons_immune_pval), "BH"),
                            nrow=nrow(cons_immune_pval))
  rownames(cons_immune_fdr) <- rownames(cons_immune_pval)
  colnames(cons_immune_fdr) <- colnames(cons_immune_pval)
  
  write.csv(cons_immune_cor, file.path(output_dir, "results", "immune_infiltration", "consensus_immune_correlation.csv"))
}

# ============ Necroptosis Score vs Immune Cells ============
score_list <- lapply(immune_cells, function(cell_id) {
  ct_res <- tryCatch(
    cor.test(necroptosis_score, as.numeric(xcell_filt[cell_id, ]), method="spearman"),
    error=function(e) list(estimate=NA_real_, p.value=1)
  )
  data.frame(CellType=cell_id, Correlation=as.numeric(ct_res$estimate), Pvalue=as.numeric(ct_res$p.value),
             stringsAsFactors=FALSE)
})
score_immune_df <- do.call(rbind, score_list)
score_immune_df$FDR <- p.adjust(score_immune_df$Pvalue, "BH")
score_immune_df <- score_immune_df[order(abs(score_immune_df$Correlation), decreasing=TRUE), ]
write.csv(score_immune_df, file.path(output_dir, "results", "immune_infiltration", "score_immune_correlation.csv"), row.names=FALSE)

# Score-immune correlation plot
p_score_imm <- ggplot(score_immune_df, aes(x=reorder(CellType, Correlation), y=Correlation, fill=FDR<0.05)) +
  geom_bar(stat="identity") + theme_bw() +
  coord_flip() + labs(title="Necroptosis Score vs Immune Cells", x="", y="Spearman Correlation") +
  scale_fill_manual(values=c("TRUE"="red", "FALSE"="grey"), name="Significant (FDR<0.05)")
save_tiff_600dpi(p_score_imm, file.path(output_dir, "results", "immune_infiltration", "score_immune_correlation_plot.tiff"), width=8, height=6)

# ============ Chord Diagram ============
if(exists("gene_immune_cor") && nrow(sig_pairs) > 0) {
  top_pairs <- sig_pairs[order(abs(gene_immune_cor[sig_pairs]), decreasing=TRUE), , drop=FALSE]
  top_pairs <- top_pairs[1:min(50, nrow(top_pairs)), , drop=FALSE]
  
  chord_data <- data.frame(
    from = rownames(gene_immune_cor)[top_pairs[,1]],
    to = colnames(gene_immune_cor)[top_pairs[,2]],
    value = gene_immune_cor[top_pairs]
  )
  
  tiff(file.path(output_dir, "results", "immune_infiltration", "gene_immune_chord.tiff"),
       width=10, height=10, units="in", res=600, compression="lzw", bg="white")
  circos.clear()
  chordDiagram(chord_data, 
               col=ifelse(chord_data$value > 0, rgb(1,0,0,0.5), rgb(0,0,1,0.5)),
               transparency=0.3,
               annotationTrack=c("grid", "axis"),
               annotationTrackHeight=c(0.03, 0.01))
  title("Gene-Immune Crosstalk")
  dev.off()
  circos.clear()
}

cat("\n\nNECROPTOSIS-IMMUNE CROSSTALK COMPLETE\n")
sink()
