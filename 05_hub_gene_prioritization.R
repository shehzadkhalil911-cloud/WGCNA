###############################################################################
# SCRIPT 5: Hub-gene Prioritization, Integration with DEGs, and Validation Prep
# Multi-Cohort WGCNA: Breast Cancer Hub Gene Discovery
###############################################################################
library(WGCNA)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(pheatmap)
library(limma)
library(pROC)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

OUT_DIR <- "F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)"
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")
RDA_DIR <- file.path(OUT_DIR, "RData")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

load(file.path(RDA_DIR, "01_processed_data_common.RData"))
load(file.path(RDA_DIR, "01_processed_data_full.RData"))
load(file.path(RDA_DIR, "02_DE_analysis.RData"))
load(file.path(RDA_DIR, "03_WGCNA_GSE42568.RData")); wgcna_42568 <- wgcna_result; rm(wgcna_result)
load(file.path(RDA_DIR, "03_WGCNA_GSE15852.RData")); wgcna_15852 <- wgcna_result; rm(wgcna_result)
load(file.path(RDA_DIR, "03_WGCNA_GSE29431.RData")); wgcna_29431 <- wgcna_result; rm(wgcna_result)
load(file.path(RDA_DIR, "04_module_matching.RData"))

cat("========== Hub-gene Prioritization ==========\n")

# Function to compute hub scores
compute_hub_score <- function(wgcna, de_results_list, dataset_name, genes) {
  mm <- abs(wgcna$MM[, grep("^ME", colnames(wgcna$MM))])
  gs <- abs(wgcna$GS)
  colnames(gs) <- "GS"
  
  # Get intramodular connectivity
  deg_table <- de_results_list$results
  deg_table <- deg_table[match(genes, deg_table$Gene), ]
  
  # Create score data frame
  score_df <- data.frame(
    Gene = genes,
    Module = wgcna$moduleColors[genes],
    MM_score = NA,
    GS_score = NA,
    DEG_logFC = deg_table$logFC,
    DEG_Pval = deg_table$adj.P.Val,
    DEG_sig = ifelse(deg_table$adj.P.Val < 0.05, 1, 0),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(genes)) {
    g <- genes[i]
    if (g %in% rownames(mm)) {
      score_df$MM_score[i] <- max(abs(mm[g, ]), na.rm = TRUE)
    }
    if (g %in% rownames(gs)) {
      score_df$GS_score[i] <- gs[g, 1]
    }
  }
  
  # Composite score (equal weights)
  score_df$MM_score_scaled <- scale(score_df$MM_score)[,1]
  score_df$GS_score_scaled <- scale(score_df$GS_score)[,1]
  score_df$HubScore <- score_df$MM_score_scaled + score_df$GS_score_scaled + 
    ifelse(score_df$DEG_sig > 0, 1, 0) * abs(score_df$DEG_logFC) / max(abs(score_df$DEG_logFC), na.rm = TRUE)
  
  score_df <- score_df[order(-score_df$HubScore), ]
  score_df
}

common_genes <- common_genes_discovery

# Get hub scores for each dataset
hub_42568 <- compute_hub_score(wgcna_42568, de_42568, "GSE42568", common_genes)
hub_15852 <- compute_hub_score(wgcna_15852, de_15852, "GSE15852", common_genes)
hub_29431 <- compute_hub_score(wgcna_29431, de_29431, "GSE29431", common_genes)

# Identify consensus hub genes (top in at least 2 datasets)
top_n <- 200
top_42568 <- hub_42568$Gene[1:min(top_n, nrow(hub_42568))]
top_15852 <- hub_15852$Gene[1:min(top_n, nrow(hub_15852))]
top_29431 <- hub_29431$Gene[1:min(top_n, nrow(hub_29431))]

consensus_hubs <- intersect(intersect(top_42568, top_15852), top_29431)
consistent_hubs <- names(which(table(c(top_42568, top_15852, top_29431)) >= 2))

cat(sprintf("Top %d hubs per dataset:\n", top_n))
cat(sprintf("  Consensus (all 3): %d\n", length(consensus_hubs)))
cat(sprintf("  Consistent (>=2): %d\n", length(consistent_hubs)))

# Intersect with DEGs for WGCNA-supported DEGs
deg_all3 <- Reduce(intersect, list(
  de_42568$results$Gene[de_42568$results$adj.P.Val < 0.05],
  de_15852$results$Gene[de_15852$results$adj.P.Val < 0.05],
  de_29431$results$Gene[de_29431$results$adj.P.Val < 0.05]
))

wgcna_supported_deg <- intersect(consistent_hubs, deg_all3)
cat(sprintf("WGCNA-supported DEGs: %d\n", length(wgcna_supported_deg)))

# Final hub gene list
final_hubs <- consistent_hubs[1:min(50, length(consistent_hubs))]
cat(sprintf("Final candidate hub genes: %d\n", length(final_hubs)))
writeLines(final_hubs, file.path(TAB_DIR, "final_hub_genes.txt"))

# Write hub gene tables
write.csv(hub_42568, file.path(TAB_DIR, "hub_scores_GSE42568.csv"), row.names = FALSE)
write.csv(hub_15852, file.path(TAB_DIR, "hub_scores_GSE15852.csv"), row.names = FALSE)
write.csv(hub_29431, file.path(TAB_DIR, "hub_scores_GSE29431.csv"), row.names = FALSE)

# Summary table of final hubs
hub_summary <- data.frame(
  Gene = final_hubs,
  Module_42568 = wgcna_42568$moduleColors[final_hubs],
  Module_15852 = wgcna_15852$moduleColors[final_hubs],
  Module_29431 = wgcna_29431$moduleColors[final_hubs],
  logFC_42568 = de_42568$results$logFC[match(final_hubs, de_42568$results$Gene)],
  logFC_15852 = de_15852$results$logFC[match(final_hubs, de_15852$results$Gene)],
  logFC_29431 = de_29431$results$logFC[match(final_hubs, de_29431$results$Gene)],
  meta_P = meta_p$meta_P[match(final_hubs, meta_p$Gene)],
  meta_FDR = meta_p$meta_FDR[match(final_hubs, meta_p$Gene)]
)
write.csv(hub_summary, file.path(TAB_DIR, "hub_gene_summary.csv"), row.names = FALSE)

# Save before consensus module part
save(hub_42568, hub_15852, hub_29431, final_hubs, consensus_hubs, consistent_hubs,
     wgcna_supported_deg, hub_summary,
     file = file.path(RDA_DIR, "05_hub_genes.RData"))

# Consensus module gene lists
for (cm in consensus_modules) {
  cm_name <- paste(cm$GSE42568, cm$GSE15852, cm$GSE29431, sep = "_")
  # Get genes in the intersection of matching modules across all 3 datasets
  mod1_genes <- names(wgcna_42568$moduleColors[wgcna_42568$moduleColors == cm$GSE42568])
  mod2_genes <- names(wgcna_15852$moduleColors[wgcna_15852$moduleColors == cm$GSE15852])
  mod3_genes <- names(wgcna_29431$moduleColors[wgcna_29431$moduleColors == cm$GSE29431])
  genes_in_mod <- Reduce(intersect, list(mod1_genes, mod2_genes, mod3_genes))
  if (length(genes_in_mod) > 0) {
    writeLines(as.character(genes_in_mod), 
               file.path(TAB_DIR, paste0("consensus_module_", cm_name, "_genes.txt")))
    cat(sprintf("Consensus module %s: %d genes\n", cm_name, length(genes_in_mod)))
  }
}

cat("\nScript 5 complete.\n")
