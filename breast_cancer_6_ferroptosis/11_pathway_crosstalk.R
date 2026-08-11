source("R/utils.R")
load_libraries()
set.seed(123)

BASE_DIR <- "F:/independent projects/cancer/breast cancer"
OUT_DIR <- file.path(BASE_DIR, "Breast cancer 6 Ferroptosis")

log_msg("=== Step 11: Ferroptosis Pathway Crosstalk Analysis ===")

deg_table <- readRDS(file.path(OUT_DIR, "results/differential_expression/combined_DEG_results.rds"))
de_frg <- readRDS(file.path(OUT_DIR, "results/ferroptosis_genes/DE_FRG_results.rds"))
module_info <- read.csv(file.path(OUT_DIR, "results/WGCNA/module_assignments.csv"), row.names = 1)
MEs <- readRDS(file.path(OUT_DIR, "results/WGCNA/module_eigengenes.rds"))
combat_expr <- readRDS(file.path(OUT_DIR, "data_processed/discovery/discovery_combat.rds"))
ferr_score <- readRDS(file.path(OUT_DIR, "results/ferroptosis_score/ferroptosis_score_vector.rds"))
metadata <- read.csv(file.path(OUT_DIR, "data_processed/discovery/discovery_metadata.csv"))

# Define pathway gene sets for crosstalk analysis
crosstalk_pathways <- list(
  LipidPeroxidation = c("ACSL4", "LPCAT3", "ALOX5", "ALOX12", "ALOX15", "ALOX15B",
                         "GPX4", "PTGS2", "CYBB", "NOX1", "NOX4", "NCF2"),
  Glutathione = c("GCLC", "GCLM", "GSR", "GPX4", "GPX2", "GSS", "SLC7A11",
                   "SLC3A2", "CHAC1", "TXNRD1", "MGST1", "MGST3"),
  IronMetabolism = c("TFRC", "TFR2", "FTH1", "FTL", "FTMT", "STEAP3",
                      "SLC40A1", "HAMP", "HFE", "HMOX1", "TF", "IREB2", "ACO1"),
  OxidativeStress = c("NFE2L2", "KEAP1", "HMOX1", "NQO1", "SQSTM1", "PRDX1",
                       "TXNRD1", "SRXN1", "SESN2", "CAT", "SOD1", "SOD2"),
  Mitochondrial = c("VDAC2", "VDAC3", "ATP5MC3", "CS", "IDH1", "IDH2",
                     "SLC25A28", "SLC25A37", "HSPA5", "LONP1", "CISD1"),
  Hypoxia = c("HIF1A", "EPAS1", "VEGFA", "SLC2A1", "HILPDA", "DDIT4",
               "BNIP3", "BNIP3L", "CA9", "LDHA"),
  Autophagy = c("ATG5", "ATG7", "ATG3", "ATG4D", "ATG13", "ATG16L1",
                 "BECN1", "MAP1LC3A", "GABARAPL1", "GABARAPL2",
                 "SQSTM1", "ULK1", "ULK2", "WIPI1", "WIPI2"),
  Apoptosis = c("TP53", "BID", "BAD", "BAX", "BAK1", "BCL2", "BCL2L1",
                 "CASP3", "CASP8", "CASP9", "CYCS", "FAS", "FASLG"),
  Necroptosis = c("RIPK1", "RIPK3", "MLKL", "CYLD", "FADD", "TRADD"),
  Inflammatory = c("IL6", "IL1B", "TNF", "CXCL2", "IFNG", "TLR4", "MYD88",
                    "NFKB1", "RELA", "TNFAIP3"),
  TNF = c("TNF", "TNFRSF1A", "TNFRSF1B", "TRADD", "FADD", "RIPK1",
           "MAP2K3", "MAP2K6", "MAPK14", "JUN", "FOS"),
  NFkB = c("NFKB1", "NFKB2", "RELA", "RELB", "IKBKB", "IKBKG", "CHUK",
            "TNFAIP3", "NFKBIA", "NFKBIE"),
  p53 = c("TP53", "CDKN1A", "MDM2", "MDM4", "BAX", "BBC3", "PMAIP1",
           "GADD45A", "SESN2", "TIGAR", "SIVA1", "SIAH1")
)

# Check which crosstalk genes are actually DEGs/DE-FRGs
crosstalk_summary <- data.frame(
  Pathway = names(crosstalk_pathways),
  N_genes = sapply(crosstalk_pathways, length),
  N_DEGs = sapply(crosstalk_pathways, function(gs) sum(gs %in% rownames(deg_table)[deg_table$Significant])),
  N_DE_FRGs = sapply(crosstalk_pathways, function(gs) sum(gs %in% rownames(de_frg))),
  stringsAsFactors = FALSE
)
crosstalk_summary$Ratio_DEGs <- crosstalk_summary$N_DEGs / crosstalk_summary$N_genes
crosstalk_summary$Ratio_DE_FRGs <- crosstalk_summary$N_DE_FRGs / crosstalk_summary$N_genes

write.csv(crosstalk_summary, file.path(OUT_DIR, "results/pathway_crosstalk/crosstalk_pathway_summary.csv"), row.names = FALSE)
log_msg(sprintf("Crosstalk pathway summary:\n%s", paste(capture.output(print(crosstalk_summary)), collapse = "\n")))

# Pathway overlap analysis
# Calculate overlap between ferroptosis genes and each crosstalk pathway
pathway_gene_overlap <- matrix(0, nrow = 1, ncol = length(crosstalk_pathways))
colnames(pathway_gene_overlap) <- names(crosstalk_pathways)
rownames(pathway_gene_overlap) <- "Ferroptosis"

for (i in seq_along(crosstalk_pathways)) {
  pathway_gene_overlap[1, i] <- length(intersect(ferr_genes <- readRDS(
    file.path(OUT_DIR, "data_processed/ferroptosis_gene_list/ferroptosis_genes_final.rds")),
    crosstalk_pathways[[i]]))
}

# Load ferr_genes again
ferr_genes <- readRDS(file.path(OUT_DIR, "data_processed/ferroptosis_gene_list/ferroptosis_genes_final.rds"))
for (i in seq_along(crosstalk_pathways)) {
  pathway_gene_overlap[1, i] <- length(intersect(ferr_genes, crosstalk_pathways[[i]]))
}

# Pathway overlay heatmap
open_tiff(file.path(OUT_DIR, "results/figures/Figure9A_pathway_overlap_heatmap.tiff"), width = 8, height = 4)
pheatmap::pheatmap(pathway_gene_overlap, display_numbers = TRUE,
                   color = colorRampPalette(c("white", "red"))(100),
                   main = "Ferroptosis Gene Overlap with Crosstalk Pathways",
                   fontsize = 8, cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

# Module-pathway heatmap
modules_top <- unique(module_info$Module[module_info$Gene %in% rownames(de_frg)])
if (length(modules_top) > 0) {
  mod_pathway_mat <- matrix(0, nrow = length(modules_top), ncol = length(crosstalk_pathways))
  rownames(mod_pathway_mat) <- modules_top
  colnames(mod_pathway_mat) <- names(crosstalk_pathways)

  for (i in seq_along(modules_top)) {
    mod_genes <- module_info$Gene[module_info$Module == modules_top[i]]
    for (j in seq_along(crosstalk_pathways)) {
      mod_pathway_mat[i, j] <- length(intersect(mod_genes, crosstalk_pathways[[j]]))
    }
  }

  if (nrow(mod_pathway_mat) > 1) {
    open_tiff(file.path(OUT_DIR, "results/figures/Figure9D_module_pathway_heatmap.tiff"), width = 9, height = 6)
    pheatmap::pheatmap(mod_pathway_mat, display_numbers = TRUE,
                       color = colorRampPalette(c("white", "red"))(100),
                       main = "Module-Pathway Crosstalk", fontsize = 7)
    dev.off()
    write.csv(mod_pathway_mat, file.path(OUT_DIR, "results/pathway_crosstalk/module_pathway_matrix.csv"))
  }
}

# Correlation-based crosstalk: module eigengenes vs pathway gene scores
pathway_scores <- matrix(NA, nrow = nrow(MEs), ncol = length(crosstalk_pathways))
colnames(pathway_scores) <- names(crosstalk_pathways)
rownames(pathway_scores) <- rownames(MEs)

for (j in seq_along(crosstalk_pathways)) {
  pw_genes <- intersect(crosstalk_pathways[[j]], rownames(combat_expr))
  if (length(pw_genes) >= 3) {
    pathway_scores[, j] <- apply(combat_expr[pw_genes, ], 2, mean)
  }
}

# Module-pathway correlations
mod_pw_cor <- cor(MEs, pathway_scores, method = "spearman", use = "pairwise")
write.csv(mod_pw_cor, file.path(OUT_DIR, "results/pathway_crosstalk/module_pathway_correlations.csv"))

# Gene-pathway chord diagram data
if (nrow(de_frg) > 1 && length(crosstalk_pathways) > 3) {
  chord_data <- data.frame()
  for (pw_name in names(crosstalk_pathways)) {
    common <- intersect(rownames(de_frg), crosstalk_pathways[[pw_name]])
    if (length(common) > 0) {
      chord_data <- rbind(chord_data,
                          data.frame(Gene = common, Pathway = pw_name, stringsAsFactors = FALSE))
    }
  }

  if (nrow(chord_data) > 2) {
    write.csv(chord_data, file.path(OUT_DIR, "results/pathway_crosstalk/gene_pathway_chord_data.csv"), row.names = FALSE)

    # Simple gene-pathway crosstalk heatmap
    chord_mat <- table(chord_data$Gene, chord_data$Pathway)
    if (nrow(chord_mat) > 1 && ncol(chord_mat) > 1) {
      open_tiff(file.path(OUT_DIR, "results/figures/Figure9C_gene_pathway_heatmap.tiff"), width = 9, height = 7)
      pheatmap::pheatmap(as.matrix(chord_mat), display_numbers = TRUE,
                         color = c("white", "red"),
                         breaks = c(-0.5, 0.5, 1.5),
                         main = "Gene-Pathway Membership",
                         cluster_rows = FALSE, cluster_cols = FALSE)
      dev.off()
    }
  }
}

# Save crosstalk pathways definition
saveRDS(crosstalk_pathways, file.path(OUT_DIR, "results/pathway_crosstalk/crosstalk_pathways.rds"))
saveRDS(pathway_scores, file.path(OUT_DIR, "results/pathway_crosstalk/pathway_scores.rds"))

log_msg("Step 11 complete: Ferroptosis pathway crosstalk analysis finished")
