# ============================================================
# 13 - Functional enrichment
# ORA (GO BP/MF/CC, KEGG, Reactome) of final stable features
# against the common-gene universe, plus GSEA on the paired-DEG
# ranked list. Table_21 + Figure_12.
# ============================================================
PROJECT_DIR <- "F:/independent projects/cancer/breast cancer/Breast cancer 8"
source(file.path(PROJECT_DIR, "config/config.R"))
source(file.path(PROJECT_DIR, "R/utils.R"))
suppressMessages({
  library(dplyr); library(ggplot2); library(clusterProfiler)
  library(org.Hs.eg.db); library(ReactomePA); library(enrichplot)
})

logmsg <- function(...) cat(sprintf(...), "\n")

nc  <- readRDS(file.path(DIR_MODELS, "nested_cv_results.rds"))
deg <- readRDS(file.path(DIR_MODELS, "deg_results.rds"))
common <- readRDS(file.path(DIR_MODELS, "common_genes.rds"))

stable <- nc$final_stable
logmsg("Stable genes: %d", length(stable))

# ---------------- map symbols to ENTREZ ----------------
map_sym <- function(syms) {
  AnnotationDbi::mapIds(org.Hs.eg.db, keys = syms, column = "ENTREZID",
                        keytype = "SYMBOL", multiVals = "first")
}
stable_entrez <- map_sym(stable)
univ_entrez   <- map_sym(common)
stable_entrez <- stable_entrez[!is.na(stable_entrez)]
univ_entrez   <- univ_entrez[!is.na(univ_entrez)]
logmsg("Mapped stable: %d, universe: %d", length(stable_entrez), length(univ_entrez))

# ---------------- ORA ----------------
run_ora <- function(term) {
  tryCatch({
    if (term == "BP")   enrichGO(stable_entrez, universe = univ_entrez, OrgDb = org.Hs.eg.db,
                                 ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
    else if (term == "MF") enrichGO(stable_entrez, universe = univ_entrez, OrgDb = org.Hs.eg.db,
                                    ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
    else if (term == "CC") enrichGO(stable_entrez, universe = univ_entrez, OrgDb = org.Hs.eg.db,
                                    ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
    else if (term == "KEGG") enrichKEGG(stable_entrez, universe = univ_entrez,
                                        organism = "hsa", pAdjustMethod = "BH",
                                        qvalueCutoff = 0.05)
    else if (term == "Reactome") enrichPathway(stable_entrez, universe = univ_entrez,
                                               organism = "human", pAdjustMethod = "BH",
                                               qvalueCutoff = 0.05, readable = TRUE)
  }, error = function(e) NULL)
}

ora <- lapply(c("BP", "MF", "CC", "KEGG", "Reactome"), function(t) {
  r <- run_ora(t)
  if (is.null(r) || nrow(as.data.frame(r)) == 0) {
    logmsg("ORA %s: no significant terms", t); return(NULL)
  }
  logmsg("ORA %s: %d significant terms", t, nrow(as.data.frame(r)))
  r
})
names(ora) <- c("BP", "MF", "CC", "KEGG", "Reactome")

# combined table
ora_tab <- do.call(rbind, lapply(names(ora), function(t) {
  r <- ora[[t]]
  if (is.null(r)) return(NULL)
  df <- as.data.frame(r)
  df <- df[order(df$p.adjust), ]
  df$source <- t
  df$GeneRatio_str <- paste0(df$GeneRatio, " (", df$BgRatio, ")")
  df[c("source", "ID", "Description", "GeneRatio_str", "pvalue", "p.adjust", "qvalue",
       "geneID", "Count")]
}))
if (!is.null(ora_tab)) {
  write.csv(ora_tab, file.path(DIR_TABLES, "Table_21_Enrichment_ORA.csv"), row.names = FALSE)
  logmsg("Table_21 written: %d rows", nrow(ora_tab))
}

# ---------------- GSEA ----------------
d <- deg$paired
d$rank_score <- -log10(d$P.Value) * sign(d$logFC)
d <- d[order(-d$rank_score), ]
gene_list <- setNames(d$rank_score, d$gene)
gene_list <- gene_list[!duplicated(names(gene_list))]
gene_list <- gene_list[is.finite(gene_list)]

gsea_go <- tryCatch(gseGO(gene_list, OrgDb = org.Hs.eg.db, ont = "BP",
                          keyType = "SYMBOL", minGSSize = 10, maxGSSize = 500,
                          pAdjustMethod = "BH", seed = SEED, verbose = FALSE),
                    error = function(e) NULL)
gene_list_entrez <- setNames(d$rank_score, map_sym(d$gene))
gene_list_entrez <- gene_list_entrez[!is.na(names(gene_list_entrez))]
gene_list_entrez <- gene_list_entrez[!duplicated(names(gene_list_entrez))]
gene_list_entrez <- gene_list_entrez[is.finite(gene_list_entrez)]
gsea_kegg <- tryCatch(gseKEGG(gene_list_entrez, organism = "hsa", minGSSize = 10,
                              maxGSSize = 500, pAdjustMethod = "BH", seed = SEED),
                      error = function(e) NULL)
gsea <- list(BP = gsea_go, KEGG = gsea_kegg)

gsea_tab <- do.call(rbind, lapply(names(gsea), function(t) {
  g <- gsea[[t]]
  if (is.null(g) || nrow(as.data.frame(g)) == 0) return(NULL)
  df <- as.data.frame(g)
  df <- df[order(df$p.adjust), ]
  df$source <- t
  df[c("source", "ID", "Description", "setSize", "enrichmentScore", "NES",
       "pvalue", "p.adjust", "qvalue", "rank", "core_enrichment")]
}))
if (!is.null(gsea_tab)) {
  write.csv(gsea_tab, file.path(DIR_TABLES, "Table_21b_Enrichment_GSEA.csv"), row.names = FALSE)
  logmsg("Table_21b written: %d rows", nrow(gsea_tab))
}

saveRDS(list(ora = ora, gsea = gsea, ora_table = ora_tab, gsea_table = gsea_tab,
             stable_entrez = stable_entrez, universe_entrez = univ_entrez),
        file.path(DIR_MODELS, "enrichment_results.rds"))

# ---------------- figures (Figure_12) ----------------
fig <- function(p, name, w = 8, h = 6) save_tiff_600dpi(p, name, w, h)

if (!is.null(ora$BP)) {
  p <- dotplot(ora$BP, showCategory = 20) + theme_bc(base_size = 9) +
    labs(title = "GO Biological Process (ORA)") +
    theme(legend.position = "right")
  fig(p, "Figure_12_Enrichment_GO_BP.tiff")
}
if (!is.null(ora$MF)) {
  p <- dotplot(ora$MF, showCategory = 15) + theme_bc(base_size = 9) +
    labs(title = "GO Molecular Function (ORA)") +
    theme(legend.position = "right")
  fig(p, "Figure_12_Enrichment_GO_MF.tiff", w = 8, h = 5.5)
}
if (!is.null(ora$CC)) {
  p <- dotplot(ora$CC, showCategory = 15) + theme_bc(base_size = 9) +
    labs(title = "GO Cellular Component (ORA)") +
    theme(legend.position = "right")
  fig(p, "Figure_12_Enrichment_GO_CC.tiff", w = 8, h = 5.5)
}
if (!is.null(ora$KEGG)) {
  p <- dotplot(ora$KEGG, showCategory = 15) + theme_bc(base_size = 9) +
    labs(title = "KEGG pathway (ORA)") +
    theme(legend.position = "right")
  fig(p, "Figure_12_Enrichment_KEGG.tiff", w = 8, h = 5.5)
}
if (!is.null(ora$Reactome)) {
  p <- dotplot(ora$Reactome, showCategory = 15) + theme_bc(base_size = 9) +
    labs(title = "Reactome pathway (ORA)") +
    theme(legend.position = "right")
  fig(p, "Figure_12_Enrichment_Reactome.tiff", w = 8, h = 5.5)
}
if (!is.null(gsea$BP)) {
  p <- gseaplot2(gsea$BP, geneSetID = head(seq_len(nrow(as.data.frame(gsea$BP))), 3),
                 pvalue_table = TRUE) + theme_bc(base_size = 9) +
    labs(title = "GSEA GO:BP top pathways") +
    theme(legend.position = "right")
  fig(p, "Figure_12_GSEA_GO_BP.tiff", w = 8, h = 6)
}
if (!is.null(gsea$KEGG)) {
  p <- gseaplot2(gsea$KEGG, geneSetID = head(seq_len(nrow(as.data.frame(gsea$KEGG))), 3),
                 pvalue_table = TRUE) + theme_bc(base_size = 9) +
    labs(title = "GSEA KEGG top pathways") +
    theme(legend.position = "right")
  fig(p, "Figure_12_GSEA_KEGG.tiff", w = 8, h = 6)
}

logmsg("Enrichment analysis complete.")
