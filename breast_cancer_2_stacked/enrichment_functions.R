# =============================================================================
# enrichment_functions.R
# Functions for GO and KEGG enrichment analysis
# =============================================================================

#' Perform GO enrichment analysis
run_go_enrichment <- function(gene_list, ontology = "BP", p_cutoff = 0.05,
                               q_cutoff = 0.2, min_gene = 5, max_gene = 500) {
    entrez_ids <- symbol_to_entrez(gene_list)
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]

    enrichGO(
        gene = as.character(entrez_ids),
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = ontology,
        pAdjustMethod = "BH",
        pvalueCutoff = p_cutoff,
        qvalueCutoff = q_cutoff,
        minGSSize = min_gene,
        maxGSSize = max_gene,
        readable = TRUE
    )
}

#' Perform KEGG pathway enrichment
run_kegg_enrichment <- function(gene_list, p_cutoff = 0.05, q_cutoff = 0.2,
                                 min_gene = 5, max_gene = 500) {
    entrez_ids <- symbol_to_entrez(gene_list)
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]

    enrichKEGG(
        gene = as.character(entrez_ids),
        organism = "hsa",
        keyType = "kegg",
        pAdjustMethod = "BH",
        pvalueCutoff = p_cutoff,
        qvalueCutoff = q_cutoff,
        minGSSize = min_gene,
        maxGSSize = max_gene
    )
}

#' Run full enrichment (GO BP, CC, MF, and KEGG)
run_full_enrichment <- function(gene_list, gene_set_name = "Gene Set") {
    results <- list()
    results$name <- gene_set_name
    results$n_genes <- length(gene_list)
    results$entrez_ids <- symbol_to_entrez(gene_list)
    results$n_mapped <- sum(!is.na(results$entrez_ids))
    results$unmapped <- gene_list[is.na(results$entrez_ids)]

    # GO enrichment
    for (ont in c("BP", "CC", "MF")) {
        results[[paste0("go_", ont)]] <- tryCatch(
            run_go_enrichment(gene_list, ontology = ont),
            error = function(e) NULL
        )
    }

    # KEGG enrichment
    results$kegg <- tryCatch(
        run_kegg_enrichment(gene_list),
        error = function(e) NULL
    )

    return(results)
}

#' Generate Pathview maps
generate_pathview_maps <- function(kegg_results, fold_changes, gene_species = "hsa",
                                    pathway_dir = PATHVIEW_DIR) {
    if (is.null(kegg_results) || nrow(kegg_results@result) == 0) {
        message("No KEGG results to map.")
        return(NULL)
    }

    entrez_ids <- symbol_to_entrez(names(fold_changes))
    fc_for_pathview <- fold_changes * 1  # ensure numeric
    names(fc_for_pathview) <- entrez_ids[names(fold_changes)]
    fc_for_pathview <- fc_for_pathview[!is.na(names(fc_for_pathview)) & 
                                          names(fc_for_pathview) != ""]

    sig_pathways <- kegg_results@result[kegg_results@result$p.adjust < 0.05, ]
    if (nrow(sig_pathways) == 0) sig_pathways <- head(kegg_results@result, 5)
    top_pathways <- head(sig_pathways$ID, 10)
    pathview_results <- list()

    for (pid in top_pathways) {
        tryCatch({
            pv_out <- pathview(
                gene.data = fc_for_pathview,
                pathway.id = pid,
                species = gene_species,
                out.suffix = "breast_cancer",
                kegg.native = TRUE,
                same.layer = TRUE,
                low = list(gene = "green", cpd = "blue"),
                mid = list(gene = "gray", cpd = "gray"),
                high = list(gene = "red", cpd = "yellow"),
                out.dir = pathway_dir
            )
            pathview_results[[pid]] <- pv_out
        }, error = function(e) {
            message("Pathview failed for ", pid, ": ", e$message)
        })
    }
    return(pathview_results)
}

#' Save enrichment tables
save_enrichment_table <- function(enrich_result, filename) {
    if (is.null(enrich_result) || nrow(enrich_result@result) == 0) {
        write.csv(data.frame(Note = "No significant enrichment found"),
                  file = filename, row.names = FALSE)
        return()
    }
    df <- enrich_result@result
    write.csv(df, file = filename, row.names = FALSE)

    # Also save as flextable Word doc
    ft <- flextable(head(df, 20))
    ft <- autofit(ft)
    word_file <- gsub("\\.csv$", ".docx", filename)
    save_as_docx(ft, path = word_file)
}
