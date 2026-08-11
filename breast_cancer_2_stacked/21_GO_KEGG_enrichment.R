setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 21_GO_KEGG_enrichment.R
# GO and KEGG enrichment analysis of consensus genes
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")
source("R/annotation_functions.R")
source("R/enrichment_functions.R")
source("R/plotting_functions.R")
load_required_packages()

consensus <- readRDS(file.path(CONSENSUS_DIR, "consensus_genes.rds"))
final_genes <- consensus$final_consensus
deg_genes <- readRDS(file.path(DE_DIR, "deg_gene_list.rds"))
de_results <- read.csv(file.path(DE_DIR, "differential_expression_full.csv"))

cat("Enrichment Analysis\n")
cat("Consensus genes:", length(final_genes), "\n")

# Run enrichment for different gene sets
gene_sets <- list(
    Consensus = final_genes,
    All_DEGs = deg_genes,
    Upregulated = de_results$Gene[de_results$DEG == "Up"],
    Downregulated = de_results$Gene[de_results$DEG == "Down"]
)

for (set_name in names(gene_sets)) {
    cat(sprintf("\n--- %s Enrichment ---\n", set_name))
    gl <- gene_sets[[set_name]]
    cat(sprintf("  Input genes: %d\n", length(gl)))
    
    # Remove NAs and empty
    gl <- gl[!is.na(gl) & gl != ""]
    
    if (length(gl) < 5) {
        cat("  Too few genes for enrichment.\n")
        next
    }
    
    enrich_results <- run_full_enrichment(gl, set_name)
    
    cat(sprintf("  Mapped to Entrez: %d/%d\n", 
                enrich_results$n_mapped, enrich_results$n_genes))
    
    # Save results
    for (ont in c("BP", "CC", "MF")) {
        res <- enrich_results[[paste0("go_", ont)]]
        if (!is.null(res)) {
            outfile <- file.path(ENRICH_DIR, 
                                 paste0(set_name, "_GO_", ont, ".csv"))
            write.csv(res@result, outfile, row.names = FALSE)
            cat(sprintf("  GO %s: %d terms\n", ont, nrow(res@result)))
            
            # Save flextable
            if (nrow(res@result) > 0) {
                ft <- flextable(head(res@result, 20))
                ft <- autofit(ft)
                save_as_docx(ft, path = gsub("\\.csv$", ".docx", outfile))
            }
        }
    }
    
    kegg_res <- enrich_results$kegg
    if (!is.null(kegg_res)) {
        outfile <- file.path(ENRICH_DIR, paste0(set_name, "_KEGG.csv"))
        write.csv(kegg_res@result, outfile, row.names = FALSE)
        cat(sprintf("  KEGG: %d pathways\n", nrow(kegg_res@result)))
        
        if (nrow(kegg_res@result) > 0) {
            ft <- flextable(head(kegg_res@result, 20))
            ft <- autofit(ft)
            save_as_docx(ft, path = gsub("\\.csv$", ".docx", outfile))
        }
    }
    
    saveRDS(enrich_results, 
            file.path(ENRICH_DIR, paste0(set_name, "_enrichment.rds")))
}

cat("\nEnrichment analysis complete.\n")

