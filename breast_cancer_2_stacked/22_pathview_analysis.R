setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 22_pathview_analysis.R
# KEGG Pathview pathway mapping
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")
source("R/annotation_functions.R")
source("R/enrichment_functions.R")
load_required_packages()

consensus <- readRDS(file.path(CONSENSUS_DIR, "consensus_genes.rds"))
final_genes <- consensus$final_consensus
de_results <- read.csv(file.path(DE_DIR, "differential_expression_full.csv"))
combat_expr <- readRDS(file.path(DISCOVERY_DIR, "discovery_combat.rds"))
meta <- readRDS(file.path(DISCOVERY_DIR, "discovery_cohort.rds"))$metadata

# Load KEGG results
kegg_file <- file.path(ENRICH_DIR, "Consensus_KEGG.csv")
if (!file.exists(kegg_file)) {
    cat("KEGG results not found. Running enrichment first...\n")
    source("scripts/21_GO_KEGG_enrichment.R")
}

enrich_rds <- file.path(ENRICH_DIR, "Consensus_enrichment.rds")
if (!file.exists(enrich_rds)) {
    cat("Enrichment RDS not found.\n")
    cat("Pathview analysis skipped.\n")
    quit(save = "no")
}
enrich_obj <- readRDS(enrich_rds)
kegg_results <- enrich_obj$kegg

if (is.null(kegg_results) || nrow(kegg_results@result) == 0) {
    cat("No KEGG enrichment results available for pathview.\n")
    write.csv(data.frame(Note = "No significant KEGG pathways found"),
              file.path(PATHVIEW_DIR, "pathview_summary.csv"), row.names = FALSE)
    cat("Pathview analysis skipped.\n")
    quit(save = "no")
}

cat("Pathview Analysis\n")
cat("Total KEGG pathways tested:", nrow(kegg_results@result), "\n")

# Calculate fold changes for consensus genes
tumor_mean <- rowMeans(combat_expr[final_genes, meta$class == "Tumor", drop = FALSE], na.rm = TRUE)
normal_mean <- rowMeans(combat_expr[final_genes, meta$class == "Normal", drop = FALSE], na.rm = TRUE)
fc <- tumor_mean - normal_mean  # Already log2 scale
names(fc) <- final_genes

# Generate Pathview maps
pv_results <- generate_pathview_maps(
    kegg_results = kegg_results,
    fold_changes = fc,
    pathway_dir = PATHVIEW_DIR
)

# Save fold change inputs
fc_df <- data.frame(
    Gene = names(fc),
    log2FC = fc,
    Tumor_Mean = tumor_mean,
    Normal_Mean = normal_mean,
    stringsAsFactors = FALSE
)
write.csv(fc_df, file.path(PATHVIEW_DIR, "fold_change_inputs.csv"), row.names = FALSE)

# Move PNGs to pathview directory
pv_files <- list.files(pattern = "^hsa.*\\.png$")
for (f in pv_files) {
    file.rename(f, file.path(PATHVIEW_DIR, f))
}
cat("Pathview maps moved to results/pathview/\n")

cat("\nPathview maps generated.\n")
cat("NOTE: KEGG pathway diagrams generated through Pathview are based on KEGG resources.\n")
cat("Publication and reuse requirements from Kanehisa Laboratories must be checked\n")
cat("before manuscript submission.\n")

