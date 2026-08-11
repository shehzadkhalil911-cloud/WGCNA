setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 13_consensus_feature_genes.R
# Identify consensus feature genes across all algorithms
# =============================================================================

source("config/analysis_config.R")
source("R/package_functions.R")
load_required_packages()

# Load stable genes from each algorithm
algo_names <- c("RF", "GB", "SVM", "KNN", "NN")
algo_genes <- list()
algo_freq <- list()

for (algo in algo_names) {
    freq_file <- file.path(FS_DIR, algo, paste0(tolower(algo), "_selection_frequency.csv"))
    if (file.exists(freq_file)) {
        freq_df <- read.csv(freq_file)
        algo_freq[[algo]] <- freq_df
        algo_genes[[algo]] <- freq_df$Gene[freq_df$Frequency >= MIN_SELECTION_FREQUENCY]
    } else {
        algo_genes[[algo]] <- character(0)
        algo_freq[[algo]] <- data.frame(Gene = character(0), Frequency = numeric(0))
    }
    cat(sprintf("%s: %d stable genes\n", algo, length(algo_genes[[algo]])))
}

# --- Consensus rules ---
# 1. Strict intersection (all 5)
strict_5 <- Reduce(intersect, algo_genes)
cat(sprintf("\nStrict 5-way intersection: %d genes\n", length(strict_5)))

# 2. At least 4 algorithms
four_of_five <- names(which(table(unlist(algo_genes)) >= 4))
cat(sprintf("At least 4 algorithms: %d genes\n", length(four_of_five)))

# 3. At least 3 algorithms
three_of_five <- names(which(table(unlist(algo_genes)) >= 3))
cat(sprintf("At least 3 algorithms: %d genes\n", length(three_of_five)))

# 4. Weighted consensus
# Combine frequency across algorithms with normalized importance
all_genes <- unique(unlist(algo_genes))
weighted_scores <- data.frame(Gene = all_genes, stringsAsFactors = FALSE)
weighted_scores$Mean_Freq <- sapply(all_genes, function(g) {
    mean(sapply(algo_freq, function(f) {
        if (g %in% f$Gene) f$Frequency[f$Gene == g] else 0
    }), na.rm = TRUE)
})
weighted_scores <- weighted_scores[order(-weighted_scores$Mean_Freq), ]
weighted_consensus <- weighted_scores$Gene[1:max(50, min(200, length(strict_5) + 50))]

cat(sprintf("Weighted consensus: %d genes\n", length(weighted_consensus)))

# --- Select final consensus rule ---
# Use at least 4-of-5 if strict intersection is too small
if (length(strict_5) >= 10) {
    final_consensus <- strict_5
    consensus_rule <- "strict_5"
    cat("Using strict 5-way intersection as final consensus.\n")
} else if (length(four_of_five) >= 10) {
    final_consensus <- four_of_five
    consensus_rule <- "four_of_five"
    cat("Using 4-of-5 consensus as final consensus.\n")
} else if (length(three_of_five) >= 10) {
    final_consensus <- three_of_five
    consensus_rule <- "three_of_five"
    cat("Using 3-of-5 consensus as final consensus.\n")
} else {
    final_consensus <- weighted_consensus[1:min(100, length(weighted_consensus))]
    consensus_rule <- "weighted"
    cat("Using weighted consensus as final consensus.\n")
}

cat("\nFinal Consensus Signature:\n")
cat("  Rule:", consensus_rule, "\n")
cat("  Genes:", length(final_consensus), "\n")
if (length(final_consensus) > 0) {
    cat("  Genes:", paste(head(final_consensus, 20), collapse = ", "), "...\n")
}

# Save all results
consensus_results <- list(
  strict_5 = strict_5,
  four_of_five = four_of_five,
  three_of_five = three_of_five,
  weighted = weighted_consensus,
  final_consensus = final_consensus,
  consensus_rule = consensus_rule,
  algorithm_genes = algo_genes,
  algorithm_frequencies = algo_freq,
  weighted_scores = weighted_scores
)
saveRDS(consensus_results, file.path(CONSENSUS_DIR, "consensus_genes.rds"))

# Write final consensus gene list
write.csv(data.frame(Gene = final_consensus), 
          file.path(CONSENSUS_DIR, "final_consensus_genes.csv"), row.names = FALSE)

# Summary table
summary_df <- data.frame(
  Algorithm = c(algo_names, "Strict_5", "4_of_5", "3_of_5", "Weighted", "Final"),
  N_Genes = c(sapply(algo_genes, length), 
              length(strict_5), length(four_of_five), 
              length(three_of_five), length(weighted_consensus),
              length(final_consensus)),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(CONSENSUS_DIR, "consensus_summary.csv"), row.names = FALSE)
print(summary_df)

cat("\nConsensus gene identification complete.\n")

