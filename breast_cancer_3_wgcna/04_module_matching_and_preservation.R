###############################################################################
# SCRIPT 4: Cross-dataset Module Matching and Preservation Analysis
# Multi-Cohort WGCNA: Breast Cancer Hub Gene Discovery
###############################################################################
library(WGCNA)
library(flashClust)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)
library(cowplot)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

OUT_DIR <- "F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)"
FIG_DIR <- file.path(OUT_DIR, "figures/preservation")
TAB_DIR <- file.path(OUT_DIR, "tables")
RDA_DIR <- file.path(OUT_DIR, "RData")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

load(file.path(RDA_DIR, "01_processed_data_common.RData"))
load(file.path(RDA_DIR, "01_processed_data_full.RData"))
load(file.path(RDA_DIR, "03_WGCNA_GSE42568.RData"))
load(file.path(RDA_DIR, "03_WGCNA_GSE15852.RData"))
load(file.path(RDA_DIR, "03_WGCNA_GSE29431.RData"))

wgcna_42568 <- wgcna_result; rm(wgcna_result)
load(file.path(RDA_DIR, "03_WGCNA_GSE15852.RData"))
wgcna_15852 <- wgcna_result; rm(wgcna_result)
load(file.path(RDA_DIR, "03_WGCNA_GSE29431.RData"))
wgcna_29431 <- wgcna_result; rm(wgcna_result)

###############################################################################
# Cross-dataset Module Matching
###############################################################################
cat("\n========== Cross-dataset Module Matching ==========\n")

module_gene_list <- function(wgcna_obj, gene_names = NULL) {
  if (is.null(gene_names)) {
    gene_names <- colnames(wgcna_obj$datExpr)
  }
  names(wgcna_obj$moduleColors) <- gene_names
  split(gene_names, wgcna_obj$moduleColors)
}

# Get common gene names from each dataset
genes_42568 <- colnames(wgcna_42568$datExpr)
genes_15852 <- colnames(wgcna_15852$datExpr)
genes_29431 <- colnames(wgcna_29431$datExpr)

modules_42568 <- module_gene_list(wgcna_42568, genes_42568)
modules_15852 <- module_gene_list(wgcna_15852, genes_15852)
modules_29431 <- module_gene_list(wgcna_29431, genes_29431)

# Convert to common gene space - find intersection of wgcna genes with common_gene_set_discovery
# The WGCNA was run on common_gene_set_discovery, so module colors align directly
common_gene_set <- common_genes_discovery

get_module_overlap <- function(mod_genes, gene_set) {
  lapply(mod_genes, function(g) intersect(g, gene_set))
}

mod42568_c <- get_module_overlap(modules_42568, common_gene_set)
mod15852_c <- get_module_overlap(modules_15852, common_gene_set)
mod29431_c <- get_module_overlap(modules_29431, common_gene_set)

# Calculate pairwise Jaccard similarities and Fisher's exact test
pairwise_module_match <- function(mod_list1, mod_list2, set_name1, set_name2) {
  n1 <- length(mod_list1)
  n2 <- length(mod_list2)
  
  jaccard_matrix <- matrix(0, n1, n2, 
                           dimnames = list(names(mod_list1), names(mod_list2)))
  fisher_pval <- matrix(1, n1, n2,
                        dimnames = list(names(mod_list1), names(mod_list2)))
  
  for (i in seq_along(mod_list1)) {
    for (j in seq_along(mod_list2)) {
      g1 <- mod_list1[[i]]
      g2 <- mod_list2[[j]]
      intersection <- length(intersect(g1, g2))
      union <- length(union(g1, g2))
      jaccard_matrix[i, j] <- ifelse(union > 0, intersection / union, 0)
      
      # Fisher's exact test
      n_total <- length(common_gene_set)
      a <- intersection
      b <- length(g1) - a
      c <- length(g2) - a
      d <- n_total - a - b - c
      ft <- fisher.test(matrix(c(a, b, c, d), 2, 2), alternative = "greater")
      fisher_pval[i, j] <- ft$p.value
    }
  }
  list(jaccard = jaccard_matrix, fisher_pval = fisher_pval)
}

match_42568_15852 <- pairwise_module_match(mod42568_c, mod15852_c, "GSE42568", "GSE15852")
match_42568_29431 <- pairwise_module_match(mod42568_c, mod29431_c, "GSE42568", "GSE29431")
match_15852_29431 <- pairwise_module_match(mod15852_c, mod29431_c, "GSE15852", "GSE29431")

# Function to find best matches
find_best_matches <- function(jaccard_mat, fisher_mat, min_jaccard = 0.05, max_pval = 0.05) {
  matches <- list()
  for (i in 1:nrow(jaccard_mat)) {
    for (j in 1:ncol(jaccard_mat)) {
      if (jaccard_mat[i, j] >= min_jaccard && fisher_mat[i, j] <= max_pval) {
        matches[[length(matches) + 1]] <- data.frame(
          Module1 = rownames(jaccard_mat)[i],
          Module2 = colnames(jaccard_mat)[j],
          Jaccard = jaccard_mat[i, j],
          Fisher_Pval = fisher_mat[i, j]
        )
      }
    }
  }
  if (length(matches) > 0) {
    do.call(rbind, matches)
  } else {
    data.frame()
  }
}

best_42568_15852 <- find_best_matches(match_42568_15852$jaccard, match_42568_15852$fisher_pval)
best_42568_29431 <- find_best_matches(match_42568_29431$jaccard, match_42568_29431$fisher_pval)
best_15852_29431 <- find_best_matches(match_15852_29431$jaccard, match_15852_29431$fisher_pval)

cat(sprintf("Matched module pairs:\n"))
cat(sprintf("  GSE42568-GSE15852: %d\n", nrow(best_42568_15852)))
cat(sprintf("  GSE42568-GSE29431: %d\n", nrow(best_42568_29431)))
cat(sprintf("  GSE15852-GSE29431: %d\n", nrow(best_15852_29431)))

write.csv(best_42568_15852, file.path(TAB_DIR, "module_match_42568_vs_15852.csv"), row.names = FALSE)
write.csv(best_42568_29431, file.path(TAB_DIR, "module_match_42568_vs_29431.csv"), row.names = FALSE)
write.csv(best_15852_29431, file.path(TAB_DIR, "module_match_15852_vs_29431.csv"), row.names = FALSE)

# Module overlap heatmaps
plot_module_overlap <- function(jaccard_mat, name1, name2, filename) {
  # Remove grey modules for clarity
  keep_rows <- !grepl("grey", rownames(jaccard_mat))
  keep_cols <- !grepl("grey", colnames(jaccard_mat))
  if (sum(keep_rows) > 0 && sum(keep_cols) > 0) {
    jaccard_sub <- jaccard_mat[keep_rows, keep_cols, drop = FALSE]
    pheatmap::pheatmap(jaccard_sub, 
                       color = colorRampPalette(c("white", "yellow", "orange", "red"))(50),
                       main = paste("Module Overlap -", name1, "vs", name2),
                       display_numbers = TRUE, number_format = "%.2f",
                       fontsize_number = 6,
                       filename = file.path(FIG_DIR, filename),
                       width = 10, height = 8)
  }
}

plot_module_overlap(match_42568_15852$jaccard, "GSE42568", "GSE15852", 
                     "module_overlap_42568_15852.png")
plot_module_overlap(match_42568_29431$jaccard, "GSE42568", "GSE29431",
                     "module_overlap_42568_29431.png")
plot_module_overlap(match_15852_29431$jaccard, "GSE15852", "GSE29431",
                     "module_overlap_15852_29431.png")

# Identify consensus modules across all 3 datasets
# A consensus module = a set of modules that match across datasets
find_consensus_modules <- function(matches12, matches13, matches23) {
  all_1 <- unique(c(as.character(matches12$Module1), as.character(matches13$Module1)))
  consensus <- list()
  for (m1 in all_1) {
    m2_candidates <- matches12$Module2[matches12$Module1 == m1]
    m3_candidates <- matches13$Module2[matches13$Module1 == m1]
    
    for (m2 in m2_candidates) {
      m3_via_12 <- matches23$Module2[matches23$Module1 == m2]
      m3_for_m1 <- intersect(m3_candidates, m3_via_12)
      if (length(m3_for_m1) > 0) {
        consensus[[length(consensus) + 1]] <- list(
          GSE42568 = m1,
          GSE15852 = m2,
          GSE29431 = m3_for_m1[1]
        )
      }
    }
  }
  # Remove duplicates
  consensus_unique <- unique(lapply(consensus, function(x) paste(sort(unlist(x)), collapse = "_")))
  consensus[!duplicated(lapply(consensus, function(x) paste(sort(unlist(x)), collapse = "_")))]
}

# Get matched pairs (best per module, with highest Jaccard)
best_match_per_module <- function(df, col1, col2) {
  best <- df %>% group_by(!!sym(col1)) %>% slice_max(Jaccard, n = 1) %>% ungroup()
  as.data.frame(best)
}

best_12 <- best_match_per_module(best_42568_15852, "Module1", "Module2")
best_13 <- best_match_per_module(best_42568_29431, "Module1", "Module2")
best_23 <- best_match_per_module(best_15852_29431, "Module1", "Module2")

consensus_modules <- find_consensus_modules(best_12, best_13, best_23)
cat(sprintf("\nConsensus modules (matched across all 3): %d\n", length(consensus_modules)))
for (cm in consensus_modules) {
  cat(sprintf("  %s - %s - %s\n", cm$GSE42568, cm$GSE15852, cm$GSE29431))
}

save(consensus_modules, match_42568_15852, match_42568_29431, match_15852_29431,
     best_42568_15852, best_42568_29431, best_15852_29431,
     modules_42568, modules_15852, modules_29431,
     file = file.path(RDA_DIR, "04_module_matching.RData"))

###############################################################################
# Module Preservation Analysis
###############################################################################
cat("\n\n========== Module Preservation Analysis ==========\n")

# Use GSE42568 as reference, test in the other two
run_preservation <- function(ref_expr, test_expr, ref_colors, ref_name, test_name) {
  cat(sprintf("Preservation: %s (ref) -> %s (test)\n", ref_name, test_name))
  
  # Multi-set format
  multiExpr <- list(
    ref = list(data = t(ref_expr)),
    test = list(data = t(test_expr))
  )
  names(multiExpr$ref$data) <- NULL
  names(multiExpr$test$data) <- NULL
  
  # Save temp files to a simple path without spaces
  old_wd <- getwd()
  setwd("C:/Users/shehz/AppData/Local/Temp/opencode/wgcna_pres")
  mp <- modulePreservation(multiExpr, list(ref = ref_colors),
                           referenceNetworks = 1,
                           nPermutations = 100,
                           networkType = "signed",
                           quickCor = 1,
                           verbose = 3)
  if (is.null(mp)) return(NULL)
  
  # Extract preservation statistics
  # For referenceNetworks=1, the structure is:
  # mp$preservation$observed[[ref_idx]][[set_idx]]
  # or mp$preservation$Z[[ref_idx]][[set_idx]]
  ref <- mp$preservation$observed[[1]][[1]]
  # For Z statistics
  z <- mp$preservation$Z[[1]][[1]]
  # Add test statistics if available
  test <- NULL
  tryCatch({
    test <- mp$preservation$observed[[1]][[2]]
  }, error = function(e) {})

  setwd(old_wd)
  
  list(ref = ref, test = test, z = z, mp = mp)
}

# Use common genes for preservation
ref_42568 <- comm_42568$expr
ref_15852 <- comm_15852$expr
ref_29431 <- comm_29431$expr

# Map module colors to common genes
# Colors are in same order as colnames(datExpr) = genes used in WGCNA
names(wgcna_42568$moduleColors) <- colnames(wgcna_42568$datExpr)
names(wgcna_15852$moduleColors) <- colnames(wgcna_15852$datExpr)
names(wgcna_29431$moduleColors) <- colnames(wgcna_29431$datExpr)

common42568_colors <- wgcna_42568$moduleColors[rownames(ref_42568)]
common42568_colors[is.na(common42568_colors)] <- "grey"

common15852_colors <- wgcna_15852$moduleColors[rownames(ref_15852)]
common15852_colors[is.na(common15852_colors)] <- "grey"

common29431_colors <- wgcna_29431$moduleColors[rownames(ref_29431)]
common29431_colors[is.na(common29431_colors)] <- "grey"

# Run preservation
pres_42568_to_15852 <- run_preservation(ref_42568, ref_15852, common42568_colors, 
                                          "GSE42568", "GSE15852")
pres_42568_to_29431 <- run_preservation(ref_42568, ref_29431, common42568_colors,
                                          "GSE42568", "GSE29431")
pres_15852_to_29431 <- run_preservation(ref_15852, ref_29431, common15852_colors,
                                          "GSE15852", "GSE29431")

# Plot preservation results
plot_preservation <- function(pres_obj, ref_name, test_name) {
  z <- pres_obj$z
  z <- z[!grepl("grey", rownames(z)), , drop = FALSE]
  
  # Zsummary plot
  df <- data.frame(Module = rownames(z), 
                   Zsummary = z[, "Zsummary.pres"],
                   medianRank = z[, "medianRank.pres"])
  df$Module <- gsub("^ME", "", df$Module)
  df <- df[order(df$Zsummary, decreasing = TRUE), ]
  
  df$color <- ifelse(df$Zsummary > 10, "Strong",
                     ifelse(df$Zsummary > 2, "Moderate", "Weak"))
  df$color <- factor(df$color, levels = c("Strong", "Moderate", "Weak"))
  
  p1 <- ggplot(df, aes(x = reorder(Module, Zsummary), y = Zsummary, fill = color)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = c(2, 10), linetype = "dashed", color = c("red", "darkgreen")) +
    scale_fill_manual(values = c("Strong" = "darkgreen", "Moderate" = "orange", "Weak" = "red")) +
    labs(title = paste("Module Preservation:", ref_name, "->", test_name),
         subtitle = "Zsummary",
         x = "Module", y = "Zsummary") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  p2 <- ggplot(df, aes(x = reorder(Module, medianRank), y = medianRank)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = paste("Module Preservation:", ref_name, "->", test_name),
         subtitle = "medianRank (lower = better preserved)",
         x = "Module", y = "medianRank") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  list(zsummary_plot = p1, medianrank_plot = p2)
}

pres_plots_1 <- plot_preservation(pres_42568_to_15852, "GSE42568", "GSE15852")
pres_plots_2 <- plot_preservation(pres_42568_to_29431, "GSE42568", "GSE29431")
pres_plots_3 <- plot_preservation(pres_15852_to_29431, "GSE15852", "GSE29431")

# Save Zsummary and medianRank tables
for (pres_name in c("pres_42568_to_15852", "pres_42568_to_29431", "pres_15852_to_29431")) {
  pres_obj <- get(pres_name)
  z_df <- as.data.frame(pres_obj$z)
  z_df$Module <- rownames(z_df)
  write.csv(z_df, file.path(TAB_DIR, paste0("preservation_", pres_name, ".csv")), row.names = FALSE)
}

# Save plots
for (i in 1:3) {
  preset <- list(pres_42568_to_15852 = pres_plots_1, 
                 pres_42568_to_29431 = pres_plots_2,
                 pres_15852_to_29431 = pres_plots_3)[[i]]
  
  ggsave(file.path(FIG_DIR, paste0("Zsummary_", names(preset)[i], ".png")), 
         preset$zsummary_plot, width = 10, height = 6, dpi = 300)
  ggsave(file.path(FIG_DIR, paste0("medianRank_", names(preset)[i], ".png")),
         preset$medianrank_plot, width = 10, height = 6, dpi = 300)
}

# Create combined preservation summary
cat("\nPreservation Summary (Zsummary):\n")
for (pres_name in c("pres_42568_to_15852", "pres_42568_to_29431", "pres_15852_to_29431")) {
  pres_obj <- get(pres_name)
  z <- pres_obj$z
  cat(sprintf("\n%s:\n", pres_name))
  z_sorted <- z[order(-z[, "Zsummary.pres"]), ]
  cat(sprintf("  Top preserved: %s (Z=%.1f)\n", 
              gsub("ME", "", rownames(z_sorted)[1]), z_sorted[1, "Zsummary.pres"]))
  cat(sprintf("  Strongly preserved (Z>10): %d\n", sum(z[, "Zsummary.pres"] > 10, na.rm=TRUE)))
  cat(sprintf("  Moderately preserved (2<Z<=10): %d\n", 
              sum(z[, "Zsummary.pres"] > 2 & z[, "Zsummary.pres"] <= 10, na.rm=TRUE)))
}

save(pres_42568_to_15852, pres_42568_to_29431, pres_15852_to_29431,
     file = file.path(RDA_DIR, "04_preservation.RData"))
cat("\nScript 4 complete.\n")
