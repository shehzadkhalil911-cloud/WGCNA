source("R/utils.R")
load_libraries()
set.seed(123)

BASE_DIR <- "F:/independent projects/cancer/breast cancer"
OUT_DIR <- file.path(BASE_DIR, "Breast cancer 6 Ferroptosis")

log_msg("=== Step 15: Module Preservation Analysis ===")

combat_expr <- readRDS(file.path(OUT_DIR, "data_processed/discovery/discovery_combat.rds"))
metadata <- read.csv(file.path(OUT_DIR, "data_processed/discovery/discovery_metadata.csv"))
module_info <- read.csv(file.path(OUT_DIR, "results/WGCNA/module_assignments.csv"))
rownames(module_info) <- module_info$Gene
discovery_list <- readRDS(file.path(OUT_DIR, "data_processed/discovery/discovery_list.rds"))

library(WGCNA)

# Per-dataset module preservation
preservation_results <- list()

for (ds in names(discovery_list)) {
  d <- discovery_list[[ds]]
  common <- intersect(rownames(d$expr), module_info$Gene)
  if (length(common) < 50) {
    log_msg(sprintf("  %s: too few common genes (%d), skipping", ds, length(common)))
    next
  }

  # Get module assignments for common genes
  mod_common <- module_info[match(common, module_info$Gene), ]
  expr_ds <- d$expr[common, ]

  # Create multiExpr
  multiExpr <- list(
    Reference = list(data = t(combat_expr[common, ])),
    Test = list(data = t(expr_ds))
  )

  # Filter zero-variance genes/samples
  multiExpr <- lapply(multiExpr, function(set) {
    gsg <- goodSamplesGenes(set$data, verbose = 0)
    if (!gsg$allOK) {
      set$data <- set$data[gsg$goodSamples, gsg$goodGenes]
    }
    set
  })

  # Ensure matching genes across reference and test
  common_genes_pres <- intersect(colnames(multiExpr$Reference$data), colnames(multiExpr$Test$data))
  multiExpr$Reference$data <- multiExpr$Reference$data[, common_genes_pres, drop = FALSE]
  multiExpr$Test$data <- multiExpr$Test$data[, common_genes_pres, drop = FALSE]
  mod_common_pres <- mod_common[match(common_genes_pres, mod_common$Gene), ]

  # Run module preservation
  pres <- tryCatch({
    modulePreservation(multiExpr, multiColor = list(mod_common_pres$Module, mod_common_pres$Module),
                       nPermutations = 30, verbose = 0, networkType = "signed",
                       referenceNetworks = 1, maxGoldModuleSize = 100)
  }, error = function(e) {
    log_msg(sprintf("  Module preservation failed for %s: %s", ds, e$message))
    NULL
  })

  if (!is.null(pres)) {
    ref <- 1
    test <- 2
    z_scores <- pres$preservation$Z[[ref]][[test]]
    median_ranks <- if (!is.null(pres$preservation$medianRank)) pres$preservation$medianRank[[ref]][[test]] else NULL

    if (!is.null(median_ranks)) {
      preservation_results[[ds]] <- list(
        Zsummary = z_scores$Zsummary.pres,
        medianRank = median_ranks$medianRank.pres,
        modules = rownames(z_scores)
      )
    } else {
      preservation_results[[ds]] <- list(
        Zsummary = z_scores$Zsummary.pres,
        medianRank = NULL,
        modules = rownames(z_scores)
      )
    }

    # Save preservation stats
    pres_df <- data.frame(
      Module = rownames(z_scores),
      Zsummary = z_scores$Zsummary.pres,
      medianRank = if (!is.null(median_ranks)) median_ranks$medianRank.pres else NA,
      stringsAsFactors = FALSE
    )
    write.csv(pres_df, file.path(OUT_DIR, sprintf("results/module_preservation/preservation_%s.csv", ds)), row.names = FALSE)

    log_msg(sprintf("  %s: %d modules assessed", ds, nrow(pres_df)))
  }
}

# Zsummary plot
zsummary_df <- do.call(rbind, lapply(names(preservation_results), function(ds) {
  if (is.null(preservation_results[[ds]])) return(NULL)
  data.frame(
    Dataset = ds,
    Module = preservation_results[[ds]]$modules,
    Zsummary = preservation_results[[ds]]$Zsummary,
    stringsAsFactors = FALSE
  )
}))

if (!is.null(zsummary_df) && nrow(zsummary_df) > 0) {
  # Thresholds: Z < 2 = no preservation, 2 <= Z < 10 = moderate, Z >= 10 = strong
  zsummary_df$Preservation <- cut(zsummary_df$Zsummary,
                                   breaks = c(-Inf, 2, 10, Inf),
                                   labels = c("No preservation", "Moderate", "Strong"))

  p_zsummary <- ggplot(zsummary_df, aes(Module, Zsummary, fill = Dataset)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = c(2, 10), linetype = "dashed", color = c("red", "green"), alpha = 0.5) +
    labs(title = "Module Preservation (Zsummary)",
         y = "Zsummary", x = "Module") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_tiff_600dpi(p_zsummary, file.path(OUT_DIR, "results/figures/Figure6A_Zsummary.tiff"),
                   width = 9, height = 6)
}

# Median rank plot
mr_df <- do.call(rbind, lapply(names(preservation_results), function(ds) {
  if (is.null(preservation_results[[ds]])) return(NULL)
  if (is.null(preservation_results[[ds]]$medianRank)) return(NULL)
  data.frame(
    Dataset = ds,
    Module = preservation_results[[ds]]$modules,
    medianRank = preservation_results[[ds]]$medianRank,
    stringsAsFactors = FALSE
  )
}))

if (!is.null(mr_df) && nrow(mr_df) > 0) {
  p_medianRank <- ggplot(mr_df, aes(Module, medianRank, color = Dataset, group = Dataset)) +
    geom_line(size = 1) + geom_point(size = 2) +
    labs(title = "Module Preservation (Median Rank)", y = "Median Rank", x = "Module") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_tiff_600dpi(p_medianRank, file.path(OUT_DIR, "results/figures/Figure6B_medianRank.tiff"),
                   width = 9, height = 6)
}

# Module overlap across datasets
mod_overlap <- matrix(0, nrow = length(unique(module_info$Module)),
                       ncol = length(names(discovery_list)))
rownames(mod_overlap) <- unique(module_info$Module)
colnames(mod_overlap) <- names(discovery_list)

for (ds in names(discovery_list)) {
  d <- discovery_list[[ds]]
  common <- intersect(rownames(d$expr), module_info$Gene)
  mod_common <- module_info[match(common, module_info$Gene), ]
  tbl <- table(mod_common$Module) / sum(table(mod_common$Module))
  mod_overlap[names(tbl), ds] <- tbl
}

open_tiff(file.path(OUT_DIR, "results/figures/Figure6C_module_overlap_heatmap.tiff"), width = 7, height = 6)
pheatmap::pheatmap(mod_overlap, display_numbers = TRUE, fontsize = 8,
                   color = colorRampPalette(c("white", "red"))(100),
                   main = "Module Overlap Across Datasets")
dev.off()
write.csv(mod_overlap, file.path(OUT_DIR, "results/module_preservation/module_overlap_matrix.csv"))

# Module size vs preservation
if (!is.null(zsummary_df) && nrow(zsummary_df) > 0) {
  module_sizes <- table(module_info$Module)
  zsummary_df$Module_Size <- module_sizes[zsummary_df$Module]

  p_size_vs_pres <- ggplot(zsummary_df, aes(Module_Size, Zsummary, color = Dataset)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_x_log10() +
    geom_hline(yintercept = c(2, 10), linetype = "dashed", alpha = 0.5) +
    labs(title = "Module Size vs Preservation",
         x = "Module Size (log scale)", y = "Zsummary") +
    theme_minimal()
  save_tiff_600dpi(p_size_vs_pres, file.path(OUT_DIR, "results/figures/Figure6D_module_size_preservation.tiff"))
}

saveRDS(preservation_results, file.path(OUT_DIR, "results/module_preservation/preservation_results.rds"))

# Try validation preservation if data available
val_data <- tryCatch(readRDS(file.path(OUT_DIR, "data_processed/validation/validation_data.rds")),
                      error = function(e) NULL)
if (!is.null(val_data)) {
  val_common <- intersect(rownames(val_data$expr), module_info$Gene)
  if (length(val_common) > 50) {
    log_msg(sprintf("Testing module preservation in validation cohort (%d genes)...", length(val_common)))
    val_mod <- module_info[match(val_common, module_info$Gene), ]
    multiExpr_val <- list(
      Reference = list(data = t(combat_expr[val_common, ])),
      Validation = list(data = t(val_data$expr[val_common, ]))
    )
    multiExpr_val <- lapply(multiExpr_val, function(set) {
      gsg <- goodSamplesGenes(set$data, verbose = 0)
      if (!gsg$allOK) {
        set$data <- set$data[gsg$goodSamples, gsg$goodGenes]
      }
      set
    })
    common_val_pres <- intersect(colnames(multiExpr_val$Reference$data), colnames(multiExpr_val$Validation$data))
    multiExpr_val$Reference$data <- multiExpr_val$Reference$data[, common_val_pres, drop = FALSE]
    multiExpr_val$Validation$data <- multiExpr_val$Validation$data[, common_val_pres, drop = FALSE]
    val_mod_pres <- val_mod[match(common_val_pres, val_mod$Gene), ]
    pres_val <- tryCatch({
      modulePreservation(multiExpr_val, multiColor = list(val_mod_pres$Module, val_mod_pres$Module),
                         nPermutations = 20, verbose = 0, networkType = "signed",
                         referenceNetworks = 1, maxGoldModuleSize = 100)
    }, error = function(e) NULL)

    if (!is.null(pres_val)) {
      z_val_scores <- pres_val$preservation$Z[[1]][[2]]
      z_val <- z_val_scores$Zsummary.pres
      mr_val <- if (!is.null(pres_val$preservation$medianRank)) pres_val$preservation$medianRank[[1]][[2]]$medianRank.pres else NULL
      pres_val_df <- data.frame(Module = rownames(z_val_scores), Zsummary = z_val, medianRank = if (is.null(mr_val)) NA else mr_val)
      write.csv(pres_val_df, file.path(OUT_DIR, "results/module_preservation/preservation_validation.csv"), row.names = FALSE)

      p_val_pres <- ggplot(pres_val_df, aes(reorder(Module, -Zsummary), Zsummary)) +
        geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
        geom_hline(yintercept = c(2, 10), linetype = "dashed", color = c("red", "green")) +
        labs(title = "Module Preservation in Validation Cohort",
             x = "Module", y = "Zsummary") +
        theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      save_tiff_600dpi(p_val_pres, file.path(OUT_DIR, "results/figures/Figure6F_validation_preservation.tiff"),
                       width = 9, height = 6)
    }
  }
}

log_msg("Step 15 complete: Module preservation analysis finished")
