###############################################################################
# PRESERVATION FIX V2: ModulePreservation + Plot + Save at each step
###############################################################################
library(WGCNA)
library(ggplot2)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

OUT_DIR <- "F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)"
FIG_DIR <- file.path(OUT_DIR, "figures/preservation")
TAB_DIR <- file.path(OUT_DIR, "tables")
RDA_DIR <- file.path(OUT_DIR, "RData")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

load(file.path(RDA_DIR, "01_processed_data_common.RData"))
load(file.path(RDA_DIR, "03_WGCNA_GSE42568.RData"))
wgcna_42568 <- wgcna_result; rm(wgcna_result)
load(file.path(RDA_DIR, "03_WGCNA_GSE15852.RData"))
wgcna_15852 <- wgcna_result; rm(wgcna_result)
load(file.path(RDA_DIR, "03_WGCNA_GSE29431.RData"))
wgcna_29431 <- wgcna_result; rm(wgcna_result)
load(file.path(RDA_DIR, "04_module_matching.RData"))

cat("========== Module Preservation Analysis ==========\n")

run_preservation <- function(ref_expr, test_expr, ref_colors, ref_name, test_name) {
  cat(sprintf("Preservation: %s (ref) -> %s (test)\n", ref_name, test_name))
  multiExpr <- list(
    ref = list(data = t(ref_expr)),
    test = list(data = t(test_expr))
  )
  names(multiExpr$ref$data) <- NULL
  names(multiExpr$test$data) <- NULL

  old_wd <- getwd()
  setwd("C:/Users/shehz/AppData/Local/Temp/opencode/wgcna_pres")
  mp <- modulePreservation(multiExpr, list(ref = ref_colors),
                           referenceNetworks = 1,
                           nPermutations = 100,
                           networkType = "signed",
                           quickCor = 1,
                           verbose = 3)
  if (is.null(mp)) return(NULL)

  ref <- mp$preservation$observed[[1]][[1]]
  z_raw <- mp$preservation$Z[[1]][[1]]
  test <- NULL
  tryCatch({ test <- mp$preservation$observed[[1]][[2]] }, error = function(e) {})

  setwd(old_wd)
  list(ref = ref, test = test, z = z_raw, mp = mp)
}

ref_42568 <- comm_42568$expr
ref_15852 <- comm_15852$expr
ref_29431 <- comm_29431$expr

names(wgcna_42568$moduleColors) <- colnames(wgcna_42568$datExpr)
names(wgcna_15852$moduleColors) <- colnames(wgcna_15852$datExpr)
names(wgcna_29431$moduleColors) <- colnames(wgcna_29431$datExpr)

common42568_colors <- wgcna_42568$moduleColors[rownames(ref_42568)]
common42568_colors[is.na(common42568_colors)] <- "grey"
common15852_colors <- wgcna_15852$moduleColors[rownames(ref_15852)]
common15852_colors[is.na(common15852_colors)] <- "grey"
common29431_colors <- wgcna_29431$moduleColors[rownames(ref_29431)]
common29431_colors[is.na(common29431_colors)] <- "grey"

# Run and save each immediately
cat("\n--- Run 1/3: GSE42568 -> GSE15852 ---\n")
pres_42568_to_15852 <- run_preservation(ref_42568, ref_15852, common42568_colors, "GSE42568", "GSE15852")
cat("  Saving intermediate...\n")
save(pres_42568_to_15852, file = file.path(RDA_DIR, "04_pres_42568_to_15852.RData"))

cat("\n--- Run 2/3: GSE42568 -> GSE29431 ---\n")
pres_42568_to_29431 <- run_preservation(ref_42568, ref_29431, common42568_colors, "GSE42568", "GSE29431")
cat("  Saving intermediate...\n")
save(pres_42568_to_29431, file = file.path(RDA_DIR, "04_pres_42568_to_29431.RData"))

cat("\n--- Run 3/3: GSE15852 -> GSE29431 ---\n")
pres_15852_to_29431 <- run_preservation(ref_15852, ref_29431, common15852_colors, "GSE15852", "GSE29431")
cat("  Saving intermediate...\n")
save(pres_15852_to_29431, file = file.path(RDA_DIR, "04_pres_15852_to_29431.RData"))

# Now plot
plot_preservation <- function(pres_obj, ref_name, test_name) {
  if (is.null(pres_obj) || is.null(pres_obj$z)) {
    cat(sprintf("  Skipping %s -> %s: NULL result\n", ref_name, test_name))
    return(NULL)
  }
  z <- pres_obj$z
  # Ensure it's a matrix
  if (is.data.frame(z)) z <- as.matrix(z)
  if (is.null(dim(z))) {
    cat(sprintf("  z is not 2D for %s -> %s, attempting fix\n", ref_name, test_name))
    z <- t(as.matrix(z))
  }
  # Check column names
  if (!("Zsummary.pres" %in% colnames(z))) {
    cat(sprintf("  Column Zsummary.pres not found for %s -> %s\n", ref_name, test_name))
    cat("  Available columns:", paste(colnames(z), collapse=", "), "\n")
    return(NULL)
  }
  # Filter grey modules if rownames exist
  if (!is.null(rownames(z))) {
    grey_idx <- grepl("grey", rownames(z), ignore.case = TRUE)
    if (sum(grey_idx) > 0) {
      z <- z[!grey_idx, , drop = FALSE]
    }
  }
  if (nrow(z) < 1) {
    cat(sprintf("  No modules left after grey filtering for %s -> %s\n", ref_name, test_name))
    return(NULL)
  }

  df <- data.frame(Module = rownames(z),
                   Zsummary = as.numeric(z[, "Zsummary.pres"]),
                   medianRank = as.numeric(z[, "medianRank.pres"]),
                   stringsAsFactors = FALSE)
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
         subtitle = "Zsummary", x = "Module", y = "Zsummary") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p2 <- ggplot(df, aes(x = reorder(Module, medianRank), y = medianRank)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = paste("Module Preservation:", ref_name, "->", test_name),
         subtitle = "medianRank (lower = better preserved)",
         x = "Module", y = "medianRank") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  list(zsummary_plot = p1, medianrank_plot = p2)
}

pres_names <- c("pres_42568_to_15852", "pres_42568_to_29431", "pres_15852_to_29431")
pres_objects <- list(pres_42568_to_15852, pres_42568_to_29431, pres_15852_to_29431)
readable_names <- list(c("GSE42568", "GSE15852"), c("GSE42568", "GSE29431"), c("GSE15852", "GSE29431"))

cat("\n--- Generating plots ---\n")
pres_plots <- list()
for (i in seq_along(pres_names)) {
  cat(sprintf("Plotting %s...\n", pres_names[i]))
  pp <- plot_preservation(pres_objects[[i]], readable_names[[i]][1], readable_names[[i]][2])
  if (!is.null(pp)) {
    ggsave(file.path(FIG_DIR, paste0("Zsummary_", pres_names[i], ".tiff")),
           pp$zsummary_plot, width = 10, height = 6, dpi = 600, compression = "lzw")
    ggsave(file.path(FIG_DIR, paste0("medianRank_", pres_names[i], ".tiff")),
           pp$medianrank_plot, width = 10, height = 6, dpi = 600, compression = "lzw")
  }
  # Save intermediate plot objects too
  pres_plots[[i]] <- pp
}

# Save tables
cat("\n--- Saving tables ---\n")
for (i in seq_along(pres_names)) {
  z_df <- as.data.frame(as.matrix(pres_objects[[i]]$z))
  z_df$Module <- rownames(pres_objects[[i]]$z)
  write.csv(z_df, file.path(TAB_DIR, paste0("preservation_", pres_names[i], ".csv")), row.names = FALSE)
}

# Summary
cat("\n--- Preservation Summary (Zsummary) ---\n")
for (i in seq_along(pres_names)) {
  z <- as.matrix(pres_objects[[i]]$z)
  if (!is.null(rownames(z))) {
    grey_idx <- grepl("grey", rownames(z), ignore.case = TRUE)
    if (sum(grey_idx) > 0) z <- z[!grey_idx, , drop = FALSE]
  }
  if (nrow(z) > 0 && "Zsummary.pres" %in% colnames(z)) {
    cat(sprintf("\n%s:\n", pres_names[i]))
    z_sorted <- z[order(-as.numeric(z[, "Zsummary.pres"])), , drop = FALSE]
    cat(sprintf("  Top preserved: %s (Z=%.1f)\n",
                gsub("ME", "", rownames(z_sorted)[1]), as.numeric(z_sorted[1, "Zsummary.pres"])))
    cat(sprintf("  Strongly preserved (Z>10): %d\n", sum(as.numeric(z[, "Zsummary.pres"]) > 10, na.rm=TRUE)))
    cat(sprintf("  Moderately preserved (2<Z<=10): %d\n",
                sum(as.numeric(z[, "Zsummary.pres"]) > 2 & as.numeric(z[, "Zsummary.pres"]) <= 10, na.rm=TRUE)))
  }
}

# Final save
save(pres_42568_to_15852, pres_42568_to_29431, pres_15852_to_29431,
     file = file.path(RDA_DIR, "04_preservation.RData"))

# Clean up intermediate files
file.remove(file.path(RDA_DIR, "04_pres_42568_to_15852.RData"))
file.remove(file.path(RDA_DIR, "04_pres_42568_to_29431.RData"))
file.remove(file.path(RDA_DIR, "04_pres_15852_to_29431.RData"))

cat("\nPreservation analysis complete.\n")
