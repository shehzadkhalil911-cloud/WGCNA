###############################################################################
# PRESERVATION FIX V3: Fix Z extraction from mp object, generate plots
###############################################################################
library(ggplot2)
options(stringsAsFactors = FALSE)

OUT_DIR <- "F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)"
FIG_DIR <- file.path(OUT_DIR, "figures/preservation")
TAB_DIR <- file.path(OUT_DIR, "tables")
RDA_DIR <- file.path(OUT_DIR, "RData")
dir.create(FIG_DIR, recursive=TRUE, showWarnings=FALSE)

load(file.path(RDA_DIR, "04_preservation.RData"))

cat("========== Fixing Preservation Z Extraction ==========\n")

fix_pres_obj <- function(pres_obj) {
  mp <- pres_obj$mp
  # Z statistics from Z$ref.ref[[2]]
  z_stats <- mp$preservation$Z[[1]][[2]]
  # Observed statistics from observed$ref.ref[[2]]
  obs_stats <- mp$preservation$observed[[1]][[2]]
  
  # Merge into single object with both Zsummary and medianRank
  z_merged <- z_stats
  z_merged$medianRank.pres <- obs_stats[, "medianRank.pres"]
  
  pres_obj$z <- z_merged
  pres_obj$ref <- obs_stats
  pres_obj$test <- mp$preservation$observed[[1]][[1]]
  pres_obj
}

pres_42568_to_15852 <- fix_pres_obj(pres_42568_to_15852)
pres_42568_to_29431 <- fix_pres_obj(pres_42568_to_29431)
pres_15852_to_29431 <- fix_pres_obj(pres_15852_to_29431)

# Verify
for (nm in c("pres_42568_to_15852", "pres_42568_to_29431", "pres_15852_to_29431")) {
  obj <- get(nm)
  z <- obj$z
  cat(sprintf("%s: Z dim=%dx%d, has Zsummary.pres: %s\n", 
              nm, nrow(z), ncol(z), "Zsummary.pres" %in% colnames(z)))
}

# Plot
plot_preservation <- function(pres_obj, ref_name, test_name) {
  z <- as.data.frame(pres_obj$z)
  if (is.null(z) || nrow(z) == 0) {
    cat(sprintf("  Skipping %s -> %s: empty Z\n", ref_name, test_name))
    return(NULL)
  }
  # Filter grey modules
  grey_idx <- grepl("grey", rownames(z), ignore.case = TRUE)
  if (sum(grey_idx) > 0) z <- z[!grey_idx, , drop = FALSE]
  if (nrow(z) < 1) return(NULL)
  
  df <- data.frame(Module = rownames(z),
                   Zsummary = z[, "Zsummary.pres"],
                   medianRank = z[, "medianRank.pres"],
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

cat("\n--- Generating preservation plots ---\n")
for (i in seq_along(pres_names)) {
  cat(sprintf("Plotting %s...\n", pres_names[i]))
  pp <- plot_preservation(pres_objects[[i]], readable_names[[i]][1], readable_names[[i]][2])
  if (!is.null(pp)) {
    ggsave(file.path(FIG_DIR, paste0("Zsummary_", pres_names[i], ".tiff")),
           pp$zsummary_plot, width = 10, height = 6, dpi = 300, compression = "lzw")
    ggsave(file.path(FIG_DIR, paste0("medianRank_", pres_names[i], ".tiff")),
           pp$medianrank_plot, width = 10, height = 6, dpi = 300, compression = "lzw")
    cat(sprintf("  Plots saved for %s\n", pres_names[i]))
  }
}

# Save tables
for (i in seq_along(pres_names)) {
  z_df <- as.data.frame(pres_objects[[i]]$z)
  z_df$Module <- rownames(z_df)
  write.csv(z_df, file.path(TAB_DIR, paste0("preservation_", pres_names[i], ".csv")), row.names = FALSE)
}

# Summary
cat("\n--- Preservation Summary ---\n")
for (i in seq_along(pres_names)) {
  z <- as.data.frame(pres_objects[[i]]$z)
  grey_idx <- grepl("grey", rownames(z), ignore.case = TRUE)
  if (sum(grey_idx) > 0) z <- z[!grey_idx, , drop = FALSE]
  if (nrow(z) > 0) {
    cat(sprintf("\n%s:\n", pres_names[i]))
    z_sorted <- z[order(-z[, "Zsummary.pres"]), ]
    cat(sprintf("  Top preserved: %s (Z=%.1f)\n", gsub("^ME", "", rownames(z_sorted)[1]), z_sorted[1, "Zsummary.pres"]))
    cat(sprintf("  Strongly preserved (Z>10): %d\n", sum(z[, "Zsummary.pres"] > 10, na.rm=TRUE)))
    cat(sprintf("  Moderately preserved (2<Z<=10): %d\n",
                sum(z[, "Zsummary.pres"] > 2 & z[, "Zsummary.pres"] <= 10, na.rm=TRUE)))
    cat(sprintf("  Not preserved (Z<=2): %d\n", sum(z[, "Zsummary.pres"] <= 2, na.rm=TRUE)))
  }
}

# Save fixed objects
save(pres_42568_to_15852, pres_42568_to_29431, pres_15852_to_29431,
     file = file.path(RDA_DIR, "04_preservation.RData"))

cat("\nPreservation fix complete.\n")
