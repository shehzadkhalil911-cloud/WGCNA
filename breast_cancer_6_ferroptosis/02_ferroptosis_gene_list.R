source("R/utils.R")
set.seed(123)

BASE_DIR <- "F:/independent projects/cancer/breast cancer"
OUT_DIR <- file.path(BASE_DIR, "Breast cancer 6 Ferroptosis")

log_msg("=== Step 2: Ferroptosis Gene List Preparation ===")

ferr_path <- file.path(BASE_DIR, "Ferroptosis Genes.xlsx")
ferr_raw <- openxlsx::read.xlsx(ferr_path, sheet = 1)

log_msg(sprintf("Raw ferroptosis file: %d rows, %d cols", nrow(ferr_raw), ncol(ferr_raw)))
log_msg(sprintf("Column names: %s", paste(colnames(ferr_raw), collapse = ", ")))

# Find Symbol column
colname_lower <- tolower(gsub("[^A-Za-z]", "", colnames(ferr_raw)))
symbol_idx <- which(colname_lower %in% c("symbol", "gene", "genesymbol", "genesymbol", "hugosymbol"))
if (length(symbol_idx) == 0) {
  # Check first column as default
  first_col <- colnames(ferr_raw)[1]
  if (grepl("symbol|gene", tolower(first_col))) {
    symbol_idx <- 1
  } else if (ncol(ferr_raw) >= 1) {
    symbol_idx <- 1
    log_msg(sprintf("WARNING: No 'Symbol' column found. Using first column: %s", first_col))
  }
}
if (length(symbol_idx) == 0) {
  stop("Cannot identify gene symbol column in ferroptosis gene list")
}

symbol_col <- colnames(ferr_raw)[symbol_idx[1]]
log_msg(sprintf("Using column: %s for gene symbols", symbol_col))

raw_genes <- as.character(ferr_raw[[symbol_col]])
raw_genes <- raw_genes[!is.na(raw_genes) & nchar(trimws(raw_genes)) > 0]
log_msg(sprintf("Raw genes extracted: %d", length(raw_genes)))

# Detect duplicates in raw list
dup_idx <- duplicated(raw_genes) | duplicated(raw_genes, fromLast = TRUE)
raw_dups <- raw_genes[dup_idx]

# Standardize using org.Hs.eg.db
library(org.Hs.eg.db)
library(AnnotationDbi)

genes_upper <- toupper(trimws(raw_genes))
keytypes <- keytypes(org.Hs.eg.db)
log_msg(sprintf("Checking %d genes against org.Hs.eg.db...", length(genes_upper)))

# Try mapping to official symbols
mapped <- tryCatch({
  suppressMessages(select(org.Hs.eg.db,
                          keys = genes_upper,
                          keytype = "SYMBOL",
                          columns = c("SYMBOL", "ENTREZID", "ENSEMBL", "GENENAME")))
}, error = function(e) {
  log_msg(sprintf("SYMBOL mapping failed: %s", e$message))
  # Try with ALIAS
  aliases <- tryCatch({
    suppressMessages(select(org.Hs.eg.db,
                            keys = genes_upper,
                            keytype = "ALIAS",
                            columns = c("SYMBOL", "ENTREZID", "ENSEMBL")))
  }, error = function(e2) {
    log_msg(sprintf("ALIAS mapping also failed: %s", e2$message))
    NULL
  })
  aliases
})

aliases_used <- character()
unmapped_genes <- character()
final_genes <- character()

if (!is.null(mapped) && nrow(mapped) > 0) {
  mapped <- mapped[!is.na(mapped$SYMBOL) & !duplicated(mapped$SYMBOL), ]

  # Check for alias-mapped genes
  direct_map <- intersect(genes_upper, mapped$SYMBOL)
  alias_map <- setdiff(mapped$SYMBOL, genes_upper)
  aliases_used <- mapped[mapped$SYMBOL %in% alias_map, ]

  # Also find unmapped
  all_mapped_symbols <- unique(mapped$SYMBOL)
  unmapped_genes <- setdiff(unique(genes_upper), all_mapped_symbols)
  final_genes <- all_mapped_symbols
} else {
  # Fallback: use original symbols
  log_msg("WARNING: Database mapping failed. Using raw gene symbols.")
  final_genes <- unique(genes_upper)
  unmapped_genes <- character()
}

# Remove known non-human artifacts
final_genes <- setdiff(final_genes, c("KIM-1", "GLUT13", "G6PDX", "MIR4715", "MIR30B",
                                       "MIR6852", "MIR137", "SNORA16A", "LOC284561",
                                       "LOC390705", "BLOC1S5-TXNDC5", "LINC00472",
                                       "YY1AP1", "LURAP1L"))

final_genes <- unique(final_genes[!is.na(final_genes) & nchar(final_genes) > 0])

# Get annotations for final genes
gene_anno <- tryCatch({
  suppressMessages(select(org.Hs.eg.db,
                          keys = final_genes,
                          keytype = "SYMBOL",
                          columns = c("ENTREZID", "ENSEMBL", "GENENAME")))
}, error = function(e) {
  data.frame(SYMBOL = final_genes, ENTREZID = NA, ENSEMBL = NA, GENENAME = NA,
             stringsAsFactors = FALSE)
})

gene_anno <- gene_anno[!duplicated(gene_anno$SYMBOL), ]
rownames(gene_anno) <- NULL

log_msg(sprintf("Final ferroptosis gene list: %d genes", length(final_genes)))
log_msg(sprintf("Unmapped genes: %d", length(unmapped_genes)))
if (length(aliases_used) > 0) {
  log_msg(sprintf("Alias-corrected genes: %d", nrow(aliases_used)))
}

# Save outputs
ferr_dir <- file.path(OUT_DIR, "data_processed/ferroptosis_gene_list")
writeLines(raw_genes, file.path(ferr_dir, "original_ferroptosis_genes.txt"))
writeLines(final_genes, file.path(ferr_dir, "final_ferroptosis_genes.txt"))
writeLines(unmapped_genes, file.path(ferr_dir, "unmapped_genes.txt"))
writeLines(raw_dups, file.path(ferr_dir, "duplicated_entries.txt"))

if (exists("aliases_used") && is.data.frame(aliases_used) && nrow(aliases_used) > 0) {
  write.csv(aliases_used, file.path(ferr_dir, "alias_corrected_genes.csv"), row.names = FALSE)
}

write.csv(gene_anno, file.path(ferr_dir, "ferroptosis_gene_annotations.csv"), row.names = FALSE)
writexl::write_xlsx(list(
  gene_annotations = gene_anno,
  original_list = data.frame(Symbol = raw_genes),
  unmapped = data.frame(Symbol = unmapped_genes)
), file.path(ferr_dir, "ferroptosis_gene_mapping.xlsx"))

# Also save for downstream use
saveRDS(final_genes, file.path(ferr_dir, "ferroptosis_genes_final.rds"))

log_msg("Step 2 complete: Ferroptosis gene list prepared")
