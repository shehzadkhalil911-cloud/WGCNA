# =============================================================================
# Script 01: Input Inventory and Gene List Preparation
# Breast Cancer 7 Cell Death Crosstalk
# =============================================================================

source(file.path("..", "config", "analysis_config.R"))

library(readxl)
library(openxlsx)
library(tidyverse)
library(data.table)
library(janitor)
library(stringr)
library(org.Hs.eg.db)

log_message("Starting input inventory and gene-list preparation.")

# --- 1. RECURSIVE FILE INVENTORY ---
inventory_files <- function(dir_path) {
  all_files <- list.files(dir_path, recursive = TRUE, full.names = TRUE)
  file_info <- data.frame(
    full_path = all_files,
    file_name = basename(all_files),
    extension = tools::file_ext(all_files),
    file_size_bytes = file.info(all_files)$size,
    stringsAsFactors = FALSE
  )
  file_info$inferred_dataset <- NA_character_
  file_info$inferred_purpose <- NA_character_
  file_info$selected <- TRUE
  file_info$reason <- NA_character_

  # Classify files
  for (i in seq_len(nrow(file_info))) {
    fn <- file_info$file_name[i]
    fp <- file_info$full_path[i]

    if (grepl("GSE42568_expression", fn)) {
      file_info$inferred_dataset[i] <- "GSE42568"
      file_info$inferred_purpose[i] <- "Discovery expression matrix"
    } else if (grepl("GSE15852_expression", fn)) {
      file_info$inferred_dataset[i] <- "GSE15852"
      file_info$inferred_purpose[i] <- "Discovery expression matrix"
    } else if (grepl("GSE29431_expression", fn)) {
      file_info$inferred_dataset[i] <- "GSE29431"
      file_info$inferred_purpose[i] <- "Discovery expression matrix"
    } else if (grepl("GSE70947_expression", fn)) {
      file_info$inferred_dataset[i] <- "GSE70947"
      file_info$inferred_purpose[i] <- "Validation expression matrix"
    } else if (grepl("GSE42568_processed", fn)) {
      file_info$inferred_dataset[i] <- "GSE42568"
      file_info$inferred_purpose[i] <- "Discovery processed R object"
    } else if (grepl("GSE15852_processed", fn)) {
      file_info$inferred_dataset[i] <- "GSE15852"
      file_info$inferred_purpose[i] <- "Discovery processed R object"
    } else if (grepl("GSE29431_processed", fn)) {
      file_info$inferred_dataset[i] <- "GSE29431"
      file_info$inferred_purpose[i] <- "Discovery processed R object"
    } else if (grepl("GSE70947_processed", fn)) {
      file_info$inferred_dataset[i] <- "GSE70947"
      file_info$inferred_purpose[i] <- "Validation processed R object"
    } else if (grepl("GSE42568_eset", fn)) {
      file_info$inferred_dataset[i] <- "GSE42568"
      file_info$inferred_purpose[i] <- "Raw ExpressionSet"
    } else if (grepl("GSE15852_eset", fn)) {
      file_info$inferred_dataset[i] <- "GSE15852"
      file_info$inferred_purpose[i] <- "Raw ExpressionSet"
    } else if (grepl("GSE29431_eset", fn)) {
      file_info$inferred_dataset[i] <- "GSE29431"
      file_info$inferred_purpose[i] <- "Raw ExpressionSet"
    } else if (grepl("GSE70947_eset", fn)) {
      file_info$inferred_dataset[i] <- "GSE70947"
      file_info$inferred_purpose[i] <- "Validation ExpressionSet"
    } else if (grepl("GSE42568_sample_classes", fn)) {
      file_info$inferred_dataset[i] <- "GSE42568"
      file_info$inferred_purpose[i] <- "Sample class assignment"
    } else if (grepl("GSE15852_sample_classes", fn)) {
      file_info$inferred_dataset[i] <- "GSE15852"
      file_info$inferred_purpose[i] <- "Sample class assignment"
    } else if (grepl("GSE29431_sample_classes", fn)) {
      file_info$inferred_dataset[i] <- "GSE29431"
      file_info$inferred_purpose[i] <- "Sample class assignment"
    } else if (grepl("GSE70947_sample_classes", fn)) {
      file_info$inferred_dataset[i] <- "GSE70947"
      file_info$inferred_purpose[i] <- "Validation sample class assignment"
    } else if (grepl("GPL570\\.soft", fn)) {
      file_info$inferred_dataset[i] <- "GSE42568/GSE15852"
      file_info$inferred_purpose[i] <- "Platform annotation (GPL570)"
    } else if (grepl("GPL96\\.soft", fn)) {
      file_info$inferred_dataset[i] <- "GSE29431"
      file_info$inferred_purpose[i] <- "Platform annotation (GPL96)"
    } else if (grepl("GPL13607\\.soft", fn)) {
      file_info$inferred_dataset[i] <- "GSE70947"
      file_info$inferred_purpose[i] <- "Platform annotation (GPL13607)"
    } else if (grepl("discovery_cohort", fn)) {
      file_info$inferred_purpose[i] <- "Pre-integrated discovery cohort"
    } else if (grepl("discovery_combat", fn)) {
      file_info$inferred_purpose[i] <- "Batch-corrected discovery cohort"
    } else if (grepl("validation_cohort", fn)) {
      file_info$inferred_purpose[i] <- "Pre-integrated validation cohort"
    } else if (grepl("dataset_inspection", fn)) {
      file_info$inferred_purpose[i] <- "Previous dataset inspection summary"
    }
  }

  return(file_info)
}

manifest <- inventory_files(DATA_DIR)
write.csv(manifest, file.path(OUTPUT$inventory, "input_file_manifest.csv"), row.names = FALSE)
saveRDS(manifest, file.path(OUTPUT$inventory, "input_file_manifest.rds"))

log_message(sprintf("Input manifest created: %d files catalogued.", nrow(manifest)))

# --- 2. GENE LIST READING ---

# 2a. Cuproptosis gene list
crg_df <- read.csv(CUPROPTOSIS_FILE, stringsAsFactors = FALSE)
crg_genes <- crg_df$GeneSymbol
crg_genes <- crg_genes[!is.na(crg_genes) & crg_genes != ""]
crg_genes <- toupper(trimws(crg_genes))
crg_genes <- unique(crg_genes)
log_message(sprintf("Cuproptosis: %d unique genes from %s", length(crg_genes), basename(CUPROPTOSIS_FILE)))

# 2b. Ferroptosis gene list
ferr_sheets <- tryCatch(excel_sheets(FERROPTOSIS_FILE), error = function(e) NULL)
log_message(sprintf("Ferroptosis sheets: %s", paste(ferr_sheets, collapse = ", ")))

ferr_df <- read_excel(FERROPTOSIS_FILE, sheet = 1)
log_message(sprintf("Ferroptosis columns: %s", paste(colnames(ferr_df), collapse = ", ")))
log_message(sprintf("Ferroptosis first few rows:"))
print(head(ferr_df, 10))

# Try to find gene column
ferr_col <- grep("gene|symbol|name|gene_symbol", colnames(ferr_df), ignore.case = TRUE, value = TRUE)
if (length(ferr_col) == 0) ferr_col <- colnames(ferr_df)[1]
log_message(sprintf("Using column: %s", ferr_col[1]))

ferr_genes <- ferr_df[[ferr_col[1]]]
ferr_genes <- ferr_genes[!is.na(ferr_genes) & ferr_genes != ""]
ferr_genes <- toupper(trimws(as.character(ferr_genes)))
ferr_genes <- unique(ferr_genes)
log_message(sprintf("Ferroptosis: %d unique genes from %s", length(ferr_genes), basename(FERROPTOSIS_FILE)))

# 2c. Necroptosis gene list
necr_sheets <- tryCatch(excel_sheets(NECROPTOSIS_FILE), error = function(e) NULL)
log_message(sprintf("Necroptosis sheets: %s", paste(necr_sheets, collapse = ", ")))

necr_df <- read_excel(NECROPTOSIS_FILE, sheet = 1)
log_message(sprintf("Necroptosis columns: %s", paste(colnames(necr_df), collapse = ", ")))
log_message(sprintf("Necroptosis first few rows:"))
print(head(necr_df, 10))

necr_col <- grep("gene|symbol|name|gene_symbol", colnames(necr_df), ignore.case = TRUE, value = TRUE)
if (length(necr_col) == 0) necr_col <- colnames(necr_df)[1]
log_message(sprintf("Using column: %s", necr_col[1]))

necr_genes <- necr_df[[necr_col[1]]]
necr_genes <- necr_genes[!is.na(necr_genes) & necr_genes != ""]
necr_genes <- toupper(trimws(as.character(necr_genes)))
necr_genes <- unique(necr_genes)
log_message(sprintf("Necroptosis: %d unique genes from %s", length(necr_genes), basename(NECROPTOSIS_FILE)))

# --- 3. MAP TO OFFICIAL SYMBOLS ---
map_to_official <- function(gene_vec, name) {
  # Query org.Hs.eg.db for official symbols
  mapped <- tryCatch({
    mapIds(org.Hs.eg.db, keys = gene_vec, keytype = "SYMBOL", column = "SYMBOL", multiVals = "first")
  }, error = function(e) {
    log_message(sprintf("Mapping issue for %s: %s", name, e$message), "WARNING")
    return(setNames(gene_vec, gene_vec))
  })
  
  result_df <- data.frame(
    input_gene = gene_vec,
    official_symbol = ifelse(is.na(mapped), gene_vec, mapped),
    mapped = !is.na(mapped),
    stringsAsFactors = FALSE
  )
  
  unmapped <- result_df$input_gene[!result_df$mapped]
  if (length(unmapped) > 0) {
    log_message(sprintf("%s: %d unmapped genes: %s", name, length(unmapped), 
                        paste(head(unmapped, 20), collapse = ", ")), "WARNING")
    write.csv(data.frame(unmapped_genes = unmapped),
              file.path(OUTPUT$inventory, sprintf("%s_unmapped_genes.csv", name)),
              row.names = FALSE)
  }
  
  return(result_df)
}

ferr_mapped <- map_to_official(ferr_genes, "ferroptosis")
necr_mapped <- map_to_official(necr_genes, "necroptosis")
crg_mapped <- map_to_official(crg_genes, "cuproptosis")

# Store standardized gene lists
ferr_genes_official <- unique(ferr_mapped$official_symbol)
necr_genes_official <- unique(necr_mapped$official_symbol)
crg_genes_official <- unique(crg_mapped$official_symbol)

# --- 4. SAVE GENE LISTS ---
gene_lists <- list(
  ferroptosis = ferr_genes_official,
  necroptosis = necr_genes_official,
  cuproptosis = crg_genes_official
)

saveRDS(gene_lists, file.path(OUTPUT$inventory, "gene_lists_raw.rds"))

# Save detailed gene list tables
write.csv(data.frame(GeneSymbol = ferr_genes_official), 
          file.path(OUTPUT$inventory, "ferroptosis_genes_standardized.csv"), row.names = FALSE)
write.csv(data.frame(GeneSymbol = necr_genes_official), 
          file.path(OUTPUT$inventory, "necroptosis_genes_standardized.csv"), row.names = FALSE)
write.csv(data.frame(GeneSymbol = crg_genes_official), 
          file.path(OUTPUT$inventory, "cuproptosis_genes_standardized.csv"), row.names = FALSE)

log_message(sprintf("Final gene counts - Ferroptosis: %d, Necroptosis: %d, Cuproptosis: %d",
                    length(ferr_genes_official), length(necr_genes_official), length(crg_genes_official)))

# --- 5. GENE LIST OVERLAP ANALYSIS ---
library(ggvenn)
library(ComplexHeatmap)
library(circlize)

# Venn diagram of three lists
venn_data <- list(
  Ferroptosis = ferr_genes_official,
  Necroptosis = necr_genes_official,
  Cuproptosis = crg_genes_official
)

p_venn <- ggvenn(venn_data, 
                 fill_color = c("#E69F00", "#56B4E9", "#009E73"),
                 stroke_size = 0.5, 
                 set_name_size = 4,
                 show_percentage = TRUE) +
  ggtitle("Overlap of Ferroptosis, Necroptosis, and Cuproptosis Gene Lists") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

save_tiff_600dpi(p_venn, file.path(OUTPUT$inventory, "gene_list_venn.tiff"))
ggsave(file.path(OUTPUT$inventory, "gene_list_venn.pdf"), plot = p_venn, width = 7, height = 6)

# Detailed overlap
ferr_set <- ferr_genes_official
necr_set <- necr_genes_official
crg_set <- crg_genes_official

all_genes <- unique(c(ferr_set, necr_set, crg_set))
overlap_matrix <- data.frame(
  Gene = all_genes,
  Ferroptosis = all_genes %in% ferr_set,
  Necroptosis = all_genes %in% necr_set,
  Cuproptosis = all_genes %in% crg_set,
  stringsAsFactors = FALSE
)

overlap_matrix$Category <- case_when(
  overlap_matrix$Ferroptosis & overlap_matrix$Necroptosis & overlap_matrix$Cuproptosis ~ "All Three",
  overlap_matrix$Ferroptosis & overlap_matrix$Necroptosis ~ "Ferroptosis + Necroptosis",
  overlap_matrix$Ferroptosis & overlap_matrix$Cuproptosis ~ "Ferroptosis + Cuproptosis",
  overlap_matrix$Necroptosis & overlap_matrix$Cuproptosis ~ "Necroptosis + Cuproptosis",
  overlap_matrix$Ferroptosis ~ "Ferroptosis Only",
  overlap_matrix$Necroptosis ~ "Necroptosis Only",
  overlap_matrix$Cuproptosis ~ "Cuproptosis Only"
)

write.csv(overlap_matrix, file.path(OUTPUT$inventory, "gene_list_overlap_matrix.csv"), row.names = FALSE)

# Shared genes
ferr_necr_crg <- intersect(intersect(ferr_set, necr_set), crg_set)
ferr_necr <- setdiff(intersect(ferr_set, necr_set), ferr_necr_crg)
ferr_crg <- setdiff(intersect(ferr_set, crg_set), ferr_necr_crg)
necr_crg <- setdiff(intersect(necr_set, crg_set), ferr_necr_crg)
ferr_only <- setdiff(ferr_set, union(necr_set, crg_set))
necr_only <- setdiff(necr_set, union(ferr_set, crg_set))
crg_only <- setdiff(crg_set, union(ferr_set, necr_set))

category_genes <- list(
  All_Three = ferr_necr_crg,
  Ferroptosis_Necroptosis = ferr_necr,
  Ferroptosis_Cuproptosis = ferr_crg,
  Necroptosis_Cuproptosis = necr_crg,
  Ferroptosis_Only = ferr_only,
  Necroptosis_Only = necr_only,
  Cuproptosis_Only = crg_only
)

for (cat_name in names(category_genes)) {
  genes <- category_genes[[cat_name]]
  if (length(genes) > 0) {
    write.csv(data.frame(GeneSymbol = genes), 
              file.path(OUTPUT$inventory, sprintf("gene_list_%s.csv", cat_name)), 
              row.names = FALSE)
  }
}

log_message(sprintf("Gene overlap: All Three=%d, Ferro+Necro=%d, Ferro+Cupro=%d, Necro+Cupro=%d",
                    length(ferr_necr_crg), length(ferr_necr), length(ferr_crg), length(necr_crg)))
log_message(sprintf("Unique: Ferroptosis=%d, Necroptosis=%d, Cuproptosis=%d",
                    length(ferr_only), length(necr_only), length(crg_only)))

# UpSet plot
upset_data <- overlap_matrix[, c("Ferroptosis", "Necroptosis", "Cuproptosis")]
upset_df <- as.data.frame(lapply(upset_data, as.integer))
rownames(upset_df) <- overlap_matrix$Gene

lib_upset <- tryCatch(find.package("UpSetR"), error = function(e) NULL)
if (!is.null(lib_upset)) {
  library(UpSetR)
  open_tiff(file.path(OUTPUT$inventory, "gene_list_upset.tiff"), width = 8, height = 6)
  upset(
    upset_df,
    sets = c("Ferroptosis", "Necroptosis", "Cuproptosis"),
    order.by = "freq",
    keep.order = TRUE,
    mainbar.y.label = "Gene Count",
    sets.x.label = "Total Genes",
    text.scale = c(1.3, 1.3, 1, 1, 1.5, 1)
  )
  dev.off()
}

# Gene count comparison bar plot
count_df <- data.frame(
  Program = c("Ferroptosis", "Necroptosis", "Cuproptosis"),
  Count = c(length(ferr_set), length(necr_set), length(crg_set)),
  stringsAsFactors = FALSE
)

p_counts <- ggplot(count_df, aes(x = Program, y = Count, fill = Program)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  geom_text(aes(label = Count), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Ferroptosis" = "#E69F00", "Necroptosis" = "#56B4E9", "Cuproptosis" = "#009E73")) +
  labs(title = "Regulated Cell Death Gene Lists", 
       subtitle = "Unique gene count per program", y = "Number of Genes") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "none")

save_tiff_600dpi(p_counts, file.path(OUTPUT$inventory, "gene_list_counts.tiff"))

# Overlap heatmap
overlap_heatmap_data <- matrix(0, nrow = 3, ncol = 3)
rownames(overlap_heatmap_data) <- c("Ferroptosis", "Necroptosis", "Cuproptosis")
colnames(overlap_heatmap_data) <- c("Ferroptosis", "Necroptosis", "Cuproptosis")
overlap_heatmap_data[1,1] <- length(ferr_set)
overlap_heatmap_data[2,2] <- length(necr_set)
overlap_heatmap_data[3,3] <- length(crg_set)
overlap_heatmap_data[1,2] <- length(intersect(ferr_set, necr_set))
overlap_heatmap_data[1,3] <- length(intersect(ferr_set, crg_set))
overlap_heatmap_data[2,1] <- length(intersect(ferr_set, necr_set))
overlap_heatmap_data[2,3] <- length(intersect(necr_set, crg_set))
overlap_heatmap_data[3,1] <- length(intersect(ferr_set, crg_set))
overlap_heatmap_data[3,2] <- length(intersect(necr_set, crg_set))

ht <- Heatmap(overlap_heatmap_data,
              name = "Overlap",
              col = colorRamp2(c(0, max(overlap_heatmap_data)), c("white", "#0072B2")),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(overlap_heatmap_data[i, j], x, y, gp = gpar(fontsize = 12, fontface = "bold"))
              },
              rect_gp = gpar(col = "grey60", lwd = 1.5),
              column_title = "Gene List Overlap",
              row_title = "Program",
              row_names_side = "left",
              column_names_side = "top",
              cluster_rows = FALSE,
              cluster_columns = FALSE)

open_tiff(file.path(OUTPUT$inventory, "gene_list_overlap_heatmap.tiff"), width = 6, height = 5)
draw(ht)
dev.off()

pdf(file.path(OUTPUT$inventory, "gene_list_overlap_heatmap.pdf"), width = 6, height = 5)
draw(ht)
dev.off()

# Map to Entrez IDs
map_to_entrez <- function(gene_vec, name) {
  entrez <- tryCatch({
    mapIds(org.Hs.eg.db, keys = gene_vec, keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
  }, error = function(e) {
    log_message(sprintf("Entrez mapping issue for %s: %s", name, e$message), "WARNING")
    return(setNames(rep(NA, length(gene_vec)), gene_vec))
  })
  return(entrez)
}

ferr_entrez <- map_to_entrez(ferr_genes_official, "ferroptosis")
necr_entrez <- map_to_entrez(necr_genes_official, "necroptosis")
crg_entrez <- map_to_entrez(crg_genes_official, "cuproptosis")

entrez_mapping <- data.frame(
  GeneSymbol = c(ferr_genes_official, necr_genes_official, crg_genes_official),
  EntrezID = c(ferr_entrez, necr_entrez, crg_entrez),
  Program = c(rep("Ferroptosis", length(ferr_genes_official)),
              rep("Necroptosis", length(necr_genes_official)),
              rep("Cuproptosis", length(crg_genes_official))),
  stringsAsFactors = FALSE
)

write.csv(entrez_mapping, file.path(OUTPUT$inventory, "gene_list_entrez_mapping.csv"), row.names = FALSE)

log_message("Input inventory and gene-list preparation complete.")
