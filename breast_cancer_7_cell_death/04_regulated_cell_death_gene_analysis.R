# =============================================================================
# Script 04: Regulated Cell Death Gene Analysis (FRG, NRG, CRG)
# Breast Cancer 7 Cell Death Crosstalk
# =============================================================================

source(file.path("..", "config", "analysis_config.R"))

library(tidyverse)
library(ggvenn)
library(ComplexHeatmap)
library(circlize)

log_message("Starting regulated cell-death gene analysis.")

# --- 1. LOAD DATA ---
gene_lists <- readRDS(file.path(OUTPUT$inventory, "gene_lists_raw.rds"))
de_results <- readRDS(file.path(OUTPUT$de, "discovery_DEG_full.rds"))
deg_sig <- read.csv(file.path(OUTPUT$de, "discovery_DEG_significant.csv"), row.names = 1)
combat_matrix <- readRDS(file.path(OUTPUT$batch_correction, "combat_corrected_matrix.rds"))
discovery_metadata <- readRDS(file.path(OUTPUT$batch_correction, "discovery_metadata.rds"))

ferr_genes <- gene_lists$ferroptosis
necr_genes <- gene_lists$necroptosis
crg_genes <- gene_lists$cuproptosis

log_message(sprintf("Loaded: %d FRGs, %d NRGs, %d CRGs", 
                    length(ferr_genes), length(necr_genes), length(crg_genes)))

# --- 2. IDENTIFY DETECTED AND DE GENES ---
all_genes_in_de <- rownames(de_results)

analyze_program <- function(program_genes, program_name, output_dir) {
  detected <- intersect(program_genes, all_genes_in_de)
  undetected <- setdiff(program_genes, all_genes_in_de)
  
  de_idx <- which(rownames(de_results) %in% program_genes & 
                    de_results$adj.P.Val < ADJ_PVAL_THRESHOLD &
                    abs(de_results$logFC) >= LOG2FC_THRESHOLD)
  de_genes <- rownames(de_results)[de_idx]
  nonsig_detected <- setdiff(detected, de_genes)
  
  de_up <- de_genes[de_results[de_idx, "logFC"] > 0]
  de_down <- de_genes[de_results[de_idx, "logFC"] < 0]
  
  # Save
  write.csv(data.frame(Gene = detected), 
            file.path(output_dir, "detected_genes.csv"), row.names = FALSE)
  write.csv(data.frame(Gene = undetected), 
            file.path(output_dir, "undetected_genes.csv"), row.names = FALSE)
  write.csv(data.frame(Gene = de_genes), 
            file.path(output_dir, "DE_genes.csv"), row.names = FALSE)
  write.csv(data.frame(Gene = de_up), 
            file.path(output_dir, "DE_genes_upregulated.csv"), row.names = FALSE)
  write.csv(data.frame(Gene = de_down), 
            file.path(output_dir, "DE_genes_downregulated.csv"), row.names = FALSE)
  write.csv(data.frame(Gene = nonsig_detected), 
            file.path(output_dir, "detected_nonsignificant.csv"), row.names = FALSE)
  
  log_message(sprintf("%s: Detected=%d, DE=%d, Up=%d, Down=%d, Non-sig=%d, Undetected=%d",
                      program_name, length(detected), length(de_genes), length(de_up), 
                      length(de_down), length(nonsig_detected), length(undetected)))
  
  # Subset DE results
  if (length(de_genes) > 0) {
    de_subset <- de_results[de_idx, ]
    write.csv(de_subset, file.path(output_dir, "DE_gene_details.csv"))
  }
  
  return(list(
    detected = detected,
    undetected = undetected,
    de_genes = de_genes,
    de_up = de_up,
    de_down = de_down,
    nonsig = nonsig_detected
  ))
}

ferr_analysis <- analyze_program(ferr_genes, "Ferroptosis", OUTPUT$ferroptosis)
necr_analysis <- analyze_program(necr_genes, "Necroptosis", OUTPUT$necroptosis)
crg_analysis <- analyze_program(crg_genes, "Cuproptosis", OUTPUT$cuproptosis)

# --- 3. VENN DIAGRAMS: DEG vs Program Genes ---
for (prog_name in c("Ferroptosis", "Necroptosis", "Cuproptosis")) {
  analysis <- list(ferr_analysis, necr_analysis, crg_analysis)[[which(c("Ferroptosis","Necroptosis","Cuproptosis") == prog_name)]]
  
  venn_data <- list(
    DEGs = rownames(deg_sig),
    Program = analysis$detected
  )
  
  p_venn <- ggvenn(venn_data,
                   fill_color = c("#D55E00", "#009E73"),
                   stroke_size = 0.5, set_name_size = 4,
                   show_percentage = TRUE) +
    ggtitle(sprintf("%s: DEG ∩ Program Genes", prog_name)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  save_tiff_600dpi(p_venn, file.path(OUTPUT$shared, sprintf("%s_DEG_overlap_venn.tiff", tolower(prog_name))))
}

# --- 4. THREE-WAY VENN: DE-FRGs, DE-NRGs, DE-CRGs ---
three_venn <- list(
  `DE-FRG` = ferr_analysis$de_genes,
  `DE-NRG` = necr_analysis$de_genes,
  `DE-CRG` = crg_analysis$de_genes
)

p_three_venn <- ggvenn(three_venn,
                       fill_color = c("#E69F00", "#56B4E9", "#009E73"),
                       stroke_size = 0.5, set_name_size = 4,
                       show_percentage = TRUE) +
  ggtitle("Overlap: DE-FRGs, DE-NRGs, DE-CRGs") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_tiff_600dpi(p_three_venn, file.path(OUTPUT$shared, "three_program_DE_venn.tiff"))

# --- 5. UPSET PLOT ---
de_frg_set <- ferr_analysis$de_genes
de_nrg_set <- necr_analysis$de_genes
de_crg_set <- crg_analysis$de_genes

all_de_genes <- unique(c(de_frg_set, de_nrg_set, de_crg_set))
upset_df <- data.frame(
  row.names = all_de_genes,
  DE_FRG = as.integer(all_de_genes %in% de_frg_set),
  DE_NRG = as.integer(all_de_genes %in% de_nrg_set),
  DE_CRG = as.integer(all_de_genes %in% de_crg_set)
)

lib_upset <- tryCatch(find.package("UpSetR"), error = function(e) NULL)
if (!is.null(lib_upset)) {
  library(UpSetR)
  open_tiff(file.path(OUTPUT$shared, "three_program_upset.tiff"), width = 8, height = 6)
  upset(upset_df,
        sets = c("DE_FRG", "DE_NRG", "DE_CRG"),
        order.by = "freq",
        keep.order = TRUE,
        mainbar.y.label = "Gene Count",
        sets.x.label = "Total DE Genes",
        text.scale = c(1.3, 1.3, 1, 1, 1.5, 1))
  dev.off()
}

# --- 6. COMBINED VOLCANO PLOT ---
de_results$Program <- "Other"
de_results$Program[rownames(de_results) %in% de_frg_set] <- "DE-FRG"
de_results$Program[rownames(de_results) %in% de_nrg_set] <- "DE-NRG"
de_results$Program[rownames(de_results) %in% de_crg_set] <- "DE-CRG"

# Three-programme volcano
p_volcano_prog <- ggplot(de_results, aes(x = logFC, y = -log10(adj.P.Val), color = Program)) +
  geom_point(data = subset(de_results, Program == "Other"), size = 0.5, alpha = 0.3, color = "grey70") +
  geom_point(data = subset(de_results, Program != "Other"), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("DE-FRG" = "#E69F00", "DE-NRG" = "#56B4E9", 
                                "DE-CRG" = "#009E73", "Other" = "grey70")) +
  geom_hline(yintercept = -log10(ADJ_PVAL_THRESHOLD), linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD), linetype = "dashed", color = "red", linewidth = 0.5) +
  labs(title = "Regulated Cell Death Genes in DEG Landscape",
       x = "log2 Fold Change", y = "-log10(Adjusted P-value)") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
save_tiff_600dpi(p_volcano_prog, file.path(OUTPUT$shared, "combined_volcano_program.tiff"))

# Programme-specific volcanos
for (prog in c("DE-FRG", "DE-NRG", "DE-CRG")) {
  de_results$Highlight <- ifelse(de_results$Program == prog, "Yes", "No")
  
  p_spec <- ggplot(de_results, aes(x = logFC, y = -log10(adj.P.Val), color = Highlight)) +
    geom_point(data = subset(de_results, Highlight == "No"), size = 0.5, alpha = 0.3, color = "grey70") +
    geom_point(data = subset(de_results, Highlight == "Yes"), size = 2.5, alpha = 0.9, color = 
                 switch(prog, "DE-FRG" = "#E69F00", "DE-NRG" = "#56B4E9", "DE-CRG" = "#009E73")) +
    geom_hline(yintercept = -log10(ADJ_PVAL_THRESHOLD), linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_vline(xintercept = c(-LOG2FC_THRESHOLD, LOG2FC_THRESHOLD), linetype = "dashed", color = "red", linewidth = 0.5) +
    labs(title = sprintf("%s Highlighted", prog), x = "log2 FC", y = "-log10(adj.P)") +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_tiff_600dpi(p_spec, file.path(OUTPUT$shared, sprintf("volcano_%s.tiff", gsub("-", "_", tolower(prog)))))
}

# --- 7. THREE-PROGRAMME EXPRESSION HEATMAP ---
all_prog_de <- unique(c(de_frg_set, de_nrg_set, de_crg_set))

if (length(all_prog_de) > 2) {
  hm_data <- combat_matrix[intersect(all_prog_de, rownames(combat_matrix)), ]
  hm_data <- t(scale(t(hm_data)))
  
  prog_ann <- rowAnnotation(
    Program = ifelse(rownames(hm_data) %in% de_frg_set, "DE-FRG", 
              ifelse(rownames(hm_data) %in% de_nrg_set, "DE-NRG", "DE-CRG")),
    col = list(Program = c("DE-FRG" = "#E69F00", "DE-NRG" = "#56B4E9", "DE-CRG" = "#009E73"))
  )
  
  col_ann <- HeatmapAnnotation(
    Group = discovery_metadata$group,
    Dataset = discovery_metadata$dataset,
    col = list(
      Group = c(Normal = "#0072B2", Tumor = "#D55E00"),
      Dataset = c(GSE42568 = unname(CB_PALETTE[1]), GSE15852 = unname(CB_PALETTE[2]), GSE29431 = unname(CB_PALETTE[3]))
    )
  )
  
  ht <- Heatmap(hm_data,
                name = "Z-score",
                col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
                top_annotation = col_ann,
                left_annotation = prog_ann,
                show_row_names = FALSE,
                show_column_names = FALSE,
                column_title = "Regulated Cell Death Gene Expression",
                clustering_distance_columns = "euclidean",
                clustering_method_columns = "ward.D2",
                row_km = 3)
  
  open_tiff(file.path(OUTPUT$shared, "programme_expression_heatmap.tiff"), width = 10, height = 8)
  draw(ht)
  dev.off()
}

# --- 8. FOLD CHANGE COMPARISON ---
fc_comparison <- data.frame(
  Gene = all_de_genes,
  logFC = de_results[all_de_genes, "logFC"],
  Program = case_when(
    all_de_genes %in% de_frg_set & all_de_genes %in% de_nrg_set & all_de_genes %in% de_crg_set ~ "All Three",
    all_de_genes %in% de_frg_set & all_de_genes %in% de_nrg_set ~ "FRG+NRG",
    all_de_genes %in% de_frg_set & all_de_genes %in% de_crg_set ~ "FRG+CRG",
    all_de_genes %in% de_nrg_set & all_de_genes %in% de_crg_set ~ "NRG+CRG",
    all_de_genes %in% de_frg_set ~ "FRG Only",
    all_de_genes %in% de_nrg_set ~ "NRG Only",
    all_de_genes %in% de_crg_set ~ "CRG Only",
    TRUE ~ "Other"
  ),
  stringsAsFactors = FALSE
)

p_fc <- ggplot(fc_comparison, aes(x = Program, y = logFC, fill = Program)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
  scale_fill_manual(values = c("All Three" = "purple", "FRG+NRG" = "#E69F00", 
                               "FRG+CRG" = "#CC79A7", "NRG+CRG" = "#56B4E9",
                               "FRG Only" = "#E69F00", "NRG Only" = "#56B4E9", 
                               "CRG Only" = "#009E73")) +
  labs(title = "Fold Change by Cross-Program Category",
       x = "", y = "log2 Fold Change") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))
save_tiff_600dpi(p_fc, file.path(OUTPUT$shared, "fold_change_by_category.tiff"))

# --- 9. EFFECT DIRECTION HEATMAP ACROSS DATASETS ---
per_dataset_deg_files <- list.files(OUTPUT$de, pattern = "DEG_GSE.*\\.csv", full.names = TRUE)
if (length(per_dataset_deg_files) > 0) {
  effect_dir <- data.frame(Gene = all_de_genes)
  for (f in per_dataset_deg_files) {
    ds_name <- gsub("DEG_|\\.csv", "", basename(f))
    ds_de <- read.csv(f, row.names = 1)
    effect_dir[[ds_name]] <- sign(ds_de[all_de_genes, "logFC"])
  }
  
  write.csv(effect_dir, file.path(OUTPUT$shared, "effect_direction_across_datasets.csv"), row.names = FALSE)
}

# --- 10. SHARED GENE ANALYSIS ---
all_program_genes <- unique(c(ferr_analysis$de_genes, necr_analysis$de_genes, crg_analysis$de_genes))

# Classify shared genes
shared_by_two <- intersect(ferr_analysis$de_genes, necr_analysis$de_genes)  # FRG ∩ NRG
shared_by_three <- Reduce(intersect, list(ferr_analysis$de_genes, necr_analysis$de_genes, crg_analysis$de_genes))

ferr_necr_shared <- setdiff(intersect(ferr_analysis$de_genes, necr_analysis$de_genes), shared_by_three)
ferr_crg_shared <- setdiff(intersect(ferr_analysis$de_genes, crg_analysis$de_genes), shared_by_three)
necr_crg_shared <- setdiff(intersect(necr_analysis$de_genes, crg_analysis$de_genes), shared_by_three)

ferr_specific <- setdiff(ferr_analysis$de_genes, union(necr_analysis$de_genes, crg_analysis$de_genes))
necr_specific <- setdiff(necr_analysis$de_genes, union(ferr_analysis$de_genes, crg_analysis$de_genes))
crg_specific <- setdiff(crg_analysis$de_genes, union(ferr_analysis$de_genes, necr_analysis$de_genes))

shared_info <- data.frame(
  Gene = all_program_genes,
  DE_FRG = all_program_genes %in% ferr_analysis$de_genes,
  DE_NRG = all_program_genes %in% necr_analysis$de_genes,
  DE_CRG = all_program_genes %in% crg_analysis$de_genes,
  Category = case_when(
    all_program_genes %in% shared_by_three ~ "All Three",
    all_program_genes %in% ferr_necr_shared ~ "Ferroptosis + Necroptosis",
    all_program_genes %in% ferr_crg_shared ~ "Ferroptosis + Cuproptosis",
    all_program_genes %in% necr_crg_shared ~ "Necroptosis + Cuproptosis",
    all_program_genes %in% ferr_specific ~ "Ferroptosis Only",
    all_program_genes %in% necr_specific ~ "Necroptosis Only",
    all_program_genes %in% crg_specific ~ "Cuproptosis Only"
  ),
  logFC = de_results[all_program_genes, "logFC"],
  adj.P.Val = de_results[all_program_genes, "adj.P.Val"],
  stringsAsFactors = FALSE
)

write.csv(shared_info, file.path(OUTPUT$shared, "shared_gene_classification.csv"), row.names = FALSE)
saveRDS(shared_info, file.path(OUTPUT$shared, "shared_gene_classification.rds"))

log_message("Regulated cell-death gene analysis complete.")
