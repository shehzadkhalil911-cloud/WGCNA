source("R/utils.R")
load_libraries()
set.seed(123)

BASE_DIR <- "F:/independent projects/cancer/breast cancer"
OUT_DIR <- file.path(BASE_DIR, "Breast cancer 6 Ferroptosis")

log_msg("Generating enrichment figures from saved CSV results...")

enrich_dir <- file.path(OUT_DIR, "results/enrichment")
fig_dir <- file.path(OUT_DIR, "results/figures")

enrich_dotplot <- function(csv_path, title, out_name, n_top = 15) {
  if (!file.exists(csv_path)) { log_msg(paste("  skip - no file:", basename(csv_path))); return() }
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  if (nrow(df) == 0) { log_msg(paste("  skip - empty:", basename(csv_path))); return() }
  
  df <- df[order(df$p.adjust), ]
  df <- head(df, n_top)
  df$Description <- factor(df$Description, levels = rev(df$Description))
  
  p <- ggplot(df, aes(Count, Description)) +
    geom_point(aes(size = Count, color = p.adjust)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(title = title, x = "Gene Count", y = "", color = "Adj. P", size = "Count") +
    theme_minimal() + theme(axis.text.y = element_text(size = 7))
  save_tiff_600dpi(p, file.path(fig_dir, out_name), width = 9, height = max(4, n_top * 0.35))
  log_msg(paste("  saved:", out_name))
}

enrich_barplot <- function(csv_path, title, out_name, n_top = 15) {
  if (!file.exists(csv_path)) { log_msg(paste("  skip - no file:", basename(csv_path))); return() }
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  if (nrow(df) == 0) { log_msg(paste("  skip - empty:", basename(csv_path))); return() }
  
  df <- df[order(df$p.adjust), ]
  df <- head(df, n_top)
  df$Description <- factor(df$Description, levels = rev(df$Description))
  df$logP <- -log10(df$p.adjust)
  
  p <- ggplot(df, aes(logP, Description)) +
    geom_col(aes(fill = p.adjust)) +
    scale_fill_gradient(low = "red", high = "blue") +
    labs(title = title, x = expression(-log[10](Adj.P)), y = "", fill = "Adj. P") +
    theme_minimal()
  save_tiff_600dpi(p, file.path(fig_dir, out_name), width = 9, height = max(4, n_top * 0.35))
  log_msg(paste("  saved:", out_name))
}

# All DEGs - GO BP
enrich_dotplot(file.path(enrich_dir, "all_DEGs_GO_BP.csv"), "GO Biological Process - All DEGs", "Figure8A_GO_BP_dot.tiff")
enrich_barplot(file.path(enrich_dir, "all_DEGs_GO_BP.csv"), "GO Biological Process - All DEGs", "Figure8B_GO_BP_bar.tiff")
# KEGG
enrich_dotplot(file.path(enrich_dir, "all_DEGs_KEGG.csv"), "KEGG Pathway - All DEGs", "Figure8C_KEGG_dot.tiff")
enrich_barplot(file.path(enrich_dir, "all_DEGs_KEGG.csv"), "KEGG Pathway - All DEGs", "Figure8D_KEGG_bar.tiff")
# Reactome
enrich_dotplot(file.path(enrich_dir, "all_DEGs_REACTOME.csv"), "Reactome - All DEGs", "Figure8E_Reactome_dot.tiff")
enrich_barplot(file.path(enrich_dir, "all_DEGs_REACTOME.csv"), "Reactome - All DEGs", "Figure8F_Reactome_bar.tiff")

# Up DEGs - GO BP
enrich_dotplot(file.path(enrich_dir, "up_DEGs_GO_BP.csv"), "GO BP - Upregulated DEGs", "Figure8G_up_GO_dot.tiff")
enrich_dotplot(file.path(enrich_dir, "up_DEGs_KEGG.csv"), "KEGG - Upregulated DEGs", "Figure8H_up_KEGG_dot.tiff")
# Down DEGs
enrich_dotplot(file.path(enrich_dir, "down_DEGs_GO_BP.csv"), "GO BP - Downregulated DEGs", "Figure8I_down_GO_dot.tiff")
enrich_dotplot(file.path(enrich_dir, "down_DEGs_KEGG.csv"), "KEGG - Downregulated DEGs", "Figure8J_down_KEGG_dot.tiff")

# DE-FRGs
enrich_dotplot(file.path(enrich_dir, "DE_FRGs_GO_BP.csv"), "GO BP - DE-FRGs", "Figure8K_DE_FRG_GO_dot.tiff")
enrich_dotplot(file.path(enrich_dir, "DE_FRGs_KEGG.csv"), "KEGG - DE-FRGs", "Figure8L_DE_FRG_KEGG_dot.tiff")
enrich_barplot(file.path(enrich_dir, "DE_FRGs_GO_BP.csv"), "GO BP - DE-FRGs", "Figure8M_DE_FRG_GO_bar.tiff")

# Hallmark
enrich_dotplot(file.path(enrich_dir, "all_DEGs_HALLMARK.csv"), "Hallmark - All DEGs", "Figure8N_Hallmark_dot.tiff")
enrich_dotplot(file.path(enrich_dir, "DE_FRGs_HALLMARK.csv"), "Hallmark - DE-FRGs", "Figure8O_DE_FRG_Hallmark_dot.tiff")

log_msg("All enrichment figures generated")
