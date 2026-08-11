# Generate Word document for selected KEGG pathways
# Fatty acid degradation, PPAR signaling, Cell cycle, Cytoskeleton in muscle cells

library(ggplot2)
library(officer)

base_dir <- normalizePath("breast_cancer_integrated_analysis")
tables_dir <- file.path(base_dir, "results", "tables")
fig_dir <- file.path(base_dir, "results", "figures", "Enrichment")
out_dir <- file.path(base_dir, "report")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

kegg <- read.csv(file.path(tables_dir, "KEGG_enrichment.csv"))

pathway_info <- list(
  list(name = "Fatty acid degradation", id = "hsa00071"),
  list(name = "PPAR signaling pathway", id = "hsa03320"),
  list(name = "Cell cycle", id = "hsa04110"),
  list(name = "Cytoskeleton in muscle cells", id = "hsa04820")
)

doc <- read_docx()
doc <- officer::body_add_par(doc, "Selected KEGG Pathway Maps", style = "heading 1")
doc <- officer::body_add_par(doc, paste("Generated:", Sys.Date()), style = "Normal")
doc <- officer::body_add_par(doc, "", style = "Normal")

for (pw in pathway_info) {
  doc <- officer::body_add_par(doc, pw$name, style = "heading 2")
  png_path <- file.path(fig_dir, paste0(pw$id, ".png"))
  if (file.exists(png_path)) {
    doc <- officer::body_add_img(doc, src = png_path, width = 6, height = 5)
  } else {
    doc <- officer::body_add_par(doc, "(Pathway map not found)", style = "Normal")
  }
  doc <- officer::body_add_break(doc)
}

doc <- officer::body_add_par(doc, "Enrichment Summary", style = "heading 1")
selected <- kegg[kegg[["Description"]] %in% sapply(pathway_info, `[[`, "name"), ]
ft <- flextable::flextable(selected[, c("Description", "Count", "pvalue", "p.adjust")])
ft <- flextable::theme_box(ft)
ft <- flextable::autofit(ft)
ft <- flextable::fontsize(ft, size = 10, part = "all")
ft <- flextable::bold(ft, part = "header")
doc <- flextable::body_add_flextable(doc, ft)

out_path <- file.path(out_dir, "Selected_KEGG_Pathways.docx")
print(doc, target = out_path)
cat("Saved:", out_path, "\n")
