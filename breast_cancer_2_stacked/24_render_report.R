setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
# =============================================================================
# 24_render_report.R
# Render the Quarto manuscript report
# =============================================================================

source("config/analysis_config.R")
cat("Rendering Quarto reports...\n")

# Check if quarto is available
has_quarto <- nchar(Sys.which("quarto")) > 0

if (has_quarto) {
    # Render main manuscript
    qmd_file <- file.path(REPORT_DIR, "manuscript_report.qmd")
    if (file.exists(qmd_file)) {
        quarto::quarto_render(
            input = qmd_file,
            output_format = c("html", "docx", "pdf"),
            output_dir = REPORT_DIR
        )
        cat("Manuscript report rendered.\n")
    } else {
        cat("manuscript_report.qmd not found. Creating skeleton...\n")
    }
    
    # Render supplementary
    supp_file <- file.path(REPORT_DIR, "supplementary_report.qmd")
    if (file.exists(supp_file)) {
        quarto::quarto_render(
            input = supp_file,
            output_format = c("html", "docx"),
            output_dir = REPORT_DIR
        )
        cat("Supplementary report rendered.\n")
    }
} else {
    cat("Quarto not found. Reports must be rendered manually.\n")
    cat("Install Quarto from https://quarto.org/docs/get-started/\n")
}

# Save session info
sink("sessionInfo.txt")
sessionInfo()
sink()
cat("Session information saved.\n")

cat("\nPipeline complete.\n")

