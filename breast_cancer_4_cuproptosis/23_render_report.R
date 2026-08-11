# 23_render_report.R
source("config/analysis_config.R")
source("R/load_functions.R")
rmd <- file.path(paths$report, "manuscript.Rmd")
out_html <- file.path(paths$report, "manuscript.html")
out_pdf <- file.path(paths$report, "manuscript.pdf")

## Prefer quarto CLI if available, else rmarkdown
qbin <- Sys.which("quarto")
if (nzchar(qbin) && file.exists(qbin)) {
  cat("Rendering with quarto CLI...\n")
  system(paste(shQuote(qbin), "render", shQuote(rmd), "--to html"))
} else {
  suppressMessages(library(rmarkdown))
  out_dir_abs <- normalizePath(paths$report, mustWork = FALSE)
  cat("Rendering with rmarkdown (html)...\n")
  rmarkdown::render(rmd, output_format = "html_document", output_file = "manuscript.html",
                    output_dir = out_dir_abs, params = list(root = "."), quiet = TRUE,
                    envir = globalenv(), knit_root_dir = getwd())
  cat("HTML written:", file.path(out_dir_abs, "manuscript.html"), "\n")
  ## attempt PDF (needs a LaTeX engine); non-fatal if unavailable
  pdf_ok <- tryCatch({
    rmarkdown::render(rmd, output_format = "pdf_document", output_file = "manuscript.pdf",
                      output_dir = out_dir_abs, params = list(root = "."), quiet = TRUE,
                      envir = globalenv(), knit_root_dir = getwd())
    TRUE
  }, error = function(e) FALSE)
  cat("PDF written:", if (pdf_ok) out_pdf else "SKIPPED (no LaTeX engine)", "\n")
}
cat("Report rendering complete.\n")
quit(save = "no", status = 0, runLast = FALSE)
