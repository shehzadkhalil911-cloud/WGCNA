suppressMessages({library(magick); library(dplyr)})
fig_dir <- "F:/independent projects/cancer/breast cancer/Breast cancer 8/figures"
png_dir <- "F:/independent projects/cancer/breast cancer/Breast cancer 8/report/figures_png"
dir.create(png_dir, showWarnings = FALSE, recursive = TRUE)

tiffs <- list.files(fig_dir, pattern = "\\.tiff$", full.names = TRUE)
tiffs <- tiffs[!grepl("^_", basename(tiffs))]
cat("Converting", length(tiffs), "figures to PNG previews...\n")
for (f in tiffs) {
  nm <- sub("\\.tiff$", ".png", basename(f))
  out <- file.path(png_dir, nm)
  tryCatch({
    img <- magick::image_read(f)
    img <- magick::image_resize(img, geometry_size_pixels(width = 1400))
    magick::image_write(img, path = out, format = "png")
  }, error = function(e) cat("  FAILED:", basename(f), "->", conditionMessage(e), "\n"))
}
cat("Done. PNG files:", length(list.files(png_dir, pattern = "\\.png$")), "\n")
