###############################################################################
# BATCH CONVERSION: All existing PNGs -> TIFF 600 dpi
###############################################################################
library(magick)

FIG_DIR <- "F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)/figures"
png_files <- list.files(FIG_DIR, pattern = "\\.png$", recursive = TRUE, full.names = TRUE)
cat(sprintf("Found %d PNG files to convert\n", length(png_files)))

for (f in png_files) {
  tiff_file <- sub("\\.png$", ".tiff", f)
  
  # Skip if TIFF already exists and is newer
  if (file.exists(tiff_file) && file.info(tiff_file)$mtime > file.info(f)$mtime) {
    next
  }
  
  img <- image_read(f)
  img <- image_convert(img, format = "tiff")
  image_write(img, tiff_file, density = 600)
  cat(sprintf("  Converted: %s -> %s\n", basename(f), basename(tiff_file)))
}

cat(sprintf("\nDone. %d files converted to TIFF at 600 dpi.\n", length(png_files)))
