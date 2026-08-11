load("F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)/RData/01_processed_data_full.RData")

cat("=== processed_42568 structure ===\n")
cat("str(processed_42568):\n")
str(processed_42568, max.level = 2, give.attr = FALSE)

cat("\n=== processed_42568[[1]] ===\n")
cat("class:", class(processed_42568[[1]]), "\n")
cat("dim:", paste(dim(processed_42568[[1]]), collapse=" x "), "\n")
cat("rownames[1:5]:", paste(rownames(processed_42568[[1]])[1:5], collapse=", "), "\n")

cat("\n=== processed_42568[[2]] ===\n")
cat("class:", class(processed_42568[[2]]), "\n")
cat("names:", paste(names(processed_42568[[2]]), collapse=", "), "\n")

cat("\n=== common_genes_discovery[1:10] ===\n")
cat(paste(common_genes_discovery[1:10], collapse=", "), "\n")

cat("\n=== hub check ===\n")
load("F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)/RData/05_hub_genes.RData")
cat("final_hubs[1:5]:", paste(final_hubs[1:5], collapse=", "), "\n")
cat("final_hubs %in% rownames(processed_42568[[1]]):", sum(final_hubs %in% rownames(processed_42568[[1]])), "/50\n")
cat("common_genes_discovery %in% rownames(processed_42568[[1]]):", sum(common_genes_discovery %in% rownames(processed_42568[[1]])), "/6932\n")
