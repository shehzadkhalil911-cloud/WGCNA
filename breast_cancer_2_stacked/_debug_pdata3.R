
setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
library(Biobase)

# GSE15852
cat('=== GSE15852 ===\n')
f <- readRDS('data/raw/GSE15852_eset.rds')
p <- pData(f)
cat('Columns:', paste(colnames(p), collapse=' | '), '\n')
cat('Title sample:', paste(head(as.character(p), 5), collapse=' | '), '\n')
cat('source_name_ch1:', paste(head(as.character(p), 5), collapse=' | '), '\n')
# Check all columns for tumor/normal
for (col in colnames(p)) {
    vals <- unique(as.character(p[[col]]))
    if (length(vals) < 10 && length(vals) > 0) {
        cat(col, ':', paste(vals, collapse=' | '), '\n')
    }
}
# Also check characteristics
char_cols <- grep('characteristics', colnames(p), value=TRUE)
for (cc in char_cols) {
    cat('---', cc, '---\n')
    cat(paste(head(as.character(p[[cc]]), 5), collapse='\n'), '\n')
}

