
setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
library(Biobase)

# GSE42568 - check tissue column
cat('=== GSE42568 tissue:ch1 ===\n')
f <- readRDS('data/raw/GSE42568_eset.rds')
p <- pData(f)
cat('tissue:ch1 values:', paste(unique(as.character(p[['tissue:ch1']])), collapse=' | '), '\n')
cat('characteristics_ch1:', paste(unique(as.character(p)), collapse=' | '), '\n')
cat('source_name_ch1:', paste(unique(as.character(p)), collapse=' | '), '\n')
cat('title first 3:', paste(head(as.character(p), 3), collapse=' || '), '\n')

# GSE15852
cat('\n=== GSE15852 ===\n')
f2 <- readRDS('data/raw/GSE15852_eset.rds')
p2 <- pData(f2)
cat('Columns:', paste(colnames(p2), collapse=' | '), '\n')
# Check all columns for tumor/normal
for (col in colnames(p2)) {
    vals <- unique(as.character(p2[[col]]))
    if (any(grepl('tumor|normal|cancer', vals, ignore.case=TRUE))) {
        cat(col, ':', paste(head(vals, 10), collapse=' | '), '\n')
    }
}

