
setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
library(Biobase)

# GSE42568
cat('=== GSE42568 phenotype ===\n')
f <- readRDS('data/raw/GSE42568_eset.rds')
p <- pData(f)
cat('Columns:', paste(colnames(p), collapse=' | '), '\n')
cat('Title sample:', paste(head(as.character(p), 5), collapse='\n  '), '\n')
cat('Source:', paste(head(as.character(p), 5), collapse='\n  '), '\n')
# Check characteristics
char_cols <- grep('characteristics', colnames(p), value=TRUE)
for (cc in char_cols) {
    cat('Characteristic:', cc, '->', paste(head(as.character(p[[cc]]), 5), collapse=' | '), '\n')
}
# Check all unique source_name values
cat('Source values:', paste(unique(as.character(p)), collapse=' | '), '\n')

cat('\n=== GSE15852 phenotype ===\n')
f2 <- readRDS('data/raw/GSE15852_eset.rds')
p2 <- pData(f2)
cat('Columns:', paste(colnames(p2), collapse=' | '), '\n')
cat('Title sample:', paste(head(as.character(p2), 5), collapse='\n  '), '\n')
cat('Source:', paste(head(as.character(p2), 5), collapse='\n  '), '\n')
char_cols2 <- grep('characteristics', colnames(p2), value=TRUE)
for (cc in char_cols2) {
    cat('Characteristic:', cc, '->', paste(head(as.character(p2[[cc]]), 5), collapse=' | '), '\n')
}

