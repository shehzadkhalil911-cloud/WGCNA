
setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
library(GEOquery)
library(Biobase)

# Check GPL570
cat('=== GPL570 ===\n')
gpl570 <- getGEO('GPL570', destdir='data/raw')
gpl_tbl <- Table(gpl570)
cat('Columns:', paste(colnames(gpl_tbl)[1:30], collapse=' | '), '\n')
cat('Has Gene.Symbol:', 'Gene.Symbol' %in% colnames(gpl_tbl), '\n')
cat('Has Symbol:', 'Symbol' %in% colnames(gpl_tbl), '\n')
cat('Has ID:', 'ID' %in% colnames(gpl_tbl), '\n')
# Check first few rows of relevant columns
for (col in intersect(colnames(gpl_tbl), c('Gene.Symbol','Symbol','Gene symbol','GeneSymbol','GENE_SYMBOL'))) {
    cat('  Column', col, ':', paste(head(gpl_tbl[[col]], 5), collapse=' | '), '\n')
}
# Search for gene-related columns
gene_cols <- grep('gene|Gene|GENE|SYMBOL|Symbol', colnames(gpl_tbl), value=TRUE)
cat('Gene-related columns:', paste(gene_cols, collapse=' | '), '\n')

