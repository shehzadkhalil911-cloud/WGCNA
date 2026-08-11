
setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
library(GEOquery)
for (gpl in c('GPL570','GPL96','GPL13607')) {
    cat('===', gpl, '===\n')
    g <- getGEO(gpl, destdir='data/raw')
    t <- Table(g)
    cat('Columns:', ncol(t), '\n')
    gene_cols <- grep('gene|Gene|GENE|SYMBOL|Symbol|symbol', colnames(t), value=TRUE)
    cat('Gene/Symbol columns:', paste(gene_cols, collapse=' | '), '\n')
    cat('ID column present:', 'ID' %in% colnames(t), '\n')
    # First value of first gene column
    if (length(gene_cols) > 0) {
        cat('  Sample values:', paste(head(t[[gene_cols[1]]], 3), collapse=' | '), '\n')
    }
    cat('\n')
}

