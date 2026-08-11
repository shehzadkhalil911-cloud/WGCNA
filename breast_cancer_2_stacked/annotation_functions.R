# =============================================================================
# annotation_functions.R
# Functions for probe-to-gene mapping, annotation, and ID conversion
# =============================================================================

#' Map probes to gene symbols using platform annotation
map_probes_to_genes <- function(probe_ids, gpl_table = NULL, platform = NULL) {
    # Find the gene symbol column (various naming conventions)
    symbol_col <- NULL
    possible_names <- c("Gene.Symbol", "Gene Symbol", "Symbol", "GENE_SYMBOL",
                        "GeneSymbol", "gene_symbol", "symbol", "GENE NAME",
                        "GeneName", "GENE_NAME", "gene_name")
    for (pn in possible_names) {
        if (!is.null(gpl_table) && pn %in% colnames(gpl_table)) {
            symbol_col <- pn
            break
        }
    }
    
    if (!is.null(symbol_col)) {
        idx <- match(probe_ids, gpl_table$ID)
        symbols <- gpl_table[[symbol_col]][idx]
    } else {
        # Fallback: check any column containing "gene" and "symbol"
        if (!is.null(gpl_table)) {
            gene_cols <- grep("gene.*symbol|symbol.*gene", 
                              colnames(gpl_table), ignore.case = TRUE, value = TRUE)
            if (length(gene_cols) > 0) {
                idx <- match(probe_ids, gpl_table$ID)
                symbols <- gpl_table[[gene_cols[1]]][idx]
            } else {
                symbols <- rep(NA, length(probe_ids))
            }
        } else {
            symbols <- rep(NA, length(probe_ids))
        }
    }
    # Clean symbols: split on "///" and take first
    symbols <- gsub("///.*", "", symbols)
    symbols <- trimws(symbols)
    symbols[symbols == ""] <- NA
    # Remove any remaining garbage (e.g., "---")
    symbols[grepl("^---+$", symbols)] <- NA
    names(symbols) <- probe_ids
    return(symbols)
}

#' Collapse duplicate probes using configurable method
collapse_duplicates <- function(expr_matrix, gene_map, method = "highest_IQR") {
    genes <- gene_map[rownames(expr_matrix)]
    valid <- !is.na(genes) & genes != ""
    expr_matrix <- expr_matrix[valid, , drop = FALSE]
    genes <- genes[valid]

    unique_genes <- unique(genes)
    collapsed <- matrix(NA, nrow = length(unique_genes), ncol = ncol(expr_matrix))
    colnames(collapsed) <- colnames(expr_matrix)
    rownames(collapsed) <- unique_genes

    for (g in unique_genes) {
        idx <- which(genes == g)
        if (length(idx) == 1) {
            collapsed[g, ] <- expr_matrix[idx, ]
        } else {
            if (method == "highest_IQR") {
                iqrs <- apply(expr_matrix[idx, , drop = FALSE], 1, IQR, na.rm = TRUE)
                best <- idx[which.max(iqrs)]
                collapsed[g, ] <- expr_matrix[best, ]
            } else if (method == "mean") {
                collapsed[g, ] <- colMeans(expr_matrix[idx, , drop = FALSE], na.rm = TRUE)
            } else {
                collapsed[g, ] <- colMeans(expr_matrix[idx, , drop = FALSE], na.rm = TRUE)
            }
        }
    }
    return(collapsed)
}

#' Map gene symbols to Entrez IDs
symbol_to_entrez <- function(symbols) {
    entrez <- mapIds(
        org.Hs.eg.db,
        keys = symbols,
        column = "ENTREZID",
        keytype = "SYMBOL",
        multiVals = "first"
    )
    return(entrez)
}

#' Map Entrez IDs to gene symbols
entrez_to_symbol <- function(entrez_ids) {
    symbols <- mapIds(
        org.Hs.eg.db,
        keys = as.character(entrez_ids),
        column = "SYMBOL",
        keytype = "ENTREZID",
        multiVals = "first"
    )
    return(symbols)
}

#' Get platform annotation from GEO
get_platform_annotation <- function(gpl_id) {
    gpl <- getGEO(gpl_id, destdir = RAW_DIR, AnnotGPL = FALSE)
    Table(gpl)
}

#' Create annotation summary
annotate_summary <- function(gene_map) {
    mapped <- sum(!is.na(gene_map) & gene_map != "")
    total <- length(gene_map)
    list(
        total_probes = total,
        mapped_genes = mapped,
        unmapped_probes = total - mapped,
        mapping_rate = mapped / total
    )
}
