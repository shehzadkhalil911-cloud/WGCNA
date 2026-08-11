# 05_integrate_discovery_datasets.R
# Integrate GSE42568, GSE15852, GSE29431 into a shared-gene discovery matrix.
source("config/analysis_config.R")
source("R/load_functions.R")
ensure_paths()
suppressMessages({ library(dplyr); library(ComplexUpset); library(ggplot2); library(UpSetR) })

## load per-dataset gene matrices produced by 03
gene_mats <- list()
meta_list <- list()
for (g in DISCOVERY_GSE) {
  m <- readRDS(file.path(paths$processed_ind, paste0(g, "_gene_expression.rds")))
  ## log2 safety: if clearly unlogged, transform
  if (max(m, na.rm = TRUE) > 100 && median(m, na.rm = TRUE) < 8) {
    m <- log2(m + 1)
    cat("Applied log2 to", g, "\n")
  }
  gene_mats[[g]] <- m
}

## shared genes across all discovery datasets
shared <- Reduce(intersect, lapply(gene_mats, rownames))
cat("Shared genes across discovery datasets:", length(shared), "\n")

## restrict + merge
mats_shared <- lapply(gene_mats, function(m) m[shared, ])
disc_expr <- do.call(cbind, mats_shared)
disc_meta <- data.frame(
  sample_id = colnames(disc_expr),
  stringsAsFactors = FALSE
)
## assign dataset by column-block order (gene_mats already in DISCOVERY_GSE order)
ncols <- sapply(mats_shared, ncol)
cs <- c(0, cumsum(ncols))
disc_meta$dataset <- ""
for (i in seq_along(DISCOVERY_GSE)) {
  disc_meta$dataset[(cs[i]+1):cs[i+1]] <- DISCOVERY_GSE[i]
}
## attach class from curated metadata
all_meta <- readRDS(file.path(paths$metadata, "sample_group_assignment.rds"))
disc_meta$class <- all_meta$assigned_group[match(disc_meta$sample_id, all_meta$sample_id)]
disc_meta$group <- disc_meta$class
disc_meta$platform <- all_meta$platform[match(disc_meta$sample_id, all_meta$sample_id)]

## per-dataset gene coverage
cov <- data.frame(
  dataset = DISCOVERY_GSE,
  n_genes = sapply(gene_mats, nrow),
  n_shared = length(shared),
  stringsAsFactors = FALSE
)
write.csv(cov, file.path(paths$discovery, "gene_coverage.csv"), row.names = FALSE)

## sample composition
comp <- disc_meta %>% group_by(dataset, group) %>% summarise(n = n(), .groups = "drop")
write.csv(comp, file.path(paths$discovery, "sample_composition.csv"), row.names = FALSE)

## before-integration density (sample genes for speed)
set.seed(GLOBAL_SEED)
sg <- sample(nrow(disc_expr), 2000)
dens_df <- stack(as.data.frame(t(disc_expr[sg, ])))
names(dens_df) <- c("expr", "sample")
dens_df$dataset <- sub("\\..*", "", as.character(dens_df$sample))
p_dens <- ggplot(dens_df, aes(x = expr, colour = dataset)) + geom_density(alpha = 0.5) +
  labs(title = "Pre-integration expression density (shared genes)") + pub_theme()

## shared-gene Upset plot
set_lists <- lapply(gene_mats, function(m) rownames(m))
names(set_lists) <- DISCOVERY_GSE
p_upset <- UpSetR::upset(fromList(set_lists), order.by = "freq",
                         main.bar.color = "#4E79A7", sets.bar.color = "#E15759",
                         text.scale = 1.1)

saveRDS(list(expression = disc_expr, metadata = disc_meta),
        file.path(paths$discovery, "discovery_cohort_raw.rds"))
saveRDS(disc_expr, file.path(paths$discovery, "discovery_expression.rds"))
write.csv(as.data.frame(disc_expr), file.path(paths$discovery, "discovery_expression.csv"))
saveRDS(shared, file.path(paths$discovery, "shared_genes.rds"))

save_publication_plot(p_dens, file.path(paths$res_pca, "preintegration_density"))
## upset plot saved directly (UpSetR returns non-ggplot object)
png(file.path(paths$res_pca, "shared_gene_upset.png"), width = 7, height = 5, units = "in", res = FIG_DPI)
UpSetR::upset(fromList(set_lists), order.by = "freq",
              main.bar.color = "#4E79A7", sets.bar.color = "#E15759", text.scale = 1.1)
dev.off()
tiff(file.path(paths$res_pca, "shared_gene_upset.tiff"), width = 7, height = 5, units = "in", res = 600, compression = "lzw")
UpSetR::upset(fromList(set_lists), order.by = "freq",
              main.bar.color = "#4E79A7", sets.bar.color = "#E15759", text.scale = 1.1)
dev.off()
pdf(file.path(paths$res_pca, "shared_gene_upset.pdf"), width = 7, height = 5)
UpSetR::upset(fromList(set_lists), order.by = "freq",
              main.bar.color = "#4E79A7", sets.bar.color = "#E15759", text.scale = 1.1)
dev.off()

cat("Discovery integration complete. Matrix:", nrow(disc_expr), "x", ncol(disc_expr), "\n")
cat("Class counts:\n"); print(table(disc_meta$group))
quit(save = "no", status = 0, runLast = FALSE)
