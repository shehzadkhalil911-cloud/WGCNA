# ============================================================
# 15a - Cross-cohort ComBat batch correction
# Combines discovery (GSE70947, Agilent) and validation
# (GSE15852, Affymetrix) on the common-gene space and applies
# ComBat (batch = cohort) with the class covariate protected,
# producing batch-corrected matrices for both cohorts.
# NOTE: ComBat is fit on both cohorts jointly; the external
# evaluation on corrected data is a harmonized-transfer test,
# not a fully independent prospective test.
# ============================================================
PROJECT_DIR <- "F:/independent projects/cancer/breast cancer/Breast cancer 8"
source(file.path(PROJECT_DIR, "config/config.R"))
suppressMessages({library(dplyr); library(sva)})

logmsg <- function(...) cat(sprintf(...), "\n")

disc <- readRDS(file.path(DIR_DISC_PROC, "discovery_processed.rds"))
val  <- readRDS(file.path(DIR_VAL_PROC, "validation_processed.rds"))
common <- readRDS(file.path(DIR_MODELS, "common_genes.rds"))

E_d <- disc$expr[common, , drop = FALSE]
E_v <- val$expr[common, , drop = FALSE]
M_d <- disc$meta
M_v <- val$meta
stopifnot(identical(rownames(E_d), common), identical(rownames(E_v), common))
stopifnot(all(!is.na(E_d)), all(!is.na(E_v)))
logmsg("Combined matrix: %d genes x %d samples", length(common), ncol(E_d) + ncol(E_v))

combat_dat <- cbind(E_d, E_v)
batch <- factor(c(rep("discovery", ncol(E_d)), rep("validation", ncol(E_v))))
class_vec <- c(M_d$class, M_v$class)
mod <- model.matrix(~ class_vec)

set.seed(SEED)
combat_out <- ComBat(dat = combat_dat, batch = batch, mod = mod,
                     par.prior = TRUE, prior.plots = FALSE)
logmsg("ComBat applied. Check: NAs in corrected data = %d", sum(is.na(combat_out)))

E_d_adj <- combat_out[, seq_len(ncol(E_d)), drop = FALSE]
E_v_adj <- combat_out[, seq_len(ncol(E_v)) + ncol(E_d), drop = FALSE]
colnames(E_d_adj) <- colnames(E_d)
colnames(E_v_adj) <- colnames(E_v)

saveRDS(list(expr = E_d_adj, meta = M_d,
             correction = "ComBat (batch=cohort, mod~class)",
             source = DATA_SOURCE),
        file.path(DIR_COMBAT_PROC, "discovery_combat.rds"))
saveRDS(list(expr = E_v_adj, meta = M_v,
             correction = "ComBat (batch=cohort, mod~class)",
             source = DATA_SOURCE),
        file.path(DIR_COMBAT_PROC, "validation_combat.rds"))
saveRDS(list(corrected = combat_out, batch = batch, mod = mod,
             common_genes = common,
             fit_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        file.path(DIR_COMBINED, "combat_fit.rds"))

# ---------------- sanity report ----------------
top <- c("TUBB3", "EZH2", "TUBA1C", "LMNB1", "MMP11")
sanity <- do.call(rbind, lapply(top, function(g) {
  d <- t.test(E_d_adj[g, M_d$class == "Tumor"], E_d_adj[g, M_d$class == "Normal"])$p.value
  v <- t.test(E_v_adj[g, M_v$class == "Tumor"], E_v_adj[g, M_v$class == "Normal"])$p.value
  data.frame(gene = g,
             discovery_tumor_vs_normal_diff = mean(E_d_adj[g, M_d$class == "Tumor"]) - mean(E_d_adj[g, M_d$class == "Normal"]),
             discovery_p = d,
             validation_tumor_vs_normal_diff = mean(E_v_adj[g, M_v$class == "Tumor"]) - mean(E_v_adj[g, M_v$class == "Normal"]),
             validation_p = v)
}))
write.csv(sanity, file.path(DIR_TABLES, "Table_05c_ComBat_Sanity_Check.csv"), row.names = FALSE)
print(sanity)

logmsg("ComBat corrected matrices saved (discovery %dx%d, validation %dx%d).",
       nrow(E_d_adj), ncol(E_d_adj), nrow(E_v_adj), ncol(E_v_adj))
