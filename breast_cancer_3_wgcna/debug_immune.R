load("F:/independent projects/cancer/breast cancer/breast cancer 3 (WCGNA)/RData/08_immune_infiltration.RData")

for (ds in c("GSE42568", "GSE15852", "GSE29431")) {
  cat("\n=== ", ds, " ===\n")
  d <- immune_results[[ds]][["diff"]]
  d <- d[order(d$FDR), ]
  cat("Top 8 by FDR:\n")
  print(d[1:8, c("CellType", "logFC", "FDR")])
  cat("Total sig (FDR<0.05):", sum(d$FDR < 0.05), "\n")
}
