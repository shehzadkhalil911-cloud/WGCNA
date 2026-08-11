setwd('F:/independent projects/cancer/breast cancer/breast cancer 2/BreastCancer_Stacked_Ensemble')
options(warn = 1)
tryCatch({
  source('R/package_functions.R')
  cat('package_functions sourced OK\n')
  cat('exists load_required_packages:', exists('load_required_packages'), '\n')
  cat('exists install_cran_packages:', exists('install_cran_packages'), '\n')
}, error = function(e) {
  cat('ERROR:', conditionMessage(e), '\n')
})
