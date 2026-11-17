# Run the whole analysis from the repository root:
#   Rscript run_all.R
# Slow steps (genome-wide screens, 100-split evaluations, forest subsampling)
# are cached in output/, so a second run only redraws figures and tables.
steps <- c("R/01_pca.R",
           "R/02_early_death_screen.R",
           "R/03_aqp1_integration.R",
           "R/04_clinical_survival.R",
           "R/05_cox_screen.R",
           "R/06_penalised_cox.R",
           "R/07_random_survival_forest.R",
           "R/08_compare_models.R")
for (f in steps) {
  message("\n==== ", f, " ====")
  t0 <- Sys.time()
  source(f, echo = FALSE, local = new.env())
  message("done in ", format(round(Sys.time() - t0, 1)))
}
