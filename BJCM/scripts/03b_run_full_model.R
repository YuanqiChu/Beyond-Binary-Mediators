# =============================================================================
# Full-scale BJCM model fit
# =============================================================================
# Runs the complete four-chain fit and saves the raw result object to
# ../results/bjcm_full_results.rds. Run this after 03a_run_reduced_check.R
# has shown healthy threshold-update acceptance rates (target ~20-50%).
#
# This will take several hours; consider running it in the background via
# run_in_background() from 00_background_run_utils.R rather than
# interactively.
# =============================================================================

library(dplyr)
library(MASS)
library(mvtnorm)
library(coda)
library(Matrix)
library(parallel)

source("00_background_run_utils.R")
source("02_bjcm_model.R")

elsa_long <- read.csv("../data/elsa_longitudinal_analysis.csv")

reset_mh_threshold_log()

cat("=== Running full-scale BJCM fit ===\n")
cat("Consider running in the background instead of interactively:\n")
cat('  job <- run_in_background(run_full_fit(elsa_long, n_chains = 4))\n')
cat('  monitor_progress("bjcm")\n')
cat("  full_result <- collect_result(job)\n")
cat("  (then skip the two lines below and go straight to saveRDS)\n\n")

t0 <- Sys.time()
full_result <- run_full_fit(elsa_long, n_chains = 4)
t1 <- Sys.time()
cat("\nFull run wall-clock time:", round(as.numeric(difftime(t1, t0, units = "hours")), 2), "hours\n")

saveRDS(full_result, "../results/bjcm_full_results.rds")

acceptance_summary <- get_mh_acceptance_summary()
cat("\n=== Threshold acceptance rates (full run) ===\n")
print(acceptance_summary)
saveRDS(acceptance_summary, "../results/mh_acceptance_summary.rds")

cat("\n=== Done. Files saved: ===\n")
cat("  ../results/bjcm_full_results.rds\n")
cat("  ../results/mh_acceptance_summary.rds\n\n")
cat("Next: run 06a_extract_final_results.R to extract the manuscript's\n")
cat("reported numbers from the saved fit.\n")
