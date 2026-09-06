# =============================================================================
# Post-hoc convergence diagnostics for the full MH-corrected BMEOP run
# (BMEOP only -- run in a separate session from any BJCM diagnostic script)
# =============================================================================
# Run this AFTER run_full_model_only_bmeop.R and compare_saved_results_bmeop.R
# have completed. Nothing here re-fits any model -- it only loads the saved
# result object and computes diagnostics, so it runs in seconds.
#
# Structure, in order:
#   1. Threshold diagnostics (SD, R-hat, ESS) -- BMEOP has a single
#      threshold set (the loneliness outcome), unlike BJCM's separate
#      mediator/outcome sets.
#   2. Reported-estimand-level diagnostic (ATE): verifies whether any
#      threshold-level mixing imperfection found in step 1 propagates to
#      the actual estimand reported in the paper.
# =============================================================================

library(coda)

mh_result <- readRDS("../results/bmeop_full_results.rds")

# -----------------------------------------------------------------------------
# 1. Threshold diagnostics
# -----------------------------------------------------------------------------
cat("=== Threshold SD, per chain ===\n")
for (chain_i in 1:4) {
  th <- mh_result$chains[[chain_i]]$thresholds
  if (is.null(th)) {
    cat(sprintf("Chain %d: could not find $thresholds -- inspect structure below\n", chain_i))
    next
  }
  sds <- apply(th[, -1, drop = FALSE], 2, sd, na.rm = TRUE)
  cat(sprintf("Chain %d: ", chain_i)); print(round(sds, 4))
}

th_list <- lapply(1:4, function(i) mh_result$chains[[i]]$thresholds[, -1])
th_mcmc <- mcmc.list(lapply(th_list, mcmc))

cat("\n=== Threshold R-hat ===\n")
print(gelman.diag(th_mcmc, multivariate = FALSE))

cat("\n=== Threshold ESS ===\n")
print(effectiveSize(th_mcmc))

# -----------------------------------------------------------------------------
# 2. Reported-estimand-level diagnostic: ATE
# -----------------------------------------------------------------------------
# If any threshold above showed elevated R-hat, this checks whether that
# imperfection propagates to the ATE actually reported in the paper, or
# washes out under averaging -- as it does for the raw regression
# coefficients (Table 2), whose individual R-hat values run higher than
# the ATE's own R-hat = 1.01.

ate_list <- lapply(1:4, function(chain_i) {
  mh_result$chains[[chain_i]]$ate
  # field name based on the general bayesian_bmeop_elsa() output format
  # (top-level object has $ate_samples pooled across chains); if this
  # per-chain field name is wrong, check str(mh_result$chains[[1]],
  # max.level = 1) first
})

if (any(sapply(ate_list, is.null))) {
  cat("\n=== Could not find per-chain $ate -- inspect structure ===\n")
  str(mh_result$chains[[1]], max.level = 1)
} else {
  ate_mcmc <- mcmc.list(lapply(ate_list, mcmc))
  
  cat("\n=== ATE (reported estimand) R-hat ===\n")
  print(gelman.diag(ate_mcmc))
  
  cat("\n=== ATE (reported estimand) ESS ===\n")
  print(effectiveSize(ate_mcmc))
}