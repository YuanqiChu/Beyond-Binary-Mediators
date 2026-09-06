# =============================================================================
# Post-hoc convergence diagnostics for the full MH-corrected BJCM run
# =============================================================================
# Run this AFTER run_full_model_only.R and compare_saved_results.R have
# completed. Nothing here re-fits any model -- it only loads the saved
# result object and computes diagnostics, so it runs in seconds.
#
# Structure, in order:
#   1. Outcome threshold diagnostics (SD, R-hat, ESS)
#   2. Mediator threshold diagnostics (SD, R-hat, ESS), for each mediator
#      that has movable thresholds (transport mobility is binary and is
#      skipped automatically)
#   3. Reported-estimand-level diagnostic (depression NIE): the raw
#      threshold parameters for depression showed elevated R-hat in step 2;
#      this final check verifies whether that imperfection propagates to
#      the actual estimand reported in the paper (it does not).
# =============================================================================

library(coda)

mh_result <- readRDS("../results/bjcm_full_results.rds")

# -----------------------------------------------------------------------------
# 1. Outcome threshold diagnostics
# -----------------------------------------------------------------------------
cat("=== Outcome threshold SD, per chain ===\n")
for (chain_i in 1:4) {
  th_out <- mh_result$results$results$chains[[chain_i]]$thresholds_out
  sds <- apply(th_out[, -1, drop = FALSE], 2, sd, na.rm = TRUE)
  cat(sprintf("Chain %d: ", chain_i)); print(round(sds, 4))
}

th_out_list <- lapply(1:4, function(i) {
  mh_result$results$results$chains[[i]]$thresholds_out[, -1]
})
th_out_mcmc <- mcmc.list(lapply(th_out_list, mcmc))

cat("\n=== Outcome threshold R-hat ===\n")
print(gelman.diag(th_out_mcmc, multivariate = FALSE))

cat("\n=== Outcome threshold ESS ===\n")
print(effectiveSize(th_out_mcmc))

# -----------------------------------------------------------------------------
# 2. Mediator threshold diagnostics
# -----------------------------------------------------------------------------
for (k in 1:3) {
  th_med_list <- lapply(1:4, function(chain_i) {
    mh_result$results$results$chains[[chain_i]]$thresholds_med[[k]]
  })
  
  # Skip mediator 2 (transport mobility) -- binary, no threshold
  if (is.null(th_med_list[[1]]) || !is.matrix(th_med_list[[1]]) || ncol(th_med_list[[1]]) <= 1) {
    cat(sprintf("\n=== Mediator %d: no threshold (binary) ===\n", k))
    next
  }
  
  cat(sprintf("\n=== Mediator %d threshold diagnostics ===\n", k))
  
  for (chain_i in 1:4) {
    sds <- apply(th_med_list[[chain_i]][, -1, drop = FALSE], 2, sd, na.rm = TRUE)
    cat(sprintf("Chain %d SD: ", chain_i)); print(round(sds, 4))
  }
  
  th_med_mcmc <- mcmc.list(lapply(th_med_list, function(m) mcmc(m[, -1, drop = FALSE])))
  
  cat("\nR-hat:\n")
  print(gelman.diag(th_med_mcmc, multivariate = FALSE))
  
  cat("\nESS:\n")
  print(effectiveSize(th_med_mcmc))
}

# -----------------------------------------------------------------------------
# 3. Reported-estimand-level diagnostic: depression NIE
# -----------------------------------------------------------------------------
# Mediator 1 (depression) showed elevated R-hat (1.04-1.11) on two of its
# seven raw threshold parameters above. This final check verifies whether
# that imperfection propagates to the estimand actually reported in the
# paper -- the depression natural indirect effect -- or whether it washes
# out under the averaging involved in computing NIE from the posterior
# draws (as it did for the analogous raw-coefficient/ATE distinction
# already documented for BMEOP).

dep_nie_list <- lapply(1:4, function(chain_i) {
  mh_result$results$results$chains[[chain_i]]$indirect_effects_individual[, 1]
})

dep_nie_mcmc <- mcmc.list(lapply(dep_nie_list, mcmc))

cat("\n=== Depression NIE (reported estimand) R-hat ===\n")
print(gelman.diag(dep_nie_mcmc))

cat("\n=== Depression NIE (reported estimand) ESS ===\n")
print(effectiveSize(dep_nie_mcmc))