# =============================================================================
# table_b2_bmeop_remaining_coefficients.R
# =============================================================================
# Appendix Table B.2: per-coefficient convergence diagnostics for the BMEOP
# covariates NOT already covered by bmeop_mh_diagnostics_summary.R (which
# handles depression, sclife, treatment, and the intercept / intercept +
# mean(U_w) pair, but only reports their point-estimate R-hat, not the
# upper bound). This script adds: (1) the remaining six covariates,
# sigma2_u, and the five outcome thresholds -- full point + upper bound;
# (2) the missing upper-bound R-hat for depression, sclife, treatment, and
# Intercept + mean(U_w), so Table B.2 can be assembled fully complete.
#
# Read-only against the already-saved MH-corrected full run -- no re-fitting.
#
# Run from MH_BMEOP/scripts/. Assumes the same data/logs/results/scripts
# layout as MH_BJCM -- adjust the path below if bmeop_full_results.rds
# is not under ../results/.
# =============================================================================

library(coda)

mh_bmeop <- readRDS("../results/bmeop_full_results.rds")

covariate_names <- rownames(mh_bmeop$beta_summary)
print(covariate_names)  # sanity check: confirms column order before indexing

beta_all  <- do.call(rbind, lapply(1:4, function(i) mh_bmeop$chains[[i]]$beta))
colnames(beta_all) <- covariate_names
chain_len <- nrow(mh_bmeop$chains[[1]]$beta)

cat("\n=============================================================\n")
cat("TABLE B.2 -- BMEOP remaining coefficients\n")
cat("=============================================================\n\n")

remaining_coefs <- c("age_gr", "dhsex2", "edqual2", "health",
                     "transport_mobility", "mobility_limitations")

for (coef_name in remaining_coefs) {
  coef_idx <- which(covariate_names == coef_name)
  chain_list <- lapply(1:4, function(i) mh_bmeop$chains[[i]]$beta[, coef_idx])
  mcmc_obj <- mcmc.list(lapply(chain_list, mcmc))
  diag_obj <- gelman.diag(mcmc_obj)
  cat(sprintf("%-22s ESS = %8.1f   R-hat = %.3f   R-hat (UB) = %.3f\n",
              coef_name, effectiveSize(mcmc_obj),
              diag_obj$psrf[1, 1], diag_obj$psrf[1, 2]))
}

# sigma2_u -- stored per-chain as a flat vector, not part of the beta matrix
sigma2u_list <- lapply(1:4, function(i) mh_bmeop$chains[[i]]$sigma2_u)
sigma2u_mcmc <- mcmc.list(lapply(sigma2u_list, mcmc))
sigma2u_diag <- gelman.diag(sigma2u_mcmc)
cat(sprintf("%-22s ESS = %8.1f   R-hat = %.3f   R-hat (UB) = %.3f\n",
            "sigma2_u", effectiveSize(sigma2u_mcmc),
            sigma2u_diag$psrf[1, 1], sigma2u_diag$psrf[1, 2]))

# thresholds (5 of them, bmeop_outcome_j2..j6)
threshold_chains <- lapply(1:4, function(i) mcmc(mh_bmeop$chains[[i]]$thresholds[, 2:6]))
threshold_mcmc <- mcmc.list(threshold_chains)
threshold_diag <- gelman.diag(threshold_mcmc)
threshold_ess  <- effectiveSize(threshold_mcmc)

cat("\nThresholds (j2-j6):\n")
for (k in 1:5) {
  cat(sprintf("  j%d  ESS = %8.1f   R-hat = %.3f   R-hat (UB) = %.3f\n",
              k + 1, threshold_ess[k],
              threshold_diag$psrf[k, 1], threshold_diag$psrf[k, 2]))
}

# ---- Fill in the upper-bound column for the coefficients already reported
#      (point-estimate only) by bmeop_mh_diagnostics_summary.R -----------
cat("\n=============================================================\n")
cat("Upper-bound R-hat for depression/sclife/treatment/Intercept+U_w_mean\n")
cat("(point estimates for these already reported by bmeop_mh_diagnostics_summary.R)\n")
cat("=============================================================\n\n")

target_coefs <- c("depression", "sclife", "livalone_TREATMENT")
for (coef_name in target_coefs) {
  coef_idx <- which(covariate_names == coef_name)
  chain_list <- lapply(1:4, function(i) mh_bmeop$chains[[i]]$beta[, coef_idx])
  mcmc_obj <- mcmc.list(lapply(chain_list, mcmc))
  diag_obj <- gelman.diag(mcmc_obj)
  cat(sprintf("%-20s R-hat = %.3f   R-hat (UB) = %.3f\n",
              coef_name, diag_obj$psrf[1, 1], diag_obj$psrf[1, 2]))
}

random_effects_all <- do.call(rbind, lapply(1:4, function(i) mh_bmeop$chains[[i]]$random_effects))
intercept_plus_remean <- beta_all[, "intercept"] + rowMeans(random_effects_all)
combined_list <- lapply(1:4, function(i) {
  idx <- ((i - 1) * chain_len + 1):(i * chain_len)
  mcmc(intercept_plus_remean[idx])
})
combined_mcmc <- mcmc.list(combined_list)
combined_diag <- gelman.diag(combined_mcmc)
cat(sprintf("%-20s R-hat = %.3f   R-hat (UB) = %.3f\n",
            "Intercept+U_w_mean", combined_diag$psrf[1, 1], combined_diag$psrf[1, 2]))

cat("\nNOTE: intercept alone (R-hat=2.463) is not re-checked here for UB --\n")
cat("its point estimate already far exceeds any convergence threshold, so\n")
cat("the upper bound is not informative and re-running it is not useful.\n")
cat("Full point estimates for depression/sclife/livalone_TREATMENT/intercept/\n")
cat("Intercept+mean(U_w) are reported by bmeop_mh_diagnostics_summary.R --\n")
cat("combine that script's Part 2/3 output with this script's output to\n")
cat("assemble the complete Table B.2.\n")