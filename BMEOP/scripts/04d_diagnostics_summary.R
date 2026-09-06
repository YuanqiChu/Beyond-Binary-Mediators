# =============================================================================
# 04d_diagnostics_summary.R
# =============================================================================
# Read-only diagnostics against the fitted BMEOP full run: coefficient and
# ATE summaries, per-coefficient convergence, and the intercept / wave-
# level random-effect non-identifiability check (Section 4.2). No model
# fitting happens in this script.
#
# Run from scripts/, so that ../results/ resolves correctly.
# =============================================================================

library(coda)

mh_bmeop <- readRDS("../results/bmeop_full_results.rds")

covariate_names <- rownames(mh_bmeop$beta_summary)

outcome_sd <- 1.634

# =============================================================================
# PART 1 -- Coefficient and ATE summaries
# =============================================================================

cat("=============================================================\n")
cat("PART 1: Coefficient and ATE summaries\n")
cat("=============================================================\n\n")

cat("--- beta_summary (all 10 covariates) ---\n")
print(mh_bmeop$beta_summary)

cat("\n--- sigma2_u summary ---\n")
sigma2_u_mean <- mean(mh_bmeop$sigma2_u_samples)
sigma2_u_sd   <- sd(mh_bmeop$sigma2_u_samples)
sigma2_u_ci   <- quantile(mh_bmeop$sigma2_u_samples, c(0.025, 0.975))
cat(sprintf("mean=%.4f  sd=%.4f  95%% CrI=[%.4f, %.4f]\n",
            sigma2_u_mean, sigma2_u_sd, sigma2_u_ci[1], sigma2_u_ci[2]))

cat("\n--- ATE full summary ---\n")
ate_mean <- mean(mh_bmeop$ate_samples)
ate_sd   <- sd(mh_bmeop$ate_samples)
ate_ci   <- quantile(mh_bmeop$ate_samples, c(0.025, 0.975))
cat(sprintf("mean=%.4f  sd=%.4f  95%% CrI=[%.4f, %.4f]  Cohen's d=%.4f\n",
            ate_mean, ate_sd, ate_ci[1], ate_ci[2], ate_mean / outcome_sd))


# =============================================================================
# PART 2 -- Per-coefficient convergence diagnostics (previously flagged ones)
# =============================================================================

cat("\n\n=============================================================\n")
cat("PART 2: Per-coefficient convergence (depression, intercept, sclife, treatment)\n")
cat("=============================================================\n\n")

beta_all <- do.call(rbind, lapply(1:4, function(i) mh_bmeop$chains[[i]]$beta))
colnames(beta_all) <- covariate_names
chain_len <- nrow(mh_bmeop$chains[[1]]$beta)

target_coefs <- c("depression", "intercept", "sclife", "livalone_TREATMENT")

for (coef_name in target_coefs) {
  coef_idx <- which(covariate_names == coef_name)
  chain_list <- lapply(1:4, function(i) mh_bmeop$chains[[i]]$beta[, coef_idx])
  mcmc_obj <- mcmc.list(lapply(chain_list, mcmc))
  cat(sprintf("%-20s ESS = %8.1f   R-hat = %.3f\n",
              coef_name, effectiveSize(mcmc_obj), gelman.diag(mcmc_obj)$psrf[1]))
}


# =============================================================================
# PART 3 -- Intercept / wave-random-effect non-identifiability check
# =============================================================================

cat("\n\n=============================================================\n")
cat("PART 3: Intercept vs. wave-level random effects\n")
cat("=============================================================\n\n")

cor_matrix <- cor(beta_all)
cat("--- Correlation of intercept with all other beta coefficients ---\n")
print(round(cor_matrix["intercept", ], 3))

random_effects_all <- do.call(rbind, lapply(1:4, function(i) mh_bmeop$chains[[i]]$random_effects))
combined <- cbind(intercept = beta_all[, "intercept"], random_effects_all)
cor_re <- cor(combined)
cat("\n--- Correlation of intercept with each wave-level random effect U_w ---\n")
print(round(cor_re["intercept", ], 3))

intercept_plus_remean <- beta_all[, "intercept"] + rowMeans(random_effects_all)
combined_list <- lapply(1:4, function(i) {
  idx <- ((i - 1) * chain_len + 1):(i * chain_len)
  mcmc(intercept_plus_remean[idx])
})
combined_mcmc <- mcmc.list(combined_list)

cat(sprintf("\nIntercept + mean(U_w): ESS = %.1f   R-hat = %.3f\n",
            effectiveSize(combined_mcmc), gelman.diag(combined_mcmc)$psrf[1]))

intercept_plus_remean_summary <- c(
  mean = mean(intercept_plus_remean),
  sd   = sd(intercept_plus_remean),
  quantile(intercept_plus_remean, c(0.025, 0.975))
)
cat("\n--- Intercept + mean(U_w): summary ---\n")
print(round(intercept_plus_remean_summary, 4))

cat("\n--- For comparison: intercept alone (unreliable per R-hat above) ---\n")
cat(sprintf("mean=%.4f  sd=%.4f  95%% CrI=[%.4f, %.4f]\n",
            mean(beta_all[, "intercept"]), sd(beta_all[, "intercept"]),
            quantile(beta_all[, "intercept"], 0.025),
            quantile(beta_all[, "intercept"], 0.975)))