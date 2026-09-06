# =============================================================================
# bjcm_intercept_wave_check.R
# =============================================================================
# BMEOP showed a near-exact additive non-identifiability between its
# intercept and the six wave-level random effects U_w (r ~ -0.999), with
# the raw intercept poorly mixed (R-hat=2.46) while the identified
# combination Intercept + mean(U_w) converged cleanly. BJCM has the same
# structural ingredient -- wave-level random effects on both the outcome
# side (U_out_wave) and each mediator side (U_med_wave) -- but this has
# never been checked. This script checks both.
#
# Read-only against the already-saved MH-corrected full run -- no re-fitting.
#
# Run from MH_BJCM/scripts/.
# =============================================================================

library(coda)

source("00_background_run_utils.R")
source("02_bjcm_model.R")

elsa_long <- read.csv("../data/elsa_longitudinal_analysis.csv")
bjcm_full <- readRDS("../results/bjcm_full_results.rds")
bjcm_raw  <- bjcm_full$results$results
chains    <- bjcm_raw$chains

inputs <- prepare_bjcm_inputs(elsa_long)
baseline_names <- colnames(inputs$X_baseline)
cat("X_baseline columns:", paste(baseline_names, collapse = ", "), "\n\n")

intercept_idx <- which(baseline_names %in% c("intercept", "(Intercept)", "Intercept"))
if (length(intercept_idx) == 0) {
  stop("Could not find an intercept column in X_baseline -- check baseline_names above ",
       "and set intercept_idx manually.")
}
cat("Intercept column index:", intercept_idx, "\n\n")

mediator_names <- c("Depression", "Transport mobility", "Mobility limitations")
n_mediators <- 3

# =============================================================================
# PART 1 -- Outcome model: gamma_out intercept vs. U_out_wave
# =============================================================================
cat("=============================================================\n")
cat("PART 1: Outcome model intercept (gamma_out) vs. U_out_wave\n")
cat("=============================================================\n\n")

gamma_out_all <- do.call(rbind, lapply(chains, function(ch) ch$gamma_out_samples))
colnames(gamma_out_all) <- baseline_names
U_out_wave_all <- do.call(rbind, lapply(chains, function(ch) ch$U_out_wave))
chain_len <- nrow(chains[[1]]$gamma_out_samples)

combined_out <- cbind(intercept = gamma_out_all[, intercept_idx], U_out_wave_all)
cor_out <- cor(combined_out)
cat("Correlation of outcome intercept with each U_out_wave column:\n")
print(round(cor_out["intercept", ], 3))

intercept_out_list <- lapply(1:4, function(i) chains[[i]]$gamma_out_samples[, intercept_idx])
intercept_out_mcmc <- mcmc.list(lapply(intercept_out_list, mcmc))
cat(sprintf("\nOutcome intercept alone:      ESS = %8.1f   R-hat = %.3f\n",
            effectiveSize(intercept_out_mcmc), gelman.diag(intercept_out_mcmc)$psrf[1]))

intercept_out_plus_remean <- gamma_out_all[, intercept_idx] + rowMeans(U_out_wave_all)
combined_out_list <- lapply(1:4, function(i) {
  idx <- ((i - 1) * chain_len + 1):(i * chain_len)
  mcmc(intercept_out_plus_remean[idx])
})
combined_out_mcmc <- mcmc.list(combined_out_list)
cat(sprintf("Outcome intercept + mean(U_w): ESS = %8.1f   R-hat = %.3f\n",
            effectiveSize(combined_out_mcmc), gelman.diag(combined_out_mcmc)$psrf[1]))

# =============================================================================
# PART 2 -- Mediator models: each mediator's own intercept vs. its own U_med_wave
# =============================================================================
cat("\n\n=============================================================\n")
cat("PART 2: Mediator model intercepts (beta_med) vs. U_med_wave\n")
cat("=============================================================\n\n")

for (k in 1:n_mediators) {
  
  beta_med_all_k <- do.call(rbind, lapply(chains, function(ch) ch$beta_med_samples[, , k]))
  colnames(beta_med_all_k) <- baseline_names
  U_med_wave_all_k <- do.call(rbind, lapply(chains, function(ch) ch$U_med_wave[, , k]))
  
  combined_med <- cbind(intercept = beta_med_all_k[, intercept_idx], U_med_wave_all_k)
  cor_med <- cor(combined_med)
  
  cat(sprintf("--- %s ---\n", mediator_names[k]))
  cat("Correlation of mediator intercept with each U_med_wave column:\n")
  print(round(cor_med["intercept", ], 3))
  
  intercept_med_list <- lapply(1:4, function(i) chains[[i]]$beta_med_samples[, intercept_idx, k])
  intercept_med_mcmc <- mcmc.list(lapply(intercept_med_list, mcmc))
  cat(sprintf("Mediator intercept alone:      ESS = %8.1f   R-hat = %.3f\n",
              effectiveSize(intercept_med_mcmc), gelman.diag(intercept_med_mcmc)$psrf[1]))
  
  intercept_med_plus_remean <- beta_med_all_k[, intercept_idx] + rowMeans(U_med_wave_all_k)
  combined_med_list <- lapply(1:4, function(i) {
    idx <- ((i - 1) * chain_len + 1):(i * chain_len)
    mcmc(intercept_med_plus_remean[idx])
  })
  combined_med_mcmc <- mcmc.list(combined_med_list)
  cat(sprintf("Mediator intercept + mean(U_w): ESS = %8.1f   R-hat = %.3f\n\n",
              effectiveSize(combined_med_mcmc), gelman.diag(combined_med_mcmc)$psrf[1]))
}

cat("=============================================================\n")
cat("Interpretation: if any 'intercept alone' R-hat is high (>1.1) while\n")
cat("the corresponding 'intercept + mean(U_w)' R-hat is clean (~1.00-1.01),\n")
cat("that model shares BMEOP's additive non-identifiability. If all\n")
cat("R-hats above are already clean, no such issue is present in BJCM.\n")