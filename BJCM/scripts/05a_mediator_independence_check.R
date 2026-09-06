# ==============================================================================
# 05a_mediator_independence_check.R
#
# The parallel mediation model treats the three mediators as conditionally
# independent given treatment and baseline covariates (Section 3.4). This
# script checks that assumption empirically.
#
# This script does NOT re-fit anything. It uses the already-fitted posterior
# means from bjcm_full_results.rds (beta_med, alpha_med, U_med_wave,
# thresholds_med) plus the analysis dataset to compute, for each person,
# a continuous "randomized quantile residual" (Dunn & Smyth, 1996) for each
# of the three mediators - the standard diagnostic for checking whether
# residuals from ordinal models are correlated after conditioning on
# treatment + covariates + wave effects (i.e. whether the conditional
# independence assumption is a reasonable approximation).
#
# If the model is correctly specified AND mediators are conditionally
# independent given (treatment, X_baseline, wave), these residuals should be
# ~ N(0,1) and uncorrelated across mediators. Non-trivial residual
# correlation would be direct empirical evidence against the independence
# assumption -- for instance, if transport difficulties causally precede
# depression within the inter-wave interval.
#
# Runtime: seconds, not hours. No MCMC involved.
# ==============================================================================

library(MASS)
library(mvtnorm)
library(coda)
library(dplyr)
library(parallel)      
if (requireNamespace("abind", quietly = TRUE)) library(abind)   

# ---- 0. Load what's needed --------------------------------------------------
source("00_background_run_utils.R")
source("02_bjcm_model.R")

elsa_long <- read.csv("../data/elsa_longitudinal_analysis.csv")
bjcm_full <- readRDS("../results/bjcm_full_results.rds")
bjcm_raw  <- bjcm_full$results$results

# ---- 1. Reconstruct the exact design matrices used for fitting -------------
# prepare_bjcm_inputs() is the same function used to build the model inputs,
# so this reconstruction is guaranteed consistent with what was actually fit.
inputs <- prepare_bjcm_inputs(elsa_long)

X_baseline <- inputs$X_baseline
treatment  <- inputs$treatment
wave_raw   <- inputs$wave
mediators_raw <- inputs$mediators   # Depression_lag, Transport_lag, Mobility_lag (raw coded)

n_mediators <- ncol(mediators_raw)
mediator_names <- c("Depression", "Transport mobility", "Mobility limitations")

# Recode each mediator to 1..K exactly as the model does internally
mediators_recoded <- apply(mediators_raw, 2, function(col) {
  match(col, sort(unique(col)))
})

unique_waves <- sort(unique(wave_raw))
wave_mapped  <- match(wave_raw, unique_waves)
n_waves      <- length(unique_waves)

cat("Reconstructed inputs: n =", nrow(X_baseline), "| n_waves =", n_waves, "\n")

# ---- 2. Pool posterior MEANS across all 4 chains for the mediator params --
# beta_med_samples: [n_store, p_baseline, n_mediators] per chain
# alpha_med_samples: [n_store, n_mediators] per chain
# U_med_wave_samples: [n_store, n_waves, n_mediators] per chain
# thresholds_med[[k]]: [n_store, n_thresh_k] per chain

chains <- bjcm_raw$chains

beta_med_mean <- apply(
  simplify2array(lapply(chains, function(ch) apply(ch$beta_med_samples, c(2, 3), mean))),
  c(1, 2), mean
)  # p_baseline x n_mediators

alpha_med_mean <- rowMeans(sapply(chains, function(ch) colMeans(ch$alpha_med_samples)))  # n_mediators

U_med_wave_mean <- apply(
  simplify2array(lapply(chains, function(ch) apply(ch$U_med_wave, c(2, 3), mean))),
  c(1, 2), mean
)  # n_waves x n_mediators

thresholds_med_mean <- lapply(1:n_mediators, function(k) {
  # NOTE: sapply()+rowMeans() breaks when a mediator has only 1 threshold
  # (binary transport mobility) - colMeans() then returns a scalar per chain,
  # so sapply() yields a plain length-4 vector instead of a matrix, and
  # rowMeans() has nothing to reduce over. Reduce(`+`, ...) works either way.
  per_chain_means <- lapply(chains, function(ch) colMeans(ch$thresholds_med[[k]]))
  Reduce(`+`, per_chain_means) / length(per_chain_means)
})

cat("Posterior means pooled across", length(chains), "chains.\n\n")

# ---- 3. Compute the latent linear predictor mu_k for each mediator ---------
mu <- matrix(NA, nrow(X_baseline), n_mediators)
for (k in 1:n_mediators) {
  mu[, k] <- as.vector(X_baseline %*% beta_med_mean[, k]) +
    treatment * alpha_med_mean[k] +
    U_med_wave_mean[wave_mapped, k]
}

# ---- 4. Randomized quantile (surrogate) residuals --------------------------
# For observed category y in {1..K}, with thresholds S = c(-Inf, tau_1..tau_{K-1}, Inf):
#   lower = Phi(S[y]   - mu)
#   upper = Phi(S[y+1] - mu)
#   u ~ Uniform(lower, upper)
#   residual = Phi^-1(u)
# Repeated R_REPS times and averaged to reduce Monte Carlo noise from the
# randomization step itself (the noise this introduces is unrelated to MCMC
# and settles down fast with more replicates - 30 is more than enough here).

set.seed(2026)
R_REPS <- 30

compute_residuals_once <- function() {
  resid_mat <- matrix(NA, nrow(X_baseline), n_mediators)
  for (k in 1:n_mediators) {
    y_k <- mediators_recoded[, k]
    S_k <- c(-Inf, thresholds_med_mean[[k]], Inf)
    valid <- !is.na(y_k)
    lower <- pnorm(S_k[y_k[valid]]     - mu[valid, k])
    upper <- pnorm(S_k[y_k[valid] + 1] - mu[valid, k])
    u <- runif(sum(valid), pmin(lower, upper), pmax(lower, upper))
    u <- pmin(pmax(u, 1e-10), 1 - 1e-10)  # avoid +-Inf from qnorm at the boundary
    resid_mat[valid, k] <- qnorm(u)
  }
  resid_mat
}

cor_reps <- array(NA, dim = c(n_mediators, n_mediators, R_REPS))
for (r in 1:R_REPS) {
  resid_mat <- compute_residuals_once()
  cor_reps[, , r] <- cor(resid_mat, use = "pairwise.complete.obs")
}

cor_mean <- apply(cor_reps, c(1, 2), mean)
cor_sd   <- apply(cor_reps, c(1, 2), sd)

dimnames(cor_mean) <- list(mediator_names, mediator_names)
dimnames(cor_sd)   <- list(mediator_names, mediator_names)

# ---- 5. Report ---------------------------------------------------------------
cat("============ CONDITIONAL RESIDUAL CORRELATION (mediator independence check) ============\n")
cat("Averaged over", R_REPS, "randomized-quantile-residual replicates.\n\n")
cat("Mean correlation matrix:\n")
print(round(cor_mean, 3))
cat("\nSD across replicates (Monte Carlo noise from randomization, not uncertainty in the estimate):\n")
print(round(cor_sd, 4))

off_diag <- cor_mean[upper.tri(cor_mean)]
cat("\nOff-diagonal pairwise correlations:\n")
pair_names <- combn(mediator_names, 2, FUN = function(x) paste(x, collapse = " vs "))
for (i in seq_along(off_diag)) {
  cat(" ", pair_names[i], ":", round(off_diag[i], 3), "\n")
}
cat("\nMax |off-diagonal correlation|:", round(max(abs(off_diag)), 3), "\n")

cat("\nInterpretation guide:\n")
cat("  |r| < 0.10 : negligible residual dependence - independence assumption well supported empirically\n")
cat("  |r| 0.10-0.20 : modest dependence - worth flagging as a quantified limitation\n")
cat("  |r| > 0.20 : non-trivial dependence - independence assumption should be revisited\n")

saveRDS(list(cor_mean = cor_mean, cor_sd = cor_sd, cor_reps = cor_reps),
        "../results/mediator_independence_check.rds")
cat("\nSaved: ../results/mediator_independence_check.rds\n")
