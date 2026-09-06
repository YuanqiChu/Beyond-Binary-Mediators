# ==============================================================================
# 05b_mediator_correlation_sensitivity.R
#
# Follow-up to 05a_mediator_independence_check.R. That diagnostic found a
# modest residual correlation between Depression and Mobility limitations
# (r = 0.165, posterior means), with the other two mediator pairs negligible
# (r = 0.045, r = 0.034). This script asks: if the true residual correlation
# between Depression and Mobility limitations were as high as (or higher
# than) what was measured, would the headline finding -- transport mobility
# dominating the NIE decomposition (~72%) -- still hold?
#
# Method: reuse simulate_parallel_mediation_effects()'s logic (Appendix A.2,
# Step 9), but replace the independent epsilon_med draw with a correlated
# draw between Depression and Mobility limitations only (Transport stays
# independent of both, consistent with the negligible correlations found in
# 05a_mediator_independence_check.R). Posterior mean parameters are reused
# as-is from bjcm_full_results.rds -- no re-fitting, no MCMC.
#
# Each mediator's path-specific NIE_k is computed as a Shapley value -- its
# marginal contribution averaged over all 3! = 6 orderings in which the
# mediators could be switched from control to treatment (Appendix A.2, Step
# 9(d)) -- matching simulate_parallel_mediation_effects() in
# 02_bjcm_model.R and removing any dependence on an arbitrary mediator
# ordering. TE, NDE, and NIE_total are unaffected, since they never involve
# an ordering choice.
#
# Runtime: seconds per rho value (deterministic simulation on posterior
# means, no MCMC).
# ==============================================================================

# ---- 0. Load what's needed (re-run 05a_mediator_independence_check.R first,
#          or re-derive the posterior means here if starting fresh) ---------
library(MASS)
library(mvtnorm)
library(coda)
library(dplyr)
library(parallel)
if (requireNamespace("abind", quietly = TRUE)) library(abind)

source("00_background_run_utils.R")
source("02_bjcm_model.R")

elsa_long <- read.csv("../data/elsa_longitudinal_analysis.csv")
bjcm_full <- readRDS("../results/bjcm_full_results.rds")
bjcm_raw  <- bjcm_full$results$results
chains    <- bjcm_raw$chains

inputs <- prepare_bjcm_inputs(elsa_long)
X_baseline <- inputs$X_baseline
treatment  <- inputs$treatment
wave_raw   <- inputs$wave
unique_waves <- sort(unique(wave_raw))
wave_mapped  <- match(wave_raw, unique_waves)

n_mediators <- 3
mediator_names <- c("Depression", "Transport mobility", "Mobility limitations")

# ---- 1. Pool posterior means across chains (mediator side, as in script 11) -
beta_med_mean <- apply(
  simplify2array(lapply(chains, function(ch) apply(ch$beta_med_samples, c(2, 3), mean))),
  c(1, 2), mean
)
alpha_med_mean <- rowMeans(sapply(chains, function(ch) colMeans(ch$alpha_med_samples)))
U_med_wave_mean <- apply(
  simplify2array(lapply(chains, function(ch) apply(ch$U_med_wave, c(2, 3), mean))),
  c(1, 2), mean
)
thresholds_med_mean <- lapply(1:n_mediators, function(k) {
  per_chain_means <- lapply(chains, function(ch) colMeans(ch$thresholds_med[[k]]))
  Reduce(`+`, per_chain_means) / length(per_chain_means)
})

# ---- 2. Pool posterior means across chains (outcome side) ------------------
beta_out_med_mean <- colMeans(do.call(rbind, lapply(chains, function(ch) colMeans(ch$beta_out_med_samples))))
lambda_out_mean   <- mean(sapply(chains, function(ch) mean(ch$lambda_out_samples)))
gamma_out_mean    <- colMeans(do.call(rbind, lapply(chains, function(ch) colMeans(ch$gamma_out_samples))))
thresholds_out_mean <- colMeans(do.call(rbind, lapply(chains, function(ch) colMeans(ch$thresholds_out))))
U_out_wave_mean   <- rowMeans(sapply(chains, function(ch) colMeans(ch$U_out_wave)))

cat("Posterior means pooled across", length(chains), "chains (mediator + outcome side).\n\n")

# ---- 3. Modified simulation: correlated Depression<->Mobility residuals ----
# Mediator order is [Depression, Transport, Mobility] (columns 1, 2, 3).
# We correlate columns 1 and 3 only, per the 05a_mediator_independence_check.R
# finding (Depression-Transport r=0.045, Transport-Mobility r=0.034, both
# negligible; Depression-Mobility r=0.165 is the one being stress-tested).
#
# NIE_k is now the Shapley value of mediator k's contribution: averaged over
# all orderings in which the K mediators could be switched from control to
# treatment (equivalently, averaged over all subsets S not containing k,
# weighted by |S|!(K-|S|-1)!/K!). This removes the dependence on any single
# fixed ordering and keeps the attribution consistent with the model's
# parallel (non-ordered) mediator structure.
simulate_with_rho <- function(rho_dep_mob, n_draws = 200, seed = 2026) {
  set.seed(seed)
  n <- nrow(X_baseline)
  
  Sigma <- diag(n_mediators)
  Sigma[1, 3] <- rho_dep_mob
  Sigma[3, 1] <- rho_dep_mob
  
  te_matrix   <- matrix(0, n, n_draws)
  nde_matrix  <- matrix(0, n, n_draws)
  nie_total_matrix <- matrix(0, n, n_draws)
  nie_individual_array <- array(0, dim = c(n, n_mediators, n_draws))
  
  S_out <- c(-Inf, thresholds_out_mean, Inf)
  
  # Shapley subset weights and subset table (K = n_mediators)
  K <- n_mediators
  shapley_weight <- sapply(0:(K - 1), function(s) {
    factorial(s) * factorial(K - s - 1) / factorial(K)
  })
  subset_membership <- as.matrix(expand.grid(rep(list(c(0, 1)), K)))
  colnames(subset_membership) <- NULL
  n_subsets <- nrow(subset_membership)
  subset_size <- rowSums(subset_membership)
  
  for (d in 1:n_draws) {
    # Correlated draw replaces the independent rnorm() draw in the original
    epsilon_med <- mvtnorm::rmvnorm(n, mean = rep(0, n_mediators), sigma = Sigma)
    epsilon_out <- rnorm(n)
    
    Z_med_treat_matrix   <- matrix(0, n, n_mediators)
    Z_med_control_matrix <- matrix(0, n, n_mediators)
    
    for (k in 1:n_mediators) {
      mu_med_treat_k <- as.vector(X_baseline %*% beta_med_mean[, k]) +
        1 * alpha_med_mean[k] + U_med_wave_mean[wave_mapped, k]
      Z_med_treat_matrix[, k] <- mu_med_treat_k + epsilon_med[, k]
      
      mu_med_control_k <- as.vector(X_baseline %*% beta_med_mean[, k]) +
        0 * alpha_med_mean[k] + U_med_wave_mean[wave_mapped, k]
      Z_med_control_matrix[, k] <- mu_med_control_k + epsilon_med[, k]
    }
    
    mediator_effects_11 <- rowSums(Z_med_treat_matrix *
                                     matrix(beta_out_med_mean, n, n_mediators, byrow = TRUE))
    mu_Y_11 <- mediator_effects_11 + 1 * lambda_out_mean +
      as.vector(X_baseline %*% gamma_out_mean) + U_out_wave_mean[wave_mapped]
    Z_Y_11 <- mu_Y_11 + epsilon_out
    
    mediator_effects_00 <- rowSums(Z_med_control_matrix *
                                     matrix(beta_out_med_mean, n, n_mediators, byrow = TRUE))
    mu_Y_00 <- mediator_effects_00 + 0 * lambda_out_mean +
      as.vector(X_baseline %*% gamma_out_mean) + U_out_wave_mean[wave_mapped]
    Z_Y_00 <- mu_Y_00 + epsilon_out
    
    mu_Y_10 <- mediator_effects_00 + 1 * lambda_out_mean +
      as.vector(X_baseline %*% gamma_out_mean) + U_out_wave_mean[wave_mapped]
    Z_Y_10 <- mu_Y_10 + epsilon_out
    
    Y_11 <- findInterval(Z_Y_11, S_out)
    Y_00 <- findInterval(Z_Y_00, S_out)
    Y_10 <- findInterval(Z_Y_10, S_out)
    
    te_matrix[, d]  <- Y_11 - Y_00
    nde_matrix[, d] <- Y_10 - Y_00
    nie_total_matrix[, d] <- Y_11 - Y_10
    
    # --- Shapley-symmetrised individual NIE_k -------------------------------
    # Step 1: evaluate v(S) = Y*(T=1, M_S=treated, M_{-S}=control) for all
    # 2^K subsets S (v(empty)=Y_10, v(full)=Y_11, already computed above).
    v_S <- matrix(0, n, n_subsets)
    for (j in 1:n_subsets) {
      membership_j <- subset_membership[j, ]
      if (all(membership_j == 0)) {
        v_S[, j] <- Y_10
      } else if (all(membership_j == 1)) {
        v_S[, j] <- Y_11
      } else {
        Z_med_S <- Z_med_control_matrix
        treated_cols <- which(membership_j == 1)
        Z_med_S[, treated_cols] <- Z_med_treat_matrix[, treated_cols]
        mediator_effects_S <- rowSums(Z_med_S *
                                        matrix(beta_out_med_mean, n, n_mediators, byrow = TRUE))
        mu_Y_S <- mediator_effects_S + 1 * lambda_out_mean +
          as.vector(X_baseline %*% gamma_out_mean) + U_out_wave_mean[wave_mapped]
        Z_Y_S <- mu_Y_S + epsilon_out
        v_S[, j] <- findInterval(Z_Y_S, S_out)
      }
    }
    
    # Step 2: for each mediator k, average the marginal contribution
    # v(S union {k}) - v(S) over all subsets S not containing k, weighted by
    # the Shapley coalition-size weight.
    for (k in 1:n_mediators) {
      contribution_k <- numeric(n)
      subsets_without_k <- which(subset_membership[, k] == 0)
      for (j in subsets_without_k) {
        membership_j_plus_k <- subset_membership[j, ]
        membership_j_plus_k[k] <- 1
        j_plus_k <- which(apply(subset_membership, 1, function(r) all(r == membership_j_plus_k)))
        w <- shapley_weight[subset_size[j] + 1]
        contribution_k <- contribution_k + w * (v_S[, j_plus_k] - v_S[, j])
      }
      nie_individual_array[, k, d] <- contribution_k
    }
  }
  
  te_mean  <- mean(rowMeans(te_matrix))
  nde_mean <- mean(rowMeans(nde_matrix))
  nie_mean <- mean(rowMeans(nie_total_matrix))
  nie_k_means <- apply(apply(nie_individual_array, c(1, 2), mean), 2, mean)
  
  list(
    rho = rho_dep_mob,
    TE = te_mean, NDE = nde_mean, NIE = nie_mean,
    NIE_depression = nie_k_means[1],
    NIE_transport  = nie_k_means[2],
    NIE_mobility   = nie_k_means[3],
    pct_depression = nie_k_means[1] / sum(nie_k_means),
    pct_transport  = nie_k_means[2] / sum(nie_k_means),
    pct_mobility   = nie_k_means[3] / sum(nie_k_means)
  )
}

# ---- 4. Sweep rho from 0 up past the measured 0.165, with margin -----------
rho_grid <- c(0, 0.05, 0.10, 0.15, 0.196, 0.20, 0.25, 0.30, 0.40, 0.50)

cat("Running sensitivity sweep over rho_dep_mob =", paste(rho_grid, collapse = ", "), "\n")
cat("(n_draws = 200 per rho value, Shapley-symmetrised NIE_k - seconds each, no MCMC)\n\n")

results_list <- lapply(rho_grid, simulate_with_rho)
results_df <- do.call(rbind, lapply(results_list, as.data.frame))

# ---- 5. Report ---------------------------------------------------------------
cat("============ SENSITIVITY OF NIE MEDIATOR BREAKDOWN TO rho(Depression, Mobility) ============\n\n")
print_df <- results_df
print_df[, -1] <- round(print_df[, -1], 4)
print(print_df, row.names = FALSE)

cat("\nTransport mobility's share of NIE, as rho increases:\n")
for (i in seq_len(nrow(results_df))) {
  cat(sprintf("  rho=%.3f : Transport = %.1f%% (Depression = %.1f%%, Mobility = %.1f%%)\n",
              results_df$rho[i], 100 * results_df$pct_transport[i],
              100 * results_df$pct_depression[i], 100 * results_df$pct_mobility[i]))
}

max_transport_drop <- max(results_df$pct_transport) - min(results_df$pct_transport)
still_dominant <- all(results_df$pct_transport > pmax(results_df$pct_depression, results_df$pct_mobility))

cat("\nMax change in Transport's share across the whole rho grid:", round(100 * max_transport_drop, 2), "percentage points\n")
cat("Does Transport remain the dominant mediator at every rho tested?", still_dominant, "\n")

saveRDS(results_df, "../results/mediator_correlation_sensitivity.rds")
cat("\nSaved: ../results/mediator_correlation_sensitivity.rds\n")

# ==============================================================================
# 6. TE / NDE / NIE convergence diagnostics (Gelman-Rubin R-hat, ESS)
# ==============================================================================
# Unrelated to the rho sensitivity sweep above -- included here only because
# this script already has the MH-corrected full result object (bjcm_full)
# loaded in memory, so no extra readRDS() is needed. Confirms the three
# main estimands actually reported in the paper (Table 4, Panel B) converged
# cleanly across the four chains, independent of the raw threshold-level
# R-hat/ESS figures reported elsewhere (see diagnose_bjcm_convergence.R).
#
# total_effects / direct_effects / indirect_effects are already
# population-average scalars per stored iteration (not per-individual), so
# no rowMeans() collapse is needed before building the mcmc.list.
# ==============================================================================

cat("\n\n============ TE / NDE / NIE convergence diagnostics (MH-corrected) ============\n\n")

te_list  <- lapply(chains, function(ch) ch$total_effects)
nde_list <- lapply(chains, function(ch) ch$direct_effects)
nie_list <- lapply(chains, function(ch) ch$indirect_effects)

te_mcmc  <- mcmc.list(lapply(te_list,  mcmc))
nde_mcmc <- mcmc.list(lapply(nde_list, mcmc))
nie_mcmc <- mcmc.list(lapply(nie_list, mcmc))

cat("TE:  ESS =", effectiveSize(te_mcmc),  " R-hat =", gelman.diag(te_mcmc)$psrf[1], "\n")
cat("NDE: ESS =", effectiveSize(nde_mcmc), " R-hat =", gelman.diag(nde_mcmc)$psrf[1], "\n")
cat("NIE: ESS =", effectiveSize(nie_mcmc), " R-hat =", gelman.diag(nie_mcmc)$psrf[1], "\n")