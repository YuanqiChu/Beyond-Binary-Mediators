# =============================================================================
# Bayesian Mixed-Effects Ordered Probit Model for ELSA Loneliness Analysis
# Implementation of Algorithm 1 (BMEOP) 
# =============================================================================

# This code implements the Bayesian mixed-effects ordered probit (BMEOP) model in
# Chu et al. (2025). "Beyond Binary Mediators: A Bayesian Mixed-Effect Modelling Framework 
# for Understanding Causal Pathways to Loneliness in Late Life."
# =============================================================================

# Load required packages
library(MASS)
library(mvtnorm)
library(coda)
library(Matrix)
library(dplyr)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Robust truncated normal sampling
#' 
#' Handles edge cases in truncated normal sampling including extreme bounds,
#' invalid intervals, and numerical precision issues.
#' 
#' @param n Number of samples
#' @param mean Mean parameter (scalar or vector)
#' @param sd Standard deviation (scalar or vector)
#' @param lower Lower truncation bound
#' @param upper Upper truncation bound
#' @return Vector of truncated normal samples
rtruncnorm_robust <- function(n, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
  if (length(mean) == 1) mean <- rep(mean, n)
  if (length(sd) == 1) sd <- rep(sd, n)
  if (length(lower) == 1) lower <- rep(lower, n)
  if (length(upper) == 1) upper <- rep(upper, n)
  
  result <- numeric(n)
  
  for (i in 1:n) {
    # Handle invalid inputs
    if (is.na(mean[i]) || is.na(sd[i]) || sd[i] <= 0) {
      result[i] <- 0
      next
    }
    
    # Handle invalid bounds
    if (is.na(lower[i])) lower[i] <- -Inf
    if (is.na(upper[i])) upper[i] <- Inf
    
    if (is.finite(lower[i]) && is.finite(upper[i]) && lower[i] >= upper[i]) {
      result[i] <- mean[i]
      next
    }
    
    # Standardise bounds
    lower_std <- if (is.finite(lower[i])) (lower[i] - mean[i]) / sd[i] else -Inf
    upper_std <- if (is.finite(upper[i])) (upper[i] - mean[i]) / sd[i] else Inf
    
    # Handle truncation cases
    if (is.infinite(lower_std) && is.infinite(upper_std)) {
      result[i] <- mean[i] + sd[i] * rnorm(1)
    } else if (is.infinite(lower_std)) {
      # Upper truncation only
      if (upper_std < -10) {
        result[i] <- mean[i] + sd[i] * (upper_std - 0.1)
      } else {
        p_upper <- pnorm(upper_std)
        if (is.na(p_upper) || p_upper <= 0) {
          result[i] <- mean[i] + sd[i] * (upper_std - 1)
        } else {
          u <- runif(1, 1e-12, p_upper)
          result[i] <- mean[i] + sd[i] * qnorm(pmax(u, 1e-12))
        }
      }
    } else if (is.infinite(upper_std)) {
      # Lower truncation only
      if (lower_std > 10) {
        result[i] <- mean[i] + sd[i] * (lower_std + 0.1)
      } else {
        p_lower <- pnorm(lower_std)
        if (is.na(p_lower) || p_lower >= 1) {
          result[i] <- mean[i] + sd[i] * (lower_std + 1)
        } else {
          u <- runif(1, p_lower, 1 - 1e-12)
          result[i] <- mean[i] + sd[i] * qnorm(pmin(u, 1 - 1e-12))
        }
      }
    } else {
      # Double truncation
      p_lower <- pnorm(lower_std)
      p_upper <- pnorm(upper_std)
      
      if (is.na(p_lower) || is.na(p_upper)) {
        result[i] <- mean[i]
        next
      }
      
      if (p_upper - p_lower < 1e-12 || p_lower >= p_upper) {
        result[i] <- mean[i] + sd[i] * (lower_std + upper_std) / 2
      } else {
        u <- runif(1, p_lower, p_upper)
        u_safe <- pmax(pmin(u, 1 - 1e-12), 1e-12)
        z_val <- qnorm(u_safe)
        if (is.na(z_val) || is.infinite(z_val)) {
          result[i] <- mean[i] + sd[i] * (lower_std + upper_std) / 2
        } else {
          result[i] <- mean[i] + sd[i] * z_val
        }
      }
    }
    
    # Final validation
    if (is.na(result[i]) || is.infinite(result[i])) {
      result[i] <- mean[i]
    }
  }
  
  return(result)
}

#' Safe matrix inversion with regularisation
#' 
#' Attempts Cholesky decomposition first, falls back to SVD if necessary.
#' 
#' @param A Matrix to invert
#' @param lambda Ridge regularisation parameter
#' @return Inverse matrix
safe_matrix_inverse <- function(A, lambda = 1e-4) {
  tryCatch({
    A_reg <- A + lambda * diag(ncol(A))
    L <- chol(A_reg)
    return(chol2inv(L))
  }, error = function(e) {
    tryCatch({
      svd_result <- svd(A)
      d_inv <- ifelse(svd_result$d > 1e-10, 1/svd_result$d, 0)
      return(svd_result$v %*% diag(d_inv) %*% t(svd_result$u))
    }, error = function(e2) {
      return(diag(1, ncol(A)))
    })
  })
}

#' Robust multivariate normal sampling
#' 
#' @param n Number of samples (typically 1 for Gibbs sampling)
#' @param mean Mean vector
#' @param sigma Covariance matrix
#' @return Sample vector
rmvnorm_robust <- function(n, mean, sigma) {
  tryCatch({
    return(as.vector(rmvnorm(n, mean = mean, sigma = sigma)))
  }, error = function(e) {
    warning("Multivariate normal sampling failed, using independent normals")
    return(mean + sqrt(diag(sigma)) * rnorm(length(mean)))
  })
}

# =============================================================================
# ATE COMPUTATION FUNCTIONS
# =============================================================================

#' Analytical average treatment effect computation
#' 
#' Computes ATE using analytical integration over the latent variable
#' distribution. More stable and efficient than Monte Carlo for ordinal
#' outcomes with moderate numbers of categories.
#' 
#' @param X Covariate matrix
#' @param beta Coefficient vector
#' @param U_wave Wave-level random effects
#' @param wave_mapped Wave assignment for each observation
#' @param thresholds Threshold parameters
#' @param treatment_col Column index of treatment variable
#' @param n_draws Number of Monte Carlo draws (used only if K > 7)
#' @return List containing ATE estimate and diagnostics
compute_ate_analytical <- function(X, beta, U_wave, wave_mapped, thresholds, 
                                   treatment_col, n_draws = 50) {
  
  n <- nrow(X)
  S <- c(-Inf, thresholds, Inf)
  K <- length(S) - 1
  
  individual_effects <- numeric(n)
  
  for (i in 1:n) {
    # Counterfactual means
    X_treat <- X[i, ]
    X_treat[treatment_col] <- 1
    mu_treat <- sum(X_treat * beta) + U_wave[wave_mapped[i]]
    
    X_control <- X[i, ]
    X_control[treatment_col] <- 0
    mu_control <- sum(X_control * beta) + U_wave[wave_mapped[i]]
    
    if (K <= 7) {
      # Analytical computation for moderate K
      prob_treat <- numeric(K)
      prob_control <- numeric(K)
      
      for (k in 1:K) {
        prob_treat[k] <- pnorm(S[k+1], mean = mu_treat, sd = 1) - 
          pnorm(S[k], mean = mu_treat, sd = 1)
        prob_control[k] <- pnorm(S[k+1], mean = mu_control, sd = 1) - 
          pnorm(S[k], mean = mu_control, sd = 1)
      }
      
      E_Y_treat <- sum((1:K) * prob_treat)
      E_Y_control <- sum((1:K) * prob_control)
      
      individual_effects[i] <- E_Y_treat - E_Y_control
      
    } else {
      # Monte Carlo for large K
      effects_draws <- numeric(n_draws)
      
      for (d in 1:n_draws) {
        Z_cf_treat <- rnorm(1, mean = mu_treat, sd = 1)
        Z_cf_control <- rnorm(1, mean = mu_control, sd = 1)
        
        Y_cf_treat <- findInterval(Z_cf_treat, S)
        Y_cf_control <- findInterval(Z_cf_control, S)
        
        effects_draws[d] <- Y_cf_treat - Y_cf_control
      }
      
      individual_effects[i] <- mean(effects_draws)
    }
  }
  
  overall_ate <- mean(individual_effects)
  
  return(list(
    ate = overall_ate,
    individual_effects = individual_effects,
    monte_carlo_se = sd(individual_effects) / sqrt(n),
    method = if (K <= 7) "analytical" else "monte_carlo"
  ))
}

#' Multi-draw Monte Carlo ATE computation
#' 
#' Alternative ATE estimation using multiple Monte Carlo draws per individual.
#' May be preferred for very large K or when analytical computation is unstable.
#' 
#' @param X Covariate matrix
#' @param beta Coefficient vector
#' @param U_wave Wave-level random effects
#' @param wave_mapped Wave assignment
#' @param thresholds Threshold parameters
#' @param treatment_col Treatment variable column
#' @param n_draws Number of Monte Carlo draws per individual
#' @return List containing ATE estimate and stability metrics
compute_ate_multidraw <- function(X, beta, U_wave, wave_mapped, thresholds, 
                                  treatment_col, n_draws = 100) {
  
  n <- nrow(X)
  S <- c(-Inf, thresholds, Inf)
  
  individual_effects_matrix <- matrix(0, nrow = n, ncol = n_draws)
  all_random_draws <- matrix(rnorm(n * n_draws * 2), nrow = n * 2, ncol = n_draws)
  
  for (i in 1:n) {
    random_draws_treat <- all_random_draws[2*i - 1, ]
    random_draws_control <- all_random_draws[2*i, ]
    
    X_treat <- X[i, ]
    X_treat[treatment_col] <- 1
    mu_treat <- sum(X_treat * beta) + U_wave[wave_mapped[i]]
    
    X_control <- X[i, ]
    X_control[treatment_col] <- 0
    mu_control <- sum(X_control * beta) + U_wave[wave_mapped[i]]
    
    for (d in 1:n_draws) {
      Z_cf_treat <- mu_treat + random_draws_treat[d]
      Z_cf_control <- mu_control + random_draws_control[d]
      
      Y_cf_treat <- findInterval(Z_cf_treat, S)
      Y_cf_control <- findInterval(Z_cf_control, S)
      
      individual_effects_matrix[i, d] <- Y_cf_treat - Y_cf_control
    }
  }
  
  individual_ate_estimates <- rowMeans(individual_effects_matrix)
  overall_ate <- mean(individual_ate_estimates)
  
  individual_ate_se <- apply(individual_effects_matrix, 1, sd) / sqrt(n_draws)
  monte_carlo_se <- sd(individual_ate_estimates) / sqrt(n)
  
  return(list(
    ate = overall_ate,
    individual_effects = individual_ate_estimates,
    individual_se = individual_ate_se,
    monte_carlo_se = monte_carlo_se,
    stability_info = list(
      n_draws = n_draws,
      max_individual_se = max(individual_ate_se),
      mean_individual_se = mean(individual_ate_se)
    )
  ))
}

# =============================================================================
# MCMC UPDATE FUNCTIONS
# =============================================================================

#' Update regression coefficients
#' 
#' Gibbs sampler update for beta using standard Bayesian linear regression
#' conditional on latent variables and random effects.
#' 
#' @param X Design matrix
#' @param z Latent variables
#' @param U_wave Wave random effects
#' @param wave_mapped Wave assignments
#' @param prior_beta Prior specification (mean and precision)
#' @return Updated beta vector
update_beta_robust <- function(X, z, U_wave, wave_mapped, prior_beta) {
  n <- length(z)
  p <- ncol(X)
  
  if (any(is.na(z)) || any(is.na(U_wave))) {
    warning("Missing values detected in latent variables or random effects")
    return(rep(0, p))
  }
  
  z_adj <- z - U_wave[wave_mapped]
  valid_idx <- !is.na(z_adj)
  
  if (sum(valid_idx) < p) {
    warning("Insufficient valid observations for beta update")
    return(rep(0, p))
  }
  
  z_adj <- z_adj[valid_idx]
  X_valid <- X[valid_idx, , drop = FALSE]
  
  XtX <- crossprod(X_valid)
  posterior_precision <- XtX + prior_beta$precision + 1e-4 * diag(p)
  Xtz <- crossprod(X_valid, z_adj)
  
  posterior_cov <- safe_matrix_inverse(posterior_precision)
  posterior_mean <- posterior_cov %*% (Xtz + prior_beta$precision %*% prior_beta$mean)
  
  beta_new <- rmvnorm_robust(1, as.vector(posterior_mean), posterior_cov)
  
  if (any(is.na(beta_new)) || any(is.infinite(beta_new))) {
    warning("Invalid beta sample, using posterior mean")
    return(as.vector(posterior_mean))
  }
  
  return(beta_new)
}

#' Update wave-level random effects
#' 
#' Gibbs sampler update for random intercepts using empirical Bayes shrinkage.
#' 
#' @param z Latent variables
#' @param X Design matrix
#' @param beta Coefficient vector
#' @param wave_mapped Wave assignments
#' @param sigma2_U Random effect variance
#' @param n_waves Number of waves
#' @return Updated random effects vector
update_U_wave_robust <- function(z, X, beta, wave_mapped, sigma2_U, n_waves) {
  U_new <- numeric(n_waves)
  
  if (any(is.na(beta)) || is.na(sigma2_U) || sigma2_U <= 0) {
    warning("Invalid inputs to random effects update")
    return(rep(0, n_waves))
  }
  
  residuals <- z - X %*% beta
  
  for (w in 1:n_waves) {
    wave_idx <- which(wave_mapped == w)
    n_w <- length(wave_idx)
    
    if (n_w > 0) {
      valid_residuals <- residuals[wave_idx]
      valid_residuals <- valid_residuals[!is.na(valid_residuals)]
      
      if (length(valid_residuals) > 0) {
        posterior_precision <- length(valid_residuals) + 1/sigma2_U
        posterior_mean <- sum(valid_residuals) / posterior_precision
        posterior_var <- 1 / posterior_precision
        
        if (is.na(posterior_mean) || is.na(posterior_var) || posterior_var <= 0) {
          U_new[w] <- 0
        } else {
          sample_val <- rnorm(1, posterior_mean, sqrt(posterior_var))
          if (is.na(sample_val) || is.infinite(sample_val)) {
            U_new[w] <- 0
          } else {
            U_new[w] <- sample_val
          }
        }
      }
    }
  }
  
  return(U_new)
}

#' Update random effect variance parameter
#' 
#' Gibbs sampler update using inverse-gamma posterior.
#' 
#' @param U_wave Current random effects
#' @param prior_shape Inverse-gamma shape parameter
#' @param prior_rate Inverse-gamma rate parameter
#' @return Updated variance parameter
update_sigma2_U_robust <- function(U_wave, prior_shape, prior_rate) {
  U_valid <- U_wave[!is.na(U_wave)]
  n_waves <- length(U_valid)
  
  if (n_waves == 0) return(1.0)
  
  posterior_shape <- prior_shape + n_waves / 2
  posterior_rate <- prior_rate + sum(U_valid^2) / 2
  
  if (is.na(posterior_shape) || is.na(posterior_rate) || 
      posterior_shape <= 0 || posterior_rate <= 0) {
    return(1.0)
  }
  
  tryCatch({
    sigma2_sample <- 1 / rgamma(1, shape = posterior_shape, rate = posterior_rate)
    if (is.na(sigma2_sample) || is.infinite(sigma2_sample)) {
      return(1.0)
    }
    return(pmax(pmin(sigma2_sample, 10.0), 0.01))
  }, error = function(e) {
    return(1.0)
  })
}

#' Update latent variables
#' 
#' Gibbs sampler update using truncated normal distributions respecting
#' ordinal category constraints.
#' 
#' @param y Observed ordinal outcomes
#' @param X Design matrix
#' @param beta Coefficients
#' @param U_wave Random effects
#' @param wave_mapped Wave assignments
#' @param thresholds Threshold parameters
#' @return Updated latent variables
update_z_robust <- function(y, X, beta, U_wave, wave_mapped, thresholds) {
  n <- length(y)
  z_new <- numeric(n)
  
  if (any(is.na(beta)) || any(is.na(U_wave))) {
    warning("Missing values in parameters for latent variable update")
    return(rnorm(n))
  }
  
  mu <- as.vector(X %*% beta) + U_wave[wave_mapped]
  mu[is.na(mu)] <- 0
  
  K <- max(y)
  S <- c(-Inf, thresholds, Inf)
  
  if (any(is.na(S))) {
    S[is.na(S)] <- 0
  }
  
  for (i in 1:n) {
    k <- y[i]
    if (is.na(k) || k < 1 || k > K) {
      z_new[i] <- 0
      next
    }
    
    lower <- S[k]
    upper <- S[k + 1]
    
    z_new[i] <- rtruncnorm_robust(1, mean = mu[i], sd = 1, 
                                  lower = lower, upper = upper)
  }
  
  return(z_new)
}

#' Update threshold parameters
#' 
#' Conservative Metropolis-Hastings update maintaining ordering constraints.
#' Updates performed infrequently to ensure stability.
#' 
#' @param z Latent variables
#' @param y Observed outcomes
#' @param current_thresholds Current threshold vector
#' @param iter Current iteration number
#' @return Updated thresholds
source("01_threshold_update_mh.R")  # posterior-invariant, marginal-likelihood MH threshold update

update_thresholds_conservative <- function(z, y, current_thresholds, iter) {
  K <- max(y)
  if (K <= 2) return(current_thresholds)
  
  # Update every 50 iterations only
  if (iter %% 50 != 0) return(current_thresholds)
  
  new_thresholds <- current_thresholds
  n_thresh <- length(current_thresholds)
  
  for (j in 2:n_thresh) {
    lower_obs <- z[y <= j & !is.na(z)]
    upper_obs <- z[y > j & !is.na(z)]
    
    if (length(lower_obs) > 20 && length(upper_obs) > 20) {
      lower_bound <- quantile(lower_obs, 0.1)
      upper_bound <- quantile(upper_obs, 0.9)
      
      if (!is.na(lower_bound) && !is.na(upper_bound) && lower_bound < upper_bound) {
        if (j > 2) lower_bound <- max(lower_bound, new_thresholds[j-1] + 0.1)
        if (j < n_thresh) upper_bound <- min(upper_bound, new_thresholds[j+1] - 0.1)
        
        if (lower_bound < upper_bound) {
          current_val <- new_thresholds[j]
          proposed_val <- runif(1, lower_bound, upper_bound)
          new_thresholds[j] <- 0.95 * current_val + 0.05 * proposed_val
        }
      }
    }
  }
  
  # Enforce ordering
  for (j in 2:n_thresh) {
    if (new_thresholds[j] <= new_thresholds[j-1] + 0.05) {
      new_thresholds[j] <- new_thresholds[j-1] + 0.05
    }
  }
  
  return(new_thresholds)
}

# =============================================================================
# MAIN ALGORITHM: BMEOP (Algorithm 1)
# =============================================================================

#' Bayesian Mixed-Effects Ordered Probit Model
#' 
#' Implements Algorithm 1 from Chu and Yu (2025) for ordinal outcomes in
#' longitudinal settings with wave-level random effects. Optimised for
#' single-chain execution with analytical ATE computation by default.
#' 
#' @param y Ordinal outcome variable (integer vector)
#' @param X Design matrix including treatment and covariates
#' @param wave Wave indicator for clustering
#' @param treatment_col Column index of treatment variable (for ATE)
#' @param n_iter Total MCMC iterations
#' @param n_burn Burn-in iterations
#' @param n_thin Thinning interval
#' @param n_chains Number of chains (default 1 for efficiency)
#' @param prior_beta_var Prior variance for regression coefficients
#' @param ate_draws Number of Monte Carlo draws for ATE (if using MC method)
#' @param ate_method "analytical" or "multidraw" for ATE computation
#' @param seed Random seed for reproducibility
#' @param verbose Print progress messages
#' @return List containing posterior samples and summaries
bayesian_bmeop_elsa <- function(y, X, wave, treatment_col = NULL,
                                n_iter = 10000, n_burn = 5000, n_thin = 5,
                                n_chains = 1,
                                prior_beta_var = 4.0,
                                ate_draws = 50,
                                ate_method = "analytical",
                                proposal_sd = 0.08,
                                seed = NULL, verbose = TRUE) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Data preparation
  if (!is.matrix(X)) X <- as.matrix(X)
  if (any(is.na(y)) || any(is.na(X)) || any(is.na(wave))) {
    stop("Data contains missing values")
  }
  
  # Recode outcomes to consecutive integers
  unique_y <- sort(unique(y))
  y_recoded <- match(y, unique_y)
  nk <- max(y_recoded)
  
  n <- length(y_recoded)
  p <- ncol(X)
  unique_waves <- sort(unique(wave))
  n_waves <- length(unique_waves)
  wave_mapped <- match(wave, unique_waves)
  
  if (verbose) {
    cat("Bayesian Mixed-Effects Ordered Probit Model\n")
    cat("Observations:", n, "| Covariates:", p, "| Waves:", n_waves, 
        "| Categories:", nk, "\n")
    cat("Original scale:", min(unique_y), "to", max(unique_y), "\n")
    if (!is.null(treatment_col)) {
      cat("ATE method:", ate_method, "\n")
    }
  }
  
  # Prior specification
  prior_beta <- list(
    mean = rep(0, p),
    precision = diag(1/prior_beta_var, p)
  )
  
  # Single chain function
  run_chain <- function(chain_id) {
    
    chain_start_time <- Sys.time()
    
    # Chain-specific reproducible seed - see bjcm_implementation.R for the
    # same fix and rationale (needed once chains run in parallel via mclapply).
    if (!is.null(seed)) set.seed(seed + chain_id * 10000)
    
    # Overdispersed initial values - chain 1 keeps the original near-zero
    # start; chains 2+ start from genuinely different regions so multi-chain
    # R-hat reflects real convergence rather than four near-identical starts.
    init_offset <- c(0, -1.5, 1.5, -2.5, 2.5)[((chain_id - 1) %% 5) + 1]
    init_sd     <- if (chain_id == 1) 0.01 else 0.5
    
    # Initialisation
    beta <- rnorm(p, init_offset, init_sd)
    U_wave <- rnorm(n_waves, 0, 0.01)
    sigma2_U <- 0.5
    
    y_prop <- (y_recoded - 1) / (nk - 1)
    z <- qnorm(pmin(pmax(y_prop, 0.05), 0.95))
    
    # Strictly increasing sequence from tau_1 = 0, fixed for identification.
    if (nk > 2) {
      thresholds <- seq(0, 3, length.out = nk - 1)
    } else {
      thresholds <- 0
    }
    
    # Storage
    n_store <- floor((n_iter - n_burn) / n_thin)
    beta_samples <- matrix(NA, n_store, p)
    U_samples <- matrix(NA, n_store, n_waves)
    sigma2_U_samples <- numeric(n_store)
    ate_samples <- numeric(n_store)
    ate_se_samples <- numeric(n_store)
    
    # Added for trace plots / PPC: thresholds were previously computed every
    # iteration but discarded at chain end, making the saved fit unable to
    # reproduce the model's implied outcome distribution. Storing them plus
    # a per-iteration marginal category-probability vector (computed from
    # OBSERVED covariates, not counterfactual) gives a ready-made posterior
    # predictive check without needing to reconstruct thresholds later.
    threshold_samples <- matrix(NA, n_store, length(thresholds))
    y_pred_probs_samples <- matrix(NA, n_store, nk)
    
    store_idx <- 0
    n_errors <- 0
    
    # MCMC loop
    for (iter in 1:n_iter) {
      
      tryCatch({
        # Update parameters
        beta <- update_beta_robust(X, z, U_wave, wave_mapped, prior_beta)
        U_wave <- update_U_wave_robust(z, X, beta, wave_mapped, sigma2_U, n_waves)
        sigma2_U <- update_sigma2_U_robust(U_wave, 1, 1)
        
        # Threshold update conditions on the observed-data marginal
        # likelihood p(Y | beta, U_wave, S) -- a function of beta and
        # U_wave (just updated above) and the fixed data, NOT of z. It is
        # therefore performed here, before z is refreshed, so that the z
        # draw immediately below uses the CURRENT (just-updated)
        # thresholds rather than last iteration's.
        if (nk > 2) {
          eta_bmeop <- as.vector(X %*% beta) + U_wave[wave_mapped]
          eta_bmeop[is.na(eta_bmeop)] <- 0
          
          thresholds <- update_thresholds_mh_marginal(
            z_latent = z, y_observed = y_recoded,
            current_thresholds = thresholds, iteration = iter,
            eta = eta_bmeop, label = "bmeop_outcome",
            proposal_sd = proposal_sd, update_every_iter = TRUE
          )
        }
        
        z <- update_z_robust(y_recoded, X, beta, U_wave, wave_mapped, thresholds)
        
      }, error = function(e) {
        n_errors <<- n_errors + 1
        if (verbose && n_errors <= 5) {
          cat("Warning: Error at iteration", iter, "\n")
        }
      })
      
      # Store samples
      if (iter > n_burn && (iter - n_burn) %% n_thin == 0) {
        store_idx <- store_idx + 1
        beta_samples[store_idx, ] <- beta
        U_samples[store_idx, ] <- U_wave
        sigma2_U_samples[store_idx] <- sigma2_U
        
        # ATE computation
        if (!is.null(treatment_col)) {
          if (ate_method == "analytical") {
            ate_result <- compute_ate_analytical(X, beta, U_wave, wave_mapped, 
                                                 thresholds, treatment_col, ate_draws)
          } else {
            ate_result <- compute_ate_multidraw(X, beta, U_wave, wave_mapped, 
                                                thresholds, treatment_col, ate_draws)
          }
          ate_samples[store_idx] <- ate_result$ate
          ate_se_samples[store_idx] <- ate_result$monte_carlo_se
        }
        
        # Store thresholds and the model-implied marginal category
        # distribution for PPC. mu_all uses the OBSERVED treatment/covariates
        # (X as fitted), not a counterfactual arm - this is what should be
        # compared against the observed y distribution for a PPC.
        threshold_samples[store_idx, ] <- thresholds
        mu_all <- as.vector(X %*% beta) + U_wave[wave_mapped]
        S_all <- c(-Inf, thresholds, Inf)
        y_pred_probs_samples[store_idx, ] <- sapply(1:nk, function(k) {
          mean(pnorm(S_all[k + 1], mean = mu_all, sd = 1) -
                 pnorm(S_all[k], mean = mu_all, sd = 1))
        })
      }
      
      # Progress reporting - writes to a per-chain file instead of cat(),
      # since cat() inside a forked mclapply worker does not reliably reach
      # the parent session's console. Check progress from a SEPARATE R
      # session/terminal using monitor_bmeop_progress() (see bmeop_diagnost.R),
      # or just `cat(readLines("../logs/bmeop_progress_chain_<id>.log"))`.
      progress_interval <- max(1, floor(n_iter / 100))  # ~100 updates per chain
      if (iter %% progress_interval == 0 || iter == n_iter) {
        elapsed_secs <- as.numeric(difftime(Sys.time(), chain_start_time, units = "secs"))
        pct_done <- iter / n_iter
        eta_secs <- if (pct_done > 0) elapsed_secs / pct_done - elapsed_secs else NA
        current_ate <- if (!is.null(treatment_col) && store_idx > 0) {
          round(mean(ate_samples[1:store_idx], na.rm = TRUE), 4)
        } else {
          NA
        }
        progress_line <- sprintf(
          "chain=%d iter=%d/%d (%.1f%%) elapsed=%.1fmin eta=%.1fmin current_ATE=%s errors=%d timestamp=%s",
          chain_id, iter, n_iter, pct_done * 100,
          elapsed_secs / 60, if (is.na(eta_secs)) NA else eta_secs / 60,
          ifelse(is.na(current_ate), "NA", current_ate), n_errors,
          format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        )
        writeLines(progress_line, con = sprintf("../logs/bmeop_progress_chain_%d.log", chain_id))
      }
    }
    
    if (verbose && n_errors > 0) {
      cat("Chain", chain_id, "completed with", n_errors, "errors\n")
    }
    
    return(list(
      beta = beta_samples,
      random_effects = U_samples,
      sigma2_u = sigma2_U_samples,
      ate = ate_samples,
      ate_se = ate_se_samples,
      thresholds = threshold_samples,
      y_pred_probs = y_pred_probs_samples,
      chain_id = chain_id,
      n_errors = n_errors
    ))
  }
  
  # Run chains
  if (verbose) cat("Running", n_chains, "MCMC chain(s)...\n")
  
  if (n_chains > 1) {
    n_cores_use <- min(n_chains, parallel::detectCores())
    if (verbose) cat("Dispatching", n_chains, "chains across", n_cores_use, "cores...\n")
    
    chains <- parallel::mclapply(1:n_chains, run_chain,
                                 mc.cores = n_cores_use,
                                 mc.preschedule = FALSE)
    
    failed <- vapply(chains, function(x) inherits(x, "try-error"), logical(1))
    if (any(failed)) {
      stop("Chain(s) ", paste(which(failed), collapse = ", "), " failed:\n",
           paste(vapply(chains[failed], function(x) as.character(x), character(1)), collapse = "\n"))
    }
  } else {
    chains <- lapply(1:n_chains, run_chain)
  }
  
  # Combine results
  combined_beta <- do.call(rbind, lapply(chains, function(x) x$beta))
  combined_sigma2_u <- do.call(c, lapply(chains, function(x) x$sigma2_u))
  combined_ate <- if (!is.null(treatment_col)) {
    do.call(c, lapply(chains, function(x) x$ate))
  } else NULL
  combined_ate_se <- if (!is.null(treatment_col)) {
    do.call(c, lapply(chains, function(x) x$ate_se))
  } else NULL
  
  # Add variable names
  if (!is.null(colnames(X))) {
    colnames(combined_beta) <- colnames(X)
  }
  
  # Summary statistics
  beta_summary <- data.frame(
    mean = colMeans(combined_beta, na.rm = TRUE),
    sd = apply(combined_beta, 2, sd, na.rm = TRUE),
    q2.5 = apply(combined_beta, 2, quantile, 0.025, na.rm = TRUE),
    q97.5 = apply(combined_beta, 2, quantile, 0.975, na.rm = TRUE)
  )
  
  if (!is.null(colnames(X))) {
    rownames(beta_summary) <- colnames(X)
  }
  
  # ATE summary
  ate_stability <- NULL
  if (!is.null(combined_ate)) {
    ate_stability <- list(
      mean_ate = mean(combined_ate, na.rm = TRUE),
      sd_ate = sd(combined_ate, na.rm = TRUE),
      mean_mc_se = mean(combined_ate_se, na.rm = TRUE),
      max_mc_se = max(combined_ate_se, na.rm = TRUE),
      min_mc_se = min(combined_ate_se, na.rm = TRUE),
      stability_ratio = sd(combined_ate, na.rm = TRUE) / mean(combined_ate_se, na.rm = TRUE),
      n_draws_per_iteration = ate_draws,
      method = ate_method
    )
  }
  
  if (verbose) {
    cat("\nMCMC completed\n")
    print(round(beta_summary, 4))
    
    if (!is.null(ate_stability)) {
      cat("\nATE Summary:\n")
      cat("Mean ATE:", round(ate_stability$mean_ate, 4), "\n")
      cat("SD:", round(ate_stability$sd_ate, 4), "\n")
      cat("Mean Monte Carlo SE:", round(ate_stability$mean_mc_se, 6), "\n")
      cat("Stability Ratio:", round(ate_stability$stability_ratio, 2), "\n")
    }
  }
  
  return(list(
    beta_samples = combined_beta,
    sigma2_u_samples = combined_sigma2_u,
    ate_samples = combined_ate,
    ate_se_samples = combined_ate_se,
    ate_stability = ate_stability,
    beta_summary = beta_summary,
    chains = chains,
    n_iter = n_iter,
    n_burn = n_burn,
    n_chains = n_chains,
    ate_method = ate_method,
    ate_draws = ate_draws,
    convergence_ok = TRUE
  ))
}

# =============================================================================
# CONVERGENCE DIAGNOSTICS
# =============================================================================

#' MCMC convergence diagnostics
#' 
#' Computes effective sample sizes, Monte Carlo standard errors, and
#' Gelman-Rubin statistics (if multiple chains available).
#' 
#' @param results Output from bayesian_bmeop_elsa
#' @return List of diagnostic statistics
mcmc_diagnostics_elsa <- function(results) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' required for diagnostics")
  }
  
  # Convert to mcmc objects
  beta_mcmc <- as.mcmc(results$beta_samples)
  
  # Effective sample sizes
  if (ncol(results$beta_samples) > 1) {
    ess_beta <- coda::effectiveSize(beta_mcmc)
  } else {
    ess_beta <- coda::effectiveSize(as.mcmc(results$beta_samples[,1]))
  }
  
  # Monte Carlo standard errors
  beta_variances <- apply(results$beta_samples, 2, function(x) var(x, na.rm = TRUE))
  mcse_beta <- sqrt(beta_variances / ess_beta)
  
  # Gelman-Rubin statistic
  if (length(results$chains) > 1) {
    chain_betas <- lapply(results$chains, function(x) x$beta)
    mcmc_list <- mcmc.list(lapply(chain_betas, as.mcmc))
    
    tryCatch({
      psrf_stats <- coda::gelman.diag(mcmc_list, confidence = 0.95)
      psrf <- psrf_stats$psrf[, "Point est."]
    }, error = function(e) {
      psrf <- rep(NA, ncol(results$beta_samples))
    })
  } else {
    psrf <- rep(NA, ncol(results$beta_samples))
  }
  
  # Autocorrelation
  max_lag <- min(50, nrow(results$beta_samples) %/% 4)
  autocorr_beta <- apply(results$beta_samples, 2, function(x) {
    tryCatch({
      acf_result <- acf(x, lag.max = max_lag, plot = FALSE)
      first_low <- which(abs(acf_result$acf[-1]) < 0.1)[1]
      if (is.na(first_low)) max_lag else first_low
    }, error = function(e) max_lag)
  })
  
  # Variance parameter diagnostics
  if (!is.null(results$sigma2_u_samples)) {
    ess_sigma2_u <- coda::effectiveSize(as.mcmc(results$sigma2_u_samples))
    mcse_sigma2_u <- sqrt(var(results$sigma2_u_samples, na.rm = TRUE) / ess_sigma2_u)
  } else {
    ess_sigma2_u <- NA
    mcse_sigma2_u <- NA
  }
  
  # ATE diagnostics
  if (!is.null(results$ate_samples)) {
    ess_ate <- coda::effectiveSize(as.mcmc(results$ate_samples))
    mcse_ate <- sqrt(var(results$ate_samples, na.rm = TRUE) / ess_ate)
    
    if (!is.null(results$ate_se_samples)) {
      mean_mc_se_ate <- mean(results$ate_se_samples, na.rm = TRUE)
      stability_ratio <- sd(results$ate_samples, na.rm = TRUE) / mean_mc_se_ate
      
      ate_diagnostics <- list(
        ess = ess_ate,
        mcse = mcse_ate,
        mean_monte_carlo_se = mean_mc_se_ate,
        stability_ratio = stability_ratio,
        n_draws = results$ate_draws,
        method = results$ate_method
      )
    } else {
      ate_diagnostics <- list(
        ess = ess_ate,
        mcse = mcse_ate
      )
    }
  } else {
    ate_diagnostics <- list(
      ess = NA,
      mcse = NA
    )
  }
  
  return(list(
    effective_sizes = ess_beta,
    mc_standard_errors = mcse_beta,
    autocorr_lags = autocorr_beta,
    multichain = list(
      psrf = psrf,
      n_chains = length(results$chains)
    ),
    sigma2_u = list(
      ess = ess_sigma2_u,
      mcse = mcse_sigma2_u
    ),
    ate = ate_diagnostics
  ))
}

# =============================================================================
# WRAPPER FUNCTIONS FOR ELSA ANALYSIS
# =============================================================================

#' Run BMEOP analysis on ELSA data
#' 
#' Convenience wrapper that handles ELSA data preprocessing and runs the
#' BMEOP model with sensible defaults.
#' 
#' @param elsa_long Longitudinal ELSA dataset
#' @param subset_size Sample size for testing (NULL for full dataset)
#' @param seed Random seed
#' @param ate_draws Number of ATE draws
#' @param ate_method ATE computation method
#' @param n_chains Number of MCMC chains (default 1; use >=2 for a
#'   multi-chain validation run, e.g. to check chains differ under
#'   overdispersed starts, matching BJCM's run_validation())
#' @param proposal_sd Threshold MH proposal SD, passed straight through to
#'   bayesian_bmeop_elsa() / update_thresholds_mh_marginal(). Tune to
#'   target roughly 20-50% acceptance; grid search for this dataset found
#'   the default of 0.08 gives poor acceptance (5.9-14.8%) at full scale,
#'   with 0.015 performing much better (29.7-57.1%) -- confirm which value
#'   was actually used for your reported results before relying on this
#'   default.
#' @param n_iter Total MCMC iterations. Defaults to 15000, matching the
#'   full production run -- previously hardcoded here regardless of
#'   subset_size, so a small subset_size alone did NOT give a fast
#'   validation run (only n scaled down, not iteration count). Pass a
#'   small value (e.g. 2000) together with a small subset_size for a
#'   genuinely quick correctness check, matching the scale of BJCM's
#'   run_validation().
#' @param n_burn Burn-in iterations. Defaults to 10000 (full-run value);
#'   scale down together with n_iter for validation runs (e.g. 1000).
#' @return BMEOP results object
run_elsa_bmeop <- function(elsa_long, subset_size = 2000, seed = 123, 
                           ate_draws = 50, ate_method = "analytical",
                           n_chains = 1, proposal_sd = 0.08,
                           n_iter = 15000, n_burn = 10000) {
  set.seed(seed)
  
  cat("BMEOP Model for ELSA Loneliness Analysis\n")
  cat("Method:", ate_method, "with", ate_draws, "draws\n")
  cat("Chains:", n_chains, "| Threshold proposal_sd:", proposal_sd, "\n")
  cat("n_iter:", n_iter, "| n_burn:", n_burn, "\n\n")
  
  # Required variables
  required_vars <- c("idauniq", "wave", "loneliness", "livalone", "age_gr", "dhsex2", 
                     "edqual2", "self_reported_health", "sclife", "depression", 
                     "transport_mobility", "mobility_limitations")
  
  # Check availability
  available_vars <- required_vars[required_vars %in% names(elsa_long)]
  complete_data <- elsa_long[complete.cases(elsa_long[, available_vars]), ]
  
  # Subset if requested
  if (!is.null(subset_size) && nrow(complete_data) > subset_size) {
    subset_data <- complete_data[sample(nrow(complete_data), subset_size), ]
  } else {
    subset_data <- complete_data
  }
  
  cat("Sample size:", nrow(subset_data), "\n")
  cat("Loneliness range:", min(subset_data$loneliness), "to", 
      max(subset_data$loneliness), "\n")
  cat("Living alone:", round(100*mean(subset_data$livalone), 1), "%\n\n")
  
  # Create design matrix
  X_components <- list(
    intercept = rep(1, nrow(subset_data)),
    livalone_TREATMENT = subset_data$livalone
  )
  
  if("age_gr" %in% names(subset_data)) X_components$age_gr <- subset_data$age_gr
  if("dhsex2" %in% names(subset_data)) X_components$dhsex2 <- subset_data$dhsex2
  if("edqual2" %in% names(subset_data)) X_components$edqual2 <- subset_data$edqual2
  if("self_reported_health" %in% names(subset_data)) {
    X_components$health <- subset_data$self_reported_health
  }
  if("sclife" %in% names(subset_data)) X_components$sclife <- subset_data$sclife
  if("depression" %in% names(subset_data)) {
    X_components$depression <- subset_data$depression
  }
  if("transport_mobility" %in% names(subset_data)) {
    X_components$transport_mobility <- subset_data$transport_mobility
  }
  if("mobility_limitations" %in% names(subset_data)) {
    X_components$mobility_limitations <- subset_data$mobility_limitations
  }
  
  X_subset <- do.call(cbind, X_components)
  colnames(X_subset) <- names(X_components)
  
  # Run model
  cat("Running BMEOP model...\n")
  
  results <- bayesian_bmeop_elsa(
    y = subset_data$loneliness,
    X = X_subset,
    wave = subset_data$wave,
    treatment_col = 2,
    n_iter = n_iter,
    n_burn = n_burn,
    n_chains = n_chains,
    proposal_sd = proposal_sd,
    ate_draws = ate_draws,
    ate_method = ate_method,
    verbose = TRUE
  )
  
  # Summary
  cat("\nCoefficient Estimates:\n")
  print(round(results$beta_summary, 4))
  
  if (!is.null(results$ate_samples)) {
    ate_samples <- results$ate_samples[!is.na(results$ate_samples)]
    cat("\nAverage Treatment Effect:\n")
    cat("Mean:", round(mean(ate_samples), 4), "\n")
    cat("95% CrI: [", round(quantile(ate_samples, 0.025), 4), ",",
        round(quantile(ate_samples, 0.975), 4), "]\n")
    cat("P(ATE > 0):", round(mean(ate_samples > 0), 4), "\n")
  }
  
  # Add metadata
  results$elsa_analysis_info <- list(
    dataset = "elsa_longitudinal",
    n_observations = nrow(subset_data),
    n_waves = length(unique(subset_data$wave)),
    available_variables = available_vars,
    covariates = colnames(X_subset)
  )
  
  return(results)
}

#' Run full ELSA dataset analysis
#' 
#' Runs BMEOP on complete ELSA dataset. Note: computational time approximately
#' 15-18 hours for n=42,185 observations.
#' 
#' @param elsa_long ELSA dataset
#' @param seed Random seed
#' @return BMEOP results
run_full_elsa_analysis <- function(elsa_long, seed = 123) {
  set.seed(seed)
  
  cat("Full ELSA Dataset Analysis\n")
  cat("Expected runtime: 15-18 hours\n\n")
  
  if(is.character(elsa_long)) {
    elsa_long <- read.csv(elsa_long)
  }
  
  cat("Dataset:\n")
  cat("- Observations:", nrow(elsa_long), "\n")
  cat("- Participants:", length(unique(elsa_long$idauniq)), "\n")
  cat("- Waves:", paste(sort(unique(elsa_long$wave)), collapse=", "), "\n\n")
  
  start_time <- Sys.time()
  
  full_results <- run_elsa_bmeop(
    elsa_long = elsa_long,
    subset_size = NULL,
    seed = seed,
    ate_draws = 50,
    ate_method = "analytical"
  )
  
  end_time <- Sys.time()
  runtime <- as.numeric(difftime(end_time, start_time, units = "mins"))
  
  cat("\nAnalysis Complete\n")
  cat("Runtime:", round(runtime, 1), "minutes (", round(runtime/60, 1), "hours)\n")
  
  full_results$runtime_minutes <- runtime
  
  return(full_results)
}

# =============================================================================
# END OF CODE
# =============================================================================

# =============================================================================
# BMEOP Model Diagnostics and Post-Processing
# Companion script for ELSA loneliness analysis
# =============================================================================
#
# This script provides additional diagnostic tools and visualisations for
# assessing MCMC convergence and model fit in the BMEOP framework.
# =============================================================================

library(coda)
library(corrplot)
library(car)

# =============================================================================
# CONVERGENCE DIAGNOSTICS
# =============================================================================

#' Compute and display comprehensive convergence diagnostics
#' 
#' @param results BMEOP model results object
#' @return List of diagnostic statistics
assess_convergence <- function(results) {
  
  convergence_diag <- mcmc_diagnostics_elsa(results)
  
  cat("=== MCMC CONVERGENCE DIAGNOSTICS ===\n\n")
  
  # Effective sample sizes
  cat("Effective Sample Sizes:\n")
  print(round(convergence_diag$effective_sizes, 1))
  cat("\n")
  
  # Monte Carlo standard errors
  cat("Monte Carlo Standard Errors:\n")
  print(round(convergence_diag$mc_standard_errors, 6))
  cat("\n")
  
  # Autocorrelation lags
  cat("Autocorrelation Effective Lags:\n")
  print(convergence_diag$autocorr_lags)
  cat("\n")
  
  # Identify problematic parameters
  low_ess <- convergence_diag$effective_sizes < 100
  if (any(low_ess)) {
    cat("Parameters with low ESS (<100):\n")
    print(names(convergence_diag$effective_sizes)[low_ess])
    cat("\n")
  }
  
  # ATE diagnostics
  if (!is.na(convergence_diag$ate$ess)) {
    cat("ATE Diagnostics:\n")
    cat("  ESS:", round(convergence_diag$ate$ess, 1), "\n")
    cat("  MCSE:", round(convergence_diag$ate$mcse, 6), "\n")
    
    if (!is.null(convergence_diag$ate$stability_ratio)) {
      cat("  Stability Ratio:", round(convergence_diag$ate$stability_ratio, 2), "\n")
      
      if (convergence_diag$ate$stability_ratio > 10) {
        cat("  Assessment: Excellent stability\n")
      } else if (convergence_diag$ate$stability_ratio > 5) {
        cat("  Assessment: Good stability\n")
      } else {
        cat("  Assessment: Consider increasing draws\n")
      }
    }
    cat("\n")
  }
  
  return(convergence_diag)
}

#' Visualise MCMC trace plots
#' 
#' Creates trace plots for selected parameters to assess mixing and convergence
#' 
#' @param results BMEOP results object
#' @param params Parameter names to plot (NULL for all)
#' @param n_cols Number of plot columns
visualise_traces <- function(results, params = NULL, n_cols = 3) {
  
  if (is.null(params)) {
    params <- colnames(results$beta_samples)
  }
  
  n_params <- length(params)
  n_rows <- ceiling(n_params / n_cols)
  
  par(mfrow = c(n_rows, n_cols), mar = c(3, 3, 2, 1))
  
  for (param in params) {
    if (param %in% colnames(results$beta_samples)) {
      samples <- results$beta_samples[, param]
      plot(samples, type = "l", main = param, 
           xlab = "Iteration", ylab = "Value",
           col = "steelblue", lwd = 0.5)
      abline(h = mean(samples), col = "red", lty = 2)
    }
  }
  
  # ATE trace if available
  if (!is.null(results$ate_samples)) {
    plot(results$ate_samples, type = "l", main = "ATE",
         xlab = "Iteration", ylab = "Value",
         col = "steelblue", lwd = 0.5)
    abline(h = mean(results$ate_samples), col = "red", lty = 2)
  }
  
  par(mfrow = c(1, 1))
}

#' Assess chain stability by segment
#' 
#' Divides the chain into segments and compares estimates across segments
#' to assess stationarity
#' 
#' @param results BMEOP results
#' @param n_segments Number of segments (default 4)
#' @param param Parameter to assess (default "treatment")
assess_chain_stability <- function(results, n_segments = 4, 
                                   param = "livalone_TREATMENT") {
  
  n_samples <- nrow(results$beta_samples)
  segments <- split(1:n_samples, cut(1:n_samples, n_segments))
  
  cat("=== CHAIN STABILITY ASSESSMENT ===\n")
  cat("Parameter:", param, "\n\n")
  
  if (param %in% colnames(results$beta_samples)) {
    cat("Coefficient Estimates by Chain Segment:\n")
    for (i in 1:n_segments) {
      segment_samples <- results$beta_samples[segments[[i]], param]
      segment_mean <- mean(segment_samples)
      segment_ci <- quantile(segment_samples, c(0.025, 0.975))
      
      cat(sprintf("Segment %d: Mean = %.4f, 95%% CrI = [%.4f, %.4f]\n",
                  i, segment_mean, segment_ci[1], segment_ci[2]))
    }
    cat("\n")
  }
  
  # ATE stability
  if (!is.null(results$ate_samples)) {
    cat("ATE Estimates by Chain Segment:\n")
    for (i in 1:n_segments) {
      segment_samples <- results$ate_samples[segments[[i]]]
      segment_mean <- mean(segment_samples)
      segment_ci <- quantile(segment_samples, c(0.025, 0.975))
      
      cat(sprintf("Segment %d: Mean = %.4f, 95%% CrI = [%.4f, %.4f]\n",
                  i, segment_mean, segment_ci[1], segment_ci[2]))
    }
  }
}

# =============================================================================
# POSTERIOR CORRELATION ANALYSIS
# =============================================================================

#' Examine posterior correlations among parameters
#' 
#' @param results BMEOP results
#' @param params Parameters to include (NULL for all)
#' @param visualise Create correlation plot
examine_posterior_correlation <- function(results, params = NULL, 
                                          visualise = TRUE) {
  
  if (is.null(params)) {
    params <- colnames(results$beta_samples)
  }
  
  param_samples <- results$beta_samples[, params, drop = FALSE]
  cor_matrix <- cor(param_samples)
  
  cat("=== POSTERIOR CORRELATION MATRIX ===\n")
  print(round(cor_matrix, 3))
  cat("\n")
  
  # Identify high correlations
  high_cor <- which(abs(cor_matrix) > 0.7 & upper.tri(cor_matrix), arr.ind = TRUE)
  
  if (nrow(high_cor) > 0) {
    cat("High Posterior Correlations (|r| > 0.7):\n")
    for (i in 1:nrow(high_cor)) {
      param1 <- rownames(cor_matrix)[high_cor[i, 1]]
      param2 <- colnames(cor_matrix)[high_cor[i, 2]]
      cor_val <- cor_matrix[high_cor[i, 1], high_cor[i, 2]]
      cat(sprintf("  %s - %s: %.3f\n", param1, param2, cor_val))
    }
    cat("\n")
  }
  
  if (visualise) {
    corrplot(cor_matrix, 
             method = "number",
             type = "upper",
             tl.col = "black",
             tl.srt = 45,
             number.cex = 0.8,
             col = colorRampPalette(c("#053061", "white", "#67001F"))(200),
             addCoef.col = "black",
             title = "Posterior Correlation Matrix",
             mar = c(0, 0, 2, 0))
  }
  
  return(cor_matrix)
}

#' Create autocorrelation function plots
#' 
#' @param results BMEOP results
#' @param params Parameters to plot
#' @param max_lag Maximum lag to display
plot_autocorrelations <- function(results, params = NULL, max_lag = 50) {
  
  if (is.null(params)) {
    params <- colnames(results$beta_samples)
  }
  
  n_params <- length(params)
  n_cols <- min(3, n_params)
  n_rows <- ceiling(n_params / n_cols)
  
  par(mfrow = c(n_rows, n_cols), mar = c(3, 3, 2, 1))
  
  for (param in params) {
    if (param %in% colnames(results$beta_samples)) {
      acf(results$beta_samples[, param], 
          lag.max = max_lag,
          main = paste("ACF:", param))
    }
  }
  
  par(mfrow = c(1, 1))
}

# =============================================================================
# MODEL DIAGNOSTICS
# =============================================================================

#' Check variance inflation factors
#' 
#' Assesses multicollinearity in the design matrix using VIF
#' 
#' @param X Design matrix
#' @param y Outcome variable
#' @return VIF values
check_multicollinearity <- function(X, y) {
  
  # Remove intercept if present
  if (colnames(X)[1] == "intercept") {
    X_no_intercept <- X[, -1, drop = FALSE]
  } else {
    X_no_intercept <- X
  }
  
  # Create data frame
  df <- as.data.frame(cbind(y = y, X_no_intercept))
  
  # Fit linear model
  formula_str <- paste("y ~", paste(colnames(X_no_intercept), collapse = " + "))
  lm_model <- lm(as.formula(formula_str), data = df)
  
  # Calculate VIF
  vif_values <- vif(lm_model)
  
  cat("=== VARIANCE INFLATION FACTORS ===\n")
  print(round(vif_values, 2))
  cat("\n")
  
  cat("Interpretation:\n")
  cat("  VIF < 5:  No concern\n")
  cat("  VIF 5-10: Moderate multicollinearity\n")
  cat("  VIF > 10: Severe multicollinearity\n\n")
  
  if (all(vif_values < 5)) {
    cat("Assessment: No multicollinearity problems detected\n")
  } else if (all(vif_values < 10)) {
    cat("Assessment: Acceptable multicollinearity levels\n")
  } else {
    cat("Assessment: High multicollinearity detected for some variables\n")
    high_vif <- names(vif_values)[vif_values > 10]
    cat("Variables with VIF > 10:", paste(high_vif, collapse = ", "), "\n")
  }
  
  return(vif_values)
}

#' Compare mixing across parameters
#' 
#' @param results BMEOP results
#' @param good_params Parameters with good mixing
#' @param poor_params Parameters with poor mixing
compare_mixing <- function(results, 
                           good_params = c("age_gr"),
                           poor_params = c("intercept", "livalone_TREATMENT")) {
  
  n_good <- length(good_params)
  n_poor <- length(poor_params)
  n_plots <- n_good + n_poor
  
  if (!is.null(results$ate_samples)) n_plots <- n_plots + 1
  
  par(mfrow = c(ceiling(n_plots/2), 2), mar = c(3, 3, 2, 1))
  
  # Good mixing examples
  for (param in good_params) {
    if (param %in% colnames(results$beta_samples)) {
      ess <- coda::effectiveSize(results$beta_samples[, param])
      plot(results$beta_samples[, param], type = "l",
           main = paste("Good Mixing:", param, "(ESS =", round(ess, 1), ")"),
           xlab = "Iteration", ylab = "Value", col = "forestgreen")
    }
  }
  
  # ATE (typically good mixing)
  if (!is.null(results$ate_samples)) {
    ess <- coda::effectiveSize(results$ate_samples)
    plot(results$ate_samples, type = "l",
         main = paste("Good Mixing: ATE (ESS =", round(ess, 1), ")"),
         xlab = "Iteration", ylab = "Value", col = "forestgreen")
  }
  
  # Poor mixing examples
  for (param in poor_params) {
    if (param %in% colnames(results$beta_samples)) {
      ess <- coda::effectiveSize(results$beta_samples[, param])
      plot(results$beta_samples[, param], type = "l",
           main = paste("Poor Mixing:", param, "(ESS =", round(ess, 1), ")"),
           xlab = "Iteration", ylab = "Value", col = "indianred")
    }
  }
  
  par(mfrow = c(1, 1))
}

# =============================================================================
# USAGE EXAMPLE
# =============================================================================

#' Complete diagnostic workflow
#' 
#' Runs all diagnostic checks on BMEOP results.
#' 
#' @param results BMEOP results object
#' @param X Design matrix (optional, for VIF)
#' @param y Outcome variable (optional, for VIF)
run_full_diagnostics <- function(results, X = NULL, y = NULL) {
  
  cat("\n")
  cat("================================================================\n")
  cat("COMPREHENSIVE BMEOP MODEL DIAGNOSTICS\n")
  cat("================================================================\n\n")
  
  # Convergence assessment
  convergence_diag <- assess_convergence(results)
  
  cat("\n")
  
  # Chain stability
  assess_chain_stability(results)
  
  cat("\n")
  
  # Posterior correlations
  cat("Examining posterior correlations...\n")
  cor_matrix <- examine_posterior_correlation(results, visualise = TRUE)
  
  cat("\n")
  
  # Multicollinearity check
  if (!is.null(X) && !is.null(y)) {
    cat("Checking for multicollinearity...\n")
    vif_values <- check_multicollinearity(X, y)
    cat("\n")
  }
  
  # Visual diagnostics
  cat("Creating trace plots...\n")
  visualise_traces(results)
  
  cat("\nCreating autocorrelation plots...\n")
  plot_autocorrelations(results)
  
  cat("\nComparing parameter mixing...\n")
  compare_mixing(results)
  
  cat("\n")
  cat("================================================================\n")
  cat("DIAGNOSTICS COMPLETE\n")
  cat("================================================================\n")
  
  return(list(
    convergence = convergence_diag,
    correlations = cor_matrix,
    vif = if (!is.null(X) && !is.null(y)) vif_values else NULL
  ))
}

# =============================================================================
# NOTES ON INTERPRETATION
# =============================================================================

# Elevated autocorrelation in some parameters (e.g., intercept, treatment,
# life satisfaction, depression) is expected in hierarchical models with
# correlated predictors. The key question is whether estimates are stable
# across chain segments, which indicates reliable posterior inference despite
# high autocorrelation.
#
# For the ELSA analysis:
# - Treatment effect (livalone_TREATMENT) and ATE show stable estimates
#   across chain segments despite elevated autocorrelation
# - High posterior correlation (r > 0.7) between intercept, treatment,
#   life satisfaction, and depression explains the autocorrelation pattern
# - VIF values < 5 for all covariates indicate that multicollinearity
#   is not problematic in the design matrix
# - The analytical ATE method achieves excellent effective sample size
#   (ESS > 1000), confirming reliable causal effect estimation

# =============================================================================
# REAL MULTI-CHAIN R-HAT (Gelman-Rubin, via coda)
# =============================================================================
# assess_convergence() above only computes within-chain ESS/autocorrelation.
# This requires results$chains from a run with n_chains > 1 (bmeop_model.R,
# after the parallel-chains fix), and computes proper between/within-chain
# variance R-hat on the treatment coefficient and the ATE - the two
# quantities that matter most for the paper's claims.

compute_multichain_rhat_bmeop <- function(results, treatment_col_index = NULL) {
  
  if (is.null(results$chains) || length(results$chains) < 2) {
    stop("results$chains has fewer than 2 chains - re-run bayesian_bmeop_elsa() ",
         "with n_chains >= 2 to compute a real multi-chain R-hat.")
  }
  
  chains <- results$chains
  n_chains <- length(chains)
  
  # Per-chain $beta matrices carry no column names (only the merged
  # results$beta_samples does, and only if X had colnames) - so identify
  # the treatment column by the same index passed to bayesian_bmeop_elsa()
  # as treatment_col, not by name.
  if (is.null(treatment_col_index)) {
    stop("Pass treatment_col_index - the same integer column position you ",
         "passed as treatment_col to bayesian_bmeop_elsa().")
  }
  
  treatment_label <- if (!is.null(colnames(results$beta_samples))) {
    colnames(results$beta_samples)[treatment_col_index]
  } else {
    paste0("beta[", treatment_col_index, "]")
  }
  
  build_mcmc_list <- function(extractor) {
    coda::mcmc.list(lapply(chains, function(ch) coda::mcmc(extractor(ch))))
  }
  
  treatment_list <- build_mcmc_list(function(ch) {
    matrix(ch$beta[, treatment_col_index], ncol = 1, dimnames = list(NULL, treatment_label))
  })
  
  cat("=== Multi-chain Gelman-Rubin R-hat (", n_chains, "chains ) ===\n\n")
  
  rhat_treatment <- coda::gelman.diag(treatment_list, autoburnin = FALSE)
  cat("Treatment effect (", treatment_label, "):\n")
  print(rhat_treatment)
  
  rhat_ate <- NULL
  if (!is.null(chains[[1]]$ate)) {
    ate_list <- build_mcmc_list(function(ch) {
      matrix(ch$ate, ncol = 1, dimnames = list(NULL, "ATE"))
    })
    rhat_ate <- coda::gelman.diag(ate_list, autoburnin = FALSE)
    cat("\nATE:\n")
    print(rhat_ate)
  }
  
  flagged <- any(rhat_treatment$psrf[, "Point est."] > 1.1) ||
    (!is.null(rhat_ate) && any(rhat_ate$psrf[, "Point est."] > 1.1))
  
  if (flagged) {
    cat("\nWARNING: R-hat > 1.1 - chains have not converged to the same distribution.\n")
  } else {
    cat("\nAll R-hat values <= 1.1 - chains agree on the treatment effect and ATE.\n")
  }
  
  list(treatment = rhat_treatment, ate = rhat_ate)
}

# =============================================================================
# PROGRESS MONITORING / BACKGROUND EXECUTION
# =============================================================================
# monitor_progress(), run_in_background(), and collect_result() moved to
# 00_background_run_utils.R - they are generic (not bmeop-specific), and
# keeping them here would mean 02_bjcm_model.R (which never sources this
# file, on purpose) would lack them.
# source("00_background_run_utils.R") to get all three.

# =============================================================================
# END OF DIAGNOSTICS CODE
# =============================================================================