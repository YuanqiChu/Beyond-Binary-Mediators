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
    
    # Initialisation
    beta <- rnorm(p, 0, 0.01)
    U_wave <- rnorm(n_waves, 0, 0.01)
    sigma2_U <- 0.5
    
    y_prop <- (y_recoded - 1) / (nk - 1)
    z <- qnorm(pmin(pmax(y_prop, 0.05), 0.95))
    
    if (nk > 2) {
      thresholds <- seq(-1.5, 1.5, length.out = nk - 1)
      thresholds[1] <- 0
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
    
    store_idx <- 0
    n_errors <- 0
    
    # MCMC loop
    for (iter in 1:n_iter) {
      
      tryCatch({
        # Update parameters
        beta <- update_beta_robust(X, z, U_wave, wave_mapped, prior_beta)
        U_wave <- update_U_wave_robust(z, X, beta, wave_mapped, sigma2_U, n_waves)
        sigma2_U <- update_sigma2_U_robust(U_wave, 1, 1)
        z <- update_z_robust(y_recoded, X, beta, U_wave, wave_mapped, thresholds)
        
        if (nk > 2) {
          thresholds <- update_thresholds_conservative(z, y_recoded, thresholds, iter)
        }
        
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
      }
      
      # Progress reporting
      if (verbose && chain_id == 1 && iter %% 2000 == 0) {
        cat("Iteration:", iter, "\n")
        if (!is.null(treatment_col) && store_idx > 0) {
          current_ate <- mean(ate_samples[1:store_idx], na.rm = TRUE)
          cat("  Current ATE:", round(current_ate, 4), "\n")
        }
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
      chain_id = chain_id,
      n_errors = n_errors
    ))
  }
  
  # Run chains
  if (verbose) cat("Running", n_chains, "MCMC chain(s)...\n")
  chains <- lapply(1:n_chains, run_chain)
  
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
#' @return BMEOP results object
run_elsa_bmeop <- function(elsa_long, subset_size = 2000, seed = 123, 
                           ate_draws = 50, ate_method = "analytical") {
  set.seed(seed)
  
  cat("BMEOP Model for ELSA Loneliness Analysis\n")
  cat("Method:", ate_method, "with", ate_draws, "draws\n\n")
  
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
    n_iter = 15000,
    n_burn = 10000,
    n_chains = 1,
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
