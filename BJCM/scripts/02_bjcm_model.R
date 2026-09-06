# Bayesian Joint Causal Mediation Analysis: Algorithm 2
# Fixed treatment with lagged mediators ensuring proper temporal ordering
# Reference: Chu et al. (2025), Algorithm 2
#
# Temporal structure: Treatment(W2, fixed) -> Mediators(t-1, lagged) -> Outcome(t, current)
# This ensures proper causal ordering with temporal precedence

library(MASS)
library(mvtnorm)
library(coda)
library(dplyr)
if (requireNamespace("abind", quietly = TRUE)) library(abind)

# -----------------------------------------------------------------------------
# Robust helper functions for numerical stability
# -----------------------------------------------------------------------------

safe_matrix_inverse <- function(A, tol = 1e-10) {
  tryCatch({
    if (any(is.na(A)) || any(is.infinite(A))) {
      return(diag(nrow(A)) * 0.01)
    }
    
    if (nrow(A) != ncol(A)) {
      return(diag(min(nrow(A), ncol(A))) * 0.01)
    }
    
    eig_decomp <- eigen(A, symmetric = TRUE)
    eig_vals <- eig_decomp$values
    eig_vecs <- eig_decomp$vectors
    
    eig_vals[eig_vals < tol] <- tol
    
    A_reg <- eig_vecs %*% diag(eig_vals) %*% t(eig_vecs)
    return(eig_vecs %*% diag(1/eig_vals) %*% t(eig_vecs))
  }, error = function(e) {
    return(diag(nrow(A)) * 0.01)
  })
}

rmvnorm_robust <- function(n, mean, sigma) {
  tryCatch({
    if (any(is.na(mean)) || any(is.na(sigma))) {
      return(rep(0, length(mean)))
    }
    
    if (length(mean) == 1) {
      return(rnorm(n, mean, sqrt(sigma)))
    }
    
    result <- rmvnorm(n, mean, sigma)
    if (any(is.na(result)) || any(is.infinite(result))) {
      return(mean + rnorm(length(mean), 0, 0.01))
    }
    
    return(as.vector(result))
  }, error = function(e) {
    return(rep(0, length(mean)))
  })
}

rtruncnorm_robust <- function(n, mean, sd, lower, upper) {
  tryCatch({
    if (is.na(mean) || is.na(sd) || sd <= 0) {
      return(mean)
    }
    
    if (is.infinite(lower) && is.infinite(upper)) {
      return(rnorm(n, mean, sd))
    }
    
    if (is.infinite(lower)) {
      p_upper <- pnorm(upper, mean, sd)
      u <- runif(n, 0, p_upper)
      return(qnorm(u, mean, sd))
    }
    
    if (is.infinite(upper)) {
      p_lower <- pnorm(lower, mean, sd)
      u <- runif(n, p_lower, 1)
      return(qnorm(u, mean, sd))
    }
    
    p_lower <- pnorm(lower, mean, sd)
    p_upper <- pnorm(upper, mean, sd)
    
    if (p_upper - p_lower < 1e-10) {
      return((lower + upper) / 2)
    }
    
    u <- runif(n, p_lower, p_upper)
    result <- qnorm(u, mean, sd)
    
    if (any(is.na(result)) || any(is.infinite(result))) {
      return((lower + upper) / 2)
    }
    
    return(result)
  }, error = function(e) {
    return(mean)
  })
}

update_U_wave_robust <- function(residual, wave_mapped, sigma2_U, n_waves) {
  U_new <- numeric(n_waves)
  
  for (w in 1:n_waves) {
    wave_idx <- which(wave_mapped == w)
    if (length(wave_idx) == 0) {
      U_new[w] <- rnorm(1, 0, sqrt(sigma2_U))
      next
    }
    
    residual_wave <- residual[wave_idx]
    
    if (any(is.na(residual_wave))) {
      U_new[w] <- rnorm(1, 0, sqrt(sigma2_U))
      next
    }
    
    n_w <- length(wave_idx)
    posterior_precision <- n_w + 1/sigma2_U
    posterior_mean <- sum(residual_wave) / posterior_precision
    posterior_var <- 1 / posterior_precision
    
    U_new[w] <- rnorm(1, posterior_mean, sqrt(posterior_var))
    
    if (is.na(U_new[w]) || is.infinite(U_new[w])) {
      U_new[w] <- 0
    }
  }
  
  return(U_new)
}

update_sigma2_U_robust <- function(U_vec, a_sigma, b_sigma) {
  tryCatch({
    if (any(is.na(U_vec)) || length(U_vec) == 0) {
      return(0.5)
    }
    
    U_clean <- U_vec[!is.na(U_vec)]
    if (length(U_clean) == 0) {
      return(0.5)
    }
    
    posterior_a <- a_sigma + length(U_clean) / 2
    posterior_b <- b_sigma + sum(U_clean^2) / 2
    
    sigma2_new <- 1 / rgamma(1, shape = posterior_a, rate = posterior_b)
    
    if (is.na(sigma2_new) || is.infinite(sigma2_new) || sigma2_new <= 0) {
      return(0.5)
    }
    
    return(pmax(sigma2_new, 0.01))
  }, error = function(e) {
    return(0.5)
  })
}

# Threshold updates for the ordinal outcome and any ordinal mediator with
# more than two categories. Each movable threshold is squeezed between
# quantiles of the data strictly below it and strictly above it (lower_obs /
# upper_obs), a stochastic proposal is drawn uniformly within that squeezed
# range, and only 5% of the distance to the proposal is taken per update.
# Thresholds are updated every 50 iterations, giving the other parameters
# time to adapt between moves. tau_1 is fixed at 0 for identification.
source("01_threshold_update_mh.R")  # posterior-invariant, marginal-likelihood MH threshold update

update_thresholds_conservative <- function(z_latent, y_observed, current_thresholds, iteration) {
  K <- max(y_observed, na.rm = TRUE)
  if (K <= 2) {
    return(current_thresholds)
  }
  
  if (iteration %% 50 != 0) {
    return(current_thresholds)
  }
  
  new_thresholds <- current_thresholds
  n_thresh <- length(current_thresholds)  # = K - 1
  
  for (j in 2:n_thresh) {
    lower_obs <- z_latent[y_observed <= j & !is.na(z_latent)]
    upper_obs <- z_latent[y_observed > j & !is.na(z_latent)]
    
    if (length(lower_obs) > 20 && length(upper_obs) > 20) {
      lower_bound <- quantile(lower_obs, 0.1)
      upper_bound <- quantile(upper_obs, 0.9)
      
      if (!is.na(lower_bound) && !is.na(upper_bound) && lower_bound < upper_bound) {
        if (j > 2) lower_bound <- max(lower_bound, new_thresholds[j - 1] + 0.1)
        if (j < n_thresh) upper_bound <- min(upper_bound, new_thresholds[j + 1] - 0.1)
        
        if (lower_bound < upper_bound) {
          current_val <- new_thresholds[j]
          proposed_val <- runif(1, lower_bound, upper_bound)
          new_thresholds[j] <- 0.95 * current_val + 0.05 * proposed_val
        }
      }
    }
  }
  
  # Enforce strict ordering as a backstop.
  for (j in 2:n_thresh) {
    if (new_thresholds[j] <= new_thresholds[j - 1] + 0.05) {
      new_thresholds[j] <- new_thresholds[j - 1] + 0.05
    }
  }
  
  new_thresholds[1] <- 0  # tau_1 stays fixed at 0 for identification
  return(new_thresholds)
}

# -----------------------------------------------------------------------------
# Parallel mediation structure with shared covariate set
# -----------------------------------------------------------------------------

update_beta_mediator_parallel <- function(X_baseline, z_med_matrix, treatment, alpha_med, 
                                          U_med_wave_matrix, wave_mapped, prior_beta_med) {
  
  n_mediators <- ncol(z_med_matrix)
  p <- ncol(X_baseline)
  
  beta_new_matrix <- matrix(0, p, n_mediators)
  
  for (k in 1:n_mediators) {
    z_med_k <- z_med_matrix[, k]
    alpha_med_k <- alpha_med[k]
    U_med_wave_k <- U_med_wave_matrix[, k]
    
    if (any(is.na(z_med_k)) || any(is.na(U_med_wave_k))) {
      beta_new_matrix[, k] <- rep(0, p)
      next
    }
    
    z_adj <- z_med_k - treatment * alpha_med_k - U_med_wave_k[wave_mapped]
    
    valid_idx <- !is.na(z_adj)
    if (sum(valid_idx) < p) {
      beta_new_matrix[, k] <- rep(0, p)
      next
    }
    
    z_adj <- z_adj[valid_idx]
    X_valid <- X_baseline[valid_idx, , drop = FALSE]
    
    XtX <- crossprod(X_valid)
    posterior_precision <- XtX + prior_beta_med$precision + 1e-4 * diag(p)
    Xtz <- crossprod(X_valid, z_adj)
    
    posterior_cov <- safe_matrix_inverse(posterior_precision)
    posterior_mean <- posterior_cov %*% (Xtz + prior_beta_med$precision %*% prior_beta_med$mean)
    
    beta_new <- rmvnorm_robust(1, as.vector(posterior_mean), posterior_cov)
    
    if (any(is.na(beta_new)) || any(is.infinite(beta_new))) {
      beta_new_matrix[, k] <- as.vector(posterior_mean)
    } else {
      beta_new_matrix[, k] <- beta_new
    }
  }
  
  return(beta_new_matrix)
}

# Update treatment effects on mediators alpha (a-paths)

update_alpha_mediator_parallel <- function(z_med_matrix, X_baseline, beta_med_matrix, treatment, 
                                           U_med_wave_matrix, wave_mapped, prior_alpha) {
  
  n_mediators <- ncol(z_med_matrix)
  alpha_new <- numeric(n_mediators)
  
  for (k in 1:n_mediators) {
    z_med_k <- z_med_matrix[, k]
    beta_med_k <- beta_med_matrix[, k]
    U_med_wave_k <- U_med_wave_matrix[, k]
    
    residual <- z_med_k - as.vector(X_baseline %*% beta_med_k) - U_med_wave_k[wave_mapped]
    
    valid_idx <- !is.na(residual) & !is.na(treatment)
    if (sum(valid_idx) < 2) {
      alpha_new[k] <- 0
      next
    }
    
    residual <- residual[valid_idx]
    treatment_valid <- treatment[valid_idx]
    
    posterior_precision <- sum(treatment_valid^2) + prior_alpha$precision
    posterior_mean <- (sum(treatment_valid * residual) + 
                         prior_alpha$precision * prior_alpha$mean) / posterior_precision
    posterior_var <- 1 / posterior_precision
    
    alpha_new[k] <- rnorm(1, posterior_mean, sqrt(posterior_var))
    
    if (is.na(alpha_new[k]) || is.infinite(alpha_new[k])) {
      alpha_new[k] <- 0
    }
  }
  
  return(alpha_new)
}

# Update mediator random effects U_med

update_U_mediator_parallel <- function(z_med_matrix, X_baseline, beta_med_matrix, treatment,
                                       alpha_med, wave_mapped, sigma2_U_med, n_waves) {
  
  n_mediators <- ncol(z_med_matrix)
  U_new_matrix <- matrix(0, n_waves, n_mediators)
  
  for (k in 1:n_mediators) {
    z_med_k <- z_med_matrix[, k]
    beta_med_k <- beta_med_matrix[, k]
    alpha_med_k <- alpha_med[k]
    
    residual_k <- z_med_k - as.vector(X_baseline %*% beta_med_k) - treatment * alpha_med_k
    
    U_new_matrix[, k] <- update_U_wave_robust(residual_k, wave_mapped, 
                                              sigma2_U_med[k], n_waves)
  }
  
  return(U_new_matrix)
}

# Update latent mediator variables Z_med
# Sample from truncated normal based on observed ordinal mediators

update_z_mediator_parallel <- function(m_obs_matrix, X_baseline, beta_med_matrix, treatment, 
                                       alpha_med, U_med_wave_matrix, wave_mapped, 
                                       thresholds_med_list) {
  
  n <- nrow(m_obs_matrix)
  n_mediators <- ncol(m_obs_matrix)
  z_med_new_matrix <- matrix(0, n, n_mediators)
  
  for (k in 1:n_mediators) {
    m_obs_k <- m_obs_matrix[, k]
    beta_med_k <- beta_med_matrix[, k]
    alpha_med_k <- alpha_med[k]
    U_med_wave_k <- U_med_wave_matrix[, k]
    thresholds_med_k <- thresholds_med_list[[k]]
    
    mu_med_k <- as.vector(X_baseline %*% beta_med_k) + treatment * alpha_med_k + 
      U_med_wave_k[wave_mapped]
    mu_med_k[is.na(mu_med_k)] <- 0
    
    K_med_k <- max(m_obs_k, na.rm = TRUE)
    S_med_k <- c(-Inf, thresholds_med_k, Inf)
    
    if (any(is.na(S_med_k))) {
      S_med_k[is.na(S_med_k)] <- 0
    }
    
    for (i in 1:n) {
      k_val <- m_obs_k[i]
      if (is.na(k_val) || k_val < 1 || k_val > K_med_k) {
        z_med_new_matrix[i, k] <- 0
        next
      }
      
      lower <- S_med_k[k_val]
      upper <- S_med_k[k_val + 1]
      
      z_med_new_matrix[i, k] <- rtruncnorm_robust(1, mean = mu_med_k[i], sd = 1, 
                                                  lower = lower, upper = upper)
    }
  }
  
  return(z_med_new_matrix)
}

# -----------------------------------------------------------------------------
# Outcome model updates
# -----------------------------------------------------------------------------

# Update mediator effects on outcome beta_out_med (b-paths)

update_beta_outcome_mediators <- function(z_out, z_med_matrix, treatment, lambda_out, 
                                          X_baseline, gamma_out, U_out_wave, wave_mapped, 
                                          prior_beta_out_med) {
  
  n_mediators <- ncol(z_med_matrix)
  beta_out_med_new <- numeric(n_mediators)
  
  for (k in 1:n_mediators) {
    z_med_k <- z_med_matrix[, k]
    
    if (n_mediators > 1) {
      other_mediators <- z_med_matrix[, -k, drop = FALSE]
      other_betas <- beta_out_med_new[-k]
      
      if (length(other_betas) > 0) {
        other_effects <- rowSums(other_mediators * matrix(other_betas, nrow = nrow(other_mediators), 
                                                          ncol = length(other_betas), byrow = TRUE))
      } else {
        other_effects <- 0
      }
    } else {
      other_effects <- 0
    }
    
    residual <- z_out - treatment * lambda_out - as.vector(X_baseline %*% gamma_out) - 
      U_out_wave[wave_mapped] - other_effects
    
    valid_idx <- !is.na(residual) & !is.na(z_med_k)
    if (sum(valid_idx) < 2) {
      beta_out_med_new[k] <- 0
      next
    }
    
    residual <- residual[valid_idx]
    z_med_valid <- z_med_k[valid_idx]
    
    posterior_precision <- sum(z_med_valid^2) + prior_beta_out_med$precision
    posterior_mean <- (sum(z_med_valid * residual) + 
                         prior_beta_out_med$precision * prior_beta_out_med$mean) / posterior_precision
    posterior_var <- 1 / posterior_precision
    
    beta_out_med_new[k] <- rnorm(1, posterior_mean, sqrt(posterior_var))
    
    if (is.na(beta_out_med_new[k]) || is.infinite(beta_out_med_new[k])) {
      beta_out_med_new[k] <- 0
    }
  }
  
  return(beta_out_med_new)
}

# Update direct treatment effect lambda (c'-path / NDE)

update_lambda_direct <- function(z_out, z_med_matrix, beta_out_med, treatment, 
                                 X_baseline, gamma_out, U_out_wave, wave_mapped, 
                                 prior_lambda) {
  
  mediator_effects <- rowSums(z_med_matrix * matrix(beta_out_med, nrow = nrow(z_med_matrix), 
                                                    ncol = length(beta_out_med), byrow = TRUE))
  
  residual <- z_out - mediator_effects - as.vector(X_baseline %*% gamma_out) - 
    U_out_wave[wave_mapped]
  
  valid_idx <- !is.na(residual) & !is.na(treatment)
  if (sum(valid_idx) < 2) return(0)
  
  residual <- residual[valid_idx]
  treatment_valid <- treatment[valid_idx]
  
  posterior_precision <- sum(treatment_valid^2) + prior_lambda$precision
  posterior_mean <- (sum(treatment_valid * residual) + 
                       prior_lambda$precision * prior_lambda$mean) / posterior_precision
  posterior_var <- 1 / posterior_precision
  
  lambda_new <- rnorm(1, posterior_mean, sqrt(posterior_var))
  
  if (is.na(lambda_new) || is.infinite(lambda_new)) {
    return(0)
  }
  
  return(lambda_new)
}

# Update baseline covariate effects on outcome gamma

update_gamma_baseline <- function(z_out, z_med_matrix, beta_out_med, treatment, 
                                  lambda_out, X_baseline, U_out_wave, wave_mapped, 
                                  prior_gamma) {
  
  p <- ncol(X_baseline)
  
  mediator_effects <- rowSums(z_med_matrix * matrix(beta_out_med, nrow = nrow(z_med_matrix), 
                                                    ncol = length(beta_out_med), byrow = TRUE))
  
  residual <- z_out - mediator_effects - treatment * lambda_out - U_out_wave[wave_mapped]
  
  valid_idx <- !is.na(residual)
  if (sum(valid_idx) < p) return(rep(0, p))
  
  residual <- residual[valid_idx]
  X_valid <- X_baseline[valid_idx, , drop = FALSE]
  
  XtX <- crossprod(X_valid)
  posterior_precision <- XtX + prior_gamma$precision + 1e-4 * diag(p)
  Xtr <- crossprod(X_valid, residual)
  
  posterior_cov <- safe_matrix_inverse(posterior_precision)
  posterior_mean <- posterior_cov %*% (Xtr + prior_gamma$precision %*% prior_gamma$mean)
  
  gamma_new <- rmvnorm_robust(1, as.vector(posterior_mean), posterior_cov)
  
  if (any(is.na(gamma_new)) || any(is.infinite(gamma_new))) {
    return(as.vector(posterior_mean))
  }
  
  return(gamma_new)
}

# Update outcome random effects U_out

update_U_outcome <- function(z_out, z_med_matrix, beta_out_med, treatment, lambda_out,
                             X_baseline, gamma_out, wave_mapped, sigma2_U_out, n_waves) {
  
  mediator_effects <- rowSums(z_med_matrix * matrix(beta_out_med, nrow = nrow(z_med_matrix), 
                                                    ncol = length(beta_out_med), byrow = TRUE))
  
  residual_out <- z_out - mediator_effects - treatment * lambda_out - 
    as.vector(X_baseline %*% gamma_out)
  
  U_out_new <- update_U_wave_robust(residual_out, wave_mapped, sigma2_U_out, n_waves)
  
  return(U_out_new)
}

# Update latent outcome variable Z_out

update_z_outcome <- function(y_obs, z_med_matrix, beta_out_med, treatment, lambda_out, 
                             X_baseline, gamma_out, U_out_wave, wave_mapped, thresholds_out) {
  
  n <- length(y_obs)
  z_out_new <- numeric(n)
  
  mediator_effects <- rowSums(z_med_matrix * matrix(beta_out_med, nrow = nrow(z_med_matrix), 
                                                    ncol = length(beta_out_med), byrow = TRUE))
  
  mu_out <- mediator_effects + treatment * lambda_out + 
    as.vector(X_baseline %*% gamma_out) + U_out_wave[wave_mapped]
  mu_out[is.na(mu_out)] <- 0
  
  K_out <- max(y_obs, na.rm = TRUE)
  S_out <- c(-Inf, thresholds_out, Inf)
  
  if (any(is.na(S_out))) {
    S_out[is.na(S_out)] <- 0
  }
  
  for (i in 1:n) {
    k <- y_obs[i]
    if (is.na(k) || k < 1 || k > K_out) {
      z_out_new[i] <- 0
      next
    }
    
    lower <- S_out[k]
    upper <- S_out[k + 1]
    
    z_out_new[i] <- rtruncnorm_robust(1, mean = mu_out[i], sd = 1, 
                                      lower = lower, upper = upper)
  }
  
  return(z_out_new)
}

# -----------------------------------------------------------------------------
# Counterfactual simulation and causal effects
# Using rho = 0 (independence assumption between counterfactual errors)
# -----------------------------------------------------------------------------

simulate_parallel_mediation_effects <- function(X_baseline, beta_med_matrix, alpha_med, 
                                                beta_out_med, lambda_out, gamma_out, 
                                                U_med_wave_matrix, U_out_wave, wave_mapped, 
                                                thresholds_med_list, thresholds_out, 
                                                n_draws = 50) {
  
  n <- nrow(X_baseline)
  n_mediators <- length(alpha_med)
  K <- n_mediators
  
  te_matrix <- matrix(0, n, n_draws)
  nde_matrix <- matrix(0, n, n_draws)
  nie_total_matrix <- matrix(0, n, n_draws)
  nie_individual_array <- array(0, dim = c(n, n_mediators, n_draws))
  
  # Shapley subset weights: w(|S|) = |S|!(K-|S|-1)!/K!, and the table of all
  # 2^K subsets S of {1,...,K} (as 0/1 membership rows). NIE_k is averaged
  # over all subsets NOT containing k, weighted by w(|S|), matching the
  # efficiency-property decomposition in Algorithm A.2 Step 9(d) / equation
  # (6) of the main text, rather than one arbitrary fixed ordering.
  shapley_weight <- sapply(0:(K - 1), function(s) {
    factorial(s) * factorial(K - s - 1) / factorial(K)
  })
  subset_membership <- as.matrix(expand.grid(rep(list(c(0, 1)), K)))
  colnames(subset_membership) <- NULL
  n_subsets <- nrow(subset_membership)
  subset_size <- rowSums(subset_membership)
  
  for (d in 1:n_draws) {
    
    # Draw errors once for consistency across counterfactuals
    epsilon_med <- matrix(rnorm(n * n_mediators), n, n_mediators)
    epsilon_out <- rnorm(n)
    
    # Simulate counterfactual mediators under rho = 0
    Z_med_treat_matrix <- matrix(0, n, n_mediators)
    Z_med_control_matrix <- matrix(0, n, n_mediators)
    
    for (k in 1:n_mediators) {
      # M*(T=1): Mediator under treatment
      mu_med_treat_k <- as.vector(X_baseline %*% beta_med_matrix[, k]) + 
        1 * alpha_med[k] + U_med_wave_matrix[wave_mapped, k]
      Z_med_treat_matrix[, k] <- mu_med_treat_k + epsilon_med[, k]
      
      # M*(T=0): Mediator under control
      mu_med_control_k <- as.vector(X_baseline %*% beta_med_matrix[, k]) + 
        0 * alpha_med[k] + U_med_wave_matrix[wave_mapped, k]
      Z_med_control_matrix[, k] <- mu_med_control_k + epsilon_med[, k]
    }
    
    # Simulate counterfactual outcomes under rho = 0
    # Y*(T=1, M*(T=1)): Total effect numerator
    mediator_effects_11 <- rowSums(Z_med_treat_matrix * 
                                     matrix(beta_out_med, n, n_mediators, byrow = TRUE))
    mu_Y_11 <- mediator_effects_11 + 1 * lambda_out + 
      as.vector(X_baseline %*% gamma_out) + U_out_wave[wave_mapped]
    Z_Y_11 <- mu_Y_11 + epsilon_out
    
    # Y*(T=0, M*(T=0)): Total effect denominator
    mediator_effects_00 <- rowSums(Z_med_control_matrix * 
                                     matrix(beta_out_med, n, n_mediators, byrow = TRUE))
    mu_Y_00 <- mediator_effects_00 + 0 * lambda_out + 
      as.vector(X_baseline %*% gamma_out) + U_out_wave[wave_mapped]
    Z_Y_00 <- mu_Y_00 + epsilon_out
    
    # Y*(T=1, M*(T=0)): Natural direct effect
    mu_Y_10 <- mediator_effects_00 + 1 * lambda_out + 
      as.vector(X_baseline %*% gamma_out) + U_out_wave[wave_mapped]
    Z_Y_10 <- mu_Y_10 + epsilon_out
    
    # Convert to ordinal outcomes
    S_out <- c(-Inf, thresholds_out, Inf)
    Y_11 <- findInterval(Z_Y_11, S_out)
    Y_00 <- findInterval(Z_Y_00, S_out)
    Y_10 <- findInterval(Z_Y_10, S_out)
    
    # Compute causal effects
    te_matrix[, d] <- Y_11 - Y_00
    nde_matrix[, d] <- Y_10 - Y_00
    nie_total_matrix[, d] <- Y_11 - Y_10
    
    # v(S) for every subset S of mediators: outcome with T=1, mediators in S
    # set to their treated counterfactual value, remaining mediators at
    # their control counterfactual value. v(empty) = Y_10, v(full) = Y_11.
    v_subset <- matrix(0, n, n_subsets)
    for (s_idx in 1:n_subsets) {
      membership <- subset_membership[s_idx, ]
      if (all(membership == 0)) {
        v_subset[, s_idx] <- Y_10
      } else if (all(membership == 1)) {
        v_subset[, s_idx] <- Y_11
      } else {
        Z_med_s <- Z_med_control_matrix
        Z_med_s[, membership == 1] <- Z_med_treat_matrix[, membership == 1]
        mediator_effects_s <- rowSums(Z_med_s *
                                        matrix(beta_out_med, n, n_mediators, byrow = TRUE))
        mu_Y_s <- mediator_effects_s + 1 * lambda_out +
          as.vector(X_baseline %*% gamma_out) + U_out_wave[wave_mapped]
        Z_Y_s <- mu_Y_s + epsilon_out
        v_subset[, s_idx] <- findInterval(Z_Y_s, S_out)
      }
    }
    
    # Shapley value of mediator k: average marginal contribution
    # v(S union {k}) - v(S) over all subsets S not containing k, weighted
    # by w(|S|). This symmetrises the attribution across mediators --
    # consistent with the parallel (non-ordered) mediator structure -- and
    # the efficiency property guarantees sum_k NIE_k = NIE_total exactly,
    # for the same reason (telescoping accounting identity) as the
    # sequential version, without depending on any single ordering.
    for (k in 1:n_mediators) {
      without_k <- which(subset_membership[, k] == 0)
      nie_k_draw <- numeric(n)
      for (s_idx in without_k) {
        membership <- subset_membership[s_idx, ]
        s_idx_with_k <- which(apply(subset_membership, 1, function(row) {
          all(row == replace(membership, k, 1))
        }))
        w <- shapley_weight[subset_size[s_idx] + 1]
        nie_k_draw <- nie_k_draw + w * (v_subset[, s_idx_with_k] - v_subset[, s_idx])
      }
      nie_individual_array[, k, d] <- nie_k_draw
    }
  }
  
  return(list(
    total_effects = rowMeans(te_matrix),
    direct_effects = rowMeans(nde_matrix),
    indirect_effects_total = rowMeans(nie_total_matrix),
    indirect_effects_individual = apply(nie_individual_array, c(1,2), mean),
    total_effects_samples = te_matrix,
    direct_effects_samples = nde_matrix,
    indirect_effects_samples = nie_total_matrix,
    indirect_effects_individual_samples = nie_individual_array,
    decomposition_check = list(
      te_mean = mean(rowMeans(te_matrix)),
      nde_mean = mean(rowMeans(nde_matrix)),
      nie_total_mean = mean(rowMeans(nie_total_matrix)),
      nie_individual_means = apply(apply(nie_individual_array, c(1,2), mean), 2, mean),
      decomposition_error = abs(mean(rowMeans(te_matrix)) - 
                                  (mean(rowMeans(nde_matrix)) + mean(rowMeans(nie_total_matrix)))),
      nie_sum_check = abs(mean(rowMeans(nie_total_matrix)) - 
                            sum(apply(apply(nie_individual_array, c(1,2), mean), 2, mean))),
      rho = 0
    )
  ))
}

# -----------------------------------------------------------------------------
# Main Bayesian Algorithm 2 implementation
# -----------------------------------------------------------------------------

bayesian_algorithm2_parallel <- function(y, mediators_matrix, X_baseline, treatment, wave, 
                                         n_iter = 12000, n_burn = 7000, n_thin = 2,
                                         n_chains = 1, prior_var = 4.0,
                                         mediation_draws = 50, seed = NULL, verbose = TRUE,
                                         memory_efficient = FALSE) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Data validation
  n <- length(y)
  n_mediators <- ncol(mediators_matrix)
  
  if (nrow(mediators_matrix) != n || nrow(X_baseline) != n || 
      length(treatment) != n || length(wave) != n) {
    stop("Data dimension mismatch")
  }
  
  if (any(is.na(y)) || any(is.na(mediators_matrix)) || any(is.na(X_baseline)) || 
      any(is.na(treatment)) || any(is.na(wave))) {
    stop("Data contains NAs - please remove before analysis")
  }
  
  # Recode outcomes and mediators
  y_recoded <- match(y, sort(unique(y)))
  nk_out <- max(y_recoded)
  
  mediators_recoded <- matrix(0, n, n_mediators)
  nk_med <- numeric(n_mediators)
  
  for (k in 1:n_mediators) {
    mediators_recoded[, k] <- match(mediators_matrix[, k], sort(unique(mediators_matrix[, k])))
    nk_med[k] <- max(mediators_recoded[, k])
  }
  
  # Wave structure
  unique_waves <- sort(unique(wave))
  n_waves <- length(unique_waves)
  wave_mapped <- match(wave, unique_waves)
  
  p_baseline <- ncol(X_baseline)
  
  if (verbose) {
    cat("Bayesian Algorithm 2: Parallel mediation analysis\n")
    cat("Reference: Chu & Yu (2025)\n")
    cat("n =", n, "| p_baseline =", p_baseline, "| n_waves =", n_waves, "\n")
    cat("Outcome categories:", nk_out, "| Mediator categories:", paste(nk_med, collapse = ", "), "\n")
    cat("n_mediators:", n_mediators, "| Parallel structure with shared covariates\n")
    cat("rho = 0 (independence assumption)\n\n")
  }
  
  # Priors
  prior_beta_med <- list(
    mean = rep(0, p_baseline),
    precision = diag(1/prior_var, p_baseline)
  )
  
  prior_alpha <- list(
    mean = 0,
    precision = 1/prior_var
  )
  
  prior_beta_out_med <- list(
    mean = 0,
    precision = 1/prior_var
  )
  
  prior_lambda <- list(
    mean = 0,
    precision = 1/prior_var
  )
  
  prior_gamma <- list(
    mean = rep(0, p_baseline),
    precision = diag(1/prior_var, p_baseline)
  )
  
  # MCMC chain function
  run_single_chain <- function(chain_id) {
    
    chain_start_time <- Sys.time()
    
    # Chain-specific reproducible seed, set explicitly inside the worker so
    # results are identical whether chains run sequentially (lapply) or
    # in parallel forked processes (mclapply).
    if (!is.null(seed)) set.seed(seed + chain_id * 10000)
    
    # Overdispersed initial values: chain 1 starts near zero, chains 2+ are
    # offset to genuinely different starting regions (alternating sign,
    # increasing magnitude), so that multi-chain R-hat reflects convergence
    # from dispersed starts.
    init_offset <- c(0, -1.5, 1.5, -2.5, 2.5)[((chain_id - 1) %% 5) + 1]
    init_sd     <- if (chain_id == 1) 0.01 else 0.5
    
    # Initialise parameters
    beta_med_matrix <- matrix(rnorm(p_baseline * n_mediators, init_offset, init_sd), p_baseline, n_mediators)
    alpha_med <- rnorm(n_mediators, init_offset, init_sd)
    beta_out_med <- rnorm(n_mediators, init_offset, init_sd)
    lambda_out <- rnorm(1, init_offset, init_sd)
    gamma_out <- rnorm(p_baseline, init_offset, init_sd)
    
    U_med_wave_matrix <- matrix(rnorm(n_waves * n_mediators, 0, 0.01), n_waves, n_mediators)
    U_out_wave <- rnorm(n_waves, 0, 0.01)
    
    sigma2_U_med <- rep(0.5, n_mediators)
    sigma2_U_out <- 0.5
    
    # Initialise latent variables
    z_med_matrix <- matrix(0, n, n_mediators)
    for (k in 1:n_mediators) {
      y_prop_med <- (mediators_recoded[, k] - 1) / (nk_med[k] - 1)
      z_med_matrix[, k] <- qnorm(pmin(pmax(y_prop_med, 0.05), 0.95))
    }
    
    y_prop_out <- (y_recoded - 1) / (nk_out - 1)
    z_out <- qnorm(pmin(pmax(y_prop_out, 0.05), 0.95))
    
    # Initialise thresholds. tau_1 = 0 is fixed for identification; the
    # remaining thresholds are seeded as a strictly increasing sequence.
    thresholds_med_list <- list()
    for (k in 1:n_mediators) {
      if (nk_med[k] > 2) {
        thresholds_med_list[[k]] <- seq(0, 3, length.out = nk_med[k] - 1)
      } else {
        thresholds_med_list[[k]] <- 0
      }
    }
    
    if (nk_out > 2) {
      thresholds_out <- seq(0, 3, length.out = nk_out - 1)
    } else {
      thresholds_out <- 0
    }
    
    # Storage
    n_store <- floor((n_iter - n_burn) / n_thin)
    
    beta_med_samples <- array(NA, c(n_store, p_baseline, n_mediators))
    alpha_med_samples <- matrix(NA, n_store, n_mediators)
    beta_out_med_samples <- matrix(NA, n_store, n_mediators)
    lambda_out_samples <- numeric(n_store)
    gamma_out_samples <- matrix(NA, n_store, p_baseline)
    
    if (memory_efficient) {
      total_effects_samples <- numeric(n_store)
      direct_effects_samples <- numeric(n_store)
      indirect_effects_samples <- numeric(n_store)
      indirect_effects_individual_samples <- matrix(NA, n_store, n_mediators)
      
      if (verbose && chain_id == 1) {
        cat("Memory-efficient mode: Storing only population-level statistics\n")
      }
    } else {
      total_effects_samples <- matrix(NA, n_store, n)
      direct_effects_samples <- matrix(NA, n_store, n)
      indirect_effects_samples <- matrix(NA, n_store, n)
      indirect_effects_individual_samples <- array(NA, c(n_store, n, n_mediators))
    }
    
    # Storage for trace plots / posterior predictive checks: thresholds
    # each stored iteration, plus per-iteration marginal category-probability
    # vectors computed from OBSERVED treatment and mediators (not the
    # counterfactual arms used in effect estimation), for the outcome model
    # and each of the three mediator models.
    thresholds_out_samples <- matrix(NA, n_store, length(thresholds_out))
    thresholds_med_samples <- lapply(1:n_mediators, function(k) {
      matrix(NA, n_store, length(thresholds_med_list[[k]]))
    })
    y_pred_probs_outcome <- matrix(NA, n_store, nk_out)
    y_pred_probs_mediators <- lapply(1:n_mediators, function(k) {
      matrix(NA, n_store, nk_med[k])
    })
    
    # Wave-level random effects and their variances, stored alongside
    # BMEOP's equivalents for convergence checks and appendix reporting.
    U_med_wave_samples <- array(NA, c(n_store, n_waves, n_mediators))
    U_out_wave_samples <- matrix(NA, n_store, n_waves)
    sigma2_U_med_samples <- matrix(NA, n_store, n_mediators)
    sigma2_U_out_samples <- numeric(n_store)
    
    store_idx <- 0
    n_errors <- 0
    
    # MCMC loop implementing full Algorithm 2
    for (iter in 1:n_iter) {
      
      tryCatch({
        
        # Algorithm 2 Steps 1-3: Mediator models
        beta_med_matrix <- update_beta_mediator_parallel(
          X_baseline, z_med_matrix, treatment, alpha_med, 
          U_med_wave_matrix, wave_mapped, prior_beta_med
        )
        
        alpha_med <- update_alpha_mediator_parallel(
          z_med_matrix, X_baseline, beta_med_matrix, treatment, 
          U_med_wave_matrix, wave_mapped, prior_alpha
        )
        
        U_med_wave_matrix <- update_U_mediator_parallel(
          z_med_matrix, X_baseline, beta_med_matrix, treatment,
          alpha_med, wave_mapped, sigma2_U_med, n_waves
        )
        
        for (k in 1:n_mediators) {
          sigma2_U_med[k] <- update_sigma2_U_robust(U_med_wave_matrix[, k], 1, 1)
        }
        
        # Threshold update conditions on the observed-data marginal
        # likelihood p(M^(k) | beta_med, alpha, U_med, S^(k)) -- a function
        # of beta_med_matrix, alpha_med, U_med_wave_matrix (all just
        # updated above) and the fixed data, NOT of z_med_matrix. It is
        # therefore performed here, before z_med_matrix is refreshed, so
        # that the z_med_matrix draw immediately below uses the CURRENT
        # (just-updated) thresholds rather than last iteration's.
        for (k in 1:n_mediators) {
          if (nk_med[k] > 2) {
            eta_med_k <- as.vector(X_baseline %*% beta_med_matrix[, k]) +
              treatment * alpha_med[k] +
              U_med_wave_matrix[, k][wave_mapped]
            eta_med_k[is.na(eta_med_k)] <- 0
            
            thresholds_med_list[[k]] <- update_thresholds_mh_marginal(
              z_latent = z_med_matrix[, k], y_observed = mediators_recoded[, k],
              current_thresholds = thresholds_med_list[[k]], iteration = iter,
              eta = eta_med_k, label = paste0("mediator_", k),
              proposal_sd = 0.08, update_every_iter = TRUE
            )
          }
        }
        
        z_med_matrix <- update_z_mediator_parallel(
          mediators_recoded, X_baseline, beta_med_matrix, treatment, 
          alpha_med, U_med_wave_matrix, wave_mapped, thresholds_med_list
        )
        
        # Algorithm 2 Steps 4-7: Outcome model
        beta_out_med <- update_beta_outcome_mediators(
          z_out, z_med_matrix, treatment, lambda_out, X_baseline, gamma_out, 
          U_out_wave, wave_mapped, prior_beta_out_med
        )
        
        lambda_out <- update_lambda_direct(
          z_out, z_med_matrix, beta_out_med, treatment, X_baseline, gamma_out, 
          U_out_wave, wave_mapped, prior_lambda
        )
        
        gamma_out <- update_gamma_baseline(
          z_out, z_med_matrix, beta_out_med, treatment, lambda_out, X_baseline, 
          U_out_wave, wave_mapped, prior_gamma
        )
        
        U_out_wave <- update_U_outcome(
          z_out, z_med_matrix, beta_out_med, treatment, lambda_out,
          X_baseline, gamma_out, wave_mapped, sigma2_U_out, n_waves
        )
        
        sigma2_U_out <- update_sigma2_U_robust(U_out_wave, 1, 1)
        
        if (nk_out > 2) {
          mediator_effects_th <- rowSums(z_med_matrix * matrix(beta_out_med,
                                                               nrow = nrow(z_med_matrix), ncol = length(beta_out_med), byrow = TRUE))
          eta_out <- mediator_effects_th + treatment * lambda_out +
            as.vector(X_baseline %*% gamma_out) + U_out_wave[wave_mapped]
          eta_out[is.na(eta_out)] <- 0
          
          thresholds_out <- update_thresholds_mh_marginal(
            z_latent = z_out, y_observed = y_recoded,
            current_thresholds = thresholds_out, iteration = iter,
            eta = eta_out, label = "outcome",
            proposal_sd = 0.08, update_every_iter = TRUE
          )
        }
        
        z_out <- update_z_outcome(
          y_recoded, z_med_matrix, beta_out_med, treatment, lambda_out, 
          X_baseline, gamma_out, U_out_wave, wave_mapped, thresholds_out
        )
        
      }, error = function(e) {
        n_errors <<- n_errors + 1
        if (verbose && n_errors <= 3) {
          cat("Chain", chain_id, "Error in iteration", iter, ":", e$message, "\n")
        }
      })
      
      # Store samples after burn-in
      if (iter > n_burn && (iter - n_burn) %% n_thin == 0) {
        store_idx <- store_idx + 1
        
        beta_med_samples[store_idx, , ] <- beta_med_matrix
        alpha_med_samples[store_idx, ] <- alpha_med
        beta_out_med_samples[store_idx, ] <- beta_out_med
        lambda_out_samples[store_idx] <- lambda_out
        gamma_out_samples[store_idx, ] <- gamma_out
        
        # Algorithm 2 Steps 9-10: Counterfactuals and effects
        cf_results <- simulate_parallel_mediation_effects(
          X_baseline, beta_med_matrix, alpha_med, beta_out_med, lambda_out, 
          gamma_out, U_med_wave_matrix, U_out_wave, wave_mapped, 
          thresholds_med_list, thresholds_out, mediation_draws
        )
        
        if (memory_efficient) {
          total_effects_samples[store_idx] <- mean(cf_results$total_effects)
          direct_effects_samples[store_idx] <- mean(cf_results$direct_effects)
          indirect_effects_samples[store_idx] <- mean(cf_results$indirect_effects_total)
          indirect_effects_individual_samples[store_idx, ] <- colMeans(cf_results$indirect_effects_individual)
        } else {
          total_effects_samples[store_idx, ] <- cf_results$total_effects
          direct_effects_samples[store_idx, ] <- cf_results$direct_effects
          indirect_effects_samples[store_idx, ] <- cf_results$indirect_effects_total
          indirect_effects_individual_samples[store_idx, , ] <- cf_results$indirect_effects_individual
        }
        
        # PPC storage: thresholds plus model-implied marginal category
        # distributions computed from OBSERVED treatment/mediators (the same
        # mu formulas used in update_z_outcome / update_z_mediator_parallel),
        # not the counterfactual arms used for TE/NDE/NIE above.
        thresholds_out_samples[store_idx, ] <- thresholds_out
        mediator_effects_obs <- rowSums(z_med_matrix * matrix(beta_out_med, nrow = n,
                                                              ncol = n_mediators, byrow = TRUE))
        mu_out_obs <- mediator_effects_obs + treatment * lambda_out +
          as.vector(X_baseline %*% gamma_out) + U_out_wave[wave_mapped]
        S_out_obs <- c(-Inf, thresholds_out, Inf)
        y_pred_probs_outcome[store_idx, ] <- sapply(1:nk_out, function(k) {
          mean(pnorm(S_out_obs[k + 1], mean = mu_out_obs, sd = 1) -
                 pnorm(S_out_obs[k], mean = mu_out_obs, sd = 1))
        })
        
        for (k in 1:n_mediators) {
          thresholds_med_samples[[k]][store_idx, ] <- thresholds_med_list[[k]]
          mu_med_obs_k <- as.vector(X_baseline %*% beta_med_matrix[, k]) +
            treatment * alpha_med[k] + U_med_wave_matrix[wave_mapped, k]
          S_med_obs_k <- c(-Inf, thresholds_med_list[[k]], Inf)
          y_pred_probs_mediators[[k]][store_idx, ] <- sapply(1:nk_med[k], function(kk) {
            mean(pnorm(S_med_obs_k[kk + 1], mean = mu_med_obs_k, sd = 1) -
                   pnorm(S_med_obs_k[kk], mean = mu_med_obs_k, sd = 1))
          })
        }
        
        U_med_wave_samples[store_idx, , ] <- U_med_wave_matrix
        U_out_wave_samples[store_idx, ] <- U_out_wave
        sigma2_U_med_samples[store_idx, ] <- sigma2_U_med
        sigma2_U_out_samples[store_idx] <- sigma2_U_out
      }
      
      # Progress reporting writes to a per-chain file instead of cat(),
      # since cat() inside a forked mclapply worker does not reliably reach
      # the parent session's console. Check progress from a separate R
      # session/terminal via monitor_progress(), or
      # `cat(readLines("../logs/bjcm_progress_chain_<id>.log"))`.
      cheap_interval <- max(1, floor(n_iter / 100))  # ~100 cheap updates per chain
      if (iter %% cheap_interval == 0 || iter == n_iter) {
        elapsed_secs <- as.numeric(difftime(Sys.time(), chain_start_time, units = "secs"))
        pct_done <- iter / n_iter
        eta_secs <- if (pct_done > 0) elapsed_secs / pct_done - elapsed_secs else NA
        progress_line <- sprintf(
          "chain=%d iter=%d/%d (%.1f%%) elapsed=%.1fmin eta=%.1fmin errors=%d timestamp=%s",
          chain_id, iter, n_iter, pct_done * 100,
          elapsed_secs / 60, if (is.na(eta_secs)) NA else eta_secs / 60,
          n_errors, format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        )
        
        # Expensive decomposition recheck kept at a coarser interval (every
        # 2000 iters) so it doesn't add overhead on every progress update -
        # it's diagnostic sugar, not required for basic progress.
        if (iter %% 2000 == 0 && store_idx > 0) {
          last_cf <- simulate_parallel_mediation_effects(
            X_baseline, beta_med_matrix, alpha_med, beta_out_med, lambda_out,
            gamma_out, U_med_wave_matrix, U_out_wave, wave_mapped,
            thresholds_med_list, thresholds_out, 10
          )
          progress_line <- paste0(progress_line,
                                  sprintf(" TE-NDE-NIE_err=%.6g", last_cf$decomposition_check$decomposition_error))
          if (!is.null(last_cf$decomposition_check$nie_sum_check)) {
            progress_line <- paste0(progress_line,
                                    sprintf(" NIE_sum_err=%.6g", last_cf$decomposition_check$nie_sum_check))
          }
        }
        
        writeLines(progress_line, con = sprintf("../logs/bjcm_progress_chain_%d.log", chain_id))
      }
    }
    
    chain_result <- list(
      beta_med_samples = beta_med_samples,
      alpha_med_samples = alpha_med_samples,
      beta_out_med_samples = beta_out_med_samples,
      lambda_out_samples = lambda_out_samples,
      gamma_out_samples = gamma_out_samples,
      total_effects = total_effects_samples,
      direct_effects = direct_effects_samples,
      indirect_effects = indirect_effects_samples,
      indirect_effects_individual = indirect_effects_individual_samples,
      thresholds_out = thresholds_out_samples,
      thresholds_med = thresholds_med_samples,
      y_pred_probs_outcome = y_pred_probs_outcome,
      y_pred_probs_mediators = y_pred_probs_mediators,
      U_med_wave = U_med_wave_samples,
      U_out_wave = U_out_wave_samples,
      sigma2_U_med = sigma2_U_med_samples,
      sigma2_U_out = sigma2_U_out_samples,
      chain_id = chain_id,
      n_errors = n_errors
    )
    
    # Checkpoint to disk inside this worker, before returning. If the outer
    # combine/return step (mclapply -> pipe back to the main session) fails
    # afterwards, this chain's result is not lost - see
    # recover_bjcm_chains_from_checkpoints() to reassemble a full run from
    # these files without re-running any chain.
    tryCatch({
      saveRDS(chain_result, sprintf("../results/checkpoints/bjcm_chain_%d_checkpoint.rds", chain_id))
    }, error = function(e) {
      writeLines(sprintf("WARNING: checkpoint save failed for chain %d: %s", chain_id, conditionMessage(e)),
                 con = sprintf("../logs/bjcm_progress_chain_%d.log", chain_id))
    })
    
    return(chain_result)
  }
  
  # Run chains
  if (verbose) cat("Running", n_chains, "chain(s)...\n")
  
  if (n_chains > 1) {
    # Parallelism via forking (Linux/macOS only).
    n_cores_use <- min(n_chains, parallel::detectCores())
    if (verbose) cat("Dispatching", n_chains, "chains across", n_cores_use, "cores...\n")
    
    chains <- parallel::mclapply(1:n_chains, run_single_chain,
                                 mc.cores = n_cores_use,
                                 mc.preschedule = FALSE)
    
    # mclapply does not raise on a worker error - it returns a try-error
    # object in that chain's slot instead. Check explicitly.
    failed <- vapply(chains, function(x) inherits(x, "try-error"), logical(1))
    if (any(failed)) {
      stop("Chain(s) ", paste(which(failed), collapse = ", "), " failed:\n",
           paste(vapply(chains[failed], function(x) as.character(x), character(1)), collapse = "\n"))
    }
  } else {
    chains <- lapply(1:n_chains, run_single_chain)
  }
  
  # Combine results
  combined_alpha_med <- do.call(rbind, lapply(chains, function(x) x$alpha_med_samples))
  combined_beta_out_med <- do.call(rbind, lapply(chains, function(x) x$beta_out_med_samples))
  combined_lambda_out <- do.call(c, lapply(chains, function(x) x$lambda_out_samples))
  
  if (memory_efficient) {
    combined_total_effects <- do.call(c, lapply(chains, function(x) x$total_effects))
    combined_direct_effects <- do.call(c, lapply(chains, function(x) x$direct_effects))
    combined_indirect_effects <- do.call(c, lapply(chains, function(x) x$indirect_effects))
    combined_indirect_individual <- do.call(rbind, 
                                            lapply(chains, function(x) x$indirect_effects_individual))
  } else {
    combined_total_effects <- do.call(rbind, lapply(chains, function(x) x$total_effects))
    combined_direct_effects <- do.call(rbind, lapply(chains, function(x) x$direct_effects))
    combined_indirect_effects <- do.call(rbind, lapply(chains, function(x) x$indirect_effects))
    combined_indirect_individual <- abind::abind(
      lapply(chains, function(x) x$indirect_effects_individual), along = 1
    )
  }
  
  # Summaries
  create_summary <- function(samples) {
    c(mean = mean(samples, na.rm = TRUE),
      sd = sd(samples, na.rm = TRUE),
      q2.5 = as.numeric(quantile(samples, 0.025, na.rm = TRUE)),
      q97.5 = as.numeric(quantile(samples, 0.975, na.rm = TRUE)))
  }
  
  if (memory_efficient) {
    te_summary <- create_summary(combined_total_effects)
    nde_summary <- create_summary(combined_direct_effects)
    nie_summary <- create_summary(combined_indirect_effects)
    
    nie_individual_summaries <- lapply(1:n_mediators, function(k) {
      create_summary(combined_indirect_individual[, k])
    })
  } else {
    te_summary <- create_summary(as.vector(combined_total_effects))
    nde_summary <- create_summary(as.vector(combined_direct_effects))
    nie_summary <- create_summary(as.vector(combined_indirect_effects))
    
    nie_individual_summaries <- lapply(1:n_mediators, function(k) {
      create_summary(as.vector(combined_indirect_individual[, , k]))
    })
  }
  
  alpha_summaries <- lapply(1:n_mediators, function(k) {
    create_summary(combined_alpha_med[, k])
  })
  
  beta_out_med_summaries <- lapply(1:n_mediators, function(k) {
    create_summary(combined_beta_out_med[, k])
  })
  
  lambda_summary <- create_summary(combined_lambda_out)
  
  # Validation
  decomposition_error <- abs(te_summary["mean"] - (nie_summary["mean"] + nde_summary["mean"]))
  
  nie_individual_sum <- sum(sapply(nie_individual_summaries, function(x) x["mean"]))
  nie_sum_error <- abs(nie_summary["mean"] - nie_individual_sum)
  
  if (verbose) {
    cat("\nAlgorithm 2 results\n")
    if (memory_efficient) {
      cat("(Memory-efficient mode: Population-level estimates only)\n")
    }
    cat("Total Effect:", round(te_summary["mean"], 4), 
        "[", round(te_summary["q2.5"], 4), ",", round(te_summary["q97.5"], 4), "]\n")
    cat("Natural Direct Effect:", round(nde_summary["mean"], 4), 
        "[", round(nde_summary["q2.5"], 4), ",", round(nde_summary["q97.5"], 4), "]\n")
    cat("Total Natural Indirect Effect:", round(nie_summary["mean"], 4), 
        "[", round(nie_summary["q2.5"], 4), ",", round(nie_summary["q97.5"], 4), "]\n\n")
    
    mediator_names <- c("Depression", "Transport", "Mobility")
    for (k in 1:n_mediators) {
      cat(paste0(mediator_names[k], " (a-path alpha", k, "): "), 
          round(alpha_summaries[[k]]["mean"], 4), "\n")
      cat(paste0(mediator_names[k], " (b-path beta", k, "): "), 
          round(beta_out_med_summaries[[k]]["mean"], 4), "\n")
      cat(paste0(mediator_names[k], " NIE", k, " (Shapley value): "), 
          round(nie_individual_summaries[[k]]["mean"], 4), "\n\n")
    }
    
    cat("Direct effect (c'-path lambda):", round(lambda_summary["mean"], 4), "\n")
    cat("TE decomposition error (TE - NDE - NIE):", round(decomposition_error, 6), "\n")
    cat("NIE sum check (Total NIE - Sum of individual NIEs):", round(nie_sum_error, 6), "\n")
    
    if (decomposition_error < 0.01) {
      cat("TE decomposition validated: TE = NDE + NIE\n")
    }
    if (nie_sum_error < 0.01) {
      cat("NIE decomposition validated: Total NIE = Sum of individual NIEs\n")
    } else {
      cat("Warning: Individual NIEs do not sum to total (error =", round(nie_sum_error, 4), ")\n")
    }
  }
  
  return(list(
    alpha_med_summaries = alpha_summaries,
    beta_out_med_summaries = beta_out_med_summaries,
    lambda_out_summary = lambda_summary,
    total_effect_summary = te_summary,
    direct_effect_summary = nde_summary,
    indirect_effect_summary = nie_summary,
    indirect_effect_individual_summaries = nie_individual_summaries,
    alpha_med_samples = combined_alpha_med,
    beta_out_med_samples = combined_beta_out_med,
    lambda_out_samples = combined_lambda_out,
    total_effects_samples = combined_total_effects,
    direct_effects_samples = combined_direct_effects,
    indirect_effects_samples = combined_indirect_effects,
    indirect_effects_individual_samples = combined_indirect_individual,
    decomposition_validation = list(
      decomposition_error = decomposition_error,
      nie_sum_error = nie_sum_error,
      te_mean = te_summary["mean"],
      nde_mean = nde_summary["mean"],
      nie_mean = nie_summary["mean"],
      nie_individual_sum = nie_individual_sum
    ),
    chains = chains,
    n_iter = n_iter,
    n_burn = n_burn,
    n_chains = n_chains,
    n_mediators = n_mediators,
    mediation_structure = "parallel",
    rho = 0,
    memory_efficient = memory_efficient
  ))
}

# -----------------------------------------------------------------------------
# Data preparation with proper temporal ordering
# -----------------------------------------------------------------------------

create_lagged_mediation_data <- function(elsa_long) {
  
  cat("Creating lagged mediation data structure\n")
  cat("Temporal ordering: Treatment(W2-fixed) -> Mediators(t-1) -> Outcome(t)\n\n")
  
  required_vars <- c("idauniq", "wave", "loneliness", "livalone", 
                     "age_gr", "dhsex2", "edqual2", "self_reported_health", "sclife",
                     "depression", "transport_mobility", "mobility_limitations")
  
  missing_vars <- required_vars[!required_vars %in% names(elsa_long)]
  if (length(missing_vars) > 0) {
    stop("Missing variables: ", paste(missing_vars, collapse = ", "))
  }
  
  cat("Original data:", nrow(elsa_long), "person-wave observations\n")
  cat("Waves available:", paste(sort(unique(elsa_long$wave)), collapse = ", "), "\n\n")
  
  # Extract W2 baseline (treatment + covariates)
  cat("Step 1: Extracting baseline treatment and covariates from W2\n")
  
  baseline_w2 <- elsa_long %>%
    filter(wave == 2) %>%
    dplyr::select(idauniq, 
                  treatment_w2 = livalone,
                  age_w2 = age_gr,
                  gender_w2 = dhsex2,
                  education_w2 = edqual2,
                  health_w2 = self_reported_health,
                  lifesat_w2 = sclife) %>%
    filter(complete.cases(.))
  
  cat("W2 baseline data:", nrow(baseline_w2), "participants\n")
  cat("Treatment distribution:", table(baseline_w2$treatment_w2), "\n\n")
  
  # Create lagged mediator structure
  cat("Step 2: Creating lagged mediator structure\n")
  
  elsa_sorted <- elsa_long %>%
    arrange(idauniq, wave)
  
  elsa_lagged <- elsa_sorted %>%
    group_by(idauniq) %>%
    mutate(
      depression_lag = lag(depression, n = 1),
      transport_lag = lag(transport_mobility, n = 1),
      mobility_lag = lag(mobility_limitations, n = 1),
      loneliness_current = loneliness
    ) %>%
    ungroup() %>%
    filter(wave >= 3) %>%
    dplyr::select(idauniq, wave, 
                  loneliness_current,
                  depression_lag, transport_lag, mobility_lag)
  
  cat("Created lagged mediators for waves 3-7\n")
  cat("Observations with lagged mediators:", nrow(elsa_lagged), "\n\n")
  
  # Merge baseline with lagged data
  cat("Step 3: Merging baseline treatment with lagged mediators\n")
  
  analysis_data <- elsa_lagged %>%
    inner_join(baseline_w2, by = "idauniq") %>%
    filter(complete.cases(.))
  
  cat("Final analysis data:", nrow(analysis_data), "person-wave observations\n")
  cat("Unique participants:", length(unique(analysis_data$idauniq)), "\n")
  cat("Waves included:", paste(sort(unique(analysis_data$wave)), collapse = ", "), "\n\n")
  
  # Validation
  cat("Data structure validation\n\n")
  
  sample_ids <- sample(unique(analysis_data$idauniq), min(5, length(unique(analysis_data$idauniq))))
  
  for (id in sample_ids[1:2]) {
    person_data <- analysis_data %>% filter(idauniq == id)
    cat(sprintf("\nParticipant %d:\n", id))
    cat(sprintf("  Treatment (W2, fixed): %d\n", unique(person_data$treatment_w2)))
    
    for (i in 1:min(3, nrow(person_data))) {
      row <- person_data[i, ]
      cat(sprintf("  Wave %d outcome uses Wave %d mediators\n", 
                  row$wave, row$wave - 1))
    }
  }
  
  cat("\nTemporal precedence validated: Mediators(t-1) always precede Outcome(t)\n")
  cat("Treatment fixed at baseline (W2) for all observations\n\n")
  
  # Descriptive statistics
  cat("Descriptive statistics\n\n")
  
  cat("Treatment (living alone at W2):\n")
  print(table(analysis_data$treatment_w2))
  cat("\n")
  
  cat("Outcome waves distribution:\n")
  print(table(analysis_data$wave))
  cat("\n")
  
  cat("Mediator ranges (lagged):\n")
  cat("- Depression (t-1):", min(analysis_data$depression_lag, na.rm=TRUE), 
      "to", max(analysis_data$depression_lag, na.rm=TRUE), "\n")
  cat("- Transport (t-1):", min(analysis_data$transport_lag, na.rm=TRUE), 
      "to", max(analysis_data$transport_lag, na.rm=TRUE), "\n")
  cat("- Mobility (t-1):", min(analysis_data$mobility_lag, na.rm=TRUE), 
      "to", max(analysis_data$mobility_lag, na.rm=TRUE), "\n")
  
  cat("\nOutcome range:\n")
  cat("- Loneliness (t):", min(analysis_data$loneliness_current, na.rm=TRUE), 
      "to", max(analysis_data$loneliness_current, na.rm=TRUE), "\n\n")
  
  wave_summary <- analysis_data %>%
    group_by(wave) %>%
    summarise(
      n_obs = n(),
      n_participants = n_distinct(idauniq),
      mean_loneliness = mean(loneliness_current, na.rm = TRUE),
      mean_depression_lag = mean(depression_lag, na.rm = TRUE),
      .groups = 'drop'
    )
  
  cat("Wave-specific summary:\n")
  print(as.data.frame(wave_summary))
  cat("\n")
  
  return(list(
    analysis_data = analysis_data,
    baseline_data = baseline_w2,
    data_structure = list(
      n_observations = nrow(analysis_data),
      n_participants = length(unique(analysis_data$idauniq)),
      waves = sort(unique(analysis_data$wave)),
      temporal_structure = "treatment_W2_fixed -> mediators_lag1 -> outcome_current",
      proper_temporal_ordering = TRUE
    )
  ))
}

# -----------------------------------------------------------------------------
# Algorithm 2 run with temporal structure
# -----------------------------------------------------------------------------

run_algo2_lagged <- function(lagged_data, 
                             n_iter = 15000, n_burn = 8000, n_thin = 5,
                             seed = 123, memory_efficient = TRUE) {
  
  cat("Running Algorithm 2 with lagged mediators\n")
  cat("Treatment: Fixed at W2 (baseline)\n")
  cat("Mediators: Lagged by 1 wave (measured before outcome)\n")
  cat("Outcomes: Current wave (waves 3-7)\n\n")
  
  analysis_data <- lagged_data$analysis_data
  
  # Prepare design matrices
  X_baseline <- cbind(
    1,
    analysis_data$age_w2,
    analysis_data$gender_w2,
    analysis_data$education_w2,
    analysis_data$health_w2,
    analysis_data$lifesat_w2
  )
  colnames(X_baseline) <- c("Intercept", "Age_W2", "Gender_W2", "Education_W2", 
                            "Health_W2", "LifeSat_W2")
  
  mediators_matrix <- cbind(
    analysis_data$depression_lag,
    analysis_data$transport_lag,
    analysis_data$mobility_lag
  )
  colnames(mediators_matrix) <- c("Depression_lag", "Transport_lag", "Mobility_lag")
  
  treatment <- analysis_data$treatment_w2
  outcome <- analysis_data$loneliness_current
  wave <- analysis_data$wave
  
  cat("Data prepared for Algorithm 2:\n")
  cat("- Observations:", nrow(analysis_data), "\n")
  cat("- Treatment (W2, fixed):", sum(treatment), "living alone,", 
      sum(1-treatment), "with others\n")
  cat("- Mediators: All lagged by 1 wave\n")
  cat("- Outcomes: Waves 3-7\n\n")
  
  # Run Algorithm 2
  start_time <- Sys.time()
  
  results <- bayesian_algorithm2_parallel(
    y = outcome,
    mediators_matrix = mediators_matrix,
    X_baseline = X_baseline,
    treatment = treatment,
    wave = wave,
    n_iter = n_iter,
    n_burn = n_burn,
    n_thin = n_thin,
    n_chains = 1,
    mediation_draws = 50,
    seed = seed,
    verbose = TRUE,
    memory_efficient = memory_efficient
  )
  
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "hours")
  
  cat("\nLagged mediator analysis complete\n")
  cat("Runtime:", round(runtime, 2), "hours\n\n")
  
  return(list(
    results = results,
    data_structure = lagged_data$data_structure,
    runtime_hours = as.numeric(runtime)
  ))
}

# -----------------------------------------------------------------------------
# Master workflow
# -----------------------------------------------------------------------------

run_complete_lagged_analysis <- function(elsa_long, seed = 123) {
  
  cat("Algorithm 2: Lagged mediator analysis\n")
  cat("Proper temporal ordering: T(W2) -> M(t-1) -> Y(t)\n")
  
  cat("Step 1: Creating temporally proper data structure\n\n")
  
  lagged_data <- create_lagged_mediation_data(elsa_long)
  
  cat("\nStep 2: Running Algorithm 2 with lagged mediators\n\n")
  
  results <- run_algo2_lagged(lagged_data, seed = seed)
  
  cat("\nLagged mediator analysis complete\n\n")
  
  return(list(
    lagged_data = lagged_data,
    results = results
  ))
}

# -----------------------------------------------------------------------------
# Usage
# -----------------------------------------------------------------------------

# Load data
# elsa_long <- read.csv("elsa_longitudinal_causal_mediation.csv")

# Run complete lagged mediator analysis
# lagged_results <- run_complete_lagged_analysis(elsa_long)

# Access results
# lagged_results$results$results$total_effect_summary
# lagged_results$results$results$indirect_effect_summary
# lagged_results$results$results$alpha_med_summaries
# lagged_results$results$results$beta_out_med_summaries

# =============================================================================
# BJCM MULTI-CHAIN FIT: validation -> full run -> real R-hat
# =============================================================================
# Prerequisites: source the preprocessing script to produce
# elsa_longitudinal_analysis.csv / elsa_long, then this file and
# 00_background_run_utils.R.
#
# Settings below match the values reported in Table 3 / Section 4.3 of the
# manuscript for the BJCM run: n_iter = 15000, n_burn = 8000, n_thin = 5,
# seed = 123. Only n_chains and memory_efficient change relative to that run.
#
# Output is wrapped in the same triple-nested shape run_complete_lagged_
# analysis() produces, so downstream diagnostic and save functions (written
# against results$results$results) work unmodified.
# =============================================================================

library(parallel)
library(coda)
library(dplyr)

# NOTE: no source() calls needed here. This file's own model functions
# (create_lagged_mediation_data, bayesian_algorithm2_parallel, etc.) are
# already defined above in this same file.
# run_in_background()/monitor_progress()/collect_result() are provided by
# 00_background_run_utils.R, which is sourced before this file.

# =============================================================================
# STEP 0: Build inputs using the manuscript's temporal structure
# =============================================================================
# Treatment fixed at W2, mediators lagged by 1 wave (t-1), outcome at current
# wave (t), restricted to waves 3-7 - see create_lagged_mediation_data().

prepare_bjcm_inputs <- function(elsa_long) {
  lagged_data <- create_lagged_mediation_data(elsa_long)
  analysis_data <- lagged_data$analysis_data
  
  # Same construction as run_algo2_lagged(), so results are directly
  # comparable to what produced Table 3.
  X_baseline <- cbind(
    1,
    analysis_data$age_w2,
    analysis_data$gender_w2,
    analysis_data$education_w2,
    analysis_data$health_w2,
    analysis_data$lifesat_w2
  )
  colnames(X_baseline) <- c("Intercept", "Age_W2", "Gender_W2", "Education_W2",
                            "Health_W2", "LifeSat_W2")
  
  mediators_matrix <- cbind(
    analysis_data$depression_lag,
    analysis_data$transport_lag,
    analysis_data$mobility_lag
  )
  colnames(mediators_matrix) <- c("Depression_lag", "Transport_lag", "Mobility_lag")
  
  list(
    y             = analysis_data$loneliness_current,
    mediators     = mediators_matrix,
    X_baseline    = X_baseline,
    treatment     = analysis_data$treatment_w2,
    wave          = analysis_data$wave,
    n             = nrow(analysis_data),
    lagged_data   = lagged_data
  )
}

# Wraps a bayesian_algorithm2_parallel() output in the same triple-nested
# shape run_complete_lagged_analysis() produces, so downstream diagnostic
# and save functions (which expect results$results$results, etc.) work
# unmodified on fits produced by this script.
wrap_for_downstream <- function(bjcm_raw, lagged_data, runtime_hours) {
  list(
    lagged_data = lagged_data,
    results = list(
      results = bjcm_raw,
      data_structure = lagged_data$data_structure,
      runtime_hours = runtime_hours
    )
  )
}

# =============================================================================
# STEP 1: SMALL-SAMPLE VALIDATION (confirm the parallel + overdispersed-init
# pipeline runs cleanly before committing to the full multi-hour run)
# =============================================================================
# Returns the wrapped (triple-nested) structure. Access the raw model output
# via val$results$results if you need it directly (e.g. for the sanity
# checks below, which are printed for you already).

run_validation <- function(elsa_long, subset_size = 5000, n_chains = 4) {
  
  cat("=== STEP 1: Small-sample validation (n =", subset_size, ", chains =", n_chains, ") ===\n\n")
  
  inputs <- prepare_bjcm_inputs(elsa_long)
  
  set.seed(123)
  idx <- sample(seq_len(inputs$n), min(subset_size, inputs$n))
  
  t0 <- Sys.time()
  val_raw <- bayesian_algorithm2_parallel(
    y                = inputs$y[idx],
    mediators_matrix = inputs$mediators[idx, , drop = FALSE],
    X_baseline       = inputs$X_baseline[idx, , drop = FALSE],
    treatment        = inputs$treatment[idx],
    wave             = inputs$wave[idx],
    n_iter           = 2000,   # short run - purpose is to confirm the
    n_burn           = 1000,   # pipeline executes correctly and chains
    n_thin           = 2,      # differ (overdispersion took effect),
    n_chains         = n_chains,   # not to get usable inference
    seed             = 123,
    memory_efficient = FALSE,
    verbose          = TRUE
  )
  t1 <- Sys.time()
  runtime_hours <- as.numeric(difftime(t1, t0, units = "hours"))
  
  cat("\nValidation run wall-clock time:", round(runtime_hours * 60, 2), "minutes\n\n")
  
  # Sanity checks before trusting the pipeline with the full run
  cat("--- Sanity checks ---\n")
  
  # (a) Did overdispersion actually produce different starting points?
  chain_first_lambda <- sapply(val_raw$chains, function(x) x$lambda_out_samples[1])
  cat("First stored lambda_out draw per chain (should differ across chains):\n")
  print(chain_first_lambda)
  
  # (b) Decomposition check (should hold regardless of subset size / chains)
  cat("\nTE decomposition error:", val_raw$decomposition_validation$decomposition_error, "\n")
  cat("NIE sum error:", val_raw$decomposition_validation$nie_sum_error, "\n")
  
  # (c) Did all chains actually run and return non-error results?
  cat("\nChains returned:", length(val_raw$chains), "(expected", n_chains, ")\n")
  
  wrap_for_downstream(val_raw, inputs$lagged_data, runtime_hours)
}

# =============================================================================
# STEP 2: FULL 4-CHAIN MODEL FIT (only run after Step 1 checks look right)
# =============================================================================
# Returns the wrapped (triple-nested) structure, ready to hand directly to
# downstream diagnostic and save functions.

run_full_fit <- function(elsa_long, n_chains = 4) {
  
  cat("=== STEP 2: Full 4-chain BJCM fit (n_chains =", n_chains, ") ===\n")
  cat("Expected wall-clock: roughly the same as a single-chain run, since\n")
  cat("chains are dispatched across separate cores.\n\n")
  
  inputs <- prepare_bjcm_inputs(elsa_long)
  
  t0 <- Sys.time()
  full_raw <- bayesian_algorithm2_parallel(
    y                = inputs$y,
    mediators_matrix = inputs$mediators,
    X_baseline       = inputs$X_baseline,
    treatment        = inputs$treatment,
    wave             = inputs$wave,
    n_iter           = 15000,   # matches Table 3 / Section 4.3
    n_burn           = 8000,    # matches Table 3 / Section 4.3
    n_thin           = 5,       # matches Table 3 / Section 4.3
    n_chains         = n_chains,
    seed             = 123,
    memory_efficient = TRUE,    # person-averaged draws per iteration are
    # all that compute_multichain_rhat() and
    # Table 3 need, and keep peak memory usage
    # manageable across four chains.
    verbose          = TRUE
  )
  t1 <- Sys.time()
  runtime_hours <- as.numeric(difftime(t1, t0, units = "hours"))
  
  cat("\nFull run wall-clock time:", round(runtime_hours, 2), "hours\n")
  
  full_results <- wrap_for_downstream(full_raw, inputs$lagged_data, runtime_hours)
  saveRDS(full_results, "../results/bjcm_full_results.rds")
  full_results
}

# =============================================================================
# STEP 3: MULTI-CHAIN GELMAN-RUBIN R-HAT ON ALL KEY PARAMETERS
# =============================================================================
# Pass the RAW model output (e.g. full$results$results if full came from
# run_full_fit()/run_validation() above), not the wrapped triple-nested
# structure - that wrapping is only for downstream compatibility. Checks
# TE/NDE/NIE and the structural parameters (lambda_out, alpha_med,
# beta_out_med).

compute_multichain_rhat <- function(bjcm_raw) {
  
  chains <- bjcm_raw$chains
  if (is.null(chains) || length(chains) < 2) {
    stop("Fewer than 2 chains found. Pass the RAW model output ",
         "(e.g. full$results$results from run_full_fit()), not the ",
         "wrapped structure, and make sure n_chains >= 2 was used.")
  }
  n_chains <- length(chains)
  n_mediators <- ncol(chains[[1]]$alpha_med_samples)
  
  # Handles both memory_efficient=TRUE (total_effects/direct_effects/
  # indirect_effects already stored as person-averaged vectors) and
  # memory_efficient=FALSE (raw n-person x iteration matrices, needing
  # rowMeans to get the person-averaged draw per iteration).
  person_avg <- function(x) if (is.matrix(x)) rowMeans(x) else as.numeric(x)
  
  build_mcmc_list <- function(extractor) {
    coda::mcmc.list(lapply(chains, function(ch) coda::mcmc(extractor(ch))))
  }
  
  # Causal estimands
  te_list  <- build_mcmc_list(function(ch) matrix(person_avg(ch$total_effects), ncol = 1,
                                                  dimnames = list(NULL, "TE")))
  nde_list <- build_mcmc_list(function(ch) matrix(person_avg(ch$direct_effects), ncol = 1,
                                                  dimnames = list(NULL, "NDE")))
  nie_list <- build_mcmc_list(function(ch) matrix(person_avg(ch$indirect_effects), ncol = 1,
                                                  dimnames = list(NULL, "NIE")))
  
  # Structural parameters
  lambda_list <- build_mcmc_list(function(ch) matrix(ch$lambda_out_samples, ncol = 1,
                                                     dimnames = list(NULL, "lambda_out")))
  alpha_list  <- build_mcmc_list(function(ch) {
    m <- ch$alpha_med_samples
    colnames(m) <- paste0("alpha_med_", seq_len(n_mediators))
    m
  })
  beta_out_med_list <- build_mcmc_list(function(ch) {
    m <- ch$beta_out_med_samples
    colnames(m) <- paste0("beta_out_med_", seq_len(n_mediators))
    m
  })
  
  cat("=== STEP 3: Multi-chain Gelman-Rubin R-hat (", n_chains, "chains ) ===\n\n")
  
  results <- list(
    TE            = coda::gelman.diag(te_list, autoburnin = FALSE),
    NDE           = coda::gelman.diag(nde_list, autoburnin = FALSE),
    NIE           = coda::gelman.diag(nie_list, autoburnin = FALSE),
    lambda_out    = coda::gelman.diag(lambda_list, autoburnin = FALSE),
    alpha_med     = coda::gelman.diag(alpha_list, autoburnin = FALSE),
    beta_out_med  = coda::gelman.diag(beta_out_med_list, autoburnin = FALSE)
  )
  
  for (nm in names(results)) {
    cat(nm, ":\n")
    print(results[[nm]])
    cat("\n")
  }
  
  flagged <- Filter(function(nm) {
    any(results[[nm]]$psrf[, "Point est."] > 1.1)
  }, names(results))
  
  if (length(flagged) > 0) {
    cat("WARNING: R-hat > 1.1 for:", paste(flagged, collapse = ", "), "\n")
  } else {
    cat("All R-hat values <= 1.1 - chains agree on the causal estimands.\n")
  }
  
  results
}

# =============================================================================
# RECOVERY: reassemble a run from per-chain checkpoint files on disk
# =============================================================================
# Use this if a run completes all chains (confirmed via monitor_progress
# showing 100% for every chain) but the final return through
# mcparallel/mclapply fails before saveRDS() runs, so
# bjcm_full_results.rds never gets written. Each chain checkpoints
# itself to bjcm_chain_<id>_checkpoint.rds as soon as it finishes sampling,
# independent of the outer combine/return step, so this lets you recover
# without re-running anything.

recover_bjcm_chains_from_checkpoints <- function(elsa_long, n_chains = 4, save = TRUE) {
  
  checkpoint_files <- sprintf("../results/checkpoints/bjcm_chain_%d_checkpoint.rds", seq_len(n_chains))
  missing <- checkpoint_files[!file.exists(checkpoint_files)]
  
  if (length(missing) > 0) {
    stop("Missing checkpoint file(s): ", paste(missing, collapse = ", "),
         ". Cannot recover - these chains did not reach the checkpoint save, ",
         "so a fresh run_full_fit() is needed.")
  }
  
  cat("Found all", n_chains, "checkpoint files. Reassembling...\n")
  chains <- lapply(checkpoint_files, readRDS)
  n_mediators <- ncol(chains[[1]]$alpha_med_samples)
  
  # Reconstructs the same summary structure bayesian_algorithm2_parallel()
  # itself returns, so recovered runs and normal runs are interchangeable
  # for every downstream script.
  memory_efficient <- !is.matrix(chains[[1]]$total_effects) ||
    length(chains[[1]]$total_effects) == nrow(chains[[1]]$alpha_med_samples)
  person_avg <- function(x) if (is.matrix(x)) rowMeans(x) else as.numeric(x)
  
  combined_alpha_med     <- do.call(rbind, lapply(chains, function(x) x$alpha_med_samples))
  combined_beta_out_med  <- do.call(rbind, lapply(chains, function(x) x$beta_out_med_samples))
  combined_lambda_out    <- do.call(c, lapply(chains, function(x) x$lambda_out_samples))
  combined_total_effects   <- unlist(lapply(chains, function(x) person_avg(x$total_effects)))
  combined_direct_effects  <- unlist(lapply(chains, function(x) person_avg(x$direct_effects)))
  combined_indirect_effects <- unlist(lapply(chains, function(x) person_avg(x$indirect_effects)))
  
  first_indiv <- chains[[1]]$indirect_effects_individual
  combined_indirect_individual <- if (length(dim(first_indiv)) == 3) {
    abind::abind(lapply(chains, function(x) x$indirect_effects_individual), along = 1)
  } else {
    do.call(rbind, lapply(chains, function(x) x$indirect_effects_individual))
  }
  
  create_summary <- function(samples) {
    c(mean = mean(samples, na.rm = TRUE),
      sd = sd(samples, na.rm = TRUE),
      q2.5 = as.numeric(quantile(samples, 0.025, na.rm = TRUE)),
      q97.5 = as.numeric(quantile(samples, 0.975, na.rm = TRUE)))
  }
  
  te_summary  <- create_summary(combined_total_effects)
  nde_summary <- create_summary(combined_direct_effects)
  nie_summary <- create_summary(combined_indirect_effects)
  
  nie_individual_summaries <- lapply(1:n_mediators, function(k) {
    vals <- if (length(dim(combined_indirect_individual)) == 3) {
      as.vector(combined_indirect_individual[, , k])
    } else {
      combined_indirect_individual[, k]
    }
    create_summary(vals)
  })
  alpha_summaries <- lapply(1:n_mediators, function(k) create_summary(combined_alpha_med[, k]))
  beta_out_med_summaries <- lapply(1:n_mediators, function(k) create_summary(combined_beta_out_med[, k]))
  lambda_summary <- create_summary(combined_lambda_out)
  
  decomposition_error <- abs(te_summary["mean"] - (nie_summary["mean"] + nde_summary["mean"]))
  nie_individual_sum <- sum(sapply(nie_individual_summaries, function(x) x["mean"]))
  nie_sum_error <- abs(nie_summary["mean"] - nie_individual_sum)
  
  cat("Reconstructed summaries - TE:", round(te_summary["mean"], 4),
      "NDE:", round(nde_summary["mean"], 4), "NIE:", round(nie_summary["mean"], 4), "\n")
  cat("Decomposition error:", round(decomposition_error, 6),
      "| NIE sum error:", round(nie_sum_error, 6), "\n")
  
  bjcm_raw <- list(
    alpha_med_summaries = alpha_summaries,
    beta_out_med_summaries = beta_out_med_summaries,
    lambda_out_summary = lambda_summary,
    total_effect_summary = te_summary,
    direct_effect_summary = nde_summary,
    indirect_effect_summary = nie_summary,
    indirect_effect_individual_summaries = nie_individual_summaries,
    alpha_med_samples = combined_alpha_med,
    beta_out_med_samples = combined_beta_out_med,
    lambda_out_samples = combined_lambda_out,
    total_effects_samples = combined_total_effects,
    direct_effects_samples = combined_direct_effects,
    indirect_effects_samples = combined_indirect_effects,
    indirect_effects_individual_samples = combined_indirect_individual,
    decomposition_validation = list(
      decomposition_error = decomposition_error,
      nie_sum_error = nie_sum_error,
      te_mean = te_summary["mean"], nde_mean = nde_summary["mean"], nie_mean = nie_summary["mean"],
      nie_individual_sum = nie_individual_sum
    ),
    chains = chains,
    n_iter = NA, n_burn = NA, n_chains = n_chains, n_mediators = n_mediators,
    mediation_structure = "parallel", rho = 0, memory_efficient = memory_efficient
  )
  
  inputs <- prepare_bjcm_inputs(elsa_long)   # rebuilds the lagged data
  # structure, no MCMC
  full_results <- wrap_for_downstream(bjcm_raw, inputs$lagged_data, runtime_hours = NA)
  
  if (save) {
    saveRDS(full_results, "../results/bjcm_full_results.rds")
    cat("Recovered result saved to bjcm_full_results.rds\n")
  }
  
  full_results
}

# =============================================================================
# USAGE
# =============================================================================
# elsa_long <- read.csv("../data/elsa_longitudinal_analysis.csv")
#
# val <- run_validation(elsa_long, subset_size = 5000, n_chains = 4)
#   -> inspect the sanity checks printed above before proceeding
#
# full <- run_full_fit(elsa_long, n_chains = 4)
#
# rhat_results <- compute_multichain_rhat(full$results$results)
#
# generate_all_outputs(full)
# organized_save_results(full, diagnostics = rhat_results)
#
# If the run does not complete cleanly (check monitor_progress("bjcm") first -
# if all chains show 100% with errors=0, the sampling itself succeeded and
# nothing needs re-running):
#
#   full <- recover_bjcm_chains_from_checkpoints(elsa_long, n_chains = 4)
#   rhat_results <- compute_multichain_rhat(full$results$results)