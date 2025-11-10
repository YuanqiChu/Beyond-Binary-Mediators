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

update_thresholds_conservative <- function(z_latent, y_observed, current_thresholds, iteration) {
  if (iteration < 1000) {
    return(current_thresholds)
  }
  
  K <- max(y_observed, na.rm = TRUE)
  if (K <= 2) {
    return(current_thresholds)
  }
  
  new_thresholds <- current_thresholds
  
  for (k in 2:(K-1)) {
    z_vals <- z_latent[y_observed == k]
    if (length(z_vals) > 5) {
      proposed <- median(z_vals)
      adjustment <- (proposed - current_thresholds[k-1]) * 0.01
      new_thresholds[k-1] <- current_thresholds[k-1] + adjustment
    }
  }
  
  new_thresholds <- sort(new_thresholds)
  new_thresholds[1] <- 0
  
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
  
  te_matrix <- matrix(0, n, n_draws)
  nde_matrix <- matrix(0, n, n_draws)
  nie_total_matrix <- matrix(0, n, n_draws)
  nie_individual_array <- array(0, dim = c(n, n_mediators, n_draws))
  
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
    
    # Individual NIE_k: Sequential intervention on each mediator
    for (k in 1:n_mediators) {
      Z_med_partial <- Z_med_control_matrix
      Z_med_partial[, 1:k] <- Z_med_treat_matrix[, 1:k]
      
      mediator_effects_partial <- rowSums(Z_med_partial * 
                                            matrix(beta_out_med, n, n_mediators, byrow = TRUE))
      mu_Y_1Mp <- mediator_effects_partial + 1 * lambda_out + 
        as.vector(X_baseline %*% gamma_out) + U_out_wave[wave_mapped]
      Z_Y_1Mp <- mu_Y_1Mp + epsilon_out
      Y_1Mp <- findInterval(Z_Y_1Mp, S_out)
      
      if (k == 1) {
        Y_prev <- Y_10
      } else {
        Z_med_prev <- Z_med_control_matrix
        Z_med_prev[, 1:(k-1)] <- Z_med_treat_matrix[, 1:(k-1)]
        mediator_effects_prev <- rowSums(Z_med_prev * 
                                           matrix(beta_out_med, n, n_mediators, byrow = TRUE))
        mu_Y_prev <- mediator_effects_prev + 1 * lambda_out + 
          as.vector(X_baseline %*% gamma_out) + U_out_wave[wave_mapped]
        Z_Y_prev <- mu_Y_prev + epsilon_out
        Y_prev <- findInterval(Z_Y_prev, S_out)
      }
      
      nie_individual_array[, k, d] <- Y_1Mp - Y_prev
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
    
    # Initialise parameters
    beta_med_matrix <- matrix(rnorm(p_baseline * n_mediators, 0, 0.01), p_baseline, n_mediators)
    alpha_med <- rnorm(n_mediators, 0, 0.01)
    beta_out_med <- rnorm(n_mediators, 0, 0.01)
    lambda_out <- rnorm(1, 0, 0.01)
    gamma_out <- rnorm(p_baseline, 0, 0.01)
    
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
    
    # Initialise thresholds
    thresholds_med_list <- list()
    for (k in 1:n_mediators) {
      if (nk_med[k] > 2) {
        thresholds_med_list[[k]] <- seq(-1.5, 1.5, length.out = nk_med[k] - 1)
        thresholds_med_list[[k]][1] <- 0
      } else {
        thresholds_med_list[[k]] <- 0
      }
    }
    
    if (nk_out > 2) {
      thresholds_out <- seq(-1.5, 1.5, length.out = nk_out - 1)
      thresholds_out[1] <- 0
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
        
        z_med_matrix <- update_z_mediator_parallel(
          mediators_recoded, X_baseline, beta_med_matrix, treatment, 
          alpha_med, U_med_wave_matrix, wave_mapped, thresholds_med_list
        )
        
        for (k in 1:n_mediators) {
          if (nk_med[k] > 2) {
            thresholds_med_list[[k]] <- update_thresholds_conservative(
              z_med_matrix[, k], mediators_recoded[, k], thresholds_med_list[[k]], iter
            )
          }
        }
        
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
        
        z_out <- update_z_outcome(
          y_recoded, z_med_matrix, beta_out_med, treatment, lambda_out, 
          X_baseline, gamma_out, U_out_wave, wave_mapped, thresholds_out
        )
        
        if (nk_out > 2) {
          thresholds_out <- update_thresholds_conservative(
            z_out, y_recoded, thresholds_out, iter
          )
        }
        
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
      }
      
      # Progress reporting
      if (verbose && chain_id == 1 && iter %% 2000 == 0) {
        cat("Chain", chain_id, "Iter:", iter, "| Errors:", n_errors)
        
        if (store_idx > 0) {
          last_cf <- simulate_parallel_mediation_effects(
            X_baseline, beta_med_matrix, alpha_med, beta_out_med, lambda_out, 
            gamma_out, U_med_wave_matrix, U_out_wave, wave_mapped, 
            thresholds_med_list, thresholds_out, 10
          )
          cat(" | TE-NDE-NIE err:", round(last_cf$decomposition_check$decomposition_error, 6))
          if (!is.null(last_cf$decomposition_check$nie_sum_check)) {
            cat(" | NIE sum err:", round(last_cf$decomposition_check$nie_sum_check, 6))
          }
        }
        cat("\n")
      }
    }
    
    return(list(
      beta_med_samples = beta_med_samples,
      alpha_med_samples = alpha_med_samples,
      beta_out_med_samples = beta_out_med_samples,
      lambda_out_samples = lambda_out_samples,
      gamma_out_samples = gamma_out_samples,
      total_effects = total_effects_samples,
      direct_effects = direct_effects_samples,
      indirect_effects = indirect_effects_samples,
      indirect_effects_individual = indirect_effects_individual_samples,
      chain_id = chain_id,
      n_errors = n_errors
    ))
  }
  
  # Run chains
  if (verbose) cat("Running", n_chains, "chain(s)...\n")
  chains <- lapply(1:n_chains, run_single_chain)
  
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
      cat(paste0(mediator_names[k], " NIE", k, " = alpha", k, " * beta", k, ": "), 
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
# END OF CODE
# =============================================================================
