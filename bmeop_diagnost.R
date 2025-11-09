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
# END OF DIAGNOSTICS CODE
# =============================================================================