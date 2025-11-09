# =============================================================================
# BMEOP Model: Usage Guide for ELSA Loneliness Analysis
# =============================================================================
#
# Files included:
# 1. bmeop_model.R - Core BMEOP algorithms
# 2. bmeop_diagnost.R - Convergence diagnostics and post-processing
# 3. bmeop_usage_guide.R - Usage examples
# =============================================================================

# =============================================================================
# SETUP
# =============================================================================

# Load the main model code
source("bmeop_model.R")

# Load diagnostic functions
source("bmeop_diagnost.R")

# Required packages (already loaded by sourcing above)
# library(MASS)
# library(mvtnorm)
# library(coda)
# library(Matrix)
# library(dplyr)
# library(corrplot)
# library(car)

# =============================================================================
# EXAMPLE 1: SMALL-SCALE TESTING
# =============================================================================

# Test the model on a subset before running on full data
# Approximate runtime: 2-3 hours for n=2,000

# Load ELSA data (adjust path as needed)
# elsa_long <- read.csv("path/to/elsa_long.csv")

# Run on subset
test_results <- run_elsa_bmeop(
  elsa_long = elsa_long,
  subset_size = 2000,           # Test on 2,000 observations
  seed = 123,
  ate_draws = 50,               # Number of Monte Carlo draws
  ate_method = "analytical",    # "analytical" or "multidraw"
  verbose = TRUE
)

# Quick assessment
print(test_results$beta_summary)
print(test_results$ate_stability)

# =============================================================================
# EXAMPLE 2: FULL ELSA ANALYSIS
# =============================================================================

# Run on complete dataset
# Approximate runtime: 15-18 hours for n=42,185
# Hardware: Apple M4 processor, 16GB RAM

full_results <- run_full_elsa_analysis(
  elsa_long = elsa_long,
  seed = 123
)

# Save results
saveRDS(full_results, "elsa_bmeop_results.rds")

# =============================================================================
# EXAMPLE 3: COMPREHENSIVE DIAGNOSTICS
# =============================================================================

# Load saved results if starting from here
# full_results <- readRDS("elsa_bmeop_results.rds")

# Extract design matrix and outcome for VIF analysis
required_vars <- c("loneliness", "livalone", "age_gr", "dhsex2", 
                   "edqual2", "self_reported_health", "sclife", 
                   "depression", "transport_mobility", "mobility_limitations")
analysis_data <- elsa_long[complete.cases(elsa_long[, required_vars]), ]

# Create design matrix (same as in model fitting)
X <- cbind(
  intercept = 1,
  livalone_TREATMENT = analysis_data$livalone,
  age_gr = analysis_data$age_gr,
  dhsex2 = analysis_data$dhsex2,
  edqual2 = analysis_data$edqual2,
  health = analysis_data$self_reported_health,
  sclife = analysis_data$sclife,
  depression = analysis_data$depression,
  transport_mobility = analysis_data$transport_mobility,
  mobility_limitations = analysis_data$mobility_limitations
)

y <- analysis_data$loneliness

# Run full diagnostic workflow
diagnostic_results <- run_full_diagnostics(
  results = full_results,
  X = X,
  y = y
)

# =============================================================================
# EXAMPLE 4: CUSTOM ANALYSIS PIPELINE
# =============================================================================

# For maximum control, use the core function directly

# Prepare data
complete_data <- elsa_long[complete.cases(elsa_long[, required_vars]), ]

# Create design matrix
X_custom <- cbind(
  intercept = 1,
  treatment = complete_data$livalone,
  # Add your covariates here
  age = complete_data$age_gr,
  female = complete_data$dhsex2,
  education = complete_data$edqual2
)

# Run BMEOP with custom settings
custom_results <- bayesian_bmeop_elsa(
  y = complete_data$loneliness,
  X = X_custom,
  wave = complete_data$wave,
  treatment_col = 2,              # Column index of treatment
  n_iter = 15000,
  n_burn = 10000,
  n_thin = 5,
  n_chains = 1,                   # Use 1 for efficiency
  prior_beta_var = 4.0,           # Prior variance
  ate_draws = 50,
  ate_method = "analytical",      # Or "multidraw"
  seed = 123,
  verbose = TRUE
)

# =============================================================================
# EXAMPLE 5: EXAMINING SPECIFIC RESULTS
# =============================================================================

# Extract treatment effect estimates
treatment_samples <- full_results$beta_samples[, "livalone_TREATMENT"]

# Summary statistics
cat("Treatment Effect (latent scale):\n")
cat("  Mean:", mean(treatment_samples), "\n")
cat("  SD:", sd(treatment_samples), "\n")
cat("  95% CrI:", quantile(treatment_samples, c(0.025, 0.975)), "\n")
cat("  P(beta > 0):", mean(treatment_samples > 0), "\n\n")

# ATE summary
if (!is.null(full_results$ate_samples)) {
  ate_samples <- full_results$ate_samples
  
  cat("Average Treatment Effect (ordinal scale):\n")
  cat("  Mean:", mean(ate_samples), "categories\n")
  cat("  SD:", sd(ate_samples), "\n")
  cat("  95% CrI:", quantile(ate_samples, c(0.025, 0.975)), "\n")
  cat("  P(ATE > 0):", mean(ate_samples > 0), "\n\n")
}

# Mediator effects
mediator_names <- c("depression", "transport_mobility", "mobility_limitations")

cat("Mediator Effects:\n")
for (med in mediator_names) {
  if (med %in% colnames(full_results$beta_samples)) {
    med_samples <- full_results$beta_samples[, med]
    cat(sprintf("  %s: %.4f [%.4f, %.4f], P(>0) = %.3f\n",
                med,
                mean(med_samples),
                quantile(med_samples, 0.025),
                quantile(med_samples, 0.975),
                mean(med_samples > 0)))
  }
}

# =============================================================================
# EXAMPLE 6: CONVERGENCE ASSESSMENT
# =============================================================================

# Quick convergence check
convergence_diag <- assess_convergence(full_results)

# Visual inspection
visualise_traces(full_results, 
                 params = c("livalone_TREATMENT", "depression", 
                            "transport_mobility", "mobility_limitations"))

# Chain stability across segments
assess_chain_stability(full_results, 
                       n_segments = 4,
                       param = "livalone_TREATMENT")

# Posterior correlations
cor_matrix <- examine_posterior_correlation(
  full_results,
  params = c("intercept", "livalone_TREATMENT", "sclife", "depression"),
  visualise = TRUE
)

# =============================================================================
# EXAMPLE 7: COMPARISON OF ATE METHODS
# =============================================================================

# Compare analytical vs multi-draw methods on a subset

# Analytical method (default, faster)
results_analytical <- run_elsa_bmeop(
  elsa_long = elsa_long,
  subset_size = 1000,
  ate_method = "analytical",
  ate_draws = 50,
  seed = 123
)

# Multi-draw method (slower but potentially more stable for large K)
results_multidraw <- run_elsa_bmeop(
  elsa_long = elsa_long,
  subset_size = 1000,
  ate_method = "multidraw",
  ate_draws = 100,
  seed = 123
)

# Compare estimates
cat("ATE Comparison:\n")
cat("  Analytical: ", mean(results_analytical$ate_samples), 
    " (SD:", sd(results_analytical$ate_samples), ")\n")
cat("  Multi-draw: ", mean(results_multidraw$ate_samples), 
    " (SD:", sd(results_multidraw$ate_samples), ")\n")

# Compare stability
cat("\nStability Ratios:\n")
cat("  Analytical:", results_analytical$ate_stability$stability_ratio, "\n")
cat("  Multi-draw:", results_multidraw$ate_stability$stability_ratio, "\n")

# =============================================================================
# TROUBLESHOOTING
# =============================================================================

# If encountering memory issues:
# - Increase thinning: n_thin = 10 or 20
# - Reduce iterations: n_iter = 10000, n_burn = 5000
# - Process in batches (subset data by wave)

# If encountering convergence warnings:
# - Check trace plots for stationarity
# - Assess stability across chain segments
# - Low ESS for individual parameters may be acceptable if
#   treatment effect and ATE are stable

# If runtime is prohibitive:
# - Use smaller subset for preliminary analysis
# - Consider parallel computation (modify code for multiple chains)
# - Reduce ate_draws to 20-30 for analytical method

# =============================================================================
# END OF USAGE GUIDE
# =============================================================================