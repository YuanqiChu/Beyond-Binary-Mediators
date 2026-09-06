# =============================================================================
# Extend the proposal_sd grid check with smaller candidates
# =============================================================================
# Follow-up to run_reduced_check_bmeop.R. That grid (0.08 down to 0.02) showed
# the MEAN acceptance rate crossing into the 20%-50% target range around
# proposal_sd=0.020-0.025, but the single HARDEST threshold (j2, historically
# the lowest-acceptance one at full scale) was still predicted to fall short
# of 20% even at proposal_sd=0.020. This script tests smaller candidates and
# re-reports both the mean AND the j2-specific predicted full-scale rate.
#
# Scaling factor: validated against the actual full run. proposal_sd=0.08 at
# n=5000 gave mean acceptance 28.89%; the actual n=50,790 full run gave 9.45%.
# Ratio = 28.89/9.45 = 3.06, matching the theoretical sqrt(50790/5000) = 3.19
# closely enough to trust this scaling for prediction.
#
# Does NOT re-run the six candidates already tested -- loads and merges with
# mh_bmeop_proposal_sd_grid_result.rds from run_reduced_check_bmeop.R.
# =============================================================================

library(dplyr)
library(MASS)
library(mvtnorm)
library(coda)
library(Matrix)
library(parallel)

source("00_background_run_utils.R")
source("02_bmeop_model.R")   # sources threshold_update_mh.R internally

elsa_long <- read.csv("../data/elsa_longitudinal_analysis.csv")

required_vars <- c("idauniq", "wave", "loneliness", "livalone", "age_gr", "dhsex2",
                   "edqual2", "self_reported_health", "sclife", "depression",
                   "transport_mobility", "mobility_limitations")

build_bmeop_X <- function(d) {
  X_components <- list(
    intercept = rep(1, nrow(d)),
    livalone_TREATMENT = d$livalone
  )
  if ("age_gr" %in% names(d)) X_components$age_gr <- d$age_gr
  if ("dhsex2" %in% names(d)) X_components$dhsex2 <- d$dhsex2
  if ("edqual2" %in% names(d)) X_components$edqual2 <- d$edqual2
  if ("self_reported_health" %in% names(d)) X_components$health <- d$self_reported_health
  if ("sclife" %in% names(d)) X_components$sclife <- d$sclife
  if ("depression" %in% names(d)) X_components$depression <- d$depression
  if ("transport_mobility" %in% names(d)) X_components$transport_mobility <- d$transport_mobility
  if ("mobility_limitations" %in% names(d)) X_components$mobility_limitations <- d$mobility_limitations
  X <- do.call(cbind, X_components)
  colnames(X) <- names(X_components)
  X
}

bmeop_complete <- elsa_long[complete.cases(elsa_long[, required_vars]), ]

run_bmeop_validation <- function(complete_data, subset_size = 5000, n_chains = 4,
                                 proposal_sd = 0.08) {
  set.seed(123)
  idx <- sample(seq_len(nrow(complete_data)), min(subset_size, nrow(complete_data)))
  d <- complete_data[idx, ]
  X <- build_bmeop_X(d)
  bayesian_bmeop_elsa(
    y = d$loneliness, X = X, wave = d$wave, treatment_col = 2,
    n_iter = 2000, n_burn = 1000, n_thin = 2,
    n_chains = n_chains, proposal_sd = proposal_sd, seed = 123, verbose = TRUE
  )
}

# ---- Load the existing grid (0.08, 0.05, 0.04, 0.03, 0.025, 0.02) ----------
grid_results <- readRDS("../results/mh_bmeop_proposal_sd_grid_result.rds")
cat("Loaded existing grid with", length(grid_results), "candidate values:",
    paste(names(grid_results), collapse = ", "), "\n\n")

# ---- Empirically validated scaling factor (n=50,790 / n=5,000) ------------
SCALE_FACTOR <- 3.06

# ---- Test smaller candidates -------------------------------------------
new_grid <- c(0.015, 0.01)

for (psd in new_grid) {
  
  psd_str <- as.character(psd)
  if (psd_str %in% names(grid_results)) {
    cat("proposal_sd =", psd, "already tested, skipping.\n\n")
    next
  }
  
  reset_mh_threshold_log()
  
  cat(sprintf("=== Running BMEOP reduced-scale check (n=5000, 4 chains, 2000 iterations, proposal_sd=%.3f) ===\n\n", psd))
  t0 <- Sys.time()
  fit_psd <- run_bmeop_validation(bmeop_complete, subset_size = 5000, n_chains = 4,
                                  proposal_sd = psd)
  t1 <- Sys.time()
  elapsed_mins <- as.numeric(difftime(t1, t0, units = "mins"))
  cat("\nRun time:", round(elapsed_mins, 2), "minutes\n\n")
  
  acc <- get_mh_acceptance_summary()
  cat("Acceptance rates:\n")
  print(acc)
  cat("\n")
  
  grid_results[[psd_str]] <- list(fit = fit_psd, acceptance = acc,
                                  elapsed_mins = elapsed_mins)
}

saveRDS(grid_results, "../results/mh_bmeop_proposal_sd_grid_result.rds")
cat("Updated: ../results/mh_bmeop_proposal_sd_grid_result.rds (now", length(grid_results), "candidates)\n\n")

# =============================================================================
# Extended summary: mean rate AND j2-specific rate (the hardest threshold),
# both raw (n=5000) and scaled prediction for the full n=50,790 run.
# =============================================================================

cat("\n\n=== Extended summary: mean and j2-specific acceptance, raw (n=5000) and predicted (n=50,790) ===\n\n")

summary_tbl <- do.call(rbind, lapply(names(grid_results), function(psd_str) {
  acc <- grid_results[[psd_str]]$acceptance
  mean_rate_n5000 <- mean(acc$rate)
  j2_rate_n5000    <- acc$rate[acc$threshold == "bmeop_outcome_j2"]
  data.frame(
    proposal_sd           = as.numeric(psd_str),
    mean_rate_n5000        = mean_rate_n5000,
    mean_rate_predicted    = mean_rate_n5000 / SCALE_FACTOR,
    j2_rate_n5000          = j2_rate_n5000,
    j2_rate_predicted      = j2_rate_n5000 / SCALE_FACTOR
  )
}))
summary_tbl <- summary_tbl[order(-summary_tbl$proposal_sd), ]
summary_tbl[ , -1] <- round(summary_tbl[ , -1], 4)
print(summary_tbl, row.names = FALSE)

cat("\nTarget for BOTH columns: 0.20-0.50 (mean_rate_predicted AND j2_rate_predicted).\n")
cat("Pick the largest proposal_sd for which j2_rate_predicted (the binding\n")
cat("constraint, historically the lowest-acceptance threshold) clears 0.20.\n\n")

cat("Reminder: predictions rely on the single validated scale factor (", SCALE_FACTOR,
    ") derived from the proposal_sd=0.08 case. Treat as a guide, not a guarantee --\n")
cat("confirm the chosen value on the full n=50,790 data before treating the\n")
cat("resulting diagnostics as final.\n")