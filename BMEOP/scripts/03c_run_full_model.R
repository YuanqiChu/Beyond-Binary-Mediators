# =============================================================================
# Full-scale BMEOP model fit
# =============================================================================
# IMPORTANT: run in a separate R session from any BJCM script (see
# 03a_run_reduced_check.R for why).
#
# Run this only after 03a_run_reduced_check.R and 03b_proposal_sd_grid_search.R
# have shown healthy threshold-update acceptance rates (~20-50%) at the
# chosen proposal_sd.
# =============================================================================

library(dplyr)
library(MASS)
library(mvtnorm)
library(coda)
library(Matrix)
library(parallel)

source("00_background_run_utils.R")
source("02_bmeop_model.R")

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

# proposal_sd = 0.015, chosen by the grid search in
# 03b_proposal_sd_grid_search.R: the largest step size in the tested grid
# whose predicted full-scale (n=50,790) acceptance rate on the hardest
# threshold (j2) cleared the 20% target (predicted ~21.1%, scaled from the
# n=5,000 reduced-scale check using the empirically validated factor
# 3.06 = [n=5,000 mean rate at 0.08] / [full-run mean rate at 0.08]).
run_bmeop_full_refit <- function(complete_data, n_chains = 4, proposal_sd = 0.015) {
  X <- build_bmeop_X(complete_data)
  bayesian_bmeop_elsa(
    y = complete_data$loneliness, X = X, wave = complete_data$wave, treatment_col = 2,
    n_iter = 15000, n_burn = 10000, n_thin = 5,
    n_chains = n_chains, proposal_sd = proposal_sd, seed = 123, verbose = TRUE
  )
}

reset_mh_threshold_log()

cat("=== Running full-scale BMEOP fit ===\n")
cat("proposal_sd = 0.015 (see comment above run_bmeop_full_refit())\n")
cat("Consider running in the background instead of interactively:\n")
cat('  job <- run_in_background(run_bmeop_full_refit(bmeop_complete, n_chains = 4, proposal_sd = 0.015))\n')
cat('  monitor_progress("bmeop")\n')
cat("  full_result <- collect_result(job)\n\n")

t0 <- Sys.time()
full_result <- run_bmeop_full_refit(bmeop_complete, n_chains = 4, proposal_sd = 0.015)
t1 <- Sys.time()
cat("\nFull run wall-clock time:", round(as.numeric(difftime(t1, t0, units = "hours")), 2), "hours\n")

saveRDS(full_result, "../results/bmeop_full_results.rds")

acceptance_summary <- get_mh_acceptance_summary()
cat("\n=== Threshold acceptance rates (full run) ===\n")
print(acceptance_summary)
saveRDS(acceptance_summary, "../results/mh_bmeop_acceptance_summary.rds")

cat("\n=== Done. Files saved: ===\n")
cat("  ../results/bmeop_full_results.rds\n")
cat("  ../results/mh_bmeop_acceptance_summary.rds\n\n")
cat("Next: run 06a_ordinal_scale_profile_table3.R and 06b_manuscript_numbers.R\n")
cat("to extract the manuscript's reported numbers.\n")
