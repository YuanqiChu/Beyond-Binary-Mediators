# =============================================================================
# BMEOP reduced-scale check: MH-corrected threshold update
# =============================================================================
# IMPORTANT: run this in a SEPARATE R session from any BJCM script. BMEOP
# and BJCM define six identically-named helper functions with DIFFERENT
# implementations (rmvnorm_robust, rtruncnorm_robust, safe_matrix_inverse,
# update_U_wave_robust, update_sigma2_U_robust, update_thresholds_conservative)
# -- sourcing both 02_bmeop_model.R (here) and 02_bjcm_model.R (MH_BJCM) in the same session
# causes silent overwrites or hard errors (per the warning in the original
# the original 08a_run_bmeop.R). Never source both in one script or one interactive session.
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
missing_vars <- required_vars[!required_vars %in% names(elsa_long)]
if (length(missing_vars) > 0) {
  stop("elsa_long is missing required columns: ", paste(missing_vars, collapse = ", "))
}
cat("elsa_long loaded -", nrow(elsa_long), "rows,",
    dplyr::n_distinct(elsa_long$idauniq), "subjects\n\n")

# BMEOP uses the full Wave 2-7 complete-case panel directly (n=50,790), not
# the lagged Wave 3-7 subset BJCM uses.
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

# -----------------------------------------------------------------------------
# Re-tuning proposal_sd: the full-scale run (n=50,790) with proposal_sd=0.08
# gave acceptance rates of 5.9%-14.8%, well below the 20%-50% target -- the
# posterior for each threshold tightens by roughly 1/sqrt(n) as n grows
# (~10x here, i.e. ~0.32x), so 0.08 is likely too wide a step now. Rather than
# guess a single replacement value, sweep a small grid at the cheap n=5000
# scale first and read off which value lands in range, before committing to
# another 4.7-hour full run. NOTE: acceptance rate at n=5000 is only a rough
# guide for n=50,790 (the true target posterior is narrower at full scale
# than what this subset shows), so treat the result as a starting point, not
# a guarantee -- confirm again on the full data if time allows.
# -----------------------------------------------------------------------------

proposal_sd_grid <- c(0.08, 0.05, 0.04, 0.03, 0.025, 0.02)

grid_results <- list()

for (psd in proposal_sd_grid) {
  
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
  
  grid_results[[as.character(psd)]] <- list(fit = fit_psd, acceptance = acc,
                                            elapsed_mins = elapsed_mins)
}

cat("\n\n=== Summary: mean acceptance rate across thresholds, by proposal_sd ===\n")
summary_tbl <- do.call(rbind, lapply(names(grid_results), function(psd_str) {
  data.frame(proposal_sd = as.numeric(psd_str),
             mean_rate = mean(grid_results[[psd_str]]$acceptance$rate))
}))
summary_tbl <- summary_tbl[order(-summary_tbl$proposal_sd), ]
print(summary_tbl, row.names = FALSE)
cat("\nTarget range: 20%-50%. Pick the largest proposal_sd whose mean_rate\n")
cat("falls in that range (larger step size = better mixing, all else equal),\n")
cat("then use that value for the full-scale re-run.\n\n")

saveRDS(grid_results, "../results/mh_bmeop_proposal_sd_grid_result.rds")
cat("Saved: ../results/mh_bmeop_proposal_sd_grid_result.rds\n\n")

# ---- Keep the rest of this script's diagnostics pointed at the ORIGINAL
#      proposal_sd=0.08 result (first grid entry), preserving prior behaviour
#      for anyone re-running this script's downstream checks unmodified ----
mh_bmeop_validation <- grid_results[["0.08"]]$fit
acceptance_summary  <- grid_results[["0.08"]]$acceptance
saveRDS(mh_bmeop_validation, "../results/mh_bmeop_validation_result.rds")
saveRDS(acceptance_summary, "../results/mh_bmeop_acceptance_summary.rds")

# BMEOP has a single threshold set (the loneliness outcome, label
# "bmeop_outcome" in 02_bmeop_model.R), unlike BJCM which has separate
# mediator and outcome thresholds.
cat("\n=== Threshold SD, per chain ===\n")
for (chain_i in 1:4) {
  th <- mh_bmeop_validation$chains[[chain_i]]$thresholds
  if (is.null(th)) {
    cat(sprintf("Chain %d: could not find $thresholds -- inspect structure below\n", chain_i))
    next
  }
  sds <- apply(th[, -1, drop = FALSE], 2, sd, na.rm = TRUE)
  cat(sprintf("Chain %d: ", chain_i)); print(round(sds, 4))
}

cat("\n=== Structure check (adjust the threshold field name above if needed) ===\n")
str(mh_bmeop_validation$chains[[1]], max.level = 1)

cat("\nInterpretation:\n")
cat("  Acceptance rate 20%-50%: healthy MH mixing.\n")
cat("  If threshold SD looks pathologically small (< 0.01, relative to the\n")
cat("  threshold spacing seen in the structure above), adjust proposal_sd\n")
cat("  in the 02_bmeop_model.R call site and re-run.\n\n")

cat("Rough extrapolated full-run time (n=50,790, 15,000 iterations, 4 chains):\n")
iter_ratio <- 15000 / 2000
data_ratio <- 50790 / 5000
cat(round((elapsed_mins / 60) * iter_ratio * data_ratio, 1), "hours\n")