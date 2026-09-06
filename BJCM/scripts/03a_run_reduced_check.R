# =============================================================================
# Reduced-scale check: MH-corrected threshold update
# =============================================================================
# Uses the existing run_validation() (already a reduced-scale run: n=5,000
# subset, 4 chains, 2,000 iterations) with the posterior-invariant,
# marginal-likelihood MH threshold update (01_threshold_update_mh.R).
# Reports acceptance rates (target ~20%-50%) and threshold SD/range so you
# can confirm the chain is actually exploring before committing to a full run.
# =============================================================================

library(dplyr)
library(MASS)
library(mvtnorm)
library(coda)
library(Matrix)
library(parallel)

source("00_background_run_utils.R")
source("02_bjcm_model.R")   # sources threshold_update_mh.R internally

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

reset_mh_threshold_log()

cat("=== Running reduced-scale check (n=5000, 4 chains, 2000 iterations) ===\n\n")
t0 <- Sys.time()
mh_validation <- run_validation(elsa_long, subset_size = 5000, n_chains = 4)
t1 <- Sys.time()
elapsed_mins <- as.numeric(difftime(t1, t0, units = "mins"))
cat("\nRun time:", round(elapsed_mins, 2), "minutes\n\n")

saveRDS(mh_validation, "../results/mh_validation_result.rds")

cat("=== MH threshold acceptance rates (target ~20%-50%) ===\n")
acceptance_summary <- get_mh_acceptance_summary()
print(acceptance_summary)
saveRDS(acceptance_summary, "../results/mh_acceptance_summary.rds")

cat("\n=== Outcome threshold SD, per chain ===\n")
for (chain_i in 1:4) {
  th_out <- mh_validation$results$results$chains[[chain_i]]$thresholds_out
  sds <- apply(th_out[, -1, drop = FALSE], 2, sd, na.rm = TRUE)
  cat(sprintf("Chain %d: ", chain_i)); print(round(sds, 4))
}

cat("\n=== Outcome threshold range (min, max), chain 1 ===\n")
th_out_c1 <- mh_validation$results$results$chains[[1]]$thresholds_out
print(apply(th_out_c1[, -1, drop = FALSE], 2, function(x) round(range(x, na.rm = TRUE), 3)))

cat("\n=== Mediator threshold SD, per chain ===\n")
for (chain_i in 1:4) {
  th_med_list <- mh_validation$results$results$chains[[chain_i]]$thresholds_med
  for (k in seq_along(th_med_list)) {
    th_med_k <- th_med_list[[k]]
    if (!is.null(th_med_k) && is.matrix(th_med_k) && ncol(th_med_k) > 1) {
      sds <- apply(th_med_k[, -1, drop = FALSE], 2, sd, na.rm = TRUE)
      cat(sprintf("Chain %d, mediator %d: ", chain_i, k)); print(round(sds, 4))
    }
  }
}

cat("\nInterpretation:\n")
cat("  Acceptance rate 20%-50%: healthy MH mixing.\n")
cat("  Threshold SD/range should now be substantially larger than the\n")
cat("  exact-conditional attempt's 1e-4 to 1e-2 (relative to ~0.6 spacing) --\n")
cat("  if not, increase proposal_sd in threshold_update_mh.R and re-run.\n\n")

cat("Rough extrapolated full-run time (n=42,325, 15,000 iterations, 4 chains):\n")
iter_ratio <- 15000 / 2000
data_ratio <- 42325 / 5000
cat(round((elapsed_mins / 60) * iter_ratio * data_ratio, 1), "hours\n")