# =============================================================================
# bjcm_ridge_detail.R
#
# Computes the ESS and intercept/wave-effect correlation numbers needed to
# update the appendix's BJCM ridge paragraph (Appendix B.2) with the actual
# numbers from this run, in the same reporting style already used there:
# "R-hat = X, ESS = Y; r ~ Z".
#
# Does not modify 02_bjcm_model.R and does not re-run the sampler -- reads
# only the per-chain draws already stored in the fitted object.
#
# Usage:
#   source("04e_ridge_detail_appendix_b2.R")
#   ridge_detail <- bjcm_ridge_detail(raw)   # raw = full_result$results$results
#
# Requires: coda
# =============================================================================

bjcm_ridge_detail <- function(bjcm_raw) {
  
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required. Install with install.packages('coda').")
  }
  
  chains <- bjcm_raw$chains
  n_mediators <- ncol(chains[[1]]$alpha_med_samples)
  
  # ---- helper: ESS and correlation for one "intercept vs mean(U_w)" pair --
  # ESS is computed with coda::effectiveSize() on the pooled mcmc.list (sums
  # the per-chain effective sample sizes, matching how ESS is reported
  # elsewhere in this paper for other multi-chain quantities). The
  # correlation r is the intercept-vs-mean(U_w) Pearson correlation, computed
  # per chain and then averaged across chains (reported range shows how
  # stable it is across chains).
  ridge_stats <- function(intercept_draws_list, meanUw_draws_list, label) {
    r_per_chain <- mapply(function(icept, muw) cor(icept, muw),
                          intercept_draws_list, meanUw_draws_list)
    
    icept_mcmc <- coda::mcmc.list(lapply(intercept_draws_list, function(x)
      coda::mcmc(matrix(x, ncol = 1))))
    ess_icept <- sum(coda::effectiveSize(icept_mcmc))
    
    rhat_icept <- coda::gelman.diag(icept_mcmc, autoburnin = FALSE)$psrf[1, "Point est."]
    
    combo_mcmc <- coda::mcmc.list(mapply(function(icept, muw)
      coda::mcmc(matrix(icept + muw, ncol = 1)),
      intercept_draws_list, meanUw_draws_list, SIMPLIFY = FALSE))
    ess_combo <- sum(coda::effectiveSize(combo_mcmc))
    rhat_combo <- coda::gelman.diag(combo_mcmc, autoburnin = FALSE)$psrf[1, "Point est."]
    
    cat(sprintf(
      "%-12s intercept: R-hat=%.3f, ESS=%.1f | combo (intercept+mean Uw): R-hat=%.3f, ESS=%.1f | r per chain: %s (range %.3f to %.3f)\n",
      label, rhat_icept, ess_icept, rhat_combo, ess_combo,
      paste(round(r_per_chain, 3), collapse = ", "),
      min(r_per_chain), max(r_per_chain)
    ))
    
    list(rhat_intercept = rhat_icept, ess_intercept = ess_icept,
         rhat_combo = rhat_combo, ess_combo = ess_combo,
         r_per_chain = r_per_chain)
  }
  
  cat("=== BJCM ridge detail: intercept vs mean(U_w), per sub-model ===\n\n")
  
  out <- list()
  
  # Outcome model
  out[["outcome"]] <- ridge_stats(
    intercept_draws_list = lapply(chains, function(ch) ch$gamma_out_samples[, 1]),
    meanUw_draws_list    = lapply(chains, function(ch) rowMeans(ch$U_out_wave)),
    label = "Outcome"
  )
  
  mediator_names <- c("Depression", "Transport", "Mobility")
  for (k in seq_len(n_mediators)) {
    out[[paste0("mediator_", k)]] <- ridge_stats(
      intercept_draws_list = lapply(chains, function(ch) ch$beta_med_samples[, 1, k]),
      meanUw_draws_list    = lapply(chains, function(ch) rowMeans(ch$U_med_wave[, , k])),
      label = mediator_names[k]
    )
  }
  
  cat("\nCopy the printed lines above (or the returned list) back to Claude\n")
  cat("to update the appendix's BJCM ridge paragraph with these exact numbers.\n")
  
  invisible(out)
}