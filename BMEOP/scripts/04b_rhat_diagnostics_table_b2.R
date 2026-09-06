# =============================================================================
# 04b_rhat_diagnostics_table_b2.R
#
# Standalone multi-chain Gelman-Rubin R-hat/ESS extraction for an already-
# fitted BMEOP result object, reproducing Table B.2's per-coefficient layout
# (each covariate individually, the identified intercept+mean(U_w)
# combination, sigma2_u, and each threshold). Does not modify
# 02_bmeop_model.R and does not re-run the sampler.
#
# IMPORTANT: unlike BJCM's run_full_fit(), run_elsa_bmeop() returns the raw
# fitted object directly -- there is NO $results$results unwrapping needed.
#
# Usage:
#   source("04b_rhat_diagnostics_table_b2.R")
#   full_result_bmeop <- readRDS("../results/bmeop_full_results.rds")
#   rhat_results <- bmeop_full_rhat(full_result_bmeop)
#
# Requires: coda
# =============================================================================

bmeop_full_rhat <- function(bmeop_raw, zero_var_tol = 1e-10) {
  
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required. Install with install.packages('coda').")
  }
  
  chains <- bmeop_raw$chains
  if (is.null(chains) || length(chains) < 2) {
    stop("Fewer than 2 chains found. Pass the object returned directly by ",
         "run_elsa_bmeop()/bayesian_bmeop_elsa() (no unwrapping needed for ",
         "BMEOP), and make sure n_chains >= 2 was used.")
  }
  n_chains <- length(chains)
  covariate_names <- colnames(chains[[1]]$beta)
  if (is.null(covariate_names)) covariate_names <- paste0("beta_", seq_len(ncol(chains[[1]]$beta)))
  n_waves <- ncol(chains[[1]]$random_effects)
  n_thresh <- ncol(chains[[1]]$thresholds)
  
  safe_gelman_uni <- function(label, extractor) {
    vecs <- tryCatch(lapply(chains, extractor), error = function(e) {
      cat(label, ": FAILED while extracting data -", conditionMessage(e), "\n")
      NULL
    })
    if (is.null(vecs)) return(NULL)
    
    min_var <- min(sapply(vecs, var, na.rm = TRUE), na.rm = TRUE)
    if (is.na(min_var) || min_var < zero_var_tol) {
      cat(label, ": SKIPPED - ~zero within-chain variance in at least one chain",
          "(fixed/never updated)\n")
      return(NULL)
    }
    
    mcmc_list <- coda::mcmc.list(lapply(vecs, function(v) coda::mcmc(matrix(v, ncol = 1))))
    res <- tryCatch(
      coda::gelman.diag(mcmc_list, autoburnin = FALSE),
      error = function(e) {
        cat(label, ": FAILED in coda::gelman.diag() -", conditionMessage(e), "\n")
        NULL
      }
    )
    if (is.null(res)) return(NULL)
    
    ess <- sum(coda::effectiveSize(mcmc_list))
    list(rhat = res$psrf[1, "Point est."], rhat_ub = res$psrf[1, "Upper C.I."], ess = ess)
  }
  
  results <- list()
  
  cat("=== BMEOP per-coefficient R-hat/ESS (", n_chains, "chains ) ===\n\n")
  
  # ---- Raw intercept (column 1) -- expected to show the additive ridge ----
  results[["intercept"]] <- safe_gelman_uni("intercept", function(ch) ch$beta[, 1])
  
  # ---- Identified combination: intercept + mean(U_w) -----------------------
  results[["intercept_plus_mean_Uw"]] <- safe_gelman_uni("intercept_plus_mean_Uw", function(ch)
    ch$beta[, 1] + rowMeans(ch$random_effects))
  
  # ---- Remaining covariates (columns 2..p), individually --------------------
  # Falls back to actual X column names when available; otherwise uses the
  # known BMEOP covariate order (living alone, age, sex, education, health,
  # life satisfaction, depression, transport, mobility) as a labelled default.
  default_names <- c("livalone_TREATMENT", "age_gr", "dhsex2", "edqual2",
                     "health", "sclife", "depression", "transport_mobility",
                     "mobility_limitations")
  if (is.null(colnames(chains[[1]]$beta))) {
    covariate_names <- c("intercept", default_names)[seq_len(ncol(chains[[1]]$beta))]
  }
  for (j in seq(2, ncol(chains[[1]]$beta))) {
    results[[covariate_names[j]]] <- safe_gelman_uni(covariate_names[j], function(ch) ch$beta[, j])
  }
  
  # ---- Wave-variance component ----------------------------------------------
  results[["sigma2_u"]] <- safe_gelman_uni("sigma2_u", function(ch) ch$sigma2_u)
  
  # ---- Thresholds -------------------------------------------------------
  # Column k of ch$thresholds IS S_k (k = 1..n_thresh) -- column 1 is the
  # fixed tau_1 = 0 (see 01_threshold_update_mh.R: "new_thresholds[1] <- 0"),
  # not an absent/movable one, so it must be labelled threshold_j1 (and will
  # be auto-skipped as constant), with columns 2..n_thresh labelled
  # threshold_j2..threshold_j{n_thresh} directly (no +1 offset).
  for (k in seq_len(n_thresh)) {
    results[[paste0("threshold_j", k)]] <- safe_gelman_uni(
      paste0("threshold_j", k), function(ch) ch$thresholds[, k])
  }
  
  # ---- ATE (the paper's primary estimand) ------------------------------------
  results[["ATE"]] <- safe_gelman_uni("ATE", function(ch) ch$ate[!is.na(ch$ate)])
  
  cat("\n--- Summary ---\n\n")
  for (nm in names(results)) {
    r <- results[[nm]]
    if (is.null(r)) next
    cat(sprintf("%-24s R-hat=%.4f (UB=%.4f)  ESS=%.1f\n", nm, r$rhat, r$rhat_ub, r$ess))
  }
  
  failed <- names(results)[sapply(results, is.null)]
  if (length(failed) > 0) {
    cat("\nCould not be computed (see messages above):", paste(failed, collapse = ", "), "\n")
  }
  
  flagged <- names(results)[sapply(results, function(r) !is.null(r) && r$rhat > 1.1)]
  if (length(flagged) > 0) {
    cat("\nWARNING: R-hat > 1.1 for:", paste(flagged, collapse = ", "), "\n")
  } else {
    cat("\nAll computed R-hat values <= 1.1.\n")
  }
  
  invisible(results)
}