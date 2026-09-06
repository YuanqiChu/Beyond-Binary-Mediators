# =============================================================================
# bjcm_rhat_diagnostics.R
#
# Standalone multi-chain Gelman-Rubin R-hat extraction for an already-fitted
# BJCM result object. Does not modify or re-source 02_bjcm_model.R and does
# not re-run the sampler -- it only reads the per-chain draws already stored
# in the fitted object.
#
# Robust version: each parameter block is computed independently. A block
# that would break coda::gelman.diag()'s internal chol() decomposition
# (because one or more columns have ~zero within-chain variance in at least
# one chain -- e.g. a parameter that barely moved, or is fixed by
# construction) has those columns automatically dropped, with a note printed
# identifying which columns and why. A block that still fails for some other
# reason is reported with its error message rather than stopping the whole
# diagnostic run.
#
# Usage:
#   source("04b_rhat_diagnostics_standalone.R")
#   full_result <- readRDS("../results/bjcm_full_results.rds")
#   raw <- full_result$results$results   # unwrap run_full_fit()'s triple nesting
#   rhat_results <- bjcm_full_rhat(raw)
#
# Requires: coda
# =============================================================================

bjcm_full_rhat <- function(bjcm_raw, zero_var_tol = 1e-10) {
  
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required. Install with install.packages('coda').")
  }
  
  chains <- bjcm_raw$chains
  if (is.null(chains) || length(chains) < 2) {
    stop("Fewer than 2 chains found. Pass the RAW model output ",
         "(e.g. full_result$results$results from run_full_fit()), not the ",
         "wrapped structure, and make sure n_chains >= 2 was used.")
  }
  n_chains    <- length(chains)
  n_mediators <- ncol(chains[[1]]$alpha_med_samples)
  p_baseline  <- dim(chains[[1]]$beta_med_samples)[2]
  n_waves     <- dim(chains[[1]]$U_out_wave)[2]
  
  person_avg <- function(x) if (is.matrix(x)) rowMeans(x) else as.numeric(x)
  
  # ---- Core helper: build an mcmc.list from a per-chain matrix extractor,
  # drop any column with ~zero within-chain variance in ANY chain (this is
  # exactly what makes coda::gelman.diag()'s internal chol(W) fail with
  # "leading minor of order k is not positive"), then run gelman.diag()
  # inside tryCatch so one bad block cannot abort the rest of the run. ------
  safe_gelman <- function(label, extractor) {
    mats <- tryCatch(lapply(chains, extractor), error = function(e) {
      cat(label, ": FAILED while extracting data -", conditionMessage(e), "\n\n")
      NULL
    })
    if (is.null(mats)) return(NULL)
    
    ncol_check <- unique(sapply(mats, ncol))
    if (length(ncol_check) != 1) {
      cat(label, ": FAILED - chains have differing numbers of columns (",
          paste(ncol_check, collapse = ", "), ")\n\n")
      return(NULL)
    }
    
    all_cols <- colnames(mats[[1]])
    if (is.null(all_cols)) all_cols <- paste0("V", seq_len(ncol(mats[[1]])))
    
    # Per-chain, per-column variance; flag a column if it is ~constant in
    # ANY chain (that alone is enough to make the pooled W singular).
    var_by_chain <- sapply(mats, function(m) apply(m, 2, var, na.rm = TRUE))
    if (is.null(dim(var_by_chain))) var_by_chain <- matrix(var_by_chain, nrow = 1)
    min_var_across_chains <- apply(var_by_chain, 1, min, na.rm = TRUE)
    
    constant_cols <- which(min_var_across_chains < zero_var_tol | is.na(min_var_across_chains))
    if (length(constant_cols) > 0) {
      cat(label, ": dropping", length(constant_cols), "column(s) with ~zero within-chain",
          "variance in at least one chain (fixed/never updated):",
          paste(all_cols[constant_cols], collapse = ", "), "\n")
      keep_cols <- setdiff(seq_along(all_cols), constant_cols)
      if (length(keep_cols) == 0) {
        cat(label, ": all columns dropped, nothing left to check.\n\n")
        return(NULL)
      }
      mats <- lapply(mats, function(m) {
        m2 <- m[, keep_cols, drop = FALSE]
        colnames(m2) <- all_cols[keep_cols]
        m2
      })
    }
    
    mcmc_list <- coda::mcmc.list(lapply(mats, coda::mcmc))
    
    res <- tryCatch(
      coda::gelman.diag(mcmc_list, autoburnin = FALSE),
      error = function(e) {
        cat(label, ": FAILED in coda::gelman.diag() -", conditionMessage(e), "\n")
        cat("  (this usually means near-collinearity between two or more\n")
        cat("  columns, not a single constant column -- e.g. thresholds\n")
        cat("  that barely separate, or two coefficients trading off\n")
        cat("  against each other across iterations)\n\n")
        NULL
      }
    )
    res
  }
  
  results <- list()
  
  # ---- Causal estimands ---------------------------------------------------
  results[["TE"]] <- safe_gelman("TE", function(ch)
    matrix(person_avg(ch$total_effects), ncol = 1, dimnames = list(NULL, "TE")))
  results[["NDE"]] <- safe_gelman("NDE", function(ch)
    matrix(person_avg(ch$direct_effects), ncol = 1, dimnames = list(NULL, "NDE")))
  results[["NIE"]] <- safe_gelman("NIE", function(ch)
    matrix(person_avg(ch$indirect_effects), ncol = 1, dimnames = list(NULL, "NIE")))
  
  # ---- Structural parameters -----------------------------------------------
  results[["lambda_out"]] <- safe_gelman("lambda_out", function(ch)
    matrix(ch$lambda_out_samples, ncol = 1, dimnames = list(NULL, "lambda_out")))
  
  results[["alpha_med"]] <- safe_gelman("alpha_med", function(ch) {
    m <- ch$alpha_med_samples
    colnames(m) <- paste0("alpha_med_", seq_len(n_mediators))
    m
  })
  
  results[["beta_out_med"]] <- safe_gelman("beta_out_med", function(ch) {
    m <- ch$beta_out_med_samples
    colnames(m) <- paste0("beta_out_med_", seq_len(n_mediators))
    m
  })
  
  results[["beta_med"]] <- safe_gelman("beta_med", function(ch) {
    m <- matrix(ch$beta_med_samples, nrow = dim(ch$beta_med_samples)[1])
    colnames(m) <- paste0("beta_med_p", rep(seq_len(p_baseline), times = n_mediators),
                          "_med", rep(seq_len(n_mediators), each = p_baseline))
    m
  })
  
  results[["gamma_out"]] <- safe_gelman("gamma_out", function(ch) {
    m <- ch$gamma_out_samples
    colnames(m) <- paste0("gamma_out_p", seq_len(ncol(m)))
    m
  })
  
  results[["thresholds_out"]] <- safe_gelman("thresholds_out", function(ch) {
    m <- ch$thresholds_out
    colnames(m) <- paste0("S_out_", seq_len(ncol(m)))
    m
  })
  
  for (k in seq_len(n_mediators)) {
    results[[paste0("thresholds_med_", k)]] <- safe_gelman(
      paste0("thresholds_med_", k),
      function(ch) {
        m <- ch$thresholds_med[[k]]
        colnames(m) <- paste0("S_med", k, "_", seq_len(ncol(m)))
        m
      }
    )
  }
  
  results[["U_out_wave"]] <- safe_gelman("U_out_wave", function(ch) {
    m <- ch$U_out_wave
    colnames(m) <- paste0("U_out_wave_", seq_len(n_waves))
    m
  })
  
  results[["U_med_wave"]] <- safe_gelman("U_med_wave", function(ch) {
    m <- matrix(ch$U_med_wave, nrow = dim(ch$U_med_wave)[1])
    colnames(m) <- paste0("U_med_wave", rep(seq_len(n_waves), times = n_mediators),
                          "_med", rep(seq_len(n_mediators), each = n_waves))
    m
  })
  
  results[["sigma2_U_out"]] <- safe_gelman("sigma2_U_out", function(ch)
    matrix(ch$sigma2_U_out, ncol = 1, dimnames = list(NULL, "sigma2_U_out")))
  
  results[["sigma2_U_med"]] <- safe_gelman("sigma2_U_med", function(ch) {
    m <- ch$sigma2_U_med
    colnames(m) <- paste0("sigma2_U_med_", seq_len(n_mediators))
    m
  })
  
  # ---- Intercept + mean(U_w) reconciled combination (per B.2 discussion) ----
  # The R-hat pattern for gamma_out_p1 / U_out_wave (and, per mediator,
  # beta_med_p1_med_k / U_med_wave*_med_k) matching almost exactly is the
  # signature of the additive non-identifiability ridge documented in
  # Appendix B.2/C.1: individually neither the raw intercept nor the
  # wave-level random effects converge, but their sum does. This block
  # computes that reconciled combination automatically for BOTH the outcome
  # model and each mediator sub-model, ASSUMING column 1 of
  # gamma_out_samples / beta_med_samples[,,k] is the intercept (the usual
  # X_baseline convention). If that assumption is wrong for your data setup,
  # the results below will look wrong (R-hat still bad) rather than silently
  # misleading -- confirm against your actual X_baseline column order if so.
  cat("--- Reconciled intercept + mean(U_w) diagnostics ---\n")
  cat("(assumes column 1 of gamma_out_samples / beta_med_samples is the\n")
  cat("intercept -- confirm against your X_baseline column order)\n\n")
  
  results[["intercept_plus_mean_Uw_out"]] <- safe_gelman("intercept_plus_mean_Uw_out", function(ch)
    matrix(ch$gamma_out_samples[, 1] + rowMeans(ch$U_out_wave), ncol = 1,
           dimnames = list(NULL, "intercept_plus_mean_Uw_out")))
  
  for (k in seq_len(n_mediators)) {
    results[[paste0("intercept_plus_mean_Uw_med_", k)]] <- safe_gelman(
      paste0("intercept_plus_mean_Uw_med_", k),
      function(ch) {
        matrix(ch$beta_med_samples[, 1, k] + rowMeans(ch$U_med_wave[, , k]), ncol = 1,
               dimnames = list(NULL, paste0("intercept_plus_mean_Uw_med_", k)))
      }
    )
  }
  cat("\n")
  
  # ---- Summary --------------------------------------------------------------
  cat("\n=== Multi-chain Gelman-Rubin R-hat (", n_chains, "chains ) ===\n\n")
  computed <- results[!sapply(results, is.null)]
  for (nm in names(computed)) {
    cat(nm, ":\n")
    print(computed[[nm]])
    cat("\n")
  }
  
  failed <- names(results)[sapply(results, is.null)]
  if (length(failed) > 0) {
    cat("Blocks that could not be computed (see messages above):",
        paste(failed, collapse = ", "), "\n\n")
  }
  
  flagged <- Filter(function(nm) {
    any(computed[[nm]]$psrf[, "Point est."] > 1.1)
  }, names(computed))
  
  if (length(flagged) > 0) {
    cat("WARNING: R-hat > 1.1 for:", paste(flagged, collapse = ", "), "\n")
  } else if (length(computed) > 0) {
    cat("All computed R-hat values <= 1.1.\n")
  }
  
  results
}