# =============================================================================
# bjcm_manuscript_numbers.R
#
# Extracts every number manuscript.tex currently hardcodes for the BJCM
# results (TE/NDE/NIE with CrI, per-mediator a-path/b-path/NIE with CrI,
# mediator % contributions, ESS, decomposition error) directly from the
# fitted object's stored summaries, at full precision -- not from the
# rounded console printout.
#
# Does not modify 02_bjcm_model.R and does not re-run the sampler.
#
# Usage:
#   source("06b_manuscript_numbers.R")
#   full_result <- readRDS("../results/bjcm_full_results.rds")
#   raw <- full_result$results$results
#   bjcm_manuscript_numbers(raw)
#
# Requires: coda
# =============================================================================

bjcm_manuscript_numbers <- function(bjcm_raw) {
  
  fmt <- function(x, d = 4) formatC(x, digits = d, format = "f")
  
  cat("=== Panel B: TE / NDE / NIE (full precision) ===\n\n")
  for (nm in c("total_effect_summary", "direct_effect_summary", "indirect_effect_summary")) {
    s <- bjcm_raw[[nm]]
    cat(sprintf("%-24s mean=%s sd=%s  95%% CrI=[%s, %s]\n",
                nm, fmt(s["mean"]), fmt(s["sd"], 5), fmt(s["q2.5"]), fmt(s["q97.5"])))
  }
  
  te_mean  <- bjcm_raw$total_effect_summary["mean"]
  nde_mean <- bjcm_raw$direct_effect_summary["mean"]
  nie_mean <- bjcm_raw$indirect_effect_summary["mean"]
  cat(sprintf("\nNDE %% of TE = %s%%   NIE %% of TE = %s%%\n",
              fmt(100 * nde_mean / te_mean, 1), fmt(100 * nie_mean / te_mean, 1)))
  
  cat("\n=== Panel C: per-mediator a-path / b-path / NIE (full precision) ===\n\n")
  mediator_names <- c("Depression", "Transport", "Mobility")
  n_mediators <- length(bjcm_raw$alpha_med_summaries)
  nie_total_mean <- nie_mean
  for (k in seq_len(n_mediators)) {
    a <- bjcm_raw$alpha_med_summaries[[k]]
    b <- bjcm_raw$beta_out_med_summaries[[k]]
    n <- bjcm_raw$indirect_effect_individual_summaries[[k]]
    cat(sprintf("%s:\n", mediator_names[k]))
    cat(sprintf("  a-path  mean=%s  95%% CrI=[%s, %s]\n", fmt(a["mean"]), fmt(a["q2.5"]), fmt(a["q97.5"])))
    cat(sprintf("  b-path  mean=%s  95%% CrI=[%s, %s]\n", fmt(b["mean"]), fmt(b["q2.5"]), fmt(b["q97.5"])))
    cat(sprintf("  NIE_k   mean=%s  95%% CrI=[%s, %s]\n", fmt(n["mean"]), fmt(n["q2.5"]), fmt(n["q97.5"])))
    cat(sprintf("  %% of NIE_total = %s%%   %% of TE = %s%%\n\n",
                fmt(100 * n["mean"] / nie_total_mean, 1), fmt(100 * n["mean"] / te_mean, 1)))
  }
  
  cat("=== Figure 2: lambda (direct effect coefficient) ===\n\n")
  l <- bjcm_raw$lambda_out_summary
  cat(sprintf("lambda mean=%s  95%% CrI=[%s, %s]\n\n", fmt(l["mean"]), fmt(l["q2.5"]), fmt(l["q97.5"])))
  
  cat("=== Decomposition validation (exact, not rounded) ===\n\n")
  dv <- bjcm_raw$decomposition_validation
  cat(sprintf("TE - NDE - NIE           = %.10g\n", dv$decomposition_error))
  cat(sprintf("Total NIE - Sum(NIE_k)   = %.10g\n\n", dv$nie_sum_error))
  
  cat("=== ESS (ATE/NDE/NIE and lambda), pooled across chains ===\n\n")
  if (requireNamespace("coda", quietly = TRUE)) {
    person_avg <- function(x) if (is.matrix(x)) rowMeans(x) else as.numeric(x)
    mk <- function(extractor) coda::mcmc.list(lapply(bjcm_raw$chains, function(ch)
      coda::mcmc(matrix(extractor(ch), ncol = 1))))
    ess_te     <- sum(coda::effectiveSize(mk(function(ch) person_avg(ch$total_effects))))
    ess_nde    <- sum(coda::effectiveSize(mk(function(ch) person_avg(ch$direct_effects))))
    ess_nie    <- sum(coda::effectiveSize(mk(function(ch) person_avg(ch$indirect_effects))))
    ess_lambda <- sum(coda::effectiveSize(mk(function(ch) ch$lambda_out_samples)))
    cat(sprintf("ESS: TE=%s  NDE=%s  NIE=%s  lambda=%s\n\n",
                fmt(ess_te, 1), fmt(ess_nde, 1), fmt(ess_nie, 1), fmt(ess_lambda, 1)))
  } else {
    cat("coda not available -- skipping ESS.\n\n")
  }
  
  cat("Copy the numbers above back to Claude to update manuscript.tex.\n")
  
  invisible(NULL)
}