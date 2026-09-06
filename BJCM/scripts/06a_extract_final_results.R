# =============================================================================
# BJCM results extraction
# =============================================================================
# Read-only extraction of the manuscript's reported numbers from the saved
# fit: TE/NDE/NIE (Table 4, Panel B), mediator-specific NIE and shares
# (Panel C), a-path/b-path coefficients (Panel A), the decomposition
# validation check, and the E-value sensitivity analysis (Section 4.6).
# Only readRDS() is used -- no model file is sourced and nothing is re-run.
# =============================================================================

results_path <- "../results/bjcm_full_results.rds"
if (!file.exists(results_path)) stop("Could not find '", results_path, "'")

full <- readRDS(results_path)
res  <- full$results$results

cat("=============================================================\n")
cat("PART 1: TE / NDE / NIE (Table 4, Panel B)\n")
cat("=============================================================\n\n")

extract_effect <- function(x, label) {
  cat(sprintf("%s: %.4f [%.4f, %.4f]\n", label, x["mean"], x["q2.5"], x["q97.5"]))
}

extract_effect(res$total_effect_summary,    "TE ")
extract_effect(res$direct_effect_summary,   "NDE")
extract_effect(res$indirect_effect_summary, "NIE")

te  <- res$total_effect_summary["mean"]
nde <- res$direct_effect_summary["mean"]
nie <- res$indirect_effect_summary["mean"]

cat(sprintf("\nNDE as %% of TE: %.1f%%\n", nde / te * 100))
cat(sprintf("NIE as %% of TE: %.1f%%\n", nie / te * 100))

outcome_sd <- 1.634  # observed loneliness SD, as used throughout the manuscript
cat(sprintf("\nTE in SD units: %.3f SD\n", te / outcome_sd))
cat(sprintf("TE as %% of scale range (3-9, range=6): %.1f%%\n", te / 6 * 100))

cat("\n\n=============================================================\n")
cat("PART 2: Mediator-specific NIE and shares (Table 4, Panel C)\n")
cat("=============================================================\n\n")

mediator_names <- c("Depression", "Transport mobility", "Mobility limitations")
nie_indiv <- res$indirect_effect_individual_summaries

for (k in 1:3) {
  m <- nie_indiv[[k]]
  share_of_nie <- m["mean"] / nie * 100
  share_of_te  <- m["mean"] / te * 100
  cat(sprintf("%-22s NIE = %.4f [%.4f, %.4f]  (%.1f%% of NIE, %.1f%% of TE)\n",
              mediator_names[k], m["mean"], m["q2.5"], m["q97.5"],
              share_of_nie, share_of_te))
}

cat("\n\n=============================================================\n")
cat("PART 3: a-path and b-path coefficients per mediator (Table 4, Panel A)\n")
cat("=============================================================\n\n")

for (k in 1:3) {
  a <- res$alpha_med_summaries[[k]]
  b <- res$beta_out_med_summaries[[k]]
  cat(sprintf("%-22s a-path (alpha%d) = %.4f [%.4f, %.4f]\n",
              mediator_names[k], k, a["mean"], a["q2.5"], a["q97.5"]))
  cat(sprintf("%-22s b-path (beta%d)  = %.4f [%.4f, %.4f]\n",
              mediator_names[k], k, b["mean"], b["q2.5"], b["q97.5"]))
}
lam <- res$lambda_out_summary
cat(sprintf("Direct effect (lambda) = %.4f [%.4f, %.4f]\n", lam["mean"], lam["q2.5"], lam["q97.5"]))

cat("\n\n=============================================================\n")
cat("PART 4: Decomposition validation (TE - NDE - sum(NIE_k) ~ 0)\n")
cat("=============================================================\n\n")

dv <- res$decomposition_validation
cat(sprintf("TE - NDE - NIE         = %.10g\n", dv$decomposition_error))
cat(sprintf("Total NIE - sum(NIE_k) = %.10g\n", dv$nie_sum_error))

cat("\n\n=============================================================\n")
cat("PART 5: E-value sensitivity analysis (Section 4.6)\n")
cat("=============================================================\n\n")
cat("RR_approx = exp(0.91 * d), d = TE / outcome_sd, E = RR + sqrt(RR*(RR-1))\n\n")

compute_e_value <- function(effect_estimate, sd_outcome) {
  d <- effect_estimate / sd_outcome
  RR_approx <- exp(0.91 * d)
  if (RR_approx >= 1) {
    e_value <- RR_approx + sqrt(RR_approx * (RR_approx - 1))
  } else {
    RR_inv <- 1 / RR_approx
    e_value <- RR_inv + sqrt(RR_inv * (RR_inv - 1))
  }
  list(d = d, RR_approx = RR_approx, e_value = e_value)
}

te_point <- res$total_effect_summary["mean"]
te_lower <- res$total_effect_summary["q2.5"]

point_result <- compute_e_value(te_point, outcome_sd)
ci_result    <- compute_e_value(te_lower, outcome_sd)

cat(sprintf("TE point estimate: %.4f\n", te_point))
cat(sprintf("  d = %.4f, RR_approx = %.4f, E-value = %.2f\n",
            point_result$d, point_result$RR_approx, point_result$e_value))
cat(sprintf("\nTE 95%% CrI lower bound: %.4f\n", te_lower))
cat(sprintf("  d = %.4f, RR_approx = %.4f, E-value (CI bound) = %.2f\n",
            ci_result$d, ci_result$RR_approx, ci_result$e_value))

cat("\n\n=== EXTRACTION COMPLETE ===\n")
