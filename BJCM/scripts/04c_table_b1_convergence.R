# =============================================================================
# table_b1_bjcm_convergence.R
# =============================================================================
# Appendix Table B.1: full Gelman-Rubin diagnostics (point estimate + 97.5%
# upper bound) for the three primary BJCM causal estimands (TE, NDE, NIE).
# Read-only against the already-saved MH-corrected full run -- no re-fitting.
#
# Run from MH_BJCM/scripts/.
# =============================================================================

library(coda)

mh <- readRDS("../results/bjcm_full_results.rds")
chains <- mh$results$results$chains

te_list  <- lapply(chains, function(ch) ch$total_effects)
nde_list <- lapply(chains, function(ch) ch$direct_effects)
nie_list <- lapply(chains, function(ch) ch$indirect_effects)

te_mcmc  <- mcmc.list(lapply(te_list,  mcmc))
nde_mcmc <- mcmc.list(lapply(nde_list, mcmc))
nie_mcmc <- mcmc.list(lapply(nie_list, mcmc))

te_diag  <- gelman.diag(te_mcmc)
nde_diag <- gelman.diag(nde_mcmc)
nie_diag <- gelman.diag(nie_mcmc)

cat("=============================================================\n")
cat("TABLE B.1 -- BJCM causal estimand convergence (point + upper CI)\n")
cat("=============================================================\n\n")

cat(sprintf("%-30s %-14s %-14s\n", "Estimand", "R-hat (point)", "R-hat (97.5% UB)"))
cat(sprintf("%-30s %-14.3f %-14.3f\n", "Total Effect (TE)",
            te_diag$psrf[1, 1], te_diag$psrf[1, 2]))
cat(sprintf("%-30s %-14.3f %-14.3f\n", "Natural Direct Effect (NDE)",
            nde_diag$psrf[1, 1], nde_diag$psrf[1, 2]))
cat(sprintf("%-30s %-14.3f %-14.3f\n", "Natural Indirect Effect (NIE)",
            nie_diag$psrf[1, 1], nie_diag$psrf[1, 2]))

cat("\nESS (for cross-check against main text):\n")
cat(sprintf("TE:  %.1f\nNDE: %.1f\nNIE: %.1f\n",
            effectiveSize(te_mcmc), effectiveSize(nde_mcmc), effectiveSize(nie_mcmc)))