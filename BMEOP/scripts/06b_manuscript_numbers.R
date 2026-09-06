# =============================================================================
# bmeop_remaining_numbers.R
#
# Pulls the two BMEOP summary values not already printed by run_elsa_bmeop()'s
# console output: the identified Intercept + mean(U_w) combination, and
# sigma2_u. Needed to finish updating Table 3 in manuscript.tex.
#
# Usage:
#   source("06b_manuscript_numbers.R")
#   full_result_bmeop <- readRDS("../results/bmeop_full_results.rds")
#   bmeop_remaining_numbers(full_result_bmeop)
# =============================================================================

bmeop_remaining_numbers <- function(bmeop_raw) {
  fmt <- function(x, d = 4) formatC(x, digits = d, format = "f")
  
  chains <- bmeop_raw$chains
  
  combo <- unlist(lapply(chains, function(ch) ch$beta[, 1] + rowMeans(ch$random_effects)))
  cat(sprintf("Intercept + mean(U_w): mean=%s  sd=%s  95%% CrI=[%s, %s]\n",
              fmt(mean(combo)), fmt(sd(combo)), fmt(quantile(combo, 0.025)), fmt(quantile(combo, 0.975))))
  
  s2u <- unlist(lapply(chains, function(ch) ch$sigma2_u))
  cat(sprintf("sigma2_u: mean=%s  sd=%s  95%% CrI=[%s, %s]\n",
              fmt(mean(s2u)), fmt(sd(s2u)), fmt(quantile(s2u, 0.025)), fmt(quantile(s2u, 0.975))))
  
  invisible(NULL)
}