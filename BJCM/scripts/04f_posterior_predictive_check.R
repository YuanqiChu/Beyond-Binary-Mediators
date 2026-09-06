# =============================================================================
# Posterior predictive check for the outcome model
# =============================================================================
# Compares the observed distribution of loneliness (categories 3-9) against
# the posterior predictive distribution implied by the fitted BJCM outcome
# model, using each individual's own observed treatment and covariates.
# Reports the full comparison table (observed, predicted, 95% posterior
# predictive interval, difference in percentage points per category) and
# regenerates Appendix Figure B.2 (ppc_check.pdf).
#
# Read-only against the saved fit: no model is re-fit and no counterfactual
# simulation is re-run.
#
# Run from scripts/, so that ../data/ and ../results/ resolve correctly.
# =============================================================================

library(ggplot2)

source("00_background_run_utils.R")
source("02_bjcm_model.R")

elsa_long <- read.csv("../data/elsa_longitudinal_analysis.csv")
inputs <- prepare_bjcm_inputs(elsa_long)

# ---- Observed category proportions -----------------------------------------
outcome_raw <- inputs$y
categories  <- sort(unique(outcome_raw))
if (!all(categories == 3:9)) {
  stop("Unexpected outcome categories: ", paste(categories, collapse = ", "),
       " -- check inputs$y coding before proceeding.")
}

observed_prop <- as.numeric(table(factor(outcome_raw, levels = categories))) /
  length(outcome_raw)

# ---- Posterior predictive category probabilities ---------------------------
full_result <- readRDS("../results/bjcm_full_results.rds")
chains <- full_result$results$results$chains

pred_probs_pooled <- do.call(rbind, lapply(chains, function(ch) ch$y_pred_probs_outcome))

if (ncol(pred_probs_pooled) != length(categories)) {
  stop("y_pred_probs_outcome has ", ncol(pred_probs_pooled),
       " columns; expected ", length(categories), " (one per category).")
}

predicted_mean <- colMeans(pred_probs_pooled)
predicted_lo   <- apply(pred_probs_pooled, 2, quantile, probs = 0.025)
predicted_hi   <- apply(pred_probs_pooled, 2, quantile, probs = 0.975)

# ---- Comparison table -------------------------------------------------------
diff_pp <- round(100 * (observed_prop - predicted_mean), 2)

comparison <- data.frame(
  category  = categories,
  observed  = round(100 * observed_prop, 2),
  predicted = round(100 * predicted_mean, 2),
  lo        = round(100 * predicted_lo, 2),
  hi        = round(100 * predicted_hi, 2),
  diff_pp   = diff_pp
)

cat("============ POSTERIOR PREDICTIVE CHECK ============\n")
print(comparison, row.names = FALSE)
cat("\nMax |difference| (percentage points):", max(abs(diff_pp)),
    "at category", categories[which.max(abs(diff_pp))], "\n")

# ---- Figure B.2 --------------------------------------------------------------
df <- data.frame(
  category  = factor(categories),
  observed  = observed_prop,
  predicted = predicted_mean,
  lo        = predicted_lo,
  hi        = predicted_hi
)

p <- ggplot(df, aes(x = category)) +
  geom_col(aes(y = observed), fill = "grey70", colour = "black",
           width = 0.65, linewidth = 0.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15,
                colour = "#D2691E", linewidth = 0.6) +
  geom_point(aes(y = predicted), colour = "#D2691E", size = 3) +
  labs(x = "Loneliness category", y = "Proportion") +
  scale_y_continuous(limits = c(0, max(df$observed, df$hi) * 1.1), expand = c(0, 0)) +
  theme_bw(base_size = 13, base_family = "serif") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")
  )

ggsave("../results/ppc_check.pdf", p, width = 6.5, height = 4.2, device = "pdf")
cat("\nSaved ppc_check.pdf\n")
