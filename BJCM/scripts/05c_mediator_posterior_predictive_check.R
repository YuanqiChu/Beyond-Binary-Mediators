# =============================================================================
# Posterior predictive check for the three mediator models (BJCM)
# =============================================================================
# Compares the observed category distribution of each lagged mediator
# (depression: CES-D 0-7, 8 categories; transport mobility: binary, 2
# categories; mobility limitations: 0-3, 4 categories) against the posterior
# predictive distribution implied by each mediator's fitted BJCM sub-model,
# using each individual's own observed treatment and baseline covariates.
#
# Mirrors 04f_posterior_predictive_check.R exactly in structure and style:
# reports the full comparison table (observed, predicted, 95% posterior
# predictive interval, difference in percentage points per category) for
# each mediator, and produces a three-panel figure (mediator_ppc_check.pdf)
# in the same visual style as Figure B.2.
#
# Read-only against the saved fit: no model is re-fit and no counterfactual
# simulation is re-run. y_pred_probs_mediators is already stored per chain
# in bjcm_full_results.rds, generated at the same time as
# y_pred_probs_outcome, so this script only pools, summarises, and plots.
#
# Unlike 04f (which re-reads the raw CSV via prepare_bjcm_inputs to get the
# observed outcome), the observed mediator values are taken directly from
# lagged_data$analysis_data saved inside bjcm_full_results.rds. This is the
# exact analysis sample the model was fit to (42,325 obs / 8,465 participants),
# so no re-reading of ../data/elsa_longitudinal_analysis.csv is needed.
#
# Run from scripts/, so that ../results/ resolves correctly.
# =============================================================================

library(ggplot2)

# ---- Mediator metadata ------------------------------------------------------
# Order must match the mediator ordering used when the model was fit
# (chains[[c]]$y_pred_probs_mediators[[k]] and thresholds_med[[k]]):
#   1 = depression (CES-D, 0-7, 8 categories)
#   2 = transport mobility (binary, 0-1, 2 categories)
#   3 = mobility limitations (0-3, 4 categories)
mediator_meta <- list(
  list(name = "Depression (CES-D)",       var = "depression_lag", categories = 0:7),
  list(name = "Transport mobility",       var = "transport_lag",  categories = 0:1),
  list(name = "Mobility limitations",     var = "mobility_lag",   categories = 0:3)
)

# ---- Load fitted results -----------------------------------------------------
full_result <- readRDS("../results/bjcm_full_results.rds")
chains <- full_result$results$results$chains
med_data <- full_result$lagged_data$analysis_data

stopifnot(length(chains[[1]]$y_pred_probs_mediators) == length(mediator_meta))

# ---- Loop over mediators: observed vs. posterior predictive -----------------
all_comparisons <- list()
all_plot_data    <- list()

for (k in seq_along(mediator_meta)) {
  
  meta       <- mediator_meta[[k]]
  categories <- meta$categories
  
  observed_raw <- med_data[[meta$var]]
  if (!all(sort(unique(observed_raw)) == categories)) {
    stop("Unexpected categories for ", meta$var, ": ",
         paste(sort(unique(observed_raw)), collapse = ", "),
         " -- check coding before proceeding.")
  }
  
  observed_prop <- as.numeric(table(factor(observed_raw, levels = categories))) /
    length(observed_raw)
  
  pred_probs_pooled <- do.call(rbind, lapply(chains, function(ch) {
    ch$y_pred_probs_mediators[[k]]
  }))
  
  if (ncol(pred_probs_pooled) != length(categories)) {
    stop("y_pred_probs_mediators[[", k, "]] has ", ncol(pred_probs_pooled),
         " columns; expected ", length(categories), " (one per category).")
  }
  
  predicted_mean <- colMeans(pred_probs_pooled)
  predicted_lo   <- apply(pred_probs_pooled, 2, quantile, probs = 0.025)
  predicted_hi   <- apply(pred_probs_pooled, 2, quantile, probs = 0.975)
  
  diff_pp <- round(100 * (observed_prop - predicted_mean), 2)
  
  comparison <- data.frame(
    mediator  = meta$name,
    category  = categories,
    observed  = round(100 * observed_prop, 2),
    predicted = round(100 * predicted_mean, 2),
    lo        = round(100 * predicted_lo, 2),
    hi        = round(100 * predicted_hi, 2),
    diff_pp   = diff_pp
  )
  all_comparisons[[k]] <- comparison
  
  cat("============ POSTERIOR PREDICTIVE CHECK:", meta$name, "============\n")
  print(comparison[, -1], row.names = FALSE)
  cat("\nMax |difference| (percentage points):", max(abs(diff_pp)),
      "at category", categories[which.max(abs(diff_pp))], "\n\n")
  
  all_plot_data[[k]] <- data.frame(
    mediator  = meta$name,
    category  = factor(categories),
    observed  = observed_prop,
    predicted = predicted_mean,
    lo        = predicted_lo,
    hi        = predicted_hi
  )
}

comparison_all <- do.call(rbind, all_comparisons)
plot_data_all  <- do.call(rbind, all_plot_data)
plot_data_all$mediator <- factor(plot_data_all$mediator,
                                 levels = sapply(mediator_meta, `[[`, "name"))

# ---- Save comparison table ---------------------------------------------------
write.csv(comparison_all, "../results/mediator_ppc_table.csv", row.names = FALSE)
cat("Saved mediator_ppc_table.csv\n")

# ---- Figure: three-panel mediator PPC, same visual style as Figure B.2 -----
p <- ggplot(plot_data_all, aes(x = category)) +
  geom_col(aes(y = observed), fill = "grey70", colour = "black",
           width = 0.65, linewidth = 0.3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15,
                colour = "#D2691E", linewidth = 0.6) +
  geom_point(aes(y = predicted), colour = "#D2691E", size = 2.6) +
  facet_wrap(~ mediator, scales = "free", nrow = 1) +
  labs(x = "Mediator category", y = "Proportion") +
  theme_bw(base_size = 12, base_family = "serif") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave("../results/mediator_ppc_check.pdf", p, width = 10, height = 3.8, device = "pdf")
cat("\nSaved mediator_ppc_check.pdf\n")

# ---- Overall summary for text --------------------------------------------
cat("\n============ OVERALL SUMMARY ============\n")
for (k in seq_along(mediator_meta)) {
  cc <- all_comparisons[[k]]
  cat(sprintf("%-22s max |diff| = %.2f pp (category %s)\n",
              mediator_meta[[k]]$name, max(abs(cc$diff_pp)),
              cc$category[which.max(abs(cc$diff_pp))]))
}