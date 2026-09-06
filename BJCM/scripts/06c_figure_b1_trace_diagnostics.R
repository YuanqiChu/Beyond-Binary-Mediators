# =============================================================================
# 06c_figure_b1_trace_diagnostics.R
#
# Regenerates Appendix Figure B.1 -- trace plots for TE, NDE, and NIE
# across the four post-burn-in MCMC chains -- directly from the saved
# BJCM result object.
#
# Nothing here re-fits any model. It only loads the saved result object
# and plots the stored per-iteration chain values, so it runs in seconds.
#
# Field names (chains[[i]]$total_effects / direct_effects /
# indirect_effects, and the rowMeans-if-matrix convention) match those
# used in 06b_manuscript_numbers.R and 04a_convergence_diagnostics.R for
# this same result object.
#
# Run from scripts/, so that ../results/ resolves correctly.
# =============================================================================

library(ggplot2)

mh_result <- readRDS("../results/bjcm_full_results.rds")
chains <- mh_result$results$results$chains
n_chains <- length(chains)

# Per-iteration, person-averaged TE/NDE/NIE for one chain. Mirrors the
# person_avg() helper in bjcm_manuscript_numbers.R: if the stored object
# is a matrix (iterations x individuals), average across individuals to
# get one value per iteration; if it is already a vector of per-iteration
# values, use it as-is.
person_avg <- function(x) if (is.matrix(x)) rowMeans(x) else as.numeric(x)

extract_estimand <- function(chains, field, label) {
  do.call(rbind, lapply(seq_along(chains), function(chain_i) {
    v <- person_avg(chains[[chain_i]][[field]])
    data.frame(
      iteration = seq_along(v),
      value     = v,
      chain     = factor(chain_i),
      estimand  = label
    )
  }))
}

df <- rbind(
  extract_estimand(chains, "total_effects",    "Total Effect (TE)"),
  extract_estimand(chains, "direct_effects",   "Natural Direct Effect (NDE)"),
  extract_estimand(chains, "indirect_effects", "Natural Indirect Effect (NIE)")
)

df$estimand <- factor(df$estimand,
                      levels = c("Total Effect (TE)",
                                 "Natural Direct Effect (NDE)",
                                 "Natural Indirect Effect (NIE)"))

cat("Chains:", n_chains, "\n")
cat("Iterations per chain (TE):",
    sum(df$estimand == "Total Effect (TE)" & df$chain == levels(df$chain)[1]),
    "\n")

chain_cols <- c("#1A3FA0", "#B32020", "#1F7A72", "#7B2D8E")[seq_len(n_chains)]

p <- ggplot(df, aes(x = iteration, y = value, colour = chain)) +
  geom_line(linewidth = 0.25, alpha = 0.85) +
  facet_wrap(~estimand, ncol = 1, scales = "free_y") +
  scale_colour_manual(values = chain_cols, name = "Chain") +
  labs(x = "Post-burn-in iteration", y = "Estimate") +
  theme_bw(base_size = 13, base_family = "serif") +
  theme(
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", hjust = 0),
    legend.position = "top"
  )

ggsave("../results/trace_diagnostics.png", p, width = 6.5, height = 7.5, dpi = 200)

cat("Saved trace_diagnostics.png\n")