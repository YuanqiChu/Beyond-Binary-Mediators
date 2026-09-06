# =============================================================================
# Figure: natural indirect effects by mediator (NIE forest plot)
# =============================================================================
# Regenerates the NIE forest plot (nie_forest.pdf / .png): the Shapley-
# symmetrised path-specific NIE for each mediator, with 95% credible
# intervals, computed directly from the saved fit.
#
# Run from scripts/, so that ../results/ resolves correctly.
# =============================================================================

library(ggplot2)

has_ggtext <- requireNamespace("ggtext", quietly = TRUE)
if (has_ggtext) library(ggtext)

full <- readRDS("../results/bjcm_full_results.rds")
res  <- full$results$results

mediator_names <- c("Depression", "Transport Mobility", "Mobility Limitations")
nie_indiv <- res$indirect_effect_individual_summaries

df <- data.frame(
  mediator = c(mediator_names, "Total NIE"),
  estimate = c(sapply(nie_indiv, function(m) m["mean"]),   res$indirect_effect_summary["mean"]),
  lo       = c(sapply(nie_indiv, function(m) m["q2.5"]),   res$indirect_effect_summary["q2.5"]),
  hi       = c(sapply(nie_indiv, function(m) m["q97.5"]),  res$indirect_effect_summary["q97.5"])
)
# Order to match the manuscript figure (smallest contributor at top)
df <- df[match(c("Mobility Limitations", "Transport Mobility", "Depression", "Total NIE"), df$mediator), ]

if (has_ggtext) {
  df$mediator_label <- ifelse(df$mediator == "Total NIE", "**Total NIE**", df$mediator)
} else {
  df$mediator_label <- df$mediator
}
df$mediator_label <- factor(df$mediator_label, levels = rev(df$mediator_label))

cols <- setNames(
  c("#7B2D8E", "#1F7A72", "#1A3FA0", "#B32020"),
  df$mediator_label
)

p <- ggplot(df, aes(x = estimate, y = mediator_label, colour = mediator_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.6) +
  geom_point(size = 3) +
  scale_colour_manual(values = cols, guide = "none") +
  labs(x = "Effect size (loneliness scale points)", y = NULL) +
  theme_bw(base_size = 13, base_family = "serif") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.ticks.y = element_blank(),
    axis.text.y = if (has_ggtext) element_markdown(size = 13) else element_text(size = 13)
  )

ggsave("../results/nie_forest.pdf", p, width = 6.5, height = 4.0, device = "pdf")
ggsave("../results/nie_forest.png", p, width = 6.5, height = 4.0, dpi = 200)

cat("Saved nie_forest.pdf / nie_forest.png\n")
