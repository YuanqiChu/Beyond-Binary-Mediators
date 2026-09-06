# =============================================================================
# 06a_ordinal_scale_profile_table3.R
# =============================================================================
# Table 3 / Figure 3 (ordinal scale profile): population-averaged predicted
# probability of each loneliness category (3-9, i.e. 7 ordinal categories)
# under the counterfactual scenarios of universal living alone (T=1) versus
# universal living with others (T=0), holding all other covariates at their
# observed values, computed directly from the fitted BMEOP posterior via
# standard posterior-predictive counterfactual simulation.
#
# For each retained posterior draw l and each individual i:
#   eta_i(t) = X_i^T beta^(l) [with treatment column set to t] + U_w(i)^(l)
#   P(Y_i = k | eta_i(t)) = Phi(S_k^(l) - eta_i(t)) - Phi(S_{k-1}^(l) - eta_i(t))
# Averaged first over individuals (population-averaged, at fixed l), then
# over posterior draws l (across all four chains).
#
# Run from MH_BMEOP/scripts/.
# =============================================================================

library(dplyr)

mh_bmeop <- readRDS("../results/bmeop_full_results.rds")
elsa_long <- read.csv("../data/elsa_longitudinal_analysis.csv")

covariate_names <- rownames(mh_bmeop$beta_summary)
treatment_idx   <- which(covariate_names == "livalone_TREATMENT")

# ---- Reconstruct the exact design matrix used for fitting ------------------
required_vars <- c("idauniq", "wave", "loneliness", "livalone", "age_gr", "dhsex2",
                   "edqual2", "self_reported_health", "sclife", "depression",
                   "transport_mobility", "mobility_limitations")
bmeop_complete <- elsa_long[complete.cases(elsa_long[, required_vars]), ]

build_bmeop_X <- function(d) {
  X_components <- list(
    intercept = rep(1, nrow(d)),
    livalone_TREATMENT = d$livalone
  )
  if ("age_gr" %in% names(d)) X_components$age_gr <- d$age_gr
  if ("dhsex2" %in% names(d)) X_components$dhsex2 <- d$dhsex2
  if ("edqual2" %in% names(d)) X_components$edqual2 <- d$edqual2
  if ("self_reported_health" %in% names(d)) X_components$health <- d$self_reported_health
  if ("sclife" %in% names(d)) X_components$sclife <- d$sclife
  if ("depression" %in% names(d)) X_components$depression <- d$depression
  if ("transport_mobility" %in% names(d)) X_components$transport_mobility <- d$transport_mobility
  if ("mobility_limitations" %in% names(d)) X_components$mobility_limitations <- d$mobility_limitations
  X <- do.call(cbind, X_components)
  colnames(X) <- names(X_components)
  X
}

X <- build_bmeop_X(bmeop_complete)
stopifnot(identical(colnames(X), covariate_names))

X_t1 <- X; X_t1[, treatment_idx] <- 1
X_t0 <- X; X_t0[, treatment_idx] <- 0

unique_waves <- sort(unique(bmeop_complete$wave))
wave_mapped  <- match(bmeop_complete$wave, unique_waves)
n_waves      <- length(unique_waves)
n            <- nrow(X)
K            <- 7  # loneliness categories 3-9

cat("n =", n, "| n_waves =", n_waves, "| K =", K, "\n\n")

# ---- Accumulate category-probability sums across all posterior draws -------
chains <- mh_bmeop$chains
n_chains <- length(chains)

prob_sum_t1 <- numeric(K)
prob_sum_t0 <- numeric(K)
n_draws_total <- 0

for (ch in chains) {
  n_iter <- nrow(ch$beta)
  for (l in 1:n_iter) {
    beta_l <- ch$beta[l, ]
    U_l    <- ch$random_effects[l, ]
    S_l    <- c(-Inf, 0, ch$thresholds[l, 2:6], Inf)  # S_1 = 0 fixed; length K+1
    
    eta_t1 <- as.vector(X_t1 %*% beta_l) + U_l[wave_mapped]
    eta_t0 <- as.vector(X_t0 %*% beta_l) + U_l[wave_mapped]
    
    probs_t1 <- matrix(NA, n, K)
    probs_t0 <- matrix(NA, n, K)
    for (k in 1:K) {
      probs_t1[, k] <- pnorm(S_l[k + 1] - eta_t1) - pnorm(S_l[k] - eta_t1)
      probs_t0[, k] <- pnorm(S_l[k + 1] - eta_t0) - pnorm(S_l[k] - eta_t0)
    }
    
    prob_sum_t1 <- prob_sum_t1 + colMeans(probs_t1)
    prob_sum_t0 <- prob_sum_t0 + colMeans(probs_t0)
    n_draws_total <- n_draws_total + 1
  }
}

prob_mean_t1 <- prob_sum_t1 / n_draws_total
prob_mean_t0 <- prob_sum_t0 / n_draws_total

cat("Total posterior draws pooled across", n_chains, "chains:", n_draws_total, "\n\n")

# ---- Report ------------------------------------------------------------------
categories <- 3:9
results_df <- data.frame(
  category          = categories,
  living_with_others = round(100 * prob_mean_t0, 1),
  living_alone        = round(100 * prob_mean_t1, 1),
  difference_pp       = round(100 * (prob_mean_t1 - prob_mean_t0), 1)
)

cat("=============================================================\n")
cat("Population-averaged predicted probability by loneliness category\n")
cat("=============================================================\n\n")
print(results_df, row.names = FALSE)

cat("\nSanity check -- each column should sum to ~100%:\n")
cat("Living with others total:", sum(results_df$living_with_others), "%\n")
cat("Living alone total:      ", sum(results_df$living_alone), "%\n")

saveRDS(results_df, "../results/ordinal_scale_profile.rds")
cat("\nSaved: ordinal_scale_profile.rds\n")