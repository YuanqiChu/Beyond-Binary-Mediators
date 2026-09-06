# =============================================================================
# Marginal-likelihood Metropolis-within-Gibbs threshold update
# =============================================================================
# Implements Step 3 of Algorithm A.1 (and, via substitution, Steps 3 and 7
# of Algorithm A.2): a Metropolis-within-Gibbs update for the ordered-probit
# thresholds, targeting the observed-data marginal likelihood
#   P(Y_i = k | eta_i, S) = Phi(S_k - eta_i) - Phi(S_{k-1} - eta_i),
#   eta_i = X_i^T beta + U_{w(i)}
# with the latent variable Z integrated out of the conditional, rather than
# conditioning on the current draw of Z_i. This avoids the poor mixing that
# direct conditioning on Z produces in ordered-probit Gibbs samplers, where
# the resulting conditional interval for S can be pathologically narrow
# (Cowles, 1996).
#
# Proposal: truncated random walk
#   S_j* = S_j + eps,  eps ~ N(0, proposal_sd_j^2)
#   truncated to (S_{j-1} + g, S_{j+1} - g)
# This is NOT symmetric in probability (the truncated-normal normalising
# constant depends on S_j itself, since the truncation bounds are fixed but
# the density is centred at S_j), so the full Hastings ratio -- including
# q(S_j | S_j*) / q(S_j* | S_j) -- is required and computed below.
#
# Only observations with Y_i in {j, j+1} are affected by S_j, so the
# log-likelihood difference is computed efficiently over those two groups
# only, not the full dataset.
#
# S_1 = 0 remains fixed for identification, matching Algorithm A.1; only
# S_2, ..., S_{K-1} are updated.
#
# Acceptance rate logging is FILE-BASED (Sys.getpid()-keyed), because chains
# run in separate forked processes and an in-memory environment's contents
# would not propagate back to the parent process.
# =============================================================================

.mh_log_dir <- "../logs/mh_threshold_logs"

reset_mh_threshold_log <- function(log_dir = .mh_log_dir) {
  if (dir.exists(log_dir)) unlink(log_dir, recursive = TRUE)
  dir.create(log_dir, showWarnings = FALSE)
  invisible(NULL)
}

log_mh_attempt <- function(label, j, accepted, log_dir = .mh_log_dir) {
  if (!dir.exists(log_dir)) dir.create(log_dir, showWarnings = FALSE)
  log_file <- file.path(log_dir, paste0("pid_", Sys.getpid(), ".csv"))
  cat(paste(label, j, as.integer(accepted), sep = ","), "\n",
      file = log_file, append = TRUE, sep = "")
  invisible(NULL)
}

get_mh_acceptance_summary <- function(log_dir = .mh_log_dir) {
  if (!dir.exists(log_dir)) {
    warning("Log directory '", log_dir, "' does not exist -- did you call ",
            "reset_mh_threshold_log() before running the model?")
    return(data.frame(threshold = character(0), proposed = integer(0),
                      accepted = integer(0), rate = numeric(0)))
  }
  log_files <- list.files(log_dir, pattern = "^pid_.*\\.csv$", full.names = TRUE)
  if (length(log_files) == 0) {
    warning("No MH threshold log files found in '", log_dir, "'.")
    return(data.frame(threshold = character(0), proposed = integer(0),
                      accepted = integer(0), rate = numeric(0)))
  }
  all_rows <- do.call(rbind, lapply(log_files, function(f) {
    lines <- readLines(f, warn = FALSE)
    if (length(lines) == 0) return(NULL)
    parts <- strsplit(lines, ",", fixed = TRUE)
    do.call(rbind, lapply(parts, function(p) {
      data.frame(threshold = paste0(p[1], "_j", p[2]), accepted = as.integer(p[3]))
    }))
  }))
  if (is.null(all_rows) || nrow(all_rows) == 0) {
    warning("Log files found but contained no logged attempts.")
    return(data.frame(threshold = character(0), proposed = integer(0),
                      accepted = integer(0), rate = numeric(0)))
  }
  agg <- aggregate(accepted ~ threshold, data = all_rows,
                   FUN = function(x) c(proposed = length(x), accepted = sum(x)))
  result <- data.frame(
    threshold = agg$threshold,
    proposed  = agg$accepted[, "proposed"],
    accepted  = agg$accepted[, "accepted"]
  )
  result$rate <- result$accepted / result$proposed
  result[order(result$threshold), ]
}

# ---- Log prior ratio for the threshold ------------------------------------
# Uniform prior on the BOUNDED ordering region 0 = S_1 < S_2 < ... <
# S_{K-1} < S_max (S_max supplied by the caller; see S_max argument of
# update_thresholds_mh_marginal() below). This is a proper prior -- an
# unbounded flat prior over the full ordered simplex would not be -- and
# because every proposal is checked against (lower_bd, upper_bd) before
# reaching this function, both S_proposed and S_current always lie inside
# the bounded region, so the ratio of two uniform densities on that region
# is always 1 regardless of S_max's value.
log_prior_ratio_threshold <- function(S_proposed, S_current, j) {
  0  # uniform prior on the bounded ordering region: log(pi(S*)/pi(S)) = 0
}

# ---- Truncated-normal proposal density (log scale) -------------------------
# log q(x | centre, sd, lower, upper) for x ~ N(centre, sd^2) truncated to
# (lower, upper). Returns -Inf if x is outside (lower, upper).
log_dtruncnorm <- function(x, centre, sd, lower, upper) {
  if (x <= lower || x >= upper) return(-Inf)
  log_dens <- dnorm(x, mean = centre, sd = sd, log = TRUE)
  log_norm_const <- log(pnorm(upper, mean = centre, sd = sd) -
                          pnorm(lower, mean = centre, sd = sd))
  log_dens - log_norm_const
}

# ---- The marginal-likelihood MH threshold update ---------------------------
# Signature matches update_thresholds_conservative(), with two
# additions (eta, proposal_sd) needed for the marginal-likelihood evaluation
# and proposal step size -- both have to be supplied at the call sites
# (see the integration notes at the bottom of this file for exactly
# what changes at each call site in the model file that sources this).
#
#   z_latent            : unused here (kept in the signature only for
#                          call-site compatibility with
#                          function; marginal-likelihood MH does not
#                          condition on Z at all)
#   y_observed          : observed category vector (Y or M^(k))
#   current_thresholds  : current threshold vector, length K-1, S[1] = 0
#   iteration           : current Gibbs iteration (kept for interface
#                          compatibility; no longer used to gate update
#                          frequency -- see note below)
#   eta                 : linear predictor vector, eta_i = X_i^T beta +
#                          U_{w(i)} (or the mediator-specific/outcome-
#                          specific equivalent), same length as y_observed
#   label               : string identifying which threshold set this is,
#                          for acceptance-rate logging
#   proposal_sd         : vector of length K-1 (or a single scalar recycled),
#                          random-walk proposal SD per threshold. Tune this
#                          to target roughly 20%-50% acceptance -- do NOT
#                          aim for acceptance = 1 (high acceptance does not
#                          imply good exploration).
#   min_spacing         : minimum allowed gap between adjacent thresholds
#                          (g in the derivation); default 0.05.
#   update_every_iter   : if TRUE (default), update every Gibbs sweep, which
#                          is the recommended setting for a valid MH kernel.
#                          Set to FALSE and pass an `iteration` divisible
#                          check yourself only if a different update
#                          cadence is needed.
#   S_max               : fixed finite upper bound on the last movable
#                          threshold, S_{K-1} < S_max. Without this the
#                          proposal region for the last threshold is
#                          unbounded and the implied uniform prior on the
#                          ordering-constrained region is improper. S_max
#                          should be chosen well outside the range the
#                          posterior ever visits (checked empirically via
#                          the trace plots / threshold posteriors) so that
#                          it has no material effect on the fitted model;
#                          default of 10 is many standard deviations beyond
#                          any threshold seen in this application (all
#                          fitted thresholds lie in roughly [0, 3]).
update_thresholds_mh_marginal <- function(z_latent, y_observed, current_thresholds,
                                          iteration, eta, label = "threshold",
                                          proposal_sd = 0.08, min_spacing = 0.05,
                                          update_every_iter = TRUE, S_max = 10) {
  K <- max(y_observed, na.rm = TRUE)
  if (K <= 2) {
    return(current_thresholds)
  }
  if (!update_every_iter && iteration %% 50 != 0) {
    return(current_thresholds)
  }
  
  new_thresholds <- current_thresholds
  n_thresh <- length(current_thresholds)  # = K - 1
  
  if (length(proposal_sd) == 1) {
    proposal_sd <- rep(proposal_sd, n_thresh)
  }
  
  for (j in 2:n_thresh) {
    
    S_j      <- new_thresholds[j]
    S_jm1    <- new_thresholds[j - 1]
    S_jp1    <- if (j < n_thresh) new_thresholds[j + 1] else S_max
    
    lower_bd <- S_jm1 + min_spacing
    upper_bd <- S_jp1 - min_spacing
    
    if (!(lower_bd < upper_bd)) {
      # No feasible region this sweep (extremely tight neighbouring
      # thresholds) -- skip and log as not attempted/rejected.
      log_mh_attempt(label, j, accepted = FALSE)
      next
    }
    
    omega_j <- proposal_sd[j]
    
    # ---- Propose: truncated-normal random walk ----
    S_j_star <- rnorm(1, mean = S_j, sd = omega_j)
    if (S_j_star <= lower_bd || S_j_star >= upper_bd) {
      # Rejected immediately: proposal fell outside the ordering-constrained
      # region. This is a legitimate MH rejection (proposal density is zero
      # there), not a bug -- log it as such.
      log_mh_attempt(label, j, accepted = FALSE)
      next
    }
    
    # ---- Observed-data log-likelihood difference (only Y_i in {j, j+1}) ----
    idx_j   <- which(y_observed == j)
    idx_jp1 <- which(y_observed == (j + 1))
    
    ll_current <- 0
    ll_proposed <- 0
    
    if (length(idx_j) > 0) {
      eta_j <- eta[idx_j]
      p_current  <- pnorm(S_j      - eta_j) - pnorm(S_jm1 - eta_j)
      p_proposed <- pnorm(S_j_star - eta_j) - pnorm(S_jm1 - eta_j)
      # Guard against numerical zero (can happen in extreme tails); such
      # observations contribute a large negative log-likelihood rather than
      # -Inf/NaN, which would otherwise silently corrupt the sum.
      p_current  <- pmax(p_current,  .Machine$double.eps)
      p_proposed <- pmax(p_proposed, .Machine$double.eps)
      ll_current  <- ll_current  + sum(log(p_current))
      ll_proposed <- ll_proposed + sum(log(p_proposed))
    }
    
    if (length(idx_jp1) > 0) {
      eta_jp1 <- eta[idx_jp1]
      p_current  <- pnorm(S_jp1 - eta_jp1) - pnorm(S_j      - eta_jp1)
      p_proposed <- pnorm(S_jp1 - eta_jp1) - pnorm(S_j_star - eta_jp1)
      p_current  <- pmax(p_current,  .Machine$double.eps)
      p_proposed <- pmax(p_proposed, .Machine$double.eps)
      ll_current  <- ll_current  + sum(log(p_current))
      ll_proposed <- ll_proposed + sum(log(p_proposed))
    }
    
    delta_ll <- ll_proposed - ll_current
    
    # ---- Prior ratio (see log_prior_ratio_threshold() above) ----
    delta_log_prior <- log_prior_ratio_threshold(S_j_star, S_j, j)
    
    # ---- Hastings correction: truncated-normal proposal is NOT symmetric,
    #      since the truncation bounds (lower_bd, upper_bd) are fixed but
    #      the proposal density is centred at the current/proposed value,
    #      so the normalising constant differs between the two directions ----
    log_q_forward  <- log_dtruncnorm(S_j_star, centre = S_j,      sd = omega_j,
                                     lower = lower_bd, upper = upper_bd)
    log_q_backward <- log_dtruncnorm(S_j,      centre = S_j_star, sd = omega_j,
                                     lower = lower_bd, upper = upper_bd)
    
    log_alpha <- delta_ll + delta_log_prior + (log_q_backward - log_q_forward)
    
    if (log(runif(1)) < log_alpha) {
      new_thresholds[j] <- S_j_star
      log_mh_attempt(label, j, accepted = TRUE)
    } else {
      log_mh_attempt(label, j, accepted = FALSE)
    }
  }
  
  new_thresholds[1] <- 0  # tau_1 stays fixed at 0 for identification
  return(new_thresholds)
}