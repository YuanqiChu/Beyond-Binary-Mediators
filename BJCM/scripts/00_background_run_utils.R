# =============================================================================
# SHARED BACKGROUND-RUN / PROGRESS-MONITORING UTILITIES
# =============================================================================
# monitor_progress(), run_in_background(), and collect_result() are generic -
# they don't touch anything model-specific - so they live here instead of
# inside 02_bmeop_model.R (MH_BMEOP) or 02_bjcm_model.R (this folder). That keeps them available to
# 03a_run_reduced_check.R, 03b_run_full_model.R, and the diagnostics/table scripts in this folder
# without any of those three ever needing to source() both model files
# together (see 03a_run_reduced_check.R in MH_BMEOP for why that's avoided).
#
# Source this file at the top of any script that calls run_in_background()
# or monitor_progress() or collect_result() - the run scripts in this folder do this.
# =============================================================================


# =============================================================================
# PROGRESS MONITORING (for long-running mclapply chains)
# =============================================================================
# bmeop_model.R and bjcm_implementation.R now write per-chain progress to
# bmeop_progress_chain_<id>.log / bjcm_progress_chain_<id>.log in the working
# directory instead of printing to console, because cat() inside a forked
# mclapply worker does not reliably reach the parent session's console.
#
# USAGE: while run_bmeop_full_refit()/run_full_refit() is running (blocking)
# in one R session/terminal, open a SECOND R session or terminal in the SAME
# working directory, source this file, and call:
#   monitor_progress("bmeop")   # or monitor_progress("bjcm")
# Re-run the call any time to get an updated snapshot; it does not loop or
# block on its own.

monitor_progress <- function(model = c("bmeop", "bjcm"), n_chains = 4) {
  model <- match.arg(model)
  prefix <- paste0(model, "_progress_chain_")
  
  found_any <- FALSE
  for (i in seq_len(n_chains)) {
    fpath <- paste0("../logs/", prefix, i, ".log")
    if (file.exists(fpath)) {
      found_any <- TRUE
      cat(readLines(fpath), "\n")
    } else {
      cat("chain=", i, " - no progress file yet (hasn't reached first checkpoint, ",
          "or hasn't started)\n", sep = "")
    }
  }
  
  if (!found_any) {
    cat("\nNo progress files found in the working directory. Check you're in\n")
    cat("the same working directory as the R session running the refit.\n")
  }
  
  invisible(NULL)
}

# =============================================================================
# RUN WITHOUT BLOCKING THE CONSOLE (no second R session needed)
# =============================================================================
# run_bmeop_full_refit()/run_full_refit() block the console via mclapply
# until all chains finish. run_in_background() wraps either call in
# parallel::mcparallel(), which forks the WHOLE call into one background
# job and returns control to your console immediately - so you can call
# monitor_progress() in the SAME session while it runs, then collect_result()
# once it's done. (Not available on Windows - mcparallel is Unix/macOS only,
# same restriction as mclapply itself.)
#
# USAGE:
#   job <- run_in_background(run_bmeop_full_refit(bmeop_complete, n_chains = 4))
#   monitor_progress("bmeop")          # check anytime, console stays free
#   monitor_progress("bmeop")          # check again later, etc.
#   bmeop_full <- collect_result(job)  # blocks ONLY once you're ready to wait
#                                       # for it to actually finish; returns
#                                       # immediately if it's already done

run_in_background <- function(expr) {
  job <- eval(substitute(parallel::mcparallel(expr)), envir = parent.frame())
  cat("Job launched in background (PID visible via job$pid). Console is free.\n")
  cat("Use monitor_progress() to check status, collect_result(job) to retrieve the result.\n")
  job
}

collect_result <- function(job, wait = TRUE) {
  result <- parallel::mccollect(job, wait = wait)
  if (is.null(result) && !wait) {
    cat("Not finished yet - call collect_result(job) again later, or with wait=TRUE to block until done.\n")
    return(invisible(NULL))
  }
  result[[as.character(job$pid)]]
}
