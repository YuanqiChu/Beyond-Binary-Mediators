# BMEOP — Bayesian Mixed-Effects Ordered Probit model

Implements Algorithm 1 (Appendix A.1) of the manuscript: the single-mediator
average-treatment-effect model for living alone → loneliness. All scripts
assume the working directory is `scripts/` and use relative paths
(`../data/`, `../results/`, `../logs/`).

**Never source anything in this folder in the same R session as anything in
`../BJCM/`** — the two models each define six identically-named helper
functions with different implementations, and sourcing both overwrites one
set with the other.

## Folder structure

```
BMEOP/
├── data/     elsa_longitudinal_analysis.csv  (not included — see root README)
├── scripts/  all R code, numbered by pipeline stage (below)
├── results/  saved fit, tables, and figures
└── logs/     background-run progress logs and MH acceptance logs
```

## Script manifest, in run order

| # | Script | Purpose |
|---|--------|---------|
| — | `00_background_run_utils.R` | Shared helpers (`run_in_background()`, `monitor_progress()`, `collect_result()`) sourced by the run scripts below. |
| 01 | `01_threshold_update_mh.R` | Marginal-likelihood Metropolis-within-Gibbs threshold update (Appendix A.1, Step 3), shared, identical copy also in `../BJCM/scripts/`. |
| 02 | `02_bmeop_model.R` | The model: Algorithm 1. |
| 03a | `03a_run_reduced_check.R` | Reduced-scale sanity check of MH acceptance rates before committing to the full run. |
| 03b | `03b_proposal_sd_grid_search.R` | Grid search over the threshold-update proposal SD, extrapolated to the full n=50,790 sample, to choose `proposal_sd = 0.015`. |
| 03c | `03c_run_full_model.R` | The full-scale run (4 chains). Saves the fit to `results/bmeop_full_results.rds`. |
| 04a | `04a_convergence_diagnostics.R` | Post-hoc R-hat/ESS diagnostics for the thresholds and the headline ATE estimand. |
| 04b | `04b_rhat_diagnostics_table_b2.R` | Standalone multi-chain Gelman–Rubin/ESS extraction, laid out per Appendix Table B.2. |
| 04c | `04c_table_b2_remaining_coefficients.R` | Completes Appendix Table B.2 with the covariates, `sigma2_u`, and thresholds not already covered by 04d. |
| 04d | `04d_diagnostics_summary.R` | Coefficient/ATE summary and the intercept / wave-level random-effect non-identifiability check (Section 4.2). |
| 06a | `06a_ordinal_scale_profile_table3.R` | Manuscript Table 3 / Figure 3: population-averaged predicted category probabilities under T=1 vs. T=0. |
| 06b | `06b_manuscript_numbers.R` | The remaining Table 3 numbers not already printed by the model's own console output (identified intercept + mean(U_w), sigma2_u). |

## Notes

BMEOP shows a near-exact additive non-identifiability between the raw
intercept and the six wave-level random effects `U_w` (r ≈ −0.999); the
*identified* combination `Intercept + mean(U_w)` is what converges cleanly
and is what's reported (see `04d_diagnostics_summary.R`).

## Results

- `bmeop_full_results.rds` — the fitted model (4 chains).
- `ordinal_scale_profile.rds` — Table 3 summary object.
