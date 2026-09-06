# BJCM — Bayesian Joint Causal Mediation model

Implements Algorithm 2 (Appendix A.2) of the manuscript: the three-mediator
(depression, transport mobility, mobility limitations) parallel mediation
model, with treatment fixed at wave 2 and mediators lagged one wave ahead
of the outcome. All scripts assume the working directory is `scripts/` and
use relative paths (`../data/`, `../results/`, `../logs/`).

**Never source anything in this folder in the same R session as anything in
`../BMEOP/`** — the two models each define six identically-named helper
functions with different implementations, and sourcing both overwrites one
set with the other.

## Folder structure

```
BJCM/
├── data/     elsa_longitudinal_analysis.csv  (not included — see root README)
├── scripts/  all R code, numbered by pipeline stage (below)
├── results/  saved fit, tables, and figures
└── logs/     background-run progress logs and MH acceptance logs
```

## Script manifest, in run order

| # | Script | Purpose |
|---|--------|---------|
| — | `00_background_run_utils.R` | Shared helpers (`run_in_background()`, `monitor_progress()`, `collect_result()`) sourced by the run scripts below. |
| 01 | `01_threshold_update_mh.R` | Marginal-likelihood Metropolis-within-Gibbs threshold update (Appendix A.1, Step 3), shared, identical copy also in `../BMEOP/scripts/`. |
| 02 | `02_bjcm_model.R` | The model: Algorithm 2, with the Shapley-symmetrised NIE decomposition (Eq. 6) embedded in `simulate_parallel_mediation_effects()`. |
| 03a | `03a_run_reduced_check.R` | Reduced-scale (n=5,000, 4 chains, 2,000 iter) sanity check of MH acceptance rates before committing to the full run. |
| 03b | `03b_run_full_model.R` | The full-scale run (4 chains, several hours). Saves the fit to `results/bjcm_full_results.rds`. |
| 04a | `04a_convergence_diagnostics.R` | Post-hoc R-hat/ESS diagnostics for the thresholds and the headline NIE estimand. |
| 04b | `04b_rhat_diagnostics_standalone.R` | Standalone multi-chain Gelman–Rubin extraction, robust to near-singular parameter blocks. |
| 04c | `04c_table_b1_convergence.R` | Reproduces Appendix Table B.1 (R-hat for TE/NDE/NIE). |
| 04d | `04d_intercept_wave_identifiability_check.R` | Checks the intercept / wave-random-effect additive non-identifiability (Appendix B.2, Table B.3). |
| 04e | `04e_ridge_detail_appendix_b2.R` | ESS / correlation numbers behind the Appendix B.2 ridge discussion. |
| 04f | `04f_posterior_predictive_check.R` | Posterior-predictive check of the outcome category distribution (Appendix Figure B.2). |
| 05a | `05a_mediator_independence_check.R` | Checks the conditional-independence-across-mediators assumption via randomized quantile residuals (Section 5, Discussion). |
| 05b | `05b_mediator_correlation_sensitivity.R` | Sensitivity of the NIE decomposition to the residual correlation found in 05a. |
| 05c | `05c_mediator_posterior_predictive_check.R` | Posterior-predictive check of each mediator's category distribution, mirroring `04f`. |
| 06a | `06a_extract_final_results.R` | Extracts every number the manuscript reports for BJCM (TE/NDE/NIE with CrI, per-mediator paths, decomposition check, E-value) from the saved fit. |
| 06b | `06b_manuscript_numbers.R` | Same extraction, packaged as a reusable function. |
| 06c | `06c_figure_b1_trace_diagnostics.R` | Regenerates Appendix Figure B.1 (TE/NDE/NIE trace plots). |
| 06d | `06d_manuscript_figures.R` | Regenerates the NIE forest plot (Figure 3), computed directly from the saved fit. |

## Results

- `bjcm_full_results.rds` — the fitted model (4 chains).
- `nie_forest.pdf` / `.png`, `ppc_check.pdf`, `mediator_ppc_check.pdf`, `trace_diagnostics.png` — manuscript and appendix figures.
- `mediator_ppc_table.csv` — the observed-vs-predicted comparison table behind `mediator_ppc_check.pdf`.
