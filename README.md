# Beyond Binary Mediators

Code accompanying:

> Chu, Y., Yu, K., Rippon, I., Victor, C. (2025). *Beyond Binary Mediators:
> A Bayesian Mixed-Effect Modelling Framework for Understanding Causal
> Pathways to Loneliness in Late Life.* Submitted to the *Journal of the
> Royal Statistical Society, Series A (JRSSA)*.

The paper fits two Bayesian ordinal models to English Longitudinal Study of
Ageing (ELSA) data to estimate the effect of living alone on loneliness in
later life:

- **BMEOP** — Bayesian Mixed-Effects Ordered Probit model (Algorithm 1,
  Appendix A.1): the single-mediator average treatment effect model.
- **BJCM** — Bayesian Joint Causal Mediation model (Algorithm 2, Appendix
  A.2): the full three-mediator (depression, transport mobility, mobility
  limitations) parallel mediation model, with a Shapley-value decomposition
  of the natural indirect effect across mediators (Eq. 6).

## Repository layout

```
Beyond-Binary-Mediators/
├── data_preparation/  — builds elsa_longitudinal_analysis.csv from raw ELSA files
├── BJCM/               — the BJCM pipeline (data, scripts, results, logs)
└── BMEOP/              — the BMEOP pipeline (data, scripts, results, logs)
```

**BJCM/ and BMEOP/ are independent.** Their model files each define six
identically-named helper functions (`rmvnorm_robust`, `rtruncnorm_robust`,
`safe_matrix_inverse`, `update_U_wave_robust`, `update_sigma2_U_robust`,
`update_thresholds_conservative`) with *different* implementations.
**Always run BJCM and BMEOP scripts in separate R sessions.** Each folder
is fully self-contained.

See `BJCM/README.md` and `BMEOP/README.md` for the run order and script
manifest within each pipeline.

## Data availability

ELSA is restricted-access microdata distributed by the UK Data Service
under an end-user licence that prohibits redistribution, so
`BJCM/data/*.csv` and `BMEOP/data/*.csv` are **not** included in this
repository. To reproduce the analysis:

1. Apply for ELSA access via the [UK Data Service](https://ukdataservice.ac.uk/)
   (study number SN 5050).
2. Run `data_preparation/build_elsa_longitudinal_analysis.R` against the raw
   wave files to build `elsa_longitudinal_analysis.csv` (see
   `data_preparation/README.md`).
3. Copy that file into both `BJCM/data/` and `BMEOP/data/`.

Posterior summaries, convergence diagnostics, and the manuscript's tables
and figures can all be reproduced from the saved fit objects in each
`results/` folder without re-running the sampler — see the per-folder
README for which script does what.

## Requirements

R (≥ 4.2 recommended) with base packages plus `ggplot2` for figures and
`coda` for Gelman–Rubin diagnostics (`ggtext` is optional, used only to
bold one label in the NIE forest plot).

## License

Code is released under the MIT License (see `LICENSE`).
