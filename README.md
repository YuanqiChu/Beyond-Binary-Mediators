# Beyond Binary Mediators

Bayesian causal mediation analysis for ordinal outcomes in longitudinal settings.

## Models

- **BMEOP** (Algorithm 1): Bayesian Mixed-Effects Ordered Probit
- **BJCM** (Algorithm 2): Bayesian Joint Causal Mediation with multiple ordinal mediators

## Application

Analysis of ELSA data (Waves 2-7, n=42,185) examining how living alone affects loneliness through depression, transport mobility, and physical mobility pathways.

## Files

- `bmeop_model.R` - Algorithm 1 implementation
- `bmeop_diagnost.R` - Algorithm 1 diagnostics
- `bmeop_usage_guide.R` - Algorithm 1 usage examples
- `bjcm_implementation.R` - Algorithm 2 implementation  

## Citation

Chu et al. (2025). "Beyond Binary Mediators: A Bayesian Mixed-Effect Modelling Framework for Understanding Causal Pathways to Loneliness in Late Life." 

## Requirements

R ≥ 4.0.0 with packages: MASS, mvtnorm, coda, Matrix, dplyr, ggplot2
