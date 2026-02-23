# R Code for Longitudinal Causal Estimation Methods

This repository provides R scripts implementing the simulation design and estimation procedures described in the accompanying manuscript.

The code reproduces a single simulation run under a longitudinal data-generating process (DGP) and implements multiple causal estimation approaches for time-varying exposures.

An R package version of these methods is currently under development.

---

## Repository Structure

- **`01_single_run_simulation.R`**  
  Generates one simulated dataset under the specified DGP and runs all estimation methods.

- **`02_methods_functions.R`**  
  Contains implementation of the following methods:
  
  - Inverse Probability of Treatment Weighting (IPTW) via CBPS  
  - Basic G-estimation (NUC-based, residualizing exposure only)  
  - Efficient G-estimation (double residualization; two-step GMM)  
  - IV-based decay model estimation (two-parameter decay model)

---

## Implemented Methods

Includes implementations of:

1. **IPTW (Marginal Structural Model)**  
   - Propensity score estimation via the `CBPS` package  
   - Weighted regression  
   - Sandwich (HC1) variance estimation  
   - Weight diagnostics

2. **Basic G-estimation (Structural Nested Mean Model)**  
   - Residualization of exposure on observed history  
   - Just-identified GMM estimation  
   - Sandwich variance estimation

3. **Efficient G-estimation**  
   - Double residualization (exposure and outcome)  
   - Two-step efficient GMM  
   - Sandwich variance estimation

4. **IV-based Decay Model**  
   - Two-parameter decay structure:  
     - `beta` (contemporaneous effect)  
     - `alpha` (decay parameter)  
   - Sandwich variance estimation  
   - Implied total longitudinal effects derived from `(beta, alpha)`

---

## Requirements

The following R packages are required:

- `MASS`
- `CBPS`
- `sandwich`
- `lmtest`
- `AER`

Install missing packages via:

```r
install.packages(c("MASS", "CBPS", "sandwich", "lmtest", "AER"))
```
