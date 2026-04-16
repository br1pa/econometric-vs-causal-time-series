# Econometric vs. Causal Structure-Learning for Time-Series Policy Decisions: UK COVID-19

This repository packages the code and data used for the paper:

> **Econometric vs. Causal Structure-Learning for Time-Series Policy Decisions: Evidence from the UK COVID-19 Policies**

The paper studies four econometric methods and eleven causal structure-learning algorithms on a UK COVID-19 time-series dataset. The attached manuscript describes a dataset with 46 columns and 866 daily observations spanning **2020-01-30 to 2022-06-13**. It uses Kalman-based imputation for missing values, k-means-based discretisation for DBN parameterisation, a knowledge-based benchmark DAG, model averaging across four econometric graphs, and causal-effect estimation using the do-operator.

## What is in this repo

- `data/raw/` – raw Excel input file
- `scripts/` – original uploaded scripts, preserved as-is

## Expected outputs

- DOT graph files for LASSO, LAR, James-Stein, and SIMONE
- the VAR lag-selection plot
- pairwise SHD comparison figure
- optional intermediate `.rds` objects for reuse

## Package dependencies

The scripts call or rely on:

`readxl`, `imputeTS`, `dplyr`, `vars`, `moments`, `tseries`, `lmtest`, `fGarch`, `lars`, `bnlearn`, `GeneNet`, `simone`, `parallel`, `igraph`, `ggplot2`, `causaleffect`, `tidyr` and `tibble`.
