# Reproducibility notes

## Original workflow dependencies

The original scripts were uploaded as separate files, but behave like one interactive R session. Key implicit dependencies include:

- `2.imputation.R` expects `df` from `1.read_data.R`
- `3.VAR_Testing.r` expects `df` and `date_column`
- `4. lars.R` expects `df` and `date_column`
- `5. james_stein.R` depends on `lars_transformations()` and `generate_graph_lars()`
- `6. simone.R` depends on `empty.graph()` from `bnlearn`, which is loaded earlier
- `8. Causal ML Algorithms.R` expects `df_full_graph` and `kb`
- `10.model_averaging.R` expects `graph_lars`, `graph_lar`, `graph_js`, `graph_simone`, `kb`, `df_full_graph`, and `dplyr`
- `12.causal_effect.R` expects `graph_js`, `df_full_graph`, `bn.fit()`, and `convert_bn_to_igraph()`

`7. parametrisation_processing.R` includes both:

1. creation of `df_full_graph`, and
2. discretisation of numeric variables for BN parameterisation.

But the causal-ML structure-learning script (`8. Causal ML Algorithms.R`) is meant to operate on the continuous/numeric version before discretisation. The code contains a comment telling the reader to 
run script 8 before finishing script 7.
