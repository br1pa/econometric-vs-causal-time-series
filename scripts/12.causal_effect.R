# 12.causal_effect.R

ds_require(c("structure_metrics_df", "ds_local_counts"))
policy_iss <- getOption("downstream.policy_iss", 10)
policy_max_unobserved_mass <- getOption("downstream.policy_max_unobserved_mass", 0.05)
policy_methods <- getOption("downstream.policy_methods", setdiff(ds_methods, "knowledge"))
policy_exposure_lags <- getOption("downstream.exposure_lags", 1L)
policy_include_max_lag <- getOption("downstream.include_max_lag_exposure", TRUE)

if (length(policy_iss) != 1L || !is.finite(policy_iss) || policy_iss <= 0 ||
    length(policy_max_unobserved_mass) != 1L || !is.finite(policy_max_unobserved_mass) ||
    policy_max_unobserved_mass < 0 || policy_max_unobserved_mass > 1) stop("Invalid settings.")

if (!length(policy_methods) || anyNA(policy_methods) || any(!policy_methods %in% ds_methods)) stop("Unknown or missing policy methods.")

policy_methods <- unique(policy_methods)
if (!is.numeric(policy_exposure_lags) || !length(policy_exposure_lags) ||
    anyNA(policy_exposure_lags) || any(!is.finite(policy_exposure_lags)) ||
    any(policy_exposure_lags < 1 | policy_exposure_lags != as.integer(policy_exposure_lags))) {
  stop("downstream.exposure_lags must contain positive integers.")
}

if (!is.logical(policy_include_max_lag) || length(policy_include_max_lag) != 1L ||
    is.na(policy_include_max_lag)) stop("downstream.include_max_lag_exposure must be TRUE/FALSE.")

policy_exposure_lags <- sort(unique(as.integer(policy_exposure_lags)))
if (any(policy_exposure_lags > max(candidate_lags))) stop("Requested exposure lag exceeds every saved model specification.")
outcome_nodes <- paste0(c("New.cases", "New.infections", "Reinfections"), "_t")
x_variables_base <- mobility_vars

policy_expected_sign <- getOption("downstream.expected_sign", NULL)
if (!is.null(policy_expected_sign)) {
  if (!is.numeric(policy_expected_sign) || is.null(names(policy_expected_sign)) ||
      anyDuplicated(names(policy_expected_sign)) ||
      any(!names(policy_expected_sign) %in% x_variables_base) ||
      anyNA(policy_expected_sign) || any(!policy_expected_sign %in% c(-1, 0, 1))) {
    stop("downstream.expected_sign must be a named numeric vector of -1, 0, +1 for mobility base names.")
  }
}

ds_policy_local <- function(Data, graph, y, iss) {
  parents <- sort(bnlearn::parents(graph, y))
  s <- ds_local_counts(Data, y, parents)
  high <- match("high", levels(Data[[y]]))
  if (is.na(high)) stop("Outcome has no 'high' level: ", y)
  alpha_config <- iss * exp(-s$log_q)
  s$posterior_high <- (s$counts[, high] + alpha_config / s$r) / (s$nj + alpha_config)
  s$default_high <- 1 / s$r
  s$root_probabilities <- lapply(parents, function(pa) {
    (tabulate(as.integer(Data[[pa]]), nbins = nlevels(Data[[pa]])) + iss / nlevels(Data[[pa]])) /
      (nrow(Data) + iss)
  })
  names(s$root_probabilities) <- parents
  s$log_weights <- rep(0, nrow(s$configurations))
  for (j in seq_along(parents)) {
    s$log_weights <- s$log_weights + log(s$root_probabilities[[j]][s$configurations[, j]])
  }
  s
}

ds_policy_probability <- function(local_fit, Data, x = NULL, value = NULL) {
  s <- local_fit
  weights_log <- s$log_weights
  take <- rep(TRUE, length(weights_log))
  if (!is.null(x) && x %in% s$parents) {
    j <- match(x, s$parents)
    code <- match(value, levels(Data[[x]]))
    if (is.na(code)) stop("Unknown intervention state.")
    take <- s$configurations[, j] == code
    weights_log <- weights_log - log(s$root_probabilities[[j]][s$configurations[, j]])
  }
  weights <- exp(weights_log[take])
  seen_mass <- sum(weights)
  if (seen_mass > 1 + 1e-8) stop("Parent-configuration probability mass exceeds one.")
  correction <- sum(weights * (s$posterior_high[take] - s$default_high))
  list(probability = min(1, max(0, s$default_high + correction)),
       correction = correction, seen_mass = min(1, max(0, seen_mass)),
       n_observed_configs = sum(take))
}

estimate_do_contrast <- function(local_fit, Data, x_var, y_var) {
  if (!all(c("low", "high") %in% levels(Data[[x_var]]))) stop("Exposure must have low & high levels.")
  low <- ds_policy_probability(local_fit, Data, x_var, "low")
  high <- ds_policy_probability(local_fit, Data, x_var, "high")
  direct <- x_var %in% local_fit$parents
  contrast <- if (direct) low$correction - high$correction else 0
  list(P_y_high_do_x_low = low$probability, P_y_high_do_x_high = high$probability,
       causal_contrast_low_minus_high = contrast, directed_path_exists = direct,
       observed_parent_mass_do_low = low$seen_mass,
       observed_parent_mass_do_high = high$seen_mass,
       unobserved_parent_mass_do_low = 1 - low$seen_mass,
       unobserved_parent_mass_do_high = 1 - high$seen_mass,
       n_observed_configs_do_low = low$n_observed_configs,
       n_observed_configs_do_high = high$n_observed_configs,
       log10_possible_parent_configurations = local_fit$log_q / log(10),
       n_outcome_parents = length(local_fit$parents))
}

policy_rows <- parent_rows <- list()
for (p in candidate_lags) {
  pn <- paste0("p", p)
  D <- df_bn_by_lag[[pn]]
  exposure_lags <- policy_exposure_lags[policy_exposure_lags <= p]
  if (policy_include_max_lag) exposure_lags <- sort(unique(c(exposure_lags, p)))
  message("Script 12: ", pn, "; exposure lags ", paste(exposure_lags, collapse = ", "))
  for (method in policy_methods) {
    entry <- ds_entry(pn, method)
    for (y in outcome_nodes) {
      local_fit <- NULL
      error <- entry$detail
      parents <- character()
      if (!is.null(entry$graph) && y %in% names(D)) {
        parents <- sort(bnlearn::parents(entry$graph, y))
        attempt <- tryCatch(ds_policy_local(D, entry$graph, y, policy_iss), error = identity)
        if (inherits(attempt, "error")) error <- conditionMessage(attempt) else local_fit <- attempt
      }
      parent_rows[[length(parent_rows) + 1L]] <- data.frame(lag = pn, model_max_lag = p,
                                                            graph_source = method, node = y, fit_status = entry$status,
                                                            n_parents = if (is.null(entry$graph) || !y %in% names(D)) NA_integer_ else length(parents),
                                                            parents = if (is.null(entry$graph)) NA_character_ else paste(parents, collapse = "; "),
                                                            ancestors = if (is.null(entry$graph)) NA_character_ else paste(parents, collapse = "; "),
                                                            stringsAsFactors = FALSE)
      for (ell in exposure_lags) for (base in x_variables_base) {
        x <- paste0(base, "_t-", ell)
        expected <- if (!is.null(policy_expected_sign) && base %in% names(policy_expected_sign)) unname(policy_expected_sign[base]) else NA_real_
        row <- data.frame(lag = pn, model_max_lag = p, exposure_lag = ell,
                          fixed_exposure_lag_comparison = ell %in% policy_exposure_lags,
                          exposure_at_model_max_lag = ell == p, graph_source = method,
                          x_base = base, x = x, y = y, x_low = "low", x_high = "high", y_high = "high",
                          fit_status = entry$status, estimation_status = "unavailable_graph", detail = error,
                          identification_status = "not_assessed", directed_path_exists = NA,
                          symbolic_effect = NA_character_,
                          contrast_definition = "P(Y=high|do(X=low))-P(Y=high|do(X=high))",
                          P_y_high_do_x_low = NA_real_, P_y_high_do_x_high = NA_real_,
                          causal_contrast_low_minus_high = NA_real_,
                          observed_parent_mass_do_low = NA_real_, observed_parent_mass_do_high = NA_real_,
                          unobserved_parent_mass_do_low = NA_real_, unobserved_parent_mass_do_high = NA_real_,
                          n_observed_configs_do_low = NA_integer_, n_observed_configs_do_high = NA_integer_,
                          log10_possible_parent_configurations = NA_real_, n_outcome_parents = NA_integer_,
                          exposure_low_count = if (x %in% names(D)) sum(D[[x]] == "low") else NA_integer_,
                          exposure_high_count = if (x %in% names(D)) sum(D[[x]] == "high") else NA_integer_,
                          outcome_high_count = if (y %in% names(D)) sum(D[[y]] == "high") else NA_integer_,
                          observed_configuration_coverage_pass = NA, support_status = "not_assessed",
                          expected_sign = expected, estimated_sign = NA_real_, direction_matches_expectation = NA,
                          iss = policy_iss, max_unobserved_mass_threshold = policy_max_unobserved_mass,
                          uncertainty = "none; exact fitted-model probability, no sampling or data uncertainty interval",
                          stringsAsFactors = FALSE)
        if (!is.null(entry$graph)) {
          if (!all(c(x, y) %in% names(D))) {
            row$estimation_status <- "missing_node"
            row$detail <- paste("Missing:", paste(setdiff(c(x, y), names(D)), collapse = "; "))
          } else if (is.null(local_fit)) {
            row$estimation_status <- "local_estimation_failed"
          } else {
            ans <- tryCatch(estimate_do_contrast(local_fit, D, x, y), error = identity)
            if (inherits(ans, "error")) {
              row$estimation_status <- "local_estimation_failed"
              row$detail <- conditionMessage(ans)
            } else {
              for (name in names(ans)) row[[name]] <- ans[[name]]
              row$detail <- NA_character_
              row$identification_status <- "identified_under_assumed_causally_sufficient_DAG"
              row$estimation_status <- if (ans$directed_path_exists) "exact_sparse_BDe" else "structural_zero"
              row$symbolic_effect <- if (ans$directed_path_exists) {
                "P(Y=high|do(X=x)) = sum_z P(Y=high|X=x,Z=z) product_j P(Z_j=z_j); Z=other parents of Y"
              } else "P(Y=high|do(X=x)) = P_G(Y=high); contrast=0 because X is not an ancestor of Y"
              row$observed_configuration_coverage_pass <-
                max(ans$unobserved_parent_mass_do_low, ans$unobserved_parent_mass_do_high) <= policy_max_unobserved_mass + 1e-12
              row$support_status <- if (row$exposure_low_count == 0L || row$exposure_high_count == 0L) {
                "exposure_state_unobserved"
              } else if (row$observed_configuration_coverage_pass) "observed_configuration_coverage_pass" else "low_observed_configuration_coverage"
              row$estimated_sign <- if (abs(ans$causal_contrast_low_minus_high) <= 1e-12) 0 else sign(ans$causal_contrast_low_minus_high)
              if (!is.na(expected)) row$direction_matches_expectation <- row$estimated_sign == expected
            }
          }
        }
        policy_rows[[length(policy_rows) + 1L]] <- row
      }
    }
  }
}

policy_contrasts_df <- ds_bind(policy_rows)
parents_ancestors_df <- ds_bind(parent_rows)
symbolic_effects_df <- policy_contrasts_df[, c("lag", "model_max_lag", "exposure_lag", "graph_source",
                                               "x", "y", "fit_status", "identification_status", "directed_path_exists", "symbolic_effect", "detail"), drop = FALSE]
policy_contrasts_by_lag <- split(policy_contrasts_df, policy_contrasts_df$lag)
parents_ancestors_by_lag <- split(parents_ancestors_df, parents_ancestors_df$lag)
symbolic_effects_by_lag <- split(symbolic_effects_df, symbolic_effects_df$lag)
conditional_contrasts_df <- policy_contrasts_df
conditional_contrasts_by_lag <- policy_contrasts_by_lag
ds_csv(policy_contrasts_df, "policy_do_contrasts_by_method_and_lag.csv")
ds_csv(policy_contrasts_df[policy_contrasts_df$lag == selected_name, , drop = FALSE],
       paste0("policy_do_contrasts_selected_", selected_name, ".csv"))
ds_csv(parents_ancestors_df, "parents_ancestors_by_method_and_lag.csv")
ds_csv(symbolic_effects_df, "symbolic_causal_effects_by_method_and_lag.csv")
ds_rds(policy_contrasts_by_lag, "policy_contrasts_by_lag.rds")
policy_settings <- list(iss = policy_iss, methods = policy_methods,
                        fixed_exposure_lags = policy_exposure_lags, include_model_max_lag_exposure = policy_include_max_lag,
                        max_unobserved_mass = policy_max_unobserved_mass, expected_sign = policy_expected_sign,
                        root_distribution = "independent Bayesian-smoothed marginals, as implied by the saved DAGs",
                        estimator = "exact sparse evaluation of uniform Dirichlet posterior means",
                        uncertainty = "No data/parameter/graph uncertainty interval is estimated.")
ds_rds(policy_settings, "policy_settings.rds")
message("Script 12 complete: ", nrow(policy_contrasts_df),
        " contrasts/status rows. Check support_status before interpreting a nonzero effect.")
