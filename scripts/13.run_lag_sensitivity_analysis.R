# 13.run_lag_sensitivity_analysis.R

ds_require(c("structure_metrics_df", "pairwise_shd_scores_df", "policy_contrasts_df",
             "edge_stability_by_method", "policy_settings"))
selected_name <- paste0("p", selected_p)
all_graph_edge_counts_df <- structure_metrics_df[, c("lag", "model_max_lag", "method",
                                                     "fit_status", "n_edges", "n_nodes", "n_rows", "detail"), drop = FALSE]
graph_edge_counts_df <- all_graph_edge_counts_df[
  all_graph_edge_counts_df$method %in% c(ds_econ_methods, "model_average"), , drop = FALSE]
causal_ml_edge_counts_df <- all_graph_edge_counts_df[
  grepl("^causal_ml_", all_graph_edge_counts_df$method), , drop = FALSE]

ds_csv(graph_edge_counts_df, "final_graph_edge_counts_by_method_and_lag.csv")
ds_csv(causal_ml_edge_counts_df, "final_causal_ml_edge_counts_by_algorithm_and_lag.csv")
ds_csv(all_graph_edge_counts_df, "final_all_graph_edge_counts_by_method_and_lag.csv")
ds_csv(structure_metrics_df, "final_structure_metrics_by_method_and_lag.csv")
ds_csv(average_bn_metrics_df, "final_average_bn_metrics_by_lag.csv")
ds_csv(pairwise_shd_scores_df, "final_pairwise_shd_by_lag.csv")
ds_csv(stable_edges_all_methods, "final_edge_stability_all_methods.csv")
ds_csv(stable_edges_all_methods[stable_edges_all_methods$stable_across_specifications, , drop = FALSE],
       "final_stable_edges_all_methods.csv")
ds_csv(parents_ancestors_df, "final_parents_ancestors_by_lag.csv")
ds_csv(symbolic_effects_df, "final_symbolic_causal_effects_by_lag.csv")
ds_csv(policy_contrasts_df, "final_policy_do_contrasts_by_method_and_lag.csv")

if (!is.null(ds_get("lag_diagnostics"))) ds_csv(lag_diagnostics, "final_lag_diagnostics.csv")

fixed_policy <- policy_contrasts_df[policy_contrasts_df$fixed_exposure_lag_comparison, , drop = FALSE]
policy_max_lag_exposures <- policy_contrasts_df[policy_contrasts_df$exposure_at_model_max_lag, , drop = FALSE]
ds_csv(fixed_policy, "policy_fixed_exposure_lag_sensitivity.csv")
ds_csv(policy_max_lag_exposures, "policy_exposure_at_model_max_lag.csv")

policy_sensitivity_rows <- list()
if (nrow(fixed_policy)) {
  keys <- do.call(paste, c(fixed_policy[c("graph_source", "x_base", "y", "exposure_lag")], sep = "\034"))
  groups <- split(seq_len(nrow(fixed_policy)), keys)
  policy_sensitivity_rows <- lapply(groups, function(ix) {
    z <- fixed_policy[ix, , drop = FALSE]
    available <- is.finite(z$causal_contrast_low_minus_high)
    covered <- available & z$support_status == "observed_configuration_coverage_pass"
    signs <- z$estimated_sign[covered]
    data.frame(graph_source = z$graph_source[1], x_base = z$x_base[1], y = z$y[1],
               exposure_lag = z$exposure_lag[1], n_eligible_specifications = nrow(z),
               n_estimated = sum(available), n_unavailable = sum(!available),
               n_structural_zero = sum(z$estimation_status == "structural_zero"),
               n_with_observed_configuration_coverage = sum(covered),
               min_contrast = if (any(available)) min(z$causal_contrast_low_minus_high[available]) else NA_real_,
               max_contrast = if (any(available)) max(z$causal_contrast_low_minus_high[available]) else NA_real_,
               same_sign_among_covered_specifications = if (sum(covered) >= 2L) length(unique(signs)) == 1L else NA,
               stringsAsFactors = FALSE)
  })
}

policy_lag_sensitivity_df <- ds_bind(policy_sensitivity_rows)
ds_csv(policy_lag_sensitivity_df, "policy_fixed_lag_sensitivity_summary.csv")
policy_support_counts <- as.data.frame(table(
  lag = factor(policy_contrasts_df$lag, levels = paste0("p", candidate_lags)),
  estimation_status = policy_contrasts_df$estimation_status,
  support_status = policy_contrasts_df$support_status), stringsAsFactors = FALSE)
policy_support_counts <- policy_support_counts[policy_support_counts$Freq > 0, , drop = FALSE]
ds_csv(policy_support_counts, "policy_support_counts.csv")

ds_make_plots <- function() {
  grDevices::pdf(ds_path("lag_sensitivity_plots.pdf"), width = 11, height = 7, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  methods <- setdiff(ds_methods, "knowledge")
  palette <- grDevices::hcl.colors(length(methods), "Dark 3")
  plot_by_lag <- function(column, title, ylab) {
    mat <- matrix(NA_real_, length(candidate_lags), length(methods), dimnames = list(NULL, methods))
    for (j in seq_along(methods)) for (i in seq_along(candidate_lags)) {
      z <- structure_metrics_df[structure_metrics_df$method == methods[j] &
                                  structure_metrics_df$model_max_lag == candidate_lags[i], column]
      if (length(z) == 1L && is.finite(z)) mat[i, j] <- z
    }
    show <- colSums(is.finite(mat)) > 0
    if (!any(show)) return(invisible(NULL))
    graphics::par(mar = c(7, 5, 4, 13), xpd = FALSE)
    graphics::matplot(candidate_lags, mat[, show, drop = FALSE], type = "b", pch = 16,
                      lty = seq_len(sum(show)), col = palette[show], xaxt = "n", xlab = "Model maximum lag (days)",
                      ylab = ylab, main = title)
    graphics::axis(1, at = candidate_lags)
    usr <- graphics::par("usr")
    graphics::legend(x = usr[2] + 0.035 * (usr[2] - usr[1]), y = usr[4],
                     xjust = 0, yjust = 1, legend = methods[show],
                     col = palette[show], lty = seq_len(sum(show)), pch = 16, bty = "n", cex = 0.75, xpd = NA)
    graphics::mtext("Unavailable fits are omitted; an empty successful graph has zero edges.", side = 1, line = 5, cex = 0.75)
  }
  plot_by_lag("n_edges", "Selected edges across maximum-lag specifications", "Number of edges")
  plot_by_lag("SHD_per_possible_edge", "Distance from the specified knowledge graph", "SHD / admissible edges")
  z <- structure_metrics_df[structure_metrics_df$lag == selected_name &
                              structure_metrics_df$method != "knowledge" & is.finite(structure_metrics_df$within_lag_BN_score_rank), , drop = FALSE]
  if (nrow(z)) {
    z <- z[order(z$within_lag_BN_score_rank, decreasing = TRUE), , drop = FALSE]
    graphics::par(mar = c(7, 13, 4, 2), xpd = FALSE)
    graphics::plot(z$within_lag_BN_score_rank, seq_len(nrow(z)), yaxt = "n", pch = 16,
                   ylim = c(0.5, nrow(z) + 0.5), ylab = "",
                   main = paste("Descriptive discrete BIC ranking at", selected_name),
                   xlab = "Within-specification rank (1 is highest score)")
    graphics::axis(2, at = seq_len(nrow(z)), labels = z$method, las = 1, cex.axis = 0.8)
    graphics::mtext("Full-sample nominal-CPT BIC; sparse configurations can make the approximation unreliable.",
                    side = 1, line = 5, cex = 0.7)
  }
}
plot_status <- "disabled"
if (isTRUE(getOption("downstream.make_plots", TRUE))) {
  plot_status <- tryCatch({ds_make_plots(); "completed"}, error = function(e) {
    warning("Plots could not be produced: ", conditionMessage(e))
    paste("failed:", conditionMessage(e))
  })
}
lag_sensitivity_summary <- list(version = downstream_version,
                                candidate_lags = candidate_lags, selected_p = selected_p, selected_name = selected_name,
                                input_audit = downstream_input_audit, model_average_pool_audit = model_average_pool_audit,
                                graph_edge_counts = graph_edge_counts_df, causal_ml_edge_counts = causal_ml_edge_counts_df,
                                all_graph_edge_counts = all_graph_edge_counts_df, average_bn_metrics = average_bn_metrics_df,
                                structure_metrics = structure_metrics_df, local_structure_metrics = local_structure_metrics_df,
                                pairwise_shd = pairwise_shd_scores_df, edge_stability = stable_edges_all_methods,
                                stable_edges = stable_edges_all_methods[stable_edges_all_methods$stable_across_specifications, , drop = FALSE],
                                parents_ancestors = parents_ancestors_df, symbolic_effects = symbolic_effects_df,
                                policy_contrasts = policy_contrasts_df, policy_fixed_lag_sensitivity = policy_lag_sensitivity_df,
                                policy_support_counts = policy_support_counts, policy_settings = policy_settings,
                                lag_diagnostics = ds_get("lag_diagnostics"), plot_status = plot_status,
                                interpretation_notes = c(
                                  "Knowledge graph is a declared reference, not a verified ground truth.",
                                  "No saved learnt edge is pruned or changed.",
                                  "Structural zero means no effect in this restricted graph, not evidence of no real-world effect.",
                                  "Exact probabilities eliminate Monte Carlo error but do not quantify data or graph uncertainty.",
                                  "The exposure-at-maximum-lag table changes the intervention horizon and is kept separate."))
ds_rds(lag_sensitivity_summary, "final_lag_sensitivity_summary.rds")
writeLines(capture.output(utils::sessionInfo()), ds_path("sessionInfo.txt"))
writeLines(lag_sensitivity_summary$interpretation_notes, ds_path("interpretation_notes.txt"))
message("Scripts 9-13 complete. Selected lag: ", selected_name,
        "; summary: ", ds_path("final_lag_sensitivity_summary.rds"))