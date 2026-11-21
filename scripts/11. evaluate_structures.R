# 11. evaluate_structures.R

ds_require(c("average_bn_fit_results_by_lag", "compute_bn_metrics"))
structure_metrics_by_lag <- local_structure_metrics_by_lag <- list()
pairwise_shd_by_lag <- list()
evaluated_graphs_by_lag <- list()

for (p in candidate_lags) {
  pn <- paste0("p", p)
  message("Script 11: evaluating all methods at ", pn)
  entries <- ds_registry(pn)
  evaluated_graphs_by_lag[[pn]] <- lapply(entries, `[[`, "graph")
  cache <- new.env(parent = emptyenv())
  rows <- local_rows <- list()
  for (method in names(entries)) {
    z <- entries[[method]]
    metrics <- NULL
    if (!is.null(z$graph)) {
      attempt <- tryCatch(compute_bn_metrics(z$graph, kb_by_lag[[pn]],
                                             df_bn_by_lag[[pn]], cache = cache), error = identity)
      if (inherits(attempt, "error")) {
        z$detail <- paste("Metric failure:", conditionMessage(attempt))
      } else {
        metrics <- attempt
        lm <- metrics$local_metrics
        local_rows[[method]] <- data.frame(lag = pn, method = method, lm,
                                           stringsAsFactors = FALSE)
      }
    }
    rows[[method]] <- ds_metric_row(metrics, pn, method, z$status, z$detail, z$warnings)
    if (!is.null(z$graph) && is.null(metrics)) rows[[method]]$score_status <- "metric_error"
  }
  
  scores <- ds_bind(rows)
  scores$within_lag_BN_score_rank <- NA_real_
  keep <- scores$method != "knowledge" & is.finite(scores$BN_Score)
  scores$within_lag_BN_score_rank[keep] <- rank(-scores$BN_Score[keep], ties.method = "min")
  structure_metrics_by_lag[[pn]] <- scores
  local_structure_metrics_by_lag[[pn]] <- ds_bind(local_rows)
  pair_methods <- setdiff(names(entries), "knowledge")
  pairs <- combn(pair_methods, 2L, simplify = FALSE)
  pairwise_shd_by_lag[[pn]] <- ds_bind(lapply(pairs, function(pair) {
    g1 <- entries[[pair[1]]]$graph
    g2 <- entries[[pair[2]]]$graph
    shd_value <- temporal_shd <- jaccard <- NA_real_
    status <- "unavailable_graph"
    if (!is.null(g1) && !is.null(g2)) {
      e1 <- ds_edge_key(ds_arcs(g1))
      e2 <- ds_edge_key(ds_arcs(g2))
      temporal_shd <- length(setdiff(e1, e2)) + length(setdiff(e2, e1))
      shd_value <- as.numeric(bnlearn::shd(g1, g2, wlbl = FALSE))
      un <- length(union(e1, e2))
      jaccard <- if (un) length(intersect(e1, e2)) / un else 1
      status <- "computed"
    }
    data.frame(lag = pn, model_max_lag = p, Method_1 = pair[1], Method_2 = pair[2],
               Pair = paste(pair, collapse = " vs "), status = status,
               fit_status_1 = entries[[pair[1]]]$status, fit_status_2 = entries[[pair[2]]]$status,
               SHD = shd_value, Temporal_SHD = temporal_shd,
               SHD_definition = "bnlearn_CPDAG_comparison_wlbl_FALSE",
               Jaccard = jaccard, stringsAsFactors = FALSE)
  }))
}

structure_metrics_df <- ds_bind(structure_metrics_by_lag)
local_structure_metrics_df <- ds_bind(local_structure_metrics_by_lag)
pairwise_shd_scores_df <- ds_bind(pairwise_shd_by_lag)
ds_csv(structure_metrics_df, "structure_metrics_by_method_and_lag.csv")
ds_csv(structure_metrics_df[structure_metrics_df$lag == selected_name, , drop = FALSE],
       paste0("structure_metrics_selected_", selected_name, ".csv"))
ds_csv(local_structure_metrics_df, "local_structure_metrics_by_method_and_lag.csv")
ds_csv(pairwise_shd_scores_df, "pairwise_shd_by_lag.csv")
ds_rds(structure_metrics_by_lag, "structure_metrics_by_lag.rds")

stability_min_frequency <- getOption("downstream.stability_min_frequency", 0.8)
stability_min_available <- getOption("downstream.stability_min_available", 2L)
if (length(stability_min_frequency) != 1L || !is.finite(stability_min_frequency) ||
    stability_min_frequency < 0 || stability_min_frequency > 1 ||
    length(stability_min_available) != 1L || !is.finite(stability_min_available) ||
    stability_min_available < 2 || stability_min_available != as.integer(stability_min_available)) {
  stop("Invalid stability settings.")
}

stability_prototype <- data.frame(method = character(), from = character(), to = character(),
                                  edge = character(), required_lag = integer(), count = integer(),
                                  n_eligible_specifications = integer(), n_available_eligible = integer(),
                                  n_unavailable_eligible = integer(), freq = numeric(),
                                  present_specifications = character(), available_eligible_specifications = character(),
                                  stable_across_specifications = logical(), stringsAsFactors = FALSE)
edge_stability_by_method <- list()

for (method in setdiff(ds_methods, "knowledge")) {
  entries <- setNames(lapply(paste0("p", candidate_lags), function(pn) ds_entry(pn, method)),
                      paste0("p", candidate_lags))
  available <- vapply(entries, function(z) !is.null(z$graph), logical(1))
  edge_sets <- lapply(entries, function(z) if (is.null(z$graph)) character() else ds_edge_key(ds_arcs(z$graph)))
  unique_arcs <- unique(ds_bind(lapply(entries, function(z) if (is.null(z$graph)) NULL else ds_arcs(z$graph)),
                                data.frame(from = character(), to = character())))
  rows <- lapply(seq_len(nrow(unique_arcs)), function(i) {
    a <- unique_arcs[i, , drop = FALSE]
    needed <- ds_lag(a$from)
    eligible <- candidate_lags >= needed
    avail <- eligible & available
    present <- vapply(edge_sets, function(e) ds_edge_key(a) %in% e, logical(1))
    n_available <- sum(avail)
    count <- sum(present & avail)
    frequency <- if (n_available) count / n_available else NA_real_
    data.frame(method = method, from = a$from, to = a$to,
               edge = paste(a$from, "->", a$to), required_lag = needed, count = count,
               n_eligible_specifications = sum(eligible), n_available_eligible = n_available,
               n_unavailable_eligible = sum(eligible & !available), freq = frequency,
               present_specifications = paste(names(entries)[present & avail], collapse = "; "),
               available_eligible_specifications = paste(names(entries)[avail], collapse = "; "),
               stable_across_specifications = n_available >= stability_min_available && frequency >= stability_min_frequency,
               stringsAsFactors = FALSE)
  })
  tab <- ds_bind(rows, stability_prototype)
  tab <- tab[order(-tab$freq, -tab$n_available_eligible, -tab$count), , drop = FALSE]
  edge_stability_by_method[[method]] <- tab
  ds_csv(tab, paste0(method, "_edge_stability_by_lag.csv"))
}
stable_edges_all_methods <- ds_bind(edge_stability_by_method, stability_prototype)
ds_csv(stable_edges_all_methods, "edge_stability_all_methods.csv")
ds_rds(edge_stability_by_method, "edge_stability_by_method.rds")
message("Script 11 complete: ", nrow(structure_metrics_df))
