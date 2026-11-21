# 10.model_averaging.R ----------------------------------------------------------

if (!exists("downstream_version", inherits = FALSE)) stop("Run Script 9.")
if (!exists("selected_p", inherits = FALSE)) stop("Run Script 3.")

ds_local_counts <- function(Data, target, parents) {
  n <- nrow(Data)
  if (length(parents)) {
    codes <- lapply(Data[parents], as.integer)
    keys <- do.call(paste, c(unname(codes), sep = ","))
  } else keys <- rep("root", n)
  unique_keys <- unique(keys)
  group <- match(keys, unique_keys)
  r <- nlevels(Data[[target]])
  counts <- matrix(tabulate((group - 1L) * r + as.integer(Data[[target]]),
                            nbins = length(unique_keys) * r), ncol = r, byrow = TRUE)
  colnames(counts) <- levels(Data[[target]])
  nj <- rowSums(counts)
  first <- match(unique_keys, keys)
  configurations <- if (length(parents)) {
    vapply(Data[parents], function(z) as.integer(z[first]), integer(length(first)))
  } else matrix(integer(0), nrow = length(first), ncol = 0L)
  configurations <- matrix(configurations, nrow = length(first), ncol = length(parents),
                           dimnames = list(NULL, parents))
  log_q <- sum(log(vapply(Data[parents], nlevels, integer(1))))
  list(counts = counts, nj = nj, configurations = configurations,
       parents = parents, target = target, r = r, log_q = log_q)
}

ds_logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

ds_local_metrics <- function(Data, target, parents) {
  s <- ds_local_counts(Data, target, parents)
  proportions <- s$counts / s$nj
  positive <- s$counts > 0
  ll <- sum(s$counts[positive] * log(proportions[positive]))
  data.frame(node = target, current = grepl("_t$", target), n_parents = length(parents),
             log_q = s$log_q, log_df = log(s$r - 1) + s$log_q,
             nominal_df = (s$r - 1) * prod(vapply(Data[parents], nlevels, integer(1))),
             n_observed_parent_configs = length(s$nj),
             observed_config_fraction = exp(log(length(s$nj)) - s$log_q),
             singleton_config_fraction = mean(s$nj == 1L), min_config_count = min(s$nj),
             LogLikelihood = ll, stringsAsFactors = FALSE)
}

ds_components <- function(nodes, a) {
  component <- seq_along(nodes)
  if (nrow(a)) for (i in seq_len(nrow(a))) {
    left <- component[match(a$from[i], nodes)]
    right <- component[match(a$to[i], nodes)]
    if (left != right) component[component == right] <- left
  }
  length(unique(component))
}

compute_bn_metrics <- function(bn_str, kb, Data, debug = FALSE, cache = NULL) {
  g <- align_bn_to_nodes(bn_str, names(Data))
  ref <- align_bn_to_nodes(kb, names(Data))
  a <- ds_arcs(g)
  parent_list <- split(a$from, a$to)
  locals <- lapply(seq_along(Data), function(i) {
    node <- names(Data)[i]
    pa <- sort(parent_list[[node]])
    cache_key <- paste(i, paste(match(pa, names(Data)), collapse = ","), sep = ":")
    if (!is.null(cache) && exists(cache_key, envir = cache, inherits = FALSE)) {
      return(get(cache_key, envir = cache, inherits = FALSE))
    }
    ans <- ds_local_metrics(Data, node, pa)
    if (!is.null(cache)) assign(cache_key, ans, envir = cache)
    ans
  })
  local_metrics <- ds_bind(locals)
  log_d <- ds_logsumexp(local_metrics$log_df)
  n_params <- sum(local_metrics$nominal_df)
  ll <- sum(local_metrics$LogLikelihood)
  penalty <- (n_params / 2) * log(nrow(Data))
  keys <- ds_edge_key(a)
  ref_keys <- ds_edge_key(ds_arcs(ref))
  shared <- length(intersect(keys, ref_keys))
  union_n <- length(union(keys, ref_keys))
  temporal_shd <- length(setdiff(keys, ref_keys)) + length(setdiff(ref_keys, keys))
  shd_value <- as.numeric(bnlearn::shd(g, ref, wlbl = FALSE))
  current <- local_metrics$current
  sparse <- any(local_metrics$observed_config_fraction[current] < 1 - 1e-12)
  summary <- list(SHD = shd_value, Temporal_SHD = temporal_shd,
                  SHD_definition = "bnlearn_CPDAG_comparison_wlbl_FALSE", BN_Score = ll - penalty,
                  n_params = n_params, log10_n_params = log_d / log(10), LogLikelihood = ll,
                  n_edges = nrow(a), n_components = ds_components(names(Data), a),
                  max_parents = max(local_metrics$n_parents), n_nodes = ncol(Data), n_rows = nrow(Data),
                  n_possible_edges = sum(current) * sum(!current),
                  SHD_per_possible_edge = shd_value / (sum(current) * sum(!current)),
                  Jaccard_with_knowledge = if (union_n) shared / union_n else 1,
                  shared_knowledge_edges = shared,
                  sparse_parent_configurations = sparse,
                  min_observed_config_fraction_current = min(local_metrics$observed_config_fraction[current]),
                  score_status = if (!is.finite(penalty)) "nominal_penalty_overflow" else if (sparse) "sparse_configurations" else "complete_configurations",
                  score_definition = "discrete_MLE_logLik_minus_half_nominal_df_log_N")
  c(summary, list(local_metrics = local_metrics))
}

ds_metric_row <- function(metrics, pn, method, status = "success", detail = NA_character_,
                          warnings = NA_character_) {
  out <- data.frame(lag = pn, model_max_lag = as.integer(sub("^p", "", pn)),
                    method = method, fit_status = status, detail = detail, warnings = warnings,
                    SHD = NA_real_, Temporal_SHD = NA_real_, SHD_definition = "bnlearn_CPDAG_comparison_wlbl_FALSE",
                    BN_Score = NA_real_, n_params = NA_real_, log10_n_params = NA_real_,
                    LogLikelihood = NA_real_, n_edges = NA_real_, n_components = NA_real_,
                    max_parents = NA_real_, n_nodes = ncol(df_bn_by_lag[[pn]]),
                    n_rows = nrow(df_bn_by_lag[[pn]]), n_possible_edges = NA_real_,
                    SHD_per_possible_edge = NA_real_, Jaccard_with_knowledge = NA_real_,
                    shared_knowledge_edges = NA_real_, sparse_parent_configurations = NA,
                    min_observed_config_fraction_current = NA_real_, score_status = "not_evaluated",
                    score_definition = "discrete_MLE_logLik_minus_half_nominal_df_log_N",
                    stringsAsFactors = FALSE)
  if (!is.null(metrics)) for (name in intersect(names(out), names(metrics))) out[[name]] <- metrics[[name]]
  out
}

ds_vote_bic <- function(g, Data, cache = NULL) {
  g <- align_bn_to_nodes(g, names(Data))
  a <- ds_arcs(g)
  parent_list <- split(a$from, a$to)
  locals <- lapply(seq_along(Data), function(i) {
    node <- names(Data)[i]
    pa <- sort(parent_list[[node]])
    key <- paste(i, paste(match(pa, names(Data)), collapse = ","), sep = ":")
    if (!is.null(cache) && exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }
    res <- ds_local_metrics(Data, node, pa)
    if (!is.null(cache)) assign(key, res, envir = cache)
    res
  })
  local_metrics <- ds_bind(locals)
  ll <- sum(local_metrics$LogLikelihood)
  n_params <- sum(local_metrics$nominal_df)
  data.frame(BN_Score = ll - n_params * log(nrow(Data)) / 2,
             LogLikelihood = ll, n_params = n_params,
             max_parents = max(local_metrics$n_parents))
}

model_graphs_by_lag <- model_edge_counts_by_lag <- average_bn_by_lag <- list()
average_bn_fit_results_by_lag <- average_bn_metrics_by_lag <- average_bn_edges_by_lag <- list()
pool_rows <- average_rows <- list()
vote_prototype <- data.frame(from = character(), to = character(), Freq = integer(),
                             available_methods = integer(), attempted_methods = integer(), support_fraction = numeric(),
                             supporting_methods = character(), selected = logical(), stringsAsFactors = FALSE)

threshold_audit_rows <- list()
for (p in candidate_lags) {
  pn <- paste0("p", p)
  entries <- ds_registry(pn, ds_econ_methods)
  available <- names(entries)[vapply(entries, function(z) !is.null(z$graph), logical(1))]
  graphs <- lapply(entries[available], `[[`, "graph")
  model_graphs_by_lag[[pn]] <- graphs
  if (!length(available)) stop("No usable econometric graphs at ", pn)
  
  edge_rows <- ds_bind(lapply(available, function(method) {
    a <- ds_arcs(graphs[[method]])
    a$method <- rep(method, nrow(a))
    a
  }), data.frame(from = character(), to = character(), method = character()))
  
  votes <- vote_prototype
  if (nrow(edge_rows)) {
    groups <- split(seq_len(nrow(edge_rows)), ds_edge_key(edge_rows))
    votes <- ds_bind(lapply(groups, function(ix) {
      supporting <- unique(edge_rows$method[ix])
      data.frame(from = edge_rows$from[ix[1]], to = edge_rows$to[ix[1]],
                 Freq = length(supporting), available_methods = length(available),
                 attempted_methods = length(ds_econ_methods),
                 support_fraction = length(supporting) / length(available),
                 supporting_methods = paste(supporting, collapse = "; "),
                 selected = FALSE, stringsAsFactors = FALSE)
    }))
    votes <- votes[order(-votes$Freq, votes$to, votes$from), , drop = FALSE]
  }
  model_edge_counts_by_lag[[pn]] <- votes
  
  score_cache <- new.env(parent = emptyenv())
  for (v in seq_len(length(available))) {
    chosen <- votes[votes$Freq >= v, , drop = FALSE]
    graph <- ds_graph(names(df_bn_by_lag[[pn]]), chosen)
    score <- ds_vote_bic(graph, df_bn_by_lag[[pn]], cache = score_cache)
    threshold_audit_rows[[paste0(pn, "_votes", v)]] <- data.frame(
      lag = pn, model_max_lag = p,
      min_votes = as.integer(v), available_methods = length(available),
      n_edges = nrow(chosen), score, stringsAsFactors = FALSE)
  }
}

model_average_threshold_bic <- ds_bind(threshold_audit_rows)
if (nrow(model_average_threshold_bic) == 0L) stop("No thresholds could be scored.")

model_average_threshold_bic$best_within_lag <- FALSE
for (pn in unique(model_average_threshold_bic$lag)) {
  idx <- which(model_average_threshold_bic$lag == pn &
                 is.finite(model_average_threshold_bic$BN_Score))
  if (!length(idx)) stop("No finite threshold BIC for ", pn)
  best <- idx[order(-model_average_threshold_bic$BN_Score[idx],
                    -model_average_threshold_bic$min_votes[idx])][1L]
  model_average_threshold_bic$best_within_lag[best] <- TRUE
}

selected_name <- paste0("p", selected_p)
selected_candidates <- model_average_threshold_bic[
  model_average_threshold_bic$lag == selected_name &
    model_average_threshold_bic$best_within_lag, , drop = FALSE]
if (nrow(selected_candidates) != 1L) {
  stop("Could not determine an optimal vote threshold for ", selected_name)
}
model_average_min_occurrences <- as.integer(selected_candidates$min_votes[1L])
model_average_threshold_bic$selected_global_threshold <-
  model_average_threshold_bic$min_votes == model_average_min_occurrences

ds_csv(model_average_threshold_bic, "model_average_threshold_BIC_by_lag.csv")
ds_csv(selected_candidates, paste0("model_average_optimal_votes_", selected_name, ".csv"))

for (p in candidate_lags) {
  pn <- paste0("p", p)
  entries <- ds_registry(pn, ds_econ_methods)
  available <- names(model_graphs_by_lag[[pn]])
  votes <- model_edge_counts_by_lag[[pn]]
  votes$selected <- votes$Freq >= model_average_min_occurrences
  model_edge_counts_by_lag[[pn]] <- votes
  
  pool_rows[[pn]] <- ds_bind(lapply(names(entries), function(method) {
    z <- entries[[method]]
    data.frame(lag = pn, method = method, status = z$status, error = z$detail,
               included_in_pool = !is.null(z$graph), n_available = length(available),
               n_attempted = length(ds_econ_methods), min_votes = model_average_min_occurrences,
               stringsAsFactors = FALSE)
  }))
  
  if (length(available) < model_average_min_occurrences) {
    g <- NULL
    status <- "insufficient_available_methods"
    detail <- paste("Available:", length(available),
                    "Required votes:", model_average_min_occurrences)
    chosen <- votes[FALSE, , drop = FALSE]
  } else {
    chosen <- votes[votes$selected, , drop = FALSE]
    g <- ds_graph(names(df_bn_by_lag[[pn]]), chosen)
    status <- if (nrow(chosen)) "success" else "success_no_edges"
    detail <- NA_character_
    ds_dot(g, paste0("average_bn_", pn, "_dot_file.dot"))
  }
  average_bn_by_lag[pn] <- list(g)
  average_bn_edges_by_lag[[pn]] <- chosen
  average_bn_fit_results_by_lag[[pn]] <- list(
    status = status, graph = g, error = detail,
    diagnostics = list(available_methods = available,
                       min_votes = model_average_min_occurrences,
                       calibrated_at = selected_name))
  metrics <- if (is.null(g)) NULL else
    compute_bn_metrics(g, kb_by_lag[[pn]], df_bn_by_lag[[pn]])
  average_bn_metrics_by_lag[pn] <- list(metrics)
  average_rows[[pn]] <- ds_metric_row(metrics, pn, "model_average", status, detail)
  ds_csv(votes, paste0("model_edge_counts_", pn, ".csv"))
  ds_csv(chosen, paste0("average_bn_edges_", pn, ".csv"))
  message("Script 10: ", pn, ": ", length(available), " available methods; ",
          status, "; consensus requires ", model_average_min_occurrences, " votes")
}
model_average_pool_audit <- ds_bind(pool_rows)
average_bn_metrics_df <- ds_bind(average_rows)
ds_csv(model_average_pool_audit, "model_average_pool_audit.csv")
ds_csv(average_bn_metrics_df, "average_bn_metrics_by_lag.csv")
ds_rds(average_bn_by_lag, "average_bn_by_lag.rds")
ds_rds(average_bn_fit_results_by_lag, "average_bn_fit_results_by_lag.rds")
message("Script 10 complete. BIC-selected vote threshold at ", selected_name,
        " = ", model_average_min_occurrences, " out of ",
        selected_candidates$available_methods[1L], " available methods.")
