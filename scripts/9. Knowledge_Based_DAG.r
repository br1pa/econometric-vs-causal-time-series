# 9. Knowledge_Based_DAG.r

if (!requireNamespace("bnlearn", quietly = TRUE)) {
  stop("Package 'bnlearn' is required. Install it before running these scripts.")
}

ds_env <- environment()
ds_get <- function(name) get0(name, envir = ds_env, inherits = FALSE)
ds_require <- function(names) {
  missing <- names[vapply(names, function(z) is.null(ds_get(z)), logical(1))]
  if (length(missing)) stop("Missing objects: ", paste(missing, collapse = ", "),
                            ". Run preceding scripts.")
}

ds_require(c("candidate_lags", "selected_p", "df_bn_by_lag",
             "empty_net_by_lag", "factor_dummy_map"))

if (!is.numeric(candidate_lags) || !length(candidate_lags) ||
    anyNA(candidate_lags) || any(!is.finite(candidate_lags)) ||
    any(candidate_lags < 1 | candidate_lags != as.integer(candidate_lags))) {
  stop("candidate_lags must contain positive integers.")
}

candidate_lags <- sort(unique(as.integer(candidate_lags)))

if (length(selected_p) != 1L || is.na(selected_p) ||
    !selected_p %in% candidate_lags) {
  stop("selected_p must be the saved selection from Script 3; no fallback is used.")
}

selected_name <- paste0("p", selected_p)
downstream_output_dir <- getOption("downstream.output_dir", "results_9_13")
dir.create(downstream_output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(downstream_output_dir)) stop("Cannot create downstream output directory.")

ds_path <- function(file) file.path(downstream_output_dir, file)
ds_csv <- function(x, file) utils::write.csv(x, ds_path(file), row.names = FALSE, na = "NA")
ds_rds <- function(x, file) saveRDS(x, ds_path(file))

ds_bind <- function(rows, prototype = data.frame()) {
  rows <- Filter(function(z) !is.null(z) && nrow(z) > 0L, rows)
  if (!length(rows)) return(prototype)
  ans <- do.call(rbind, rows)
  rownames(ans) <- NULL
  ans
}

ds_arcs <- function(graph) {
  a <- bnlearn::arcs(graph)
  data.frame(from = as.character(a[, "from"]), to = as.character(a[, "to"]),
             stringsAsFactors = FALSE, check.names = FALSE)
}

ds_edge_key <- function(arcs) paste(arcs$from, arcs$to, sep = "\034")

ds_lag <- function(nodes) {
  out <- rep(NA_integer_, length(nodes))
  out[grepl("_t$", nodes)] <- 0L
  past <- grepl("_t-[1-9][0-9]*$", nodes)
  out[past] <- as.integer(sub(".*_t-([0-9]+)$", "\\1", nodes[past]))
  out
}

align_bn_to_nodes <- function(bn_str, required_nodes) {
  if (!inherits(bn_str, "bn")) stop("Expected a successful bnlearn graph.")
  ns <- bnlearn::nodes(bn_str)
  if (anyDuplicated(ns) || !setequal(ns, required_nodes)) stop("Graph/data node sets differ.")
  a <- ds_arcs(bn_str)
  if (nrow(a) && (anyNA(ds_lag(a$from)) || anyNA(ds_lag(a$to)) ||
                  any(ds_lag(a$from) < 1L | ds_lag(a$to) != 0L))) {
    stop("Graph contains an edge other than lagged-to-current.")
  }
  if (anyDuplicated(ds_edge_key(a))) stop("Graph contains duplicated arcs.")
  bn_str
}

ds_graph <- function(nodes, arc_df) {
  g <- bnlearn::empty.graph(nodes)
  if (nrow(arc_df)) bnlearn::arcs(g) <- as.matrix(arc_df[, c("from", "to"), drop = FALSE])
  align_bn_to_nodes(g, nodes)
}

ds_dot <- function(g, file) {
  quote_dot <- function(x) encodeString(x, quote = '"')
  a <- ds_arcs(g)
  writeLines(c("digraph {", paste0("  ", quote_dot(bnlearn::nodes(g)), ";"),
               if (nrow(a)) paste0("  ", quote_dot(a$from), " -> ", quote_dot(a$to), ";"),
               "}"), ds_path(file))
}

# Validate the discrete data without changing factor levels or column names.
input_rows <- list()
base_nodes <- NULL
base_levels <- NULL
for (p in candidate_lags) {
  pn <- paste0("p", p)
  D <- df_bn_by_lag[[pn]]
  if (!is.data.frame(D) || nrow(D) < 2L || !ncol(D)) stop("Missing/invalid BN data: ", pn)
  if (anyDuplicated(names(D)) || anyNA(names(D)) || anyNA(D) ||
      !all(vapply(D, is.factor, logical(1))) || any(vapply(D, nlevels, integer(1)) < 2L)) {
    stop("BN data must have unique names, factor columns with retained levels, and no NAs: ", pn)
  }
  current <- names(D)[grepl("_t$", names(D))]
  bases <- sub("_t$", "", current)
  if (is.null(base_nodes)) {
    base_nodes <- bases
    base_levels <- lapply(D[current], levels)
  }
  if (!identical(bases, base_nodes)) stop("Current node sets differ between lag specifications.")
  expected <- c(current, unlist(lapply(seq_len(p), function(k) paste0(bases, "_t-", k))))
  if (!setequal(names(D), expected)) stop("Incomplete or unexpected lag-expanded columns: ", pn)
  for (j in seq_along(bases)) {
    nn <- c(paste0(bases[j], "_t"), paste0(bases[j], "_t-", seq_len(p)))
    if (!all(vapply(D[nn], function(x) identical(levels(x), base_levels[[j]]), logical(1)))) {
      stop("Factor levels differ across lagged copies: ", bases[j])
    }
  }
  align_bn_to_nodes(empty_net_by_lag[[pn]], names(D))
  if (nrow(ds_arcs(empty_net_by_lag[[pn]]))) stop("empty_net_by_lag is not empty: ", pn)
  input_rows[[pn]] <- data.frame(lag = pn, model_max_lag = p, n_rows = nrow(D),
                                 n_nodes = ncol(D), n_current = length(current),
                                 n_lagged = length(expected) - length(current),
                                 selected = p == selected_p)
}
downstream_input_audit <- ds_bind(input_rows)
ds_csv(downstream_input_audit, "downstream_input_audit.csv")

infection_vars <- c("Reinfections", "New.infections", "New.cases")
mobility_vars <- c("Flights.7.day.moving.average",
                   "OpenTable.restaurant.bookings.London.index",
                   "Google.homeworking.Greater.London.mobility.index",
                   "Google.workplace.Greater.London.mobility.index", "Apple.walking.London.mobility.index",
                   "Google.parks.Greater.London.mobility.index",
                   "Google.retail.recreation.Greater.London.mobility.index",
                   "Google.grocery.pharmacy.Greater.London.mobility.index",
                   "Google.transit.stations.mobility.index", "TfL.Tube.mobility.index",
                   "TfL.Bus.mobility.index", "Citymapper.journeys.mobility.index")
make_arcs <- function(from, to) expand.grid(from = from, to = to, stringsAsFactors = FALSE)
knowledge_arcs <- unique(ds_bind(list(
  make_arcs("Majority.COVID.19.variant", c("Reinfections", "New.infections",
                                           "Hospital.admissions", "Patients.in.MVBs", "COVID.19.deaths.on.certificate",
                                           "Lockdown.severity", "Second.dose.uptake")),
  make_arcs("First.dose.uptake", "Second.dose.uptake"),
  make_arcs("Second.dose.uptake", "Third.dose.uptake"),
  make_arcs(c("Reinfections", "New.infections", "Tests.across.all.4.Pillars"), "New.cases"),
  make_arcs("New.cases", c("Hospital.admissions", "Second.dose.uptake")),
  make_arcs("Hospital.admissions", "Patients.in.hospital"),
  make_arcs("Patients.in.hospital", c("COVID.19.deaths.on.certificate",
                                      "Patients.in.MVBs", "Face.masks", "Lockdown.severity")),
  make_arcs("Patients.in.MVBs", c("COVID.19.deaths.on.certificate", "Face.masks", "Lockdown.severity")),
  make_arcs("COVID.19.deaths.on.certificate", c("Excess.mortality", "Face.masks", "Lockdown.severity")),
  make_arcs("Lockdown.severity", c("Schools", mobility_vars)),
  make_arcs("Season", c("Reinfections", "New.infections", "Lockdown.severity", mobility_vars))
)))

knowledge_arcs$min_lag <- 1L
knowledge_arcs$max_lag <- 21L

expand_prior_variable <- function(variable_name, factor_dummy_map) {
  ans <- if (variable_name %in% names(factor_dummy_map)) factor_dummy_map[[variable_name]] else variable_name
  if (!length(ans)) stop("No encoded columns for prior variable: ", variable_name)
  if (any(!ans %in% base_nodes)) stop("Unknown encoded prior variable: ", paste(ans, collapse = ", "))
  as.character(ans)
}

kb_by_lag <- list()
knowledge_expansion_rows <- list()

for (p in candidate_lags) {
  pn <- paste0("p", p)
  expanded <- lapply(seq_len(nrow(knowledge_arcs)), function(i) {
    z <- knowledge_arcs[i, ]
    hi <- min(p, z$max_lag)
    if (z$min_lag > hi) return(NULL)
    f <- expand_prior_variable(z$from, factor_dummy_map)
    t <- expand_prior_variable(z$to, factor_dummy_map)
    grid <- expand.grid(from_base = f, to_base = t, edge_lag = seq.int(z$min_lag, hi),
                        stringsAsFactors = FALSE)
    data.frame(lag = pn, prior_from = z$from, prior_to = z$to, grid,
               from = paste0(grid$from_base, "_t-", grid$edge_lag),
               to = paste0(grid$to_base, "_t"), stringsAsFactors = FALSE)
  })
  expanded <- ds_bind(expanded)
  a <- unique(expanded[, c("from", "to"), drop = FALSE])
  kb_by_lag[[pn]] <- ds_graph(names(df_bn_by_lag[[pn]]), a)
  knowledge_expansion_rows[[pn]] <- expanded
  ds_dot(kb_by_lag[[pn]], paste0("knowledge_graph_", pn, "_dot_file.dot"))
}
knowledge_expansion_audit <- ds_bind(knowledge_expansion_rows)
ds_csv(knowledge_arcs, "knowledge_base_arcs.csv")
ds_csv(knowledge_expansion_audit, "knowledge_expansion_audit.csv")
ds_rds(kb_by_lag, "kb_by_lag.rds")

ds_econ_methods <- c("lasso", "lar", "james_stein", "simone")
ds_ml_algorithms <- unique(c("pc.stable", "gs", "iamb", "fast.iamb", "inter.iamb",
                             "iamb.fdr", "hc", "tabu", "mmhc", "rsmax2", "h2pc", ds_get("causal_ml_algorithms")))
ds_methods <- c(ds_econ_methods, paste0("causal_ml_", ds_ml_algorithms), "model_average", "knowledge")
ds_entry <- function(pn, method) {
  fit <- NULL
  if (method %in% ds_econ_methods) {
    prefix <- if (method == "james_stein") "js" else method
    fit <- ds_get(paste0(prefix, "_fit_results_by_lag"))[[pn]]
    g <- ds_get(paste0(prefix, "_graphs_by_lag"))[[pn]]
  } else if (grepl("^causal_ml_", method)) {
    alg <- sub("^causal_ml_", "", method)
    fit <- ds_get("causal_ml_fit_results_by_lag")[[pn]][[alg]]
    source_graphs <- ds_get("causal_ml_graphs_by_lag")
    if (is.null(source_graphs)) source_graphs <- ds_get("converted_dags_by_lag")
    g <- source_graphs[[pn]][[alg]]
  } else if (method == "model_average") {
    g <- ds_get("average_bn_by_lag")[[pn]]
    fit <- ds_get("average_bn_fit_results_by_lag")[[pn]]
  } else if (method == "knowledge") {
    g <- ds_get("kb_by_lag")[[pn]]
  } else stop("Unknown method: ", method)
  detail <- if (!is.null(fit$error)) paste(fit$error, collapse = "; ") else NA_character_
  warns <- if (!is.null(fit$diagnostics$warnings)) paste(fit$diagnostics$warnings, collapse = "; ") else NA_character_
  if (!is.null(fit$status) && !fit$status %in% c("success", "success_no_edges")) {
    return(list(graph = NULL, status = fit$status, detail = detail, warnings = warns))
  }
  if (!is.null(fit) && "graph" %in% names(fit)) g <- fit$graph
  if (is.null(g)) return(list(graph = NULL, status = "missing_graph", detail = detail, warnings = warns))
  checked <- tryCatch(align_bn_to_nodes(g, names(df_bn_by_lag[[pn]])), error = identity)
  if (inherits(checked, "error")) return(list(graph = NULL, status = "invalid_graph",
                                              detail = conditionMessage(checked), warnings = warns))
  list(graph = checked, status = if (nrow(ds_arcs(checked))) "success" else "success_no_edges",
       detail = detail, warnings = warns)
}
ds_registry <- function(pn, methods = ds_methods) setNames(lapply(methods, function(m) ds_entry(pn, m)), methods)
downstream_version <- "2026-10-14-v1"
message("Script 9 complete. Selected lag retained: ", selected_name,
        ". Outputs: ", downstream_output_dir)
