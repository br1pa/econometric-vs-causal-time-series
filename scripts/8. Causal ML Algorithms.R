# 8. Causal ML Algorithms.R -----------------------------------------------------

suppressPackageStartupMessages({
  library(bnlearn)
  library(dplyr)
})

set.seed(123)


# ------------------------------------------------------------------------------
# 0. Configuration and checks
# ------------------------------------------------------------------------------

causal_ml_algorithms <- c(
  "pc.stable",
  "gs",
  "iamb",
  "fast.iamb",
  "inter.iamb",
  "iamb.fdr",
  "h2pc",
  "hc",
  "tabu",
  "mmhc",
  "rsmax2"
)

causal_ml_alpha <- 0.05
causal_ml_test <- "mi-g-sh"
causal_ml_score <- "ebic-g"
causal_ml_seed <- 123L


if (!exists("df_structure_by_lag")) {
  stop(
    "df_structure_by_lag does not exist. Run ",
    "4. parametrisation_processing.R before ",
    "8. Causal ML Algorithms.R."
  )
}

if (!exists("candidate_lags")) {
  stop("candidate_lags does not exist. Run 3. lag_strategy.R first.")
}

candidate_lags <- sort(unique(as.integer(candidate_lags)))

if (length(candidate_lags) == 0L || anyNA(candidate_lags) ||
    any(candidate_lags < 1L)) {
  stop("candidate_lags must contain positive integers.")
}


supported_causal_ml_algorithms <- c(
  "pc.stable",
  "gs",
  "iamb",
  "fast.iamb",
  "inter.iamb",
  "iamb.fdr",
  "h2pc",
  "hc",
  "tabu",
  "mmhc",
  "rsmax2"
)

unknown_algorithms <- setdiff(
  causal_ml_algorithms,
  supported_causal_ml_algorithms
)

if (length(unknown_algorithms) > 0L) {
  stop(
    "Unsupported causal-ML algorithms: ",
    paste(unknown_algorithms, collapse = ", ")
  )
}


if (!exists("new_fit_result", mode = "function") ||
    !exists("fit_succeeded", mode = "function")) {
  helper_candidates <- c(
    "4a.fit_result_helpers.R",
    "4a. fit_result_helpers.R"
  )
  
  helper_file <- helper_candidates[file.exists(helper_candidates)][1]
  
  if (is.na(helper_file)) {
    stop(
      "The common fit-result helpers are unavailable. Create and run ",
      "4a.fit_result_helpers.R before 8. Causal ML Algorithms.R."
    )
  }
  
  source(helper_file)
}


# ------------------------------------------------------------------------------
# 1. Utility functions
# ------------------------------------------------------------------------------

clean_condition_text <- function(x) {
  if (inherits(x, "condition")) {
    x <- conditionMessage(x)
  }
  
  x <- as.character(x)
  
  if (length(x) == 0L || is.na(x[1])) {
    return(NA_character_)
  }
  
  x <- gsub("[\r\n\t]+", " ", x[1])
  trimws(x)
}


causal_ml_scalar_or_default <- function(x, default) {
  if (is.null(x) || length(x) == 0L || is.na(x[1])) {
    return(default)
  }
  
  x[1]
}


node_time <- function(node) {
  node <- as.character(node)
  
  current <- grepl("_t$", node)
  lagged <- grepl("_t-[0-9]+$", node)
  
  if (any(!current & !lagged)) {
    stop(
      "Cannot identify the time slice for node(s): ",
      paste(node[!current & !lagged], collapse = ", ")
    )
  }
  
  out <- integer(length(node))
  out[current] <- 0L
  out[lagged] <- -as.integer(
    sub(".*_t-([0-9]+)$", "\\1", node[lagged])
  )
  
  out
}


causal_ml_algorithm_class <- function(algorithm) {
  if (algorithm %in% c(
    "pc.stable",
    "gs",
    "iamb",
    "fast.iamb",
    "inter.iamb",
    "iamb.fdr"
  )) {
    return("constraint_based")
  }
  
  if (algorithm %in% c("hc", "tabu")) {
    return("score_based")
  }
  
  if (algorithm %in% c("h2pc", "mmhc", "rsmax2")) {
    return("hybrid")
  }
  
  stop("Unknown algorithm: ", algorithm)
}


# ------------------------------------------------------------------------------
# 2. Prepare the continuous causal-ML data
# ------------------------------------------------------------------------------

prepare_causal_ml_data <- function(df_structure_p) {
  data <- as.data.frame(
    df_structure_p,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (nrow(data) < 2L) {
    stop("The causal-ML dataset contains fewer than two observations.")
  }
  
  if (ncol(data) < 2L) {
    stop("The causal-ML dataset contains fewer than two nodes.")
  }
  
  if (is.null(names(data)) || anyNA(names(data)) ||
      any(names(data) == "")) {
    stop("Every causal-ML column must have a non-missing name.")
  }
  
  if (anyDuplicated(names(data))) {
    stop("The causal-ML dataset contains duplicated node names.")
  }
  
  valid_temporal_name <- grepl(
    "(_t$)|(_t-[0-9]+$)",
    names(data)
  )
  
  if (!all(valid_temporal_name)) {
    stop(
      "Invalid temporal node names: ",
      paste(names(data)[!valid_temporal_name], collapse = ", ")
    )
  }
  
  current_nodes <- names(data)[grepl("_t$", names(data))]
  lagged_nodes <- names(data)[grepl("_t-[0-9]+$", names(data))]
  
  if (length(current_nodes) == 0L) {
    stop("No current-time nodes were found.")
  }
  
  if (length(lagged_nodes) == 0L) {
    stop("No lagged nodes were found.")
  }
  
  numeric_columns <- vapply(data, is.numeric, logical(1))
  
  if (!all(numeric_columns)) {
    stop(
      "All causal-ML columns must already be numeric. Non-numeric columns: ",
      paste(names(data)[!numeric_columns], collapse = ", ")
    )
  }
  
  matrix_data <- data.matrix(data)
  storage.mode(matrix_data) <- "double"
  
  if (anyNA(matrix_data) || any(!is.finite(matrix_data))) {
    stop("The causal-ML input contains missing or non-finite values.")
  }
  
  minimum <- apply(matrix_data, 2L, min)
  maximum <- apply(matrix_data, 2L, max)
  range_value <- maximum - minimum
  constant <- !is.finite(range_value) | range_value <= 0
  
  safe_range <- range_value
  safe_range[constant] <- 1
  
  scaled_matrix <- sweep(matrix_data, 2L, minimum, FUN = "-")
  scaled_matrix <- sweep(scaled_matrix, 2L, safe_range, FUN = "/")
  
  if (any(constant)) {
    scaled_matrix[, constant] <- 0
  }
  
  if (anyNA(scaled_matrix) || any(!is.finite(scaled_matrix))) {
    stop("Min-max scaling produced missing or non-finite values.")
  }
  
  scaled_data <- as.data.frame(
    scaled_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  names(scaled_data) <- names(data)
  
  node_audit <- data.frame(
    node = names(data),
    time = node_time(names(data)),
    minimum = as.numeric(minimum),
    maximum = as.numeric(maximum),
    range = as.numeric(range_value),
    constant = as.logical(constant),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  list(
    data = scaled_data,
    node_audit = node_audit,
    n_observations = nrow(scaled_data),
    n_nodes = ncol(scaled_data),
    n_current_nodes = length(current_nodes),
    n_lagged_nodes = length(lagged_nodes),
    n_constant_nodes = sum(constant),
    constant_nodes = names(data)[constant]
  )
}


# ------------------------------------------------------------------------------
# 3. Temporal blacklist
# ------------------------------------------------------------------------------

make_temporal_blacklist <- function(nodes) {
  nodes <- as.character(nodes)
  
  if (length(nodes) < 2L) {
    stop("At least two nodes are required to construct a blacklist.")
  }
  
  if (anyDuplicated(nodes)) {
    stop("Cannot construct a blacklist from duplicated node names.")
  }
  
  node_times <- node_time(nodes)
  current_indicator <- node_times == 0L
  lagged_indicator <- node_times < 0L
  
  if (!all(current_indicator | lagged_indicator)) {
    stop("Every node must be either lagged or current.")
  }
  
  n_nodes <- length(nodes)
  
  from_index <- rep(seq_len(n_nodes), each = n_nodes)
  to_index <- rep(seq_len(n_nodes), times = n_nodes)
  
  non_self <- from_index != to_index
  allowed <- lagged_indicator[from_index] & current_indicator[to_index]
  blacklist_indicator <- non_self & !allowed
  
  blacklist <- cbind(
    from = nodes[from_index[blacklist_indicator]],
    to = nodes[to_index[blacklist_indicator]]
  )
  
  storage.mode(blacklist) <- "character"
  
  expected_allowed <- sum(lagged_indicator) * sum(current_indicator)
  expected_blacklisted <- n_nodes * (n_nodes - 1L) - expected_allowed
  
  if (nrow(blacklist) != expected_blacklisted) {
    stop("The temporal blacklist has an unexpected number of rows.")
  }
  
  if (any(blacklist[, "from"] == blacklist[, "to"])) {
    stop("The temporal blacklist unexpectedly contains self-arcs.")
  }
  
  blacklist
}


# ------------------------------------------------------------------------------
# 4. bnlearn algorithm calls
# ------------------------------------------------------------------------------

run_bnlearn_causal_algorithm <- function(algorithm,
                                         dataframe,
                                         blacklist,
                                         alpha_value,
                                         test_name,
                                         score_name) {
  switch(
    algorithm,
    
    "pc.stable" = bnlearn::pc.stable(
      x = dataframe,
      alpha = alpha_value,
      test = test_name,
      blacklist = blacklist,
      debug = FALSE
    ),
    
    "gs" = bnlearn::gs(
      x = dataframe,
      alpha = alpha_value,
      test = test_name,
      blacklist = blacklist,
      debug = FALSE
    ),
    
    "iamb" = bnlearn::iamb(
      x = dataframe,
      alpha = alpha_value,
      test = test_name,
      blacklist = blacklist,
      debug = FALSE
    ),
    
    "fast.iamb" = bnlearn::fast.iamb(
      x = dataframe,
      alpha = alpha_value,
      test = test_name,
      blacklist = blacklist,
      debug = FALSE
    ),
    
    "inter.iamb" = bnlearn::inter.iamb(
      x = dataframe,
      alpha = alpha_value,
      test = test_name,
      blacklist = blacklist,
      debug = FALSE
    ),
    
    "iamb.fdr" = bnlearn::iamb.fdr(
      x = dataframe,
      alpha = alpha_value,
      test = test_name,
      blacklist = blacklist,
      debug = FALSE
    ),
    
    "h2pc" = bnlearn::h2pc(
      x = dataframe,
      blacklist = blacklist,
      restrict.args = list(
        alpha = alpha_value,
        test = test_name
      ),
      maximize.args = list(
        score = score_name
      ),
      debug = FALSE
    ),
    
    "hc" = bnlearn::hc(
      x = dataframe,
      score = score_name,
      blacklist = blacklist,
      debug = FALSE
    ),
    
    "tabu" = bnlearn::tabu(
      x = dataframe,
      score = score_name,
      blacklist = blacklist,
      debug = FALSE
    ),
    
    "mmhc" = bnlearn::mmhc(
      x = dataframe,
      blacklist = blacklist,
      restrict.args = list(
        alpha = alpha_value,
        test = test_name
      ),
      maximize.args = list(
        score = score_name
      ),
      debug = FALSE
    ),
    
    "rsmax2" = bnlearn::rsmax2(
      x = dataframe,
      blacklist = blacklist,
      restrict.args = list(
        alpha = alpha_value,
        test = test_name
      ),
      maximize.args = list(
        score = score_name
      ),
      debug = FALSE
    ),
    
    stop("Unknown algorithm: ", algorithm)
  )
}


# ------------------------------------------------------------------------------
# 5. Convert a learnt structure into the common temporal DAG
# ------------------------------------------------------------------------------

convert_raw_bn_to_temporal_dag <- function(raw_graph, required_nodes) {
  if (is.null(raw_graph) || !inherits(raw_graph, "bn")) {
    stop("The learned object is not a bnlearn graph.")
  }
  
  raw_nodes <- bnlearn::nodes(raw_graph)
  
  if (length(raw_nodes) != length(required_nodes) ||
      !setequal(raw_nodes, required_nodes)) {
    stop("The learned graph does not contain the required node set.")
  }
  
  raw_arc_matrix <- bnlearn::arcs(raw_graph)
  
  if (is.null(raw_arc_matrix) || nrow(raw_arc_matrix) == 0L) {
    raw_arcs <- data.frame(
      from = character(0),
      to = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    raw_arcs <- as.data.frame(
      raw_arc_matrix,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    names(raw_arcs) <- c("from", "to")
  }
  
  if (nrow(raw_arcs) == 0L) {
    return(list(
      graph = bnlearn::empty.graph(required_nodes),
      raw_arcs = raw_arcs,
      n_raw_arc_rows = 0L,
      n_direct_temporal_arc_rows = 0L,
      n_reverse_temporal_arc_rows = 0L,
      n_non_temporal_arc_rows = 0L,
      n_duplicate_temporal_rows_collapsed = 0L,
      n_final_edges = 0L
    ))
  }
  
  if (any(!raw_arcs$from %in% required_nodes) ||
      any(!raw_arcs$to %in% required_nodes)) {
    stop("The learned graph contains arcs involving unknown nodes.")
  }
  
  from_time <- node_time(raw_arcs$from)
  to_time <- node_time(raw_arcs$to)
  
  direct_temporal <- from_time < 0L & to_time == 0L
  reverse_temporal <- from_time == 0L & to_time < 0L
  temporal_adjacency <- direct_temporal | reverse_temporal
  
  temporal_arcs <- data.frame(
    from = ifelse(
      direct_temporal[temporal_adjacency],
      raw_arcs$from[temporal_adjacency],
      raw_arcs$to[temporal_adjacency]
    ),
    to = ifelse(
      direct_temporal[temporal_adjacency],
      raw_arcs$to[temporal_adjacency],
      raw_arcs$from[temporal_adjacency]
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  n_temporal_rows_before_deduplication <- nrow(temporal_arcs)
  temporal_arcs <- unique(temporal_arcs)
  
  final_graph <- bnlearn::empty.graph(required_nodes)
  
  if (nrow(temporal_arcs) > 0L) {
    temporal_arc_matrix <- as.matrix(
      temporal_arcs[, c("from", "to"), drop = FALSE]
    )
    colnames(temporal_arc_matrix) <- c("from", "to")
    bnlearn::arcs(final_graph) <- temporal_arc_matrix
  }
  
  final_arc_matrix <- bnlearn::arcs(final_graph)
  
  if (!is.null(final_arc_matrix) && nrow(final_arc_matrix) > 0L) {
    valid_order <- node_time(final_arc_matrix[, "from"]) < 0L &
      node_time(final_arc_matrix[, "to"]) == 0L
    
    if (!all(valid_order)) {
      stop("The converted graph violates the required temporal ordering.")
    }
  }
  
  list(
    graph = final_graph,
    raw_arcs = raw_arcs,
    n_raw_arc_rows = nrow(raw_arcs),
    n_direct_temporal_arc_rows = sum(direct_temporal),
    n_reverse_temporal_arc_rows = sum(reverse_temporal),
    n_non_temporal_arc_rows = sum(!temporal_adjacency),
    n_duplicate_temporal_rows_collapsed =
      n_temporal_rows_before_deduplication - nrow(temporal_arcs),
    n_final_edges = nrow(temporal_arcs)
  )
}


clean_temporal_dag <- function(bn_structure, nodes) {
  convert_raw_bn_to_temporal_dag(
    raw_graph = bn_structure,
    required_nodes = nodes
  )$graph
}


convert_to_temporal_dag <- function(bn_structure, nodes) {
  if (is.null(bn_structure)) {
    return(NULL)
  }
  
  clean_temporal_dag(
    bn_structure = bn_structure,
    nodes = nodes
  )
}


has_any_edges <- function(bn_structure) {
  if (is.null(bn_structure) || !inherits(bn_structure, "bn")) {
    return(FALSE)
  }
  
  arc_matrix <- bnlearn::arcs(bn_structure)
  !is.null(arc_matrix) && nrow(arc_matrix) > 0L
}


# ------------------------------------------------------------------------------
# 6. Fit-result construction
# ------------------------------------------------------------------------------

failed_causal_ml_outcome <- function(algorithm,
                                     p,
                                     stage,
                                     error,
                                     rng_seed = NA_integer_,
                                     elapsed_seconds = NA_real_,
                                     warnings = character(0),
                                     preparation = NULL,
                                     n_blacklisted_arcs = NA_integer_) {
  n_observations <- if (is.null(preparation)) {
    NA_integer_
  } else {
    preparation$n_observations
  }
  
  n_nodes <- if (is.null(preparation)) {
    NA_integer_
  } else {
    preparation$n_nodes
  }
  
  n_constant_nodes <- if (is.null(preparation)) {
    NA_integer_
  } else {
    preparation$n_constant_nodes
  }
  
  result <- new_fit_result(
    status = "failed",
    error = clean_condition_text(error),
    graph = NULL,
    coefficients = NULL,
    diagnostics = list(
      method = "causal_ml",
      algorithm = algorithm,
      algorithm_class = causal_ml_algorithm_class(algorithm),
      p = as.integer(p),
      failed_stage = stage,
      rng_seed = as.integer(rng_seed),
      elapsed_seconds = as.numeric(elapsed_seconds),
      warnings = unique(as.character(warnings)),
      bnlearn_version = as.character(utils::packageVersion("bnlearn")),
      test = if (causal_ml_algorithm_class(algorithm) == "score_based") {
        NA_character_
      } else {
        causal_ml_test
      },
      alpha = if (causal_ml_algorithm_class(algorithm) == "score_based") {
        NA_real_
      } else {
        causal_ml_alpha
      },
      score = if (causal_ml_algorithm_class(algorithm) ==
                  "constraint_based") {
        NA_character_
      } else {
        causal_ml_score
      },
      n_observations = as.integer(n_observations),
      n_nodes = as.integer(n_nodes),
      n_constant_nodes = as.integer(n_constant_nodes),
      n_blacklisted_arcs = as.integer(n_blacklisted_arcs),
      n_raw_arc_rows = NA_integer_,
      n_direct_temporal_arc_rows = NA_integer_,
      n_reverse_temporal_arc_rows = NA_integer_,
      n_non_temporal_arc_rows = NA_integer_,
      n_duplicate_temporal_rows_collapsed = NA_integer_,
      n_final_edges = NA_integer_
    )
  )
  
  list(
    result = result,
    raw_arcs = NULL
  )
}


fit_causal_ml_algorithm <- function(algorithm,
                                    dataframe,
                                    blacklist,
                                    p,
                                    preparation,
                                    alpha_value,
                                    test_name,
                                    score_name,
                                    rng_seed) {
  stage <- "initialisation"
  captured_warnings <- character(0)
  start_time <- unname(proc.time()[3L])
  
  tryCatch({
    set.seed(as.integer(rng_seed))
    
    stage <- "structure_learning"
    
    raw_graph <- withCallingHandlers(
      run_bnlearn_causal_algorithm(
        algorithm = algorithm,
        dataframe = dataframe,
        blacklist = blacklist,
        alpha_value = alpha_value,
        test_name = test_name,
        score_name = score_name
      ),
      warning = function(w) {
        captured_warnings <<- c(
          captured_warnings,
          clean_condition_text(w)
        )
        tryInvokeRestart("muffleWarning")
      }
    )
    
    stage <- "learned_object_validation"
    
    if (is.null(raw_graph) || !inherits(raw_graph, "bn")) {
      stop("The algorithm did not return a bnlearn graph.")
    }
    
    stage <- "temporal_conversion"
    
    conversion <- convert_raw_bn_to_temporal_dag(
      raw_graph = raw_graph,
      required_nodes = colnames(dataframe)
    )
    
    stage <- "final_graph_validation"
    
    final_graph <- conversion$graph
    
    if (is.null(final_graph) || !inherits(final_graph, "bn")) {
      stop("Temporal conversion did not return a bnlearn graph.")
    }
    
    final_nodes <- bnlearn::nodes(final_graph)
    
    if (length(final_nodes) != ncol(dataframe) ||
        !setequal(final_nodes, colnames(dataframe))) {
      stop("The final temporal graph has an incorrect node set.")
    }
    
    elapsed_seconds <- unname(proc.time()[3L]) - start_time
    
    status <- if (conversion$n_final_edges == 0L) {
      "success_no_edges"
    } else {
      "success"
    }
    
    result <- new_fit_result(
      status = status,
      error = NA_character_,
      graph = final_graph,
      coefficients = NULL,
      diagnostics = list(
        method = "causal_ml",
        algorithm = algorithm,
        algorithm_class = causal_ml_algorithm_class(algorithm),
        p = as.integer(p),
        failed_stage = NA_character_,
        rng_seed = as.integer(rng_seed),
        elapsed_seconds = as.numeric(elapsed_seconds),
        warnings = unique(captured_warnings),
        bnlearn_version = as.character(utils::packageVersion("bnlearn")),
        test = if (causal_ml_algorithm_class(algorithm) == "score_based") {
          NA_character_
        } else {
          test_name
        },
        alpha = if (causal_ml_algorithm_class(algorithm) == "score_based") {
          NA_real_
        } else {
          alpha_value
        },
        score = if (causal_ml_algorithm_class(algorithm) ==
                    "constraint_based") {
          NA_character_
        } else {
          score_name
        },
        n_observations = preparation$n_observations,
        n_nodes = preparation$n_nodes,
        n_constant_nodes = preparation$n_constant_nodes,
        constant_nodes = preparation$constant_nodes,
        n_blacklisted_arcs = nrow(blacklist),
        n_raw_arc_rows = conversion$n_raw_arc_rows,
        n_direct_temporal_arc_rows =
          conversion$n_direct_temporal_arc_rows,
        n_reverse_temporal_arc_rows =
          conversion$n_reverse_temporal_arc_rows,
        n_non_temporal_arc_rows =
          conversion$n_non_temporal_arc_rows,
        n_duplicate_temporal_rows_collapsed =
          conversion$n_duplicate_temporal_rows_collapsed,
        n_final_edges = conversion$n_final_edges
      )
    )
    
    list(
      result = result,
      raw_arcs = conversion$raw_arcs
    )
  }, error = function(e) {
    elapsed_seconds <- unname(proc.time()[3L]) - start_time
    
    failed_causal_ml_outcome(
      algorithm = algorithm,
      p = p,
      stage = stage,
      error = e,
      rng_seed = rng_seed,
      elapsed_seconds = elapsed_seconds,
      warnings = captured_warnings,
      preparation = preparation,
      n_blacklisted_arcs = nrow(blacklist)
    )
  })
}


# ------------------------------------------------------------------------------
# 7. Status and output helpers
# ------------------------------------------------------------------------------

causal_ml_status_row <- function(result, lag_name) {
  diagnostics <- result$diagnostics
  warning_values <- diagnostics$warnings
  
  warning_text <- if (is.null(warning_values) ||
                      length(warning_values) == 0L) {
    NA_character_
  } else {
    paste(unique(warning_values), collapse = " | ")
  }
  
  data.frame(
    lag = lag_name,
    p = as.integer(
      causal_ml_scalar_or_default(diagnostics$p, NA_integer_)
    ),
    method = "causal_ml",
    algorithm = causal_ml_scalar_or_default(
      diagnostics$algorithm,
      NA_character_
    ),
    algorithm_class = causal_ml_scalar_or_default(
      diagnostics$algorithm_class,
      NA_character_
    ),
    status = causal_ml_scalar_or_default(result$status, "failed"),
    failed_stage = causal_ml_scalar_or_default(
      diagnostics$failed_stage,
      NA_character_
    ),
    error = clean_condition_text(result$error),
    warnings = warning_text,
    rng_seed = as.integer(
      causal_ml_scalar_or_default(diagnostics$rng_seed, NA_integer_)
    ),
    elapsed_seconds = as.numeric(
      causal_ml_scalar_or_default(
        diagnostics$elapsed_seconds,
        NA_real_
      )
    ),
    n_observations = as.integer(
      causal_ml_scalar_or_default(
        diagnostics$n_observations,
        NA_integer_
      )
    ),
    n_nodes = as.integer(
      causal_ml_scalar_or_default(diagnostics$n_nodes, NA_integer_)
    ),
    n_constant_nodes = as.integer(
      causal_ml_scalar_or_default(
        diagnostics$n_constant_nodes,
        NA_integer_
      )
    ),
    n_blacklisted_arcs = as.integer(
      causal_ml_scalar_or_default(
        diagnostics$n_blacklisted_arcs,
        NA_integer_
      )
    ),
    n_raw_arc_rows = as.integer(
      causal_ml_scalar_or_default(
        diagnostics$n_raw_arc_rows,
        NA_integer_
      )
    ),
    n_direct_temporal_arc_rows = as.integer(
      causal_ml_scalar_or_default(
        diagnostics$n_direct_temporal_arc_rows,
        NA_integer_
      )
    ),
    n_reverse_temporal_arc_rows = as.integer(
      causal_ml_scalar_or_default(
        diagnostics$n_reverse_temporal_arc_rows,
        NA_integer_
      )
    ),
    n_non_temporal_arc_rows = as.integer(
      causal_ml_scalar_or_default(
        diagnostics$n_non_temporal_arc_rows,
        NA_integer_
      )
    ),
    n_duplicate_temporal_rows_collapsed = as.integer(
      causal_ml_scalar_or_default(
        diagnostics$n_duplicate_temporal_rows_collapsed,
        NA_integer_
      )
    ),
    n_final_edges = as.integer(
      causal_ml_scalar_or_default(
        diagnostics$n_final_edges,
        NA_integer_
      )
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


write_causal_ml_algorithm_outputs <- function(outcome,
                                              lag_name,
                                              algorithm) {
  result <- outcome$result
  
  if (!is.null(outcome$raw_arcs)) {
    utils::write.csv(
      outcome$raw_arcs,
      file = paste0(
        "causal_ml_raw_arcs_",
        algorithm,
        "_",
        lag_name,
        ".csv"
      ),
      row.names = FALSE,
      na = "NA"
    )
  }
  
  if (fit_succeeded(result)) {
    bnlearn::write.dot(
      result$graph,
      file = paste0(
        "causal_ml_",
        algorithm,
        "_",
        lag_name,
        "_dot_file.dot"
      )
    )
  }
}


# ------------------------------------------------------------------------------
# 8. Run every algorithm for every lag specification
# ------------------------------------------------------------------------------

causal_ml_fit_results_by_lag <- list()

causal_ml_graphs_by_lag <- list()
converted_dags_by_lag <- list()

causal_ml_status_rows <- list()
causal_ml_preprocessing_rows <- list()
causal_ml_node_audit_by_lag <- list()
causal_ml_blacklist_summary_rows <- list()


for (p in candidate_lags) {
  p_name <- paste0("p", p)
  
  message("============================================================")
  message("Running causal ML for lag specification ", p_name)
  message("============================================================")
  
  results_p <- setNames(
    vector("list", length(causal_ml_algorithms)),
    causal_ml_algorithms
  )
  
  graphs_p <- setNames(
    vector("list", length(causal_ml_algorithms)),
    causal_ml_algorithms
  )
  
  preparation <- NULL
  blacklist <- NULL
  
  preparation_error <- NULL
  constraint_error <- NULL
  
  if (!p_name %in% names(df_structure_by_lag)) {
    preparation_error <- simpleError(
      paste0("No structure-learning dataset exists for ", p_name, ".")
    )
  } else {
    preparation <- tryCatch(
      prepare_causal_ml_data(df_structure_by_lag[[p_name]]),
      error = function(e) {
        preparation_error <<- e
        NULL
      }
    )
  }
  
  if (is.null(preparation_error) && !is.null(preparation)) {
    blacklist <- tryCatch(
      make_temporal_blacklist(colnames(preparation$data)),
      error = function(e) {
        constraint_error <<- e
        NULL
      }
    )
  }
  
  if (!is.null(preparation)) {
    node_audit_p <- preparation$node_audit
    node_audit_p$lag <- p_name
    node_audit_p$p <- as.integer(p)
    node_audit_p <- node_audit_p[
      ,
      c(
        "lag",
        "p",
        "node",
        "time",
        "minimum",
        "maximum",
        "range",
        "constant"
      ),
      drop = FALSE
    ]
    
    causal_ml_node_audit_by_lag[[p_name]] <- node_audit_p
    
    utils::write.csv(
      node_audit_p,
      file = paste0("causal_ml_node_audit_", p_name, ".csv"),
      row.names = FALSE,
      na = "NA"
    )
    
    causal_ml_preprocessing_rows[[p_name]] <- data.frame(
      lag = p_name,
      p = as.integer(p),
      n_observations = preparation$n_observations,
      n_nodes = preparation$n_nodes,
      n_current_nodes = preparation$n_current_nodes,
      n_lagged_nodes = preparation$n_lagged_nodes,
      n_constant_nodes = preparation$n_constant_nodes,
      constant_nodes = if (length(preparation$constant_nodes) == 0L) {
        NA_character_
      } else {
        paste(preparation$constant_nodes, collapse = " | ")
      },
      stringsAsFactors = FALSE
    )
  } else {
    causal_ml_node_audit_by_lag[[p_name]] <- data.frame()
    
    causal_ml_preprocessing_rows[[p_name]] <- data.frame(
      lag = p_name,
      p = as.integer(p),
      n_observations = NA_integer_,
      n_nodes = NA_integer_,
      n_current_nodes = NA_integer_,
      n_lagged_nodes = NA_integer_,
      n_constant_nodes = NA_integer_,
      constant_nodes = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  
  causal_ml_blacklist_summary_rows[[p_name]] <- data.frame(
    lag = p_name,
    p = as.integer(p),
    n_blacklisted_arcs = if (is.null(blacklist)) {
      NA_integer_
    } else {
      nrow(blacklist)
    },
    n_allowed_lagged_to_current_arcs = if (is.null(preparation)) {
      NA_integer_
    } else {
      preparation$n_lagged_nodes * preparation$n_current_nodes
    },
    stringsAsFactors = FALSE
  )
  
  for (algorithm_index in seq_along(causal_ml_algorithms)) {
    algorithm <- causal_ml_algorithms[algorithm_index]
    algorithm_seed <- as.integer(
      causal_ml_seed + p * 100L + algorithm_index
    )
    
    message("Learning structure using ", algorithm, " at ", p_name)
    
    if (!is.null(preparation_error)) {
      outcome <- failed_causal_ml_outcome(
        algorithm = algorithm,
        p = p,
        stage = "data_preparation",
        error = preparation_error,
        rng_seed = algorithm_seed,
        preparation = preparation
      )
    } else if (!is.null(constraint_error)) {
      outcome <- failed_causal_ml_outcome(
        algorithm = algorithm,
        p = p,
        stage = "temporal_blacklist",
        error = constraint_error,
        rng_seed = algorithm_seed,
        preparation = preparation
      )
    } else {
      outcome <- fit_causal_ml_algorithm(
        algorithm = algorithm,
        dataframe = preparation$data,
        blacklist = blacklist,
        p = p,
        preparation = preparation,
        alpha_value = causal_ml_alpha,
        test_name = causal_ml_test,
        score_name = causal_ml_score,
        rng_seed = algorithm_seed
      )
    }
    
    results_p[[algorithm]] <- outcome$result
    
    graphs_p[algorithm] <- list(
      if (fit_succeeded(outcome$result)) {
        outcome$result$graph
      } else {
        NULL
      }
    )
    
    status_key <- paste(p_name, algorithm, sep = "__")
    causal_ml_status_rows[[status_key]] <- causal_ml_status_row(
      result = outcome$result,
      lag_name = p_name
    )
    
    write_causal_ml_algorithm_outputs(
      outcome = outcome,
      lag_name = p_name,
      algorithm = algorithm
    )
    
    message(
      "Completed ",
      algorithm,
      " at ",
      p_name,
      " with status: ",
      outcome$result$status
    )
    
    rm(outcome)
    invisible(gc(verbose = FALSE))
  }
  
  causal_ml_fit_results_by_lag[[p_name]] <- results_p
  converted_dags_by_lag[[p_name]] <- graphs_p
  causal_ml_graphs_by_lag[[p_name]] <- graphs_p
  
  status_p <- dplyr::bind_rows(
    causal_ml_status_rows[
      grepl(paste0("^", p_name, "__"), names(causal_ml_status_rows))
    ]
  )
  
  utils::write.csv(
    status_p,
    file = paste0("causal_ml_status_", p_name, ".csv"),
    row.names = FALSE,
    na = "NA"
  )
  
  rm(blacklist)
  invisible(gc(verbose = FALSE))
}


# ------------------------------------------------------------------------------
# 9. Combined audit files and saved R objects
# ------------------------------------------------------------------------------

causal_ml_status_by_lag <- dplyr::bind_rows(causal_ml_status_rows)
causal_ml_preprocessing_by_lag <- dplyr::bind_rows(
  causal_ml_preprocessing_rows
)
causal_ml_blacklist_summary_by_lag <- dplyr::bind_rows(
  causal_ml_blacklist_summary_rows
)


utils::write.csv(
  causal_ml_status_by_lag,
  file = "causal_ml_status_by_lag.csv",
  row.names = FALSE,
  na = "NA"
)

utils::write.csv(
  causal_ml_preprocessing_by_lag,
  file = "causal_ml_preprocessing_by_lag.csv",
  row.names = FALSE,
  na = "NA"
)

utils::write.csv(
  causal_ml_blacklist_summary_by_lag,
  file = "causal_ml_blacklist_summary_by_lag.csv",
  row.names = FALSE,
  na = "NA"
)


causal_ml_run_configuration <- data.frame(
  setting = c(
    "algorithms",
    "alpha",
    "conditional_independence_test",
    "score",
    "base_seed",
    "scaling",
    "permitted_arcs",
    "failed_graph_representation",
    "constraint_output_orientation"
  ),
  value = c(
    paste(causal_ml_algorithms, collapse = " | "),
    as.character(causal_ml_alpha),
    causal_ml_test,
    causal_ml_score,
    as.character(causal_ml_seed),
    "full-sample min-max scaling",
    "lagged nodes to current nodes only",
    "NULL",
    "lagged-to-current temporal orientation"
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  causal_ml_run_configuration,
  file = "causal_ml_run_configuration.csv",
  row.names = FALSE
)


saveRDS(
  causal_ml_fit_results_by_lag,
  file = "causal_ml_fit_results_by_lag.rds"
)

saveRDS(
  converted_dags_by_lag,
  file = "converted_dags_by_lag.rds"
)

saveRDS(
  causal_ml_graphs_by_lag,
  file = "causal_ml_graphs_by_lag.rds"
)


message("Causal-ML structure learning is complete.")
print(causal_ml_status_by_lag)
