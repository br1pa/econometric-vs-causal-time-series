# 7. simone.R -------------------------------------------------------------------

library(igraph)

if (!exists("is.directed", mode = "function")) {
  is.directed <- igraph::is_directed
}

library(simone)
library(bnlearn)
library(dplyr)

if (!exists("new_fit_result", mode = "function") ||
    !exists("failed_fit", mode = "function") ||
    !exists("fit_succeeded", mode = "function")) {
  
  source("4a. Helpers.R")
}

set.seed(123)


# ------------------------------------------------------------------------------
# 0. Checks and settings
# ------------------------------------------------------------------------------

if (!exists("df_structure_by_lag")) {
  stop(
    "df_structure_by_lag does not exist. ",
    "Run 4. parametrisation_processing.R before 7. simone.R."
  )
}

if (!exists("candidate_lags")) {
  stop(
    "candidate_lags does not exist. ",
    "Run 3. lag_strategy.R before 7. simone.R."
  )
}

if (!exists("fit_minmax_scaler", mode = "function") ||
    !exists("apply_minmax_scaler", mode = "function")) {
  
  stop(
    "Min-max scaling functions do not exist. ",
    "Run 5. lars.R first."
  )
}

if (!exists("selected_p")) {
  selected_p <- 7L
}

simone_results_threshold <- 0

# ------------------------------------------------------------------------------
# 1. Convert a coefficient matrix into a temporal DAG
# ------------------------------------------------------------------------------

generate_simone_graph_from_lagged_coefs <- function(
    coef_df,
    threshold = 0
) {
  
  if (is.null(coef_df)) {
    stop("coef_df is NULL.")
  }
  
  if (!is.data.frame(coef_df) &&
      !is.matrix(coef_df)) {
    stop("coef_df must be a data frame or matrix.")
  }
  
  if (is.null(rownames(coef_df)) ||
      is.null(colnames(coef_df))) {
    stop(
      "coef_df must have target row names ",
      "and predictor column names."
    )
  }
  
  if (!is.numeric(threshold) ||
      length(threshold) != 1L ||
      !is.finite(threshold) ||
      threshold < 0) {
    stop("threshold must be one non-negative number.")
  }
  
  coefficient_matrix <- as.matrix(coef_df)
  storage.mode(coefficient_matrix) <- "double"
  
  predictor_nodes <- colnames(coefficient_matrix)
  target_nodes <- rownames(coefficient_matrix)
  all_nodes <- unique(c(predictor_nodes, target_nodes))
  
  graph <- bnlearn::empty.graph(all_nodes)
  
  selected_positions <- which(
    is.finite(coefficient_matrix) &
      abs(coefficient_matrix) > threshold,
    arr.ind = TRUE
  )
  
  if (nrow(selected_positions) == 0) {
    return(graph)
  }
  
  selected_arcs <- cbind(
    from = predictor_nodes[selected_positions[, "col"]],
    to = target_nodes[selected_positions[, "row"]]
  )
  
  selected_arcs <- unique(selected_arcs)
  
  arcs(graph) <- selected_arcs
  
  return(graph)
}


# ------------------------------------------------------------------------------
# 2. Map numeric SIMONE node names to input column names
# ------------------------------------------------------------------------------

map_simone_node <- function(
    node,
    available_cols
) {
  
  node <- as.character(node)
  
  node_index <- suppressWarnings(
    as.integer(node)
  )
  
  if (!is.na(node_index) &&
      node_index >= 1L &&
      node_index <= length(available_cols)) {
    
    return(
      available_cols[node_index]
    )
  }
  
  return(node)
}


# ------------------------------------------------------------------------------
# 3. Apply SIMONE to one lag-expanded DBN dataset
# ------------------------------------------------------------------------------

simone_apply_structure_data <- function(
    df_structure_p,
    p,
    clustering = TRUE
) {
  
  stage <- "input_validation"
  
  tryCatch({
    
    # --------------------------------------------------------------------------
    # 3.1 Validate inputs
    # --------------------------------------------------------------------------
    
    if (!is.data.frame(df_structure_p)) {
      stop("df_structure_p must be a data frame.")
    }
    
    if (nrow(df_structure_p) == 0 ||
        ncol(df_structure_p) == 0) {
      stop("df_structure_p is empty.")
    }
    
    if (anyDuplicated(names(df_structure_p))) {
      stop("df_structure_p contains duplicated column names.")
    }
    
    if (!is.numeric(p) ||
        length(p) != 1L ||
        !is.finite(p) ||
        p < 1) {
      stop("p must be one positive integer.")
    }
    
    p <- as.integer(p)
    
    if (!is.logical(clustering) ||
        length(clustering) != 1L ||
        is.na(clustering)) {
      stop("clustering must be either TRUE or FALSE.")
    }
    
    # --------------------------------------------------------------------------
    # 3.2 Identify the lagged predictors and current outcomes
    # --------------------------------------------------------------------------
    
    stage <- "identify_nodes"
    
    x_cols <- grep(
      "_t-[0-9]+$",
      names(df_structure_p),
      value = TRUE
    )
    
    y_cols <- grep(
      "_t$",
      names(df_structure_p),
      value = TRUE
    )
    
    if (length(x_cols) == 0) {
      stop(
        "No lagged predictor columns found. ",
        "Expected names ending in _t-1, _t-2, etc."
      )
    }
    
    if (length(y_cols) == 0) {
      stop(
        "No current outcome columns found. ",
        "Expected names ending in _t."
      )
    }
    
    if (length(intersect(x_cols, y_cols)) > 0) {
      stop(
        "Some variables were classified as both lagged ",
        "predictors and current outcomes."
      )
    }
    
    all_cols <- c(x_cols, y_cols)
    
    # --------------------------------------------------------------------------
    # 3.3 Prepare and scale the data
    # --------------------------------------------------------------------------
    
    stage <- "prepare_data"
    
    df_simone <- df_structure_p[
      ,
      all_cols,
      drop = FALSE
    ]
    
    simone_matrix <- data.matrix(df_simone)
    storage.mode(simone_matrix) <- "double"
    
    if (anyNA(simone_matrix)) {
      stop(
        "The SIMONE input contains missing values ",
        "before preprocessing."
      )
    }
    
    if (any(!is.finite(simone_matrix))) {
      stop(
        "The SIMONE input contains non-finite values ",
        "before preprocessing."
      )
    }
    
    simone_scaler <- fit_minmax_scaler(
      simone_matrix
    )
    
    simone_matrix <- apply_minmax_scaler(
      simone_matrix,
      simone_scaler
    )
    
    df_simone <- as.data.frame(
      simone_matrix
    )
    
    names(df_simone) <- all_cols
    
    if (anyNA(simone_matrix)) {
      stop(
        "The SIMONE input contains missing values ",
        "after preprocessing."
      )
    }
    
    if (any(!is.finite(simone_matrix))) {
      stop(
        "The SIMONE input contains non-finite values ",
        "after preprocessing."
      )
    }
    
    nonconstant_flag <- vapply(
      df_simone,
      function(variable) {
        length(unique(variable)) > 1L
      },
      logical(1)
    )
    
    nonconstant_cols <- names(df_simone)[
      nonconstant_flag
    ]
    
    removed_constant_cols <- names(df_simone)[
      !nonconstant_flag
    ]
    
    df_simone_fit <- df_simone[
      ,
      nonconstant_cols,
      drop = FALSE
    ]
    
    if (ncol(df_simone_fit) < 2L) {
      stop(
        "SIMONE cannot be fitted because fewer than ",
        "two nonconstant variables remain."
      )
    }
    
    coef_df <- data.frame(
      matrix(
        0,
        nrow = length(y_cols),
        ncol = length(x_cols),
        dimnames = list(
          y_cols,
          x_cols
        )
      ),
      check.names = FALSE
    )
    
    # --------------------------------------------------------------------------
    # 3.4 Configure and fit SIMONE
    # --------------------------------------------------------------------------
    
    stage <- "set_options"
    
    simone_control <- setOptions(
      clusters.crit = "BIC",
      verbose = FALSE
    )
    
    stage <- "simone_fit"
    
    fit_simone <- simone::simone(
      X = df_simone_fit,
      type = "steady-state",
      clustering = clustering,
      control = simone_control
    )
    
    if (is.null(fit_simone)) {
      stop("simone() returned NULL.")
    }
    
    # --------------------------------------------------------------------------
    # 3.5 Extract the BIC-selected network
    # --------------------------------------------------------------------------
    
    stage <- "get_network"
    
    simone_network <- getNetwork(
      fit_simone,
      selection = "BIC"
    )
    
    if (is.null(simone_network)) {
      stop("getNetwork() returned NULL.")
    }
    
    if (is.null(simone_network$A)) {
      stop("The SIMONE network does not contain an adjacency matrix A.")
    }
    
    edge_matrix <- as.matrix(
      simone_network$A
    )
    
    storage.mode(edge_matrix) <- "double"
    
    expected_dimension <- ncol(df_simone_fit)
    
    if (nrow(edge_matrix) != expected_dimension ||
        ncol(edge_matrix) != expected_dimension) {
      stop(
        "The dimensions of the SIMONE adjacency matrix do not ",
        "match the number of variables supplied to SIMONE."
      )
    }
    
    if (anyNA(edge_matrix) ||
        any(!is.finite(edge_matrix))) {
      stop(
        "The SIMONE adjacency matrix contains missing ",
        "or non-finite values."
      )
    }
    
    if (is.null(rownames(edge_matrix))) {
      rownames(edge_matrix) <- colnames(df_simone_fit)
    }
    
    if (is.null(colnames(edge_matrix))) {
      colnames(edge_matrix) <- colnames(df_simone_fit)
    }
    
    matrix_is_symmetric <- isTRUE(
      all.equal(
        unname(edge_matrix),
        unname(t(edge_matrix)),
        tolerance = 1e-8,
        check.attributes = FALSE
      )
    )
    
    if (!matrix_is_symmetric) {
      stop(
        "The SIMONE adjacency matrix is not symmetric, ",
        "although an undirected network was expected."
      )
    }
    
    # --------------------------------------------------------------------------
    # 3.6 Locate the unique nonzero undirected edges
    # --------------------------------------------------------------------------
    
    stage <- "identify_network_edges"
    
    nonzero_mask <- is.finite(edge_matrix) &
      abs(edge_matrix) > 0
    
    # Exclude diagonal entries.
    diag(nonzero_mask) <- FALSE
    
    nonzero_entries <- which(
      nonzero_mask &
        upper.tri(nonzero_mask),
      arr.ind = TRUE
    )
    
    n_network_edges <- nrow(nonzero_entries)
    
    if (n_network_edges == 0L) {
      
      return(
        new_fit_result(
          status = "success_no_edges",
          coefficients = coef_df,
          diagnostics = list(
            method = "simone",
            p = p,
            clustering = clustering,
            selection = "BIC",
            n_observations = nrow(df_simone_fit),
            n_input_nodes = length(all_cols),
            n_nonconstant_nodes = length(nonconstant_cols),
            n_removed_constant_nodes =
              length(removed_constant_cols),
            removed_constant_nodes =
              removed_constant_cols,
            n_selected_network_edges = 0L,
            n_selected_edges = 0L
          )
        )
      )
    }
    
    # --------------------------------------------------------------------------
    # 3.7 Retain only lagged-to-current temporal edges
    # --------------------------------------------------------------------------
    
    stage <- "map_temporal_edges"
    
    for (i in seq_len(nrow(nonzero_entries))) {
      
      row_index <- nonzero_entries[i, "row"]
      column_index <- nonzero_entries[i, "col"]
      
      node_row <- map_simone_node(
        rownames(edge_matrix)[row_index],
        colnames(df_simone_fit)
      )
      
      node_column <- map_simone_node(
        colnames(edge_matrix)[column_index],
        colnames(df_simone_fit)
      )
      
      edge_value <- edge_matrix[
        row_index,
        column_index
      ]
      
      if (!is.finite(edge_value)) {
        next
      }
      
      # If row is lagged and column is current, orient row -> column.
      if (node_row %in% x_cols &&
          node_column %in% y_cols) {
        
        coef_df[
          node_column,
          node_row
        ] <- edge_value
      }
      
      # If column is lagged and row is current, orient column -> row.
      if (node_column %in% x_cols &&
          node_row %in% y_cols) {
        
        coef_df[
          node_row,
          node_column
        ] <- edge_value
      }
    }
    
    coefficient_matrix <- as.matrix(coef_df)
    
    n_temporal_edges <- sum(
      is.finite(coefficient_matrix) &
        abs(coefficient_matrix) > 0,
      na.rm = TRUE
    )
    
    fit_status <- if (n_temporal_edges == 0L) {
      "success_no_edges"
    } else {
      "success"
    }
    
    # --------------------------------------------------------------------------
    # 3.8 Return the completed fit result
    # --------------------------------------------------------------------------
    
    new_fit_result(
      status = fit_status,
      coefficients = coef_df,
      diagnostics = list(
        method = "simone",
        p = p,
        clustering = clustering,
        selection = "BIC",
        n_observations = nrow(df_simone_fit),
        n_input_nodes = length(all_cols),
        n_nonconstant_nodes = length(nonconstant_cols),
        n_removed_constant_nodes =
          length(removed_constant_cols),
        removed_constant_nodes =
          removed_constant_cols,
        n_selected_network_edges = n_network_edges,
        n_selected_edges = n_temporal_edges
      )
    )
    
  }, error = function(e) {
    
    failed_fit(
      error = e,
      method = "simone",
      p = p,
      stage = stage
    )
  })
}


# ------------------------------------------------------------------------------
# 4. Run SIMONE for all candidate lag specifications
# ------------------------------------------------------------------------------

simone_fit_results_by_lag <- list()
simone_results_by_lag <- list()
simone_graphs_by_lag <- list()

for (p in candidate_lags) {
  
  p_name <- paste0("p", p)
  
  message("=================================================")
  message("Running SIMONE for lag specification ", p_name)
  message("=================================================")
  
  if (!p_name %in% names(df_structure_by_lag)) {
    warning(
      "No structure-learning dataset found for ",
      p_name,
      ". Skipping."
    )
    next
  }
  
  df_structure_p <- df_structure_by_lag[[p_name]]
  
  simone_fit_p <- simone_apply_structure_data(
    df_structure_p = df_structure_p,
    p = p,
    clustering = TRUE
  )
  
  if (fit_succeeded(simone_fit_p)) {
    
    graph_attempt <- tryCatch(
      generate_simone_graph_from_lagged_coefs(
        simone_fit_p$coefficients,
        threshold = simone_results_threshold
      ),
      error = function(e) e
    )
    
    if (inherits(graph_attempt, "condition")) {
      
      simone_fit_p <- failed_fit(
        error = graph_attempt,
        method = "simone",
        p = p,
        stage = "graph_construction"
      )
      
    } else {
      
      simone_fit_p$graph <- graph_attempt
      
      graph_arcs <- bnlearn::arcs(
        simone_fit_p$graph
      )
      
      n_graph_edges <- if (is.null(graph_arcs)) {
        0L
      } else {
        nrow(graph_arcs)
      }
      
      simone_fit_p$diagnostics$n_graph_edges <-
        n_graph_edges
      
      simone_fit_p$status <- if (n_graph_edges == 0L) {
        "success_no_edges"
      } else {
        "success"
      }
    }
  }
  
  simone_fit_results_by_lag[[p_name]] <- simone_fit_p
  
  simone_results_by_lag[p_name] <- list(
    simone_fit_p$coefficients
  )
  
  simone_graphs_by_lag[p_name] <- list(
    simone_fit_p$graph
  )
  
  if (!is.null(simone_fit_p$coefficients)) {
    
    write.csv(
      simone_fit_p$coefficients,
      file = paste0(
        "simone_coefficients_",
        p_name,
        ".csv"
      ),
      row.names = TRUE
    )
  }
  
  if (!is.null(simone_fit_p$graph)) {
    
    bnlearn::write.dot(
      simone_fit_p$graph,
      file = paste0(
        "simone_",
        p_name,
        "_dot_file.dot"
      )
    )
  }
  
  if (identical(simone_fit_p$status, "failed")) {
    
    warning(
      "SIMONE failed for ",
      p_name,
      " during ",
      simone_fit_p$diagnostics$failed_stage,
      ": ",
      simone_fit_p$error
    )
  }
}

saveRDS(
  simone_fit_results_by_lag,
  file = "simone_fit_results_by_lag.rds"
)