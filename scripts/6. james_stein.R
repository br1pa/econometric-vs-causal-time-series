# 6. james_stein.R --------------------------------------------------------------

library(GeneNet)
library(bnlearn)
library(dplyr)

set.seed(123)

min_max_norm <- function(x) {
  lo <- min(x)
  hi <- max(x)
  
  if (!is.finite(lo) || !is.finite(hi)) {
    stop("Non-finite values in input.")
  }
  
  if (hi == lo) {
    return(rep(0, length(x)))
  }
  
  (x - lo) / (hi - lo)
}

# ------------------------------------------------------------------------------
# 0. Checks
# ------------------------------------------------------------------------------

if (!exists("df_structure_by_lag")) {
  stop("df_structure_by_lag does not exist. 
       Run 4.parametrisation_processing.R before 6.james_stein.R.")
}

if (!exists("selected_p")) {
  selected_p <- 7
}

if (!exists("new_fit_result", mode = "function") ||
    !exists("failed_fit", mode = "function") ||
    !exists("fit_succeeded", mode = "function")) {
  source("4a. Helpers.R")
}

js_ggm_method <- "prob"
js_ggm_cutoff <- 0.95
js_results_threshold <- 0

# ------------------------------------------------------------------------------
# 1. Convert GeneNet output into a lagged coefficient matrix
# ------------------------------------------------------------------------------

james_stein_apply_structure_data <- function(
    df_structure_p,
    p,
    ggm_cutoff = js_ggm_cutoff
) {
  
  stage <- "input_validation"
  
  tryCatch({
    
    # --------------------------------------------------------------------------
    # 1. Validate inputs
    # --------------------------------------------------------------------------
    
    if (!is.data.frame(df_structure_p)) {
      stop("df_structure_p must be a data frame.")
    }
    
    if (nrow(df_structure_p) == 0 || ncol(df_structure_p) == 0) {
      stop("df_structure_p is empty.")
    }
    
    if (anyDuplicated(names(df_structure_p))) {
      stop("df_structure_p contains duplicated column names.")
    }
    
    if (!is.numeric(p) || length(p) != 1L ||
        !is.finite(p) || p < 1) {
      stop("p must be one positive integer.")
    }
    
    p <- as.integer(p)
    
    if (!is.numeric(ggm_cutoff) ||
        length(ggm_cutoff) != 1L ||
        !is.finite(ggm_cutoff) ||
        ggm_cutoff < 0 ||
        ggm_cutoff > 1) {
      stop("ggm_cutoff must be a number between zero and one.")
    }
    
    if (!exists("min_max_norm", mode = "function")) {
      stop(
        "min_max_norm() is unavailable. ",
        "Run 5. lars.R first."
      )
    }
    
    # --------------------------------------------------------------------------
    # 2. Identify lagged predictors and current outcomes
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
        "Some columns were classified as both lagged ",
        "predictors and current outcomes."
      )
    }
    
    # --------------------------------------------------------------------------
    # 3. Prepare and scale the GeneNet input
    # --------------------------------------------------------------------------
    
    stage <- "prepare_data"
    
    selected_cols <- c(x_cols, y_cols)
    
    df_js <- df_structure_p[
      ,
      selected_cols,
      drop = FALSE
    ]
    
    df_js[] <- lapply(
      df_js,
      min_max_norm
    )
    
    x_mat <- data.matrix(df_js)
    storage.mode(x_mat) <- "double"
    
    if (anyNA(x_mat)) {
      stop(
        "The James-Stein input contains missing values ",
        "after preprocessing."
      )
    }
    
    if (any(!is.finite(x_mat))) {
      stop(
        "The James-Stein input contains non-finite values ",
        "after preprocessing."
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
    # 4. Estimate the shrinkage partial-correlation matrix
    # --------------------------------------------------------------------------
    
    stage <- "ggm.estimate.pcor"
    
    dyn <- GeneNet::ggm.estimate.pcor(
      x_mat,
      method = "static"
    )
    
    if (is.null(dyn)) {
      stop("ggm.estimate.pcor() returned NULL.")
    }
    
    # --------------------------------------------------------------------------
    # 5. Test the candidate network edges
    # --------------------------------------------------------------------------
    
    stage <- "network.test.edges"
    
    tested_edges <- GeneNet::network.test.edges(
      dyn,
      plot = FALSE
    )
    
    if (is.null(tested_edges)) {
      stop("network.test.edges() returned NULL.")
    }
    
    # --------------------------------------------------------------------------
    # 6. Select edges using the specified posterior-probability cutoff
    # --------------------------------------------------------------------------
    
    stage <- "extract.network"
    
    js_net <- GeneNet::extract.network(
      tested_edges,
      method.ggm = js_ggm_method,
      cutoff.ggm = ggm_cutoff
    )
    
    if (is.null(js_net)) {
      
      return(
        new_fit_result(
          status = "success_no_edges",
          coefficients = coef_df,
          diagnostics = list(
            method = "james_stein",
            p = p,
            n_observations = nrow(x_mat),
            n_input_nodes = ncol(x_mat),
            n_lagged_predictors = length(x_cols),
            n_current_targets = length(y_cols),
            n_selected_network_edges = 0L,
            n_selected_edges = 0L,
            ggm_method = js_ggm_method,
            ggm_cutoff = ggm_cutoff
          )
        )
      )
    }
    
    js_net <- as.data.frame(
      js_net,
      stringsAsFactors = FALSE
    )
    
    if (nrow(js_net) == 0) {
      
      return(
        new_fit_result(
          status = "success_no_edges",
          coefficients = coef_df,
          diagnostics = list(
            method = "james_stein",
            p = p,
            n_observations = nrow(x_mat),
            n_input_nodes = ncol(x_mat),
            n_lagged_predictors = length(x_cols),
            n_current_targets = length(y_cols),
            n_selected_network_edges = 0L,
            n_selected_edges = 0L,
            ggm_method = js_ggm_method,
            ggm_cutoff = ggm_cutoff
          )
        )
      )
    }
    
    # --------------------------------------------------------------------------
    # 7. Validate the GeneNet output
    # --------------------------------------------------------------------------
    
    stage <- "validate_genenet_output"
    
    required_js_columns <- c(
      "node1",
      "node2",
      "pcor"
    )
    
    missing_js_columns <- setdiff(
      required_js_columns,
      names(js_net)
    )
    
    if (length(missing_js_columns) > 0) {
      stop(
        "GeneNet output is missing required column(s): ",
        paste(
          missing_js_columns,
          collapse = ", "
        ),
        ". Edge weights must come from 'pcor'."
      )
    }
    
    map_genenet_node <- function(node,
                                 available_columns) {
      
      node <- as.character(node)
      
      node_index <- suppressWarnings(
        as.integer(node)
      )
      
      if (!is.na(node_index) &&
          node_index >= 1L &&
          node_index <= length(available_columns)) {
        
        return(
          available_columns[node_index]
        )
      }
      
      return(node)
    }
    
    # --------------------------------------------------------------------------
    # 8. Retain only lagged-to-current temporal edges
    # --------------------------------------------------------------------------
    
    stage <- "map_temporal_edges"
    
    for (i in seq_len(nrow(js_net))) {
      
      node1 <- map_genenet_node(
        js_net[["node1"]][i],
        colnames(df_js)
      )
      
      node2 <- map_genenet_node(
        js_net[["node2"]][i],
        colnames(df_js)
      )
      
      pcor_value <- suppressWarnings(
        as.numeric(
          as.character(
            js_net[["pcor"]][i]
          )
        )
      )
      
      if (!is.finite(pcor_value)) {
        next
      }
      
      # GeneNet returns undirected partial-correlation edges.
      # Retain and orient only lagged-current adjacencies.
      
      if (node1 %in% x_cols &&
          node2 %in% y_cols) {
        
        coef_df[node2, node1] <- pcor_value
      }
      
      if (node2 %in% x_cols &&
          node1 %in% y_cols) {
        
        coef_df[node1, node2] <- pcor_value
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
    # 9. Return the successful result
    # --------------------------------------------------------------------------
    
    new_fit_result(
      status = fit_status,
      coefficients = coef_df,
      diagnostics = list(
        method = "james_stein",
        p = p,
        n_observations = nrow(x_mat),
        n_input_nodes = ncol(x_mat),
        n_lagged_predictors = length(x_cols),
        n_current_targets = length(y_cols),
        n_selected_network_edges = nrow(js_net),
        n_selected_edges = n_temporal_edges,
        ggm_method = js_ggm_method,
        ggm_cutoff = ggm_cutoff
      )
    )
    
  }, error = function(e) {
    
    failed_fit(
      error = e,
      method = "james_stein",
      p = p,
      stage = stage
    )
  })
}


# ------------------------------------------------------------------------------
# 2. Convert coefficient matrix into a temporal DAG
# ------------------------------------------------------------------------------

generate_graph_from_lagged_coefs <- function(coef_df, threshold = 0) {
  
  predictor_nodes <- colnames(coef_df)
  target_nodes <- rownames(coef_df)
  all_nodes <- unique(c(predictor_nodes, target_nodes))
  
  g <- empty.graph(all_nodes)
  
  for (target in target_nodes) {
    for (pred in predictor_nodes) {
      
      val <- coef_df[target, pred]
      
      if (is.finite(val) && abs(val) > threshold) {
        g <- set.arc(
          g,
          from = pred,
          to = target,
          check.cycles = TRUE
        )
      }
    }
  }
  
  return(g)
}


# ------------------------------------------------------------------------------
# 3. Run James-Stein shrinkage for all candidate lag specifications
# ------------------------------------------------------------------------------

js_fit_results_by_lag <- list()
js_results_by_lag <- list()
js_graphs_by_lag <- list()

for (p in candidate_lags) {
  
  p_name <- paste0("p", p)
  
  message("=================================================")
  message("Running James-Stein for lag specification ", p_name)
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
  
  js_fit_p <- james_stein_apply_structure_data(
    df_structure_p = df_structure_p,
    p = p,
    ggm_cutoff = js_ggm_cutoff
  )
  
  if (fit_succeeded(js_fit_p)) {
    
    graph_attempt <- tryCatch(
      generate_graph_from_lagged_coefs(
        js_fit_p$coefficients,
        threshold = js_results_threshold
      ),
      error = function(e) e
    )
    
    if (inherits(graph_attempt, "condition")) {
      
      js_fit_p <- failed_fit(
        error = graph_attempt,
        method = "james_stein",
        p = p,
        stage = "graph_construction"
      )
      
    } else {
      
      js_fit_p$graph <- graph_attempt
      
      graph_arcs <- bnlearn::arcs(
        js_fit_p$graph
      )
      
      n_graph_edges <- if (is.null(graph_arcs)) {
        0L
      } else {
        nrow(graph_arcs)
      }
      
      js_fit_p$diagnostics$n_graph_edges <- n_graph_edges
      
      js_fit_p$status <- if (n_graph_edges == 0L) {
        "success_no_edges"
      } else {
        "success"
      }
    }
  }
  
  js_fit_results_by_lag[[p_name]] <- js_fit_p
  
  js_results_by_lag[p_name] <- list(
    js_fit_p$coefficients
  )
  
  js_graphs_by_lag[p_name] <- list(
    js_fit_p$graph
  )
  
  if (!is.null(js_fit_p$coefficients)) {
    
    write.csv(
      js_fit_p$coefficients,
      file = paste0(
        "js_coefficients_",
        p_name,
        ".csv"
      ),
      row.names = TRUE
    )
  }
  
  if (!is.null(js_fit_p$graph)) {
    
    bnlearn::write.dot(
      js_fit_p$graph,
      file = paste0(
        "js_",
        p_name,
        "_dot_file.dot"
      )
    )
  }
}