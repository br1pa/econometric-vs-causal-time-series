library(igraph)

compute_bn_metrics <- function(bn_str, kb, Data, debug = FALSE) {
  # Ensure that bn_str is a bn object
  if (!inherits(bn_str, "bn")) {
    stop("The bn_str argument must be a bn object.")
  }
  
  # Ensure that kb is a bn object
  if (!inherits(kb, "bn")) {
    stop("The kb argument must be a bn object.")
  }
  
  # Ensure that Data is a data frame
  if (!is.data.frame(Data)) {
    stop("The Data argument must be a data frame.")
  }
  
  # Compute SHD (Structural Hamming Distance)
  shd_val <- shd(bn_str, kb, debug = debug)
  
  # Compute number of parameters
  nparams_val <- nparams(bn_str, data = Data, effective = TRUE)
  
  # Compute log-likelihood
  logLik_val <- logLik(bn_str, data = Data)
  
  # Compute BN score
  bn_score <- logLik(bn_str, data = Data) - ((nparams(bn_str, data = Data, effective = TRUE) * log(nrow(Data)))/2)
  
  # Compute number of edges
  arc_set <- arcs(bn_str)
  number_of_edges <- nrow(arc_set)
  
  # Compute number of independent graphical fragments
  graph <- graph_from_data_frame(arc_set, directed = TRUE)
  number_of_components <- count_components(graph)
  
  # Create a list to store the results
  results <- list(SHD = shd_val,
                  BN_Score = bn_score,
                  n_params = nparams_val,
                  LogLikelihood = as.numeric(logLik_val),
                  n_edges = number_of_edges,
                  n_components = number_of_components)
  
  return(results)
}

edge_selector <- function(edge_matrix, n_occurences) {
  
  arc_table <- as.data.frame(edge_matrix)
  arc_counts <- arc_table %>%
    count(from, to) %>%
    rename(Freq = n) %>%
    arrange(desc(Freq)) %>%
    filter(Freq >= n_occurences)
  
  return(arc_counts)
}

extract_edges <- function(bn) {
  edges <- arcs(bn)
  return(edges)
}

all_edges <- do.call(rbind, lapply(list(graph_lars, graph_lar, graph_js, graph_simone), extract_edges))
arcs_selected <- edge_selector(all_edges, n_occurences = 2) #2 is the optimal choice

arcs_to_add <- as.matrix(arcs_selected[,1:2])
average_bn <- empty.graph(nodes(empty_net))
arcs(average_bn) <- arcs_to_add

compute_bn_metrics(average_bn, kb, df_full_graph)
