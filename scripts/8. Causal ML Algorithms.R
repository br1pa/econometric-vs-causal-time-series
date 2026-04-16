library(parallel)
library(bnlearn)
library(igraph)

convert_bn_to_igraph <- function(bn) {
  # Get arc list
  arc_mat <- bnlearn::arcs(bn)  # matrix with columns "from", "to"
  
  # If there are no edges, build an empty graph with the right node names
  if (is.null(arc_mat) || nrow(arc_mat) == 0) {
    node_names <- bnlearn::nodes(bn)
    g <- igraph::make_empty_graph(n = length(node_names), directed = TRUE)
    igraph::V(g)$name <- node_names
    return(g)
  }
  
  # Convert arcs matrix to a data.frame for igraph
  arc_df <- as.data.frame(arc_mat, stringsAsFactors = FALSE)
  colnames(arc_df) <- c("from", "to")
  
  # Build directed igraph
  g <- igraph::graph_from_data_frame(arc_df, directed = TRUE, vertices = unique(c(arc_df$from, arc_df$to)))
  
  return(g)
}

learn_bn_structures_continuous <- function(dataframe) {

  # Algorithms compatible with continuous data
  algorithms <- c("pc.stable",
                  "gs",  
                  "iamb", 
                  "fast.iamb", 
                  "inter.iamb",
                  "iamb.fdr",
                  "h2pc",
                  "hc", 
                  "tabu", 
                  "mmhc", 
                  "rsmax2")
  
  # Initialize a list to store the learned structures
  learned_structures <- list()
  
  # Iterate over the algorithms and learn the structures
  for (alg in algorithms) {
    cat("Learning structure using", alg, "\n")
    tryCatch({
      # Call the correct bnlearn function based on the algorithm name
      bn_structure <- switch(alg,
                             "pc.stable" = pc.stable(dataframe, alpha = 0.05, test = 'mi-g-sh', cluster = makeCluster(detectCores() - 1)),
                             "gs" = gs(dataframe, alpha = 0.05, test = 'mi-g-sh', cluster = makeCluster(detectCores() - 1)),
                             "iamb" = iamb(dataframe, alpha = 0.05, test = 'mi-g-sh', cluster = makeCluster(detectCores() - 1)),
                             "fast.iamb" = fast.iamb(dataframe, alpha = 0.05, test = 'mi-g-sh', cluster = makeCluster(detectCores() - 1)),
                             "inter.iamb" = inter.iamb(dataframe, alpha = 0.05, test = 'mi-g-sh', cluster = makeCluster(detectCores() - 1)),
                             "iamb.fdr" = iamb.fdr(dataframe, alpha = 0.05, test = 'mi-g-sh', cluster = makeCluster(detectCores() - 1)),
                             "h2pc" = h2pc(dataframe, restrict.args = list (alpha = 0.05, test = 'mi-g-sh'), maximize.args = list(score = "ebic-g")),
                             "hc" =  hc(dataframe, score = "ebic-g"),
                             "tabu" = tabu(dataframe, score = "ebic-g"),
                             "mmhc" = mmhc(dataframe, restrict.args = list (alpha = 0.05, test = 'mi-g-sh'), maximize.args = list(score = "ebic-g")),
                             "rsmax2" = rsmax2(dataframe, restrict.args = list (alpha = 0.05, test = 'mi-g-sh'), maximize.args = list(score = "ebic-g")),
                             stop(paste("Unknown algorithm:", alg))
      )
      # Store the learned structure in the list
      learned_structures[[alg]] <- bn_structure
    }, error = function(e) {
      cat("Error in", alg, ":", e$message, "\n")
      learned_structures[[alg]] <- NULL
    })
  }
  
  # Return the list of learned structures
  return(learned_structures)
}

min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

df_learnt_bns <- df_full_graph %>% dplyr::mutate(across(where(is.numeric), min_max_norm))

colnames(df_learnt_bns) <- colnames(df_full_graph)

learnt_bns <- learn_bn_structures_continuous(df_learnt_bns)

# GO BACK TO 7. PARAMETRISATION_PROCESSING.R AND RUN THE REST -------------------------------------------------------------

has_directed_edges <- function(bn_structure) {
  # Extract the arc set from the bn structure
  arcs <- arcs(bn_structure)
  
  # Check if there are any arcs defined
  if (nrow(arcs) > 0) {
    return(TRUE)  # The structure has directed edges
  } else {
    return(FALSE)  # The structure is a skeleton without directed edges
  }
}

# Filter learnt_bns to include only structures with directed edges
filtered_learnt_bns <- Filter(has_directed_edges, learnt_bns)

# Initialize a list for errors
errors <- list()

# Apply cextend with error handling
converted_dags <- lapply(seq_along(filtered_learnt_bns), function(i) {
  tryCatch(
    {
      # Apply cextend
      cextend(filtered_learnt_bns[[i]])
    },
    error = function(e) {
      # Log the error with index
      errors[[i]] <- e$message
      return(NULL) # Return NULL for failed structures
    }
  )
})

# Print errors
if (length(errors) > 0) {
  message("Structures with errors: ", paste(names(errors), collapse = ", "))
  print(errors)
}

# Assign names to the list if missing
algorithm_names <- c("pc.stable", "gs", "iamb", "fast.iamb", "inter.iamb", 
                     "iamb.fdr", "h2pc", "hc", "tabu", "mmhc", "rsmax2")
names(converted_dags) <- algorithm_names

# Iterate through the algorithms in converted_dags
for (alg in names(converted_dags)) {
  print(paste("Processing:", alg))
  
  # Extract the Bayesian network structure
  bn_structure <- converted_dags[[alg]]
  
  if (is.null(bn_structure)) {
    print(paste("Skipping", alg, "because it is NULL"))
    next
  }
  
  # Initialize variables to store the results
  number_of_edges <- NA
  shd_val <- NA
  nparams_val <- NA
  logLik_val <- NA
  bn_score <- NA
  
  tryCatch(
    {
      # Compute arcs and number of edges
      arc_set <- arcs(bn_structure)
      number_of_edges <- nrow(arc_set)
      
      # Compute SHD
      shd_val <- shd(bn_structure, kb)
      
      # Compute number of parameters
      nparams_val <- nparams(bn_structure, data = df_full_graph, effective = TRUE)
      
      # Compute log-likelihood
      logLik_val <- logLik(bn_structure, data = df_full_graph)
      
      # Compute Bayesian network score
      bn_score <- logLik_val - ((nparams_val * log(nrow(df_full_graph))) / 2)
      
      # Print results for this algorithm
      print(paste("Number of edges for", alg, ":", number_of_edges))
      print(paste("SHD for", alg, ":", shd_val))
      print(paste("Number of parameters for", alg, ":", nparams_val))
      print(paste("Log-likelihood for", alg, ":", logLik_val))
      print(paste("Bayesian network score for", alg, ":", bn_score))
    },
    error = function(e) {
      print(paste("Error processing", alg, ":", e$message))
    }
  )
}

# --------------------------------------------------------------------------------

# List of nodes for which you want to find parents and ancestors
nodes <- c("New.cases_t", 
           "New.infections_t",
           "Reinfections_t")

# Initialize a list to store the results
results <- list()

# Iterate over the nodes
for (node in nodes) {
  # Find parents and ancestors
  parents <- parents(converted_dags$tabu, node)
  ancestors <- ancestors(converted_dags$tabu, node)
  
  # Store the results in the list
  results[[node]] <- list("Parents" = parents, "Ancestors" = ancestors)
}

print(results)

igraph_network <- convert_bn_to_igraph(converted_dags$tabu)

cat("Node Names: ", toString(V(igraph_network)$name), "\n")

effect <- causal.effect(x = c("Flights.7.day.moving.average_t-1"), y = c("New.cases_t"), simp = TRUE, prune = TRUE, G = igraph_network)
print(effect)

lar_bn_fitted <- bn.fit(converted_dags$tabu, df_full_graph)

# Define your variables
x_variables <- c("Citymapper.journeys.mobility.index_t-1",
                 "TfL.Bus.mobility.index_t-1",
                 "Flights.7.day.moving.average_t-1",
                 "TfL.Tube.mobility.index_t-1",
                 "Google.grocery.pharmacy.Greater.London.mobility.index_t-1",
                 "Apple.walking.London.mobility.index_t-1",
                 "Google.transit.stations.mobility.index_t-1",
                 "Google.retail.recreation.Greater.London.mobility.index_t-1",
                 "OpenTable.restaurant.bookings.London.index_t-1")

y_variables <- c("New.cases_t", "New.infections_t", "Reinfections_t")

# Initialize a list to store ACE values for each (X, Y) pair
results <- list()

# Loop through each combination of X and Y variables
for (y_var in y_variables) {
  for (x_var in x_variables) {
    
    # Compute E_Y_do_X for different values of x_var
    E_Y_do_X_2 <- cpquery(
      lar_bn_fitted,
      event = (get(y_var) == "3"),  # Dynamically get the Y variable
      evidence = (get(x_var) == "2"),  # Dynamically get the X variable
      n = 100000
    )
    
    E_Y_do_X_3 <- cpquery(
      lar_bn_fitted,
      event = (get(y_var) == "3"),  # Dynamically get the Y variable
      evidence = (get(x_var) == "3"),  # Dynamically get the X variable
      n = 100000
    )
    
    E_Y_do_X_1 <- cpquery(
      lar_bn_fitted,
      event = (get(y_var) == "3"),  # Dynamically get the Y variable
      evidence = (get(x_var) == "1"),  # Dynamically get the X variable
      n = 100000
    )
    
    # Calculate ACE for the pair (X, Y)
    ACE_ij <- 0.5 * (E_Y_do_X_2 - E_Y_do_X_3 + E_Y_do_X_1 - E_Y_do_X_3)
    
    # Print or store the result
    print(paste("ACE for", y_var, "and", x_var, "is", ACE_ij))
    results[[paste(y_var, x_var, sep = "_")]] <- ACE_ij
  }
}

# View results
results
