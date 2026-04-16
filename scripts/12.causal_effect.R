library(causaleffect)


nodes <- c("New.cases_t", 
           "New.infections_t",
           "Reinfections_t")

# Initialize a list to store the results
results <- list()

# Iterate over the nodes
for (node in nodes) {
  # Find parents and ancestors
  parents <- parents(graph_js, node)
  ancestors <- ancestors(graph_js, node)
  
  # Store the results in the list
  results[[node]] <- list("Parents" = parents, "Ancestors" = ancestors)
}

print(results)

igraph_network <- convert_bn_to_igraph(graph_js)

cat("Node Names: ", toString(V(igraph_network)$name), "\n")

effect <- causal.effect(x = c("OpenTable.restaurant.bookings.London.index_t-1"), y = c("Reinfections_t"), simp = TRUE, prune = TRUE, G = igraph_network)
print(effect)

lar_bn_fitted <- bn.fit(graph_js, df_full_graph)

# Define variables
x_variables <- c("Citymapper.journeys.mobility.index_t-1", 
                 "OpenTable.restaurant.bookings.London.index_t-1")

y_variables <- c("Reinfections_t")

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
