library(simone)

simone_apply <- function(df){
  
  # Parameters
  # ----------
  # df : data.frame
  #   Dataframe to apply SIMONE models to
  # 
  # Returns
  # -------
  # total_edge_matrix : matrix
  #   A matrix where the rows and columns represent nodes, and the values represent how many times the edge between those two nodes occurs in 
  #   the submodels
  
  ctrl = setOptions(clusters.crit = 'BIC', verbose = FALSE)
  
  total_edge_matrix <- matrix(0, nrow = (ncol(df) - 1), ncol = (ncol(df) - 1))
  
  ordered_col_names_df <- subset(df, select=-c(Date))
  colnames(total_edge_matrix) <- colnames(ordered_col_names_df)
  rownames(total_edge_matrix) <- colnames(ordered_col_names_df)
  
  
  fit_simone <- simone(subset(df, select=-c(Date)), type ='time-course', clustering = TRUE, control = ctrl)
  var <- getNetwork(fit_simone, selection = 'BIC')
  edge_matrix <- var$A
  
  total_edge_matrix[rownames(edge_matrix), colnames(edge_matrix)] <- edge_matrix
  
  total_edge_matrix[,]<- edge_matrix
  
  
  return(total_edge_matrix)
}

generate_graph_simone <- function(simone_results, simone_results_threshold = 1){
  
  nodest1 <- paste0(colnames(simone_results), "_t-1")
  nodest2 <- paste0(colnames(simone_results), "_t")
  nodes <- c(nodest1,nodest2)
  
  dot_file <- empty.graph(nodes)
  
  for (j in 1:ncol(simone_results)) {
    yt <- paste0(colnames(simone_results)[[j]], "_t")
    for (i in 1:nrow(simone_results)){
      if (simone_results[i,j] != 0 && is.na(simone_results[i,j]) == FALSE && abs(simone_results[i,j]) > simone_results_threshold) {
        xtminus1 <- paste0(colnames(simone_results)[[i]], "_t-1")
        dot_file <- set.arc(dot_file, from = xtminus1, to = yt, check.cycles = TRUE)
        
        
      }
    }
  }
  return(dot_file)
}

simone_results <- simone_apply(df)

graph_simone <- generate_graph_simone(simone_results, simone_results_threshold = 0)

print(graph_simone)

write.dot(graph_simone, file = "simone_dot_file.dot")