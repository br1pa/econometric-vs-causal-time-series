library(GeneNet)

james_stein <- function(df, ggm_cutoff){
  
  df2 <- df
  df2 <- df2 %>% dplyr::select(-c(Date))
  n <- ncol(df2)
  coef_df <- data.frame(matrix(ncol = n, nrow = n))
  names(coef_df) <- names(df2)
  rownames(coef_df) <- names(df2)
  
  
  x <- lars_transformations(df, type ='x')
  x <- data.matrix(x)
  
  dyn <- ggm.estimate.pcor(x, method = "dynamic")
  arth.arcs <- network.test.edges(dyn, plot = FALSE)
  print(arth.arcs)
  arth.net <- extract.network(arth.arcs, method.ggm = 'prob', cutoff.ggm = ggm_cutoff)
  
  for (i in 1:nrow(arth.net)){
    col <- arth.net[i, 'node1']
    row <- arth.net[i, 'node2']
    
    coef_df[row, col] <- arth.net[i, 'pcor']
  }
  
  return(coef_df)
}

js_results <- james_stein(df, ggm_cutoff = 0.05)

graph_js <- generate_graph_lars(js_results, lars_results_threshold =  0)

print(graph_js)

write.dot(graph_js, file = "js_dot_file.dot")