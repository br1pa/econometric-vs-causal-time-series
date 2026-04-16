discretize_kmeans <- function(data, num_breaks) {
  cluster <- kmeans(data, centers = num_breaks)
  binned_data <- cut(data, breaks = c(-Inf, sort(cluster$centers), Inf), include.lowest = TRUE, labels = FALSE)
  return(binned_data)
}


for (col_name in names(df)) {
  if (col_name %in% names(factor_levels)) {
    df[[col_name]] <- factor(df[[col_name]], labels = factor_levels[[col_name]])
  }
}

add_suffix_and_duplicate <- function(df) {
  # Adding the t-1 suffix to original columns
  names(df) <- paste0(names(df), "_t-1")
  
  # Duplicating the columns with a t suffix
  df_t = df %>% 
    rename_with(~paste0(sub("_t-1$", "", .), "_t"), everything())
  
  # Combining both sets of columns
  return(bind_cols(df, df_t))
}

df_full_graph <- add_suffix_and_duplicate(df)

df_full_graph <- df_full_graph %>% dplyr::select(-Date_t, -`Date_t-1`)
# run 8. Causal ML Algorithms.R -----------------------------------------------------------------------
numeric_columns <- sapply(df_full_graph, is.numeric)
df_full_graph[numeric_columns] <- lapply(df_full_graph[numeric_columns], discretize_kmeans, 2)
any(is.na(df_full_graph))
df_full_graph[numeric_columns] <- lapply(df_full_graph[numeric_columns], factor)

df_full_graph <- as.data.frame(df_full_graph)
nodes <- colnames(df_full_graph)          # Get the column names of your data as nodes
empty_net <- empty.graph(nodes)  # Create an empty graph with these nodes

for (i in nodes) {
  for (j in nodes) {
    if (i != j) {
      empty_net <- set.arc(empty_net, from = i, to = j)
    }
  }
}
