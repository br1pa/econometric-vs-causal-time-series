library(lars)
library(bnlearn)
library(dplyr)

lars_transformations <- function(df, type){
  
  # This function transforms a data frame into a format required to run graphical lasso. It orders the data frame by date, 
  # after which it will remove the first or last observation based on if the data frame will be the x or y input for graphical lasso.
  
  #  Parameters
  #  ----------
  #  df : data.frame
  #    The dataframe which will be used for graphical lasso
  #  
  #  type : string
  #    "x" or "y" based on whether we are treating the dataframe as the x or y arguements for graphical lasso
  
  #  Returns
  #  -------
      
  #  df : dataframe
  #    Returns the same dataframe ordered by date with rows dropped.
  
  df <- df %>% 
    arrange(Date) 
  
  if (type == "x"){
    
    df <- df %>% 
      slice_tail(n=-1)%>%                          
      tidyr::unite(rowname, Date) %>%
      tibble::column_to_rownames() 
  }
  
  
  if (type == "y"){
    
    df <- df %>% 
      slice_head(n=-1)%>%                          
      tidyr::unite(rowname, Date) %>%
      tibble::column_to_rownames() 
  }
  
  return(df)
}

lars_apply <- function(df){
  
  #  This function enables the application of graphical lasso on each column of a dataframe, and binds these
  #  results in a dataframe.

  #  Parameters
  #  ----------
  
  #  df : data.frame
  #    The dataframe which will be used for graphical lasso  
  
  #  Returns
  #  -------
      
  #  coef_df : dataframe
  #    A dataframe with the results of graphical lasso running on each column of x appended to it.
  
  
  n <- ncol(df)
  coef_df <- data.frame(matrix(ncol = n, nrow = 0))
  names(coef_df) <- names(df)
  coef_df <- coef_df %>% dplyr::select(-c('Date'))
  
  
  # Normalizing results 
  min_max_norm <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }
  
  df_index <- df %>% dplyr::select(Date)
  df_vals <- df %>% dplyr::select(-c(Date))
  
  vals_names <- names(df_vals)
  
  df_vals<- as.data.frame(lapply(df_vals, min_max_norm))
  names(df_vals) <- vals_names
  
  df <- cbind(df_index, df_vals)
  
  
  for (column_name in names(df)) {
    
    # Skip date columns in coef_df
    
    if((column_name=='Date') ) next
    
    # If that column has been dropped from the previous steps
    
    if(!(column_name %in% names(df))){
      coeff <- rep(NA, n)
      names(coeff) <- names(df_final)
      coef_df <- bind_rows(coef_df, coeff)
      
    } else{
      
      x <- lars_transformations(df, type = 'x')
      x <- data.matrix(x)
      
      y <- lars_transformations(df, type = 'y')
      y <- y %>% dplyr::select(column_name)
      y <- data.matrix(y)
      
      lasso_cv <- cv.lars(y = y , x = x, mode = "fraction", plot.it = FALSE)
      frac <- lasso_cv$index[which.min(lasso_cv$cv)]
      lars_pred <- predict(lars(y = y, x = x, type = 'lasso'), s = frac, type = 'coef', mode = 'fraction')
      
      coeffs <- lars_pred$coefficients
      coeffs <- as.data.frame(t(coeffs))
      
      coef_df <- bind_rows(coef_df, coeffs)
    }
    
    
  }
  
  return(coef_df)
  
}


generate_graph_lars <- function(lars_results, lars_results_threshold = 0){
  
  nodest1 <- paste0(colnames(lars_results), "_t-1")
  nodest2 <- paste0(colnames(lars_results), "_t")
  nodes <- c(nodest1,nodest2)
  
  dot_file <- empty.graph(nodes)
  
  for (i in 1:nrow(lars_results)) {
    yt <- paste0(colnames(lars_results)[[i]], "_t")
    for (j in 1:ncol(lars_results)){
      if (lars_results[i,j] != 0 && is.na(lars_results[i,j]) == FALSE && abs(lars_results[i,j]) > lars_results_threshold) {
        xtminus1 <- paste0(colnames(lars_results)[[j]], "_t-1")
        dot_file <- set.arc(dot_file, from = xtminus1, to = yt, check.cycles = TRUE)
        
        
      }
    }
  }
  return(dot_file)
}

lar <- function(df) {
  
  n <- ncol(df)
  coef_df <- data.frame(matrix(ncol = n, nrow = 0))
  names(coef_df) <- names(df)
  coef_df <- coef_df %>% dplyr::select(-c('Date'))
  
  min_max_norm <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }
  
  df_index <- df %>% dplyr::select(c(Date))
  df_vals <- df %>% dplyr::select(-c(Date))
  vals_names <- names(df_vals)
  
  df_vals <- as.data.frame(lapply(df_vals, min_max_norm))
  names(df_vals) <- vals_names
  df <- cbind(df_index, df_vals)
  
  for (column_name in names(df)) {
    
    if (column_name == 'Date') next
    print(paste("Processing column:", column_name))
    
    x <- lars_transformations(df, type = 'x')
    x <- data.matrix(x)
    
    y <- lars_transformations(df, type = 'y')
    y <- y %>% dplyr::select(column_name)
    y <- data.matrix(y)
    
    lasso_cv <- cv.lars(y = y , x = x, type = 'lar', plot.it = FALSE)
    
    lars_pred <- predict(lars(y = y, x = x, type = 'lar'), type = 'coef')
    best_coef_index <- which.min(lasso_cv$cv)
    
    raw   <- lars_pred$coefficients[best_coef_index, , drop = TRUE]    
    coeffs <- as.data.frame(t(raw))                                    
    names(coeffs) <- vals_names      

    for (col in names(coef_df)) {
      if (!col %in% names(coeffs)) {
        coeffs[[col]] <- NA
      }
    }
    
    if (ncol(coeffs) == ncol(coef_df)) {
      coef_df <- rbind(coef_df, coeffs)
    } else {
      warning(paste("Mismatch in columns for", column_name, "Expected:", ncol(coef_df), "Got:", ncol(coeffs)))
    }
    
    print(coeffs)
  }
  
  return(coef_df)
}


# Initialize a list to store factor levels
factor_levels <- list()

# Loop through columns and convert factors to numeric
for (col_name in names(df)) {
  if (is.factor(df[[col_name]])) {
    # Store levels
    factor_levels[[col_name]] <- levels(df[[col_name]])
    
    # Convert factor to numeric
    df[[col_name]] <- as.numeric(df[[col_name]])
  }
}

df <- cbind(Date = date_column, df)


df <- df %>%
  mutate(across(where(is.factor), as.integer))

lars_results <- lars_apply(df)

lars_results_threshold = 0.4

graph_lars <- generate_graph_lars(lars_results, lars_results_threshold =  lars_results_threshold)

print(graph_lars)

write.dot(graph_lars, file = "lars_dot_file.dot")

lar_results <- lar(df)

lar_results_threshold = 0.4

graph_lar <- generate_graph_lars(lar_results, lars_results_threshold =  lar_results_threshold)

print(graph_lar)

write.dot(graph_lar, file = "lar_dot_file.dot")