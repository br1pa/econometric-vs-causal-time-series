library(imputeTS)

continuous_vars <- c()

for (var_name in names(df)) {
  if (is.numeric(df[[var_name]])) {
    continuous_vars <- c(continuous_vars, var_name)
  }
}

# Impute missing values in continuous columns using state space models
for (var_name in continuous_vars) {
  df[[var_name]] <- na_kalman(df[[var_name]], model = "auto.arima")
}


# Extracting categorical columns
categorical_vars <- c()

for (var_name in names(df)) {
  if (is.factor(df[[var_name]])) {
    categorical_vars <- c(categorical_vars, var_name)
  }
}

# Convert categorical columns to dummy variables
formula_str <- paste(" ~ ", paste(categorical_vars, collapse = " + "), " - 1")
df_dummy <- as.data.frame(model.matrix(as.formula(formula_str), data = df))

# Impute missing values in the dummy variables
df_dummy[] <- lapply(df_dummy, function(column) {
  return(na_kalman(column, model = "auto.arima"))
})

# Convert back to categorical from dummy
for (var_name in categorical_vars) {
  # Extract columns related to current variable
  related_cols <- grep(paste0("^", var_name, "\\."), names(df_dummy), value = TRUE)
  
  if (length(related_cols) == 0) {
    next # if there are no related columns, skip this iteration
  }
  
  # Find the max dummy column for each row
  max_col <- apply(df_dummy[related_cols], 1, which.max)
  
  # Check if any of the rows don't have a max value
  if (any(is.na(max_col))) {
    cat("No max column for some rows in variable:", var_name, "\n")
    next
  }
  
  # Convert back to factor
  df[[var_name]] <- factor(sapply(max_col, function(idx) {
    sub(paste0("^", var_name, "\\."), "", related_cols[idx])
  }))
}

df <- na.omit(df)

