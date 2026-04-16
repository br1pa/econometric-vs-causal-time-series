# VAR testing -------------------------------------------------------------------

library(dplyr)
library(vars)
library(moments) 
library(tseries)  
library(lmtest)
library(fGarch)

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

date_column <- date_column[-(1:5)]

df <- cbind(Date = date_column, df)

df <- df %>% dplyr::select(-c('Date'))

df <- df %>%
  mutate(across(where(is.factor), as.integer))

lag_selection <- VARselect(df, lag.max = 18, type = "const")

lag_selection$selection[['AIC(n)']]
lag_selection$selection[["HQ(n)"]]
lag_selection$selection[['SC(n)']]
lag_selection$selection[['FPE(n)']]

var_model <- VAR(df, p = 1) 

residuals_var <- residuals(var_model)

# Get the number of variables in the VAR model
num_vars <- ncol(residuals_var)

# Perform Jarque-Bera test on the residuals of each variable in the VAR
for (i in seq_len(num_vars)) {
  series_resid <- residuals_var[, i]
  
  jb_test <- jarque.bera.test(series_resid)
  
  cat("Jarque-Bera test for", colnames(residuals_var)[i], ":\n")
  print(jb_test)
  cat("\n")
}


for (i in seq_len(num_vars)) {
  print(paste("Tests for", colnames(df)[i]))
  
  # Breusch-Godfrey Test
  fit <- lm(var_model$varresult[[i]])
  print("Breusch-Godfrey Test:")
  print(bgtest(fit))
  
  cat("\n")
}

# Initialize vectors to store metrics for each lag.max value
max_lag <- 18
AIC_values <- numeric(max_lag)
HQ_values <- numeric(max_lag)
SC_values <- numeric(max_lag)
FPE_values <- numeric(max_lag)

# Loop through lag.max values from 1 to max_lag
for (lag in 1:max_lag) {
  lag_selection <- VARselect(df, lag.max = lag, type = "const")
  
  # Store the selection metrics for each lag
  AIC_values[lag] <- lag_selection$selection[['AIC(n)']]
  HQ_values[lag] <- lag_selection$selection[['HQ(n)']]
  SC_values[lag] <- lag_selection$selection[['SC(n)']]
  FPE_values[lag] <- lag_selection$selection[['FPE(n)']]
}

# Plotting
plot(1:max_lag, AIC_values, type = "b", col = "blue", xlab = "lag.max", ylab = "Value", 
     main = "Selection Metrics for Different lag.max Values")
lines(1:max_lag, HQ_values, type = "b", col = "red")
lines(1:max_lag, SC_values, type = "b", col = "green")
lines(1:max_lag, FPE_values, type = "b", col = "purple")
legend("topright", legend = c("AIC", "HQ", "SC", "FPE"), col = c("blue", "red", "green", "purple"), lty = 1)

pacf(df$New.infections, main = "Partial Autocorrelation of New Infections", 
     xlab = "Lags", ylab = "PACF")