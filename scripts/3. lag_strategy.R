# 3. lag_strategy.R -------------------------------------------------------------

library(dplyr)
library(vars)
library(lmtest)
library(moments)


# ------------------------------------------------------------------------------
# 0. Checks and lag settings
# ------------------------------------------------------------------------------

if (!exists("df")) {
  stop(
    "The object df does not exist. ",
    "Run 1. read_data.R and 2. imputation.R first."
  )
}

if (!is.data.frame(df)) {
  stop("df must be a data frame.")
}

if (!"Date" %in% names(df)) {
  stop(
    "Date must remain in df so chronological ordering can be verified."
  )
}

if (!inherits(df$Date, "Date")) {
  stop("df$Date must have class Date.")
}

if (anyNA(df)) {
  stop(
    "The input data contain missing values. ",
    "Complete 2. imputation.R before running this script."
  )
}

if (anyDuplicated(df$Date)) {
  stop("The data contain duplicated dates.")
}

if (any(as.integer(diff(df$Date)) != 1L)) {
  stop(
    "The observations are not consecutive daily observations."
  )
}

sensitivity_lags <- c(
  1L,
  7L,
  14L,
  17L,
  21L
)

# Diagnostic-test settings.
diagnostic_alpha <- 0.05
bg_order <- 7L
pacf_lag_max <- 30L


# ------------------------------------------------------------------------------
# 1. Prepare the dataset
# ------------------------------------------------------------------------------

prepare_var_data <- function(data) {
  
  if (!"Date" %in% names(data)) {
    stop("The data do not contain Date.")
  }
  
  data_ordered <- data[
    order(data$Date),
    ,
    drop = FALSE
  ]
  
  predictor_data <- data_ordered[
    ,
    setdiff(
      names(data_ordered),
      "Date"
    ),
    drop = FALSE
  ]
  
  # Identify constant variables before factor encoding.
  original_constant_flag <- vapply(
    predictor_data,
    function(variable) {
      length(unique(variable)) <= 1L
    },
    logical(1)
  )
  
  original_constant_variables <- names(
    predictor_data
  )[original_constant_flag]
  
  if (length(original_constant_variables) > 0L) {
    
    warning(
      "Dropping constant variables from the VAR diagnostics: ",
      paste(
        original_constant_variables,
        collapse = ", "
      )
    )
    
    predictor_data <- predictor_data[
      ,
      !original_constant_flag,
      drop = FALSE
    ]
  }
  
  if (ncol(predictor_data) == 0L) {
    stop(
      "All variables are constant; VAR diagnostics cannot be fitted."
    )
  }
  
  model_formula <- stats::reformulate(
    names(predictor_data)
  )
  
  model_terms <- stats::terms(
    model_formula,
    data = predictor_data
  )
  
  encoded_matrix <- stats::model.matrix(
    model_formula,
    data = predictor_data
  )
  
  encoded_assignment <- attr(
    encoded_matrix,
    "assign"
  )
  
  keep_nonintercept <- colnames(
    encoded_matrix
  ) != "(Intercept)"
  
  encoded_matrix <- encoded_matrix[
    ,
    keep_nonintercept,
    drop = FALSE
  ]
  
  encoded_assignment <- encoded_assignment[
    keep_nonintercept
  ]
  
  storage.mode(encoded_matrix) <- "double"
  
  term_labels <- attr(
    model_terms,
    "term.labels"
  )
  
  source_variables <- term_labels[
    encoded_assignment
  ]
  
  if (anyNA(encoded_matrix) ||
      any(!is.finite(encoded_matrix))) {
    
    problematic_columns <- colnames(
      encoded_matrix
    )[
      colSums(
        is.na(encoded_matrix) |
          !is.finite(encoded_matrix)
      ) > 0
    ]
    
    stop(
      "Missing or non-finite values were found in encoded VAR columns: ",
      paste(
        problematic_columns,
        collapse = ", "
      )
    )
  }
  
  encoded_variable_flag <- apply(
    encoded_matrix,
    MARGIN = 2,
    FUN = function(variable) {
      length(unique(variable)) > 1L
    }
  )
  
  encoding_audit <- data.frame(
    encoded_variable = colnames(
      encoded_matrix
    ),
    source_variable = source_variables,
    retained = encoded_variable_flag,
    stringsAsFactors = FALSE
  )
  
  removed_encoded_variables <- colnames(
    encoded_matrix
  )[!encoded_variable_flag]
  
  if (length(removed_encoded_variables) > 0L) {
    
    warning(
      "Dropping constant encoded variables: ",
      paste(
        removed_encoded_variables,
        collapse = ", "
      )
    )
    
    encoded_matrix <- encoded_matrix[
      ,
      encoded_variable_flag,
      drop = FALSE
    ]
  }
  
  if (ncol(encoded_matrix) == 0L) {
    stop(
      "No nonconstant encoded variables remain for VAR diagnostics."
    )
  }
  
  var_data <- as.data.frame(
    encoded_matrix,
    check.names = FALSE
  )
  
  constant_variable_audit <- data.frame(
    variable = original_constant_variables,
    reason = rep(
      "constant_before_encoding",
      length(original_constant_variables)
    ),
    stringsAsFactors = FALSE
  )
  
  list(
    data = var_data,
    encoding_audit = encoding_audit,
    constant_variable_audit =
      constant_variable_audit
  )
}


prepared_var_data <- prepare_var_data(
  df
)

df_var <- prepared_var_data$data

var_encoding_audit <-
  prepared_var_data$encoding_audit

var_constant_variable_audit <-
  prepared_var_data$constant_variable_audit

var_data_summary <- data.frame(
  n_daily_observations = nrow(df_var),
  n_original_non_date_variables =
    ncol(df) - 1L,
  n_encoded_var_series = ncol(df_var),
  n_constant_original_variables_removed =
    nrow(var_constant_variable_audit),
  n_constant_encoded_variables_removed =
    sum(!var_encoding_audit$retained),
  first_date = as.character(min(df$Date)),
  last_date = as.character(max(df$Date)),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 2. Calculate feasible unrestricted VAR lag limits
# ------------------------------------------------------------------------------

number_of_observations <- nrow(df_var)
number_of_series <- ncol(df_var)
number_of_deterministic_terms <- 1L

# For a VAR(p), the number of usable outcome observations is n-p,
# and each equation has k*p+1 regressors.
#
# Feasibility requires:
#
# n - p - (k*p + 1) > 0.

max_equation_feasible_p <- floor(
  (
    number_of_observations -
      number_of_deterministic_terms -
      1L
  ) /
    (
      number_of_series + 1L
    )
)

# A joint Gaussian VAR likelihood also requires sufficient residual
# degrees of freedom to estimate a full-rank k-by-k covariance matrix:
#
# n - p - (k*p + 1) >= k.

max_joint_likelihood_p <- floor(
  (
    number_of_observations -
      number_of_deterministic_terms -
      number_of_series
  ) /
    (
      number_of_series + 1L
    )
)

max_equation_feasible_p <- max(
  0L,
  as.integer(max_equation_feasible_p)
)

max_joint_likelihood_p <- max(
  0L,
  as.integer(max_joint_likelihood_p)
)

if (max_equation_feasible_p < 1L) {
  stop(
    "There are not enough observations to estimate even ",
    "a VAR(1)."
  )
}

if (max_joint_likelihood_p < 1L) {
  stop(
    "There are not enough observations to calculate ",
    "a full joint VAR likelihood even for VAR(1)."
  )
}

message(
  "Maximum equation-level OLS lag: ",
  max_equation_feasible_p
)

message(
  "Maximum lag with enough residual degrees of freedom ",
  "for a joint VAR likelihood: ",
  max_joint_likelihood_p
)


# ------------------------------------------------------------------------------
# 3. Select a lag using unrestricted joint-VAR information criteria
# ------------------------------------------------------------------------------

max_lag_for_VARselect <- min(
  21L,
  max_joint_likelihood_p
)

lag_selection_attempt <- tryCatch(
  list(
    result = vars::VARselect(
      df_var,
      lag.max = max_lag_for_VARselect,
      type = "const"
    ),
    error = NA_character_
  ),
  error = function(e) {
    list(
      result = NULL,
      error = conditionMessage(e)
    )
  }
)

if (is.null(lag_selection_attempt$result)) {
  stop(
    "VARselect() failed: ",
    lag_selection_attempt$error
  )
}

lag_selection <- lag_selection_attempt$result

lag_selection_criterion <- "SC(n)"

if (!lag_selection_criterion %in%
    names(lag_selection$selection)) {
  
  stop(
    "The requested lag-selection criterion ",
    lag_selection_criterion,
    " was not returned by VARselect()."
  )
}

selected_p <- as.integer(
  lag_selection$selection[[lag_selection_criterion]]
)

if (!is.finite(selected_p) ||
    selected_p < 1L ||
    selected_p > max_lag_for_VARselect) {
  
  stop(
    "VARselect() returned an invalid selected lag: ",
    selected_p
  )
}

varselect_summary <- data.frame(
  criterion = names(
    lag_selection$selection
  ),
  selected_lag = as.integer(
    lag_selection$selection
  ),
  search_maximum_lag =
    max_lag_for_VARselect,
  selected_at_search_boundary =
    as.integer(lag_selection$selection) ==
    max_lag_for_VARselect,
  stringsAsFactors = FALSE
)

selected_p_is_boundary <-
  selected_p == max_lag_for_VARselect

if (selected_p_is_boundary) {
  
  warning(
    lag_selection_criterion,
    " selected the largest jointly feasible lag, p = ",
    selected_p,
    ". "
  )
}


# ------------------------------------------------------------------------------
# 4. Define the graph-learning lag specifications
# ------------------------------------------------------------------------------

candidate_lags <- sort(
  unique(
    c(
      sensitivity_lags,
      selected_p
    )
  )
)

common_max_p <- max(
  candidate_lags
)

message(
  "Graph-learning candidate lags: ",
  paste(
    candidate_lags,
    collapse = ", "
  )
)

message(
  "All diagnostic models will use outcomes beginning after p = ",
  common_max_p,
  " so that their effective samples are aligned."
)


# ------------------------------------------------------------------------------
# 5. Helpers for residual diagnostics
# ------------------------------------------------------------------------------

empty_equation_diagnostics <- function() {
  
  data.frame(
    p = integer(0),
    equation = character(0),
    BG_p_value = numeric(0),
    BG_reject = logical(0),
    BG_error = character(0),
    JB_p_value = numeric(0),
    JB_reject = logical(0),
    JB_error = character(0),
    max_abs_standardised_residual_PACF =
      numeric(0),
    PACF_error = character(0),
    stringsAsFactors = FALSE
  )
}


safe_bg_test <- function(
    equation_fit,
    order
) {
  
  tryCatch(
    list(
      p_value = as.numeric(
        lmtest::bgtest(
          equation_fit,
          order = order
        )$p.value
      ),
      error = NA_character_
    ),
    error = function(e) {
      list(
        p_value = NA_real_,
        error = conditionMessage(e)
      )
    }
  )
}


safe_jb_test <- function(residual_vector) {
  
  residual_vector <- as.numeric(
    residual_vector
  )
  
  residual_vector <- residual_vector[
    is.finite(residual_vector)
  ]
  
  if (length(residual_vector) < 8L) {
    return(
      list(
        p_value = NA_real_,
        error = paste(
          "Fewer than eight finite residuals",
          "were available."
        )
      )
    )
  }
  
  if (stats::sd(residual_vector) == 0) {
    return(
      list(
        p_value = NA_real_,
        error = "The residual series is constant."
      )
    )
  }
  
  tryCatch(
    list(
      p_value = as.numeric(
        moments::jarque.test(
          residual_vector
        )$p.value
      ),
      error = NA_character_
    ),
    error = function(e) {
      list(
        p_value = NA_real_,
        error = conditionMessage(e)
      )
    }
  )
}


safe_residual_pacf <- function(
    residual_vector,
    requested_lag_max
) {
  
  residual_vector <- as.numeric(
    residual_vector
  )
  
  residual_vector <- residual_vector[
    is.finite(residual_vector)
  ]
  
  if (length(residual_vector) < 3L) {
    return(
      list(
        maximum = NA_real_,
        error = paste(
          "Fewer than three finite residuals",
          "were available."
        )
      )
    )
  }
  
  residual_standard_deviation <- stats::sd(
    residual_vector
  )
  
  if (!is.finite(
    residual_standard_deviation
  ) ||
  residual_standard_deviation <= 0) {
    
    return(
      list(
        maximum = NA_real_,
        error = paste(
          "The residual series cannot be",
          "standardised."
        )
      )
    )
  }
  
  standardised_residuals <- (
    residual_vector -
      mean(residual_vector)
  ) / residual_standard_deviation
  
  usable_lag_max <- min(
    as.integer(requested_lag_max),
    length(standardised_residuals) - 1L
  )
  
  if (usable_lag_max < 1L) {
    return(
      list(
        maximum = NA_real_,
        error = "No positive PACF lag was available."
      )
    )
  }
  
  tryCatch({
    
    pacf_values <- stats::pacf(
      standardised_residuals,
      plot = FALSE,
      lag.max = usable_lag_max
    )$acf
    
    pacf_values <- as.numeric(
      pacf_values
    )
    
    pacf_values <- pacf_values[
      is.finite(pacf_values)
    ]
    
    if (length(pacf_values) == 0L) {
      
      list(
        maximum = NA_real_,
        error = paste(
          "No finite PACF estimates were",
          "returned."
        )
      )
      
    } else {
      
      list(
        maximum = max(
          abs(pacf_values)
        ),
        error = NA_character_
      )
    }
    
  }, error = function(e) {
    
    list(
      maximum = NA_real_,
      error = conditionMessage(e)
    )
  })
}


# ------------------------------------------------------------------------------
# 6. Assess one candidate lag
# ------------------------------------------------------------------------------

assess_var_lag <- function(
    data,
    dates,
    p,
    common_max_p,
    bg_order = 7L,
    pacf_lag_max = 30L,
    alpha = 0.05
) {
  
  p <- as.integer(p)
  common_max_p <- as.integer(
    common_max_p
  )
  
  if (p < 1L) {
    stop("p must be positive.")
  }
  
  if (p > common_max_p) {
    stop("p cannot exceed common_max_p.")
  }
  
  if (length(dates) != nrow(data)) {
    stop(
      "The date vector and VAR data have different lengths."
    )
  }
  
  first_input_row <-
    common_max_p - p + 1L
  
  if (first_input_row < 1L ||
      first_input_row > nrow(data)) {
    stop(
      "The aligned first input row is outside the data."
    )
  }
  
  data_fit <- data[
    seq.int(
      first_input_row,
      nrow(data)
    ),
    ,
    drop = FALSE
  ]
  
  dates_fit <- dates[
    seq.int(
      first_input_row,
      length(dates)
    )
  ]
  
  k <- ncol(data_fit)
  n_input <- nrow(data_fit)
  n_effective <- n_input - p
  
  # Each equation contains k*p lagged coefficients and one intercept.
  n_regressors <- k * p + 1L
  residual_df <- n_effective -
    n_regressors
  
  equation_ols_feasible <-
    residual_df > 0L
  
  joint_likelihood_dimension_feasible <-
    residual_df >= k
  
  first_outcome_date <- if (
    n_effective > 0L
  ) {
    as.character(
      dates_fit[p + 1L]
    )
  } else {
    NA_character_
  }
  
  last_outcome_date <- if (
    n_effective > 0L
  ) {
    as.character(
      dates_fit[length(dates_fit)]
    )
  } else {
    NA_character_
  }
  
  make_summary <- function(
    fit_status,
    fit_error = NA_character_,
    likelihood_error = NA_character_,
    residual_covariance_rank = NA_integer_,
    residual_covariance_full_rank = NA,
    joint_likelihood_available = FALSE,
    n_rank_deficient_equations = NA_integer_,
    log_likelihood = NA_real_,
    likelihood_parameters = NA_real_,
    AIC = NA_real_,
    BIC = NA_real_,
    HQ = NA_real_,
    max_root = NA_real_,
    n_BG_tested = 0L,
    n_BG_reject = 0L,
    n_BG_test_errors = 0L,
    n_JB_tested = 0L,
    n_JB_reject = 0L,
    n_JB_test_errors = 0L,
    max_abs_standardised_residual_PACF =
      NA_real_,
    n_PACF_tested = 0L,
    n_PACF_test_errors = 0L,
    pacf_bound = NA_real_) {
    
    data.frame(
      p = p,
      k = k,
      n_input = n_input,
      n_effective = n_effective,
      n_regressors_per_equation =
        n_regressors,
      residual_df_per_equation =
        residual_df,
      equation_ols_feasible =
        equation_ols_feasible,
      joint_likelihood_dimension_feasible =
        joint_likelihood_dimension_feasible,
      residual_covariance_rank =
        residual_covariance_rank,
      residual_covariance_full_rank =
        residual_covariance_full_rank,
      joint_likelihood_available =
        joint_likelihood_available,
      n_rank_deficient_equations =
        n_rank_deficient_equations,
      first_outcome_date =
        first_outcome_date,
      last_outcome_date =
        last_outcome_date,
      log_likelihood = log_likelihood,
      likelihood_parameters =
        likelihood_parameters,
      AIC = AIC,
      BIC = BIC,
      HQ = HQ,
      max_root = max_root,
      n_BG_tested = n_BG_tested,
      n_BG_reject = n_BG_reject,
      n_BG_test_errors =
        n_BG_test_errors,
      n_JB_tested = n_JB_tested,
      n_JB_reject = n_JB_reject,
      n_JB_test_errors =
        n_JB_test_errors,
      max_abs_std_resid_pacf =
        max_abs_standardised_residual_PACF,
      n_PACF_tested =
        n_PACF_tested,
      n_PACF_test_errors =
        n_PACF_test_errors,
      pacf_bound = pacf_bound,
      fit_status = fit_status,
      fit_error = fit_error,
      likelihood_error =
        likelihood_error,
      stringsAsFactors = FALSE
    )
  }
  
  if (!equation_ols_feasible) {
    
    return(
      list(
        summary = make_summary(
          fit_status =
            "not_equation_ols_feasible",
          fit_error = paste0(
            "Equation-level OLS is infeasible: ",
            n_effective,
            " effective observations, ",
            n_regressors,
            " regressors per equation, and ",
            residual_df,
            " residual degrees of freedom."
          )
        ),
        equation_diagnostics =
          empty_equation_diagnostics()
      )
    )
  }
  
  fit_attempt <- tryCatch(
    list(
      fit = vars::VAR(
        data_fit,
        p = p,
        type = "const"
      ),
      error = NA_character_
    ),
    error = function(e) {
      list(
        fit = NULL,
        error = conditionMessage(e)
      )
    }
  )
  
  if (is.null(fit_attempt$fit)) {
    
    return(
      list(
        summary = make_summary(
          fit_status = "fit_failed",
          fit_error = fit_attempt$error
        ),
        equation_diagnostics =
          empty_equation_diagnostics()
      )
    )
  }
  
  fit <- fit_attempt$fit
  
  residual_matrix <- as.matrix(
    residuals(fit)
  )
  
  storage.mode(residual_matrix) <-
    "double"
  
  if (ncol(residual_matrix) != k) {
    stop(
      "The VAR residual matrix has an unexpected number of columns."
    )
  }
  
  residual_covariance_rank <- qr(
    residual_matrix
  )$rank
  
  residual_covariance_full_rank <-
    residual_covariance_rank == k
  
  equation_ranks <- vapply(
    fit$varresult,
    function(equation_fit) {
      equation_fit$rank
    },
    integer(1)
  )
  
  n_rank_deficient_equations <- sum(
    equation_ranks < n_regressors
  )
  
  # --------------------------------------------------------------------------
  # Residual diagnostic tests by equation
  # --------------------------------------------------------------------------
  
  equation_names <- names(
    fit$varresult
  )
  
  if (is.null(equation_names)) {
    equation_names <- colnames(
      residual_matrix
    )
  }
  
  bg_results <- lapply(
    fit$varresult,
    safe_bg_test,
    order = bg_order
  )
  
  jb_results <- lapply(
    seq_len(ncol(residual_matrix)),
    function(column_index) {
      safe_jb_test(
        residual_matrix[
          ,
          column_index
        ]
      )
    }
  )
  
  pacf_results <- lapply(
    seq_len(ncol(residual_matrix)),
    function(column_index) {
      safe_residual_pacf(
        residual_matrix[
          ,
          column_index
        ],
        requested_lag_max =
          pacf_lag_max
      )
    }
  )
  
  bg_p_values <- vapply(
    bg_results,
    function(result) {
      result$p_value
    },
    numeric(1)
  )
  
  bg_errors <- vapply(
    bg_results,
    function(result) {
      result$error
    },
    character(1)
  )
  
  jb_p_values <- vapply(
    jb_results,
    function(result) {
      result$p_value
    },
    numeric(1)
  )
  
  jb_errors <- vapply(
    jb_results,
    function(result) {
      result$error
    },
    character(1)
  )
  
  pacf_maximum_values <- vapply(
    pacf_results,
    function(result) {
      result$maximum
    },
    numeric(1)
  )
  
  pacf_errors <- vapply(
    pacf_results,
    function(result) {
      result$error
    },
    character(1)
  )
  
  equation_diagnostics <- data.frame(
    p = rep(
      p,
      length(equation_names)
    ),
    equation = equation_names,
    BG_p_value = bg_p_values,
    BG_reject =
      is.finite(bg_p_values) &
      bg_p_values < alpha,
    BG_error = bg_errors,
    JB_p_value = jb_p_values,
    JB_reject =
      is.finite(jb_p_values) &
      jb_p_values < alpha,
    JB_error = jb_errors,
    max_abs_standardised_residual_PACF =
      pacf_maximum_values,
    PACF_error = pacf_errors,
    stringsAsFactors = FALSE
  )
  
  n_BG_tested <- sum(
    is.finite(bg_p_values)
  )
  
  n_BG_reject <- sum(
    is.finite(bg_p_values) &
      bg_p_values < alpha
  )
  
  n_BG_test_errors <- sum(
    !is.na(bg_errors)
  )
  
  n_JB_tested <- sum(
    is.finite(jb_p_values)
  )
  
  n_JB_reject <- sum(
    is.finite(jb_p_values) &
      jb_p_values < alpha
  )
  
  n_JB_test_errors <- sum(
    !is.na(jb_errors)
  )
  
  n_PACF_tested <- sum(
    is.finite(pacf_maximum_values)
  )
  
  n_PACF_test_errors <- sum(
    !is.na(pacf_errors)
  )
  
  maximum_pacf <- if (
    n_PACF_tested == 0L
  ) {
    NA_real_
  } else {
    max(
      pacf_maximum_values,
      na.rm = TRUE
    )
  }
  
  pacf_bound <- 1.96 /
    sqrt(nrow(residual_matrix))
  
  # --------------------------------------------------------------------------
  # VAR stability roots
  # --------------------------------------------------------------------------
  
  root_values <- tryCatch(
    vars::roots(fit),
    error = function(e) {
      numeric(0)
    }
  )
  
  maximum_root <- if (
    length(root_values) == 0L ||
    all(!is.finite(root_values))
  ) {
    NA_real_
  } else {
    max(
      Mod(root_values),
      na.rm = TRUE
    )
  }
  
  # --------------------------------------------------------------------------
  # Joint likelihood and information criteria
  # --------------------------------------------------------------------------
  
  log_likelihood_value <- NA_real_
  likelihood_parameter_count <- NA_real_
  aic_value <- NA_real_
  bic_value <- NA_real_
  hq_value <- NA_real_
  likelihood_error <- NA_character_
  joint_likelihood_available <- FALSE
  
  if (!joint_likelihood_dimension_feasible) {
    
    likelihood_error <- paste(
      "Only",
      residual_df,
      "residual degrees of freedom are available for",
      k,
      "series."
    )
    
  } else if (!residual_covariance_full_rank) {
    
    likelihood_error <- paste(
      "The residual matrix has rank",
      residual_covariance_rank,
      "but rank",
      k,
      "is required for a nonsingular covariance matrix."
    )
    
  } else {
    
    likelihood_attempt <- tryCatch(
      list(
        value = stats::logLik(fit),
        error = NA_character_
      ),
      error = function(e) {
        list(
          value = NULL,
          error = conditionMessage(e)
        )
      }
    )
    
    if (is.null(likelihood_attempt$value)) {
      
      likelihood_error <-
        likelihood_attempt$error
      
    } else {
      
      likelihood_object <-
        likelihood_attempt$value
      
      log_likelihood_value <- as.numeric(
        likelihood_object
      )
      
      likelihood_parameter_count <- attr(
        likelihood_object,
        "df"
      )
      
      if (is.null(
        likelihood_parameter_count
      ) ||
      !is.finite(
        likelihood_parameter_count
      )) {
        
        likelihood_parameter_count <-
          k * n_regressors +
          k * (k + 1L) / 2L
      }
      
      likelihood_n <- attr(
        likelihood_object,
        "nobs"
      )
      
      if (is.null(likelihood_n) ||
          !is.finite(likelihood_n)) {
        likelihood_n <-
          nrow(residual_matrix)
      }
      
      if (is.finite(
        log_likelihood_value
      )) {
        
        aic_value <-
          -2 * log_likelihood_value +
          2 * likelihood_parameter_count
        
        bic_value <-
          -2 * log_likelihood_value +
          log(likelihood_n) *
          likelihood_parameter_count
        
        hq_value <-
          -2 * log_likelihood_value +
          2 * log(log(likelihood_n)) *
          likelihood_parameter_count
        
        joint_likelihood_available <-
          TRUE
        
      } else {
        
        likelihood_error <-
          "logLik() returned a non-finite value."
      }
    }
  }
  
  fit_status <- if (
    n_rank_deficient_equations > 0L
  ) {
    "success_with_rank_deficient_equations"
  } else if (
    !joint_likelihood_available
  ) {
    "equation_fit_only"
  } else {
    "success"
  }
  
  list(
    summary = make_summary(
      fit_status = fit_status,
      likelihood_error =
        likelihood_error,
      residual_covariance_rank =
        residual_covariance_rank,
      residual_covariance_full_rank =
        residual_covariance_full_rank,
      joint_likelihood_available =
        joint_likelihood_available,
      n_rank_deficient_equations =
        n_rank_deficient_equations,
      log_likelihood =
        log_likelihood_value,
      likelihood_parameters =
        likelihood_parameter_count,
      AIC = aic_value,
      BIC = bic_value,
      HQ = hq_value,
      max_root = maximum_root,
      n_BG_tested = n_BG_tested,
      n_BG_reject = n_BG_reject,
      n_BG_test_errors =
        n_BG_test_errors,
      n_JB_tested = n_JB_tested,
      n_JB_reject = n_JB_reject,
      n_JB_test_errors =
        n_JB_test_errors,
      max_abs_standardised_residual_PACF =
        maximum_pacf,
      n_PACF_tested =
        n_PACF_tested,
      n_PACF_test_errors =
        n_PACF_test_errors,
      pacf_bound = pacf_bound
    ),
    equation_diagnostics =
      equation_diagnostics
  )
}


# ------------------------------------------------------------------------------
# 7. Run aligned diagnostics for every graph-learning lag
# ------------------------------------------------------------------------------

lag_assessments <- lapply(
  candidate_lags,
  function(p) {
    
    message(
      "Assessing unrestricted VAR diagnostics for p = ",
      p
    )
    
    assess_var_lag(
      data = df_var,
      dates = df$Date,
      p = p,
      common_max_p =
        common_max_p,
      bg_order = bg_order,
      pacf_lag_max =
        pacf_lag_max,
      alpha = diagnostic_alpha
    )
  }
)

names(lag_assessments) <- paste0(
  "p",
  candidate_lags
)

lag_diagnostics <- dplyr::bind_rows(
  lapply(
    lag_assessments,
    function(result) {
      result$summary
    }
  )
)

lag_equation_diagnostics <-
  dplyr::bind_rows(
    lapply(
      lag_assessments,
      function(result) {
        result$equation_diagnostics
      }
    )
  )

write.csv(
  lag_diagnostics,
  file = "lag_diagnostics.csv",
  row.names = FALSE
)

write.csv(
  lag_equation_diagnostics,
  file = "lag_equation_diagnostics.csv",
  row.names = FALSE
)