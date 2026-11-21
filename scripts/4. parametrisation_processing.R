# 4. parametrisation_processing.R ----------------------------------------------

library(dplyr)
library(bnlearn)

set.seed(123)


# ------------------------------------------------------------------------------
# 0. Checks
# ------------------------------------------------------------------------------

if (!exists("df")) {
  stop(
    "The object df does not exist. ",
    "Run scripts 1 to 3 before this script."
  )
}

if (!is.data.frame(df)) {
  stop("df must be a data frame.")
}

if (!"Date" %in% names(df)) {
  stop("df does not contain a Date column.")
}

if (!inherits(df$Date, "Date")) {
  stop("df$Date must have class Date.")
}

if (anyNA(df)) {
  stop(
    "df contains missing values. ",
    "Complete 2. imputation.R before running this script."
  )
}

if (anyDuplicated(df$Date)) {
  stop("df contains duplicated dates.")
}

if (any(as.integer(diff(df$Date)) != 1L)) {
  stop(
    "The observations are not consecutive daily observations."
  )
}

if (!exists("candidate_lags")) {
  stop(
    "candidate_lags does not exist. ",
    "Run 3. lag_strategy.R first."
  )
}

if (length(candidate_lags) == 0L ||
    any(!is.finite(candidate_lags)) ||
    any(candidate_lags < 1L)) {
  stop("candidate_lags contains invalid values.")
}

candidate_lags <- sort(
  unique(
    as.integer(candidate_lags)
  )
)

if (!exists("selected_p")) {
  selected_p <- candidate_lags[1]
}

date_column <- df$Date


# ------------------------------------------------------------------------------
# 1. Remove any constant base variables
# ------------------------------------------------------------------------------

base_variable_names <- setdiff(
  names(df),
  "Date"
)

constant_base_flag <- vapply(
  df[
    ,
    base_variable_names,
    drop = FALSE
  ],
  function(variable) {
    length(unique(variable)) <= 1L
  },
  logical(1)
)

constant_base_variables <-
  base_variable_names[
    constant_base_flag
  ]

df_preprocessing_base <- df

if (length(constant_base_variables) > 0L) {
  
  warning(
    "Removing constant variables before structure learning: ",
    paste(
      constant_base_variables,
      collapse = ", "
    )
  )
  
  df_preprocessing_base <-
    df_preprocessing_base[
      ,
      !names(df_preprocessing_base) %in%
        constant_base_variables,
      drop = FALSE
    ]
}

constant_base_variable_audit <- data.frame(
  variable = constant_base_variables,
  reason = rep(
    "constant_in_completed_dataset",
    length(constant_base_variables)
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 2. Encode factor variables using k-1 indicators
# ------------------------------------------------------------------------------

encode_factor_variables <- function(data) {
  
  if (!"Date" %in% names(data)) {
    stop(
      "The data supplied to encode_factor_variables() ",
      "must contain Date."
    )
  }
  
  factor_variables <- names(data)[
    vapply(
      data,
      is.factor,
      logical(1)
    )
  ]
  
  if (length(factor_variables) == 0L) {
    
    empty_mapping_audit <- data.frame(
      source_variable = character(0),
      factor_level = character(0),
      reference_level = logical(0),
      dummy_variable = character(0),
      stringsAsFactors = FALSE
    )
    
    return(
      list(
        data = data,
        factor_dummy_map =
          setNames(
            list(),
            character(0)
          ),
        factor_dummy_audit =
          empty_mapping_audit
      )
    )
  }
  
  for (variable_name in factor_variables) {
    
    if (anyNA(data[[variable_name]])) {
      stop(
        "Factor variable ",
        variable_name,
        " contains missing values."
      )
    }
    
    if (nlevels(data[[variable_name]]) < 2L) {
      stop(
        "Factor variable ",
        variable_name,
        " contains fewer than two levels."
      )
    }
  }
  
  factor_formula <- stats::reformulate(
    factor_variables
  )
  
  # First factor level serves as the reference category.
  
  factor_contrasts <- setNames(
    lapply(
      factor_variables,
      function(variable_name) {
        
        stats::contr.treatment(
          n = levels(
            data[[variable_name]]
          ),
          base = 1L
        )
      }
    ),
    factor_variables
  )
  
  dummy_matrix_full <- stats::model.matrix(
    factor_formula,
    data = data,
    contrasts.arg = factor_contrasts
  )
  
  column_assignment <- attr(
    dummy_matrix_full,
    "assign"
  )
  
  keep_nonintercept <- colnames(
    dummy_matrix_full
  ) != "(Intercept)"
  
  dummy_matrix <- dummy_matrix_full[
    ,
    keep_nonintercept,
    drop = FALSE
  ]
  
  dummy_assignment <- column_assignment[
    keep_nonintercept
  ]
  
  storage.mode(dummy_matrix) <- "double"
  
  if (anyNA(dummy_matrix) ||
      any(!is.finite(dummy_matrix))) {
    stop(
      "The factor dummy matrix contains missing or non-finite values."
    )
  }
  
  factor_dummy_map <- setNames(
    lapply(
      seq_along(factor_variables),
      function(variable_index) {
        
        colnames(dummy_matrix)[
          dummy_assignment ==
            variable_index
        ]
      }
    ),
    factor_variables
  )
  
  # Confirm that k-1 indicators were produced.
  for (variable_name in factor_variables) {
    
    expected_dummy_count <-
      nlevels(data[[variable_name]]) -
      1L
    
    actual_dummy_count <- length(
      factor_dummy_map[[variable_name]]
    )
    
    if (actual_dummy_count !=
        expected_dummy_count) {
      
      stop(
        "Expected ",
        expected_dummy_count,
        " dummy variables for ",
        variable_name,
        " but model.matrix() produced ",
        actual_dummy_count,
        "."
      )
    }
  }
  
  factor_dummy_audit <- dplyr::bind_rows(
    lapply(
      factor_variables,
      function(variable_name) {
        
        factor_levels <- levels(
          data[[variable_name]]
        )
        
        dummy_variables <- factor_dummy_map[[variable_name]]
        
        data.frame(
          source_variable = variable_name,
          factor_level = factor_levels,
          reference_level =
            seq_along(factor_levels) == 1L,
          dummy_variable = c(
            NA_character_,
            dummy_variables
          ),
          stringsAsFactors = FALSE
        )
      }
    )
  )
  
  nonfactor_data <- data[
    ,
    setdiff(
      names(data),
      factor_variables
    ),
    drop = FALSE
  ]
  
  encoded_data <- cbind(
    nonfactor_data,
    as.data.frame(
      dummy_matrix,
      check.names = FALSE
    )
  )
  
  encoded_data <- as.data.frame(
    encoded_data,
    check.names = FALSE
  )
  
  if (anyDuplicated(names(encoded_data))) {
    stop(
      "Duplicated variable names were produced during factor encoding."
    )
  }
  
  if (!inherits(encoded_data$Date, "Date")) {
    stop(
      "The Date class was lost during factor encoding."
    )
  }
  
  analysis_variables <- setdiff(
    names(encoded_data),
    "Date"
  )
  
  numeric_check <- vapply(
    encoded_data[
      ,
      analysis_variables,
      drop = FALSE
    ],
    is.numeric,
    logical(1)
  )
  
  if (!all(numeric_check)) {
    stop(
      "The following encoded analysis variables are not numeric: ",
      paste(
        names(numeric_check)[
          !numeric_check
        ],
        collapse = ", "
      )
    )
  }
  
  list(
    data = encoded_data,
    factor_dummy_map =
      factor_dummy_map,
    factor_dummy_audit =
      factor_dummy_audit
  )
}


encoded_factor_result <-
  encode_factor_variables(
    df_preprocessing_base
  )

df_analysis <-
  encoded_factor_result$data

factor_dummy_map <-
  encoded_factor_result$factor_dummy_map

factor_dummy_audit <-
  encoded_factor_result$factor_dummy_audit

# ------------------------------------------------------------------------------
# 3. Check consistency with the VAR representation
# ------------------------------------------------------------------------------

analysis_node_names <- setdiff(
  names(df_analysis),
  "Date"
)

if (exists("df_var")) {
  
  if (nrow(df_analysis) != nrow(df_var)) {
    stop(
      "df_analysis and df_var contain different numbers of observations."
    )
  }
  
  only_in_structure_data <- setdiff(
    analysis_node_names,
    names(df_var)
  )
  
  only_in_var_data <- setdiff(
    names(df_var),
    analysis_node_names
  )
  
  if (length(only_in_structure_data) > 0L ||
      length(only_in_var_data) > 0L) {
    
    stop(
      "The encoded variables used for structure learning do not ",
      "match those used for the VAR diagnostics. ",
      "Only in structure data: ",
      paste(
        only_in_structure_data,
        collapse = ", "
      ),
      ". Only in VAR data: ",
      paste(
        only_in_var_data,
        collapse = ", "
      ),
      "."
    )
  }
}


# ------------------------------------------------------------------------------
# 4. Fit reproducible two-category k-means discretisers
# ------------------------------------------------------------------------------

fit_kmeans_discretiser <- function(
    x,
    variable_name,
    nstart = 50L
) {
  
  x <- as.numeric(x)
  
  if (anyNA(x)) {
    stop(
      variable_name,
      " contains missing values before discretisation."
    )
  }
  
  if (any(!is.finite(x))) {
    stop(
      variable_name,
      " contains non-finite values before discretisation."
    )
  }
  
  unique_values <- unique(x)
  
  if (length(unique_values) == 1L) {
    
    return(
      list(
        status = "constant",
        lower_centre =
          unique_values[1],
        upper_centre =
          unique_values[1],
        cutpoint = NA_real_
      )
    )
  }
  
  fit_attempt <- tryCatch(
    list(
      fit = stats::kmeans(
        x = x,
        centers = 2L,
        nstart = as.integer(nstart)
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
    stop(
      "k-means discretisation failed for ",
      variable_name,
      ": ",
      fit_attempt$error
    )
  }
  
  centres <- sort(
    as.numeric(
      fit_attempt$fit$centers
    )
  )
  
  if (length(centres) != 2L ||
      any(!is.finite(centres))) {
    stop(
      "Invalid k-means centres were returned for ",
      variable_name,
      "."
    )
  }
  
  cutpoint <- mean(centres)
  
  list(
    status = "two_categories",
    lower_centre = centres[1],
    upper_centre = centres[2],
    cutpoint = cutpoint
  )
}


apply_discretiser <- function(
    x,
    discretiser
) {
  
  x <- as.numeric(x)
  
  if (identical(
    discretiser$status,
    "constant"
  )) {
    
    return(
      factor(
        rep(
          "constant",
          length(x)
        ),
        levels = "constant"
      )
    )
  }
  
  factor(
    ifelse(
      x <= discretiser$cutpoint,
      "low",
      "high"
    ),
    levels = c(
      "low",
      "high"
    )
  )
}


numeric_base_variables <- setdiff(
  names(df_analysis)[
    vapply(
      df_analysis,
      is.numeric,
      logical(1)
    )
  ],
  "Date"
)

if (length(numeric_base_variables) == 0L) {
  stop(
    "No numeric variables are available for discretisation."
  )
}

bn_discretisers <- setNames(
  lapply(
    numeric_base_variables,
    function(variable_name) {
      
      fit_kmeans_discretiser(x = df_analysis[[variable_name]],
        variable_name = variable_name,
        nstart = 50L
      )
    }
  ),
  numeric_base_variables
)

bn_cutpoints <- setNames(
  lapply(
    bn_discretisers,
    function(discretiser) {
      discretiser$cutpoint
    }
  ),
  numeric_base_variables
)

df_bn_base <- df_analysis

for (variable_name in numeric_base_variables) {
  
  df_bn_base[[variable_name]] <-
    apply_discretiser(x = df_analysis[[variable_name]],
      discretiser =
        bn_discretisers[[variable_name]]
    )
}

if (!inherits(df_bn_base$Date, "Date")) {
  stop(
    "The Date class was lost during discretisation."
  )
}

bn_factor_check <- vapply(
  df_bn_base[
    ,
    numeric_base_variables,
    drop = FALSE
  ],
  is.factor,
  logical(1)
)

if (!all(bn_factor_check)) {
  stop(
    "Some discretised variables are not factors."
  )
}


# ------------------------------------------------------------------------------
# 5. Save the discretisation audit
# ------------------------------------------------------------------------------

bn_discretisation_audit <- dplyr::bind_rows(
  lapply(
    numeric_base_variables,
    function(variable_name) {
      
      discretiser <- bn_discretisers[[variable_name]]
      
      discretised_values <-
        df_bn_base[[variable_name]]
      
      category_counts <- table(
        discretised_values
      )
      
      get_category_count <- function(
    category_name
      ) {
        
        if (category_name %in%
            names(category_counts)) {
          
          return(
            as.integer(
              category_counts[[category_name]]
            )
          )
        }
        
        return(0L)
      }
      
      data.frame(
        variable = variable_name,
        status = discretiser$status,
        lower_centre =
          discretiser$lower_centre,
        upper_centre =
          discretiser$upper_centre,
        cutpoint =
          discretiser$cutpoint,
        n_low =
          get_category_count("low"),
        n_high =
          get_category_count("high"),
        n_constant =
          get_category_count("constant"),
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  bn_discretisation_audit,
  file = "bn_discretisation_audit.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Create one lag-expanded DBN dataset
# ------------------------------------------------------------------------------

make_dbn_data <- function(
    data,
    p
) {
  
  if (!is.data.frame(data)) {
    stop("data must be a data frame.")
  }
  
  if (!"Date" %in% names(data)) {
    stop("data must contain a Date column.")
  }
  
  if (!inherits(data$Date, "Date")) {
    stop("data$Date must have class Date.")
  }
  
  if (anyNA(data)) {
    stop(
      "The data supplied to make_dbn_data() contain missing values."
    )
  }
  
  if (anyDuplicated(data$Date)) {
    stop(
      "The data supplied to make_dbn_data() contain duplicated dates."
    )
  }
  
  data <- data[
    order(data$Date),
    ,
    drop = FALSE
  ]
  
  if (any(as.integer(diff(data$Date)) != 1L)) {
    stop(
      "The data supplied to make_dbn_data() are not consecutive daily observations."
    )
  }
  
  if (!is.numeric(p) ||
      length(p) != 1L ||
      !is.finite(p) ||
      p < 1L ||
      p != as.integer(p)) {
    stop("p must be one positive integer.")
  }
  
  p <- as.integer(p)
  
  number_of_rows <- nrow(data)
  
  if (p >= number_of_rows) {
    stop(
      "p must be smaller than the number of observations."
    )
  }
  
  variable_names <- setdiff(
    names(data),
    "Date"
  )
  
  if (length(variable_names) == 0L) {
    stop(
      "No analysis variables remain after removing Date."
    )
  }
  
  outcome_indices <- seq.int(
    from = p + 1L,
    to = number_of_rows
  )
  
  lagged_data <- data.frame(
    Date = data$Date[
      outcome_indices
    ],
    check.names = FALSE
  )
  
  # Add lagged predictors in the order:
  # all variables at t-1, all variables at t-2, ..., all variables at t-p.
  for (lag_value in seq_len(p)) {
    
    predictor_indices <- seq.int(
      from = p + 1L - lag_value,
      to = number_of_rows -
        lag_value
    )
    
    for (variable_name in variable_names) {
      
      lagged_node_name <- paste0(
        variable_name,
        "_t-",
        lag_value
      )
      
      lagged_data[[lagged_node_name]] <- data[[variable_name]][predictor_indices]
    }
  }
  
  # Add current outcomes at time t.
  for (variable_name in variable_names) {
    
    current_node_name <- paste0(
      variable_name,
      "_t"
    )
    
    lagged_data[[current_node_name]] <- data[[variable_name]][outcome_indices]
  }
  
  expected_rows <- number_of_rows - p
  
  expected_nodes <- length(
    variable_names
  ) * (p + 1L)
  
  if (nrow(lagged_data) !=
      expected_rows) {
    stop(
      "Unexpected number of lag-expanded rows for p = ",
      p,
      "."
    )
  }
  
  if ((ncol(lagged_data) - 1L) !=
      expected_nodes) {
    stop(
      "Unexpected number of lag-expanded nodes for p = ",
      p,
      "."
    )
  }
  
  if (anyNA(lagged_data)) {
    
    affected_columns <- names(
      lagged_data
    )[
      colSums(
        is.na(lagged_data)
      ) > 0
    ]
    
    stop(
      "Lag expansion produced missing values for p = ",
      p,
      " in: ",
      paste(
        affected_columns,
        collapse = ", "
      ),
      ". No rows have been silently removed."
    )
  }
  
  if (anyDuplicated(names(lagged_data))) {
    stop(
      "Lag expansion produced duplicated node names for p = ",
      p,
      "."
    )
  }
  
  rownames(lagged_data) <- NULL
  
  return(lagged_data)
}


# ------------------------------------------------------------------------------
# 7. Prepare numeric data for structure learning
# ------------------------------------------------------------------------------

prepare_structure_data <- function(
    lagged_data
) {
  
  if (!"Date" %in% names(lagged_data)) {
    stop(
      "lagged_data must contain Date."
    )
  }
  
  structure_data <- lagged_data[
    ,
    setdiff(
      names(lagged_data),
      "Date"
    ),
    drop = FALSE
  ]
  
  numeric_check <- vapply(
    structure_data,
    is.numeric,
    logical(1)
  )
  
  if (!all(numeric_check)) {
    stop(
      "The following structure-learning nodes are not numeric: ",
      paste(
        names(numeric_check)[
          !numeric_check
        ],
        collapse = ", "
      )
    )
  }
  
  if (anyNA(structure_data) ||
      any(
        !is.finite(
          data.matrix(structure_data)
        )
      )) {
    
    stop(
      "The structure-learning data contain missing or non-finite values."
    )
  }
  
  as.data.frame(
    structure_data,
    check.names = FALSE
  )
}


# ------------------------------------------------------------------------------
# 8. Create numeric and discretised DBN data for every lag
# ------------------------------------------------------------------------------

df_lagged_raw_by_lag <- list()
df_bn_lagged_by_lag <- list()

df_structure_by_lag <- list()
df_bn_by_lag <- list()

empty_net_by_lag <- list()

lagged_data_audit_rows <- list()

for (p in candidate_lags) {
  
  p_name <- paste0(
    "p",
    p
  )
  
  message(
    "Creating lag-expanded datasets for ",
    p_name
  )
  
  # Numeric structure-learning version.
  df_lagged_p <- make_dbn_data(
    data = df_analysis,
    p = p
  )
  
  df_structure_p <- prepare_structure_data(
    df_lagged_p
  )
  
  # Discretised BN parameterisation version.
  df_bn_lagged_p <- make_dbn_data(
    data = df_bn_base,
    p = p
  )
  
  df_bn_p <- df_bn_lagged_p[
    ,
    setdiff(
      names(df_bn_lagged_p),
      "Date"
    ),
    drop = FALSE
  ]
  
  if (!all(
    vapply(
      df_bn_p,
      is.factor,
      logical(1)
    )
  )) {
    
    nonfactor_nodes <- names(
      df_bn_p
    )[
      !vapply(
        df_bn_p,
        is.factor,
        logical(1)
      )
    ]
    
    stop(
      "The following BN nodes are not factors for ",
      p_name,
      ": ",
      paste(
        nonfactor_nodes,
        collapse = ", "
      )
    )
  }
  
  if (anyNA(df_bn_p)) {
    stop(
      "The discretised BN data contain missing values for ",
      p_name,
      "."
    )
  }
  
  if (!identical(
    colnames(df_structure_p),
    colnames(df_bn_p)
  )) {
    
    stop(
      "The numeric and discretised node sets differ for ",
      p_name,
      "."
    )
  }
  
  if (!identical(
    df_lagged_p$Date,
    df_bn_lagged_p$Date
  )) {
    
    stop(
      "The numeric and discretised outcome dates differ for ",
      p_name,
      "."
    )
  }
  
  expected_rows <- nrow(df_analysis) -
    p
  
  if (nrow(df_structure_p) !=
      expected_rows ||
      nrow(df_bn_p) !=
      expected_rows) {
    
    stop(
      "Unexpected number of rows for ",
      p_name,
      "."
    )
  }
  
  node_names <- colnames(
    df_structure_p
  )
  
  lagged_predictor_nodes <- grep(
    "_t-[0-9]+$",
    node_names,
    value = TRUE
  )
  
  current_outcome_nodes <- grep(
    "_t$",
    node_names,
    value = TRUE
  )
  
  expected_lagged_nodes <-
    length(analysis_node_names) * p
  
  expected_current_nodes <-
    length(analysis_node_names)
  
  if (length(lagged_predictor_nodes) !=
      expected_lagged_nodes) {
    
    stop(
      "Unexpected number of lagged predictor nodes for ",
      p_name,
      "."
    )
  }
  
  if (length(current_outcome_nodes) !=
      expected_current_nodes) {
    
    stop(
      "Unexpected number of current outcome nodes for ",
      p_name,
      "."
    )
  }
  
  df_lagged_raw_by_lag[[p_name]] <-
    df_lagged_p
  
  df_bn_lagged_by_lag[[p_name]] <-
    df_bn_lagged_p
  
  df_structure_by_lag[[p_name]] <-
    df_structure_p
  
  df_bn_by_lag[[p_name]] <-
    as.data.frame(
      df_bn_p,
      check.names = FALSE
    )
  
  empty_net_by_lag[[p_name]] <-
    bnlearn::empty.graph(
      node_names
    )
  
  lagged_data_audit_rows[[p_name]] <-
    data.frame(
      lag_name = p_name,
      maximum_lag = p,
      n_input_rows =
        nrow(df_analysis),
      n_output_rows =
        nrow(df_structure_p),
      n_base_nodes =
        length(analysis_node_names),
      n_lagged_predictor_nodes =
        length(lagged_predictor_nodes),
      n_current_outcome_nodes =
        length(current_outcome_nodes),
      n_total_nodes =
        length(node_names),
      numeric_missing_values =
        sum(is.na(df_structure_p)),
      discrete_missing_values =
        sum(is.na(df_bn_p)),
      first_outcome_date =
        as.character(
          df_lagged_p$Date[1]
        ),
      last_outcome_date =
        as.character(
          df_lagged_p$Date[
            nrow(df_lagged_p)
          ]
        ),
      stringsAsFactors = FALSE
    )
}