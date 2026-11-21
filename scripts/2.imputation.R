# 2. imputation.R ---------------------------------------------------------------

library(imputeTS)

# ------------------------------------------------------------------------------
# 0. Validate the imported dataset
# ------------------------------------------------------------------------------

if (!exists("df")) {
  stop(
    "The object df does not exist. ",
    "Run 1. read_data.R before 2. imputation.R."
  )
}

if (!is.data.frame(df)) {
  stop("df must be a data frame.")
}

if (!"Date" %in% names(df)) {
  stop("The imported dataset does not contain a Date variable.")
}

if (!inherits(df$Date, "Date")) {
  stop(
    "df$Date must have class Date. ",
    "Check 1. read_data.R."
  )
}

if (anyNA(df$Date)) {
  stop("The Date variable contains missing values.")
}

if (anyDuplicated(df$Date)) {
  stop("The Date variable contains duplicated dates.")
}

if (any(as.integer(diff(df$Date)) != 1L)) {
  stop(
    "The imported observations are not consecutive daily observations."
  )
}

initial_number_of_rows <- nrow(df)
initial_number_of_columns <- ncol(df)
initial_dates <- df$Date
initial_missing_count <- sum(is.na(df))

numeric_variables <- setdiff(
  names(df)[
    vapply(
      df,
      is.numeric,
      logical(1)
    )
  ],
  "Date"
)

categorical_variables <- names(df)[
  vapply(
    df,
    is.factor,
    logical(1)
  )
]

expected_variables <- c(
  "Date",
  numeric_variables,
  categorical_variables
)

unexpected_variables <- setdiff(
  names(df),
  expected_variables
)

if (length(unexpected_variables) > 0L) {
  stop(
    "The following variables are neither numeric nor factors: ",
    paste(
      unexpected_variables,
      collapse = ", "
    )
  )
}

missing_before_imputation <- vapply(
  df,
  function(x) {
    sum(is.na(x))
  },
  integer(1)
)


# ------------------------------------------------------------------------------
# 1. Helper for recording imputation status
# ------------------------------------------------------------------------------

make_imputation_status_row <- function(
    variable,
    variable_type,
    status,
    missing_before,
    imputed_values,
    remaining_missing,
    error = NA_character_,
    number_of_levels = NA_integer_,
    tied_reconstructions = NA_integer_
) {
  
  data.frame(
    variable = variable,
    variable_type = variable_type,
    status = status,
    error = error,
    missing_before = as.integer(missing_before),
    imputed_values = as.integer(imputed_values),
    remaining_missing = as.integer(remaining_missing),
    number_of_levels = as.integer(number_of_levels),
    tied_reconstructions =
      as.integer(tied_reconstructions),
    stringsAsFactors = FALSE
  )
}


imputation_status_rows <- list()


# ------------------------------------------------------------------------------
# 2. Impute numeric variables
# ------------------------------------------------------------------------------

for (variable_name in numeric_variables) {
  
  message(
    "Imputing numeric variable: ",
    variable_name
  )
  
  original_values <- as.numeric(
    df[[variable_name]]
  )
  
  missing_rows <- is.na(
    original_values
  )
  
  number_missing_before <- sum(
    missing_rows
  )
  
  if (number_missing_before == 0L) {
    
    imputation_status_rows[[variable_name]] <-
      make_imputation_status_row(
        variable = variable_name,
        variable_type = "numeric",
        status = "no_missing_values",
        missing_before = 0L,
        imputed_values = 0L,
        remaining_missing = 0L
      )
    
    next
  }
  
  number_observed <- sum(
    !missing_rows
  )
  
  if (number_observed == 0L) {
    
    imputation_status_rows[[variable_name]] <-
      make_imputation_status_row(
        variable = variable_name,
        variable_type = "numeric",
        status = "failed",
        missing_before =
          number_missing_before,
        imputed_values = 0L,
        remaining_missing =
          number_missing_before,
        error = "The variable contains no observed values."
      )
    
    next
  }
  
  # If all observed values are identical, an ARIMA state-space model
  # is unnecessary and may fail. Fill missing values with that value.
  observed_values <- original_values[
    !missing_rows
  ]
  
  if (length(unique(observed_values)) == 1L) {
    
    imputed_values <- original_values
    
    imputed_values[missing_rows] <-
      unique(observed_values)
    
    df[[variable_name]] <-
      imputed_values
    
    imputation_status_rows[[variable_name]] <-
      make_imputation_status_row(
        variable = variable_name,
        variable_type = "numeric",
        status = "success_constant_series",
        missing_before =
          number_missing_before,
        imputed_values =
          number_missing_before,
        remaining_missing = 0L
      )
    
    next
  }
  
  imputation_attempt <- tryCatch(
    list(
      values = imputeTS::na_kalman(
        original_values,
        model = "auto.arima",
        smooth = FALSE
      ),
      error = NA_character_
    ),
    error = function(e) {
      list(
        values = NULL,
        error = conditionMessage(e)
      )
    }
  )
  
  if (is.null(imputation_attempt$values)) {
    
    imputation_status_rows[[variable_name]] <-
      make_imputation_status_row(
        variable = variable_name,
        variable_type = "numeric",
        status = "failed",
        missing_before =
          number_missing_before,
        imputed_values = 0L,
        remaining_missing =
          number_missing_before,
        error = imputation_attempt$error
      )
    
    next
  }
  
  imputed_values <- as.numeric(
    imputation_attempt$values
  )
  
  if (length(imputed_values) !=
      length(original_values)) {
    
    imputation_status_rows[[variable_name]] <-
      make_imputation_status_row(
        variable = variable_name,
        variable_type = "numeric",
        status = "failed",
        missing_before =
          number_missing_before,
        imputed_values = 0L,
        remaining_missing =
          number_missing_before,
        error = paste(
          "na_kalman() returned",
          length(imputed_values),
          "values instead of",
          length(original_values),
          "."
        )
      )
    
    next
  }
  
  # Confirm that originally observed values were not altered.
  observed_values_preserved <- isTRUE(
    all.equal(
      imputed_values[!missing_rows],
      original_values[!missing_rows],
      tolerance = sqrt(.Machine$double.eps),
      check.attributes = FALSE
    )
  )
  
  if (!observed_values_preserved) {
    
    imputation_status_rows[[variable_name]] <-
      make_imputation_status_row(
        variable = variable_name,
        variable_type = "numeric",
        status = "failed",
        missing_before =
          number_missing_before,
        imputed_values = 0L,
        remaining_missing =
          number_missing_before,
        error = paste(
          "na_kalman() altered values that were",
          "originally observed."
        )
      )
    
    next
  }
  
  nonfinite_imputations <- missing_rows &
    !is.na(imputed_values) &
    !is.finite(imputed_values)
  
  if (any(nonfinite_imputations)) {
    
    imputation_status_rows[[variable_name]] <-
      make_imputation_status_row(
        variable = variable_name,
        variable_type = "numeric",
        status = "failed",
        missing_before =
          number_missing_before,
        imputed_values = 0L,
        remaining_missing =
          number_missing_before,
        error = paste(
          "na_kalman() returned non-finite",
          "imputed values."
        )
      )
    
    next
  }
  
  df[[variable_name]] <- imputed_values
  
  number_imputed <- sum(
    missing_rows &
      !is.na(imputed_values)
  )
  
  number_remaining <- sum(
    is.na(imputed_values)
  )
  
  numeric_status <- if (
    number_remaining == 0L
  ) {
    "success"
  } else {
    "success_with_remaining_missing"
  }
  
  imputation_status_rows[[variable_name]] <-
    make_imputation_status_row(
      variable = variable_name,
      variable_type = "numeric",
      status = numeric_status,
      missing_before =
        number_missing_before,
      imputed_values =
        number_imputed,
      remaining_missing =
        number_remaining
    )
}


# ------------------------------------------------------------------------------
# 3. Impute categorical variables
# ------------------------------------------------------------------------------

for (variable_name in categorical_variables) {
  
  message(
    "Imputing categorical variable: ",
    variable_name
  )
  
  original_factor <- df[[variable_name]]
  missing_rows <- is.na(original_factor)
  
  number_missing_before <- sum(
    missing_rows
  )
  
  factor_levels <- levels(
    original_factor
  )
  
  number_of_levels <- length(
    factor_levels
  )
  
  if (number_missing_before == 0L) {
    
    imputation_status_rows[[variable_name]] <-
      make_imputation_status_row(
        variable = variable_name,
        variable_type = "categorical",
        status = "no_missing_values",
        missing_before = 0L,
        imputed_values = 0L,
        remaining_missing = 0L,
        number_of_levels =
          number_of_levels,
        tied_reconstructions = 0L
      )
    
    next
  }
  
  if (number_of_levels == 0L) {
    
    imputation_status_rows[[variable_name]] <-
      make_imputation_status_row(
        variable = variable_name,
        variable_type = "categorical",
        status = "failed",
        missing_before =
          number_missing_before,
        imputed_values = 0L,
        remaining_missing =
          number_missing_before,
        error = paste(
          "The factor has no observed levels."
        ),
        number_of_levels = 0L
      )
    
    next
  }
  
  # A one-level factor can be imputed deterministically.
  if (number_of_levels == 1L) {
    
    reconstructed_values <- as.character(
      original_factor
    )
    
    reconstructed_values[missing_rows] <-
      factor_levels[1]
    
    df[[variable_name]] <- factor(
      reconstructed_values,
      levels = factor_levels
    )
    
    imputation_status_rows[[variable_name]] <-
      make_imputation_status_row(
        variable = variable_name,
        variable_type = "categorical",
        status = "success_single_level",
        missing_before =
          number_missing_before,
        imputed_values =
          number_missing_before,
        remaining_missing = 0L,
        number_of_levels = 1L,
        tied_reconstructions = 0L
      )
    
    next
  }
  
  # Create one indicator for every observed factor level.
  dummy_matrix <- do.call(
    cbind,
    lapply(
      factor_levels,
      function(level_name) {
        
        ifelse(
          is.na(original_factor),
          NA_real_,
          as.numeric(
            original_factor ==
              level_name
          )
        )
      }
    )
  )
  
  dummy_matrix <- as.matrix(
    dummy_matrix
  )
  
  storage.mode(dummy_matrix) <- "double"
  
  colnames(dummy_matrix) <- factor_levels
  
  indicator_results <- lapply(
    seq_len(ncol(dummy_matrix)),
    function(column_index) {
      
      level_name <- factor_levels[
        column_index
      ]
      
      indicator_values <- dummy_matrix[
        ,
        column_index
      ]
      
      observed_indicator_values <-
        indicator_values[
          !is.na(indicator_values)
        ]
      
      # A constant indicator can be filled directly.
      if (length(unique(
        observed_indicator_values
      )) == 1L) {
        
        completed_indicator <-
          indicator_values
        
        completed_indicator[
          is.na(completed_indicator)
        ] <- unique(
          observed_indicator_values
        )
        
        return(
          list(
            success = TRUE,
            values = completed_indicator,
            error = NA_character_,
            level = level_name
          )
        )
      }
      
      tryCatch(
        list(
          success = TRUE,
          values = imputeTS::na_kalman(
            indicator_values,
            model = "auto.arima",
            smooth = FALSE
          ),
          error = NA_character_,
          level = level_name
        ),
        error = function(e) {
          list(
            success = FALSE,
            values = NULL,
            error = conditionMessage(e),
            level = level_name
          )
        }
      )
    }
  )
  
  failed_indicators <- !vapply(
    indicator_results,
    function(result) {
      result$success
    },
    logical(1)
  )
  
  if (any(failed_indicators)) {
    
    failure_messages <- vapply(
      indicator_results[
        failed_indicators
      ],
      function(result) {
        paste0(
          result$level,
          ": ",
          result$error
        )
      },
      character(1)
    )
    
    imputation_status_rows[[variable_name]] <-
      make_imputation_status_row(
        variable = variable_name,
        variable_type = "categorical",
        status = "failed",
        missing_before =
          number_missing_before,
        imputed_values = 0L,
        remaining_missing =
          number_missing_before,
        error = paste(
          failure_messages,
          collapse = " | "
        ),
        number_of_levels =
          number_of_levels
      )
    
    next
  }
  
  imputed_dummy_matrix <- do.call(
    cbind,
    lapply(
      indicator_results,
      function(result) {
        as.numeric(result$values)
      }
    )
  )
  
  imputed_dummy_matrix <- as.matrix(
    imputed_dummy_matrix
  )
  
  storage.mode(imputed_dummy_matrix) <-
    "double"
  
  colnames(imputed_dummy_matrix) <-
    factor_levels
  
  if (nrow(imputed_dummy_matrix) !=
      nrow(df) ||
      ncol(imputed_dummy_matrix) !=
      number_of_levels) {
    
    imputation_status_rows[[variable_name]] <-
      make_imputation_status_row(
        variable = variable_name,
        variable_type = "categorical",
        status = "failed",
        missing_before =
          number_missing_before,
        imputed_values = 0L,
        remaining_missing =
          number_missing_before,
        error = paste(
          "The imputed indicator matrix has",
          "unexpected dimensions."
        ),
        number_of_levels =
          number_of_levels
      )
    
    next
  }
  
  nonfinite_indicator_values <-
    !is.na(imputed_dummy_matrix) &
    !is.finite(imputed_dummy_matrix)
  
  if (any(nonfinite_indicator_values)) {
    
    imputation_status_rows[[variable_name]] <-
      make_imputation_status_row(
        variable = variable_name,
        variable_type = "categorical",
        status = "failed",
        missing_before =
          number_missing_before,
        imputed_values = 0L,
        remaining_missing =
          number_missing_before,
        error = paste(
          "The imputed indicator matrix contains",
          "non-finite values."
        ),
        number_of_levels =
          number_of_levels
      )
    
    next
  }

  finite_indicator_rows <- apply(
    imputed_dummy_matrix,
    MARGIN = 1,
    FUN = function(values) {
      all(is.finite(values))
    }
  )
  
  imputable_rows <- missing_rows &
    finite_indicator_rows
  
  reconstructed_values <- as.character(
    original_factor
  )
  
  number_of_ties <- 0L
  
  if (any(imputable_rows)) {
    
    reconstruction_scores <-
      imputed_dummy_matrix[
        imputable_rows,
        ,
        drop = FALSE
      ]
    
    row_maximum <- apply(
      reconstruction_scores,
      MARGIN = 1,
      FUN = max
    )
    
    tie_tolerance <- sqrt(
      .Machine$double.eps
    )
    
    number_at_maximum <- rowSums(
      abs(
        reconstruction_scores -
          row_maximum
      ) <= tie_tolerance
    )
    
    number_of_ties <- sum(
      number_at_maximum > 1L
    )
    
    selected_level_index <- max.col(
      reconstruction_scores,
      ties.method = "first"
    )
    
    reconstructed_values[
      imputable_rows
    ] <- factor_levels[
      selected_level_index
    ]
  }
  
  reconstructed_factor <- factor(
    reconstructed_values,
    levels = factor_levels
  )
  
  # Confirm that originally observed categories remain unchanged.
  observed_categories_preserved <- identical(
    as.character(
      reconstructed_factor[
        !missing_rows
      ]
    ),
    as.character(
      original_factor[
        !missing_rows
      ]
    )
  )
  
  if (!observed_categories_preserved) {
    
    imputation_status_rows[[variable_name]] <-
      make_imputation_status_row(
        variable = variable_name,
        variable_type = "categorical",
        status = "failed",
        missing_before =
          number_missing_before,
        imputed_values = 0L,
        remaining_missing =
          number_missing_before,
        error = paste(
          "Categorical reconstruction altered",
          "originally observed values."
        ),
        number_of_levels =
          number_of_levels,
        tied_reconstructions =
          number_of_ties
      )
    
    next
  }
  
  df[[variable_name]] <- reconstructed_factor
  
  number_imputed <- sum(
    missing_rows &
      !is.na(reconstructed_factor)
  )
  
  number_remaining <- sum(
    is.na(reconstructed_factor)
  )
  
  categorical_status <- if (
    number_remaining == 0L
  ) {
    "success"
  } else {
    "success_with_remaining_missing"
  }
  
  imputation_status_rows[[variable_name]] <-
    make_imputation_status_row(
      variable = variable_name,
      variable_type = "categorical",
      status = categorical_status,
      missing_before =
        number_missing_before,
      imputed_values =
        number_imputed,
      remaining_missing =
        number_remaining,
      number_of_levels =
        number_of_levels,
      tied_reconstructions =
        number_of_ties
    )
}


# ------------------------------------------------------------------------------
# 4. Save the variable-level imputation audit
# ------------------------------------------------------------------------------

imputation_variable_audit <- do.call(
  rbind,
  imputation_status_rows
)

rownames(imputation_variable_audit) <-
  NULL

variable_order <- match(
  imputation_variable_audit$variable,
  names(df)
)

imputation_variable_audit <-
  imputation_variable_audit[
    order(variable_order),
    ,
    drop = FALSE
  ]

write.csv(
  imputation_variable_audit,
  file = "imputation_variable_audit.csv",
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 5. Stop if an imputation method failed
# ------------------------------------------------------------------------------

failed_imputations <-
  imputation_variable_audit$status ==
  "failed"

if (any(failed_imputations)) {
  
  failed_variables <-
    imputation_variable_audit$variable[
      failed_imputations
    ]
  
  stop(
    "Imputation failed for: ",
    paste(
      failed_variables,
      collapse = ", "
    ),
    ". See imputation_variable_audit.csv for details. ",
    "No rows have been removed."
  )
}


# ------------------------------------------------------------------------------
# 6. Identify any rows that remain incomplete
# ------------------------------------------------------------------------------

missing_after_imputation <- vapply(
  df,
  function(x) {
    sum(is.na(x))
  },
  integer(1)
)

remaining_incomplete_rows <- which(
  !complete.cases(df)
)

if (length(remaining_incomplete_rows) > 0L) {
  
  removed_rows_audit <- data.frame(
    row_in_imported_data =
      remaining_incomplete_rows,
    Date = as.character(
      df$Date[
        remaining_incomplete_rows
      ]
    ),
    missing_variables = vapply(
      remaining_incomplete_rows,
      function(row_index) {
        
        paste(
          names(df)[
            is.na(
              df[row_index, , drop = TRUE]
            )
          ],
          collapse = " | "
        )
      },
      character(1)
    ),
    stringsAsFactors = FALSE
  )
  
} else {
  
  removed_rows_audit <- data.frame(
    row_in_imported_data = integer(0),
    Date = character(0),
    missing_variables = character(0),
    stringsAsFactors = FALSE
  )
}

write.csv(
  removed_rows_audit,
  file = "imputation_removed_rows.csv",
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 7. Remove only genuinely unresolved rows
# ------------------------------------------------------------------------------

if (length(remaining_incomplete_rows) > 0L) {
  
  message(
    "Removing ",
    length(remaining_incomplete_rows),
    " unresolved row(s) on date(s): ",
    paste(
      as.character(
        df$Date[
          remaining_incomplete_rows
        ]
      ),
      collapse = ", "
    )
  )
  
  df <- df[
    -remaining_incomplete_rows,
    ,
    drop = FALSE
  ]
}

rownames(df) <- NULL
df <- droplevels(df)


# ------------------------------------------------------------------------------
# 8. Validate the completed dataset
# ------------------------------------------------------------------------------

if (nrow(df) == 0L) {
  stop("No observations remain after imputation.")
}

if (!all(complete.cases(df))) {
  stop(
    "The final dataset still contains missing values."
  )
}

if (anyDuplicated(df$Date)) {
  stop(
    "The final dataset contains duplicated dates."
  )
}

if (nrow(df) > 1L &&
    any(as.integer(diff(df$Date)) != 1L)) {
  
  stop(
    "Removing unresolved rows created a gap in the daily time series. ",
    "Row-based lags would no longer correspond to calendar-day lags."
  )
}

if (ncol(df) != initial_number_of_columns) {
  stop(
    "The number of columns changed during imputation."
  )
}

date_column <- df$Date
