# 5. lars.R --------------------------------------------------------------------
#
# Fits LASSO and Least Angle Regression (LAR) temporal graphs for every
# candidate maximum-lag specification created by
# 4. parametrisation_processing.R.

suppressPackageStartupMessages({
  library(glmnet)
  library(lars)
  library(bnlearn)
  library(dplyr)
})

set.seed(123)


# ------------------------------------------------------------------------------
# 0. Configuration and checks
# ------------------------------------------------------------------------------

validation_fraction <- 0.20
minimum_training_rows <- 30L
minimum_validation_rows <- 10L
lasso_nlambda <- 100L

# Coefficients smaller than this are treated as numerical zero.
coefficient_zero_tolerance <- 1e-8


if (!exists("df_structure_by_lag")) {
  stop(
    "df_structure_by_lag does not exist. Run ",
    "4. parametrisation_processing.R before 5. lars.R."
  )
}

if (!exists("candidate_lags")) {
  stop("candidate_lags does not exist. Run 3. lag_strategy.R first.")
}

candidate_lags <- sort(unique(as.integer(candidate_lags)))

if (length(candidate_lags) == 0L || anyNA(candidate_lags) ||
    any(candidate_lags < 1L)) {
  stop("candidate_lags must contain positive integers.")
}

# ------------------------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------------------------

scalar_or_default <- function(x, default) {
  if (is.null(x) || length(x) == 0L || is.na(x[1])) {
    return(default)
  }
  
  x[1]
}


fit_minmax_scaler <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  
  if (ncol(x) == 0L) {
    stop("The predictor matrix has no columns.")
  }
  
  mins <- apply(x, 2L, min)
  maxs <- apply(x, 2L, max)
  ranges <- maxs - mins
  
  usable <- is.finite(mins) & is.finite(maxs) & ranges > 0
  
  safe_ranges <- ranges
  safe_ranges[!usable] <- 1
  
  list(
    mins = mins,
    ranges = safe_ranges,
    usable = usable
  )
}


apply_minmax_scaler <- function(x, scaler) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  
  if (ncol(x) != length(scaler$mins)) {
    stop("The matrix and scaler have incompatible numbers of columns.")
  }
  
  scaled <- sweep(x, 2L, scaler$mins, FUN = "-")
  scaled <- sweep(scaled, 2L, scaler$ranges, FUN = "/")
  
  storage.mode(scaled) <- "double"
  scaled
}


find_duplicate_columns <- function(x, digits = 14L) {
  x <- as.matrix(x)
  
  if (ncol(x) <= 1L) {
    return(rep(FALSE, ncol(x)))
  }
  
  duplicated(
    as.data.frame(
      t(signif(x, digits = digits)),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  )
}


prepare_blocked_design <- function(x_raw,
                                   validation_fraction,
                                   minimum_training_rows,
                                   minimum_validation_rows) {
  x_raw <- as.matrix(x_raw)
  storage.mode(x_raw) <- "double"
  
  if (anyNA(x_raw) || any(!is.finite(x_raw))) {
    stop("The predictor matrix contains missing or non-finite values.")
  }
  
  n <- nrow(x_raw)
  train_end <- floor(n * (1 - validation_fraction))
  validation_start <- train_end + 1L
  n_validation <- n - train_end
  
  if (train_end < minimum_training_rows) {
    stop("There are too few rows in the chronological training block.")
  }
  
  if (n_validation < minimum_validation_rows) {
    stop("There are too few rows in the chronological validation block.")
  }
  
  train_index <- seq_len(train_end)
  validation_index <- seq.int(validation_start, n)
  
  x_train_raw <- x_raw[train_index, , drop = FALSE]
  scaler <- fit_minmax_scaler(x_train_raw)
  
  zero_variance_names <- colnames(x_raw)[!scaler$usable]
  usable_names <- colnames(x_raw)[scaler$usable]
  
  if (length(usable_names) == 0L) {
    stop("No predictor varies in the chronological training block.")
  }
  
  x_all_scaled <- apply_minmax_scaler(x_raw, scaler)
  x_candidate <- x_all_scaled[, usable_names, drop = FALSE]
  x_train_candidate <- x_candidate[train_index, , drop = FALSE]
  
  duplicate_indicator <- find_duplicate_columns(x_train_candidate)
  duplicate_names <- colnames(x_train_candidate)[duplicate_indicator]
  eligible_names <- colnames(x_train_candidate)[!duplicate_indicator]
  
  if (length(eligible_names) == 0L) {
    stop("No eligible predictors remain after preprocessing.")
  }
  
  x_full <- x_all_scaled[, eligible_names, drop = FALSE]
  x_train <- x_full[train_index, , drop = FALSE]
  x_validation <- x_full[validation_index, , drop = FALSE]
  
  if (anyNA(x_full) || any(!is.finite(x_full))) {
    stop("Scaling produced missing or non-finite predictor values.")
  }
  
  training_rank <- qr(x_train, tol = 1e-7, LAPACK = FALSE)$rank
  
  list(
    x_full = x_full,
    x_train = x_train,
    x_validation = x_validation,
    train_index = train_index,
    validation_index = validation_index,
    train_end = train_end,
    validation_start = validation_start,
    eligible_names = eligible_names,
    zero_variance_names = zero_variance_names,
    duplicate_names = duplicate_names,
    training_rank = as.integer(training_rank),
    n_original_predictors = ncol(x_raw)
  )
}


validation_summary <- function(prediction_matrix, y_validation) {
  prediction_matrix <- as.matrix(prediction_matrix)
  storage.mode(prediction_matrix) <- "double"
  y_validation <- as.numeric(y_validation)
  
  if (nrow(prediction_matrix) != length(y_validation)) {
    stop("The validation predictions and outcome have different lengths.")
  }
  
  if (anyNA(prediction_matrix) || any(!is.finite(prediction_matrix))) {
    stop("Validation predictions contain missing or non-finite values.")
  }
  
  squared_errors <- sweep(
    prediction_matrix,
    1L,
    y_validation,
    FUN = "-"
  )^2
  
  mse <- colMeans(squared_errors)
  
  standard_error <- apply(
    squared_errors,
    2L,
    function(z) {
      if (length(z) <= 1L) {
        return(0)
      }
      
      stats::sd(z) / sqrt(length(z))
    }
  )
  
  if (anyNA(mse) || any(!is.finite(mse))) {
    stop("At least one validation MSE is missing or non-finite.")
  }
  
  standard_error[!is.finite(standard_error)] <- 0
  
  list(
    squared_errors = squared_errors,
    mse = mse,
    standard_error = standard_error
  )
}


select_one_standard_error <- function(mse,
                                      standard_error,
                                      sparse_order) {
  mse <- as.numeric(mse)
  standard_error <- as.numeric(standard_error)
  sparse_order <- as.numeric(sparse_order)
  
  if (length(mse) == 0L ||
      length(standard_error) != length(mse) ||
      length(sparse_order) != length(mse)) {
    stop("Invalid inputs to the one-standard-error selector.")
  }
  
  minimum_index <- which.min(mse)
  cutoff <- mse[minimum_index] + standard_error[minimum_index]
  admissible <- which(mse <= cutoff + sqrt(.Machine$double.eps))
  
  if (length(admissible) == 0L) {
    stop("No path position satisfies the one-standard-error rule.")
  }
  
  selected_index <- admissible[which.min(sparse_order[admissible])]
  
  list(
    selected_index = selected_index,
    minimum_index = minimum_index,
    selected_mse = mse[selected_index],
    minimum_mse = mse[minimum_index],
    standard_error_at_minimum = standard_error[minimum_index],
    one_se_cutoff = cutoff
  )
}


target_status_row <- function(target,
                              p,
                              method,
                              status,
                              failed_stage = NA_character_,
                              error = NA_character_,
                              n_eligible_predictors = NA_integer_,
                              n_selected_predictors = NA_integer_,
                              max_absolute_coefficient = NA_real_,
                              target_handling = "regular_fit") {
  data.frame(
    target = target,
    p = as.integer(p),
    method = method,
    status = status,
    failed_stage = failed_stage,
    error = error,
    n_eligible_predictors = as.integer(n_eligible_predictors),
    n_selected_predictors = as.integer(n_selected_predictors),
    max_absolute_coefficient = as.numeric(max_absolute_coefficient),
    target_handling = target_handling,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


tuning_row <- function(target,
                       p,
                       method,
                       engine,
                       path_mode,
                       selected_path_value,
                       validation_mse,
                       minimum_validation_mse,
                       validation_se_at_minimum,
                       one_se_cutoff,
                       design,
                       n_selected_predictors,
                       max_absolute_coefficient,
                       target_handling = "regular_fit") {
  data.frame(
    target = target,
    p = as.integer(p),
    method = method,
    engine = engine,
    path_mode = path_mode,
    selection_rule = "one_standard_error",
    selected_path_value = as.numeric(selected_path_value),
    validation_mse = as.numeric(validation_mse),
    minimum_validation_mse = as.numeric(minimum_validation_mse),
    validation_se_at_minimum = as.numeric(validation_se_at_minimum),
    one_se_cutoff = as.numeric(one_se_cutoff),
    train_end = as.integer(design$train_end),
    validation_start = as.integer(design$validation_start),
    n_eligible_predictors = length(design$eligible_names),
    n_removed_zero_variance = length(design$zero_variance_names),
    n_removed_duplicates = length(design$duplicate_names),
    training_rank = as.integer(design$training_rank),
    n_selected_predictors = as.integer(n_selected_predictors),
    max_absolute_coefficient = as.numeric(max_absolute_coefficient),
    target_handling = target_handling,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


failed_target_result <- function(target,
                                 p,
                                 method,
                                 stage,
                                 error,
                                 design = NULL) {
  if (inherits(error, "condition")) {
    error <- conditionMessage(error)
  }
  
  n_eligible <- if (is.null(design)) {
    NA_integer_
  } else {
    length(design$eligible_names)
  }
  
  list(
    status = "failed",
    coefficients = NULL,
    tuning = NULL,
    status_row = target_status_row(
      target = target,
      p = p,
      method = method,
      status = "failed",
      failed_stage = stage,
      error = as.character(error),
      n_eligible_predictors = n_eligible,
      target_handling = "failed"
    )
  )
}


constant_training_target_result <- function(y_raw,
                                            target,
                                            p,
                                            method,
                                            design) {
  y_training <- y_raw[design$train_index]
  y_validation <- y_raw[design$validation_index]
  
  intercept_prediction <- rep(mean(y_training), length(y_validation))
  validation_mse <- mean((y_validation - intercept_prediction)^2)
  
  coefficients <- stats::setNames(
    rep(0, length(design$eligible_names)),
    design$eligible_names
  )
  
  list(
    status = "success_no_edges",
    coefficients = coefficients,
    tuning = tuning_row(
      target = target,
      p = p,
      method = method,
      engine = "intercept_only",
      path_mode = "none",
      selected_path_value = 0,
      validation_mse = validation_mse,
      minimum_validation_mse = validation_mse,
      validation_se_at_minimum = 0,
      one_se_cutoff = validation_mse,
      design = design,
      n_selected_predictors = 0L,
      max_absolute_coefficient = 0,
      target_handling = "constant_training_intercept_only"
    ),
    status_row = target_status_row(
      target = target,
      p = p,
      method = method,
      status = "success_no_edges",
      n_eligible_predictors = length(design$eligible_names),
      n_selected_predictors = 0L,
      max_absolute_coefficient = 0,
      target_handling = "constant_training_intercept_only"
    )
  )
}


# ------------------------------------------------------------------------------
# 2. LASSO target fit
# ------------------------------------------------------------------------------

fit_lasso_target <- function(y_raw,
                             target,
                             p,
                             design,
                             coefficient_tolerance,
                             nlambda) {
  stage <- "target_validation"
  
  tryCatch({
    y_raw <- as.numeric(y_raw)
    
    if (length(y_raw) != nrow(design$x_full)) {
      stop("The target and predictor matrix have different row counts.")
    }
    
    if (anyNA(y_raw) || any(!is.finite(y_raw))) {
      stop("The target contains missing or non-finite values.")
    }
    
    y_training_raw <- y_raw[design$train_index]
    
    if (length(unique(y_training_raw)) <= 1L) {
      return(
        constant_training_target_result(
          y_raw = y_raw,
          target = target,
          p = p,
          method = "lasso",
          design = design
        )
      )
    }
    
    y_min <- min(y_training_raw)
    y_range <- max(y_training_raw) - y_min
    
    if (!is.finite(y_range) || y_range <= 0) {
      stop("The target has no usable training-period range.")
    }
    
    y_scaled <- (y_raw - y_min) / y_range
    y_training <- y_scaled[design$train_index]
    y_validation <- y_scaled[design$validation_index]
    
    stage <- "training_path"
    
    lambda_min_ratio <- if (nrow(design$x_train) < ncol(design$x_train)) {
      1e-2
    } else {
      1e-4
    }
    
    training_fit <- glmnet::glmnet(
      x = design$x_train,
      y = y_training,
      family = "gaussian",
      alpha = 1,
      standardize = FALSE,
      intercept = TRUE,
      nlambda = nlambda,
      lambda.min.ratio = lambda_min_ratio,
      thresh = 1e-8,
      maxit = 100000L
    )
    
    lambdas <- as.numeric(training_fit$lambda)
    
    if (length(lambdas) == 0L || anyNA(lambdas) ||
        any(!is.finite(lambdas))) {
      stop("glmnet did not return a valid lambda path.")
    }
    
    stage <- "blocked_validation"
    
    validation_predictions <- predict(
      training_fit,
      newx = design$x_validation,
      s = lambdas,
      type = "response",
      exact = FALSE
    )
    
    losses <- validation_summary(
      prediction_matrix = validation_predictions,
      y_validation = y_validation
    )
    
    # glmnet returns lambda from largest (sparsest) to smallest (densest).
    sparse_order <- rank(-lambdas, ties.method = "first")
    
    selection <- select_one_standard_error(
      mse = losses$mse,
      standard_error = losses$standard_error,
      sparse_order = sparse_order
    )
    
    selected_lambda <- lambdas[selection$selected_index]
    
    stage <- "full_refit"
    
    full_fit <- glmnet::glmnet(
      x = design$x_full,
      y = y_scaled,
      family = "gaussian",
      alpha = 1,
      standardize = FALSE,
      intercept = TRUE,
      lambda = selected_lambda,
      thresh = 1e-8,
      maxit = 100000L
    )
    
    stage <- "coefficient_extraction"
    
    coefficient_matrix <- as.matrix(
      stats::coef(
        full_fit,
        s = selected_lambda,
        exact = FALSE
      )
    )
    
    coefficient_vector <- coefficient_matrix[, 1L]
    coefficient_vector <- coefficient_vector[
      names(coefficient_vector) != "(Intercept)"
    ]
    
    coefficients <- stats::setNames(
      rep(0, length(design$eligible_names)),
      design$eligible_names
    )
    
    matched_names <- intersect(
      names(coefficient_vector),
      names(coefficients)
    )
    
    coefficients[matched_names] <- coefficient_vector[matched_names]
    
    if (anyNA(coefficients) || any(!is.finite(coefficients))) {
      stop("The fitted LASSO coefficients are missing or non-finite.")
    }
    
    coefficients[abs(coefficients) <= coefficient_tolerance] <- 0
    
    n_selected <- sum(abs(coefficients) > coefficient_tolerance)
    max_coefficient <- if (length(coefficients) == 0L) {
      0
    } else {
      max(abs(coefficients))
    }
    
    target_status <- if (n_selected == 0L) {
      "success_no_edges"
    } else {
      "success"
    }
    
    list(
      status = target_status,
      coefficients = coefficients,
      tuning = tuning_row(
        target = target,
        p = p,
        method = "lasso",
        engine = "glmnet",
        path_mode = "lambda",
        selected_path_value = selected_lambda,
        validation_mse = selection$selected_mse,
        minimum_validation_mse = selection$minimum_mse,
        validation_se_at_minimum = selection$standard_error_at_minimum,
        one_se_cutoff = selection$one_se_cutoff,
        design = design,
        n_selected_predictors = n_selected,
        max_absolute_coefficient = max_coefficient
      ),
      status_row = target_status_row(
        target = target,
        p = p,
        method = "lasso",
        status = target_status,
        n_eligible_predictors = length(design$eligible_names),
        n_selected_predictors = n_selected,
        max_absolute_coefficient = max_coefficient
      )
    )
  }, error = function(e) {
    failed_target_result(
      target = target,
      p = p,
      method = "lasso",
      stage = stage,
      error = e,
      design = design
    )
  })
}


# ------------------------------------------------------------------------------
# 3. LAR target fit
# ------------------------------------------------------------------------------

lar_path_predictions <- function(fit, newx) {
  newx <- as.matrix(newx)
  storage.mode(newx) <- "double"
  
  beta_path <- as.matrix(fit$beta)
  storage.mode(beta_path) <- "double"
  
  if (ncol(newx) != ncol(beta_path)) {
    stop("The LAR coefficient path and new predictor matrix are incompatible.")
  }
  
  if (length(fit$meanx) != ncol(newx)) {
    stop("The LAR centring vector has the wrong length.")
  }
  
  centred_newx <- sweep(newx, 2L, fit$meanx, FUN = "-")
  predictions <- centred_newx %*% t(beta_path)
  predictions <- predictions + as.numeric(fit$mu)
  
  predictions
}


lar_step_numbers <- function(fit) {
  step_numbers <- suppressWarnings(as.integer(rownames(fit$beta)))
  
  if (length(step_numbers) != nrow(fit$beta) || anyNA(step_numbers)) {
    step_numbers <- seq_len(nrow(fit$beta)) - 1L
  }
  
  step_numbers
}


fit_lar_target <- function(y_raw,
                           target,
                           p,
                           design,
                           coefficient_tolerance) {
  stage <- "target_validation"
  
  tryCatch({
    y_raw <- as.numeric(y_raw)
    
    if (length(y_raw) != nrow(design$x_full)) {
      stop("The target and predictor matrix have different row counts.")
    }
    
    if (anyNA(y_raw) || any(!is.finite(y_raw))) {
      stop("The target contains missing or non-finite values.")
    }
    
    y_training_raw <- y_raw[design$train_index]
    
    if (length(unique(y_training_raw)) <= 1L) {
      return(
        constant_training_target_result(
          y_raw = y_raw,
          target = target,
          p = p,
          method = "lar",
          design = design
        )
      )
    }
    
    y_min <- min(y_training_raw)
    y_range <- max(y_training_raw) - y_min
    
    if (!is.finite(y_range) || y_range <= 0) {
      stop("The target has no usable training-period range.")
    }
    
    y_scaled <- (y_raw - y_min) / y_range
    y_training <- y_scaled[design$train_index]
    y_validation <- y_scaled[design$validation_index]
    
    stage <- "training_path"
    
    training_max_steps <- min(
      ncol(design$x_train),
      nrow(design$x_train) - 1L
    )
    
    training_fit <- lars::lars(
      x = design$x_train,
      y = y_training,
      type = "lar",
      normalize = FALSE,
      intercept = TRUE,
      use.Gram = FALSE,
      max.steps = training_max_steps
    )
    
    beta_training <- as.matrix(training_fit$beta)
    
    if (nrow(beta_training) == 0L ||
        ncol(beta_training) != length(design$eligible_names)) {
      stop("lars did not return a valid coefficient path.")
    }
    
    stage <- "blocked_validation"
    
    validation_predictions <- lar_path_predictions(
      fit = training_fit,
      newx = design$x_validation
    )
    
    step_numbers <- lar_step_numbers(training_fit)
    
    finite_path_rows <- apply(beta_training, 1L, function(z) {
      all(is.finite(z))
    }) & apply(validation_predictions, 2L, function(z) {
      all(is.finite(z))
    })
    
    if (!any(finite_path_rows)) {
      stop("No finite LAR path positions are available for validation.")
    }
    
    beta_training <- beta_training[finite_path_rows, , drop = FALSE]
    validation_predictions <- validation_predictions[
      ,
      finite_path_rows,
      drop = FALSE
    ]
    step_numbers <- step_numbers[finite_path_rows]
    
    path_order <- order(step_numbers)
    beta_training <- beta_training[path_order, , drop = FALSE]
    validation_predictions <- validation_predictions[
      ,
      path_order,
      drop = FALSE
    ]
    step_numbers <- step_numbers[path_order]
    
    losses <- validation_summary(
      prediction_matrix = validation_predictions,
      y_validation = y_validation
    )
    
    selection <- select_one_standard_error(
      mse = losses$mse,
      standard_error = losses$standard_error,
      sparse_order = step_numbers
    )
    
    selected_step <- step_numbers[selection$selected_index]
    
    stage <- "full_refit"
    
    full_max_steps <- min(
      ncol(design$x_full),
      nrow(design$x_full) - 1L
    )
    
    full_fit <- lars::lars(
      x = design$x_full,
      y = y_scaled,
      type = "lar",
      normalize = FALSE,
      intercept = TRUE,
      use.Gram = FALSE,
      max.steps = full_max_steps
    )
    
    full_beta_path <- as.matrix(full_fit$beta)
    full_step_numbers <- lar_step_numbers(full_fit)
    
    eligible_full_rows <- which(
      full_step_numbers <= selected_step &
        apply(full_beta_path, 1L, function(z) all(is.finite(z)))
    )
    
    if (length(eligible_full_rows) == 0L) {
      stop("The selected LAR step is unavailable in the full-sample refit.")
    }
    
    selected_full_row <- eligible_full_rows[
      which.max(full_step_numbers[eligible_full_rows])
    ]
    
    used_full_step <- full_step_numbers[selected_full_row]
    coefficients <- as.numeric(full_beta_path[selected_full_row, ])
    names(coefficients) <- design$eligible_names
    
    if (anyNA(coefficients) || any(!is.finite(coefficients))) {
      stop("The fitted LAR coefficients are missing or non-finite.")
    }
    
    coefficients[abs(coefficients) <= coefficient_tolerance] <- 0
    
    n_selected <- sum(abs(coefficients) > coefficient_tolerance)
    max_coefficient <- if (length(coefficients) == 0L) {
      0
    } else {
      max(abs(coefficients))
    }
    
    target_status <- if (n_selected == 0L) {
      "success_no_edges"
    } else {
      "success"
    }
    
    list(
      status = target_status,
      coefficients = coefficients,
      tuning = tuning_row(
        target = target,
        p = p,
        method = "lar",
        engine = "lars",
        path_mode = "step",
        selected_path_value = used_full_step,
        validation_mse = selection$selected_mse,
        minimum_validation_mse = selection$minimum_mse,
        validation_se_at_minimum = selection$standard_error_at_minimum,
        one_se_cutoff = selection$one_se_cutoff,
        design = design,
        n_selected_predictors = n_selected,
        max_absolute_coefficient = max_coefficient
      ),
      status_row = target_status_row(
        target = target,
        p = p,
        method = "lar",
        status = target_status,
        n_eligible_predictors = length(design$eligible_names),
        n_selected_predictors = n_selected,
        max_absolute_coefficient = max_coefficient
      )
    )
  }, error = function(e) {
    failed_target_result(
      target = target,
      p = p,
      method = "lar",
      stage = stage,
      error = e,
      design = design
    )
  })
}


# ------------------------------------------------------------------------------
# 4. Coefficient matrix to temporal graph
# ------------------------------------------------------------------------------

generate_graph_from_lagged_coefs <- function(coef_df,
                                             threshold = coefficient_zero_tolerance) {
  coefficient_matrix <- as.matrix(coef_df)
  storage.mode(coefficient_matrix) <- "double"
  
  predictor_nodes <- colnames(coefficient_matrix)
  target_nodes <- rownames(coefficient_matrix)
  
  if (is.null(predictor_nodes) || is.null(target_nodes)) {
    stop("The coefficient matrix must have row and column names.")
  }
  
  if (anyNA(coefficient_matrix) || any(!is.finite(coefficient_matrix))) {
    stop("The coefficient matrix contains missing or non-finite values.")
  }
  
  graph <- bnlearn::empty.graph(unique(c(predictor_nodes, target_nodes)))
  
  selected_locations <- which(
    abs(coefficient_matrix) > threshold,
    arr.ind = TRUE
  )
  
  if (nrow(selected_locations) > 0L) {
    for (i in seq_len(nrow(selected_locations))) {
      graph <- bnlearn::set.arc(
        graph,
        from = predictor_nodes[selected_locations[i, "col"]],
        to = target_nodes[selected_locations[i, "row"]],
        check.cycles = FALSE
      )
    }
  }
  
  graph
}


# ------------------------------------------------------------------------------
# 5. Fit one method at one maximum lag
# ------------------------------------------------------------------------------

fit_method_at_lag <- function(df_structure_p,
                              p,
                              method,
                              validation_fraction = 0.20,
                              coefficient_tolerance = coefficient_zero_tolerance) {
  if (!method %in% c("lasso", "lar")) {
    stop("method must be either 'lasso' or 'lar'.")
  }
  
  stage <- "input_validation"
  
  tryCatch({
    x_cols <- grep("_t-[0-9]+$", names(df_structure_p), value = TRUE)
    y_cols <- grep("_t$", names(df_structure_p), value = TRUE)
    
    if (length(x_cols) == 0L) {
      stop("No lagged predictor columns were found.")
    }
    
    if (length(y_cols) == 0L) {
      stop("No current target columns were found.")
    }
    
    x_raw <- data.matrix(df_structure_p[, x_cols, drop = FALSE])
    y_raw <- df_structure_p[, y_cols, drop = FALSE]
    
    if (anyNA(x_raw) || anyNA(y_raw)) {
      stop("The LASSO/LAR input contains missing values.")
    }
    
    if (any(!is.finite(x_raw))) {
      stop("The LASSO/LAR predictor matrix contains non-finite values.")
    }
    
    stage <- "design_preparation"
    
    design <- prepare_blocked_design(
      x_raw = x_raw,
      validation_fraction = validation_fraction,
      minimum_training_rows = minimum_training_rows,
      minimum_validation_rows = minimum_validation_rows
    )
    
    coefficient_matrix <- matrix(
      0,
      nrow = length(y_cols),
      ncol = length(x_cols),
      dimnames = list(y_cols, x_cols)
    )
    
    target_results <- vector("list", length(y_cols))
    names(target_results) <- y_cols
    
    stage <- "target_fitting"
    
    for (target in y_cols) {
      message(
        "Fitting ",
        toupper(method),
        " target ",
        target,
        " at p = ",
        p
      )
      
      if (method == "lasso") {
        target_result <- fit_lasso_target(
          y_raw = y_raw[[target]],
          target = target,
          p = p,
          design = design,
          coefficient_tolerance = coefficient_tolerance,
          nlambda = lasso_nlambda
        )
      } else {
        target_result <- fit_lar_target(
          y_raw = y_raw[[target]],
          target = target,
          p = p,
          design = design,
          coefficient_tolerance = coefficient_tolerance
        )
      }
      
      target_results[[target]] <- target_result
      
      if (target_result$status %in% c("success", "success_no_edges")) {
        coefficient_matrix[
          target,
          names(target_result$coefficients)
        ] <- target_result$coefficients
      }
    }
    
    target_status <- dplyr::bind_rows(
      lapply(target_results, function(x) x$status_row)
    )
    
    tuning <- dplyr::bind_rows(
      lapply(target_results, function(x) x$tuning)
    )
    
    n_failed_targets <- sum(target_status$status == "failed")
    n_successful_targets <- sum(
      target_status$status %in% c("success", "success_no_edges")
    )
    
    if (n_failed_targets > 0L) {
      failed_rows <- target_status[target_status$status == "failed", ]
      
      error_text <- paste0(
        n_failed_targets,
        " of ",
        length(y_cols),
        " target models failed. First failure: ",
        failed_rows$target[1],
        " [",
        failed_rows$failed_stage[1],
        "]: ",
        failed_rows$error[1]
      )
      
      return(
        new_fit_result(
          status = "failed",
          error = error_text,
          graph = NULL,
          coefficients = NULL,
          diagnostics = list(
            method = method,
            p = as.integer(p),
            failed_stage = "target_fitting",
            n_targets = length(y_cols),
            n_successful_targets = n_successful_targets,
            n_failed_targets = n_failed_targets,
            n_selected_edges = NA_integer_,
            n_graph_edges = NA_integer_,
            target_status = target_status,
            tuning = tuning,
            design = list(
              train_end = design$train_end,
              validation_start = design$validation_start,
              n_original_predictors = design$n_original_predictors,
              n_eligible_predictors = length(design$eligible_names),
              zero_variance_predictors = design$zero_variance_names,
              duplicate_predictors = design$duplicate_names,
              training_rank = design$training_rank
            )
          )
        )
      )
    }
    
    stage <- "graph_construction"
    
    coefficient_df <- as.data.frame(
      coefficient_matrix,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    
    attr(coefficient_df, "tuning") <- tuning
    attr(coefficient_df, "target_status") <- target_status
    
    graph <- generate_graph_from_lagged_coefs(
      coefficient_df,
      threshold = coefficient_tolerance
    )
    
    n_selected_edges <- sum(
      abs(coefficient_matrix) > coefficient_tolerance
    )
    n_graph_edges <- nrow(bnlearn::arcs(graph))
    
    if (n_selected_edges != n_graph_edges) {
      stop(
        "The number of selected coefficients does not equal the number of ",
        "graph arcs."
      )
    }
    
    run_status <- if (n_graph_edges == 0L) {
      "success_no_edges"
    } else {
      "success"
    }
    
    new_fit_result(
      status = run_status,
      error = NA_character_,
      graph = graph,
      coefficients = coefficient_df,
      diagnostics = list(
        method = method,
        p = as.integer(p),
        failed_stage = NA_character_,
        n_targets = length(y_cols),
        n_successful_targets = n_successful_targets,
        n_failed_targets = 0L,
        n_constant_training_targets = sum(
          target_status$target_handling ==
            "constant_training_intercept_only"
        ),
        n_selected_edges = as.integer(n_selected_edges),
        n_graph_edges = as.integer(n_graph_edges),
        target_status = target_status,
        tuning = tuning,
        design = list(
          train_end = design$train_end,
          validation_start = design$validation_start,
          n_original_predictors = design$n_original_predictors,
          n_eligible_predictors = length(design$eligible_names),
          zero_variance_predictors = design$zero_variance_names,
          duplicate_predictors = design$duplicate_names,
          training_rank = design$training_rank
        )
      )
    )
  }, error = function(e) {
    new_fit_result(
      status = "failed",
      error = conditionMessage(e),
      graph = NULL,
      coefficients = NULL,
      diagnostics = list(
        method = method,
        p = as.integer(p),
        failed_stage = stage,
        n_targets = NA_integer_,
        n_successful_targets = NA_integer_,
        n_failed_targets = NA_integer_,
        n_selected_edges = NA_integer_,
        n_graph_edges = NA_integer_,
        target_status = data.frame(),
        tuning = data.frame()
      )
    )
  })
}


# ------------------------------------------------------------------------------
# 6. Output helpers
# ------------------------------------------------------------------------------

fit_status_row <- function(result, lag_name, method) {
  diagnostics <- result$diagnostics
  
  data.frame(
    lag = lag_name,
    method = method,
    status = scalar_or_default(result$status, "failed"),
    failed_stage = scalar_or_default(
      diagnostics$failed_stage,
      NA_character_
    ),
    error = scalar_or_default(result$error, NA_character_),
    n_targets = as.integer(
      scalar_or_default(diagnostics$n_targets, NA_integer_)
    ),
    n_successful_targets = as.integer(
      scalar_or_default(diagnostics$n_successful_targets, NA_integer_)
    ),
    n_failed_targets = as.integer(
      scalar_or_default(diagnostics$n_failed_targets, NA_integer_)
    ),
    n_constant_training_targets = as.integer(
      scalar_or_default(
        diagnostics$n_constant_training_targets,
        NA_integer_
      )
    ),
    n_selected_edges = as.integer(
      scalar_or_default(diagnostics$n_selected_edges, NA_integer_)
    ),
    n_graph_edges = as.integer(
      scalar_or_default(diagnostics$n_graph_edges, NA_integer_)
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}


write_successful_fit_outputs <- function(result, lag_name, method) {
  target_status_file <- paste0(method, "_target_status_", lag_name, ".csv")
  tuning_file <- paste0(method, "_tuning_", lag_name, ".csv")
  
  utils::write.csv(
    result$diagnostics$target_status,
    file = target_status_file,
    row.names = FALSE,
    na = "NA"
  )
  
  utils::write.csv(
    result$diagnostics$tuning,
    file = tuning_file,
    row.names = FALSE,
    na = "NA"
  )
  
  if (fit_succeeded(result)) {
    utils::write.csv(
      result$coefficients,
      file = paste0(method, "_coefficients_", lag_name, ".csv"),
      row.names = TRUE,
      na = "NA"
    )
    
    bnlearn::write.dot(
      result$graph,
      file = paste0(method, "_", lag_name, "_dot_file.dot")
    )
  }
}


# ------------------------------------------------------------------------------
# 7. Run all candidate lag specifications
# ------------------------------------------------------------------------------

lasso_fit_results_by_lag <- list()
lar_fit_results_by_lag <- list()
lasso_results_by_lag <- list()
lar_results_by_lag <- list()
lasso_graphs_by_lag <- list()
lar_graphs_by_lag <- list()
lasso_tuning_by_lag <- list()
lar_tuning_by_lag <- list()


for (p in candidate_lags) {
  p_name <- paste0("p", p)
  
  message("============================================================")
  message("Running LASSO and LAR for ", p_name)
  message("============================================================")
  
  if (!p_name %in% names(df_structure_by_lag)) {
    missing_error <- paste0(
      "No structure-learning dataset exists for ",
      p_name,
      "."
    )
    
    lasso_result <- new_fit_result(
      status = "failed",
      error = missing_error,
      graph = NULL,
      coefficients = NULL,
      diagnostics = list(
        method = "lasso",
        p = as.integer(p),
        failed_stage = "input_lookup",
        n_targets = NA_integer_,
        n_successful_targets = NA_integer_,
        n_failed_targets = NA_integer_,
        n_selected_edges = NA_integer_,
        n_graph_edges = NA_integer_,
        target_status = data.frame(),
        tuning = data.frame()
      )
    )
    
    lar_result <- new_fit_result(
      status = "failed",
      error = missing_error,
      graph = NULL,
      coefficients = NULL,
      diagnostics = list(
        method = "lar",
        p = as.integer(p),
        failed_stage = "input_lookup",
        n_targets = NA_integer_,
        n_successful_targets = NA_integer_,
        n_failed_targets = NA_integer_,
        n_selected_edges = NA_integer_,
        n_graph_edges = NA_integer_,
        target_status = data.frame(),
        tuning = data.frame()
      )
    )
  } else {
    df_structure_p <- df_structure_by_lag[[p_name]]
    
    lasso_result <- fit_method_at_lag(
      df_structure_p = df_structure_p,
      p = p,
      method = "lasso",
      validation_fraction = validation_fraction
    )
    
    lar_result <- fit_method_at_lag(
      df_structure_p = df_structure_p,
      p = p,
      method = "lar",
      validation_fraction = validation_fraction
    )
  }
  
  lasso_fit_results_by_lag[[p_name]] <- lasso_result
  lar_fit_results_by_lag[[p_name]] <- lar_result
  
  write_successful_fit_outputs(
    result = lasso_result,
    lag_name = p_name,
    method = "lasso"
  )
  
  write_successful_fit_outputs(
    result = lar_result,
    lag_name = p_name,
    method = "lar"
  )
  
  lasso_results_by_lag[[p_name]] <- if (fit_succeeded(lasso_result)) {
    lasso_result$coefficients
  } else {
    NULL
  }
  
  lar_results_by_lag[[p_name]] <- if (fit_succeeded(lar_result)) {
    lar_result$coefficients
  } else {
    NULL
  }
  
  lasso_graphs_by_lag[[p_name]] <- if (fit_succeeded(lasso_result)) {
    lasso_result$graph
  } else {
    NULL
  }
  
  lar_graphs_by_lag[[p_name]] <- if (fit_succeeded(lar_result)) {
    lar_result$graph
  } else {
    NULL
  }
  
  lasso_tuning_by_lag[[p_name]] <- lasso_result$diagnostics$tuning
  lar_tuning_by_lag[[p_name]] <- lar_result$diagnostics$tuning
}


# ------------------------------------------------------------------------------
# 8. Combined status files and saved result objects
# ------------------------------------------------------------------------------

lasso_status_by_lag <- dplyr::bind_rows(
  lapply(names(lasso_fit_results_by_lag), function(p_name) {
    fit_status_row(
      result = lasso_fit_results_by_lag[[p_name]],
      lag_name = p_name,
      method = "lasso"
    )
  })
)

lar_status_by_lag <- dplyr::bind_rows(
  lapply(names(lar_fit_results_by_lag), function(p_name) {
    fit_status_row(
      result = lar_fit_results_by_lag[[p_name]],
      lag_name = p_name,
      method = "lar"
    )
  })
)

lasso_target_status <- dplyr::bind_rows(
  lapply(lasso_fit_results_by_lag, function(result) {
    result$diagnostics$target_status
  })
)

lar_target_status <- dplyr::bind_rows(
  lapply(lar_fit_results_by_lag, function(result) {
    result$diagnostics$target_status
  })
)


utils::write.csv(
  lasso_status_by_lag,
  file = "lasso_status_by_lag.csv",
  row.names = FALSE,
  na = "NA"
)

utils::write.csv(
  lar_status_by_lag,
  file = "lar_status_by_lag.csv",
  row.names = FALSE,
  na = "NA"
)

utils::write.csv(
  lasso_target_status,
  file = "lasso_target_status.csv",
  row.names = FALSE,
  na = "NA"
)

utils::write.csv(
  lar_target_status,
  file = "lar_target_status.csv",
  row.names = FALSE,
  na = "NA"
)


run_configuration <- data.frame(
  setting = c(
    "validation_fraction",
    "minimum_training_rows",
    "minimum_validation_rows",
    "lasso_nlambda",
    "coefficient_zero_tolerance",
    "selection_rule",
    "constant_training_target_handling",
    "lar_validation_prediction"
  ),
  value = c(
    as.character(validation_fraction),
    as.character(minimum_training_rows),
    as.character(minimum_validation_rows),
    as.character(lasso_nlambda),
    format(coefficient_zero_tolerance, scientific = TRUE),
    "one_standard_error",
    "intercept_only_zero_coefficients",
    "direct_from_beta_path"
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  run_configuration,
  file = "lars_run_configuration.csv",
  row.names = FALSE
)


saveRDS(
  lasso_fit_results_by_lag,
  file = "lasso_fit_results_by_lag.rds"
)

saveRDS(
  lar_fit_results_by_lag,
  file = "lar_fit_results_by_lag.rds"
)


message("LASSO and LAR processing is complete.")
print(lasso_status_by_lag)
print(lar_status_by_lag)
