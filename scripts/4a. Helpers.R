# 4a.helpers.R -------------------------------------------------------

new_fit_result <- function(status,
                           graph = NULL,
                           coefficients = NULL,
                           diagnostics = list(),
                           error = NA_character_) {
  
  valid_statuses <- c(
    "success",
    "success_no_edges",
    "failed"
  )
  
  if (length(status) != 1L || !status %in% valid_statuses) {
    stop(
      "status must be one of: ",
      paste(valid_statuses, collapse = ", ")
    )
  }
  
  if (inherits(error, "condition")) {
    error <- conditionMessage(error)
  }
  
  structure(
    list(
      status = status,
      error = error,
      graph = graph,
      coefficients = coefficients,
      diagnostics = diagnostics
    ),
    class = c("fit_result", "list")
  )
}


failed_fit <- function(error,
                       method,
                       p,
                       stage) {
  
  new_fit_result(
    status = "failed",
    error = error,
    graph = NULL,
    coefficients = NULL,
    diagnostics = list(
      method = method,
      p = p,
      failed_stage = stage
    )
  )
}


fit_succeeded <- function(result) {
  
  inherits(result, "fit_result") &&
    result$status %in% c(
      "success",
      "success_no_edges"
    )
}