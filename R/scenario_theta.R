#' Treatment effects of the scenarios
#'
#' Checks the argument \code{theta} of \code{qdm_scenarios()}.
#'
#' @param design An object of class \code{"qdm_design"}.
#' @param theta A numeric vector (one endpoint), a list of two numeric
#'   vectors (two endpoints) or \code{"thresholds"}.
#'
#' @return A list with one numeric vector of effects per endpoint.
#' @noRd
scenario_theta <- function(design, theta) {
  if (is.character(theta)) {
    if (!identical(theta, "thresholds")) {
      stop(paste("'theta' must be numeric (a list of two numeric vectors",
                 "for two endpoints) or \"thresholds\"."), call. = FALSE)
    }
    return(threshold_effects(design$criterion, design$effect))
  }
  if (design$criterion$n_endpoints == 1L) {
    check_numeric(theta, "theta", len = NULL)
    return(list(as.numeric(theta)))
  }
  if (!is.list(theta) || length(theta) != 2L) {
    stop(paste("For two endpoints, 'theta' must be a list of two numeric",
               "vectors (the effects on endpoints 1 and 2) or",
               "\"thresholds\"."), call. = FALSE)
  }
  for (j in 1:2) {
    check_numeric(theta[[j]], sprintf("theta[[%d]]", j), len = NULL)
  }
  lapply(unname(theta), as.numeric)
}
