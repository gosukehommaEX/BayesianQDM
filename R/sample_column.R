#' Probabilities of a region for one scenario
#'
#' @param p A vector or matrix of probabilities, or a list of two of them
#'   (marginal rules).
#' @param j Index of the scenario.
#' @param per_scenario Logical; whether \code{p} has one column per scenario
#'   (otherwise \code{p} is the same for all scenarios).
#'
#' @return \code{p} for scenario \code{j}, of the same form.
#' @noRd
sample_column <- function(p, j, per_scenario = is.matrix(p)) {
  if (is.list(p)) {
    return(lapply(p, sample_column, j = j, per_scenario = per_scenario))
  }
  if (per_scenario) p[, j] else p
}
