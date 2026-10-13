#' Probabilities of a region for one scenario
#'
#' @param p A vector (the same for all scenarios) or a matrix with one column
#'   per scenario.
#' @param j Index of the scenario.
#'
#' @return A numeric vector.
#' @noRd
sample_column <- function(p, j) {
  if (is.matrix(p)) p[, j] else p
}
