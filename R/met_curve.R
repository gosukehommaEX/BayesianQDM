#' Probability that a criterion is met, as a function of the threshold
#'
#' @param p Numeric vector of probabilities of a region, one per outcome.
#' @param w Numeric vector of the probabilities of the outcomes.
#' @param grid Numeric vector of thresholds.
#'
#' @return A numeric vector of the same length as \code{grid}.
#' @noRd
met_curve <- function(p, w, grid) {
  vapply(grid, function(g) sum(w[meets_threshold(p, g)]), numeric(1))
}
