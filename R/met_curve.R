#' Probability that a criterion is met, as a function of the threshold
#'
#' @param p Probabilities of a region, one per outcome: a numeric vector, a
#'   matrix over pairs of treatment and control outcomes, or a list of two
#'   (one per endpoint) for the marginal rules.
#' @param w Probabilities of the outcomes: a numeric vector, or a list with
#'   elements \code{t} and \code{c} for a matrix \code{p} over pairs (see
#'   \code{sample_weights()}).
#' @param grid Numeric vector of thresholds.
#' @param rule One of \code{"joint"}, \code{"and"} or \code{"or"}; for the
#'   marginal rules, the same threshold is used for both endpoints.
#'
#' @return A numeric vector of the same length as \code{grid}.
#' @noRd
met_curve <- function(p, w, grid, rule = "joint") {
  k <- if (rule == "joint") 1L else 2L
  vapply(grid, function(g) {
    ind <- meets_criterion(p, rep(g, k), rule)
    if (is.list(w)) {
      sum(w$t * (ind %*% w$c))
    } else {
      sum(w[ind])
    }
  }, numeric(1))
}
