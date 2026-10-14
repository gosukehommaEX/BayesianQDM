#' Integer keys of count vectors
#'
#' @param x Integer matrix of counts with columns \code{"00"}, \code{"01"},
#'   \code{"10"} and \code{"11"} and a common total.
#'
#' @return A numeric vector of keys that differ between different rows: the
#'   counts of the last three columns as digits in base \code{b}, where
#'   \code{b} exceeds the largest count (exact for totals below about
#'   200000).
#' @noRd
outcome_key <- function(x) {
  b <- max(x, 0) + 1
  x[, 2L] + b * x[, 3L] + b ^ 2 * x[, 4L]
}
