#' Quantile of the posterior distribution of the treatment effect
#'
#' @param cdf A distribution function (from \code{decision_cdf()}).
#' @param p A probability strictly between 0 and 1.
#' @param support The range of the effect (from \code{effect_support()}).
#'
#' @return The quantile, found by \code{stats::uniroot()} after the search
#'   of a finite interval when the range is infinite.
#' @noRd
effect_quantile <- function(cdf, p, support) {
  lower <- if (is.finite(support[1L])) support[1L] else -1
  upper <- if (is.finite(support[2L])) support[2L] else 1
  k <- 0L
  while (!is.finite(support[1L]) && cdf(lower) > p && k < 100L) {
    lower <- 2 * lower - 1
    k <- k + 1L
  }
  k <- 0L
  while (!is.finite(support[2L]) && cdf(upper) < p && k < 100L) {
    upper <- 2 * upper + 1
    k <- k + 1L
  }
  stats::uniroot(function(q) cdf(q) - p, c(lower, upper),
                 tol = 1e-8 * (upper - lower))$root
}
