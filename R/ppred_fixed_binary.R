#' Posterior predictive distribution of a future effect with a fixed control
#'
#' Probability that the future effect \eqn{Y / m - v} exceeds (or does not
#' exceed) \code{q}, where \eqn{Y} has the Beta-binomial distribution with
#' \code{m} trials and parameters \code{a} and \code{b}, and \eqn{v} is the
#' fixed control response rate. A future effect that differs from \code{q} by
#' less than \eqn{10^{-10}} is treated as equal to \code{q}, so that ties
#' such as \eqn{3 / 10 - 0.3 = 0} are not decided by rounding errors.
#'
#' @param q Null value of the future effect.
#' @param m Number of future subjects.
#' @param a,b Numeric vectors of Beta parameters.
#' @param value Fixed control response rate.
#' @param lower.tail Logical; if \code{TRUE}, \eqn{P(Y / m - v \le q)},
#'   otherwise \eqn{P(Y / m - v > q)}.
#'
#' @return A numeric vector of probabilities, one per element of \code{a}.
#' @noRd
ppred_fixed_binary <- function(q, m, a, b, value, lower.tail = TRUE) {
  k <- 0:m
  above <- (k / m - value) - q > 1e-10
  args <- recycle_common(list(a = a, b = b))
  out <- vapply(seq_along(args$a), function(i) {
    p <- dbetabinom(k, m, args$a[i], args$b[i])
    sum(p[if (lower.tail) !above else above])
  }, numeric(1))
  pmin(pmax(out, 0), 1)
}
