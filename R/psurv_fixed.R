#' Distribution function of a milestone survival probability
#'
#' \eqn{\Pr(S \le x)} (or \eqn{\Pr(S > x)}) for \eqn{S = \exp(-\lambda t_0)}
#' with \eqn{\lambda \sim \textrm{Gamma}(a, b)} (shape and rate). For
#' \eqn{0 < x < 1}, \eqn{S > x} if and only if \eqn{\lambda < -\log(x) / t_0};
#' \eqn{\Pr(S > x)} is 1 for \eqn{x \le 0} and 0 for \eqn{x \ge 1}. Used for
#' a fixed control with \eqn{x} the threshold plus the fixed value.
#'
#' @param x Numeric vector of values.
#' @param milestone Milestone time \eqn{t_0}.
#' @param shape,rate Numeric vectors of Gamma parameters.
#' @param lower.tail Logical; if \code{TRUE}, \eqn{\Pr(S \le x)}.
#'
#' @return A numeric vector of probabilities.
#' @noRd
psurv_fixed <- function(x, milestone, shape, rate, lower.tail = TRUE) {
  args <- recycle_common(list(x = x, shape = shape, rate = rate))
  inside <- args$x > 0 & args$x < 1
  # P(S > x) outside (0, 1)
  upper <- ifelse(args$x <= 0, 1, 0)
  # P(lambda < -log(x) / t0) = P(S > x); the lower tail is the upper tail of
  # the hazard
  lam <- -log(args$x[inside]) / milestone
  if (lower.tail) {
    out <- 1 - upper
    out[inside] <- stats::pgamma(lam, args$shape[inside], args$rate[inside],
                                 lower.tail = FALSE)
  } else {
    out <- upper
    out[inside] <- stats::pgamma(lam, args$shape[inside], args$rate[inside])
  }
  out
}
