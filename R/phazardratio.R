#' Distribution function of the ratio of two hazards
#'
#' Computes \eqn{\Pr(\lambda_t / \lambda_c \le q)} (or the upper tail
#' \eqn{\Pr(\lambda_t / \lambda_c > q)}) for hazards with independent Gamma
#' distributions, \eqn{\lambda_t \sim \textrm{Gamma}(a_t, b_t)} and
#' \eqn{\lambda_c \sim \textrm{Gamma}(a_c, b_c)} (shape and rate). This is the
#' posterior probability used for a time-to-event endpoint evaluated by the
#' hazard ratio with a concurrent, external or hypothetical control.
#'
#' Because \eqn{2 b_k \lambda_k} has the chi-squared distribution with
#' \eqn{2 a_k} degrees of freedom, \eqn{(b_t \lambda_t / a_t) /
#' (b_c \lambda_c / a_c)} has the F distribution with \eqn{2 a_t} and
#' \eqn{2 a_c} degrees of freedom, and
#' \deqn{\Pr(\lambda_t / \lambda_c \le q) =
#'   F_{2 a_t, 2 a_c}\left(\frac{q \, a_c \, b_t}{a_t \, b_c}\right),}
#' which is computed exactly by \code{\link[stats]{pf}}.
#'
#' @param q Numeric vector of quantiles (hazard ratios).
#' @param shape_t,rate_t Positive shape and rate of the treatment hazard.
#' @param shape_c,rate_c Positive shape and rate of the control hazard.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are
#'   \eqn{\Pr(\lambda_t / \lambda_c \le q)}, otherwise
#'   \eqn{\Pr(\lambda_t / \lambda_c > q)}.
#'
#' @return A numeric vector of probabilities. All arguments except
#'   \code{lower.tail} are recycled to a common length; each must have length
#'   1 or the common length.
#'
#' @family distributions
#' @examples
#' # Posterior probability that the hazard ratio is below 0.7 with 20 events
#' # in 300 months of follow-up versus 30 events in 250 months (Gamma(0.1, 1)
#' # priors)
#' phazardratio(0.7, shape_t = 20.1, rate_t = 301, shape_c = 30.1,
#'              rate_c = 251)
#' @export
phazardratio <- function(q, shape_t, rate_t, shape_c, rate_c,
                         lower.tail = TRUE) {
  check_numeric(q, "q", len = NULL)
  check_numeric(shape_t, "shape_t", len = NULL, lower = 0, lower_open = TRUE)
  check_numeric(rate_t, "rate_t", len = NULL, lower = 0, lower_open = TRUE)
  check_numeric(shape_c, "shape_c", len = NULL, lower = 0, lower_open = TRUE)
  check_numeric(rate_c, "rate_c", len = NULL, lower = 0, lower_open = TRUE)
  check_flag(lower.tail, "lower.tail")
  args <- recycle_common(list(q = q, shape_t = shape_t, rate_t = rate_t,
                              shape_c = shape_c, rate_c = rate_c))
  x <- pmax(args$q, 0) * args$shape_c * args$rate_t /
    (args$shape_t * args$rate_c)
  stats::pf(x, 2 * args$shape_t, 2 * args$shape_c, lower.tail = lower.tail)
}
