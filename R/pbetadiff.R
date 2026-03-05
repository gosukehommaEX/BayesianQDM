#' Cumulative Distribution Function of the Difference Between Two Beta Variables
#'
#' Calculates the cumulative distribution function (CDF) of the difference
#' between two independent Beta-distributed random variables. Specifically,
#' computes \eqn{P(X_t - X_c \le q)} or \eqn{P(X_t - X_c > q)}, where
#' \eqn{X_t \sim \mathrm{Beta}(\alpha_t, \beta_t)} and
#' \eqn{X_c \sim \mathrm{Beta}(\alpha_c, \beta_c)}.
#'
#' @param q A numeric scalar in \code{[-1, 1]} representing the quantile threshold.
#' @param alpha_t A positive numeric scalar giving the first shape parameter of
#'        the treatment Beta distribution.
#' @param alpha_c A positive numeric scalar giving the first shape parameter of
#'        the control Beta distribution.
#' @param beta_t A positive numeric scalar giving the second shape parameter of
#'        the treatment Beta distribution.
#' @param beta_c A positive numeric scalar giving the second shape parameter of
#'        the control Beta distribution.
#' @param lower.tail A logical scalar; if \code{TRUE} (default), the function
#'        returns \eqn{P(X_t - X_c \le q)}, otherwise \eqn{P(X_t - X_c > q)}.
#'
#' @return A numeric scalar in \code{[0, 1]}.
#'
#' @details
#' The upper-tail probability is obtained via the convolution formula:
#' \deqn{P(X_t - X_c > q)
#'   = \int_0^1 P(X_c < x - q)\, f_t(x)\, dx
#'   = \int_0^1 F_{\mathrm{Beta}(\alpha_c,\beta_c)}(x - q)\,
#'              f_{\mathrm{Beta}(\alpha_t,\beta_t)}(x)\, dx}
#' where \eqn{f_t} is the density of \eqn{X_t} and \eqn{F_c} is the CDF of
#' \eqn{X_c}. Boundary cases are handled automatically by \code{pbeta}
#' (\eqn{0} for \eqn{x - q \le 0}, \eqn{1} for \eqn{x - q \ge 1}), so
#' integrating over \code{[0, 1]} is safe.
#'
#' This single-integral convolution replaces an equivalent double-integral
#' formulation based on Appell's F1 hypergeometric function, yielding the
#' same result with lower computational cost because both \code{pbeta} and
#' \code{dbeta} are implemented in compiled C code.
#'
#' @examples
#' # P(Beta(0.5, 0.5) - Beta(0.5, 0.5) > 0.2)
#' pbetadiff(0.2, 0.5, 0.5, 0.5, 0.5, lower.tail = FALSE)
#'
#' # P(Beta(2, 3) - Beta(1, 4) > -0.1)
#' pbetadiff(-0.1, 2, 1, 3, 4, lower.tail = FALSE)
#'
#' # P(Beta(1, 1) - Beta(1, 1) > 0) - should be approximately 0.5
#' pbetadiff(0, 1, 1, 1, 1, lower.tail = FALSE)
#'
#' # Lower tail: P(Beta(2, 2) - Beta(2, 2) <= 0.1)
#' pbetadiff(0.1, 2, 2, 2, 2, lower.tail = TRUE)
#'
#' @importFrom stats dbeta integrate pbeta
#' @export
pbetadiff <- function(q, alpha_t, alpha_c, beta_t, beta_c, lower.tail = TRUE) {

  # --- Input validation ---
  if (!is.numeric(q) || length(q) != 1L || is.na(q) || q < -1 || q > 1) {
    stop("'q' must be a single numeric value in [-1, 1]")
  }

  if (!is.numeric(alpha_t) || length(alpha_t) != 1L || is.na(alpha_t) || alpha_t <= 0) {
    stop("'alpha_t' must be a single positive numeric value")
  }

  if (!is.numeric(alpha_c) || length(alpha_c) != 1L || is.na(alpha_c) || alpha_c <= 0) {
    stop("'alpha_c' must be a single positive numeric value")
  }

  if (!is.numeric(beta_t) || length(beta_t) != 1L || is.na(beta_t) || beta_t <= 0) {
    stop("'beta_t' must be a single positive numeric value")
  }

  if (!is.numeric(beta_c) || length(beta_c) != 1L || is.na(beta_c) || beta_c <= 0) {
    stop("'beta_c' must be a single positive numeric value")
  }

  if (!is.logical(lower.tail) || length(lower.tail) != 1L || is.na(lower.tail)) {
    stop("'lower.tail' must be a single logical value (TRUE or FALSE)")
  }

  # --- Compute P(X_t - X_c > q) via convolution ---
  # Integrand: f_Beta(alpha_t, beta_t)(x) * F_Beta(alpha_c, beta_c)(x - q)
  # pbeta handles x - q outside [0, 1] automatically.
  result <- integrate(
    function(x) pbeta(x - q, alpha_c, beta_c) * dbeta(x, alpha_t, beta_t),
    lower = 0,
    upper = 1
  )[["value"]]

  # Convert to lower tail if requested:
  #   lower.tail = TRUE  -> P(X_t - X_c <= q) = 1 - P(X_t - X_c > q)
  #   lower.tail = FALSE -> P(X_t - X_c > q)
  result <- lower.tail + c(1, -1)[lower.tail + 1] * result

  return(result)
}
