#' Cumulative Distribution Function of the Difference Between Two Beta Variables
#'
#' Calculates the cumulative distribution function (CDF) of the difference
#' between two independent Beta-distributed random variables. Specifically,
#' computes \eqn{P(X_1 - X_2 \le q)} or \eqn{P(X_1 - X_2 > q)}, where
#' \eqn{X_1 \sim \mathrm{Beta}(\alpha_1, \beta_1)} and
#' \eqn{X_2 \sim \mathrm{Beta}(\alpha_2, \beta_2)}.
#'
#' @param q A numeric scalar in \code{[-1, 1]} representing the quantile threshold.
#' @param alpha1 A positive numeric scalar giving the first shape parameter of
#'        the first Beta distribution.
#' @param alpha2 A positive numeric scalar giving the first shape parameter of
#'        the second Beta distribution.
#' @param beta1 A positive numeric scalar giving the second shape parameter of
#'        the first Beta distribution.
#' @param beta2 A positive numeric scalar giving the second shape parameter of
#'        the second Beta distribution.
#' @param lower.tail A logical scalar; if \code{TRUE} (default), the function
#'        returns \eqn{P(X_1 - X_2 \le q)}, otherwise \eqn{P(X_1 - X_2 > q)}.
#'
#' @return A numeric scalar in \code{[0, 1]}.
#'
#' @details
#' The upper-tail probability is obtained via the convolution formula:
#' \deqn{P(X_1 - X_2 > q)
#'   = \int_0^1 P(X_2 < x - q)\, f_1(x)\, dx
#'   = \int_0^1 F_{\mathrm{Beta}(\alpha_2,\beta_2)}(x - q)\,
#'              f_{\mathrm{Beta}(\alpha_1,\beta_1)}(x)\, dx}
#' where \eqn{f_1} is the density of \eqn{X_1} and \eqn{F_2} is the CDF of
#' \eqn{X_2}. Boundary cases are handled automatically by \code{pbeta}
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
#' p2betadiff(0.2, 0.5, 0.5, 0.5, 0.5, lower.tail = FALSE)
#'
#' # P(Beta(2, 3) - Beta(1, 4) > -0.1)
#' p2betadiff(-0.1, 2, 3, 1, 4, lower.tail = FALSE)
#'
#' # P(Beta(1, 1) - Beta(1, 1) > 0) - should be approximately 0.5
#' p2betadiff(0, 1, 1, 1, 1, lower.tail = FALSE)
#'
#' # Lower tail: P(Beta(2, 2) - Beta(2, 2) <= 0.1)
#' p2betadiff(0.1, 2, 2, 2, 2, lower.tail = TRUE)
#'
#' @importFrom stats dbeta integrate pbeta
#' @export
p2betadiff <- function(q, alpha1, alpha2, beta1, beta2, lower.tail = TRUE) {

  # --- Input validation ---
  if (!is.numeric(q) || length(q) != 1L || is.na(q) || q < -1 || q > 1) {
    stop("'q' must be a single numeric value in [-1, 1]")
  }

  if (!is.numeric(alpha1) || length(alpha1) != 1L || is.na(alpha1) || alpha1 <= 0) {
    stop("'alpha1' must be a single positive numeric value")
  }

  if (!is.numeric(alpha2) || length(alpha2) != 1L || is.na(alpha2) || alpha2 <= 0) {
    stop("'alpha2' must be a single positive numeric value")
  }

  if (!is.numeric(beta1) || length(beta1) != 1L || is.na(beta1) || beta1 <= 0) {
    stop("'beta1' must be a single positive numeric value")
  }

  if (!is.numeric(beta2) || length(beta2) != 1L || is.na(beta2) || beta2 <= 0) {
    stop("'beta2' must be a single positive numeric value")
  }

  if (!is.logical(lower.tail) || length(lower.tail) != 1L || is.na(lower.tail)) {
    stop("'lower.tail' must be a single logical value (TRUE or FALSE)")
  }

  # --- Compute P(X1 - X2 > q) via convolution ---
  # Integrand: f_Beta(alpha1,beta1)(x) * F_Beta(alpha2,beta2)(x - q)
  # pbeta handles x - q outside [0, 1] automatically.
  result <- integrate(
    function(x) pbeta(x - q, alpha2, beta2) * dbeta(x, alpha1, beta1),
    lower = 0,
    upper = 1
  )[["value"]]

  # Convert to lower tail if requested:
  #   lower.tail = TRUE  -> P(X1 - X2 <= q) = 1 - P(X1 - X2 > q)
  #   lower.tail = FALSE -> P(X1 - X2 > q)
  result <- lower.tail + c(1, -1)[lower.tail + 1] * result

  return(result)
}
