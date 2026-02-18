#' Cumulative Distribution Function of the Difference of Two t-Distributed Variables
#' by Numerical Integration
#'
#' Calculates the cumulative distribution function (CDF) of the difference between
#' two independent non-standardised t-distributed random variables using exact
#' numerical integration via convolution. Specifically, computes
#' \eqn{P(T_1 - T_2 \le q)} or \eqn{P(T_1 - T_2 > q)}, where
#' \eqn{T_i \sim t(\mu_i, \sigma_i^2, \nu_i)}.
#'
#' @param q A numeric scalar representing the quantile threshold.
#' @param mu.t1 A numeric scalar or vector giving the location parameter of the
#'        first t-distribution.
#' @param mu.t2 A numeric scalar or vector giving the location parameter of the
#'        second t-distribution.
#' @param sd.t1 A positive numeric scalar or vector giving the scale parameter of
#'        the first t-distribution.
#' @param sd.t2 A positive numeric scalar or vector giving the scale parameter of
#'        the second t-distribution.
#' @param nu.t1 A numeric scalar giving the degrees of freedom of the first
#'        t-distribution. Must be greater than 2 for finite variance.
#' @param nu.t2 A numeric scalar giving the degrees of freedom of the second
#'        t-distribution. Must be greater than 2 for finite variance.
#' @param lower.tail A logical scalar; if \code{TRUE} (default), the function
#'        returns \eqn{P(T_1 - T_2 \le q)}, otherwise \eqn{P(T_1 - T_2 > q)}.
#'
#' @return A numeric scalar or vector in \code{[0, 1]}.  When \code{mu.t1},
#'         \code{mu.t2}, \code{sd.t1}, or \code{sd.t2} are vectors of length
#'         \eqn{n}, a vector of length \eqn{n} is returned.
#'
#' @details
#' The upper-tail probability is obtained via the convolution formula:
#' \deqn{P(T_1 - T_2 > q)
#'   = \int_{-\infty}^{\infty} f_1(x)\, F_2(x - q)\, dx}
#' where \eqn{f_1} is the density of \eqn{T_1} and \eqn{F_2} is the CDF of
#' \eqn{T_2}. Integration bounds are set to \eqn{\mu_1 \pm 8\sigma^*_1}, where
#' \eqn{\sigma^*_1 = \sigma_1 \sqrt{(\nu_1 + 1)/(\nu_1 - 2)}} is the standard
#' deviation of \eqn{T_1}. Relative and absolute tolerances are \code{1e-6} and
#' \code{1e-8}, respectively.
#'
#' When the input parameters are vectors, \code{mapply} applies the scalar
#' integration function across all parameter sets.
#'
#' Advantages:
#' \itemize{
#'   \item Exact within numerical precision.
#'   \item Handles arbitrary parameter combinations.
#' }
#' Computational note: this method is substantially slower than the
#' Moment-Matching approximation (\code{\link{pMM2tdiff}}) because it calls
#' \code{integrate()} once per parameter set. For large-scale simulation
#' (many parameter sets), prefer \code{CalcMethod = 'MM'} in
#' \code{\link{pGNGsinglecontinuous}}.
#'
#' @examples
#' # P(T1 - T2 > 3) with equal parameters
#' pNI2tdiff(q = 3, mu.t1 = 2, mu.t2 = 0, sd.t1 = 1, sd.t2 = 1,
#'           nu.t1 = 17, nu.t2 = 17, lower.tail = FALSE)
#'
#' # P(T1 - T2 > 1) with unequal variances
#' pNI2tdiff(q = 1, mu.t1 = 5, mu.t2 = 3, sd.t1 = 2, sd.t2 = 1.5,
#'           nu.t1 = 10, nu.t2 = 15, lower.tail = FALSE)
#'
#' # P(T1 - T2 > 0) with different degrees of freedom
#' pNI2tdiff(q = 0, mu.t1 = 1, mu.t2 = 1, sd.t1 = 1, sd.t2 = 1,
#'           nu.t1 = 5, nu.t2 = 20, lower.tail = FALSE)
#'
#' # Lower tail: P(T1 - T2 <= 2)
#' pNI2tdiff(q = 2, mu.t1 = 3, mu.t2 = 0, sd.t1 = 1.5, sd.t2 = 1.2,
#'           nu.t1 = 12, nu.t2 = 15, lower.tail = TRUE)
#'
#' # Vectorised usage
#' pNI2tdiff(q = 1, mu.t1 = c(2, 3, 4), mu.t2 = c(0, 1, 2),
#'           sd.t1 = c(1, 1.2, 1.5), sd.t2 = c(1, 1.1, 1.3),
#'           nu.t1 = 10, nu.t2 = 10, lower.tail = FALSE)
#'
#' @importFrom stats dt integrate pt
#' @export
pNI2tdiff <- function(q, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2,
                      lower.tail = TRUE) {

  # --- Input validation ---
  if (!is.numeric(q) || length(q) != 1L || is.na(q)) {
    stop("'q' must be a single numeric value")
  }

  for (nm in c("mu.t1", "mu.t2", "sd.t1", "sd.t2")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) < 1L || any(is.na(val))) {
      stop(paste0("'", nm, "' must be a numeric scalar or vector with no missing values"))
    }
  }

  if (any(sd.t1 <= 0)) stop("'sd.t1' must contain only positive values")
  if (any(sd.t2 <= 0)) stop("'sd.t2' must contain only positive values")

  if (!is.numeric(nu.t1) || length(nu.t1) != 1L || is.na(nu.t1) || nu.t1 <= 2) {
    stop("'nu.t1' must be a single numeric value greater than 2")
  }

  if (!is.numeric(nu.t2) || length(nu.t2) != 1L || is.na(nu.t2) || nu.t2 <= 2) {
    stop("'nu.t2' must be a single numeric value greater than 2")
  }

  if (!is.logical(lower.tail) || length(lower.tail) != 1L || is.na(lower.tail)) {
    stop("'lower.tail' must be a single logical value (TRUE or FALSE)")
  }

  # --- Scalar integration function for one parameter set ---
  .integrate_one <- function(mu1, mu2, sd1, sd2) {
    # Integration bounds: mu1 +/- 8 * SD(T1)
    spread1 <- sd1 * sqrt((nu.t1 + 1) / (nu.t1 - 2))
    integrate(
      function(x) {
        # PDF of T1 evaluated at x
        f1 <- dt((x - mu1) / sd1, df = nu.t1) / sd1
        # CDF of T2 evaluated at x - q
        F2 <- pt(((x - q) - mu2) / sd2, df = nu.t2)
        f1 * F2
      },
      lower   = mu1 - 8 * spread1,
      upper   = mu1 + 8 * spread1,
      rel.tol = 1e-6,
      abs.tol = 1e-8
    )[["value"]]
  }

  # --- Apply integration across all parameter sets ---
  results <- mapply(.integrate_one,
                    mu1 = mu.t1, mu2 = mu.t2, sd1 = sd.t1, sd2 = sd.t2)

  # Convert to lower tail if requested:
  #   lower.tail = TRUE  -> P(T1 - T2 <= q) = 1 - P(T1 - T2 > q)
  #   lower.tail = FALSE -> P(T1 - T2 > q)
  results <- lower.tail + c(1, -1)[lower.tail + 1] * results

  return(results)
}
