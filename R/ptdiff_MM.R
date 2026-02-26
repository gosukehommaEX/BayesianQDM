#' Cumulative Distribution Function of the Difference of Two t-Distributed Variables
#' by the Moment-Matching Approximation
#'
#' Calculates the cumulative distribution function (CDF) of the difference between
#' two independent non-standardised t-distributed random variables using the
#' Moment-Matching (MM) approximation. Specifically, computes
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
#'        t-distribution. Must be greater than 4 for the MM approximation.
#' @param nu.t2 A numeric scalar giving the degrees of freedom of the second
#'        t-distribution. Must be greater than 4 for the MM approximation.
#' @param lower.tail A logical scalar; if \code{TRUE} (default), the function
#'        returns \eqn{P(T_1 - T_2 \le q)}, otherwise \eqn{P(T_1 - T_2 > q)}.
#'
#' @return A numeric scalar or vector in \code{[0, 1]}.  When \code{mu.t1},
#'         \code{mu.t2}, \code{sd.t1}, or \code{sd.t2} are vectors of length
#'         \eqn{n}, a vector of length \eqn{n} is returned.
#'
#' @details
#' The difference \eqn{Z = T_1 - T_2} is approximated by a single t-distribution
#' \eqn{Z \approx t(\mu^*, {\sigma^*}^2, \nu^*)} whose parameters are determined
#' by matching the first two even moments of \eqn{Z}:
#' \itemize{
#'   \item \eqn{\mu^* = \mu_1 - \mu_2}
#'   \item \eqn{{\sigma^*}^2 = (\mathrm{Var}_1 + \mathrm{Var}_2)(\nu^* - 2)/\nu^*},
#'         where \eqn{\mathrm{Var}_i = \sigma_i^2 \nu_i / (\nu_i - 2)}.
#'   \item \eqn{\nu^*} is obtained from the fourth-moment equation.
#' }
#'
#' The approximation requires \eqn{\nu_1 > 4} and \eqn{\nu_2 > 4} for finite
#' fourth moments. It is exact in the normal limit (\eqn{\nu \to \infty}) and
#' works well in practice when \eqn{\nu > 10}. Because it reduces to a single
#' call to \code{pt()}, it is orders of magnitude faster than the numerical
#' integration method (\code{\link{ptdiff_NI}}) and fully vectorised.
#'
#' @examples
#' # P(T1 - T2 > 3) with equal parameters
#' ptdiff_MM(q = 3, mu.t1 = 2, mu.t2 = 0, sd.t1 = 1, sd.t2 = 1,
#'           nu.t1 = 17, nu.t2 = 17, lower.tail = FALSE)
#'
#' # P(T1 - T2 > 1) with unequal variances
#' ptdiff_MM(q = 1, mu.t1 = 5, mu.t2 = 3, sd.t1 = 2, sd.t2 = 1.5,
#'           nu.t1 = 10, nu.t2 = 15, lower.tail = FALSE)
#'
#' # P(T1 - T2 > 0) with different degrees of freedom
#' ptdiff_MM(q = 0, mu.t1 = 1, mu.t2 = 1, sd.t1 = 1, sd.t2 = 1,
#'           nu.t1 = 5, nu.t2 = 20, lower.tail = FALSE)
#'
#' # Lower tail: P(T1 - T2 <= 2)
#' ptdiff_MM(q = 2, mu.t1 = 3, mu.t2 = 0, sd.t1 = 1.5, sd.t2 = 1.2,
#'           nu.t1 = 12, nu.t2 = 15, lower.tail = TRUE)
#'
#' # Vectorised usage
#' ptdiff_MM(q = 1, mu.t1 = c(2, 3, 4), mu.t2 = c(0, 1, 2),
#'           sd.t1 = c(1, 1.2, 1.5), sd.t2 = c(1, 1.1, 1.3),
#'           nu.t1 = 10, nu.t2 = 10, lower.tail = FALSE)
#'
#' @importFrom stats pt
#' @export
ptdiff_MM <- function(q, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2,
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

  if (!is.numeric(nu.t1) || length(nu.t1) != 1L || is.na(nu.t1) || nu.t1 <= 4) {
    stop("'nu.t1' must be a single numeric value greater than 4")
  }

  if (!is.numeric(nu.t2) || length(nu.t2) != 1L || is.na(nu.t2) || nu.t2 <= 4) {
    stop("'nu.t2' must be a single numeric value greater than 4")
  }

  if (!is.logical(lower.tail) || length(lower.tail) != 1L || is.na(lower.tail)) {
    stop("'lower.tail' must be a single logical value (TRUE or FALSE)")
  }

  # --- Moment-Matching approximation ---
  # Step 1: Mean of the difference
  mu.diff <- mu.t1 - mu.t2

  # Step 2: Effective degrees of freedom (nu.star) from fourth-moment matching
  var1    <- sd.t1 ^ 2 * nu.t1 / (nu.t1 - 2)
  var2    <- sd.t2 ^ 2 * nu.t2 / (nu.t2 - 2)
  Qu_star <- (var1 + var2) ^ 2
  term1   <- (sd.t1 ^ 2) ^ 2 * nu.t1 ^ 2 / ((nu.t1 - 2) * (nu.t1 - 4))
  term2   <- (sd.t2 ^ 2) ^ 2 * nu.t2 ^ 2 / ((nu.t2 - 2) * (nu.t2 - 4))
  term3   <- 2 * var1 * var2
  Qu      <- term1 + term2 + term3
  nu.star <- (2 * Qu_star - 4 * Qu) / (Qu_star - Qu)

  # Step 3: Scale parameter (sigma.star) from second-moment matching
  sigma.star <- sqrt((var1 + var2) * (nu.star - 2) / nu.star)

  # Step 4: Standardised statistic and CDF evaluation
  standardized <- (q - mu.diff) / sigma.star
  results      <- pt(standardized, df = nu.star, lower.tail = FALSE)

  # Convert to lower tail if requested:
  #   lower.tail = TRUE  -> P(T1 - T2 <= q) = 1 - P(T1 - T2 > q)
  #   lower.tail = FALSE -> P(T1 - T2 > q)
  results <- lower.tail + c(1, -1)[lower.tail + 1] * results

  return(results)
}
