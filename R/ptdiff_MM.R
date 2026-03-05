#' Cumulative Distribution Function of the Difference of Two t-Distributed Variables
#' by the Moment-Matching Approximation
#'
#' Calculates the cumulative distribution function (CDF) of the difference between
#' two independent non-standardised t-distributed random variables using the
#' Moment-Matching (MM) approximation. Specifically, computes
#' \eqn{P(T_t - T_c \le q)} or \eqn{P(T_t - T_c > q)}, where
#' \eqn{T_k \sim t(\mu_k, \sigma_k^2, \nu_k)} for \eqn{k \in \{t, c\}}.
#'
#' @param q A numeric scalar representing the quantile threshold.
#' @param mu_t A numeric scalar or vector giving the location parameter of the
#'        first t-distribution (treatment).
#' @param mu_c A numeric scalar or vector giving the location parameter of the
#'        second t-distribution (control).
#' @param sd_t A positive numeric scalar or vector giving the scale parameter of
#'        the first t-distribution (treatment).
#' @param sd_c A positive numeric scalar or vector giving the scale parameter of
#'        the second t-distribution (control).
#' @param nu_t A numeric scalar giving the degrees of freedom of the first
#'        t-distribution (treatment). Must be greater than 4 for the MM
#'        approximation.
#' @param nu_c A numeric scalar giving the degrees of freedom of the second
#'        t-distribution (control). Must be greater than 4 for the MM
#'        approximation.
#' @param lower.tail A logical scalar; if \code{TRUE} (default), the function
#'        returns \eqn{P(T_t - T_c \le q)}, otherwise \eqn{P(T_t - T_c > q)}.
#'
#' @return A numeric scalar or vector in \code{[0, 1]}.  When \code{mu_t},
#'         \code{mu_c}, \code{sd_t}, or \code{sd_c} are vectors of length
#'         \eqn{n}, a vector of length \eqn{n} is returned.
#'
#' @details
#' The difference \eqn{Z = T_t - T_c} is approximated by a single t-distribution
#' \eqn{Z \approx t(\mu^*, {\sigma^*}^2, \nu^*)} whose parameters are determined
#' by matching the first two even moments of \eqn{Z}:
#' \itemize{
#'   \item \eqn{\mu^* = \mu_t - \mu_c}
#'   \item \eqn{{\sigma^*}^2 = (\mathrm{Var}_t + \mathrm{Var}_c)(\nu^* - 2)/\nu^*},
#'         where \eqn{\mathrm{Var}_k = \sigma_k^2 \nu_k / (\nu_k - 2)}.
#'   \item \eqn{\nu^*} is obtained from the fourth-moment equation.
#' }
#'
#' The approximation requires \eqn{\nu_t > 4} and \eqn{\nu_c > 4} for finite
#' fourth moments. It is exact in the normal limit (\eqn{\nu \to \infty}) and
#' works well in practice when \eqn{\nu > 10}. Because it reduces to a single
#' call to \code{pt()}, it is orders of magnitude faster than the numerical
#' integration method (\code{\link{ptdiff_NI}}) and fully vectorised.
#'
#' @examples
#' # P(T_t - T_c > 3) with equal parameters
#' ptdiff_MM(q = 3, mu_t = 2, mu_c = 0, sd_t = 1, sd_c = 1,
#'           nu_t = 17, nu_c = 17, lower.tail = FALSE)
#'
#' # P(T_t - T_c > 1) with unequal variances
#' ptdiff_MM(q = 1, mu_t = 5, mu_c = 3, sd_t = 2, sd_c = 1.5,
#'           nu_t = 10, nu_c = 15, lower.tail = FALSE)
#'
#' # P(T_t - T_c > 0) with different degrees of freedom
#' ptdiff_MM(q = 0, mu_t = 1, mu_c = 1, sd_t = 1, sd_c = 1,
#'           nu_t = 5, nu_c = 20, lower.tail = FALSE)
#'
#' # Lower tail: P(T_t - T_c <= 2)
#' ptdiff_MM(q = 2, mu_t = 3, mu_c = 0, sd_t = 1.5, sd_c = 1.2,
#'           nu_t = 12, nu_c = 15, lower.tail = TRUE)
#'
#' # Vectorised usage
#' ptdiff_MM(q = 1, mu_t = c(2, 3, 4), mu_c = c(0, 1, 2),
#'           sd_t = c(1, 1.2, 1.5), sd_c = c(1, 1.1, 1.3),
#'           nu_t = 10, nu_c = 10, lower.tail = FALSE)
#'
#' @importFrom stats pt
#' @export
ptdiff_MM <- function(q, mu_t, mu_c, sd_t, sd_c, nu_t, nu_c,
                      lower.tail = TRUE) {

  # --- Input validation ---
  if (!is.numeric(q) || length(q) != 1L || is.na(q)) {
    stop("'q' must be a single numeric value")
  }

  for (nm in c("mu_t", "mu_c", "sd_t", "sd_c")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) < 1L || any(is.na(val))) {
      stop(paste0("'", nm, "' must be a numeric scalar or vector with no missing values"))
    }
  }

  if (any(sd_t <= 0)) stop("'sd_t' must contain only positive values")
  if (any(sd_c <= 0)) stop("'sd_c' must contain only positive values")

  if (!is.numeric(nu_t) || length(nu_t) != 1L || is.na(nu_t) || nu_t <= 4) {
    stop("'nu_t' must be a single numeric value greater than 4")
  }

  if (!is.numeric(nu_c) || length(nu_c) != 1L || is.na(nu_c) || nu_c <= 4) {
    stop("'nu_c' must be a single numeric value greater than 4")
  }

  if (!is.logical(lower.tail) || length(lower.tail) != 1L || is.na(lower.tail)) {
    stop("'lower.tail' must be a single logical value (TRUE or FALSE)")
  }

  # --- Moment-Matching approximation ---
  # Step 1: Mean of the difference
  mu.diff <- mu_t - mu_c

  # Step 2: Effective degrees of freedom (nu.star) from fourth-moment matching
  var_t   <- sd_t ^ 2 * nu_t / (nu_t - 2)
  var_c   <- sd_c ^ 2 * nu_c / (nu_c - 2)
  Qu_star <- (var_t + var_c) ^ 2
  term1   <- (sd_t ^ 2) ^ 2 * nu_t ^ 2 / ((nu_t - 2) * (nu_t - 4))
  term2   <- (sd_c ^ 2) ^ 2 * nu_c ^ 2 / ((nu_c - 2) * (nu_c - 4))
  term3   <- 2 * var_t * var_c
  Qu      <- term1 + term2 + term3
  nu.star <- (2 * Qu_star - 4 * Qu) / (Qu_star - Qu)

  # Step 3: Scale parameter (sigma.star) from second-moment matching
  sigma.star <- sqrt((var_t + var_c) * (nu.star - 2) / nu.star)

  # Step 4: Standardised statistic and CDF evaluation
  standardized <- (q - mu.diff) / sigma.star
  results      <- pt(standardized, df = nu.star, lower.tail = FALSE)

  # Convert to lower tail if requested:
  #   lower.tail = TRUE  -> P(T_t - T_c <= q) = 1 - P(T_t - T_c > q)
  #   lower.tail = FALSE -> P(T_t - T_c > q)
  results <- lower.tail + c(1, -1)[lower.tail + 1] * results

  return(results)
}
