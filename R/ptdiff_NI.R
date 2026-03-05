#' Cumulative Distribution Function of the Difference of Two t-Distributed Variables
#' by Numerical Integration
#'
#' Calculates the cumulative distribution function (CDF) of the difference between
#' two independent non-standardised t-distributed random variables using exact
#' numerical integration via convolution. Specifically, computes
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
#'        t-distribution (treatment). Must be greater than 2 for finite variance.
#' @param nu_c A numeric scalar giving the degrees of freedom of the second
#'        t-distribution (control). Must be greater than 2 for finite variance.
#' @param lower.tail A logical scalar; if \code{TRUE} (default), the function
#'        returns \eqn{P(T_t - T_c \le q)}, otherwise \eqn{P(T_t - T_c > q)}.
#'
#' @return A numeric scalar or vector in \code{[0, 1]}.  When \code{mu_t},
#'         \code{mu_c}, \code{sd_t}, or \code{sd_c} are vectors of length
#'         \eqn{n}, a vector of length \eqn{n} is returned.
#'
#' @details
#' The upper-tail probability is obtained via the convolution formula:
#' \deqn{P(T_t - T_c > q)
#'   = \int_{-\infty}^{\infty} f_t(x)\, F_c(x - q)\, dx}
#' where \eqn{f_t} is the density of \eqn{T_t} and \eqn{F_c} is the CDF of
#' \eqn{T_c}. Integration bounds are set to \eqn{\mu_t \pm 8\sigma^*_t}, where
#' \eqn{\sigma^*_t = \sigma_t \sqrt{(\nu_t + 1)/(\nu_t - 2)}} is the standard
#' deviation of \eqn{T_t}. Relative and absolute tolerances are \code{1e-6} and
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
#' Moment-Matching approximation (\code{\link{ptdiff_MM}}) because it calls
#' \code{integrate()} once per parameter set. For large-scale simulation
#' (many parameter sets), prefer \code{CalcMethod = 'MM'} in
#' \code{\link{pbayesdecisionprob1cont}}.
#'
#' @examples
#' # P(T_t - T_c > 3) with equal parameters
#' ptdiff_NI(q = 3, mu_t = 2, mu_c = 0, sd_t = 1, sd_c = 1,
#'           nu_t = 17, nu_c = 17, lower.tail = FALSE)
#'
#' # P(T_t - T_c > 1) with unequal variances
#' ptdiff_NI(q = 1, mu_t = 5, mu_c = 3, sd_t = 2, sd_c = 1.5,
#'           nu_t = 10, nu_c = 15, lower.tail = FALSE)
#'
#' # P(T_t - T_c > 0) with different degrees of freedom
#' ptdiff_NI(q = 0, mu_t = 1, mu_c = 1, sd_t = 1, sd_c = 1,
#'           nu_t = 5, nu_c = 20, lower.tail = FALSE)
#'
#' # Lower tail: P(T_t - T_c <= 2)
#' ptdiff_NI(q = 2, mu_t = 3, mu_c = 0, sd_t = 1.5, sd_c = 1.2,
#'           nu_t = 12, nu_c = 15, lower.tail = TRUE)
#'
#' # Vectorised usage
#' ptdiff_NI(q = 1, mu_t = c(2, 3, 4), mu_c = c(0, 1, 2),
#'           sd_t = c(1, 1.2, 1.5), sd_c = c(1, 1.1, 1.3),
#'           nu_t = 10, nu_c = 10, lower.tail = FALSE)
#'
#' @importFrom stats dt integrate pt
#' @export
ptdiff_NI <- function(q, mu_t, mu_c, sd_t, sd_c, nu_t, nu_c,
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

  if (!is.numeric(nu_t) || length(nu_t) != 1L || is.na(nu_t) || nu_t <= 2) {
    stop("'nu_t' must be a single numeric value greater than 2")
  }

  if (!is.numeric(nu_c) || length(nu_c) != 1L || is.na(nu_c) || nu_c <= 2) {
    stop("'nu_c' must be a single numeric value greater than 2")
  }

  if (!is.logical(lower.tail) || length(lower.tail) != 1L || is.na(lower.tail)) {
    stop("'lower.tail' must be a single logical value (TRUE or FALSE)")
  }

  # --- Core scalar integration function ---
  # Computes P(T_t - T_c > q) for a single parameter set via convolution.
  integrand_upper <- function(x, q, mu_t, mu_c, sd_t, sd_c, nu_t, nu_c) {
    # f_t(x) * F_c(x - q): density of T_t times CDF of T_c evaluated at (x - q)
    dt((x - mu_t) / sd_t, df = nu_t) / sd_t *
      pt((x - q - mu_c) / sd_c, df = nu_c)
  }

  # Integration bounds: mu_t +/- 8 * sd*(T_t), where sd*(T_t) = sd_t * sqrt((nu_t+1)/(nu_t-2))
  scalar_upper <- function(q, mu_t, mu_c, sd_t, sd_c, nu_t, nu_c) {
    spread <- sd_t * sqrt((nu_t + 1) / (nu_t - 2))
    lo     <- mu_t - 8 * spread
    hi     <- mu_t + 8 * spread
    integrate(
      integrand_upper,
      lower     = lo,
      upper     = hi,
      q         = q,
      mu_t      = mu_t,
      mu_c      = mu_c,
      sd_t      = sd_t,
      sd_c      = sd_c,
      nu_t      = nu_t,
      nu_c      = nu_c,
      rel.tol   = 1e-6,
      abs.tol   = 1e-8
    )$value
  }

  # --- Vectorised computation ---
  # Use mapply to apply the scalar function across all parameter combinations
  p_upper <- mapply(
    scalar_upper,
    q    = q,
    mu_t = mu_t,
    mu_c = mu_c,
    sd_t = sd_t,
    sd_c = sd_c,
    MoreArgs = list(nu_t = nu_t, nu_c = nu_c)
  )

  # Clamp to [0, 1] to guard against minor numerical overshoot
  p_upper <- pmin(pmax(p_upper, 0), 1)

  if (lower.tail) 1 - p_upper else p_upper
}
