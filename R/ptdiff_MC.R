#' Cumulative Distribution Function of the Difference of Two
#' Independent t-Distributed Variables via Monte Carlo Simulation
#'
#' Calculates the cumulative distribution function (CDF) of the difference
#' between two independent non-standardised t-distributed random variables
#' using Monte Carlo simulation. Specifically, computes
#' \eqn{P(T_t - T_c \le q)} or \eqn{P(T_t - T_c > q)}, where
#' \eqn{T_k \sim t(\mu_k, \sigma_k^2, \nu_k)} for \eqn{k \in \{t, c\}}.
#'
#' @param nMC A positive integer giving the number of Monte Carlo draws.
#'        Typical values range from 10,000 (quick estimates) to 100,000
#'        or more (high precision). Larger values reduce Monte Carlo error
#'        at the cost of computation time.
#' @param q A numeric scalar representing the quantile threshold.
#' @param mu_t A numeric scalar or vector giving the location parameter of
#'        the t-distribution for the treatment group.
#' @param mu_c A numeric scalar or vector giving the location parameter of
#'        the t-distribution for the control group.
#' @param sd_t A positive numeric scalar or vector giving the scale parameter
#'        of the t-distribution for the treatment group.
#' @param sd_c A positive numeric scalar or vector giving the scale parameter
#'        of the t-distribution for the control group.
#' @param nu_t A numeric scalar giving the degrees of freedom of the
#'        t-distribution for the treatment group.
#'        Must be greater than 2 for finite variance.
#' @param nu_c A numeric scalar giving the degrees of freedom of the
#'        t-distribution for the control group.
#'        Must be greater than 2 for finite variance.
#' @param lower.tail A logical scalar; if \code{TRUE} (default), the function
#'        returns \eqn{P(T_t - T_c \le q)}, otherwise
#'        \eqn{P(T_t - T_c > q)}.
#'
#' @return A numeric scalar or vector in \code{[0, 1]}.  When \code{mu_t},
#'         \code{mu_c}, \code{sd_t}, or \code{sd_c} are vectors of length
#'         \eqn{n}, a vector of length \eqn{n} is returned.
#'
#' @details
#' The algorithm proceeds as follows:
#' \enumerate{
#'   \item Generate an \code{nMC x n} matrix of standard t draws for
#'         \eqn{T_t} (\eqn{\nu_t} degrees of freedom), then scale and shift
#'         each column by \code{sd_t[i]} and \code{mu_t[i]}.
#'   \item Repeat for \eqn{T_c} with \eqn{\nu_c} degrees of freedom.
#'   \item Compute the \code{nMC x n} difference matrix
#'         \eqn{D = T_t - T_c}.
#'   \item Return \code{colMeans(D > q)} as the estimated
#'         \eqn{P(T_t - T_c > q)} for each parameter set.
#' }
#'
#' All operations are matrix-based (no R-level loop over parameter sets), so
#' performance scales well with \eqn{n}. However, note that the matrix size
#' is \code{nMC x n}, so memory usage grows linearly with both \code{nMC}
#' and \eqn{n}. When \eqn{n} is large (e.g., \code{nsim x n_scenarios} in
#' \code{\link{pbayesdecisionprob1cont}}), memory requirements can become
#' prohibitive; in such cases prefer \code{CalcMethod = 'MM'}.
#'
#' Monte Carlo error is approximately
#' \eqn{\sqrt{p(1 - p) / \mathrm{nMC}}}; near \eqn{p = 0.5} this is
#' roughly \eqn{0.5 / \sqrt{\mathrm{nMC}}}.
#'
#' @examples
#' # P(T_t - T_c > 3) with equal parameters
#' ptdiff_MC(nMC = 1e5, q = 3, mu_t = 2, mu_c = 0, sd_t = 1, sd_c = 1,
#'           nu_t = 17, nu_c = 17, lower.tail = FALSE)
#'
#' # P(T_t - T_c > 1) with unequal scales
#' ptdiff_MC(nMC = 1e5, q = 1, mu_t = 5, mu_c = 3, sd_t = 2, sd_c = 1.5,
#'           nu_t = 10, nu_c = 15, lower.tail = FALSE)
#'
#' # P(T_t - T_c > 0) with different degrees of freedom
#' ptdiff_MC(nMC = 1e5, q = 0, mu_t = 1, mu_c = 1, sd_t = 1, sd_c = 1,
#'           nu_t = 5, nu_c = 20, lower.tail = FALSE)
#'
#' # Lower tail: P(T_t - T_c <= 2)
#' ptdiff_MC(nMC = 1e5, q = 2, mu_t = 3, mu_c = 0, sd_t = 1.5, sd_c = 1.2,
#'           nu_t = 12, nu_c = 15, lower.tail = TRUE)
#'
#' # Vectorised usage
#' ptdiff_MC(nMC = 1e5, q = 1, mu_t = c(2, 3, 4), mu_c = c(0, 1, 2),
#'           sd_t = c(1, 1.2, 1.5), sd_c = c(1, 1.1, 1.3),
#'           nu_t = 10, nu_c = 10, lower.tail = FALSE)
#'
#' @importFrom stats rt
#' @export
ptdiff_MC <- function(nMC, q, mu_t, mu_c, sd_t, sd_c, nu_t, nu_c,
                      lower.tail = TRUE) {

  # --- Input validation ---
  if (!is.numeric(nMC) || length(nMC) != 1L || is.na(nMC) ||
      nMC != floor(nMC) || nMC < 1L) {
    stop("'nMC' must be a single positive integer")
  }

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

  # --- Monte Carlo estimation ---
  # Number of parameter sets
  n <- max(length(mu_t), length(mu_c), length(sd_t), length(sd_c))

  # Recycle location/scale parameters to length n
  mu_t <- rep(mu_t, length.out = n)
  mu_c <- rep(mu_c, length.out = n)
  sd_t <- rep(sd_t, length.out = n)
  sd_c <- rep(sd_c, length.out = n)

  # Generate nMC x n matrices of scaled and shifted t draws.
  # rep(..., each = nMC) broadcasts column-wise so that column j uses
  # scale sd_t[j] and location mu_t[j].
  rand.t <- matrix(rt(nMC * n, df = nu_t), nrow = nMC, ncol = n) *
    rep(sd_t, each = nMC) + rep(mu_t, each = nMC)

  rand.c <- matrix(rt(nMC * n, df = nu_c), nrow = nMC, ncol = n) *
    rep(sd_c, each = nMC) + rep(mu_c, each = nMC)

  # Column-wise proportion of differences exceeding q: P(T_t - T_c > q)
  results <- colMeans(rand.t - rand.c > q)

  # Convert to lower tail if requested:
  #   lower.tail = TRUE  -> P(T_t - T_c <= q) = 1 - P(T_t - T_c > q)
  #   lower.tail = FALSE -> P(T_t - T_c > q)
  results <- lower.tail + c(1, -1)[lower.tail + 1] * results

  return(results)
}
