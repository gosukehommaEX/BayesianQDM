#' Cumulative Distribution Function of the Difference of Two t-Distributed Variables
#' by Monte Carlo Simulation
#'
#' Calculates the cumulative distribution function (CDF) of the difference between
#' two independent non-standardised t-distributed random variables using Monte Carlo
#' simulation. Specifically, computes \eqn{P(T_1 - T_2 \le q)} or
#' \eqn{P(T_1 - T_2 > q)}, where \eqn{T_i \sim t(\mu_i, \sigma_i^2, \nu_i)}.
#'
#' @param nMC A positive integer giving the number of Monte Carlo draws. Typical
#'        values range from 10,000 (quick estimates) to 100,000 or more (high
#'        precision). Larger values reduce Monte Carlo error at the cost of
#'        computation time.
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
#' The algorithm proceeds as follows:
#' \enumerate{
#'   \item Generate an \code{nMC x n} matrix of standard t draws for \eqn{T_1}
#'         (\eqn{\nu_1} degrees of freedom), then scale and shift each column
#'         by \code{sd.t1[i]} and \code{mu.t1[i]}.
#'   \item Repeat for \eqn{T_2} with \eqn{\nu_2} degrees of freedom.
#'   \item Compute the \code{nMC x n} difference matrix \eqn{D = T_1 - T_2}.
#'   \item Return \code{colMeans(D > q)} as the estimated \eqn{P(T_1 - T_2 > q)}
#'         for each parameter set.
#' }
#'
#' All operations are matrix-based (no R-level loop over parameter sets), so
#' performance scales well with \eqn{n}. However, note that the matrix size is
#' \code{nMC x n}, so memory usage grows linearly with both \code{nMC} and
#' \eqn{n}. When \eqn{n} is large (e.g., \code{nsim x n_scenarios} in
#' \code{\link{pGNGsinglecontinuous}}), memory requirements can become prohibitive;
#' in such cases prefer \code{CalcMethod = 'MM'}.
#'
#' Monte Carlo error is approximately \eqn{\sqrt{p(1-p)/\mathrm{nMC}}};
#' near \eqn{p = 0.5} this is roughly \eqn{0.5/\sqrt{\mathrm{nMC}}}.
#'
#' @examples
#' # P(T1 - T2 > 3) with equal parameters
#' ptdiff_MC(nMC = 1e5, q = 3, mu.t1 = 2, mu.t2 = 0, sd.t1 = 1, sd.t2 = 1,
#'           nu.t1 = 17, nu.t2 = 17, lower.tail = FALSE)
#'
#' # P(T1 - T2 > 1) with unequal variances
#' ptdiff_MC(nMC = 1e5, q = 1, mu.t1 = 5, mu.t2 = 3, sd.t1 = 2, sd.t2 = 1.5,
#'           nu.t1 = 10, nu.t2 = 15, lower.tail = FALSE)
#'
#' # P(T1 - T2 > 0) with different degrees of freedom
#' ptdiff_MC(nMC = 1e5, q = 0, mu.t1 = 1, mu.t2 = 1, sd.t1 = 1, sd.t2 = 1,
#'           nu.t1 = 5, nu.t2 = 20, lower.tail = FALSE)
#'
#' # Lower tail: P(T1 - T2 <= 2)
#' ptdiff_MC(nMC = 1e5, q = 2, mu.t1 = 3, mu.t2 = 0, sd.t1 = 1.5, sd.t2 = 1.2,
#'           nu.t1 = 12, nu.t2 = 15, lower.tail = TRUE)
#'
#' # Vectorised usage
#' ptdiff_MC(nMC = 1e5, q = 1, mu.t1 = c(2, 3, 4), mu.t2 = c(0, 1, 2),
#'           sd.t1 = c(1, 1.2, 1.5), sd.t2 = c(1, 1.1, 1.3),
#'           nu.t1 = 10, nu.t2 = 10, lower.tail = FALSE)
#'
#' @importFrom stats rt
#' @export
ptdiff_MC <- function(nMC, q, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2,
                      lower.tail = TRUE) {

  # --- Input validation ---
  if (!is.numeric(nMC) || length(nMC) != 1L || is.na(nMC) ||
      nMC != floor(nMC) || nMC < 1L) {
    stop("'nMC' must be a single positive integer")
  }

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

  # --- Monte Carlo estimation ---
  # Number of parameter sets
  n <- max(length(mu.t1), length(mu.t2), length(sd.t1), length(sd.t2))

  # Recycle location/scale parameters to length n
  mu.t1 <- rep(mu.t1, length.out = n)
  mu.t2 <- rep(mu.t2, length.out = n)
  sd.t1 <- rep(sd.t1, length.out = n)
  sd.t2 <- rep(sd.t2, length.out = n)

  # Generate nMC x n matrices of scaled and shifted t draws.
  # rep(..., each = nMC) broadcasts column-wise so that column j uses
  # scale sd.t1[j] and location mu.t1[j].
  rand.t1 <- matrix(rt(nMC * n, df = nu.t1), nrow = nMC, ncol = n) *
    rep(sd.t1, each = nMC) + rep(mu.t1, each = nMC)

  rand.t2 <- matrix(rt(nMC * n, df = nu.t2), nrow = nMC, ncol = n) *
    rep(sd.t2, each = nMC) + rep(mu.t2, each = nMC)

  # Column-wise proportion of differences exceeding q: P(T1 - T2 > q)
  results <- colMeans(rand.t1 - rand.t2 > q)

  # Convert to lower tail if requested:
  #   lower.tail = TRUE  -> P(T1 - T2 <= q) = 1 - P(T1 - T2 > q)
  #   lower.tail = FALSE -> P(T1 - T2 > q)
  results <- lower.tail + c(1, -1)[lower.tail + 1] * results

  return(results)
}
