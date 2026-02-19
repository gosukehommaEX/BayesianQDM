#' Compute Joint Bivariate Binary Probabilities from Marginal Rates and Correlation
#'
#' Converts the parameterization \eqn{(\pi_{k1}, \pi_{k2}, \rho_k)} into the
#' four-cell joint probability vector \eqn{(p_{k,00}, p_{k,01}, p_{k,10}, p_{k,11})}
#' for a bivariate binary outcome in arm \eqn{k}.  The conversion uses the
#' standard moment-matching identity for the Pearson correlation of two Bernoulli
#' variables, and checks that the requested correlation is feasible for the
#' supplied marginal rates.
#'
#' @param pi1 A single numeric value in \eqn{(0, 1)} giving the marginal
#'        response probability for Endpoint 1.
#' @param pi2 A single numeric value in \eqn{(0, 1)} giving the marginal
#'        response probability for Endpoint 2.
#' @param rho A single numeric value giving the Pearson correlation between
#'        Endpoint 1 and Endpoint 2.  Must lie within the feasible range
#'        \eqn{[\rho_{\min}, \rho_{\max}]} determined by \code{pi1} and
#'        \code{pi2} (see Details).  Use \code{rho = 0} for independence.
#' @param tol A single positive numeric value specifying the tolerance used
#'        when checking whether \code{rho} lies within its feasible range and
#'        whether the resulting probabilities are non-negative.  Default is
#'        \code{1e-10}.
#'
#' @return A named numeric vector of length 4 with elements
#'         \code{p00}, \code{p01}, \code{p10}, \code{p11}, where
#'         \code{p_lm = Pr(Endpoint 1 = l, Endpoint 2 = m)}.
#'         All elements are non-negative and sum to 1.
#'
#' @details
#' For a bivariate binary outcome \eqn{(Y_1, Y_2)} with marginal success
#' probabilities \eqn{\pi_1 = \Pr(Y_1 = 1)} and \eqn{\pi_2 = \Pr(Y_2 = 1)},
#' the Pearson correlation is
#' \deqn{
#'   \rho = \frac{p_{11} - \pi_1 \pi_2}
#'              {\sqrt{\pi_1 (1-\pi_1) \pi_2 (1-\pi_2)}}.
#' }
#' Solving for \eqn{p_{11}} gives
#' \deqn{p_{11} = \rho \sqrt{\pi_1(1-\pi_1)\pi_2(1-\pi_2)} + \pi_1 \pi_2,}
#' from which the remaining probabilities follow:
#' \deqn{p_{10} = \pi_1 - p_{11}, \quad
#'       p_{01} = \pi_2 - p_{11}, \quad
#'       p_{00} = 1 - p_{10} - p_{01} - p_{11}.}
#'
#' For \eqn{p_{11}} to keep all four cell probabilities in \eqn{[0, 1]},
#' the correlation must satisfy
#' \deqn{
#'   \rho_{\min} = \frac{\max(0,\, \pi_1 + \pi_2 - 1) - \pi_1 \pi_2}
#'                     {\sqrt{\pi_1(1-\pi_1)\pi_2(1-\pi_2)}}
#'   \le \rho \le
#'   \frac{\min(\pi_1, \pi_2) - \pi_1 \pi_2}
#'        {\sqrt{\pi_1(1-\pi_1)\pi_2(1-\pi_2)}}
#'   = \rho_{\max}.
#' }
#' The function raises an error if \code{rho} falls outside this range
#' (subject to \code{tol}).
#'
#' @examples
#' # Example 1: Independent endpoints (rho = 0)
#' getjointbin(pi1 = 0.3, pi2 = 0.4, rho = 0.0)
#'
#' # Example 2: Positive correlation
#' getjointbin(pi1 = 0.3, pi2 = 0.4, rho = 0.3)
#'
#' # Example 3: Negative correlation
#' getjointbin(pi1 = 0.3, pi2 = 0.4, rho = -0.2)
#'
#' # Example 4: Verify cell probabilities sum to 1
#' p <- getjointbin(pi1 = 0.25, pi2 = 0.35, rho = 0.1)
#' sum(p)  # Should be 1
#'
#' # Example 5: Verify marginal recovery
#' p <- getjointbin(pi1 = 0.25, pi2 = 0.35, rho = 0.1)
#' p["p10"] + p["p11"]  # Should equal 0.25 (pi1)
#' p["p01"] + p["p11"]  # Should equal 0.35 (pi2)
#'
#' @export
getjointbin <- function(pi1, pi2, rho, tol = 1e-10) {

  # --- Input validation ---
  if (!is.numeric(pi1) || length(pi1) != 1L || is.na(pi1) ||
      pi1 <= 0 || pi1 >= 1) {
    stop("'pi1' must be a single numeric value in (0, 1)")
  }

  if (!is.numeric(pi2) || length(pi2) != 1L || is.na(pi2) ||
      pi2 <= 0 || pi2 >= 1) {
    stop("'pi2' must be a single numeric value in (0, 1)")
  }

  if (!is.numeric(rho) || length(rho) != 1L || is.na(rho)) {
    stop("'rho' must be a single numeric value")
  }

  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("'tol' must be a single positive numeric value")
  }

  # --- Compute the standard deviation product (denominator of rho formula) ---
  # denom = sqrt(pi1 * (1 - pi1) * pi2 * (1 - pi2))
  denom <- sqrt(pi1 * (1 - pi1) * pi2 * (1 - pi2))

  # --- Compute feasible range of rho given (pi1, pi2) ---
  # Lower bound: p11 >= max(0, pi1 + pi2 - 1)  =>  rho >= rho_min
  # Upper bound: p11 <= min(pi1, pi2)           =>  rho <= rho_max
  rho_min <- (max(0, pi1 + pi2 - 1) - pi1 * pi2) / denom
  rho_max <- (min(pi1, pi2)         - pi1 * pi2) / denom

  # --- Check feasibility of rho ---
  if (rho < rho_min - tol || rho > rho_max + tol) {
    stop(sprintf(
      paste0("'rho' = %.6f is infeasible for pi1 = %.4f, pi2 = %.4f.\n",
             "  Feasible range: [%.6f, %.6f]"),
      rho, pi1, pi2, rho_min, rho_max
    ))
  }

  # Clamp rho to the feasible range to avoid tiny numerical violations
  rho <- min(max(rho, rho_min), rho_max)

  # --- Recover the four joint probabilities ---
  # p11 is derived from the rho formula rearranged for p11
  p11 <- rho * denom + pi1 * pi2

  # Clamp p11 to its valid interval [max(0, pi1+pi2-1), min(pi1, pi2)]
  # to absorb any residual floating-point error
  p11 <- min(max(p11, max(0, pi1 + pi2 - 1)), min(pi1, pi2))

  p10 <- pi1 - p11
  p01 <- pi2 - p11
  p00 <- 1   - p10 - p01 - p11

  # Assemble the probability vector (order: 00, 01, 10, 11)
  p <- c(p00 = p00, p01 = p01, p10 = p10, p11 = p11)

  # --- Sanity checks on the resulting probabilities ---
  if (any(p < -tol)) {
    stop(sprintf(
      paste0("Computed joint probabilities contain a negative value ",
             "(min = %.2e). This should not occur; please check inputs."),
      min(p)
    ))
  }

  # Replace tiny negative values caused by floating-point arithmetic with 0
  p[p < 0] <- 0

  # Re-normalise to ensure exact summation to 1
  p <- p / sum(p)

  return(p)
}
