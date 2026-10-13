#' Joint cell probabilities of two binary endpoints
#'
#' Converts the marginal response rates of two binary endpoints and their
#' correlation into the probabilities of the four joint response patterns.
#'
#' For two Bernoulli variables with response rates \eqn{\pi_1} and
#' \eqn{\pi_2}, the Pearson correlation \eqn{\rho} determines
#' \deqn{p_{11} = \pi_1 \pi_2 + \rho \sqrt{\pi_1 (1 - \pi_1) \pi_2 (1 - \pi_2)},}
#' and \eqn{p_{10} = \pi_1 - p_{11}}, \eqn{p_{01} = \pi_2 - p_{11}},
#' \eqn{p_{00} = 1 - p_{10} - p_{01} - p_{11}}. All four probabilities are
#' non-negative only if
#' \eqn{\max(0, \pi_1 + \pi_2 - 1) \le p_{11} \le \min(\pi_1, \pi_2)}, which
#' bounds the feasible range of \eqn{\rho}. An error is raised if \code{rho}
#' is outside this range by more than \code{tol}.
#'
#' @param pi1 Response rate of endpoint 1, strictly between 0 and 1.
#' @param pi2 Response rate of endpoint 2, strictly between 0 and 1.
#' @param rho Correlation between the two endpoints.
#' @param tol Tolerance for the feasibility check.
#'
#' @return A named numeric vector of the probabilities of the cells
#'   \code{"00"}, \code{"01"}, \code{"10"} and \code{"11"} (first digit:
#'   endpoint 1; second digit: endpoint 2), summing to 1.
#'
#' @family distributions
#' @examples
#' cellprob_binary2(pi1 = 0.3, pi2 = 0.4, rho = 0.3)
#' @export
cellprob_binary2 <- function(pi1, pi2, rho, tol = 1e-10) {
  check_numeric(pi1, "pi1", lower = 0, upper = 1, lower_open = TRUE,
                upper_open = TRUE)
  check_numeric(pi2, "pi2", lower = 0, upper = 1, lower_open = TRUE,
                upper_open = TRUE)
  check_numeric(rho, "rho")
  check_numeric(tol, "tol", lower = 0)
  sd_prod <- sqrt(pi1 * (1 - pi1) * pi2 * (1 - pi2))
  p11_min <- max(0, pi1 + pi2 - 1)
  p11_max <- min(pi1, pi2)
  rho_min <- (p11_min - pi1 * pi2) / sd_prod
  rho_max <- (p11_max - pi1 * pi2) / sd_prod
  if (rho < rho_min - tol || rho > rho_max + tol) {
    stop(sprintf(paste("'rho' = %s is infeasible for pi1 = %s and pi2 = %s;",
                       "the feasible range is from %s to %s."),
                 format(rho, digits = 6), format(pi1, digits = 6),
                 format(pi2, digits = 6), format(rho_min, digits = 6),
                 format(rho_max, digits = 6)),
         call. = FALSE)
  }
  # Clamp to the feasible range to absorb rounding errors
  p11 <- min(max(pi1 * pi2 + rho * sd_prod, p11_min), p11_max)
  p10 <- pi1 - p11
  p01 <- pi2 - p11
  p00 <- 1 - p10 - p01 - p11
  p <- pmax(c("00" = p00, "01" = p01, "10" = p10, "11" = p11), 0)
  p / sum(p)
}
