#' Calculate Appell's First Hypergeometric Function of Two Variables
#'
#' This function computes Appell's first hypergeometric function F1(a; b1, b2; c; x, y)
#' using numerical integration of its integral representation. The Appell F1 function
#' is a generalization of the Gauss hypergeometric function to two variables and
#' appears in various statistical applications including beta distribution differences.
#'
#' @param a A numeric value representing the first parameter (must satisfy Re(a) > 0
#'        and Re(c - a) > 0 for convergence).
#' @param b1 A numeric value representing the second parameter.
#' @param b2 A numeric value representing the third parameter.
#' @param c A numeric value representing the fourth parameter (must be greater than a).
#' @param x A numeric value representing the first variable (typically |x| < 1 for convergence).
#' @param y A numeric value representing the second variable (typically |y| < 1 for convergence).
#'
#' @return A numeric value of Appell's first hypergeometric function F1(a; b1, b2; c; x, y).
#'
#' @details
#' The function evaluates Appell's F1 hypergeometric function using the integral representation:
#' \deqn{F_1(a; b_1, b_2; c; x, y) = \frac{\Gamma(c)}{\Gamma(a)\Gamma(c-a)} \int_0^1 t^{a-1}(1-t)^{c-a-1}(1-xt)^{-b_1}(1-yt)^{-b_2} dt}
#'
#' where \eqn{\Gamma(\cdot)} is the gamma function. This implementation uses numerical
#' integration via the \code{integrate} function from the stats package.
#'
#' The Appell F1 function satisfies the following properties:
#' \itemize{
#'   \item Symmetry: \eqn{F_1(a; b_1, b_2; c; x, y) = F_1(a; b_2, b_1; c; y, x)}
#'   \item Reduction to Gauss hypergeometric: \eqn{F_1(a; b, 0; c; x, y) = {}_2F_1(a, b; c; x)}
#'   \item Special case: \eqn{F_1(a; 0, 0; c; x, y) = 1}
#' }
#'
#' @examples
#' # Calculate F1(0.5; 0.5, 0; 1; 0.96, 1.2) - reduces to hypergeometric function
#' AppellsF1(0.5, 0.5, 0, 1, 0.96, 1.2)
#'
#' # Calculate F1(1; 1, 1; 2; 0.3, 0.4)
#' AppellsF1(1, 1, 1, 2, 0.3, 0.4)
#'
#' # Verify symmetry property: F1(a; b1, b2; c; x, y) = F1(a; b2, b1; c; y, x)
#' result1 <- AppellsF1(1, 2, 3, 4, 0.2, 0.3)
#' result2 <- AppellsF1(1, 3, 2, 4, 0.3, 0.2)
#' all.equal(result1, result2)
#'
#' @importFrom stats integrate
#' @export
AppellsF1 <- function(a, b1, b2, c, x, y) {
  # Calculate the normalization coefficient using gamma functions
  coefficient <- gamma(c) / (gamma(a) * gamma(c - a))

  # Define the integrand function for numerical integration
  integrand <- function(t) {
    t ^ (a - 1) * (1 - t) ^ (c - a - 1) * (1 - x * t) ^ (-b1) * (1 - y * t) ^ (-b2)
  }

  # Perform numerical integration from 0 to 1 using adaptive quadrature
  integral_result <- integrate(integrand, lower = 0, upper = 1)[['value']]

  # Return the final result as coefficient times integral
  return(coefficient * integral_result)
}
