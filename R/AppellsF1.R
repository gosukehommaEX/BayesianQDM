#' Calculate Appell's First Hypergeometric Function of Two Variables
#'
#' This function computes the numeric value of Appell's first hypergeometric function
#' F1(a; b1, b2; c; x, y) using integral representation. The Appell hypergeometric
#' function is a generalization of the hypergeometric function to two variables.
#'
#' @param a A numeric value representing the first parameter.
#' @param b1 A numeric value representing the second parameter.
#' @param b2 A numeric value representing the third parameter.
#' @param c A numeric value representing the fourth parameter.
#' @param x A numeric value representing the first variable.
#' @param y A numeric value representing the second variable.
#'
#' @return A numeric value of Appell's first hypergeometric function F1(a; b1, b2; c; x, y).
#'
#' @details
#' The function uses the integral representation:
#' \deqn{F_1(a; b_1, b_2; c; x, y) = \frac{\Gamma(c)}{\Gamma(a)\Gamma(c-a)} \int_0^1 t^{a-1}(1-t)^{c-a-1}(1-xt)^{-b_1}(1-yt)^{-b_2} dt}
#'
#' @examples
#' # Calculate F1(0.5; 0.5, 0; 1; 0.96, 1.2)
#' AppellsF1(0.5, 0.5, 0, 1, 0.96, 1.2)
#'
#' # Calculate F1(1; 1, 1; 2; 0.3, 0.4)
#' AppellsF1(1, 1, 1, 2, 0.3, 0.4)
#'
#' @importFrom stats integrate
#' @export
AppellsF1 <- function(a, b1, b2, c, x, y) {
  # Calculate the coefficient using gamma functions
  coefficient <- gamma(c) / (gamma(a) * gamma(c - a))

  # Define the integrand function
  integrand <- function(t) {
    t ^ (a - 1) * (1 - t) ^ (c - a - 1) * (1 - x * t) ^ (-b1) * (1 - y * t) ^ (-b2)
  }

  # Perform numerical integration from 0 to 1
  integral_result <- integrate(integrand, lower = 0, upper = 1)[['value']]

  # Return the final result
  return(coefficient * integral_result)
}
