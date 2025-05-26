#' The Appell's first hypergeometric function of two variables x and y.
#'
#' This function calculates the numeric value of the Appell's first hypergeometric function.
#'
#' @param a  A numeric value.
#' @param b1 A numeric value.
#' @param b2 A numeric value.
#' @param c A numeric value.
#' @param x A numeric value.
#' @param y A numeric value.
#'
#' @return The \code{AppellsF1} gives the numeric value of the Appell's first hypergeometric function.
#'
#' @examples
#' AppellsF1(0.5, 0.5, 0, 1, 0.96, 1.2)
#'
#' @importFrom stats integrate
#' @export
AppellsF1 = function(a, b1, b2, c, x, y) {
  '*'(
    gamma(c) / (gamma(a) * gamma(c - a)),
    integrate(
      function(t) {
        t ^ (a - 1) * (1 - t) ^ (c - a - 1) * (1 - x * t) ^ (-b1) * (1 - y * t) ^ (-b2)
      },
      lower = 0,
      upper = 1
    )[['value']]
  )
}
