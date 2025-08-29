#' Density of difference between two beta variables
#'
#' This function calculates density of difference between two beta variables
#' using Appell’s first hypergeometric function of two variables.
#'
#' @param pi A numeric value representing the difference of the two beta variables.
#' @param alpha1 The first non-negative numeric value of beta function for group 1.
#' @param alpha2 The first non-negative numeric value of beta function for group 2.
#' @param beta1 The second non-negative numeric value of beta function for group 1.
#' @param beta2 The second non-negative numeric value of beta function for group 2.
#'
#' @return The \code{ddiff2beta} gives the numeric value of the density.
#'
#' @examples
#' ddiff2beta(0.2, 0.5, 0.5, 0.5, 0.5)
#'
#' @export
ddiff2beta = function(pi, alpha1, alpha2, beta1, beta2) {
  k = 1 / (beta(alpha1, beta1) * beta(alpha2, beta2))
  gamma = alpha1 + alpha2 + beta1 + beta2 - 2
  if((pi >= -1) & (pi < 0)) {
    density = '*'(
      k * beta(alpha1, beta2) * (-pi) ^ (beta1 + beta2 - 1) * (1 + pi) ^ (alpha1 + beta2 - 1),
      AppellsF1(beta2, 1 - alpha2, gamma, alpha1 + beta2, 1 - pi ^ 2, 1 + pi)
    )
  } else if((pi >= 0) & (pi < 1)) {
    if ((pi == 0) & (alpha1 + alpha2 > 1) & (beta1 + beta2 > 1)) {
      density = k * beta(alpha1 + alpha2 - 1, beta1 + beta2 - 1)
    } else {
      density = '*'(
        k * beta(alpha2, beta1) * pi ^ (beta1 + beta2 - 1) * (1 - pi) ^ (alpha2 + beta1 - 1),
        AppellsF1(beta1, gamma, 1 - alpha1, alpha2 + beta1, 1 - pi, 1 - pi ^ 2)
      )
    }
  } else {
    density = 0
  }
  return(density)
}
