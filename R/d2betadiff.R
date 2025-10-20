#' Probability Density Function of the Difference Between Two Beta Variables
#'
#' This function calculates the probability density function (PDF) of the difference
#' between two independent beta-distributed random variables. The computation uses
#' Appell's first hypergeometric function for numerical evaluation. Specifically,
#' it computes the density f(π) where π = X1 - X2, with X1 ~ Beta(alpha1, beta1)
#' and X2 ~ Beta(alpha2, beta2).
#'
#' @param pi A numeric value in the range \code{[-1, 1]} representing the difference
#'        between the two beta variables.
#' @param alpha1 A positive numeric value representing the first shape parameter
#'        of the first beta distribution.
#' @param alpha2 A positive numeric value representing the first shape parameter
#'        of the second beta distribution.
#' @param beta1 A positive numeric value representing the second shape parameter
#'        of the first beta distribution.
#' @param beta2 A positive numeric value representing the second shape parameter
#'        of the second beta distribution.
#'
#' @return A numeric value representing the probability density at the specified
#'         difference value. Returns 0 if pi is outside the range \code{[-1, 1]}.
#'
#' @details
#' The function evaluates the density using different formulations depending on
#' the range of π:
#' \itemize{
#'   \item For π ∈ \code{[-1, 0)}: Uses Appell's F1 function with parameters adjusted
#'         for negative differences
#'   \item For π ∈ \code{[0, 1)}: Uses Appell's F1 function with parameters adjusted
#'         for positive differences
#'   \item Special case at π = 0: When alpha1 + alpha2 > 1 and beta1 + beta2 > 1,
#'         uses a simplified beta function form
#'   \item For π ∉ \code{[-1, 1]}: Returns 0 (outside support of the distribution)
#' }
#'
#' The normalization constant k is computed as:
#' \deqn{k = \frac{1}{B(\alpha_1, \beta_1) B(\alpha_2, \beta_2)}}
#'
#' where B(·,·) is the beta function.
#'
#' @examples
#' # Calculate density at π = 0.2 for symmetric beta distributions
#' d2betadiff(0.2, 0.5, 0.5, 0.5, 0.5)
#'
#' # Calculate density at π = 0 for Beta(2,3) - Beta(1,4)
#' d2betadiff(0, 2, 3, 1, 4)
#'
#' # Calculate density at π = -0.3 for informative priors
#' d2betadiff(-0.3, 3, 2, 2, 3)
#'
#' # Verify that density is 0 outside [-1, 1]
#' d2betadiff(1.5, 1, 1, 1, 1)
#'
#' @export
d2betadiff <- function(pi, alpha1, alpha2, beta1, beta2) {
  # Calculate the normalization constant using beta functions
  k <- 1 / (beta(alpha1, beta1) * beta(alpha2, beta2))

  # Calculate gamma parameter for Appell's function
  gamma <- alpha1 + alpha2 + beta1 + beta2 - 2

  # Case 1: pi in [-1, 0) - negative difference region
  if((pi >= -1) & (pi < 0)) {
    density <- '*'(
      k * beta(alpha1, beta2) * (-pi) ^ (beta1 + beta2 - 1) * (1 + pi) ^ (alpha1 + beta2 - 1),
      AppellsF1(beta2, 1 - alpha2, gamma, alpha1 + beta2, 1 - pi ^ 2, 1 + pi)
    )
  }
  # Case 2: pi in [0, 1) - non-negative difference region
  else if((pi >= 0) & (pi < 1)) {
    # Special case at pi = 0 when certain parameter conditions are met
    if ((pi == 0) & (alpha1 + alpha2 > 1) & (beta1 + beta2 > 1)) {
      density <- k * beta(alpha1 + alpha2 - 1, beta1 + beta2 - 1)
    } else {
      density <- '*'(
        k * beta(alpha2, beta1) * pi ^ (beta1 + beta2 - 1) * (1 - pi) ^ (alpha2 + beta1 - 1),
        AppellsF1(beta1, gamma, 1 - alpha1, alpha2 + beta1, 1 - pi, 1 - pi ^ 2)
      )
    }
  }
  # Case 3: pi outside [-1, 1] - outside support of distribution
  else {
    density <- 0
  }

  return(density)
}
