#' Cumulative Distribution Function of the Difference Between Two Beta Variables
#'
#' This function calculates the cumulative distribution function (CDF) of the difference
#' between two independent beta-distributed random variables using Appell's first
#' hypergeometric function. Specifically, it computes P(X1 - X2 > q) where X1 ~ Beta(alpha1, beta1)
#' and X2 ~ Beta(alpha2, beta2).
#'
#' @param q A numeric value representing the quantile threshold.
#' @param alpha1 A positive numeric value representing the first shape parameter of the first beta distribution.
#' @param alpha2 A positive numeric value representing the first shape parameter of the second beta distribution.
#' @param beta1 A positive numeric value representing the second shape parameter of the first beta distribution.
#' @param beta2 A positive numeric value representing the second shape parameter of the second beta distribution.
#'
#' @return A numeric value representing P(X1 - X2 > q), the probability that the difference
#'         between the two beta variables exceeds the quantile q.
#'
#' @details
#' The function uses a piecewise approach based on the range of the difference π = X1 - X2:
#' \itemize{
#'   \item For π ∈ [-1, 0): Uses Appell's F1 function with specific parameters
#'   \item For π ∈ [0, 1): Uses a different parameterization of Appell's F1 function
#'   \item Special case handling when π = 0 and certain parameter conditions are met
#' }
#'
#' @examples
#' # Calculate P(Beta(0.5, 0.5) - Beta(0.5, 0.5) > 0.2)
#' pBetadiff(0.2, 0.5, 0.5, 0.5, 0.5)
#'
#' # Calculate P(Beta(2, 3) - Beta(1, 4) > -0.1)
#' pBetadiff(-0.1, 2, 3, 1, 4)
#'
#' # Calculate P(Beta(1, 1) - Beta(1, 1) > 0)
#' pBetadiff(0, 1, 1, 1, 1)
#'
#' @importFrom stats integrate
#' @export
pBetadiff <- function(q, alpha1, alpha2, beta1, beta2) {
  # Calculate the normalization constant using beta functions
  k <- 1 / (beta(alpha1, beta1) * beta(alpha2, beta2))

  # Calculate gamma parameter for Appell's function
  gamma <- alpha1 + alpha2 + beta1 + beta2 - 2

  # Compute the posterior probability using numerical integration
  results <- integrate(
    Vectorize(function(pi) {
      # Case 1: pi in [-1, 0)
      if((pi >= -1) & (pi < 0)) {
        density <- '*'(
          k * beta(alpha1, beta2) * (-pi) ^ (beta1 + beta2 - 1) * (1 + pi) ^ (alpha1 + beta2 - 1),
          AppellsF1(beta2, 1 - alpha2, gamma, alpha1 + beta2, 1 - pi ^ 2, 1 + pi)
        )
      }
      # Case 2: pi in [0, 1)
      else if((pi >= 0) & (pi < 1)) {
        # Special case: pi = 0 with specific parameter conditions
        if ((pi == 0) & (alpha1 + alpha2 > 1) & (beta1 + beta2 > 1)) {
          density <- k * beta(alpha1 + alpha2 - 1, beta1 + beta2 - 1)
        } else {
          # General case for pi in (0, 1)
          density <- '*'(
            k * beta(alpha2, beta1) * pi ^ (beta1 + beta2 - 1) * (1 - pi) ^ (alpha2 + beta1 - 1),
            AppellsF1(beta1, gamma, 1 - alpha1, alpha2 + beta1, 1 - pi, 1 - pi ^ 2)
          )
        }
      }
      # Case 3: pi outside [-1, 1]
      else {
        density <- 0
      }

      return(density)
    }),
    # Integration bounds: from q to 1 (upper tail probability)
    lower = q,
    upper = 1
  )[['value']]

  return(results)
}
