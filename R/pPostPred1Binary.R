#' Calculate Bayesian Posterior Probability or Bayesian Posterior Predictive Probability
#' for a Clinical Trial When Outcome is Single Binary
#'
#' This function computes Bayesian posterior probability or posterior predictive probability
#' for a binary outcome clinical trials. The function can handle controlled, uncontrolled, and
#' external control designs, using beta-binomial conjugate priors.
#'
#' @param prob A character string specifying the type of probability to calculate
#'        (\code{prob = 'posterior'} or \code{prob = 'predictive'}).
#' @param design A character string specifying the type of design
#'        (\code{design = 'controlled'}, \code{design = 'uncontrolled'}, or \code{design = 'external'}).
#' @param theta0 A numeric value representing the pre-specified threshold value.
#' @param n1 A positive integer representing the number of patients in group 1 for a proof-of-concept (PoC) trial.
#' @param n2 A positive integer representing the number of patients in group 2 for the PoC trial.
#' @param y1 A non-negative integer representing the observed number of responders in group 1 for the PoC trial.
#' @param y2 A non-negative integer representing the observed number of responders in group 2 for the PoC trial.
#' @param a1 A positive numeric value representing the first shape parameter of the prior distribution for group 1.
#' @param a2 A positive numeric value representing the first shape parameter of the prior distribution for group 2.
#' @param b1 A positive numeric value representing the second shape parameter of the prior distribution for group 1.
#' @param b2 A positive numeric value representing the second shape parameter of the prior distribution for group 2.
#' @param m1 A positive integer representing the number of patients in group 1 for the future trial data.
#' @param m2 A positive integer representing the number of patients in group 2 for the future trial data.
#' @param ne1 A positive integer representing the number of patients in group 1 for the external data.
#' @param ne2 A positive integer representing the number of patients in group 2 for the external data.
#' @param ye1 A non-negative integer representing the observed number of responders in group 1 for the external data.
#' @param ye2 A non-negative integer representing the observed number of responders in group 2 for the external data.
#' @param ae1 A positive numeric value representing the scale parameter (power parameter) for group 1.
#' @param ae2 A positive numeric value representing the scale parameter (power parameter) for group 2.
#' @param lower.tail logical; if TRUE (default), probabilities are P(theta <= theta0), otherwise, P(theta > theta0)
#'
#' @return A numeric value representing the Bayesian posterior probability or Bayesian posterior
#'         predictive probability.
#'
#' @details
#' The function can obtain:
#' \itemize{
#'   \item Bayesian posterior probability
#'   \item Bayesian posterior predictive probability
#' }
#'
#' The prior distribution of the proportion of responders (πij) for each treatment group (j=1,2) follows
#' a beta distribution. For posterior probability, the posterior distribution of πij follows a beta distribution.
#' For posterior predictive probability, the predictive distribution of future trial data follows
#' a beta-binomial distribution. The function can account for external (historical) data through
#' power priors.
#'
#' @examples
#' # Calculate posterior probability with external control
#' pPostPred1Binary(
#'   prob = 'posterior', design = 'external', theta0 = 0.15,
#'   n1 = 12, n2 = 15, y1 = 7, y2 = 9, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   m1 = NULL, m2 = NULL, ne1 = 12, ne2 = 12, ye1 = 6, ye2 = 6, ae1 = 0.5, ae2 = 0.5, lower.tail = FALSE
#' )
#'
#' # Calculate posterior predictive probability with external control
#' pPostPred1Binary(
#'   prob = 'predictive', design = 'external', theta0 = 0.5,
#'   n1 = 12, n2 = 15, y1 = 7, y2 = 7, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   m1 = 12, m2 = 12, ne1 = 12, ne2 = 12, ye1 = 6, ye2 = 6, ae1 = 0.5, ae2 = 0.5, lower.tail = FALSE
#' )
#'
#' @importFrom stats integrate
#' @export
pPostPred1Binary <- function(prob = 'posterior', design = 'controlled', theta0,
                             n1, n2, y1, y2, a1, a2, b1, b2,
                             m1, m2, ne1, ne2, ye1, ye2, ae1, ae2, lower.tail = TRUE) {
  # Check parameter sets for posterior predictive probability
  if((prob == 'predictive') & (sum(sapply(list(m1, m2), is.null)) > 0)) {
    stop('If you calculate posterior predictive probability, m1 and m2 should be non-null')
  }

  # Check parameter sets for external design
  if((design == 'external') & (sum(sapply(list(ne1, ne2, ye1, ye2, ae1, ae2), is.null)) > 0)) {
    stop('If you use the external data, ne1, ne2, ye1, ye2, ae1 and ae2 should be non-null')
  }

  # Set shape parameters for posterior distributions
  # Incorporate external data if design is 'external'
  s11 <- y1 + a1 + (design == 'external') * c(ye1 * ae1, 0)[1]
  s12 <- y2 + a2 + (design == 'external') * c(ye2 * ae2, 0)[1]
  s21 <- n1 - y1 + b1 + (design == 'external') * c((ne1 - ye1) * ae1, 0)[1]
  s22 <- n2 - y2 + b2 + (design == 'external') * c((ne2 - ye2) * ae2, 0)[1]

  if(prob == 'posterior') {
    # Calculate posterior probability using beta distributions
    results <- pBetadiff(theta0, s11, s12, s21, s22, lower.tail)
  } else if(prob == 'predictive') {
    # Calculate posterior predictive probability using beta-binomial distributions
    results <- pBetaBinomdiff(theta0, m1, m2, s11, s12, s21, s22, lower.tail)
  }

  return(results)
}
