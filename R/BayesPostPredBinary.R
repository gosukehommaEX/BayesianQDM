#' Calculate bayesian posterior probability or bayesian posterior predictive probability
#' for a clinical trial when outcome is binary.
#'
#' The function can obtain:
#' \itemize{
#'   \item Bayesian posterior probability
#'   \item Bayesian posterior predictive probability
#' }
#' Prior distribution of proportion of responders (pi_{j}) for each treatment group (j=1,2) is
#' following beta distribution.
#' For posterior probability, posterior distribution of pi_{j} is following beta distribution.
#' For posterior predictive probability, predictive distribution of future trial data is
#' following beta-binomial distribution.
#' The function can account for external (historical) data.
#'
#' @param prob A type of probability (i.e., \code{prob='posterior'} or \code{prob='predictive'}).
#' @param external When external data is used, \code{external=TRUE}.
#' @param theta0 A numeric pre-specified threshold value.
#' @param n1 A number of patients in group 1 for a proof-of-concept (PoC) trial.
#' @param n2 A number of patients in group 2 for the PoC trial.
#' @param y1 An observed number of responders in group 1 for the PoC trial.
#' @param y2 An observed number of responders in group 2 for the PoC trial.
#' @param a1 The first non-negative numeric value of the prior distribution for group 1.
#' @param a2 The first non-negative numeric value of the prior distribution for group 2.
#' @param b1 The second non-negative numeric value of the prior distribution for group 1.
#' @param b2 The second non-negative numeric value of the prior distribution for group 2.
#' @param m1 A number of patients in group 1 for the future trial data.
#' @param m2 A number of patients in group 2 for the future trial data.
#' @param ne1 A number of patients in group 1 for the external data.
#' @param ne2 A number of patients in group 2 for the external data.
#' @param ye1 An observed number of responders in group 1 for the external data.
#' @param ye2 An observed number of responders in group 2 for the external data.
#' @param ae1 A scale parameter (power parameter) for group 1.
#' @param ae2 A scale parameter (power parameter) for group 2.
#'
#' @return The \code{BayesPostPredBinary} gives the numeric value of the
#' bayesian posterior probability or bayesian posterior predictive probability.
#'
#' @examples
#' BayesPostPredBinary(
#'   prob = 'posterior', external = TRUE, theta0 = 0.15,
#'   n1 = 12, n2 = 15, y1 = 7, y2 = 9, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   m1 = NULL, m2 = NULL, ne1 = 12, ne2 = 12, ye1 = 6, ye2 = 6, ae1 = 0.5, ae2 = 0.5
#' )
#' BayesPostPredBinary(
#'   prob = 'predictive', external = TRUE, theta0 = 0.5,
#'   n1 = 12, n2 = 15, y1 = 7, y2 = 7, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
#'   m1 = 12, m2 = 12, ne1 = 12, ne2 = 12, ye1 = 6, ye2 = 6, ae1 = 0.5, ae2 = 0.5
#' )
#'
#' @importFrom stats integrate
#' @importFrom VGAM dbetabinom.ab
#' @export
BayesPostPredBinary = function(prob, external, theta0,
                               n1, n2, y1, y2, a1, a2, b1, b2,
                               m1, m2, ne1, ne2, ye1, ye2, ae1, ae2) {
  # Check parameter sets
  if((prob == 'predictive') & (sum(sapply(list(m1, m2), is.null)) > 0)) {
    stop('If you calculate posterior predictive probability, m1 and m2 should be non-null')
  }
  if((external == TRUE) & (sum(sapply(list(ne1, ne2, ye1, ye2, ae1, ae2), is.null)) > 0)) {
    stop('If you use the external data, ne1, ne2, ye1, ye2, ae1 and ae2 should be non-null')
  }
  # Set shape parameters
  s11 = y1 + a1 + external * c(ye1 * ae1, 0)[1]
  s12 = y2 + a2 + external * c(ye2 * ae2, 0)[1]
  s21 = n1 - y1 + b1 + external * c((ne1 - ye1) * ae1, 0)[1]
  s22 = n2 - y2 + b2 + external * c((ne2 - ye2) * ae2, 0)[1]
  if(prob == 'posterior') {
    # A posterior probability
    g = integrate(Vectorize(function(theta) ddiff2beta(theta, s11, s12, s21, s22)), theta0, 1)[['value']]
  } else if(prob == 'predictive') {
    # Probability mass functions of beta-binomial distribution
    dbetabinom1 = VGAM::dbetabinom.ab(0:m1, m1, s11, s21)
    dbetabinom2 = VGAM::dbetabinom.ab(0:m2, m2, s12, s22)
    # A posterior predictive probability
    I = (outer(m2 * 0:m1, m1 * 0:m2, '-') / (m1 * m2) > theta0)
    g = as.numeric(crossprod(dbetabinom1[row(I)[I]], dbetabinom2[col(I)[I]]))
  }
  return(g)
}
