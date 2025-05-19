#' Calculate bayesian posterior probability or bayesian posterior predictive probability
#' for a clinical trial when outcome is continuous.
#'
#' The function can obtain:
#' \itemize{
#'   \item Bayesian posterior probability
#'   \item Bayesian posterior predictive probability
#' }
#' Prior distribution of mean and variance of outcomes for each treatment group (k=1,2) is (1) Normal-Inverse-X^2 or (2) Vague.
#' Posterior distribution or posterior predictive distribution of outcome for each treatment group is following t-distribution.
#'
#' @param msim A numeric value of a number of data generation in PoC trial.
#' @param design A type of trial design (i.e., \code{design='controlled'} or \code{prob='uncontrolled'}).
#' @param prob A type of probability (i.e., \code{prob='posterior'} or \code{prob='predictive'}).
#' @param prior A prior distribution (i.e., \code{prior='N-Inv-Chisq'} or \code{prior='vague'}).
#' @param theta0 Numeric pre-specified threshold value(s) (user can set multiple values).
#' @param n1 A number of patients in group 1 for a proof-of-concept (PoC) trial.
#' @param n2 A number of patients in group 2 for the PoC trial.
#' @param m1 A number of patients in group 1 for the future trial data.
#' @param m2 A number of patients in group 2 for the future trial data.
#' @param kappa01 A true value related to the mean for conjugate prior of Normal-Inverse-X^2 in group 1.
#' @param kappa02 A true value related to the mean for conjugate prior of Normal-Inverse-X^2 in group 2.
#' @param nu01 A true value related to the variance for conjugate prior of Normal-Inverse-X^2 in group 1.
#' @param nu02 A true value related to the variance for conjugate prior of Normal-Inverse-X^2 in group 2.
#' @param mu01 A true mean value of outcomes in group 1 for the PoC trial.
#' @param mu02 A true mean value of outcomes in group 2 for the PoC trial.
#' @param sigma01 A true standard deviation value of outcomes in group 1 for the PoC trial.
#' @param sigma02 A true standard deviation value of outcomes in group 2 for the PoC trial.
#' @param bar.y1 A sample mean of group 1.
#' @param bar.y2 A sample mean of group 2.
#' @param s1 A sample standard deviation of group 1.
#' @param s2 A sample standard deviation of group 2.
#' @param r A parameter value associated with the distribution of mean for group 2 When \code{prob='uncontrolled'}.
#' @param seed A numeric value of a seed number for data generation.
#'
#' @return The \code{BayesPostPredContinuous} gives the numeric value of the
#' bayesian posterior probability or bayesian posterior predictive probability.
#'
#' @examples
#' BayesPostPredContinuous(
#'   msim = 10000, design = 'controlled', prob = 'posterior', prior = 'N-Inv-Chisq', theta0 = c(2, 0), n1 = 12, n2 = 12, m1 = NULL, m2 = NULL,
#'   kappa01 = 5, kappa02 = 5, nu01 = 5, nu02 = 5, mu01 = 5, mu02 = 5, sigma01 = sqrt(5), sigma02 = sqrt(5),
#'   bar.y1 = 2, bar.y2 = 0, s1 = 1, s2 = 1, r = NULL, seed = 1
#' )
#' BayesPostPredContinuous(
#'   msim = 10000, design = 'uncontrolled', prob = 'predictive', prior = 'vague', theta0 = 0.5, n1 = 12, n2 = NULL, m1 = 120, m2 = 120,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL, mu01 = NULL, mu02 = 0, sigma01 = NULL, sigma02 = NULL,
#'   bar.y1 = 2, bar.y2 = NULL, s1 = 1, s2 = NULL, r = 12, seed = 1
#' )
#'
#' @importFrom stats rt
#' @export
BayesPostPredContinuous = function(msim, design, prob, prior, theta0, n1, n2, m1, m2,
                                   kappa01, kappa02, nu01, nu02, mu01, mu02, sigma01, sigma02,
                                   bar.y1, bar.y2, s1, s2, r, seed) {
  # Check parameter sets
  if((design == 'uncontrolled') & (is.null(r))) {
    stop('If you consider uncontrolled design, r should be non-null')
  }
  if((prob == 'predictive') & (sum(sapply(list(m1, m2), is.null)) > 0)) {
    stop('If you calculate posterior predictive probability, m1 and m2 should be non-null')
  }
  if((prior == 'N-Inv-Chisq') & (sum(sapply(list(kappa01, nu01, sigma01), is.null)) > 0)) {
    stop('If you use the N-Inv-Chisq prior, kappa01, nu01 and sigma01 should be non-null')
  }
  # Set seed number
  set.seed(seed)
  # Indicator function that prior distribution is an N-Inv-Chisq
  I.NIC = rep(as.logical(prior == 'N-Inv-Chisq'), 2)
  # Indicator function that the type of probability is a posterior
  I.POS = as.logical(prob == 'posterior')
  # Renaming parameters
  nk = c(n1, n2); mu0k = c(mu01, mu02); kappa0k = c(kappa01, kappa02); nu0k = c(nu01, nu02);
  sigma0k = c(sigma01, sigma02); bar.yk = c(bar.y1, bar.y2); sk = c(s1, s2)
  # Sample size
  kappa.nk = nk + ifelse(I.NIC, kappa0k, 0)
  # Degree of freedom
  nu.nk = nk - 1 + ifelse(I.NIC, nu0k + 1, 0)
  # Mean of a posterior distribution
  mu.nk = bar.yk * ifelse(I.NIC, nk / kappa.nk, 1) + ifelse(I.NIC, kappa0k * mu0k / kappa.nk, 0)
  # Variance of a posterior distribution
  var.nk = '+'(
    sk ^ 2 * ifelse(I.NIC, (nk - 1) / nu.nk, 1),
    ifelse(I.NIC, (nu0k * sigma0k ^ 2 + nk * kappa0k / kappa.nk * (mu0k - bar.yk) ^ 2) / nu.nk, 0)
  )
  ## Generate t-distributed random variables
  # Group 1
  rand.t1 = sqrt(var.nk[1] / kappa.nk[1]) * ifelse(I.POS, 1, sqrt((1 + kappa.nk[1]) / m1)) * rt(msim, df = nu.nk[1], ncp = 0)
  # Group 2
  if(design == 'controlled') {
    rand.t2 = sqrt(var.nk[2] / kappa.nk[2]) * ifelse(I.POS, 1, sqrt((1 + kappa.nk[2]) / m2)) * rt(msim, df = nu.nk[2], ncp = 0)
    # Calculate theta values (= t_{df_{1}} - t_{df_{2}})
    theta =  rand.t1 - rand.t2 + diff(rev(mu.nk))
  } else if(design == 'uncontrolled') {
    # Calculate theta values (= t_{df_{1}} - t_{df_{2}})
    theta = ifelse(I.POS, sqrt(r + 1), sqrt(r * m1 / m2 + 1)) * rand.t1 + mu.nk[1] - mu02
  }
  # The probability of exceeding \theta_{0}
  g = sapply(theta0, function(i) sum(c(theta) > i) / msim)
  return(g)
}
