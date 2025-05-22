#' Calculate the Go, NoGo and Gray probabilities for a clinical trial when outcome is continuous
#' under the Bayesian framework using two metrics:
#' (i)  posterior probability for the treatment effect to be greater than a threshold, and
#' (ii) posterior predictive probability of phase III study success.
#'
#' The function can obtain:
#' \itemize{
#'   \item Go probability
#'   \item NoGo probability
#'   \item Gray probability
#' }
#' The function can be used for controlled design and uncontrolled design.
#' @param nsim A numeric value of the number of iterations.
#' @param design A type of trial design (i.e., \code{design='controlled'} or \code{prob='uncontrolled'}).
#' @param prob A type of probability (i.e., \code{prob='posterior'} or \code{prob='predictive'}).
#' @param prior A prior distribution (i.e., \code{prior='N-Inv-Chisq'} or \code{prior='vague'}).
#' @param theta0 Numeric pre-specified threshold value(s) (user can set multiple values).
#' @param gamma1 A numeric value of a minimum probability to declare success.
#' @param gamma2 A numeric value of a futility threshold.
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
#' @param mu1 A true mean of group 1 for PoC trial.
#' @param mu2 A true mean of group 2 for PoC trial.
#' @param sigma1 A true standard deviation of group 1 for PoC trial.
#' @param sigma2 A true standard deviation of group 2 for PoC trial.
#' @param r A parameter value associated with the distribution of mean for group 2 When \code{prob='uncontrolled'}.
#' @param seed A numeric value of the seed number.
#'
#' @return The \code{BayesDecisionProbContinuous} gives the dataframe including the
#'  true means for both groups, and the Go, NoGo and Gray probabilities.
#'
#' @examples
#' BayesDecisionProbContinuous(
#'   nsim = 1000, design = 'controlled', prob = 'posterior', prior = 'N-Inv-Chisq', theta0 = c(2, 0), gamma1 = 0.8, gamma2 = 0.3,
#'   n1 = 12, n2 = 12, m1 = NULL, m2 = NULL, kappa01 = 5, kappa02 = 5, nu01 = 5, nu02 = 5, mu01 = 5, mu02 = 5, sigma01 = sqrt(5), sigma02 = sqrt(5),
#'   mu1 = 4, mu2 = 0, sigma1 = 1, sigma2 = 1, r = NULL, seed = 1
#' )
#'
#' @importFrom stats rnorm
#' @export
BayesDecisionProbContinuous = function(nsim, design, prob, prior, theta0, gamma1, gamma2,
                                       n1, n2, m1, m2, kappa01, kappa02, nu01, nu02, mu01, mu02, sigma01, sigma02,
                                       mu1, mu2, sigma1, sigma2, r, seed) {
  # Set seed number
  set.seed(seed)
  # Random numbers of outcomes for group 1 in PoC study
  y1 = matrix(rnorm(nsim * n1, mu1, sigma1), nrow = nsim)
  # Sample mean for group 1
  bar.y1 = rowSums(y1) / n1
  # Sample standard deviation for group 1
  s1 = sqrt(rowSums((y1 - bar.y1) ^ 2) / (n1 - 1))
  if(design == 'controlled') {
    # Random numbers of outcomes for group 2 in PoC study
    y2 = matrix(rnorm(nsim * n2, mu2, sigma2), nrow = nsim)
    # Sample mean for group 2
    bar.y2 = rowSums(y2) / n2
    # Sample standard deviation for group 2
    s2 = sqrt(rowSums((y2 - bar.y2) ^ 2) / (n2 - 1))
  } else if(design == 'uncontrolled') {
    # Sample mean for group 2
    bar.y2 = NULL
    # Sample standard deviation for group 2
    s2 = NULL
  }
  # Probability of success
  prob.success = BayesPostPredContinuous(
    design, prob, prior, theta0, n1, n2, m1, m2,
    kappa01, kappa02, nu01, nu02, mu01, mu02, sigma01, sigma02,
    bar.y1, bar.y2, s1, s2, r
  )
  if(prob == 'posterior') {
    # Go, NoGo and Gray probabilities
    Go = sum(prob.success[, 1] >= gamma1) / nsim
    NoGo = sum(prob.success[, 2] <= gamma2) / nsim
    Gray = 1 - Go - NoGo
  } else if(prob == 'predictive') {
    # Go, NoGo and Gray probabilities
    Go = sum(prob.success >= gamma1) / nsim
    NoGo = sum(prob.success <= gamma2) / nsim
    Gray = 1 - Go - NoGo
  }
  if(sum(Gray < 0) > 0) {
    stop('Because negative gray probability(s) is obtained, re-consider appropriate threshold')
  }
  # Return output
  if(is.null(mu2)) {
    results = data.frame(mu1, Go, NoGo, Gray)
  } else {
    results = data.frame(mu1, mu2, Go, NoGo, Gray)
  }
  return(results)
}
