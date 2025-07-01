#' Calculate the Go, NoGo and Gray Probabilities for a Clinical Trial When Outcome is Continuous
#' Under the Bayesian Framework Using Two Metrics
#'
#' This function calculates Go, NoGo, and Gray probabilities for continuous outcome clinical trials
#' under the Bayesian framework using two metrics: (i) posterior probability for the treatment
#' effect to be greater than a threshold, and (ii) posterior predictive probability of phase III
#' study success. The function supports controlled and uncontrolled designs with multiple calculation methods.
#'
#' @param nsim A positive integer representing the number of iterations for calculating posterior/posterior predictive probability.
#' @param prob A character string specifying the type of probability to use
#'        (\code{prob = 'posterior'} or \code{prob = 'predictive'}).
#' @param design A character string specifying the type of trial design
#'        (\code{design = 'controlled'} or \code{design = 'uncontrolled'}).
#' @param prior A character string specifying the prior distribution
#'        (\code{prior = 'N-Inv-Chisq'} or \code{prior = 'vague'}).
#' @param CalcMethod A character string specifying the calculation method
#'        (\code{CalcMethod = 'NI'} for numerical integration, \code{CalcMethod = 'MC'} for Monte Carlo method,
#'        or \code{CalcMethod = 'WS'} for Welch-Satterthwaite approximation).
#' @param theta.TV A numeric value representing the pre-specified threshold value for calculating
#'        Go probability when \code{prob = 'posterior'}.
#' @param theta.MAV A numeric value representing the pre-specified threshold value for calculating
#'        NoGo probability when \code{prob = 'posterior'}.
#' @param theta.NULL A numeric value representing the pre-specified threshold value for calculating
#'        Go/NoGo probabilities when \code{prob = 'predictive'}.
#' @param nMC A positive integer representing the number of iterations for Monte Carlo simulation
#'        (required only if \code{CalcMethod = 'MC'}).
#' @param gamma1 A numeric value between 0 and 1 representing the minimum probability to declare success.
#' @param gamma2 A numeric value between 0 and 1 representing the futility threshold.
#' @param n1 A positive integer representing the number of patients in group 1 for a proof-of-concept (PoC) trial.
#' @param n2 A positive integer representing the number of patients in group 2 for the PoC trial.
#' @param m1 A positive integer representing the number of patients in group 1 for the future trial data.
#' @param m2 A positive integer representing the number of patients in group 2 for the future trial data.
#' @param kappa01 A positive numeric value representing the prior precision parameter related to the mean
#'        for conjugate prior of Normal-Inverse-Chi-squared in group 1.
#' @param kappa02 A positive numeric value representing the prior precision parameter related to the mean
#'        for conjugate prior of Normal-Inverse-Chi-squared in group 2.
#' @param nu01 A positive numeric value representing the prior degrees of freedom related to the variance
#'        for conjugate prior of Normal-Inverse-Chi-squared in group 1.
#' @param nu02 A positive numeric value representing the prior degrees of freedom related to the variance
#'        for conjugate prior of Normal-Inverse-Chi-squared in group 2.
#' @param mu01 A numeric value representing the prior mean value of outcomes in group 1 for the PoC trial.
#' @param mu02 A numeric value representing the prior mean value of outcomes in group 2 for the PoC trial.
#' @param sigma01 A positive numeric value representing the prior standard deviation of outcomes in group 1 for the PoC trial.
#' @param sigma02 A positive numeric value representing the prior standard deviation of outcomes in group 2 for the PoC trial.
#' @param mu1 A numeric value representing the true mean of group 1 for PoC trial.
#' @param mu2 A numeric value representing the true mean of group 2 for PoC trial.
#' @param sigma1 A positive numeric value representing the true standard deviation of group 1 for PoC trial.
#' @param sigma2 A positive numeric value representing the true standard deviation of group 2 for PoC trial.
#' @param r A positive numeric value representing the parameter value associated with the distribution
#'        of mean for group 2 when \code{design = 'uncontrolled'}.
#' @param seed A numeric value representing the seed number for reproducible random number generation.
#'
#' @return A data frame containing the true means for both groups, and the Go, NoGo, and Gray probabilities.
#'
#' @details
#' The function can obtain:
#' \itemize{
#'   \item Go probability
#'   \item NoGo probability
#'   \item Gray probability
#' }
#'
#' The function can be used for controlled design and uncontrolled design. The decision framework is based on:
#' \itemize{
#'   \item Go: Probability that the treatment effect exceeds the efficacy threshold
#'   \item NoGo: Probability that the treatment effect is below the futility threshold
#'   \item Gray: Intermediate zone where neither Go nor NoGo criteria are met
#' }
#'
#' The function uses simulation to generate observed data and then applies Bayesian methods to calculate
#' decision probabilities. Three calculation methods are available for computing the underlying probabilities:
#' numerical integration (NI), Monte Carlo simulation (MC), and Welch-Satterthwaite approximation (WS).
#'
#' @examples
#' # Calculate Go/NoGo/Gray probabilities using posterior probability for controlled design
#' BayesDecisionProbContinuous(
#'   nsim = 1000, prob = 'posterior', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'NI',
#'   theta.TV = 2, theta.MAV = 0, theta.NULL = NULL, nMC = NULL, gamma1 = 0.8, gamma2 = 0.3,
#'   n1 = 12, n2 = 12, m1 = NULL, m2 = NULL, kappa01 = 5, kappa02 = 5, nu01 = 5, nu02 = 5,
#'   mu01 = 5, mu02 = 5, sigma01 = sqrt(5), sigma02 = sqrt(5), mu1 = 4, mu2 = 0,
#'   sigma1 = 1, sigma2 = 1, r = NULL, seed = 1
#' )
#'
#' @importFrom stats rnorm
#' @export
BayesDecisionProbContinuous <- function(nsim, prob, design, prior, CalcMethod, theta.TV, theta.MAV, theta.NULL, nMC,
                                        gamma1, gamma2, n1, n2, m1, m2, kappa01, kappa02, nu01, nu02, mu01, mu02,
                                        sigma01, sigma02, mu1, mu2, sigma1, sigma2, r, seed) {
  # Check parameter sets for posterior probability
  if((prob == 'posterior') & (sum(sapply(list(theta.TV, theta.MAV), is.null)) > 0)) {
    stop('If you calculate the Go, NoGo and Gray probabilities using posterior probability, theta.TV and theta.MAV should be non-null')
  }

  # Check parameter sets for posterior predictive probability
  if((prob == 'predictive') & (sum(sapply(list(m1, m2), is.null)) > 0)) {
    stop('If you calculate the Go, NoGo and Gray probabilities using posterior predictive probability, m1 and m2 should be non-null')
  }

  # Check parameter sets for posterior predictive probability threshold
  if((prob == 'predictive') & (is.null(theta.NULL))) {
    stop('If you calculate the Go, NoGo and Gray probabilities using posterior predictive probability, theta.NULL should be non-null')
  }

  # Check parameter sets for uncontrolled design
  if((design == 'uncontrolled') & (is.null(r))) {
    stop('If you consider uncontrolled design, r should be non-null')
  }

  # Set seed number for reproducible results
  set.seed(seed)

  # Generate random numbers of outcomes for group 1 in PoC study
  y1 <- matrix(rnorm(nsim * n1, mu1, sigma1), nrow = nsim)

  # Calculate sample mean for group 1
  bar.y1 <- rowSums(y1) / n1

  # Calculate sample standard deviation for group 1
  s1 <- sqrt(rowSums((y1 - bar.y1) ^ 2) / (n1 - 1))

  if(design == 'controlled') {
    # Generate random numbers of outcomes for group 2 in PoC study
    y2 <- matrix(rnorm(nsim * n2, mu2, sigma2), nrow = nsim)

    # Calculate sample mean for group 2
    bar.y2 <- rowSums(y2) / n2

    # Calculate sample standard deviation for group 2
    s2 <- sqrt(rowSums((y2 - bar.y2) ^ 2) / (n2 - 1))
  } else if(design == 'uncontrolled') {
    # Set sample mean and standard deviation for group 2 to NULL for uncontrolled design
    bar.y2 <- NULL
    s2 <- NULL
  }

  # Set values of theta0 based on probability type
  if(prob == 'posterior') {
    theta0 <- c(theta.TV, theta.MAV)
  } else {
    theta0 <- theta.NULL
  }

  # Calculate Go, NoGo and Gray probabilities
  list.Go.and.NoGo.probs <- lapply(seq(theta0), function(i) {
    # Calculate probability of success for each simulated dataset
    prob.success <- BayesPostPredContinuous(
      prob, design, prior, CalcMethod, theta0[i], nMC, n1, n2, m1, m2,
      kappa01, kappa02, nu01, nu02, mu01, mu02, sigma01, sigma02,
      bar.y1, bar.y2, s1, s2, r
    )

    if(prob == 'posterior') {
      # Calculate Go probability (for theta.TV) or NoGo probability (for theta.MAV)
      Go   <- ifelse(i == 1, 1, NA) * sum(prob.success >= gamma1) / nsim
      NoGo <- ifelse(i == 1, NA, 1) * sum(prob.success <= gamma2) / nsim
    } else {
      # Calculate Go and NoGo probabilities for posterior predictive
      Go <- sum(prob.success >= gamma1) / nsim
      NoGo <- sum(prob.success <= gamma2) / nsim
    }

    return(c(Go, NoGo))
  })

  # Combine Go and NoGo probabilities across thresholds
  Go.and.NoGo.probs <- do.call(pmax, c(list(na.rm = TRUE), list.Go.and.NoGo.probs))

  # Calculate Gray probability (complement of Go and NoGo)
  Gray.prob <- 1 - sum(Go.and.NoGo.probs)

  # Check for negative Gray probabilities
  if(Gray.prob < 0) {
    print('Because negative gray probability(s) is obtained, re-consider appropriate threshold')
  }

  # Prepare results data frame
  if(is.null(mu2)) {
    # For uncontrolled design, only include mu1
    results <- data.frame(mu1, Go = Go.and.NoGo.probs[1], NoGo = Go.and.NoGo.probs[2], Gray = Gray.prob)
  } else {
    # For controlled design, include both mu1 and mu2
    results <- data.frame(mu1, mu2, Go = Go.and.NoGo.probs[1], NoGo = Go.and.NoGo.probs[2], Gray = Gray.prob)
  }

  return(results)
}
