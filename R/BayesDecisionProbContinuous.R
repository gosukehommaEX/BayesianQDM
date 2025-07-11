#' Calculate the Go, NoGo and Gray Probabilities for a Clinical Trial When Outcome is Continuous
#' Under the Bayesian Framework Using Two Metrics
#'
#' @description
#' This function calculates Go, NoGo, and Gray probabilities for continuous outcome clinical trials
#' under the Bayesian framework using two metrics: (i) posterior probability for the treatment
#' effect to be greater than a threshold, and (ii) posterior predictive probability of phase III
#' study success. The function supports controlled, uncontrolled, and external control designs with
#' multiple calculation methods.
#'
#' @param nsim A positive integer representing the number of simulations.
#' @param prob A character string specifying the type of probability to calculate
#' ('posterior' or 'predictive').
#' @param design A character string specifying the type of trial design
#' ('controlled', 'uncontrolled', or 'external').
#' @param prior A character string specifying the prior distribution
#' ('N-Inv-Chisq' or 'vague').
#' @param CalcMethod A character string specifying the calculation method
#' ('NI' for numerical integration, 'MC' for Monte Carlo method,
#' 'WS' for Welch-Satterthwaite approximation, or 'MCMC' for MCMC).
#' @param theta.TV A numeric value representing the target value threshold.
#' @param theta.MAV A numeric value representing the minimum acceptable value threshold.
#' @param theta.NULL A numeric value representing the null hypothesis value (default: NULL).
#' @param nMC A positive integer representing the number of iterations for Monte Carlo simulation
#' (required only if CalcMethod = 'MC').
#' @param nMCMCsample A positive integer representing the number of iterations for MCMC sampling
#' (required only if CalcMethod = 'MCMC').
#' @param gamma1 A numeric value representing the go threshold.
#' @param gamma2 A numeric value representing the no-go threshold.
#' @param n1 A positive integer representing the number of patients in group 1.
#' @param n2 A positive integer representing the number of patients in group 2.
#' @param m1 A positive integer representing the number of patients in group 1 for future trial.
#' @param m2 A positive integer representing the number of patients in group 2 for future trial.
#' @param kappa01 A positive numeric value representing the prior sample size for group 1 (N-Inv-Chisq prior).
#' @param kappa02 A positive numeric value representing the prior sample size for group 2 (N-Inv-Chisq prior).
#' @param nu01 A positive numeric value representing the prior degrees of freedom for group 1 (N-Inv-Chisq prior).
#' @param nu02 A positive numeric value representing the prior degrees of freedom for group 2 (N-Inv-Chisq prior).
#' @param mu01 A numeric value representing the prior mean for group 1 (N-Inv-Chisq prior).
#' @param mu02 A numeric value representing the prior mean for group 2 (N-Inv-Chisq prior).
#' @param sigma01 A positive numeric value representing the prior standard deviation for group 1 (N-Inv-Chisq prior).
#' @param sigma02 A positive numeric value representing the prior standard deviation for group 2 (N-Inv-Chisq prior).
#' @param mu1 A numeric value representing the true mean of group 1.
#' @param mu2 A numeric value representing the true mean of group 2.
#' @param sigma1 A positive numeric value representing the true standard deviation of group 1.
#' @param sigma2 A positive numeric value representing the true standard deviation of group 2.
#' @param r A positive numeric value for uncontrolled design (default: NULL).
#' @param ne1 A positive integer representing the sample size for group 1 in external trial (default: NULL).
#' @param ne2 A positive integer representing the sample size for group 2 in external trial (default: NULL).
#' @param alpha01 A positive numeric value representing the scale parameter of the power prior for group 1
#' (required for external design, can be NULL if no external treatment data).
#' @param alpha02 A positive numeric value representing the scale parameter of the power prior for group 2
#' (required for external design, can be NULL if no external control data).
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
#' The function can be used for controlled design, uncontrolled design, and external control design.
#' The decision framework is based on:
#' \itemize{
#'   \item Go: Probability that the treatment effect exceeds the efficacy threshold
#'   \item NoGo: Probability that the treatment effect is below the futility threshold
#'   \item Gray: Intermediate zone where neither Go nor NoGo criteria are met
#' }
#'
#' The function uses simulation to generate observed data and then applies Bayesian methods to calculate
#' decision probabilities. Four calculation methods are available for computing the underlying probabilities:
#' numerical integration (NI), Monte Carlo simulation (MC), Welch-Satterthwaite approximation (WS),
#' and MCMC for external data incorporation.
#'
#' @examples
#' # Example 1: Numerical Integration (NI) method
#' BayesDecisionProbContinuous(
#'   nsim = 100, prob = 'posterior', design = 'controlled',
#'   prior = 'N-Inv-Chisq', CalcMethod = 'NI',
#'   theta.TV = 2, theta.MAV = 0, theta.NULL = NULL,
#'   nMC = NULL, nMCMCsample = NULL,
#'   gamma1 = 0.8, gamma2 = 0.3,
#'   n1 = 12, n2 = 12, m1 = NULL, m2 = NULL,
#'   kappa01 = 5, kappa02 = 5, nu01 = 5, nu02 = 5,
#'   mu01 = 5, mu02 = 5, sigma01 = sqrt(5), sigma02 = sqrt(5),
#'   mu1 = 4, mu2 = 0, sigma1 = 1, sigma2 = 1,
#'   r = NULL, ne1 = NULL, ne2 = NULL,
#'   alpha01 = NULL, alpha02 = NULL, seed = 1
#' )
#'
#' # Example 2: Monte Carlo (MC) method
#' BayesDecisionProbContinuous(
#'   nsim = 100, prob = 'posterior', design = 'controlled',
#'   prior = 'vague', CalcMethod = 'MC',
#'   theta.TV = 1.5, theta.MAV = -1, theta.NULL = NULL,
#'   nMC = 1e+4, nMCMCsample = NULL,
#'   gamma1 = 0.8, gamma2 = 0.2,
#'   n1 = 12, n2 = 12, m1 = NULL, m2 = NULL,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
#'   mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
#'   mu1 = 2, mu2 = 0, sigma1 = 1, sigma2 = 1,
#'   r = NULL, ne1 = NULL, ne2 = NULL,
#'   alpha01 = NULL, alpha02 = NULL, seed = 2
#' )
#'
#' # Example 3: Welch-Satterthwaite (WS) method
#' BayesDecisionProbContinuous(
#'   nsim = 100, prob = 'predictive', design = 'controlled',
#'   prior = 'vague', CalcMethod = 'WS',
#'   theta.TV = NULL, theta.MAV = NULL, theta.NULL = 1,
#'   nMC = NULL, nMCMCsample = NULL,
#'   gamma1 = 0.8, gamma2 = 0.2,
#'   n1 = 12, n2 = 12, m1 = 24, m2 = 24,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
#'   mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
#'   mu1 = 2.5, mu2 = 1.2, sigma1 = 1, sigma2 = 1,
#'   r = NULL, ne1 = NULL, ne2 = NULL,
#'   alpha01 = NULL, alpha02 = NULL, seed = 3
#' )
#'
#' \dontrun{
#' # Example 4: MCMC method with external control data
#' BayesDecisionProbContinuous(
#'   nsim = 100, prob = 'posterior', design = 'external',
#'   prior = 'vague', CalcMethod = 'MCMC',
#'   theta.TV = 1, theta.MAV = -1, theta.NULL = NULL,
#'   nMC = NULL, nMCMCsample = 3000,
#'   gamma1 = 0.8, gamma2 = 0.2,
#'   n1 = 12, n2 = 12, m1 = NULL, m2 = NULL,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
#'   mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
#'   mu1 = 2, mu2 = 0, sigma1 = 1, sigma2 = 1,
#'   r = NULL, ne1 = NULL, ne2 = 20,
#'   alpha01 = NULL, alpha02 = 0.5, seed = 4
#' )
#' }
#'
#' @importFrom stats rnorm
#' @export
BayesDecisionProbContinuous <- function(nsim, prob, design, prior, CalcMethod, theta.TV, theta.MAV, theta.NULL,
                                        nMC = NULL, nMCMCsample = NULL,
                                        gamma1, gamma2, n1, n2, m1, m2, kappa01, kappa02, nu01, nu02,
                                        mu01, mu02, sigma01, sigma02, mu1, mu2, sigma1, sigma2,
                                        r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL, seed) {

  # Set seed for reproducibility
  set.seed(seed)

  # Check parameter sets
  if (!prob %in% c("posterior", "predictive")) {
    stop("prob must be either 'posterior' or 'predictive'")
  }

  if (!design %in% c("controlled", "uncontrolled", "external")) {
    stop("design must be 'controlled', 'uncontrolled', or 'external'")
  }

  if (!prior %in% c("vague", "N-Inv-Chisq")) {
    stop("prior must be either 'vague' or 'N-Inv-Chisq'")
  }

  if (!CalcMethod %in% c("NI", "MC", "WS", "MCMC")) {
    stop("CalcMethod must be 'NI', 'MC', 'WS', or 'MCMC'")
  }

  # For external design, only MCMC is supported
  if (design == "external" && CalcMethod != "MCMC") {
    stop("For external design, CalcMethod must be 'MCMC'. Other methods are not supported for external data incorporation.")
  }

  # For non-external designs, MCMC is not supported
  if (design != "external" && CalcMethod == "MCMC") {
    stop("MCMC method is only available for external design.")
  }

  # Validate required parameters for each method
  if (CalcMethod == "MC" && is.null(nMC)) {
    stop("nMC must be specified for Monte Carlo method")
  }

  if (CalcMethod == "MCMC" && is.null(nMCMCsample)) {
    stop("nMCMCsample must be specified for MCMC method")
  }

  if (prob == "predictive" && (missing(m1) || missing(m2))) {
    stop("m1 and m2 must be specified for predictive probability")
  }

  # Set seed number for reproducible results
  set.seed(seed)

  # Generate random numbers of outcomes for group 1 in PoC study
  y1 <- matrix(rnorm(nsim * n1, mu1, sigma1), nrow = nsim)

  # Calculate sample mean for group 1
  bar.y1 <- rowSums(y1) / n1

  # Calculate sample standard deviation for group 1
  s1 <- sqrt(rowSums((y1 - bar.y1) ^ 2) / (n1 - 1))

  if(design == 'controlled' | design == 'external') {
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
      prob, design, prior, CalcMethod, theta0[i], nMC, nMCMCsample, n1, n2, m1, m2,
      kappa01, kappa02, nu01, nu02, mu01, mu02, sigma01, sigma02,
      bar.y1, bar.y2, s1, s2, r, ne1, ne2, alpha01, alpha02
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
  if(design == 'uncontrolled') {
    # For uncontrolled design, only include mu1
    results <- data.frame(mu1, Go = Go.and.NoGo.probs[1], NoGo = Go.and.NoGo.probs[2], Gray = Gray.prob)
  } else {
    # For controlled and external designs, include both mu1 and mu2
    results <- data.frame(mu1, mu2, Go = Go.and.NoGo.probs[1], NoGo = Go.and.NoGo.probs[2], Gray = Gray.prob)
  }

  return(results)
}
