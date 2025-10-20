#' Calculate the Go, NoGo and Gray Probabilities for a Clinical Trial When Outcome is Continuous Under the Bayesian Framework Using Two Metrics
#'
#' @description
#' This function computes the Go, NoGo, and Gray probabilities for continuous outcome clinical trials
#' using the Bayesian framework. The function supports controlled, uncontrolled, and external control designs with
#' multiple calculation methods including numerical integration, Monte Carlo simulation, and Welch-Satterthwaite approximation.
#' For external control designs, power priors are incorporated using exact conjugate representation.
#'
#' @param nsim A positive integer representing the number of iterations for calculating posterior/posterior predictive probability.
#' @param prob A character string specifying the type of probability to use
#'        (\code{prob = 'posterior'} or \code{prob = 'predictive'}).
#' @param design A character string specifying the type of trial design
#'        (\code{design = 'controlled'}, \code{design = 'uncontrolled'}, or \code{design = 'external'}).
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
#'        Go probability when \code{prob = 'predictive'}.
#' @param nMC A positive integer representing the number of iterations for Monte Carlo simulation
#'        (required if \code{CalcMethod = 'MC'}).
#' @param gamma1 A numeric value between 0 and 1 representing the minimum probability to declare success.
#' @param gamma2 A numeric value between 0 and 1 representing the futility threshold.
#' @param n1 A positive integer representing the number of patients in group 1 for the proof-of-concept (PoC) trial.
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
#'        mean of group 2 for \code{design = 'uncontrolled'}.
#' @param ne1 A positive integer representing the sample size for group 1 in external trial
#'        (required if external design, can be NULL if no external treatment data).
#' @param ne2 A positive integer representing the sample size for group 2 in external trial
#'        (required if external design, can be NULL if no external control data).
#' @param alpha01 A positive numeric value representing the scale parameter (power prior) for group 1
#'        (required if external design, can be NULL if no external treatment data).
#' @param alpha02 A positive numeric value representing the scale parameter (power prior) for group 2
#'        (required if external design, can be NULL if no external control data).
#' @param bar.ye1 A numeric value representing the external sample mean of group 1 (required if external treatment data available).
#' @param bar.ye2 A numeric value representing the external sample mean of group 2 (required if external control data available).
#' @param se1 A positive numeric value representing the external sample standard deviation of group 1 (required if external treatment data available).
#' @param se2 A positive numeric value representing the external sample standard deviation of group 2 (required if external control data available).
#' @param Gray_inc_Miss A logical value representing that Miss probability is included into Gray probability if TRUE, not otherwise.
#' @param seed A numeric value representing the seed number for reproducible random number generation.
#'
#' @return A data frame containing the true means for both groups and the Go, NoGo, and Gray probabilities.
#'
#' @details
#' The function performs Monte Carlo simulation to evaluate operating characteristics by:
#' \itemize{
#'   \item Generating random trial data based on specified true parameters
#'   \item Computing posterior or predictive probabilities for each simulated trial
#'   \item Classifying each trial as Go, NoGo, or Gray based on decision thresholds
#' }
#'
#' For external control designs, power priors are incorporated using exact conjugate representation:
#' \itemize{
#'   \item Power priors for normal data are mathematically equivalent to Normal-Inverse-Chi-squared distributions
#'   \item This enables closed-form computation without MCMC sampling
#'   \item Alpha parameters control the degree of borrowing (0 = no borrowing, 1 = full borrowing)
#' }
#'
#' Decision rules:
#' \itemize{
#'   \item **Go**: P(treatment effect > threshold) ≥ γ₁
#'   \item **NoGo**: P(treatment effect > threshold) ≤ γ₂
#'   \item **Gray**: γ₂ < P(treatment effect > threshold) < γ₁
#' }
#'
#' @examples
#' # Example 1: Controlled design with vague prior and NI method
#' pGoNogoGray1Continuous(
#'   nsim = 100, prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
#'   theta.TV = 1.5, theta.MAV = -0.5, theta.NULL = NULL,
#'   nMC = NULL, gamma1 = 0.7, gamma2 = 0.2,
#'   n1 = 15, n2 = 15, m1 = NULL, m2 = NULL,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
#'   mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
#'   mu1 = 3, mu2 = 1, sigma1 = 1.2, sigma2 = 1.1,
#'   r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
#'   bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL, seed = 2
#' )
#'
#' # Example 2: External design with control data
#' \dontrun{
#' pGoNogoGray1Continuous(
#'   nsim = 100, prob = 'posterior', design = 'external', prior = 'vague', CalcMethod = 'WS',
#'   theta.TV = 1.0, theta.MAV = 0.0, theta.NULL = NULL,
#'   nMC = NULL, gamma1 = 0.8, gamma2 = 0.2,
#'   n1 = 12, n2 = 12, m1 = NULL, m2 = NULL,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
#'   mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
#'   mu1 = 2, mu2 = 0, sigma1 = 1, sigma2 = 1,
#'   r = NULL, ne1 = NULL, ne2 = 20, alpha01 = NULL, alpha02 = 0.5,
#'   bar.ye1 = NULL, bar.ye2 = 0, se1 = NULL, se2 = 1, seed = 4
#' )
#' }
#'
#' # Example 3: Controlled design with predictive probability
#' pGoNogoGray1Continuous(
#'   nsim = 100, prob = 'predictive', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'NI',
#'   theta.TV = NULL, theta.MAV = NULL, theta.NULL = 2.0,
#'   nMC = NULL, gamma1 = 0.75, gamma2 = 0.15,
#'   n1 = 15, n2 = 15, m1 = 50, m2 = 50,
#'   kappa01 = 3, kappa02 = 3, nu01 = 4, nu02 = 4,
#'   mu01 = 3.5, mu02 = 1.5, sigma01 = 1.5, sigma02 = 1.5,
#'   mu1 = 3.2, mu2 = 1.3, sigma1 = 1.4, sigma2 = 1.2,
#'   r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
#'   bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL, seed = 3
#' )
#'
#' # Example 4: External design with predictive probability
#' \dontrun{
#' pGoNogoGray1Continuous(
#'   nsim = 100, prob = 'predictive', design = 'external', prior = 'vague', CalcMethod = 'MC',
#'   theta.TV = NULL, theta.MAV = NULL, theta.NULL = 1.5,
#'   nMC = 5000, gamma1 = 0.7, gamma2 = 0.2,
#'   n1 = 12, n2 = 12, m1 = 30, m2 = 30,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
#'   mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
#'   mu1 = 2.5, mu2 = 1.0, sigma1 = 1.3, sigma2 = 1.1,
#'   r = NULL, ne1 = 15, ne2 = 18, alpha01 = 0.6, alpha02 = 0.7,
#'   bar.ye1 = 2.3, bar.ye2 = 0.9, se1 = 1.2, se2 = 1.0, seed = 5
#' )
#' }
#'
#' @importFrom stats rnorm
#' @export
pGoNogoGray1Continuous <- function(nsim, prob, design, prior, CalcMethod, theta.TV, theta.MAV, theta.NULL,
                                   nMC = NULL, gamma1, gamma2, n1, n2, m1, m2, kappa01, kappa02, nu01, nu02,
                                   mu01, mu02, sigma01, sigma02, mu1, mu2, sigma1, sigma2,
                                   r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
                                   bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL, Gray_inc_Miss = FALSE, seed) {

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

  if (!CalcMethod %in% c("NI", "MC", "WS")) {
    stop("CalcMethod must be 'NI', 'MC', or 'WS'")
  }

  # Validate required parameters for each method
  if (CalcMethod == "MC" && is.null(nMC)) {
    stop("nMC must be specified for Monte Carlo method")
  }

  if (prob == "predictive" && (missing(m1) || missing(m2))) {
    stop("m1 and m2 must be specified for predictive probability")
  }

  # Validate external design parameters
  if (design == "external") {
    if (is.null(ne1) && is.null(ne2)) {
      stop("For external design, at least one of ne1 or ne2 must be specified")
    }
    if (!is.null(ne1) && (is.null(alpha01) || is.null(bar.ye1) || is.null(se1))) {
      stop("For external treatment data, alpha01, bar.ye1, and se1 must be specified")
    }
    if (!is.null(ne2) && (is.null(alpha02) || is.null(bar.ye2) || is.null(se2))) {
      stop("For external control data, alpha02, bar.ye2, and se2 must be specified")
    }
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
    # For uncontrolled design, set bar.y2 and s2 to NULL
    bar.y2 <- NULL
    s2 <- NULL
  }

  # External data are provided as fixed historical sample statistics
  # No data generation needed - these are already observed historical values

  # Initialize vectors to store Go and NoGo probabilities
  if(prob == 'posterior') {
    # Define thresholds: theta.TV for Go, theta.MAV for NoGo
    theta_values <- c(theta.TV, theta.MAV)
  } else {
    # For predictive probability, use theta.NULL (single threshold)
    theta_values <- c(theta.NULL, theta.NULL)
  }

  # Calculate posterior/posterior predictive probabilities for each threshold
  gPost <- lapply(seq_along(theta_values), function(i) {
    pPostPred1Continuous(
      prob = prob, design = design, prior = prior, CalcMethod = CalcMethod,
      theta0 = theta_values[i], nMC = nMC, n1 = n1, n2 = n2, m1 = m1, m2 = m2,
      kappa01 = kappa01, kappa02 = kappa02, nu01 = nu01, nu02 = nu02,
      mu01 = mu01, mu02 = mu02, sigma01 = sigma01, sigma02 = sigma02,
      bar.y1 = bar.y1, bar.y2 = bar.y2, s1 = s1, s2 = s2,
      r = r, ne1 = ne1, ne2 = ne2, alpha01 = alpha01, alpha02 = alpha02,
      bar.ye1 = bar.ye1, bar.ye2 = bar.ye2, se1 = se1, se2 = se2, lower.tail = c(FALSE, TRUE)[i]
    )
  })

  # Calculate Go, NoGo and Miss probabilities based on decision criteria
  GoNogoProb <- sapply(seq(3), function(j) {
    # Create indicator matrix for Go (j=1) or NoGo (j=2) decisions
    if(j == 1) {
      I <- as.numeric((gPost[[1]] >= gamma1) & (gPost[[2]] < gamma2))
    } else if(j == 2) {
      I <- as.numeric((gPost[[1]] < gamma1) & (gPost[[2]] >= gamma2))
    } else {
      I <- as.numeric((gPost[[1]] >= gamma1) & (gPost[[2]] >= gamma2))
    }
    sum(I) / nsim
  })

  # Check for positive Miss probabilities
  if(sum(GoNogoProb[3]) > 0) {
    stop('Because positive Miss probability(s) is obtained, re-consider appropriate threshold')
  }

  # Calculate Gray probability (complement of Go and NoGo)
  GrayProb <- 1 - sum(GoNogoProb[-3])
  if(Gray_inc_Miss) {
    GrayProb <- GrayProb + GoNogoProb[3]
  }

  # Create results data frame
  if(design == 'uncontrolled') {
    # For uncontrolled design, only include mu1
    results <- data.frame(
      mu1, Go = GoNogoProb[1], Gray = GrayProb, NoGo = GoNogoProb[2], Miss = GoNogoProb[3]
    )
  } else {
    # For controlled and external designs, include both mu1 and mu2
    results <- data.frame(
      mu1, mu2, Go = GoNogoProb[1], Gray = GrayProb, NoGo = GoNogoProb[2], Miss = GoNogoProb[3]
    )
  }

  return(results)
}
