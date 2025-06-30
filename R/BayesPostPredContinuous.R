#' Calculate Bayesian Posterior Probability or Bayesian Posterior Predictive Probability
#' for a Clinical Trial When Outcome is Continuous
#'
#' This function computes Bayesian posterior probability or posterior predictive probability
#' for continuous outcome clinical trials. The function supports controlled and uncontrolled designs
#' with Normal-Inverse-Chi-squared or vague priors, using three calculation methods: numerical
#' integration, Monte Carlo simulation, and Welch-Satterthwaite approximation.
#'
#' @param prob A character string specifying the type of probability to calculate
#'        (\code{prob = 'posterior'} or \code{prob = 'predictive'}).
#' @param design A character string specifying the type of trial design
#'        (\code{design = 'controlled'} or \code{design = 'uncontrolled'}).
#' @param prior A character string specifying the prior distribution
#'        (\code{prior = 'N-Inv-Chisq'} or \code{prior = 'vague'}).
#' @param CalcMethod A character string specifying the calculation method
#'        (\code{CalcMethod = 'NI'} for numerical integration, \code{CalcMethod = 'MC'} for Monte Carlo method,
#'        or \code{CalcMethod = 'WS'} for Welch-Satterthwaite approximation).
#' @param theta0 A numeric value representing the pre-specified threshold value.
#' @param nMC A positive integer representing the number of iterations for Monte Carlo simulation
#'        (required only if \code{CalcMethod = 'MC'}).
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
#' @param bar.y1 A numeric value representing the sample mean of group 1.
#' @param bar.y2 A numeric value representing the sample mean of group 2.
#' @param s1 A positive numeric value representing the sample standard deviation of group 1.
#' @param s2 A positive numeric value representing the sample standard deviation of group 2.
#' @param r A positive numeric value representing the parameter value associated with the distribution
#'        of mean for group 2 when \code{design = 'uncontrolled'}.
#'
#' @return A numeric vector representing the Bayesian posterior probability or Bayesian posterior
#'         predictive probability. The function can handle vectorized inputs.
#'
#' @details
#' The function can obtain:
#' \itemize{
#'   \item Bayesian posterior probability
#'   \item Bayesian posterior predictive probability
#' }
#'
#' Prior distribution of mean and variance of outcomes for each treatment group (k=1,2) can be either
#' (1) Normal-Inverse-Chi-squared or (2) Vague. The posterior distribution or posterior predictive
#' distribution of outcome for each treatment group follows a t-distribution.
#'
#' Three calculation methods are available:
#' \itemize{
#'   \item NI: Numerical integration method for exact computation
#'   \item MC: Monte Carlo simulation for flexible approximation
#'   \item WS: Welch-Satterthwaite approximation for computational efficiency
#' }
#'
#' @examples
#' # Calculate posterior probability using Normal-Inverse-Chi-squared prior with numerical integration
#' BayesPostPredContinuous(
#'   prob = 'posterior', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'NI',
#'   nMC = NULL, theta0 = 2, n1 = 12, n2 = 12, m1 = NULL, m2 = NULL,
#'   kappa01 = 5, kappa02 = 5, nu01 = 5, nu02 = 5, mu01 = 5, mu02 = 5,
#'   sigma01 = sqrt(5), sigma02 = sqrt(5), bar.y1 = 2, bar.y2 = 0, s1 = 1, s2 = 1, r = NULL
#' )
#'
#' # Calculate posterior predictive probability for uncontrolled design using vague prior
#' BayesPostPredContinuous(
#'   prob = 'predictive', design = 'uncontrolled', prior = 'vague', CalcMethod = 'NI',
#'   nMC = NULL, theta0 = 0.5, n1 = 12, n2 = NULL, m1 = 120, m2 = 120,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL, mu01 = NULL, mu02 = 0,
#'   sigma01 = NULL, sigma02 = NULL, bar.y1 = 2, bar.y2 = NULL, s1 = 1, s2 = NULL, r = 12
#' )
#'
#' @export
BayesPostPredContinuous <- function(prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
                                    theta0, nMC = 1e+4, n1, n2, m1, m2,
                                    kappa01, kappa02, nu01, nu02, mu01, mu02, sigma01, sigma02,
                                    bar.y1, bar.y2, s1, s2, r) {
  # Check parameter sets for uncontrolled design
  if((design == 'uncontrolled') & (is.null(r))) {
    stop('If you consider uncontrolled design, r should be non-null')
  }

  # Check parameter sets for posterior predictive probability
  if((prob == 'predictive') & (sum(sapply(list(m1, m2), is.null)) > 0)) {
    stop('If you calculate posterior predictive probability, m1 and m2 should be non-null')
  }

  # Check parameter sets for Normal-Inverse-Chi-squared prior
  if((prior == 'N-Inv-Chisq') & (sum(sapply(list(kappa01, nu01, sigma01), is.null)) > 0)) {
    stop('If you use the N-Inv-Chisq prior, kappa01, nu01 and sigma01 should be non-null')
  }

  # Define parameters for calculating posterior/posterior predictive probabilities
  if(prior == 'N-Inv-Chisq') {
    # Calculate updated precision parameters
    kappa.n1 <- kappa01 + n1
    kappa.n2 <- kappa02 + n2

    # Calculate updated degrees of freedom
    nu.t1 <- nu01 + n1
    if(design == 'controlled') {
      nu.t2 <- nu02 + n2
    } else if(design == 'uncontrolled') {
      nu.t2 <- nu.t1
    }

    # Calculate posterior means of t-distributions
    mu.t1 <- (kappa01 * mu01 + n1 * bar.y1) / kappa.n1
    if(design == 'controlled') {
      mu.t2 <- (kappa02 * mu02 + n2 * bar.y2) / kappa.n2
    } else if(design == 'uncontrolled') {
      mu.t2 <- mu02
    }

    # Calculate posterior variance for group 1
    var.n1 <- (nu01 * sigma01 ^ 2 + (n1 - 1) * s1 ^ 2 + n1 * kappa01 / (kappa01 + n1) * (mu01 - bar.y1) ^ 2) / nu.t1

    # Calculate posterior variance for group 2 (controlled design only)
    if(design == 'controlled') {
      var.n2 <- (nu02 * sigma02 ^ 2 + (n2 - 1) * s2 ^ 2 + n2 * kappa02 / (kappa02 + n2) * (mu02 - bar.y2) ^ 2) / nu.t2
    } else if(design == 'uncontrolled') {
      var.n2 <- NULL
    }

    # Calculate standard deviations of t-distributions based on probability type
    if(prob == 'posterior') {
      sd.t1 <- sqrt(var.n1 / kappa.n1)
      if(design == 'controlled') {
        sd.t2 <- sqrt(var.n2 / kappa.n2)
      } else if(design == 'uncontrolled') {
        sd.t2 <- sqrt(r) * sd.t1
      }
    } else if(prob == 'predictive') {
      sd.t1 <- sqrt((1 + kappa.n1) * var.n1 / (kappa.n1 * m1))
      if(design == 'controlled') {
        sd.t2 <- sqrt((1 + kappa.n2) * var.n2 / (kappa.n2 * m2))
      } else if(design == 'uncontrolled') {
        sd.t2 <- sqrt(r) * sd.t1
      }
    }
  } else if(prior == 'vague') {
    # Calculate degrees of freedom for vague priors
    nu.t1 <- n1 - 1
    if(design == 'controlled') {
      nu.t2 <- n2 - 1
    } else if(design == 'uncontrolled') {
      nu.t2 <- nu.t1
    }

    # Set means of t-distributions to sample means
    mu.t1 <- bar.y1
    if(design == 'controlled') {
      mu.t2 <- bar.y2
    } else if(design == 'uncontrolled') {
      mu.t2 <- mu02
    }

    # Calculate standard deviations of t-distributions based on probability type
    if(prob == 'posterior') {
      sd.t1 <- sqrt(s1 ^ 2 / n1)
      if(design == 'controlled') {
        sd.t2 <- sqrt(s2 ^ 2 / n2)
      } else if(design == 'uncontrolled') {
        sd.t2 <- sqrt(r) * sd.t1
      }
    } else if(prob == 'predictive') {
      sd.t1 <- sqrt((1 + n1) * s1 ^ 2 / (n1 * m1))
      if(design == 'controlled') {
        sd.t2 <- sqrt((1 + n2) * s2 ^ 2 / (n2 * m2))
      } else if(design == 'uncontrolled') {
        sd.t2 <- sqrt(r) * sd.t1
      }
    }
  }

  # Calculate the probability of exceeding θ₀ using the specified method
  if(CalcMethod == 'NI') {
    results <- pNIdifft(theta0, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2)
  } else if(CalcMethod == 'MC') {
    results <- pMCdifft(nMC, theta0, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2)
  } else if(CalcMethod == 'WS') {
    results <- pWSdifft(theta0, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2)
  }

  # Return results
  return(results)
}
