#' Calculate Bayesian Posterior Probability or Bayesian Posterior Predictive Probability for a Clinical Trial When Outcome is Continuous
#'
#' This function computes Bayesian posterior probability or posterior predictive probability
#' for continuous outcome clinical trials. The function supports controlled, uncontrolled, and
#' external control designs with Normal-Inverse-Chi-squared or vague priors, using three calculation
#' methods: numerical integration, Monte Carlo simulation, and Welch-Satterthwaite approximation.
#' For external control designs, MCMC sampling is used to incorporate historical data through power priors.
#'
#' @param prob A character string specifying the type of probability to calculate.
#' Options are 'posterior' (default) for posterior probability or 'predictive' for posterior predictive probability.
#' @param design A character string specifying the trial design.
#' Options are 'controlled' (default), 'uncontrolled', or 'external'.
#' @param prior A character string specifying the prior distribution.
#' Options are 'vague' (default) or 'N-Inv-Chisq' for Normal-Inverse-Chi-squared.
#' @param CalcMethod A character string specifying the calculation method.
#' Options are 'NI' (numerical integration, default), 'MC' (Monte Carlo), 'WS' (Welch-Satterthwaite),
#' or 'MCMC' (MCMC sampling for external design).
#' @param theta0 A numeric value representing the threshold for the treatment effect.
#' @param nMC A positive integer representing the number of Monte Carlo iterations for MC method (default: NULL).
#' @param nMCMCsample A positive integer representing the number of MCMC iterations for external design (default: NULL).
#' @param n1 A positive integer representing the sample size for group 1 in PoC trial.
#' @param n2 A positive integer representing the sample size for group 2 in PoC trial.
#' @param m1 A positive integer representing the sample size for group 1 in future trial (for predictive probability).
#' @param m2 A positive integer representing the sample size for group 2 in future trial (for predictive probability).
#' @param kappa01 A positive numeric value representing the prior precision parameter for group 1 (N-Inv-Chisq prior).
#' @param kappa02 A positive numeric value representing the prior precision parameter for group 2 (N-Inv-Chisq prior).
#' @param nu01 A positive numeric value representing the prior degrees of freedom for group 1 (N-Inv-Chisq prior).
#' @param nu02 A positive numeric value representing the prior degrees of freedom for group 2 (N-Inv-Chisq prior).
#' @param mu01 A numeric value representing the prior mean for group 1 (N-Inv-Chisq prior).
#' @param mu02 A numeric value representing the prior mean for group 2 (N-Inv-Chisq prior).
#' @param sigma01 A positive numeric value representing the prior standard deviation for group 1 (N-Inv-Chisq prior).
#' @param sigma02 A positive numeric value representing the prior standard deviation for group 2 (N-Inv-Chisq prior).
#' @param bar.y1 A numeric value representing the sample mean of group 1.
#' @param bar.y2 A numeric value representing the sample mean of group 2.
#' @param s1 A positive numeric value representing the sample standard deviation of group 1.
#' @param s2 A positive numeric value representing the sample standard deviation of group 2.
#' @param r A positive numeric value for uncontrolled design (default: NULL).
#' @param ne1 A positive integer representing the sample size for group 1 in external trial (default: NULL).
#' @param ne2 A positive integer representing the sample size for group 2 in external trial (default: NULL).
#' @param alpha01 A positive numeric value representing the power prior scale parameter for group 1 (default: NULL).
#' @param alpha02 A positive numeric value representing the power prior scale parameter for group 2 (default: NULL).
#'
#' @return A numeric vector representing the Bayesian posterior probability or Bayesian posterior
#' predictive probability. The function can handle vectorized inputs.
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
#' For controlled and uncontrolled designs, three calculation methods are available:
#' \itemize{
#'   \item NI: Numerical integration method for exact computation
#'   \item MC: Monte Carlo simulation for flexible approximation
#'   \item WS: Welch-Satterthwaite approximation for computational efficiency
#' }
#'
#' For external control designs, the function uses MCMC sampling to incorporate historical
#' data through power prior methodology:
#' \itemize{
#'   \item MCMC: Markov Chain Monte Carlo sampling for posterior inference with external data
#'   \item Power priors allow controlled borrowing from historical data
#'   \item alpha parameters control the degree of borrowing (0 = no borrowing, 1 = full borrowing)
#' }
#'
#' The external design supports:
#' \itemize{
#'   \item External control data only (ne2, alpha02)
#'   \item External treatment data only (ne1, alpha01)
#'   \item Both external control and treatment data
#' }
#'
#' @examples
#' # Example 1: Numerical Integration (NI) method with N-Inv-Chisq prior
#' BayesPostPredContinuous(
#'   prob = 'posterior', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'NI',
#'   theta0 = 2, n1 = 12, n2 = 12, kappa01 = 5, kappa02 = 5, nu01 = 5, nu02 = 5,
#'   mu01 = 5, mu02 = 5, sigma01 = sqrt(5), sigma02 = sqrt(5),
#'   bar.y1 = 2, bar.y2 = 0, s1 = 1, s2 = 1
#' )
#'
#' # Example 2: Monte Carlo (MC) method with vague prior
#' BayesPostPredContinuous(
#'   prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'MC',
#'   theta0 = 1, nMC = 10000, n1 = 12, n2 = 12,
#'   bar.y1 = 3, bar.y2 = 1, s1 = 1.5, s2 = 1.2
#' )
#'
#' # Example 3: Welch-Satterthwaite (WS) approximation with N-Inv-Chisq prior
#' BayesPostPredContinuous(
#'   prob = 'predictive', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'WS',
#'   theta0 = 0.5, n1 = 15, n2 = 15, m1 = 100, m2 = 100,
#'   kappa01 = 3, kappa02 = 3, nu01 = 4, nu02 = 4, mu01 = 2, mu02 = 2,
#'   sigma01 = 2, sigma02 = 2, bar.y1 = 2.5, bar.y2 = 1.8, s1 = 1.8, s2 = 1.6
#' )
#'
#' # Example 2: External control design with MCMC method
#' BayesPostPredContinuous(
#'   prob = 'posterior', design = 'external', CalcMethod = 'MCMC',
#'   theta0 = 1.5, nMCMCsample = 5000, n1 = 12, n2 = 12,
#'   bar.y1 = 4, bar.y2 = 2, s1 = 1.2, s2 = 1.1,
#'   ne2 = 20, alpha02 = 0.5
#' )
#'
#' @export
BayesPostPredContinuous <- function(prob = "posterior", design = "controlled", prior = "vague", CalcMethod = "NI",
                                    theta0, nMC = NULL, nMCMCsample = NULL,
                                    n1, n2, m1, m2, kappa01, kappa02, nu01, nu02,
                                    mu01, mu02, sigma01, sigma02, bar.y1, bar.y2, s1, s2,
                                    r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL) {

  # Input validation
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

  # Define parameters for calculating posterior/posterior predictive probabilities
  if(!is.null(prior) && prior == 'N-Inv-Chisq') {
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
  } else if(!is.null(prior) && prior == 'vague') {
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
  } else if((design == 'external') && (CalcMethod == 'MCMC')) {
    # Use MCMC sampling for external design with power prior
    results <- pMCMCdiff(nMCMCsample, theta0, bar.y1, bar.y2, s1, s2, n1, n2, ne1, ne2, alpha01, alpha02)
  } else {
    stop("Invalid combination of design and CalcMethod.")
  }

  # Return results
  return(results)
}
