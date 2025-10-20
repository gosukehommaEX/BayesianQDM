#' Calculate Bayesian Posterior Probability or Bayesian Posterior Predictive Probability for a Clinical Trial When Outcome is Continuous
#'
#' @description
#' This function computes Bayesian posterior probability or posterior predictive probability
#' for continuous outcome clinical trials. The function supports controlled, uncontrolled, and
#' external control designs with Normal-Inverse-Chi-squared or vague priors, using three calculation
#' methods: numerical integration, Monte Carlo simulation, and Welch-Satterthwaite approximation.
#' For external control designs, power priors are incorporated using exact conjugate representation
#' as Normal-Inverse-Chi-squared distributions, enabling closed-form computation without MCMC sampling.
#'
#' @param prob A character string specifying the type of probability to use
#'        (\code{prob = 'posterior'} or \code{prob = 'predictive'}).
#' @param design A character string specifying the type of trial design
#'        (\code{design = 'controlled'}, \code{design = 'uncontrolled'}, or \code{design = 'external'}).
#' @param prior A character string specifying the prior distribution
#'        (\code{prior = 'N-Inv-Chisq'} or \code{prior = 'vague'}).
#' @param CalcMethod A character string specifying the calculation method
#'        (\code{CalcMethod = 'NI'} for numerical integration, \code{CalcMethod = 'MC'} for Monte Carlo method,
#'        or \code{CalcMethod = 'WS'} for Welch-Satterthwaite approximation).
#' @param theta0 A numeric value representing the pre-specified threshold value.
#' @param nMC A positive integer representing the number of iterations for Monte Carlo simulation
#'        (required if \code{CalcMethod = 'MC'}).
#' @param n1 A positive integer representing the number of patients in group 1 for the PoC trial.
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
#'        mean of group 2 for \code{design = 'uncontrolled'}.
#' @param ne1 A positive integer representing the sample size for group 1 in external trial (can be NULL if no external treatment data).
#' @param ne2 A positive integer representing the sample size for group 2 in external trial (can be NULL if no external control data).
#' @param alpha01 A positive numeric value representing the power prior scale parameter for group 1 (can be NULL if no external treatment data).
#' @param alpha02 A positive numeric value representing the power prior scale parameter for group 2 (can be NULL if no external control data).
#' @param bar.ye1 A numeric value representing the external sample mean of group 1 (required if external treatment data available).
#' @param bar.ye2 A numeric value representing the external sample mean of group 2 (required if external control data available).
#' @param se1 A positive numeric value representing the external sample standard deviation of group 1 (required if external treatment data available).
#' @param se2 A positive numeric value representing the external sample standard deviation of group 2 (required if external control data available).
#' @param lower.tail logical; if TRUE (default), probabilities are P(theta <= theta0), otherwise, P(theta > theta0)
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
#' For controlled and uncontrolled designs, three calculation methods are available:
#' \itemize{
#'   \item NI: Numerical integration method for exact computation
#'   \item MC: Monte Carlo simulation for flexible approximation
#'   \item WS: Welch-Satterthwaite approximation for computational efficiency
#' }
#'
#' For external control designs, power priors are incorporated using exact conjugate representation:
#' \itemize{
#'   \item Power priors for normal data are mathematically equivalent to Normal-Inverse-Chi-squared distributions
#'   \item This enables closed-form posterior computation without MCMC sampling
#'   \item Alpha parameters control the degree of borrowing (0 = no borrowing, 1 = full borrowing)
#'   \item The method preserves complete Bayesian rigor with no approximation
#' }
#'
#' The external design supports:
#' \itemize{
#'   \item External control data only (ne2, alpha02, bar.ye2, se2)
#'   \item External treatment data only (ne1, alpha01, bar.ye1, se1)
#'   \item Both external control and treatment data
#' }
#'
#' @examples
#' # Example 1: Numerical Integration (NI) method with N-Inv-Chisq prior
#' pPostPred1Continuous(
#'   prob = 'posterior', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'NI',
#'   theta0 = 2, n1 = 12, n2 = 12, kappa01 = 5, kappa02 = 5, nu01 = 5, nu02 = 5,
#'   mu01 = 5, mu02 = 5, sigma01 = sqrt(5), sigma02 = sqrt(5),
#'   bar.y1 = 2, bar.y2 = 0, s1 = 1, s2 = 1, lower.tail = FALSE
#' )
#'
#' # Example 2: Monte Carlo (MC) method with vague prior
#' pPostPred1Continuous(
#'   prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'MC',
#'   theta0 = 1, nMC = 10000, n1 = 12, n2 = 12,
#'   bar.y1 = 3, bar.y2 = 1, s1 = 1.5, s2 = 1.2, lower.tail = FALSE
#' )
#'
#' # Example 3: Welch-Satterthwaite (WS) approximation with N-Inv-Chisq prior
#' pPostPred1Continuous(
#'   prob = 'predictive', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'WS',
#'   theta0 = 0.5, n1 = 15, n2 = 15, m1 = 100, m2 = 100,
#'   kappa01 = 3, kappa02 = 3, nu01 = 4, nu02 = 4, mu01 = 2, mu02 = 2,
#'   sigma01 = 2, sigma02 = 2, bar.y1 = 2.5, bar.y2 = 1.8, s1 = 1.8, s2 = 1.6, lower.tail = FALSE
#' )
#'
#' # Example 4: External control design with power prior (NI method)
#' pPostPred1Continuous(
#'   prob = 'posterior', design = 'external', prior = 'vague', CalcMethod = 'NI',
#'   theta0 = 1.5, n1 = 12, n2 = 12, bar.y1 = 4, bar.y2 = 2, s1 = 1.2, s2 = 1.1,
#'   ne2 = 20, alpha02 = 0.5, bar.ye2 = 1.8, se2 = 1.0, lower.tail = FALSE
#' )
#'
#' # Example 5: External design with both treatment and control data
#' pPostPred1Continuous(
#'   prob = 'posterior', design = 'external', prior = 'N-Inv-Chisq', CalcMethod = 'WS',
#'   theta0 = 1.0, n1 = 15, n2 = 15, bar.y1 = 3.5, bar.y2 = 2.0, s1 = 1.3, s2 = 1.1,
#'   kappa01 = 2, kappa02 = 2, nu01 = 3, nu02 = 3, mu01 = 3, mu02 = 2,
#'   sigma01 = 1.5, sigma02 = 1.5,
#'   ne1 = 25, ne2 = 25, alpha01 = 0.7, alpha02 = 0.7,
#'   bar.ye1 = 3.2, bar.ye2 = 1.9, se1 = 1.4, se2 = 1.2, lower.tail = FALSE
#' )
#'
#' @export
pPostPred1Continuous <- function(prob = "posterior", design = "controlled", prior = "vague", CalcMethod = "NI",
                                 theta0, nMC = NULL, n1, n2, m1 = NULL, m2 = NULL, kappa01 = NULL, kappa02 = NULL,
                                 nu01 = NULL, nu02 = NULL, mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
                                 bar.y1, bar.y2, s1, s2, r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
                                 bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL, lower.tail = TRUE) {

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

  if (!CalcMethod %in% c("NI", "MC", "WS")) {
    stop("CalcMethod must be 'NI', 'MC', or 'WS'")
  }

  # Validate required parameters for each method
  if (CalcMethod == "MC" && is.null(nMC)) {
    stop("nMC must be specified for Monte Carlo method")
  }

  if (prob == "predictive" && (is.null(m1) || is.null(m2))) {
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

  # Validate N-Inv-Chisq prior parameters
  if (!is.null(prior) && prior == 'N-Inv-Chisq' && design != 'external') {
    required_params <- c("kappa01", "kappa02", "nu01", "nu02", "mu01", "mu02", "sigma01", "sigma02")
    missing_params <- required_params[sapply(required_params, function(x) is.null(get(x, envir = environment())))]
    if (length(missing_params) > 0) {
      stop(paste("For N-Inv-Chisq prior, the following parameters are required:", paste(missing_params, collapse = ", ")))
    }
  }

  # Calculate hyperparameters for different designs
  if (design == "external") {
    # Power prior implementation using exact conjugate representation
    if (!is.null(prior) && prior == 'N-Inv-Chisq') {
      ## Informative prior case
      # Group 1 (Treatment) - Apply power prior if external data available
      if (!is.null(ne1) && !is.null(alpha01)) {
        # Power prior parameters for group 1
        mu.n1 <- (alpha01 * ne1 * bar.ye1 + kappa01 * mu01) / (alpha01 * ne1 + kappa01)
        # Posterior parameters after current data
        mu.t1 <- '/'(
          alpha01 * ne1 * bar.ye1 + kappa01 * mu01 + n1 * bar.y1,
          alpha01 * ne1 + kappa01 + n1
        )
        kappa.star.n1 <- alpha01 * ne1 + kappa01 + n1
        nu.t1 <- alpha01 * ne1 + nu01 + n1
        sigma2.star.n1 <- '/'(
          '+'(
            '+'(
              alpha01 * (ne1 - 1) * se1 ^ 2 + nu01 * sigma01 ^ 2,
              (alpha01 * ne1 * kappa01 * (bar.ye1 - mu01) ^ 2) / (alpha01 * ne1 + kappa01)
            ),
            '+'(
              (n1 - 1) * s1 ^ 2,
              n1 * (alpha01 * ne1 + kappa01) / (alpha01 * ne1 + kappa01 + n1) * (mu.n1 - bar.y1) ^ 2
            )
          ),
          nu.t1
        )
      } else {
        # No external treatment data - use original prior
        mu.t1 <- (kappa01 * mu01 + n1 * bar.y1) / (kappa01 + n1)
        kappa.star.n1 <- kappa01 + n1
        nu.t1 <- nu01 + n1
        sigma2.star.n1 <- '+'(
          nu01 * sigma01 ^ 2 + (n1 - 1) * s1 ^ 2,
          n1 * kappa01 * (mu01 - bar.y1) ^ 2 / (kappa01 + n1)
        ) / nu.t1
      }
      # Group 2 (Control) - Apply power prior if external data available
      if (!is.null(ne2) && !is.null(alpha02)) {
        # Power prior parameters for group 1
        mu.n2 <- (alpha02 * ne2 * bar.ye2 + kappa02 * mu02) / (alpha02 * ne2 + kappa02)
        # Posterior parameters after current data
        mu.t2 <- '/'(
          alpha02 * ne2 * bar.ye2 + kappa02 * mu02 + n2 * bar.y2,
          alpha02 * ne2 + kappa02 + n2
        )
        kappa.star.n2 <- alpha02 * ne2 + kappa02 + n2
        nu.t2 <- alpha02 * ne2 + nu02 + n2
        sigma2.star.n2 <- '/'(
          '+'(
            '+'(
              alpha02 * (ne2 - 2) * se2 ^ 2 + nu02 * sigma02 ^ 2,
              (alpha02 * ne2 * kappa02 * (bar.ye2 - mu02) ^ 2) / (alpha02 * ne2 + kappa02)
            ),
            '+'(
              (n2 - 2) * s2 ^ 2,
              n2 * (alpha02 * ne2 + kappa02) / (alpha02 * ne2 + kappa02 + n2) * (mu.n2 - bar.y2) ^ 2
            )
          ),
          nu.t2
        )
      } else {
        # No external control data - use original prior
        mu.t2 <- (kappa02 * mu02 + n2 * bar.y2) / (kappa02 + n2)
        kappa.star.n2 <- kappa02 + n2
        nu.t2 <- nu02 + n2
        Sy2 <-
          sigma2.star.n2 <- '+'(
            nu02 * sigma02 ^ 2 + (n2 - 1) * s2 ^ 2,
            n2 * kappa02 * (mu02 - bar.y2) ^ 2 / (kappa02 + n2)
          ) / nu.t2
      }
    } else {
      ## Vague prior case
      # Group 1 (Treatment) - Apply power prior if external data available
      if (!is.null(ne1) && !is.null(alpha01)) {
        # Posterior parameters with power prior
        mu.t1 <- (alpha01 * ne1 * bar.ye1 + n1 * bar.y1) / (alpha01 * ne1 + n1)
        kappa.star.n1 <- alpha01 * ne1 + n1
        nu.t1 <- alpha01 * ne1 + n1 - 1
        sigma2.star.n1 <- '/'(
          '+'(
            alpha01 * (ne1 - 1) * se1 ^ 2 + (n1 - 1) * s1 ^ 2,
            (alpha01 * ne1 * n1 * (bar.ye1 - bar.y1) ^ 2) / (alpha01 * ne1 + n1)
          ),
          alpha01 * ne1 + n1
        )
      } else {
        # No external treatment data - use vague prior
        mu.t1 <- bar.y1
        kappa.star.n1 <- n1
        nu.t1 <- n1 - 1
        sigma2.star.n1 <- s1 ^ 2
      }
      # Group 2 (Control) - Apply power prior if external data available
      if (!is.null(ne2) && !is.null(alpha02)) {
        # Posterior parameters with power prior
        mu.t2 <- (alpha02 * ne2 * bar.ye2 + n2 * bar.y2) / (alpha02 * ne2 + n2)
        kappa.star.n2 <- alpha02 * ne2 + n2
        nu.t2 <- alpha02 * ne2 + n2 - 2
        sigma2.star.n2 <- '/'(
          '+'(
            alpha02 * (ne2 - 2) * se2 ^ 2 + (n2 - 2) * s2 ^ 2,
            (alpha02 * ne2 * n2 * (bar.ye2 - bar.y2) ^ 2) / (alpha02 * ne2 + n2)
          ),
          alpha02 * ne2 + n2
        )
      } else {
        # No external control data - use vague prior
        mu.t2 <- bar.y2
        kappa.star.n2 <- n2
        nu.t2 <- n2 - 1
        sigma2.star.n2 <- s2 ^ 2
      }
    }
    # Calculate standard deviations of t-distributions based on probability type
    if (prob == 'posterior') {
      sd.t1 <- sqrt(sigma2.star.n1 / kappa.star.n1)
      sd.t2 <- sqrt(sigma2.star.n2 / kappa.star.n2)
    } else if (prob == 'predictive') {
      sd.t1 <- sqrt((1 + 1 / kappa.star.n1) * sigma2.star.n1 / m1)
      sd.t2 <- sqrt((1 + 1 / kappa.star.n2) * sigma2.star.n2 / m2)
    }
  } else {
    # For controlled and uncontrolled designs
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
      var.n1 <- '+'(
        nu01 * sigma01 ^ 2 + (n1 - 1) * s1 ^ 2,
        n1 * kappa01 * (mu01 - bar.y1) ^ 2 / (kappa01 + n1)
      ) / nu.t1
      # Calculate posterior variance for group 2 (controlled design only)
      if(design == 'controlled') {
        var.n2 <- '+'(
          nu02 * sigma02 ^ 2 + (n2 - 1) * s2 ^ 2,
          n2 * kappa02 * (mu02 - bar.y2) ^ 2 / (kappa02 + n2)
        )
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
        sd.t1 <- sqrt((1 + 1 / kappa.n1) * var.n1 / m1)
        if(design == 'controlled') {
          sd.t2 <- sqrt((1 + 1 / kappa.n2) * var.n2 / m2)
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
        sd.t1 <- sqrt((1 + 1 / n1) * s1 ^ 2 / m1)
        if(design == 'controlled') {
          sd.t2 <- sqrt((1 + 1 / n2) * s2 ^ 2 / m2)
        } else if(design == 'uncontrolled') {
          sd.t2 <- sqrt(r) * sd.t1
        }
      }
    }
  }

  # Calculate the probability of below or exceeding θ₀ using the specified method
  if(CalcMethod == 'NI') {
    results <- pNIdifft(theta0, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2, lower.tail)
  } else if(CalcMethod == 'MC') {
    results <- pMCdifft(nMC, theta0, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2, lower.tail)
  } else if(CalcMethod == 'WS') {
    results <- pWSdifft(theta0, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2, lower.tail)
  } else {
    stop("Invalid CalcMethod.")
  }

  # Return results
  return(results)
}
