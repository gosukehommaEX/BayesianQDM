#' Calculate Bayesian Posterior Probability or Posterior Predictive Probability
#' for a Clinical Trial with a Single Continuous Endpoint
#'
#' This function computes Bayesian posterior probability or posterior predictive
#' probability for continuous outcome clinical trials. The function supports
#' controlled, uncontrolled, and external control designs, with Normal-Inverse-Chi-squared
#' conjugate priors or vague priors. External data can be incorporated through
#' power priors.
#'
#' @param prob A character string specifying the type of probability to calculate.
#'        Options are \code{'posterior'} for posterior probability or \code{'predictive'}
#'        for posterior predictive probability.
#' @param design A character string specifying the type of trial design. Options are
#'        \code{'controlled'} for randomized controlled trials, \code{'uncontrolled'}
#'        for single-arm studies, or \code{'external'} for designs incorporating external data.
#' @param prior A character string specifying the prior distribution type. Options are
#'        \code{'vague'} for Jeffreys prior or \code{'N-Inv-Chisq'} for Normal-Inverse-Chi-squared
#'        conjugate prior.
#' @param CalcMethod A character string specifying the calculation method. Options are
#'        \code{'NI'} for numerical integration, \code{'MC'} for Monte Carlo simulation,
#'        or \code{'WS'} for Welch-Satterthwaite approximation.
#' @param theta0 A numeric value representing the pre-specified threshold value for
#'        the treatment effect (difference in means).
#' @param nMC A positive integer representing the number of Monte Carlo samples
#'        (required if \code{CalcMethod = 'MC'}, otherwise set to NULL).
#' @param n1 A positive integer representing the number of patients in group 1
#'        (treatment) for the proof-of-concept (PoC) trial.
#' @param n2 A positive integer representing the number of patients in group 2.
#'        For \code{design = 'controlled'} or \code{'external'}: sample size of the
#'        control group. For \code{design = 'uncontrolled'}: set equal to n1
#'        (not actually used in calculations, but required for consistency).
#' @param m1 A positive integer representing the number of patients in group 1 for
#'        the future trial (required if \code{prob = 'predictive'}, otherwise set to NULL).
#' @param m2 A positive integer representing the number of patients in group 2 for
#'        the future trial (required if \code{prob = 'predictive'}, otherwise set to NULL).
#' @param kappa01 A positive numeric value representing the prior precision parameter
#'        for group 1 (required if \code{prior = 'N-Inv-Chisq'}, otherwise set to NULL).
#' @param kappa02 A positive numeric value representing the prior precision parameter
#'        for group 2 (required if \code{prior = 'N-Inv-Chisq'} and \code{design = 'controlled'},
#'        otherwise set to NULL).
#' @param nu01 A positive numeric value representing the prior degrees of freedom
#'        for group 1 (required if \code{prior = 'N-Inv-Chisq'}, otherwise set to NULL).
#' @param nu02 A positive numeric value representing the prior degrees of freedom
#'        for group 2 (required if \code{prior = 'N-Inv-Chisq'} and \code{design = 'controlled'},
#'        otherwise set to NULL).
#' @param mu01 A numeric value representing the prior mean for group 1
#'        (required if \code{prior = 'N-Inv-Chisq'}, otherwise set to NULL).
#' @param mu02 A numeric value representing the prior mean or hypothetical control mean.
#'        For \code{design = 'controlled'} with \code{prior = 'N-Inv-Chisq'}: prior mean
#'        for group 2. For \code{design = 'uncontrolled'}: hypothetical control mean
#'        based on historical data or prior knowledge (required for uncontrolled design).
#'        Otherwise set to NULL.
#' @param sigma01 A positive numeric value representing the prior standard deviation
#'        for group 1 (required if \code{prior = 'N-Inv-Chisq'}, otherwise set to NULL).
#' @param sigma02 A positive numeric value representing the prior standard deviation
#'        for group 2 (required if \code{prior = 'N-Inv-Chisq'} and \code{design = 'controlled'},
#'        otherwise set to NULL).
#' @param bar.y1 A numeric value representing the sample mean for group 1.
#' @param bar.y2 A numeric value representing the sample mean for group 2
#'        (required if \code{design = 'controlled'} or \code{'external'}, otherwise set to NULL).
#' @param s1 A positive numeric value representing the sample standard deviation for group 1.
#' @param s2 A positive numeric value representing the sample standard deviation for group 2
#'        (required if \code{design = 'controlled'} or \code{'external'}, otherwise set to NULL).
#' @param r A positive numeric value representing the variance scaling factor for
#'        hypothetical control. For \code{design = 'uncontrolled'}: specifies how the
#'        variance of the hypothetical control relates to the treatment variance
#'        (e.g., r = 1 means equal variances). Set to NULL for controlled or external designs.
#' @param ne1 A positive integer representing the number of patients in group 1 for
#'        the external data (required if \code{design = 'external'} and external treatment
#'        data used, otherwise set to NULL).
#' @param ne2 A positive integer representing the number of patients in group 2 for
#'        the external data (required if \code{design = 'external'} and external control
#'        data used, otherwise set to NULL).
#' @param alpha01 A numeric value in (0, 1] representing the power prior scale parameter
#'        for group 1 (required if \code{design = 'external'} and external treatment data
#'        used, otherwise set to NULL). Controls the degree of borrowing.
#' @param alpha02 A numeric value in (0, 1] representing the power prior scale parameter
#'        for group 2 (required if \code{design = 'external'} and external control data
#'        used, otherwise set to NULL). Controls the degree of borrowing.
#' @param bar.ye1 A numeric value representing the sample mean for external group 1
#'        (required if \code{design = 'external'} and external treatment data used,
#'        otherwise set to NULL).
#' @param bar.ye2 A numeric value representing the sample mean for external group 2
#'        (required if \code{design = 'external'} and external control data used,
#'        otherwise set to NULL).
#' @param se1 A positive numeric value representing the sample standard deviation for
#'        external group 1 (required if \code{design = 'external'} and external treatment
#'        data used, otherwise set to NULL).
#' @param se2 A positive numeric value representing the sample standard deviation for
#'        external group 2 (required if \code{design = 'external'} and external control
#'        data used, otherwise set to NULL).
#' @param lower.tail A logical value; if TRUE (default), probabilities are
#'        P(treatment effect ≤ theta0), otherwise P(treatment effect > theta0).
#'
#' @return A numeric value in \code{[0, 1]} representing the Bayesian posterior probability
#'         or posterior predictive probability that the treatment effect exceeds
#'         (or is below) the threshold theta0.
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
#' **Design-specific handling**:
#' \itemize{
#'   \item **Controlled design**: Uses observed control data (bar.y2, s2, n2) directly
#'   \item **Uncontrolled design**: Uses hypothetical control specified by \code{mu02}
#'         (hypothetical mean) and \code{r} (variance scaling factor). The treatment
#'         variance is scaled by r to obtain the hypothetical control variance.
#'         Parameters bar.y2, s2 should be set to NULL.
#'   \item **External design**: Incorporates external data through power priors with
#'         exact conjugate representation. For vague prior, the effective posterior is
#'         computed using sufficient statistics from both current and external data.
#' }
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
#' # Example 1: Controlled design with posterior probability
#' pPPsinglecontinuous(
#'   prob = 'posterior', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'NI',
#'   theta0 = 2, n1 = 12, n2 = 12,
#'   kappa01 = 5, kappa02 = 5, nu01 = 5, nu02 = 5,
#'   mu01 = 5, mu02 = 5, sigma01 = sqrt(5), sigma02 = sqrt(5),
#'   bar.y1 = 2, bar.y2 = 0, s1 = 1, s2 = 1,
#'   m1 = NULL, m2 = NULL, r = NULL,
#'   ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
#'   bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 2: Uncontrolled design with hypothetical control
#' # mu02 = 1.5 is the hypothetical control mean
#' # r = 1.2 means hypothetical control variance is 1.2 times treatment variance
#' pPPsinglecontinuous(
#'   prob = 'posterior', design = 'uncontrolled', prior = 'vague', CalcMethod = 'WS',
#'   theta0 = 1.5, n1 = 15, n2 = 15,
#'   bar.y1 = 3.5, bar.y2 = NULL, s1 = 1.2, s2 = NULL,
#'   mu02 = 1.5, r = 1.2,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
#'   mu01 = NULL, sigma01 = NULL, sigma02 = NULL,
#'   m1 = NULL, m2 = NULL,
#'   ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
#'   bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 3: Posterior predictive probability for controlled design
#' pPPsinglecontinuous(
#'   prob = 'predictive', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'WS',
#'   theta0 = 0.5, n1 = 15, n2 = 15, m1 = 100, m2 = 100,
#'   kappa01 = 3, kappa02 = 3, nu01 = 4, nu02 = 4,
#'   mu01 = 2, mu02 = 2, sigma01 = 2, sigma02 = 2,
#'   bar.y1 = 2.5, bar.y2 = 1.8, s1 = 1.8, s2 = 1.6,
#'   r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
#'   bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' # Example 4: External control design with power prior
#' pPPsinglecontinuous(
#'   prob = 'posterior', design = 'external', prior = 'vague', CalcMethod = 'NI',
#'   theta0 = 1.5, n1 = 12, n2 = 12,
#'   bar.y1 = 4, bar.y2 = 2, s1 = 1.2, s2 = 1.1,
#'   ne2 = 20, alpha02 = 0.5, bar.ye2 = 1.8, se2 = 1.0,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
#'   mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
#'   m1 = NULL, m2 = NULL, r = NULL,
#'   ne1 = NULL, alpha01 = NULL, bar.ye1 = NULL, se1 = NULL,
#'   lower.tail = FALSE
#' )
#'
#' @export
pPPsinglecontinuous <- function(prob = "posterior", design = "controlled", prior = "vague", CalcMethod = "NI",
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

  # Validate uncontrolled design parameters
  if (design == "uncontrolled") {
    if (is.null(mu02)) {
      stop("For uncontrolled design, mu02 (hypothetical control mean) must be specified")
    }
    if (is.null(r)) {
      stop("For uncontrolled design, r (variance scaling factor) must be specified")
    }
  }

  # For external design with vague prior, use special handling
  if (design == 'external' && prior == 'vague') {
    # Group 1 (Treatment) - Apply power prior if external data available
    if (!is.null(ne1) && !is.null(alpha01)) {
      # Posterior parameters with power prior for treatment
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
      # Posterior parameters with power prior for control
      mu.t2 <- (alpha02 * ne2 * bar.ye2 + n2 * bar.y2) / (alpha02 * ne2 + n2)
      kappa.star.n2 <- alpha02 * ne2 + n2
      nu.t2 <- alpha02 * ne2 + n2 - 1
      sigma2.star.n2 <- '/'(
        '+'(
          alpha02 * (ne2 - 1) * se2 ^ 2 + (n2 - 1) * s2 ^ 2,
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
        # For uncontrolled: use same df as treatment
        nu.t2 <- nu.t1
      }

      # Calculate posterior means of t-distributions
      mu.t1 <- (kappa01 * mu01 + n1 * bar.y1) / kappa.n1
      if(design == 'controlled') {
        mu.t2 <- (kappa02 * mu02 + n2 * bar.y2) / kappa.n2
      } else if(design == 'uncontrolled') {
        # For uncontrolled: use hypothetical control mean
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
          # For uncontrolled: scale treatment sd by r
          sd.t2 <- sqrt(r) * sd.t1
        }
      } else if(prob == 'predictive') {
        sd.t1 <- sqrt((1 + 1 / kappa.n1) * var.n1 / m1)
        if(design == 'controlled') {
          sd.t2 <- sqrt((1 + 1 / kappa.n2) * var.n2 / m2)
        } else if(design == 'uncontrolled') {
          # For uncontrolled: scale treatment sd by r
          sd.t2 <- sqrt(r) * sd.t1
        }
      }
    } else if(!is.null(prior) && prior == 'vague') {
      # Calculate degrees of freedom for vague priors
      nu.t1 <- n1 - 1
      if(design == 'controlled') {
        nu.t2 <- n2 - 1
      } else if(design == 'uncontrolled') {
        # For uncontrolled: use same df as treatment
        nu.t2 <- nu.t1
      }

      # Set means of t-distributions to sample means
      mu.t1 <- bar.y1
      if(design == 'controlled') {
        mu.t2 <- bar.y2
      } else if(design == 'uncontrolled') {
        # For uncontrolled: use hypothetical control mean
        mu.t2 <- mu02
      }

      # Calculate standard deviations of t-distributions based on probability type
      if(prob == 'posterior') {
        sd.t1 <- sqrt(s1 ^ 2 / n1)
        if(design == 'controlled') {
          sd.t2 <- sqrt(s2 ^ 2 / n2)
        } else if(design == 'uncontrolled') {
          # For uncontrolled: scale treatment sd by r
          sd.t2 <- sqrt(r) * sd.t1
        }
      } else if(prob == 'predictive') {
        sd.t1 <- sqrt((1 + 1 / n1) * s1 ^ 2 / m1)
        if(design == 'controlled') {
          sd.t2 <- sqrt((1 + 1 / n2) * s2 ^ 2 / m2)
        } else if(design == 'uncontrolled') {
          # For uncontrolled: scale treatment sd by r
          sd.t2 <- sqrt(r) * sd.t1
        }
      }
    }
  }

  # Calculate the probability of below or exceeding θ₀ using the specified method
  # All three designs (controlled, uncontrolled, external) use the same calculation
  if(CalcMethod == 'NI') {
    results <- pNI2tdiff(theta0, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2, lower.tail)
  } else if(CalcMethod == 'MC') {
    results <- pMC2tdiff(nMC, theta0, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2, lower.tail)
  } else if(CalcMethod == 'WS') {
    results <- pWS2tdiff(theta0, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2, lower.tail)
  } else {
    stop("Invalid CalcMethod.")
  }

  # Return results
  return(results)
}
