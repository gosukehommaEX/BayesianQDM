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
#'
#' @return The \code{BayesPostPredContinuous} gives the numeric value(s) of the
#' bayesian posterior probability or bayesian posterior predictive probability.
#'
#' @examples
#' BayesPostPredContinuous(
#'   design = 'controlled', prob = 'posterior', prior = 'N-Inv-Chisq', theta0 = c(2, 0), n1 = 12, n2 = 12, m1 = NULL, m2 = NULL,
#'   kappa01 = 5, kappa02 = 5, nu01 = 5, nu02 = 5, mu01 = 5, mu02 = 5, sigma01 = sqrt(5), sigma02 = sqrt(5),
#'   bar.y1 = 2, bar.y2 = 0, s1 = 1, s2 = 1, r = NULL
#' )
#' BayesPostPredContinuous(
#'   design = 'uncontrolled', prob = 'predictive', prior = 'vague', theta0 = 0.5, n1 = 12, n2 = NULL, m1 = 120, m2 = 120,
#'   kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL, mu01 = NULL, mu02 = 0, sigma01 = NULL, sigma02 = NULL,
#'   bar.y1 = 2, bar.y2 = NULL, s1 = 1, s2 = NULL, r = 12
#' )
#'
#' @importFrom fGarch dstd pstd
#' @importFrom stats integrate
#' @export
BayesPostPredContinuous = function(design, prob, prior, theta0, n1, n2, m1, m2,
                                   kappa01, kappa02, nu01, nu02, mu01, mu02, sigma01, sigma02,
                                   bar.y1, bar.y2, s1, s2, r) {
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
  # Define parameters for calculating posterior/posterior predictive probabilities
  if(prior == 'N-Inv-Chisq') {
    # Sample size
    kappa.n1 = kappa01 + n1
    kappa.n2 = kappa02 + n2
    # Degree of freedom
    nu.n1 = nu01 + n1
    if(design == 'controlled') {
      nu.n2 = nu02 + n2
    } else if(design == 'uncontrolled') {
      nu.n2 = nu.n1
    }
    # Means of t-distributions
    mu.t1 = (kappa01 * mu01 + n1 * bar.y1) / kappa.n1
    if(design == 'controlled') {
      mu.t2 = (kappa02 * mu02 + n2 * bar.y2) / kappa.n2
    } else if(design == 'uncontrolled') {
      mu.t2 = mu02
    }
    # Variance of a posterior distribution
    var.n1 = (nu01 * sigma01 ^ 2 + (n1 - 1) * s1 ^ 2 + n1 * kappa01 / (kappa01 + n1) * (mu01 - bar.y1) ^ 2) / nu.n1
    if(design == 'controlled') {
      var.n2 = (nu02 * sigma02 ^ 2 + (n2 - 1) * s2 ^ 2 + n2 * kappa02 / (kappa02 + n2) * (mu02 - bar.y2) ^ 2) / nu.n2
    } else if(design == 'uncontrolled') {
      var.n2 = NULL
    }
    # Standard deviations of t-distributions
    if(prob == 'posterior') {
      sd.t1 = sqrt(var.n1 / kappa.n1)
      if(design == 'controlled') {
        sd.t2 = sqrt(var.n2 / kappa.n2)
      } else if(design == 'uncontrolled') {
        sd.t2 = sqrt(r) * sd.t1
      }
    } else if(prob == 'predictive') {
      sd.t1 = sqrt((1 + kappa.n1) * var.n1 / (kappa.n1 * m1))
      if(design == 'controlled') {
        sd.t2 = sqrt((1 + kappa.n2) * var.n2 / (kappa.n2 * m2))
      } else if(design == 'uncontrolled') {
        sd.t2 = sqrt(r) * sd.t1
      }
    }
  } else if(prior == 'vague') {
    # Degree of freedom
    nu.n1 = n1 - 1
    if(design == 'controlled') {
      nu.n2 = n2 - 1
    } else if(design == 'uncontrolled') {
      nu.n2 = nu.n1
    }
    # Means of t-distributions
    mu.t1 = bar.y1
    if(design == 'controlled') {
      mu.t2 = bar.y2
    } else if(design == 'uncontrolled') {
      mu.t2 = mu02
    }
    # Standard deviations of t-distributions
    if(prob == 'posterior') {
      sd.t1 = sqrt(s1 ^ 2 / n1)
      if(design == 'controlled') {
        sd.t2 = sqrt(s2 ^ 2 / n2)
      } else if(design == 'uncontrolled') {
        sd.t2 = sqrt(r) * sd.t1
      }
    } else if(prob == 'predictive') {
      sd.t1 = sqrt((1 + n1) * s1 ^ 2 / (n1 * m1))
      if(design == 'controlled') {
        sd.t2 = sqrt((1 + n2) * s2 ^ 2 / (n2 * m2))
      } else if(design == 'uncontrolled') {
        sd.t2 = sqrt(r) * sd.t1
      }
    }
  }
  # The probability of exceeding \theta_{0}
  g = Vectorize(
    function(mu.t1, s.t1, mu.t2, s.t2, Theta0) {
      integrate(
        function(x, mu.t1, s.t1, mu.t2, s.t2, Theta0) {
          '*'(
            fGarch::dstd(x,          mean = mu.t1, sd = s.t1, nu = nu.n1),
            fGarch::pstd(x - Theta0, mean = mu.t2, sd = s.t2, nu = nu.n2)
          )
        },
        -Inf, Inf, mu.t1 = mu.t1, s.t1 = s.t1, mu.t2 = mu.t2, s.t2 = s.t2, Theta0 = Theta0
      )[['value']]
    }
  )
  # Result
  result = data.frame(matrix(sapply(theta0, function(i) g(mu.t1, sd.t1, mu.t2, sd.t2, i)), ncol = length(theta0)))
  if(length(theta0) == 2) {
    names(result) = paste0(ifelse(prob == 'posterior', 'g.post.', 'g.pred.'), c('MAV', 'TV'))
  } else if(length(theta0) == 1) {
    names(result) = paste0(ifelse(prob == 'posterior', 'g.post.', 'g.pred.'), 'NULL')
  }
  return(result)
}
