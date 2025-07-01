#' Calculate the Go, NoGo and Gray Probabilities for a Clinical Trial When Outcome is Binary
#' Under the Bayesian Framework Using Two Metrics
#'
#' This function calculates Go, NoGo, and Gray probabilities for binary outcome clinical trials
#' under the Bayesian framework using two metrics: (i) posterior probability for the treatment
#' effect to be greater than a threshold, and (ii) posterior predictive probability of phase III
#' study success. The function supports controlled, uncontrolled, and external control designs.
#'
#' @param prob A character string specifying the type of probability to use
#'        (\code{prob = 'posterior'} or \code{prob = 'predictive'}).
#' @param design A character string specifying the type of design
#'        (\code{design = 'controlled'}, \code{design = 'uncontrolled'}, or \code{design = 'external'}).
#' @param theta.TV A numeric value representing the pre-specified threshold value for calculating
#'        Go probability when \code{prob = 'posterior'}.
#' @param theta.MAV A numeric value representing the pre-specified threshold value for calculating
#'        NoGo probability when \code{prob = 'posterior'}.
#' @param theta.NULL A numeric value representing the pre-specified threshold value for calculating
#'        Go/NoGo probabilities when \code{prob = 'predictive'}.
#' @param gamma1 A numeric value between 0 and 1 representing the minimum probability to declare success.
#' @param gamma2 A numeric value between 0 and 1 representing the futility threshold.
#' @param pi1 A numeric value or vector representing true response probability(s) for group 1.
#' @param pi2 A numeric value or vector representing true response probability(s) for group 2.
#' @param n1 A positive integer representing the number of patients in group 1 for a proof-of-concept (PoC) trial.
#' @param n2 A positive integer representing the number of patients in group 2 for the PoC trial.
#' @param a1 A positive numeric value representing the first shape parameter of the prior distribution for group 1.
#' @param a2 A positive numeric value representing the first shape parameter of the prior distribution for group 2.
#' @param b1 A positive numeric value representing the second shape parameter of the prior distribution for group 1.
#' @param b2 A positive numeric value representing the second shape parameter of the prior distribution for group 2.
#' @param z A non-negative integer representing the hypothetical observed number of responders in group 2 for an uncontrolled design.
#' @param m1 A positive integer representing the number of patients in group 1 for the future trial data.
#' @param m2 A positive integer representing the number of patients in group 2 for the future trial data.
#' @param ne1 A positive integer representing the number of patients in group 1 for the external data.
#' @param ne2 A positive integer representing the number of patients in group 2 for the external data.
#' @param ye1 A non-negative integer representing the observed number of responders in group 1 for the external data.
#' @param ye2 A non-negative integer representing the observed number of responders in group 2 for the external data.
#' @param ae1 A positive numeric value representing the scale parameter (power parameter) for group 1.
#' @param ae2 A positive numeric value representing the scale parameter (power parameter) for group 2.
#'
#' @return A data frame containing the true response probabilities for both groups, and the
#'         Go, NoGo, and Gray probabilities.
#'
#' @details
#' The function can obtain:
#' \itemize{
#'   \item Go probability
#'   \item NoGo probability
#'   \item Gray probability
#' }
#'
#' The function can be used for controlled design, uncontrolled design, and design using
#' external (historical) data. The decision framework is based on:
#' \itemize{
#'   \item Go: Probability that the treatment effect exceeds the efficacy threshold
#'   \item NoGo: Probability that the treatment effect is below the futility threshold
#'   \item Gray: Intermediate zone where neither Go nor NoGo criteria are met
#' }
#'
#' @examples
#' # Calculate Go/NoGo/Gray probabilities using posterior probability for controlled design
#' BayesDecisionProbBinary(
#'   prob = 'posterior', design = 'controlled', theta.TV = 0.4, theta.MAV = 0.2, theta.NULL = NULL,
#'   gamma1 = 0.5, gamma2 = 0.2, pi1 = c(0.2, 0.4, 0.6, 0.8), pi2 = rep(0.2, 4), n1 = 12, n2 = 12,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5, z = NULL, m1 = NULL, m2 = NULL, ne1 = NULL, ne2 = NULL,
#'   ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
#' )
#'
#' # Calculate Go/NoGo/Gray probabilities using posterior predictive probability for controlled design
#' BayesDecisionProbBinary(
#'   prob = 'predictive', design = 'controlled', theta.TV = NULL, theta.MAV = NULL, theta.NULL = 0,
#'   gamma1 = 0.9, gamma2 = 0.3, pi1 = c(0.2, 0.4, 0.6, 0.8), pi2 = rep(0.2, 4), n1 = 12, n2 = 12,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5, z = NULL, m1 = 30, m2 = 30, ne1 = NULL, ne2 = NULL,
#'   ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
#' )
#'
#' @importFrom stats dbinom
#' @export
BayesDecisionProbBinary <- function(prob = 'posterior', design = 'controlled', theta.TV, theta.MAV, theta.NULL = NULL, gamma1, gamma2,
                                    pi1, pi2, n1, n2, a1, a2, b1, b2, z = NULL,
                                    m1, m2, ne1, ne2, ye1, ye2, ae1, ae2) {
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
  if((design == 'uncontrolled') & (is.null(z))) {
    stop('If you consider uncontrolled design, z should be non-null')
  }

  # Check parameter sets for external design
  if((design == 'external') & (sum(sapply(list(ne1, ne2, ye1, ye2, ae1, ae2), is.null)) > 0)) {
    stop('If you use the external data, ne1, ne2, ye1, ye2, ae1 and ae2 should be non-null')
  }

  # Set values of theta0 based on probability type
  if(prob == 'posterior') {
    theta0 <- c(theta.TV, theta.MAV)
  } else {
    theta0 <- theta.NULL
  }

  # Calculate Bayesian posterior probability or Bayesian posterior predictive probability
  # Define possible outcomes for group 1
  Y1 <- 0:n1
  # Define possible outcomes for group 2 (fixed value for uncontrolled design)
  if(design == 'uncontrolled') { Y2 <- z } else { Y2 <- 0:n2 }

  # Calculate posterior/posterior predictive probabilities for each threshold
  gPost <- lapply(seq(length(theta0)), function(i) {
    sapply(Y1, function(y1) {
      sapply(Y2, function(y2) {
        BayesPostPredBinary(
          prob, design, theta0[i],
          n1, n2, y1, y2, a1, a2, b1, b2,
          m1, m2, ne1, ne2, ye1, ye2, ae1, ae2
        )
      })
    })
  })

  # Calculate Go/NoGo probabilities based on decision criteria
  GoNogoProb <- matrix(
    sapply(seq(2), function(j) {
      # Create indicator matrix for Go (j=1) or NoGo (j=2) decisions
      I <- matrix((c(1, -1)[j] * gPost[[ifelse(prob == 'posterior', j, 1)]] >= c(gamma1, -gamma2)[j]), nrow = length(Y2))

      if(design == 'uncontrolled') {
        # For uncontrolled design, sum over group 1 outcomes only
        colSums(outer(col(I)[I] - 1, pi1, FUN = function(X, Y) dbinom(X, n1, Y)))
      } else {
        # For controlled design, sum over both group outcomes
        diag(crossprod(
          outer(col(I)[I] - 1, pi1, FUN = function(X, Y) dbinom(X, n1, Y)),
          outer(row(I)[I] - 1, pi2, FUN = function(X, Y) dbinom(X, n2, Y))
        ))
      }
    }),
    ncol = 2
  )

  # Calculate Gray probability (complement of Go and NoGo)
  GrayProb <- 1 - rowSums(GoNogoProb)

  # Check for negative Gray probabilities
  if(sum(GrayProb < 0) > 0) {
    print('Because negative gray probability(s) is obtained, re-consider appropriate threshold')
  }

  # Prepare results data frame
  results <- data.frame(
    pi1, pi2, Go = GoNogoProb[, 1], Gray = GrayProb, NoGo = GoNogoProb[, 2]
  )

  return(results)
}
