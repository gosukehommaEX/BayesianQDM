#' Calculate the Go, NoGo and Gray probabilities for a clinical trial when outcome is binary
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
#' The function can be used for controlled design, uncontrolled design and design using an external (historical) data.
#'
#' @param prob A type of probability (i.e., \code{prob='posterior'} or \code{prob='predictive'}).
#' @param design A type of design (i.e., \code{design='controlled'}, \code{design='uncontrolled'} or \code{design='external'}).
#' @param theta.TV A numeric pre-specified threshold value for calculating Go probability when \code{prob='posterior'}.
#' @param theta.MAV A numeric pre-specified threshold value for calculating NoGo probability when \code{prob='posterior'}.
#' @param theta.NULL A numeric pre-specified threshold value for calculating Go/NoGo probabilities when \code{prob='predictive'}.
#' @param gamma1 A numeric value of a minimum probability to declare success.
#' @param gamma2 A numeric value of a futility threshold.
#' @param pi1 A numeric value or a vector representing true response probability(s) for group 1.
#' @param pi2 A numeric value or a vector representing true response probability(s) for group 2.
#' @param n1 A number of patients in group 1 for a proof-of-concept (PoC) trial.
#' @param n2 A number of patients in group 2 for the PoC trial.
#' @param a1 The first non-negative numeric value of the prior distribution for group 1.
#' @param a2 The first non-negative numeric value of the prior distribution for group 2.
#' @param b1 The second non-negative numeric value of the prior distribution for group 1.
#' @param b2 The second non-negative numeric value of the prior distribution for group 2.
#' @param z A hypothetical observed number of responders in group 2 for an uncontrolled design.
#' @param m1 A number of patients in group 1 for the future trial data.
#' @param m2 A number of patients in group 2 for the future trial data.
#' @param ne1 A number of patients in group 1 for the external data.
#' @param ne2 A number of patients in group 2 for the external data.
#' @param ye1 An observed number of responders in group 1 for the external data.
#' @param ye2 An observed number of responders in group 2 for the external data.
#' @param ae1 A scale parameter (power parameter) for group 1.
#' @param ae2 A scale parameter (power parameter) for group 2.
#'
#' @return The \code{BayesDecisionProbBinary} gives the dataframe including the
#'  true response probabilities for both groups, and the Go, NoGo and Gray probabilities.
#'
#' @examples
#' BayesDecisionProbBinary(
#'   prob = 'posterior', design = 'controlled', theta.TV = 0.4, theta.MAV = 0.2, theta.NULL = NULL,
#'   gamma1 = 0.5, gamma2 = 0.2, pi1 = c(0.2, 0.4, 0.6, 0.8), pi2 = rep(0.2, 4), n1 = 12, n2 = 12,
#'   a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5, z = NULL, m1 = NULL, m2 = NULL, ne1 = NULL, ne2 = NULL,
#'   ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
#' )
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
  # Check parameter sets
  if((prob == 'posterior') & (sum(sapply(list(theta.TV, theta.MAV), is.null)) > 0)) {
    stop('If you calculate the Go, NoGo and Gray probabilities using posterior probability, theta.TV and theta.MAV should be non-null')
  }
  if((prob == 'predictive') & (sum(sapply(list(m1, m2), is.null)) > 0)) {
    stop('If you calculate the Go, NoGo and Gray probabilities using posterior predictive probability, m1 and m2 should be non-null')
  }
  if((prob == 'predictive') & (is.null(theta.NULL))) {
    stop('If you calculate the Go, NoGo and Gray probabilities using posterior predictive probability, theta.NULL should be non-null')
  }
  if((design == 'uncontrolled') & (is.null(z))) {
    stop('If you consider uncontrolled design, z should be non-null')
  }
  if((design == 'external') & (sum(sapply(list(ne1, ne2, ye1, ye2, ae1, ae2), is.null)) > 0)) {
    stop('If you use the external data, ne1, ne2, ye1, ye2, ae1 and ae2 should be non-null')
  }
  # Set values of theta0
  if(prob == 'posterior') {
    theta0 <- c(theta.TV, theta.MAV)
  } else {
    theta0 <- theta.NULL
  }
  # Calculate bayesian posterior probability or bayesian posterior predictive probability
  Y1 <- 0:n1
  if(design == 'uncontrolled') { Y2 <- z } else { Y2 <- 0:n2 }
  # Posterior/Posterior predictive probabilities
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
  # Go/NoGo probabilities
  GoNogoProb <- matrix(
    sapply(seq(2), function(j) {
      I <- matrix((c(1, -1)[j] * gPost[[ifelse(prob == 'posterior', j, 1)]] >= c(gamma1, -gamma2)[j]), nrow = length(Y2))
      if(design == 'uncontrolled') {
        colSums(outer(col(I)[I] - 1, pi1, FUN = function(X, Y) dbinom(X, n1, Y)))
      } else {
        diag(crossprod(
          outer(col(I)[I] - 1, pi1, FUN = function(X, Y) dbinom(X, n1, Y)),
          outer(row(I)[I] - 1, pi2, FUN = function(X, Y) dbinom(X, n2, Y))
        ))
      }
    }),
    ncol = 2
  )
  # Gray probability
  GrayProb <- 1 - rowSums(GoNogoProb)
  if(sum(GrayProb < 0) > 0) {
    print('Because negative gray probability(s) is obtained, re-consider appropriate threshold')
  }
  # Results
  results <- data.frame(
    pi1, pi2, Go = GoNogoProb[, 1], Gray = GrayProb, NoGo = GoNogoProb[, 2]
  )
  return(results)
}
