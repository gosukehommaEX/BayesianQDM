% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayesDecisionProbBinary.R
\name{BayesDecisionProbBinary}
\alias{BayesDecisionProbBinary}
\title{Calculate the Go, NoGo and Gray Probabilities for a Clinical Trial When Outcome is Binary
Under the Bayesian Framework Using Two Metrics}
\usage{
BayesDecisionProbBinary(
  prob = "posterior",
  design = "controlled",
  theta.TV,
  theta.MAV,
  theta.NULL = NULL,
  gamma1,
  gamma2,
  pi1,
  pi2,
  n1,
  n2,
  a1,
  a2,
  b1,
  b2,
  z = NULL,
  m1,
  m2,
  ne1,
  ne2,
  ye1,
  ye2,
  ae1,
  ae2
)
}
\arguments{
\item{prob}{A character string specifying the type of probability to use
(\code{prob = 'posterior'} or \code{prob = 'predictive'}).}

\item{design}{A character string specifying the type of design
(\code{design = 'controlled'}, \code{design = 'uncontrolled'}, or \code{design = 'external'}).}

\item{theta.TV}{A numeric value representing the pre-specified threshold value for calculating
Go probability when \code{prob = 'posterior'}.}

\item{theta.MAV}{A numeric value representing the pre-specified threshold value for calculating
NoGo probability when \code{prob = 'posterior'}.}

\item{theta.NULL}{A numeric value representing the pre-specified threshold value for calculating
Go/NoGo probabilities when \code{prob = 'predictive'}.}

\item{gamma1}{A numeric value between 0 and 1 representing the minimum probability to declare success.}

\item{gamma2}{A numeric value between 0 and 1 representing the futility threshold.}

\item{pi1}{A numeric value or vector representing true response probability(s) for group 1.}

\item{pi2}{A numeric value or vector representing true response probability(s) for group 2.}

\item{n1}{A positive integer representing the number of patients in group 1 for a proof-of-concept (PoC) trial.}

\item{n2}{A positive integer representing the number of patients in group 2 for the PoC trial.}

\item{a1}{A positive numeric value representing the first shape parameter of the prior distribution for group 1.}

\item{a2}{A positive numeric value representing the first shape parameter of the prior distribution for group 2.}

\item{b1}{A positive numeric value representing the second shape parameter of the prior distribution for group 1.}

\item{b2}{A positive numeric value representing the second shape parameter of the prior distribution for group 2.}

\item{z}{A non-negative integer representing the hypothetical observed number of responders in group 2 for an uncontrolled design.}

\item{m1}{A positive integer representing the number of patients in group 1 for the future trial data.}

\item{m2}{A positive integer representing the number of patients in group 2 for the future trial data.}

\item{ne1}{A positive integer representing the number of patients in group 1 for the external data.}

\item{ne2}{A positive integer representing the number of patients in group 2 for the external data.}

\item{ye1}{A non-negative integer representing the observed number of responders in group 1 for the external data.}

\item{ye2}{A non-negative integer representing the observed number of responders in group 2 for the external data.}

\item{ae1}{A positive numeric value representing the scale parameter (power parameter) for group 1.}

\item{ae2}{A positive numeric value representing the scale parameter (power parameter) for group 2.}
}
\value{
A data frame containing the true response probabilities for both groups, and the
Go, NoGo, and Gray probabilities.
}
\description{
This function calculates Go, NoGo, and Gray probabilities for binary outcome clinical trials
under the Bayesian framework using two metrics: (i) posterior probability for the treatment
effect to be greater than a threshold, and (ii) posterior predictive probability of phase III
study success. The function supports controlled, uncontrolled, and external control designs.
}
\details{
The function can obtain:
\itemize{
\item Go probability
\item NoGo probability
\item Gray probability
}

The function can be used for controlled design, uncontrolled design, and design using
external (historical) data. The decision framework is based on:
\itemize{
\item Go: Probability that the treatment effect exceeds the efficacy threshold
\item NoGo: Probability that the treatment effect is below the futility threshold
\item Gray: Intermediate zone where neither Go nor NoGo criteria are met
}
}
\examples{
# Calculate Go/NoGo/Gray probabilities using posterior probability for controlled design
BayesDecisionProbBinary(
  prob = 'posterior', design = 'controlled', theta.TV = 0.4, theta.MAV = 0.2, theta.NULL = NULL,
  gamma1 = 0.5, gamma2 = 0.2, pi1 = c(0.2, 0.4, 0.6, 0.8), pi2 = rep(0.2, 4), n1 = 12, n2 = 12,
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5, z = NULL, m1 = NULL, m2 = NULL, ne1 = NULL, ne2 = NULL,
  ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
)

# Calculate Go/NoGo/Gray probabilities using posterior predictive probability for controlled design
BayesDecisionProbBinary(
  prob = 'predictive', design = 'controlled', theta.TV = NULL, theta.MAV = NULL, theta.NULL = 0,
  gamma1 = 0.9, gamma2 = 0.3, pi1 = c(0.2, 0.4, 0.6, 0.8), pi2 = rep(0.2, 4), n1 = 12, n2 = 12,
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5, z = NULL, m1 = 30, m2 = 30, ne1 = NULL, ne2 = NULL,
  ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
)

}
