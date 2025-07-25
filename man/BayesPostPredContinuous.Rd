% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BayesPostPredContinuous.R
\name{BayesPostPredContinuous}
\alias{BayesPostPredContinuous}
\title{Calculate Bayesian Posterior Probability or Bayesian Posterior Predictive Probability for a Clinical Trial When Outcome is Continuous}
\usage{
BayesPostPredContinuous(
  prob = "posterior",
  design = "controlled",
  prior = "vague",
  CalcMethod = "NI",
  theta0,
  nMC = NULL,
  nMCMCsample = NULL,
  n1,
  n2,
  m1,
  m2,
  kappa01,
  kappa02,
  nu01,
  nu02,
  mu01,
  mu02,
  sigma01,
  sigma02,
  bar.y1,
  bar.y2,
  s1,
  s2,
  r = NULL,
  ne1 = NULL,
  ne2 = NULL,
  alpha01 = NULL,
  alpha02 = NULL
)
}
\arguments{
\item{prob}{A character string specifying the type of probability to calculate.
Options are 'posterior' (default) for posterior probability or 'predictive' for posterior predictive probability.}

\item{design}{A character string specifying the trial design.
Options are 'controlled' (default), 'uncontrolled', or 'external'.}

\item{prior}{A character string specifying the prior distribution.
Options are 'vague' (default) or 'N-Inv-Chisq' for Normal-Inverse-Chi-squared.}

\item{CalcMethod}{A character string specifying the calculation method.
Options are 'NI' (numerical integration, default), 'MC' (Monte Carlo), 'WS' (Welch-Satterthwaite),
or 'MCMC' (MCMC sampling for external design).}

\item{theta0}{A numeric value representing the threshold for the treatment effect.}

\item{nMC}{A positive integer representing the number of Monte Carlo iterations for MC method (default: NULL).}

\item{nMCMCsample}{A positive integer representing the number of MCMC iterations for external design (default: NULL).}

\item{n1}{A positive integer representing the sample size for group 1 in PoC trial.}

\item{n2}{A positive integer representing the sample size for group 2 in PoC trial.}

\item{m1}{A positive integer representing the sample size for group 1 in future trial (for predictive probability).}

\item{m2}{A positive integer representing the sample size for group 2 in future trial (for predictive probability).}

\item{kappa01}{A positive numeric value representing the prior precision parameter for group 1 (N-Inv-Chisq prior).}

\item{kappa02}{A positive numeric value representing the prior precision parameter for group 2 (N-Inv-Chisq prior).}

\item{nu01}{A positive numeric value representing the prior degrees of freedom for group 1 (N-Inv-Chisq prior).}

\item{nu02}{A positive numeric value representing the prior degrees of freedom for group 2 (N-Inv-Chisq prior).}

\item{mu01}{A numeric value representing the prior mean for group 1 (N-Inv-Chisq prior).}

\item{mu02}{A numeric value representing the prior mean for group 2 (N-Inv-Chisq prior).}

\item{sigma01}{A positive numeric value representing the prior standard deviation for group 1 (N-Inv-Chisq prior).}

\item{sigma02}{A positive numeric value representing the prior standard deviation for group 2 (N-Inv-Chisq prior).}

\item{bar.y1}{A numeric value representing the sample mean of group 1.}

\item{bar.y2}{A numeric value representing the sample mean of group 2.}

\item{s1}{A positive numeric value representing the sample standard deviation of group 1.}

\item{s2}{A positive numeric value representing the sample standard deviation of group 2.}

\item{r}{A positive numeric value for uncontrolled design (default: NULL).}

\item{ne1}{A positive integer representing the sample size for group 1 in external trial (default: NULL).}

\item{ne2}{A positive integer representing the sample size for group 2 in external trial (default: NULL).}

\item{alpha01}{A positive numeric value representing the power prior scale parameter for group 1 (default: NULL).}

\item{alpha02}{A positive numeric value representing the power prior scale parameter for group 2 (default: NULL).}
}
\value{
A numeric vector representing the Bayesian posterior probability or Bayesian posterior
predictive probability. The function can handle vectorized inputs.
}
\description{
This function computes Bayesian posterior probability or posterior predictive probability
for continuous outcome clinical trials. The function supports controlled, uncontrolled, and
external control designs with Normal-Inverse-Chi-squared or vague priors, using three calculation
methods: numerical integration, Monte Carlo simulation, and Welch-Satterthwaite approximation.
For external control designs, MCMC sampling is used to incorporate historical data through power priors.
}
\details{
The function can obtain:
\itemize{
\item Bayesian posterior probability
\item Bayesian posterior predictive probability
}

Prior distribution of mean and variance of outcomes for each treatment group (k=1,2) can be either
(1) Normal-Inverse-Chi-squared or (2) Vague. The posterior distribution or posterior predictive
distribution of outcome for each treatment group follows a t-distribution.

For controlled and uncontrolled designs, three calculation methods are available:
\itemize{
\item NI: Numerical integration method for exact computation
\item MC: Monte Carlo simulation for flexible approximation
\item WS: Welch-Satterthwaite approximation for computational efficiency
}

For external control designs, the function uses MCMC sampling to incorporate historical
data through power prior methodology:
\itemize{
\item MCMC: Markov Chain Monte Carlo sampling for posterior inference with external data
\item Power priors allow controlled borrowing from historical data
\item alpha parameters control the degree of borrowing (0 = no borrowing, 1 = full borrowing)
}

The external design supports:
\itemize{
\item External control data only (ne2, alpha02)
\item External treatment data only (ne1, alpha01)
\item Both external control and treatment data
}
}
\examples{
# Example 1: Numerical Integration (NI) method with N-Inv-Chisq prior
BayesPostPredContinuous(
  prob = 'posterior', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'NI',
  theta0 = 2, n1 = 12, n2 = 12, kappa01 = 5, kappa02 = 5, nu01 = 5, nu02 = 5,
  mu01 = 5, mu02 = 5, sigma01 = sqrt(5), sigma02 = sqrt(5),
  bar.y1 = 2, bar.y2 = 0, s1 = 1, s2 = 1
)

# Example 2: Monte Carlo (MC) method with vague prior
BayesPostPredContinuous(
  prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'MC',
  theta0 = 1, nMC = 10000, n1 = 12, n2 = 12,
  bar.y1 = 3, bar.y2 = 1, s1 = 1.5, s2 = 1.2
)

# Example 3: Welch-Satterthwaite (WS) approximation with N-Inv-Chisq prior
BayesPostPredContinuous(
  prob = 'predictive', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'WS',
  theta0 = 0.5, n1 = 15, n2 = 15, m1 = 100, m2 = 100,
  kappa01 = 3, kappa02 = 3, nu01 = 4, nu02 = 4, mu01 = 2, mu02 = 2,
  sigma01 = 2, sigma02 = 2, bar.y1 = 2.5, bar.y2 = 1.8, s1 = 1.8, s2 = 1.6
)

# Example 4: External control design with MCMC method
BayesPostPredContinuous(
  prob = 'posterior', design = 'external', CalcMethod = 'MCMC',
  theta0 = 1.5, nMCMCsample = 5000, n1 = 12, n2 = 12,
  bar.y1 = 4, bar.y2 = 2, s1 = 1.2, s2 = 1.1,
  ne2 = 20, alpha02 = 0.5
)

}
