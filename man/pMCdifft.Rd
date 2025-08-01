% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pMCdifft.R
\name{pMCdifft}
\alias{pMCdifft}
\title{Cumulative Distribution Function of the Difference of Two t-Distributed Variables by Monte Carlo Simulation}
\usage{
pMCdifft(nMC, q, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2)
}
\arguments{
\item{nMC}{A positive integer representing the number of Monte Carlo iterations for simulation.}

\item{q}{A numeric value representing the quantile threshold.}

\item{mu.t1}{A numeric value representing the location parameter of the first t-distribution.}

\item{mu.t2}{A numeric value representing the location parameter of the second t-distribution.}

\item{sd.t1}{A positive numeric value representing the scale parameter of the first t-distribution.}

\item{sd.t2}{A positive numeric value representing the scale parameter of the second t-distribution.}

\item{nu.t1}{A positive numeric value representing the degrees of freedom of the first t-distribution.}

\item{nu.t2}{A positive numeric value representing the degrees of freedom of the second t-distribution.}
}
\value{
A numeric vector representing P(T1 - T2 > q), the probability that the difference
between the two t-distributed variables exceeds the quantile q.
}
\description{
This function calculates the cumulative distribution function (CDF) of the difference
between two independent t-distributed random variables using Monte Carlo simulation.
Specifically, it computes P(T1 - T2 > q) where T1 and T2 follow t-distributions
with potentially different location, scale, and degrees of freedom parameters.
}
\details{
This function uses Monte Carlo simulation to approximate the distribution of the
difference between two t-distributed variables. The method involves:
\itemize{
\item Generating nMC random samples from each t-distribution
\item Computing the difference for each pair of samples
\item Calculating the proportion of differences that exceed the threshold q
\item The accuracy increases with larger nMC values
}

Monte Carlo simulation provides a flexible and intuitive approach that can handle
any combination of parameters, but the accuracy depends on the number of simulations.
Larger nMC values yield more precise estimates but require more computational time.
}
\examples{
# Calculate P(t1 - t2 > 3) for equal parameters using 10,000 simulations
pMCdifft(nMC = 1e+4, q = 3, mu.t1 = 2, mu.t2 = 0, sd.t1 = 1, sd.t2 = 1, nu.t1 = 17, nu.t2 = 17)

# Calculate P(t1 - t2 > 1) for unequal variances using 50,000 simulations
pMCdifft(nMC = 5e+4, q = 1, mu.t1 = 5, mu.t2 = 3, sd.t1 = 2, sd.t2 = 1.5, nu.t1 = 10, nu.t2 = 15)

# Calculate P(t1 - t2 > 0) for high precision with 100,000 simulations
pMCdifft(nMC = 1e+5, q = 0, mu.t1 = 1, mu.t2 = 0, sd.t1 = 1, sd.t2 = 1, nu.t1 = 5, nu.t2 = 20)

}
