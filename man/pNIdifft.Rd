% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pNIdifft.R
\name{pNIdifft}
\alias{pNIdifft}
\title{Cumulative distribution function (CDF) of the difference of two t-distributed outcomes by convolution of a difference of the t-distributions.}
\usage{
pNIdifft(q, mu.t1, mu.t2, sd.t1, sd.t2, nu.t1, nu.t2)
}
\arguments{
\item{q}{Quantile value.}

\item{mu.t1}{Location parameter of t-distribution for group 1.}

\item{mu.t2}{Location parameter of t-distribution for group 2.}

\item{sd.t1}{Scale parameter of t-distribution for group 1.}

\item{sd.t2}{Scale parameter of t-distribution for group 2.}

\item{nu.t1}{Degrees of freedom of t-distribution for group 1.}

\item{nu.t2}{Degrees of freedom of t-distribution for group 2.}
}
\value{
Cumulative probability(s)
}
\description{
This function calculates CDF(s) of difference between two variables following t-distribution
using a numerical integration.
}
\examples{
pNIdifft(q = 3, mu.t1 = 2, mu.t2 = 0, sd.t1 = 1, sd.t2 = 1, nu.t1 = 17, nu.t2 = 17)

}
