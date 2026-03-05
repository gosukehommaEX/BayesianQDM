# CRAN Submission Comments

## Test Environments

* local R installation, R 4.5.1
* ubuntu-latest (on GitHub Actions), R-release
* windows-latest (on GitHub Actions), R-release
* macOS-latest (on GitHub Actions), R-release
* ubuntu-latest (on GitHub Actions), R-devel

## R CMD check results

0 errors | 0 warnings | 0 notes

## Package Summary

BayesianQDM provides a comprehensive Bayesian quantitative decision-making framework
for clinical trials, supporting both single and two-endpoint analyses for binary and
continuous outcomes. The package implements Go/NoGo/Gray decision criteria using
posterior probabilities and posterior predictive probabilities under conjugate Bayesian
models, with support for controlled, uncontrolled, and external designs.
The package also provides grid-search functions to find optimal Go and NoGo thresholds
that satisfy user-specified operating characteristic criteria for all supported
endpoint types and study designs.

## Dependencies

* Imports: mvtnorm, stats
* Suggests: testthat (>= 3.0.0), knitr, rmarkdown, dplyr, ggplot2, tidyr, purrr
* All dependencies are well-established CRAN packages

## Notes

* No MCMC sampling is used; all computations rely on exact conjugate posterior
  representations or closed-form approximations, keeping run times short
* Monte Carlo functions have an adjustable iteration count (`nMC`) to balance
  speed and precision; default values are set conservatively for CRAN checks
* All examples in the documentation are designed to complete within a few seconds

## Previous CRAN Submissions

This is the initial submission to CRAN.
