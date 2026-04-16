# CRAN Submission Comments

## Test Environments

* local R installation, R 4.5.1
* ubuntu-latest (on GitHub Actions), R-release
* windows-latest (on GitHub Actions), R-release
* macOS-latest (on GitHub Actions), R-release
* ubuntu-latest (on GitHub Actions), R-devel

## R CMD check results

0 errors | 0 warnings | 1 note

* checking for future file timestamps ... NOTE
  unable to verify current time
  This is a transient network issue in the check environment and is not
  related to the package itself.

## Package Summary

BayesianQDM provides a comprehensive Bayesian quantitative decision-making framework
for clinical trials, supporting both single and two-endpoint analyses for binary and
continuous outcomes. The package implements Go/NoGo/Gray decision criteria using
posterior probabilities and posterior predictive probabilities under conjugate Bayesian
models, with support for controlled, uncontrolled, and external designs.
The package also provides grid-search functions to find optimal Go and NoGo thresholds
that satisfy user-specified operating characteristic criteria for all supported
endpoint types and study designs. S3 print() and plot() methods are provided for all
four decision probability classes, enabling formatted display and visualisation of
operating characteristics across treatment scenarios.

## Dependencies

* Imports: ggplot2 (>= 3.4.0), gridExtra, mvtnorm, stats
* Suggests: testthat (>= 3.0.0), knitr, rmarkdown, dplyr, tidyr, purrr
* All dependencies are well-established CRAN packages

## Notes

* No MCMC sampling is used; all computations rely on exact conjugate posterior
  representations or closed-form approximations, keeping run times short
* Monte Carlo functions have an adjustable iteration count (`nMC`) to balance
  speed and precision; default values are set conservatively for CRAN checks
* All examples in the documentation are designed to complete within a few seconds

## Resubmission

This is a resubmission following reviewer feedback. The following changes were made:

* Expanded the MCMC acronym to "Markov chain Monte Carlo (MCMC)" in DESCRIPTION.
* Added a reference to Kang, Yamaguchi, and Han (2026) in DESCRIPTION.
* Replaced `\dontrun{}` with `\donttest{}` in examples for `getgamma2bin()` and
  `getgamma2cont()`.
