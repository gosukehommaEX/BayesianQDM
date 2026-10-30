## Submission

This is a major update of BayesianQDM from version 0.1.0 to 1.0.0. The
interface has been redesigned and is not backward compatible: a study is
specified by small objects combined into a design, and three functions
compute the decision for observed data, the operating characteristics and
calibrated decision thresholds for any design. Formulas and numerical methods
that were wrong in version 0.1.0 have been corrected. NEWS.md lists the
changes and the correspondence between the old and new functions, and a new
vignette explains the migration. Version 1.0.0 also adds a time-to-event
endpoint evaluated by the milestone survival probability or the hazard
ratio, single-arm designs compared with a fixed control value, decision
criteria for endpoints on which smaller effects are favorable, and tables
and plots of the observed data that lead to each decision.

## Test environments

* Local: Windows 11, R 4.6.0
* win-builder: R-devel (2026-10-28 r90783) and R-release (4.6.2)
* GitHub Actions: ubuntu-latest (R-devel, R-release and R-oldrel-1),
  windows-latest (R-release) and macos-latest (R-release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies on CRAN.

## Notes

* The heavy computations (numerical integration and Monte Carlo loops) are
  written in C++ with Rcpp.
* The examples, tests and vignettes use small numbers of simulated studies
  and Monte Carlo draws; the tests take less than two minutes in R CMD check.
