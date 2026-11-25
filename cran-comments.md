## Submission

This is a minor release, version 1.1.0. It adds functions that build a
design step by step with the native pipe, reports of a design and its
results, the R code of a design, the comparison of designs, the search for
the smallest sample size and two vignettes, and the package now depends on
R (>= 4.1.0) for the native pipe. Apart from the default grid of
qdm_calibrate() and the legends of the plots of its results, the functions
of version 1.0.1 give the same results. The changes are listed in NEWS.md.

Version 1.0.1 (2026-10-31) came one day after version 1.0.0 because the
CRAN team asked us to fix a compiler warning with GCC 16; this is the first
planned update after version 1.0.0.

## Test environments

* Local: Windows 11, R 4.6.1
* win-builder: R-devel (2026-11-17 r90926)
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
