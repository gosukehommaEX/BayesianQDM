## Submission

This is a patch release of BayesianQDM 1.0.0. It fixes the installation
warning on r-devel-linux-x86_64-fedora-gcc (-Warray-bounds with GCC 16 in
src/regions_continuous2.cpp), which the CRAN team asked us to correct before
2026-11-20. The warning came from std::sort() on a fixed-size array in which
at most two elements are sorted; the call is replaced by a comparison and a
swap, so the results are unchanged. We do not have access to GCC 16; with
the compilers used in the test environments below, the package installs
without warnings.

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
