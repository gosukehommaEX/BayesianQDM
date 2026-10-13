# BayesianQDM (development version)

This version is a complete redesign of the package towards version 1.0.0.
The functions of version 0.1.0 are being replaced by a new interface in which
a trial is described by small specification objects. A table of
correspondence between the old and new functions will be added before the
release.

## Specification functions

* Priors: `prior_beta()`, `prior_dirichlet()`, `prior_nix()`, `prior_niw()`
  and `prior_vague()`.
* Data summaries: `data_binary()`, `data_continuous()`, `data_binary2()` and
  `data_continuous2()`.
* External data with a power prior weight from 0 to 1: `borrow()`.
* Hypothetical control for uncontrolled studies: `hypothetical()`.
* Decision criteria: `criterion_posterior()` and `criterion_predictive()`.
  For all endpoints, the Go criterion is that the probability of the Go
  region is at least `gamma_go` and the NoGo criterion is that the
  probability of the NoGo region is at least `gamma_nogo`.
* `qdm_regions()` shows the numbering of the regions for two endpoints.

## Distribution functions

* `pbetadiff()` has the argument order `(q, a_t, b_t, a_c, b_c)`. The range
  of integration is restricted to the support and split at quantiles of both
  Beta distributions, so concentrated distributions are integrated
  accurately. The computation is done in C++.
* `pbetabinomdiff()` has the argument order
  `(q, m_t, a_t, b_t, m_c, a_c, b_c)` and computes the probability mass
  functions on the log scale, so large numbers of future subjects no longer
  cause overflow.
* `ptdiff()` replaces `ptdiff_NI()`, `ptdiff_MC()` and `ptdiff_MM()`. The
  numerical integration now covers the whole real line (version 0.1.0
  truncated the range, with errors of up to about 3e-4 for small degrees of
  freedom).
* `rdirichlet()` always returns a matrix.
* `cellprob_binary2()` replaces `getjointbin()`.

# BayesianQDM 0.1.0

## Resubmission

* Expanded the MCMC acronym to "Markov chain Monte Carlo (MCMC)" in DESCRIPTION.
* Added a reference to Kang, Yamaguchi, and Han (2026) in DESCRIPTION.
* Replaced `\dontrun{}` with `\donttest{}` in examples for `getgamma2bin()` and
  `getgamma2cont()`.

## Initial Release

Initial release providing a comprehensive Bayesian quantitative decision-making
framework for clinical trials with single and two-endpoint analyses for binary
and continuous outcomes.

### Decision Probability Functions

* `pbayesdecisionprob1bin()` - Go/NoGo/Gray decision probabilities for a single
  binary endpoint
* `pbayesdecisionprob1cont()` - Go/NoGo/Gray decision probabilities for a single
  continuous endpoint
* `pbayesdecisionprob2bin()` - Go/NoGo/Gray decision probabilities for two binary
  endpoints
* `pbayesdecisionprob2cont()` - Go/NoGo/Gray decision probabilities for two
  continuous endpoints

### S3 Methods for Decision Probability Classes

* `print.pbayesdecisionprob1bin()` - Print method for `pbayesdecisionprob1bin` objects
* `print.pbayesdecisionprob1cont()` - Print method for `pbayesdecisionprob1cont` objects
* `print.pbayesdecisionprob2bin()` - Print method for `pbayesdecisionprob2bin` objects
* `print.pbayesdecisionprob2cont()` - Print method for `pbayesdecisionprob2cont` objects
* `plot.pbayesdecisionprob1bin()` - Plot method for `pbayesdecisionprob1bin` objects;
  displays Go/NoGo/Gray decision probabilities as a line plot across treatment scenarios
* `plot.pbayesdecisionprob1cont()` - Plot method for `pbayesdecisionprob1cont` objects;
  displays Go/NoGo/Gray decision probabilities as a line plot across treatment scenarios
* `plot.pbayesdecisionprob2bin()` - Plot method for `pbayesdecisionprob2bin` objects;
  displays Go/NoGo/Gray decision probabilities as a tile or scatter plot over a
  two-dimensional treatment scenario grid
* `plot.pbayesdecisionprob2cont()` - Plot method for `pbayesdecisionprob2cont` objects;
  displays Go/NoGo/Gray decision probabilities as a tile or scatter plot over a
  two-dimensional treatment scenario grid

### S3 Methods for Optimal Threshold Search Classes

* `plot.getgamma1bin()` - Plot method for `getgamma1bin` objects; displays
  calibration curves of marginal Go and NoGo probabilities against the threshold
  grid for a single binary endpoint
* `plot.getgamma1cont()` - Plot method for `getgamma1cont` objects; displays
  calibration curves of marginal Go and NoGo probabilities against the threshold
  grid for a single continuous endpoint
* `plot.getgamma2bin()` - Plot method for `getgamma2bin` objects; displays
  calibration curves of marginal Go and NoGo probabilities against the threshold
  grid for two binary endpoints
* `plot.getgamma2cont()` - Plot method for `getgamma2cont` objects; displays
  calibration curves of marginal Go and NoGo probabilities against the threshold
  grid for two continuous endpoints

### Posterior and Predictive Probability Functions

* `pbayespostpred1bin()` - Posterior or predictive probability for a single binary
  endpoint
* `pbayespostpred1cont()` - Posterior or predictive probability for a single
  continuous endpoint
* `pbayespostpred2bin()` - Joint region probabilities for two binary endpoints
* `pbayespostpred2cont()` - Joint region probabilities for two continuous endpoints

### Optimal Threshold Search Functions

* `getgamma1bin()` - Find optimal Go/NoGo thresholds for a single binary endpoint
  by grid search over candidate gamma values
* `getgamma1cont()` - Find optimal Go/NoGo thresholds for a single continuous
  endpoint by grid search over candidate gamma values
* `getgamma2bin()` - Find optimal Go/NoGo thresholds for two binary endpoints by
  grid search over candidate gamma value pairs
* `getgamma2cont()` - Find optimal Go/NoGo thresholds for two continuous endpoints
  by grid search over candidate gamma value pairs

### Distribution Functions for Continuous Endpoints

* `ptdiff_NI()` - CDF for the difference of two t-distributions via numerical
  integration
* `ptdiff_MC()` - CDF for the difference of two t-distributions via Monte Carlo
  simulation
* `ptdiff_MM()` - CDF for the difference of two t-distributions via
  Moment-Matching approximation

### Distribution Functions for Binary Endpoints

* `pbetadiff()` - CDF for the difference of two independent Beta distributions
* `pbetabinomdiff()` - Beta-binomial posterior predictive probability

### Sampling Functions

* `rdirichlet()` - Random sampler for the Dirichlet distribution

### Utility Functions

* `getjointbin()` - Joint binary probability from marginals and a correlation
  parameter
* `allmultinom()` - Enumerate all multinomial outcome combinations

### Study Designs

* Controlled design
* Uncontrolled design (hypothetical control)
* External design with power priors

### Prior Distributions

* Binary endpoints: Beta and Dirichlet conjugate priors
* Single continuous endpoints: Normal-Inverse-Chi-squared conjugate prior and
  vague (Jeffreys) prior
* Two continuous endpoints: Normal-Inverse-Wishart conjugate prior and vague prior

### Calculation Methods for Continuous Endpoints

* NI (Numerical Integration): exact computation via adaptive quadrature
* MC (Monte Carlo): simulation-based estimation
* MM (Moment-Matching): closed-form approximation; fully vectorised and
  recommended for large-scale simulation studies
