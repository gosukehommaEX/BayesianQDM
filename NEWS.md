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

## Designs and operations: one binary endpoint

* `design_binary()` specifies a study with one binary endpoint: sample sizes,
  Beta priors, external data for either group (`borrow()`), and a concurrent,
  external-only (`n_c = 0` with `external_c`) or hypothetical (`n_c = 0` with
  `hypothetical_c`) control group.
* `qdm_decide()` computes the probabilities of the Go and NoGo regions for
  observed data and returns the decision: Go, NoGo, Gray, or Miss (both
  criteria met).
* `qdm_oc()` computes the probabilities of the four decisions under
  scenarios for the true response rates by exact enumeration of the
  outcomes. Miss decisions are reported in their own column with a warning;
  the arguments `error_if_Miss` and `Gray_inc_Miss` of version 0.1.0 are
  removed.
* `qdm_calibrate()` chooses `gamma_go` and `gamma_nogo` as the smallest grid
  values for which the probabilities that the Go criterion is met under one
  scenario and that the NoGo criterion is met under another are at most
  their targets, and returns the design with these thresholds and its
  decision probabilities.
* `print()` methods for designs, decisions, operating characteristics and
  calibrations, and `plot()` methods (ggplot2) for operating characteristics
  and calibrations. Without a concurrent control group, the thresholds in
  the plot of operating characteristics are shifted by the mean of the
  distribution of the control group.
* Probabilities are no longer set to 0 below `.Machine$double.eps^0.25`. A
  probability that is below a threshold by less than 1e-10 is treated as
  equal to the threshold, so that exact ties are not decided by rounding
  errors.

## Designs and operations: one continuous endpoint

* `design_continuous()` specifies a study with one normally distributed
  endpoint: sample sizes, Normal-Inverse-chi-squared (`prior_nix()`) or vague
  (`prior_vague()`) priors, external data for either group (`borrow()`), a
  concurrent, external-only or hypothetical control group, and the method
  used to compute the probabilities (`"NI"`, `"MM"` or `"MC"`).
* External data enter through the conjugate representation of the power
  prior. With the vague prior, the posterior distribution of the mean has
  `a0 * n_e + n - 1` degrees of freedom, and the posterior scale uses the
  same denominator (version 0.1.0 used `a0 * n_e + n`, so the posterior did
  not reduce to the one without external data as `a0` tends to 0).
* The posterior predictive distribution of the mean of `m` future subjects
  is `t(nu_n, mu_n, sigma_n^2 * (1 / kappa_n + 1 / m))`. Version 0.1.0 used
  the squared scale `sigma_n^2 * (1 + 1 / kappa_n) / m`, which treats the `m`
  future subjects as independent draws from the predictive distribution of
  one subject and understates the uncertainty for `m` of 2 or more.
* A hypothetical control (`hypothetical(mean = mu0c, r = r)`) gives the
  control mean the distribution `t(nu_nt, mu0c, r * sigma_nt^2 / kappa_nt)`
  (Kang et al., 2026, Section 4.2.4). For the predictive criterion its
  squared scale is `r * sigma_nt^2 * (1 / kappa_nt + 1 / m_c)`.
* `qdm_oc()` and `qdm_calibrate()` simulate the sample means and standard
  deviations of `nsim` studies, with the same random numbers for all
  scenarios. With `method = "MC"`, the probabilities of the Go and NoGo
  regions are estimated from the same draws, so that they add up to 1 for
  the predictive criterion.
* `qdm_decide()` has a `seed` argument for designs with `method = "MC"`.

## Designs and operations: two binary endpoints

* `design_binary2()` specifies a study with two binary endpoints: sample
  sizes, Dirichlet priors for the probabilities of the four combinations of
  responses (`prior_dirichlet()`), external data for either group
  (`borrow(data_binary2(...), a0)`), a concurrent, external-only or
  hypothetical control group, the number of Monte Carlo draws `nMC`, and
  whether the operating characteristics are computed by enumerating all
  outcomes (`enumerate = TRUE`, the default) or by simulation.
* For two endpoints, `criterion_posterior()` and `criterion_predictive()`
  have the argument `rule`. With `rule = "joint"` (the default), the
  criteria are on the probabilities of unions of regions (`go_regions` and
  `nogo_regions`, numbered as in `qdm_regions()`). With `rule = "and"` or
  `rule = "or"`, the criteria are on each endpoint separately (the
  probability that the effect exceeds its TV, or is at most its MAV), with
  one threshold for both endpoints or one per endpoint.
* The posterior predictive probabilities of the regions are computed
  exactly from the joint distribution of the numbers of future responders
  on the two endpoints (Dirichlet-multinomial distribution) instead of by
  Monte Carlo simulation. A future effect equal to the null value is not
  counted as larger; version 0.1.0 compared differences of rounded
  proportions, so some exact ties were counted as larger.
* The posterior probabilities of the regions are estimated with `nMC`
  draws (default 10000) from the Dirichlet posterior distributions, computed
  in C++. The probabilities of all regions, and of the Go and NoGo regions,
  are estimated from the same draws, and `qdm_decide()` reports all of them.
  `qdm_decide()`, `qdm_oc()` and `qdm_calibrate()` take `seed` for these
  draws.
* `qdm_oc()` takes scenarios with the columns `pi_t1`, `pi_t2`, `rho_t`,
  `pi_c1`, `pi_c2` and `rho_c`. With enumeration, the region probabilities
  are computed once for every pair of treatment and control outcomes, with
  the same Monte Carlo draws for all pairs, and weighted by the multinomial
  probabilities under each scenario. With simulation, all scenarios use the
  same uniform random numbers. Enumeration is limited to 1e7 pairs of
  outcomes (for example, 24 subjects in each group). A scenario whose
  correlation is outside the feasible range for its response rates is not
  evaluable: its probabilities are `NA` with a warning (version 0.1.0
  stopped with an error in its Monte Carlo path). `qdm_calibrate()` requires
  feasible scenarios.
* `plot()` of operating characteristics for two endpoints draws a heat map
  of the probability of each decision over the effects on the two endpoints
  and marks scenarios that are not evaluable with "NE".

## Designs and operations: two continuous endpoints

* `design_continuous2()` specifies a study with two continuous endpoints:
  sample sizes, Normal-Inverse-Wishart (`prior_niw()`) or vague
  (`prior_vague()`) priors, external data for either group
  (`borrow(data_continuous2(...), a0)`), a concurrent, external-only or
  hypothetical (`hypothetical(mean = c(m1, m2), r = r)`) control group, the
  decision criterion for two endpoints (with `rule` as for two binary
  endpoints), and the method used to compute the probabilities (`"MM"` or
  `"MC"`).
* The posterior parameters are updated with the external data (power prior)
  and then with the current data; the vague prior corresponds to
  `kappa0 = 0`, `nu0 = -1` and `Lambda0 = 0`.
* The posterior predictive distribution of the mean of `m` future subjects
  has the scale matrix `(1 / kappa_n + 1 / m) * Lambda_n / (nu_n - 1)`.
  Version 0.1.0 used `(1 + kappa_n) / kappa_n * Lambda_n / (nu_n - 1) / m`,
  which treats the `m` future subjects as independent draws from the
  predictive distribution of one subject.
* With `method = "MM"`, the moment-matching approximation of the difference
  of two bivariate t distributions is used, and its rectangle probabilities
  are computed by numerical integration in C++ for real degrees of freedom
  (version 0.1.0 rounded the degrees of freedom to an integer for
  `mvtnorm::pmvt()`). With `method = "MC"`, the probabilities of all
  regions are estimated from the same draws, and every simulated study uses
  its own draws (version 0.1.0 used the same draws for all simulated
  studies, so the Monte Carlo error did not average out).
* `qdm_oc()` takes scenarios with the columns `mu_t1`, `mu_t2`, `sd_t1`,
  `sd_t2`, `rho_t`, `mu_c1`, `mu_c2`, `sd_c1`, `sd_c2` and `rho_c`, and
  simulates the sample means and the sums of squares and cross-products
  (Bartlett decomposition of the Wishart distribution) of `nsim` studies
  with the same random numbers for all scenarios.

## Distribution functions

* `pbetadiff()` has the argument order `(q, a_t, b_t, a_c, b_c)`. The range
  of integration is restricted to the support and split at quantiles of both
  Beta distributions, so concentrated distributions are integrated
  accurately. When a shape parameter of the control group is below 1, the
  integral near the unbounded end of its density is computed after a change
  of variables that removes the singularity. The computation is done in C++.
* `pbetabinomdiff()` has the argument order
  `(q, m_t, a_t, b_t, m_c, a_c, b_c)` and computes the probability mass
  functions on the log scale, so large numbers of future subjects no longer
  cause overflow. Values of `q` far outside [-1, 1] are handled exactly.
* `ptdiff()` replaces `ptdiff_NI()`, `ptdiff_MC()` and `ptdiff_MM()`. The
  numerical integration now covers the whole real line (version 0.1.0
  truncated the range, with errors of up to about 3e-4 for small degrees of
  freedom). It integrates over the group with the smaller scale, splits the
  range around the peak of the density and the change of the other factor,
  and integrates the parts far from the origin on a logarithmic scale, so
  that scale ratios from 1e-9 to 1e3 and degrees of freedom from 0.5 to
  `Inf` are handled.
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
