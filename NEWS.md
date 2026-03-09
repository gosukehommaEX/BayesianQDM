# BayesianQDM 0.1.0

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
