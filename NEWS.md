# BayesianQDM 0.1.0

## Initial Release

The initial release of BayesianQDM provides a comprehensive framework for Bayesian decision-making in clinical trials with support for both binary and continuous endpoints.

### Core Functions

#### Binary Endpoints

* `BayesPostPredBinary()` - Posterior and posterior predictive probability calculation
* `BayesDecisionProbBinary()` - Go/NoGo/Gray decision probability framework
* `pBetadiff()` - Beta distribution differences
* `pBetaBinomdiff()` - Beta-Binomial distribution differences

#### Continuous Endpoints

* `BayesPostPredContinuous()` - Posterior and posterior predictive probability calculation
* `BayesDecisionProbContinuous()` - Go/NoGo/Gray decision probability framework
* Distribution difference functions:
  - `pNIdifft()` - Numerical integration method
  - `pMCdifft()` - Monte Carlo simulation method
  - `pWSdifft()` - Welch-Satterthwaite approximation method
  - `pMCMCdiff()` - MCMC sampling method with power prior for external data

#### Utility Functions

* `AppellsF1()` - Appell's hypergeometric function F1 for numerical calculations

### Study Designs

* **Controlled design** - Standard randomized controlled trials
* **Uncontrolled design** - Single-arm studies with historical controls
* **External control design** - Power prior incorporation of historical data

### Prior Distributions

* **Binary endpoints** - Beta priors with flexible parameterization
* **Continuous endpoints** - Normal-Inverse-Chi-squared conjugate priors and vague priors

### Calculation Methods

* **NI** (Numerical Integration) - Exact calculation using convolution
* **WS** (Welch-Satterthwaite) - Fast approximation for unequal variances
* **MC** (Monte Carlo) - Simulation-based flexible approach
* **MCMC** - Markov Chain Monte Carlo sampling for external data incorporation

### Documentation

* Comprehensive function documentation with examples
* Three detailed vignettes:
  - Introduction to BayesianQDM
  - Binary endpoints analysis
  - Continuous endpoints analysis
* Complete test suite using testthat

### Dependencies

* Base R stats functions
* bayesDP package (suggested) for external data incorporation via MCMC
* Compatible with tidyverse ecosystem for data manipulation and visualization

### Quality Assurance

* Extensive test coverage for all core functions
* Input validation and informative error messages
* Numerical precision checks and method comparisons
* Operating characteristics validation

## Future Development

Planned features for future releases:
* Additional calculation methods for complex scenarios
* Enhanced visualization functions
* Adaptive design utilities
* Extended external data incorporation methods
* Performance optimizations for large-scale simulations
