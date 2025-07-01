# BayesianQDM 0.1.0

## Initial Release

This is the initial release of BayesianQDM, providing comprehensive methods for Bayesian quantitative decision-making in clinical trials.

### Core Features

#### Binary Endpoints
* `BayesDecisionProbBinary()` - Calculate Go/NoGo/Gray probabilities for binary outcomes
* `BayesPostPredBinary()` - Calculate posterior and posterior predictive probabilities
* `pBetadiff()` - Cumulative distribution function for difference of beta variables
* `pBetaBinomdiff()` - CDF for difference of beta-binomial variables
* `AppellsF1()` - Appell's first hypergeometric function implementation

#### Continuous Endpoints
* `BayesDecisionProbContinuous()` - Calculate Go/NoGo/Gray probabilities for continuous outcomes
* `BayesPostPredContinuous()` - Calculate posterior and posterior predictive probabilities
* `pNIdifft()` - Exact numerical integration method for t-distribution differences
* `pWSdifft()` - Welch-Satterthwaite approximation method
* `pMCdifft()` - Monte Carlo simulation method
* `pINLAdifft()` - INLA-based method for external data incorporation

### Study Designs
* **Controlled design** - Standard randomized controlled trials
* **Uncontrolled design** - Single-arm studies with historical controls
* **External control design** - Power prior incorporation of historical data

### Prior Distributions
* **Binary endpoints** - Beta priors with flexible parameterization
* **Continuous endpoints** - Normal-Inverse-Chi-squared conjugate priors and vague priors

### Calculation Methods
* **NI (Numerical Integration)** - Exact calculation using convolution
* **WS (Welch-Satterthwaite)** - Fast approximation for unequal variances
* **MC (Monte Carlo)** - Simulation-based flexible approach
* **INLA** - Integrated Nested Laplace Approximation for external data

### Documentation
* Comprehensive function documentation with examples
* Three detailed vignettes:
  - Introduction to BayesianQDM
  - Binary endpoints analysis
  - Continuous endpoints analysis
* Complete test suite using testthat

### Dependencies
* Base R stats functions
* INLA package (suggested) for external data incorporation
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
