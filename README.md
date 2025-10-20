# BayesianQDM <img src="man/figures/badge-BayesianQDM.png" align="right" height="139" alt="BayesianQDM package logo" />

<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/BayesianQDM)](http://cran.r-project.org/package=BayesianQDM)
[![R-CMD-check](https://github.com/gosukehommaEX/BayesianQDM/workflows/R-CMD-check/badge.svg)](https://github.com/gosukehommaEX/BayesianQDM/actions)
[![Codecov test coverage](https://codecov.io/gh/gosukehommaEX/BayesianQDM/branch/main/graph/badge.svg)](https://codecov.io/gh/gosukehommaEX/BayesianQDM?branch=main)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

## Overview
  
`BayesianQDM` is an R package that provides comprehensive methods for Bayesian quantitative decision-making in clinical trials. The package supports both binary and continuous endpoints with various study designs including controlled, uncontrolled, and external control designs.

### Key Features

- **Endpoint Types**: Binary (response/non-response) and continuous (biomarker changes, scores)
- **Study Designs**: Controlled, uncontrolled, and external control designs
- **Decision Framework**: Go/NoGo/Gray zones for evidence-based decision making
- **Calculation Methods**: Multiple approaches optimized for accuracy and speed
- **Prior Integration**: Flexible prior specifications including efficient power priors for external data

### Decision Framework

The package implements a three-zone Bayesian decision framework:
  
- **Go**: Sufficient evidence to proceed to next phase (high posterior probability of efficacy)
- **NoGo**: Insufficient evidence to proceed (low posterior probability of efficacy)  
- **Gray**: Inconclusive evidence (intermediate posterior probability requiring additional data)

For technical details about the methodology, please refer to Kang et al. (20XX).

## Installation

### From CRAN (when available)
``` r
install.packages("BayesianQDM")
```

### Development Version
``` r
# Install from GitHub
devtools::install_github("gosukehommaEX/BayesianQDM")

# Or with dependencies
devtools::install_github("gosukehommaEX/BayesianQDM", dependencies = TRUE)
```

## Quick Start

### Binary Endpoints Example

``` r
library(BayesianQDM)

# Calculate Go/NoGo/Gray probabilities for binary endpoint
result_binary <- pGNGsinglebinary(
  prob = 'posterior', design = 'controlled', 
  theta.TV = 0.3, theta.MAV = 0.1, theta.NULL = NULL,
  gamma1 = 0.8, gamma2 = 0.2,
  pi1 = c(0.3, 0.5, 0.7), pi2 = rep(0.2, 3), 
  n1 = 15, n2 = 15,
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  z = NULL, m1 = NULL, m2 = NULL,
  ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
)

print(result_binary)
```

### Continuous Endpoints Example

``` r
# Calculate decision probabilities for continuous endpoint  
result_continuous <- pGNGsinglecontinuous(
  nsim = 100, prob = 'posterior', design = 'controlled', 
  prior = 'vague', CalcMethod = 'WS',
  theta.TV = 1.5, theta.MAV = 0.5, theta.NULL = NULL,
  nMC = NULL, gamma1 = 0.8, gamma2 = 0.3,
  n1 = 12, n2 = 12, m1 = NULL, m2 = NULL,
  kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
  mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
  mu1 = 4.5, mu2 = 2.0, sigma1 = 1.5, sigma2 = 1.3,
  r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
  bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
  seed = 123
)

print(result_continuous)
```

### External Control Design Example

``` r
# External control design with power priors
result_power_prior <- pGNGsinglecontinuous(
  nsim = 100, prob = 'posterior', design = 'external', 
  prior = 'vague', CalcMethod = 'WS',
  theta.TV = 1.0, theta.MAV = 0.0, theta.NULL = NULL,
  nMC = NULL, gamma1 = 0.8, gamma2 = 0.2,
  n1 = 15, n2 = 15, m1 = NULL, m2 = NULL,
  kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
  mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
  mu1 = 3.0, mu2 = 1.5, sigma1 = 1.2, sigma2 = 1.1,
  r = NULL,
  ne1 = 25, ne2 = 25,                    # External sample sizes
  alpha01 = 0.6, alpha02 = 0.7,          # Power prior parameters
  bar.ye1 = 2.8, bar.ye2 = 1.4,          # Historical sample means
  se1 = 1.3, se2 = 1.2,                  # Historical sample SDs
  seed = 123
)

print(result_power_prior)
```

### Operating Characteristics
Evaluate decision framework performance:
``` r
# Simulate trials across different scenarios
scenarios <- expand.grid(
  true_effect = c(0, 0.5, 1.0, 1.5, 2.0),
  sample_size = c(15, 25, 35)
)

# Calculate Go probabilities for each scenario
# (implementation details in vignettes)
```

## Documentation

### 📚 Comprehensive Vignettes

The package includes detailed vignettes with practical examples:
  
- **[Introduction to BayesianQDM](https://gosukehommaEX.github.io/BayesianQDM/articles/BayesianQDM.html)**: Overview and quick start guide
- **[Binary Endpoints](https://gosukehommaEX.github.io/BayesianQDM/articles/binary-endpoints.html)**: Detailed examples for binary outcome analysis  
- **[Continuous Endpoints](https://gosukehommaEX.github.io/BayesianQDM/articles/continuous-endpoints.html)**: Comprehensive guide for continuous outcome analysis

### 📖 Function Documentation

- **[Full Reference](https://gosukehommaEX.github.io/BayesianQDM/reference/index.html)**: Complete function documentation with examples

### 🌐 Package Website

Visit our [pkgdown website](https://gosukehommaEX.github.io/BayesianQDM/) for:
- Interactive documentation
- Downloadable examples
- Method comparisons
- Best practice guides

### Local Access

Access vignettes locally after installation:
``` r
# View available vignettes
vignette(package = "BayesianQDM")

# Open specific vignettes
vignette("BayesianQDM")
vignette("binary-endpoints", package = "BayesianQDM")  
vignette("continuous-endpoints", package = "BayesianQDM")
```

## Core Functions

### Decision Making Functions
- `pGNGsinglebinary()` - Go/NoGo/Gray probabilities for binary endpoints
- `pGNGsinglecontinuous()` - Go/NoGo/Gray probabilities for continuous endpoints

### Probability Calculation Functions  
- `pPPsinglebinary()` - Posterior/predictive probabilities for binary endpoints
- `pPPsinglecontinuous()` - Posterior/predictive probabilities for continuous endpoints

### Utility Functions
- `p2betadiff()` - CDF for difference of two beta distributions
- `p2betabinomdiff()` - Beta-binomial posterior predictive probability
- `d2betadiff()` - PDF for difference of two beta distributions
- `pNI2tdiff()` - Numerical integration for t-distribution differences
- `pMC2tdiff()` - Monte Carlo for t-distribution differences
- `pWS2tdiff()` - Welch-Satterthwaite approximation for t-distribution differences
- `AppellsF1()` - Appell's hypergeometric function F1

## Performance Considerations

### Computational Efficiency
The package now uses **exact conjugate representation** for power priors, providing:

- **No MCMC required**: Closed-form computation for external data incorporation
- **Significant speed improvement**: Orders of magnitude faster than traditional power prior methods
- **Complete Bayesian rigor**: No approximation involved in power prior calculations

### Method Selection
- **NI method**: Most accurate, recommended for final analyses
- **WS method**: Fast approximation, excellent for large simulation studies
- **MC method**: Flexible, handles complex scenarios with adjustable precision

For large simulation studies, consider:
- Using WS method for efficient screening
- NI method for final analyses
- Parallel processing for Monte Carlo methods
- Leveraging exact power prior computation for external data scenarios

## Technical Innovation

### Exact Conjugate Power Priors
BayesianQDM implements a breakthrough approach to power prior computation:

- **Mathematical equivalence**: Power priors for normal data are exactly equivalent to Normal-Inverse-Chi-squared distributions
- **Closed-form solutions**: Eliminates computational burden of MCMC sampling
- **Theoretical rigor**: Preserves complete Bayesian framework without approximation
- **Practical benefits**: Enables efficient simulation studies and real-time analysis

This innovation makes power prior analysis computationally feasible for large-scale simulation studies while maintaining full statistical rigor.

## Citation

If you use BayesianQDM in your research, please cite:

```
Homma G, Yamaguchi Y (2025). BayesianQDM: Bayesian Quantitative
Decision-Making Framework for Binary and Continuous Endpoints.
R package version 0.1.0.
```

## References

Kang et al. (20XX). [Title]. [Journal].

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.html) file for details.
