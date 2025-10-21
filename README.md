# BayesianQDM <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/BayesianQDM)](http://cran.r-project.org/package=BayesianQDM)
[![R-CMD-check](https://github.com/gosukehommaEX/BayesianQDM/workflows/R-CMD-check/badge.svg)](https://github.com/gosukehommaEX/BayesianQDM/actions)
[![Codecov test coverage](https://codecov.io/gh/gosukehommaEX/BayesianQDM/branch/main/graph/badge.svg)](https://codecov.io/gh/gosukehommaEX/BayesianQDM?branch=main)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

## Overview

BayesianQDM provides a comprehensive framework for **Bayesian Quantitative Decision-Making** in clinical trials. The package enables researchers to make **Go/NoGo/Gray** decisions using posterior and predictive probabilities for both binary and continuous endpoints.

## Key Features

- ✅ **Multiple Study Designs**: Controlled, uncontrolled, and external control designs
- 📊 **Flexible Endpoints**: Support for both binary and continuous outcomes
- 🔢 **Multiple Calculation Methods**: For continuous endpoints (Numerical Integration, Monte Carlo, Welch-Satterthwaite)
- 🎯 **Three-Zone Decision Framework**: Go/NoGo/Gray probability calculation
- 📈 **Power Priors**: Incorporation of historical/external data

## Installation

### From CRAN

```r
install.packages("BayesianQDM")
```

### Development Version

```r
# install.packages("devtools")
devtools::install_github("gosukehommaEX/BayesianQDM")
```

## Quick Start

### Binary Endpoint Example

```r
library(BayesianQDM)

# Calculate Go/NoGo/Gray probabilities for binary endpoint
result <- pGNGsinglebinary(
  prob = 'posterior', 
  design = 'controlled',
  theta.TV = 0.3, theta.MAV = 0.1, theta.NULL = NULL,
  gamma1 = 0.8, gamma2 = 0.2,
  pi1 = c(0.3, 0.5, 0.7), 
  pi2 = rep(0.2, 3),
  n1 = 15, n2 = 15,
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  z = NULL, m1 = NULL, m2 = NULL,
  ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
)

print(result)
```

### Continuous Endpoint Example

```r
# Calculate Go/NoGo/Gray probabilities for continuous endpoint
result <- pGNGsinglecontinuous(
  nsim = 100, 
  prob = 'posterior', 
  design = 'controlled', 
  prior = 'vague', 
  CalcMethod = 'NI',
  theta.TV = 1.5, theta.MAV = 0.5, theta.NULL = NULL,
  nMC = NULL, gamma1 = 0.8, gamma2 = 0.3,
  n1 = 15, n2 = 15, m1 = NULL, m2 = NULL,
  kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
  mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
  mu1 = 4.0, mu2 = 2.0, sigma1 = 1.5, sigma2 = 1.3,
  r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
  bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
  seed = 123
)

print(result)
```

## Documentation

### 📚 Vignettes

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

### Distribution Functions
- `p2betadiff()` - CDF for difference of two beta distributions
- `p2betabinomdiff()` - Beta-binomial posterior predictive probability
- `pNI2tdiff()` - Numerical integration for t-distribution differences
- `pMC2tdiff()` - Monte Carlo for t-distribution differences
- `pWS2tdiff()` - Welch-Satterthwaite approximation for t-distribution differences
- `d2betadiff()` - Density function for difference of two beta distributions

### Utility Functions
- `AppellsF1()` - Appell's hypergeometric function F1

## Study Designs

### Controlled Design
Standard randomized controlled trials with treatment and control groups.

### Uncontrolled Design
Single-arm studies using historical controls with beta priors (binary) or informative priors (continuous).

### External Control Design
Incorporating historical/external data through power priors for enhanced decision making. For binary endpoints, use `ne1`, `ne2`, `ye1`, `ye2`, `ae1`, `ae2` parameters. For continuous endpoints, use `ne1`, `ne2`, `alpha01`, `alpha02`, `bar.ye1`, `bar.ye2`, `se1`, `se2` parameters.

## Calculation Methods

### For Continuous Endpoints

| Method | Description | Use Case |
|--------|-------------|----------|
| **NI** | Numerical Integration | Most accurate, recommended for final analyses |
| **WS** | Welch-Satterthwaite | Fast approximation for unequal variances |
| **MC** | Monte Carlo | Flexible simulation-based approach |

### Method Comparison Example

``` r
# Compare calculation methods
mu1 <- 3.5; mu2 <- 1.8; sd1 <- 1.3; sd2 <- 1.1; nu1 <- 14; nu2 <- 16

# Numerical integration (exact)
prob_ni <- pNI2tdiff(q = 1.5, mu.t1 = mu1, mu.t2 = mu2, 
                     sd.t1 = sd1, sd.t2 = sd2, nu.t1 = nu1, nu.t2 = nu2)

# Welch-Satterthwaite approximation (fast)  
prob_ws <- pWS2tdiff(q = 1.5, mu.t1 = mu1, mu.t2 = mu2, 
                     sd.t1 = sd1, sd.t2 = sd2, nu.t1 = nu1, nu.t2 = nu2)

cat("NI:", round(prob_ni, 4), "WS:", round(prob_ws, 4), 
    "Diff:", round(abs(prob_ni - prob_ws), 4))
```

## Study Design Examples

### Controlled Design
Standard two-arm randomized controlled trial.

```r
# Binary endpoint
result <- pPPsinglebinary(
  prob = 'posterior', 
  design = 'controlled', 
  theta0 = 0.15,
  n1 = 12, n2 = 15, y1 = 7, y2 = 5, 
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  m1 = NULL, m2 = NULL,
  ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
)
```

### Uncontrolled Design
Single-arm study with historical control assumption.

```r
result <- pPPsinglebinary(
  prob = 'posterior', 
  design = 'uncontrolled', 
  theta0 = 0.15,
  n1 = 20, n2 = 20, y1 = 12, y2 = 3.5,
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  m1 = NULL, m2 = NULL,
  ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
)
```

### External Control Design
Incorporating historical data using power priors.

```r
result <- pPPsinglebinary(
  prob = 'posterior', 
  design = 'external', 
  theta0 = 0.15,
  n1 = 20, n2 = 20, y1 = 12, y2 = 8,
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  m1 = NULL, m2 = NULL,
  ne1 = 15, ne2 = 25, ye1 = 9, ye2 = 10, ae1 = 0.5, ae2 = 0.5
)
```

## References

For theoretical background and methodology details, please refer to the package vignettes and function documentation.

## Citation

To cite BayesianQDM in publications, please use:

```
Homma, G., Yamaguchi, Y. (2025). BayesianQDM: Bayesian Quantitative 
Decision-Making Framework for Binary and Continuous Endpoints. 
R package version 0.1.0.
```

## License

This package is licensed under GPL (>= 3).

## Authors

- Gosuke Homma
- Yusuke Yamaguchi

## Issues and Contributions

To report bugs or request features, please visit our [GitHub repository](https://github.com/gosukehommaEX/BayesianQDM/issues).
