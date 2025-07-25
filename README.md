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
- **Prior Integration**: Flexible prior specifications including power priors for external data

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

# Or with dependencies for external data features
devtools::install_github("gosukehommaEX/BayesianQDM", dependencies = TRUE)
```

### Optional Dependencies
For external data incorporation with MCMC sampling:
``` r
install.packages("bayesDP")
```

## Quick Start

### Binary Endpoints Example

``` r
library(BayesianQDM)

# Calculate Go/NoGo/Gray probabilities for binary endpoint
result_binary <- BayesDecisionProbBinary(
  prob = 'posterior', design = 'controlled', 
  theta.TV = 0.3, theta.MAV = 0.1, gamma1 = 0.8, gamma2 = 0.2,
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
result_continuous <- BayesDecisionProbContinuous(
  nsim = 100, prob = 'posterior', design = 'controlled', 
  prior = 'vague', CalcMethod = 'WS',
  theta.TV = 1.5, theta.MAV = 0.5, theta.NULL = NULL,
  nMC = NULL, nMCMCsample = NULL,
  gamma1 = 0.8, gamma2 = 0.3,
  n1 = 12, n2 = 12, m1 = NULL, m2 = NULL,
  mu1 = 4.5, mu2 = 2.0, sigma1 = 1.5, sigma2 = 1.3,
  seed = 123
)

print(result_continuous)
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
- `BayesDecisionProbBinary()` - Go/NoGo/Gray probabilities for binary endpoints
- `BayesDecisionProbContinuous()` - Go/NoGo/Gray probabilities for continuous endpoints

### Probability Calculation Functions  
- `BayesPostPredBinary()` - Posterior/predictive probabilities for binary data
- `BayesPostPredContinuous()` - Posterior/predictive probabilities for continuous data

### Distribution Functions
- `pBetadiff()`, `pBetaBinomdiff()` - Beta distribution differences
- `pNIdifft()`, `pWSdifft()`, `pMCdifft()`, `pMCMCdiff()` - t-distribution differences  
- `AppellsF1()` - Appell's hypergeometric function

## Calculation Methods

### For Continuous Endpoints

| Method | Description | Use Case |
|--------|-------------|----------|
| **NI** | Numerical Integration | Most accurate, recommended for final analyses |
| **WS** | Welch-Satterthwaite | Fast approximation, good for simulations |
| **MC** | Monte Carlo | Flexible, handles complex scenarios |
| **MCMC** | Markov Chain Monte Carlo | External data incorporation with power priors |

### Method Comparison Example

``` r
# Compare calculation methods
mu1 <- 3.5; mu2 <- 1.8; sd1 <- 1.3; sd2 <- 1.1; nu1 <- 14; nu2 <- 16

# Numerical integration (exact)
prob_ni <- pNIdifft(q = 1.5, mu.t1 = mu1, mu.t2 = mu2, 
                   sd.t1 = sd1, sd.t2 = sd2, nu.t1 = nu1, nu.t2 = nu2)

# Welch-Satterthwaite approximation (fast)  
prob_ws <- pWSdifft(q = 1.5, mu.t1 = mu1, mu.t2 = mu2,
                   sd.t1 = sd1, sd.t2 = sd2, nu.t1 = nu1, nu.t2 = nu2)

cat("NI:", round(prob_ni, 4), "WS:", round(prob_ws, 4), 
    "Diff:", round(abs(prob_ni - prob_ws), 4))
```

## Study Design Examples

### Controlled Design
Standard two-arm randomized controlled trial.

### Uncontrolled Design  
Single-arm study with comparison to historical control.

### External Control Design
Incorporating historical data using power priors.

## Prior Distributions

### Binary Endpoints
- **Beta priors**: Beta(a, b) for response probabilities
- **Common choices**: Beta(0.5, 0.5) for Jeffreys prior, Beta(1, 1) for uniform

### Continuous Endpoints  
- **Vague priors**: Non-informative approach letting data drive conclusions
- **Normal-Inverse-Chi-squared**: Conjugate priors for incorporating prior knowledge
- **Power priors**: For external data incorporation with controlled borrowing

## Advanced Features

### Power Prior Integration
Control the degree of borrowing from external data:
``` r
# α = 0: No borrowing
# α = 1: Full borrowing  
# 0 < α < 1: Partial borrowing

alpha_values <- c(0, 0.5, 1.0)
comparison <- sapply(alpha_values, function(alpha) {
  BayesPostPredBinary(
    prob = 'posterior', design = 'external', theta0 = 0.2,
    n1 = 15, n2 = 15, y1 = 8, y2 = 6,
    a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
    m1 = NULL, m2 = NULL,
    ne1 = 20, ne2 = 20, ye1 = 10, ye2 = 4,
    ae1 = alpha, ae2 = alpha
  )
})
```

### External Data with MCMC
For continuous endpoints with historical data:
``` r
# Example with external control data (requires bayesDP package)
if (requireNamespace("bayesDP", quietly = TRUE)) {
  result_mcmc <- pMCMCdiff(
    nMCMCsample = 5000,
    q = 1.5,
    mu.n1 = 4.2, mu.n2 = 2.8,
    sd.n1 = 1.5, sd.n2 = 1.4,
    n1 = 12, n2 = 12,
    ne1 = NULL, ne2 = 24,  # External control data only
    alpha01 = NULL, alpha02 = 0.5  # 50% borrowing from control
  )
  
  cat("Posterior probability with external data:", round(result_mcmc, 4))
}
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

## Performance Considerations

- **NI method**: Most accurate but moderate computational cost
- **WS method**: Fast approximation suitable for simulations
- **MC method**: Flexible but computationally intensive (use larger nMC for precision)
- **MCMC method**: Efficient for external data scenarios with power priors

For large simulation studies, consider:
- Using WS method for screening
- NI method for final analyses
- Parallel processing for Monte Carlo methods
- Caching results for repeated analyses

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

## Support

- **Documentation**: See package vignettes and function help
- **Issues**: Report bugs or request features on [GitHub](https://github.com/gosukehommaEX/BayesianQDM/issues)
- **Questions**: Use GitHub Discussions for methodology questions

---
