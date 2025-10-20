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

### Development version from GitHub

``` r
# install.packages("devtools")
devtools::install_github("gosukehommaEX/BayesianQDM")
```

## Quick Start

### Binary Endpoints Example

``` r
library(BayesianQDM)

# Calculate decision probabilities for binary endpoint
result_binary <- pGNGsinglebinary(
  prob = 'posterior', 
  design = 'controlled', 
  theta.TV = 0.3, 
  theta.MAV = 0.1, 
  theta.NULL = NULL,
  gamma1 = 0.8, 
  gamma2 = 0.2,
  pi1 = c(0.2, 0.4, 0.6, 0.8), 
  pi2 = rep(0.2, 4), 
  n1 = 12, 
  n2 = 12,
  a1 = 0.5, 
  a2 = 0.5, 
  b1 = 0.5, 
  b2 = 0.5,
  z = NULL, 
  m1 = NULL, 
  m2 = NULL,
  ne1 = NULL, 
  ne2 = NULL, 
  ye1 = NULL, 
  ye2 = NULL, 
  ae1 = NULL, 
  ae2 = NULL
)

print(result_binary)

# Visualize results
library(ggplot2)
library(tidyr)

result_long <- pivot_longer(result_binary, 
                            cols = c(Go, Gray, NoGo), 
                            names_to = "Decision", 
                            values_to = "Probability")

ggplot(result_long, aes(x = pi1 - pi2, y = Probability, fill = Decision)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("Go" = "#2E7D32", "Gray" = "#BABABA", "NoGo" = "#D32F2F")) +
  labs(title = 'Decision Probabilities for Binary Endpoint',
       x = 'Treatment Effect (π1 - π2)',
       y = 'Probability') +
  theme_minimal()
```

### Continuous Endpoints Example

``` r
# Calculate decision probabilities for continuous endpoint  
result_continuous <- pGNGsinglecontinuous(
  nsim = 100,
  prob = 'posterior', 
  design = 'controlled', 
  prior = 'vague', 
  CalcMethod = 'WS',
  theta.TV = 1.5, 
  theta.MAV = 0.5, 
  theta.NULL = NULL,
  nMC = NULL, 
  gamma1 = 0.8, 
  gamma2 = 0.3,
  n1 = 15, 
  n2 = 15, 
  m1 = NULL, 
  m2 = NULL,
  kappa01 = NULL, 
  kappa02 = NULL, 
  nu01 = NULL, 
  nu02 = NULL,
  mu01 = NULL, 
  mu02 = NULL, 
  sigma01 = NULL, 
  sigma02 = NULL,
  mu1 = 3.0, 
  mu2 = 1.2, 
  sigma1 = 1.5, 
  sigma2 = 1.5,
  r = NULL, 
  ne1 = NULL, 
  ne2 = NULL, 
  alpha01 = NULL, 
  alpha02 = NULL,
  bar.ye1 = NULL,
  bar.ye2 = NULL,
  se1 = NULL,
  se2 = NULL,
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
Incorporating historical/external data through power priors for enhanced decision making.

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

``` r
result <- pPPsinglebinary(
  prob = 'posterior', design = 'controlled', theta0 = 0.15,
  n1 = 20, n2 = 20, y1 = 12, y2 = 8,
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  m1 = NULL, m2 = NULL,
  ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL,
  ae1 = NULL, ae2 = NULL
)
```

### Uncontrolled Design
Single-arm study with historical control.

``` r
result <- pPPsinglebinary(
    prob = 'posterior', design = 'uncontrolled', theta0 = 0.15,
    n1 = 20, n2 = 20, y1 = 12, y2 = 3.5,
    a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
    m1 = NULL, m2 = NULL,
    ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL,
    ae1 = NULL, ae2 = NULL
)
```

### External Control Design
Incorporating historical data using power priors.

``` r
result <- pPPsinglebinary(
  prob = 'posterior', design = 'external', theta0 = 0.15,
  n1 = 20, n2 = 20, y1 = 12, y2 = 8,
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  m1 = NULL, m2 = NULL,
  ne1 = 30, ne2 = 30, ye1 = 18, ye2 = 12,
  ae1 = 0.5, ae2 = 0.5
)
```

## Best Practices

### Parameter Selection
- **Thresholds**: Align with clinical meaningfulness (e.g., MCID)
- **Decision criteria**: Balance Type I/II errors with study objectives
- **Priors**: Document rationale and conduct sensitivity analyses

### Validation
- **Operating characteristics**: Evaluate across relevant scenarios
- **Sensitivity analysis**: Test robustness to assumptions
- **Method comparison**: Verify consistency across calculation approaches

### Reporting
- **Transparency**: Report all inputs (priors, thresholds, methods)
- **Uncertainty**: Include credible intervals and sensitivity analyses
- **Interpretation**: Provide clinical context for decision probabilities

## Integration with R Ecosystem

The package integrates seamlessly with:
- **tidyverse**: Data manipulation and visualization workflows
- **ggplot2**: Creating publication-quality decision probability plots
- **knitr/rmarkdown**: Reproducible reporting and documentation

## Citation

To cite BayesianQDM in publications, please use:

```
Homma, G., Yamaguchi, Y. (2025). BayesianQDM: Bayesian Quantitative 
Decision-Making Framework for Binary and Continuous Endpoints. 
R package version 0.1.0.
```

## References

Kang et al. (20XX). [Title]. [Journal].

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Support

- **Documentation**: See package vignettes and function help
- **Issues**: Report bugs or request features on [GitHub](https://github.com/gosukehommaEX/BayesianQDM/issues)
- **Questions**: Use GitHub Discussions for methodology questions
