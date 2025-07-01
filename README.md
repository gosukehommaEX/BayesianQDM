# BayesianQDM <img src="man/figures/badge-BayesianQDM.png" align="right" height="139" /></a>

<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/BayesianQDM)](http://cran.r-project.org/package=BayesianQDM)
[![R-CMD-check](https://github.com/gosukehommaEX/BayesianQDM/workflows/R-CMD-check/badge.svg)](https://github.com/gosukehommaEX/BayesianQDM/actions)
[![Codecov test coverage](https://codecov.io/gh/gosukehommaEX/BayesianQDM/branch/main/graph/badge.svg)](https://codecov.io/gh/gosukehommaEX/BayesianQDM?branch=main)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/gosukehommaEX/BayesianQDM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/gosukehommaEX/BayesianQDM/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/gosukehommaEX/BayesianQDM/graph/badge.svg)](https://app.codecov.io/gh/gosukehommaEX/BayesianQDM)
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

You can install the development version of BayesianQDM from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("gosukehommaEX/BayesianQDM")
```

Once submitted to CRAN, you will be able to install the stable version with:

``` r
install.packages("BayesianQDM")
```

## Quick Start

### Binary Endpoint Example

``` r
library(BayesianQDM)

# Calculate Go/NoGo/Gray probabilities for binary endpoint
result_binary <- BayesDecisionProbBinary(
  prob = 'posterior', design = 'controlled', 
  theta.TV = 0.3, theta.MAV = 0.1, 
  gamma1 = 0.8, gamma2 = 0.2,
  pi1 = c(0.2, 0.4, 0.6, 0.8), pi2 = rep(0.2, 4), 
  n1 = 20, n2 = 20,
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5
)

print(result_binary)
```

### Continuous Endpoint Example

``` r
# Calculate Go/NoGo/Gray probabilities for continuous endpoint  
result_continuous <- BayesDecisionProbContinuous(
  nsim = 1000, prob = 'posterior', design = 'controlled', 
  prior = 'vague', CalcMethod = 'NI',
  theta.TV = 1.5, theta.MAV = 0.5, 
  gamma1 = 0.8, gamma2 = 0.3,
  n1 = 15, n2 = 15, mu1 = 3.0, mu2 = 1.2,
  sigma1 = 1.5, sigma2 = 1.5, seed = 123
)

print(result_continuous)
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
- `pNIdifft()`, `pWSdifft()`, `pMCdifft()` - t-distribution differences  
- `AppellsF1()` - Appell's hypergeometric function

## Calculation Methods

### For Continuous Endpoints

| Method | Description | Use Case |
|--------|-------------|----------|
| **NI** | Numerical Integration | Most accurate, recommended for final analyses |
| **WS** | Welch-Satterthwaite | Fast approximation, good for simulations |
| **MC** | Monte Carlo | Flexible, handles complex scenarios |
| **INLA** | Integrated Nested Laplace | External data incorporation |

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

``` r
# External control example
external_result <- BayesPostPredBinary(
  prob = 'posterior', design = 'external', theta0 = 0.15,
  n1 = 20, n2 = 20, y1 = 12, y2 = 8,
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  ne1 = 30, ne2 = 30, ye1 = 15, ye2 = 6,  
  ae1 = 0.5, ae2 = 0.5  # Power prior weights
)
```

## Visualization Example

``` r
library(dplyr)
library(tidyr) 
library(ggplot2

# Create decision probability plot
results <- BayesDecisionProbBinary(
  prob = 'posterior', design = 'controlled',
  theta.TV = 0.3, theta.MAV = 0.1, gamma1 = 0.8, gamma2 = 0.2,
  pi1 = seq(0.1, 0.9, by = 0.05), pi2 = rep(0.2, length(seq(0.1, 0.9, by = 0.05))),
  n1 = 20, n2 = 20, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5
) %>%
  mutate(theta = pi1 - pi2)

# Plot decision probabilities
results %>%
  pivot_longer(cols = c(Go, NoGo, Gray), names_to = 'Decision', values_to = 'Probability') %>%
  mutate(Decision = factor(Decision, levels = c('Go', 'Gray', 'NoGo'))) %>%
  ggplot(aes(x = theta, y = Probability)) +
  geom_line(aes(colour = Decision, linetype = Decision), linewidth = 1.2) +
  scale_color_manual(values = c('Go' = '#2E8B57', 'Gray' = '#808080', 'NoGo' = '#DC143C')) +
  theme_bw() +
  labs(title = 'Go/Gray/NoGo Decision Probabilities',
       x = 'Treatment Effect (π₁ - π₂)', y = 'Probability')
```

## Documentation

The package includes comprehensive documentation:

### Vignettes
- **Introduction to BayesianQDM**: Overview and quick start guide
- **Binary Endpoints**: Detailed examples for binary outcome analysis  
- **Continuous Endpoints**: Comprehensive guide for continuous outcome analysis

Access vignettes with:
``` r
vignette("BayesianQDM")
vignette("binary-endpoints", package = "BayesianQDM")  
vignette("continuous-endpoints", package = "BayesianQDM")
```

### Function Help
``` r
?BayesDecisionProbBinary
?BayesDecisionProbContinuous
?BayesPostPredBinary
?BayesPostPredContinuous
```

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
    ne1 = 20, ne2 = 20, ye1 = 10, ye2 = 4,
    ae1 = alpha, ae2 = alpha
  )
})
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
- **INLA method**: Efficient for external data scenarios

For large simulation studies, consider:
- Using WS method for screening
- NI method for final analyses
- Parallel processing for Monte Carlo approaches

## Common Use Cases

### Phase II Proof-of-Concept
- **Objective**: Determine if treatment shows efficacy signal
- **Typical settings**: θ_TV = 1.5 × MCID, γ₁ = 0.8

### Dose-Finding Studies
- **Objective**: Select optimal dose for Phase III
- **Considerations**: Multiple comparisons, dose-response modeling

### Biomarker-Driven Trials  
- **Objective**: Evaluate treatment in biomarker-defined populations
- **Challenges**: Smaller samples, incorporation of biomarker information

## Best Practices

### Parameter Selection
- **Thresholds**: Align with clinical meaningfulness (e.g., MCID)
- **Decision criteria**: Balance Type I/II errors for study objectives  
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
- **parallel**: High-performance computing for large simulations

## Contributing

Contributions are welcome! Please see our contributing guidelines and code of conduct. For bugs or feature requests, please open an issue on GitHub.

## Citation

If you use BayesianQDM in your research, please cite:

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

---

**Note**: This package is designed for research and educational purposes. For regulatory submissions, please ensure compliance with applicable guidelines and consider consulting with biostatisticians experienced in Bayesian methods.
