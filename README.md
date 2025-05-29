# BayesianQDM <img src="man/figures/badge-BayesianQDM.png" align="right" height="139" /></a>

<!-- badges: start -->
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/BayesianQDM)](http://cran.r-project.org/package=BayesianQDM)
<!-- badges: end -->

## Overview

`BayesianQDM` is an R package that provides methods to calculate posterior probabilities and posterior predictive probabilities, and Go, NoGo and Gray probabilities for quantitative decision-making framework under bayesian paradigm.

For technical details about the methodology, please refer to Kang et al.(20XX).

## Installation

You can install the development version of BayesianQDM from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("gosukehommaEX/BayesianQDM")
```

## Usage

### Binary endpoint

``` r
# Load packages
library(BayesianQDM)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

# Calculate bayesian decision probabilities
results = tibble(
  design = rep(c('controlled', 'uncontrolled', 'external'), each = 2),
  prob = rep(c('posterior', 'predictive'), 3),
  theta.TV = c(0.5, NA, 0.5, NA, 0.5, NA),
  theta.MAV = c(0.15, NA, 0.15, NA, 0.15, NA),
  theta.NULL = c(NA, 0.5, NA, 0.5, NA, 0.5),
  gamma1 = c(0.8, 0.8, 0.8, 0.8, 0.6, 0.6),
  gamma2 = c(0.3, 0.3, 0.3, 0.3, 0.3, 0.3)
) %>%
  group_by_all() %>%
  reframe(
    GoNoGoGray = map(list(1), ~ {
      BayesDecisionProbBinary(
        prob = prob, design = design, theta.TV = theta.TV, theta.MAV = theta.MAV, theta.NULL = theta.NULL, gamma1 = gamma1, gamma2 = gamma2,
        pi1 = seq(0, 1, l = 101), pi2 = rep(0.2, 101), n1 = 12, n2 = 15, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5, z = 0.5,
        m1 = 12, m2 = 12, ne1 = 12, ne2 = 12, ye1 = 6, ye2 = 6, ae1 = 0.5, ae2 = 0.5
      )
    })
  ) %>%
  unnest(GoNoGoGray) %>%
  mutate(
    Gray = if_else(Gray < 0, 0, Gray)
  ) %>% 
  mutate(theta = pi1 - pi2)

# Display a figure showing bayesian decision making
figure = results %>%
  pivot_longer(
    cols = c(Go, NoGo, Gray), names_to = 'Decision', values_to = 'Prob'
  ) %>%
  mutate(
    design = factor(design, levels = c('controlled', 'uncontrolled', 'external')),
    prob = factor(prob, levels = c('posterior', 'predictive')),
    Decision = factor(Decision, levels = c('Go', 'Gray', 'NoGo'))
  ) %>%
  ggplot(aes(x = theta, y = Prob)) +
  facet_grid(
    prob ~ design
  ) +
  geom_line(aes(colour = Decision, linetype = Decision), linewidth = 1) +
  theme_bw() +
  scale_color_manual(
    values = c('Go' = '#658D1B', 'Gray' = '#939597', 'NoGo' = '#D91E49'),
    labels =  c('Go', 'Gray', 'NoGo')
  ) +
  scale_x_continuous(
    breaks = seq(0 - 0.2, 1 - 0.2, l = 6)
  ) +
  scale_y_continuous(
    breaks = seq(0, 1, l = 11)
  ) +
  labs(
    title = 'Probability of making Go/Gray/NoGo decision',
    x = expression(theta),
    y = 'Probability'
  ) +
  theme(
    text = element_text(size = 40),
    panel.spacing = unit(1.5, 'lines'),
    legend.background = element_rect(fill = 'white', color = 'black'),
    legend.key.width = unit(4, 'cm'),
    legend.text = element_text(size = 40),
    legend.title = element_blank(), 
    legend.position = 'bottom'
  )
```

<img src="man/figures/README-figure-binary.png" width="100%"/>

### Continuous endpoint

``` r
# Load packages
library(BayesianQDM)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggh4x)

# Calculate bayesian decision probabilities
results = tibble(
  CalcMethod = c('NI', 'MC', 'WS')
) %>%
  group_by_all() %>%
  reframe(
    design = rep(c('controlled', 'uncontrolled'), each = 4),
    prob = rep(rep(c('posterior', 'predictive'), each = 2), 2),
    prior = rep(c('N-Inv-Chisq', 'vague'), 4)
  ) %>%
  group_by_all() %>%
  mutate(
    muk = list(
      tibble(
        mu1 = seq(0, 8, by = 0.2),
        mu2 = 0
      )
    )
  ) %>%
  reframe(
    GoNoGoGray = map(list(1), ~ {
      muk[[1]] %>%
        group_by_all() %>%
        reframe(
          BayesDecisionProbContinuous(
            nsim = 10000, prob = prob, design = design, prior = prior, CalcMethod = CalcMethod,
            theta.TV = 2, theta.MAV = 0, theta.NULL = 0.5, nMC = 10000, gamma1 = 0.8, gamma2 = 0.3,
            n1 = 12, n2 = 12, m1 = 120, m2 = 120, kappa01 = 5, kappa02 = 5,
            nu01 = 5, nu02 = 5, mu01 = 5, mu02 = 5, sigma01 = sqrt(5), sigma02 = sqrt(5),
            mu1 = mu1, mu2 = mu2, sigma1 = 1, sigma2 = 1, r = 12, seed = 1
          )
        )
    })
  ) %>%
  unnest(GoNoGoGray) %>%
  mutate(
    theta = mu1 - mu2
  )

# Display a figure showing bayesian decision making
figure = results %>%
  pivot_longer(
    cols = c(Go, NoGo, Gray), names_to = 'Decision', values_to = 'Prob'
  ) %>%
  mutate(
    CalcMethod = factor(CalcMethod, levels = c('NI', 'MC', 'WS')),
    design = factor(design, levels = c('controlled', 'uncontrolled')),
    prob = factor(prob, levels = c('posterior', 'predictive')),
    prior = factor(prior, levels = c('N-Inv-Chisq', 'vague')),
    Decision = factor(Decision, levels = c('Go', 'Gray', 'NoGo'))
  ) %>%
  ggplot(aes(x = theta, y = Prob)) +
  facet_nested(
    prob ~ design + prior,
    nest_line = element_line(colour = 'black')
  ) +
  geom_line(aes(colour = Decision, linetype = CalcMethod), linewidth = 1) +
  theme_bw() +
  scale_color_manual(
    values = c('Go' = '#658D1B', 'Gray' = '#939597', 'NoGo' = '#D91E49'),
    labels =  c('Go', 'Gray', 'NoGo')
  ) +
  scale_x_continuous(
    #expand=c(0, 0),
    limits = c(0, 8),
    breaks = seq(0, 8, l = 9)
  ) +
  scale_y_continuous(
    #expand=c(0, 0),
    limits = c(0, 1),
    breaks = seq(0, 1, l = 11)
  ) +
  labs(
    title = 'Probability of making Go/Gray/NoGo decision',
    x = expression(theta),
    y = 'Probability'
  ) +
  theme(
    text = element_text(size = 40),
    legend.background = element_rect(fill = 'white', color = 'black'),
    legend.key.width = unit(2, 'cm'),
    legend.text = element_text(size = 40),
    legend.title = element_blank(),
    legend.position = 'bottom'
  )
```
<img src="man/figures/README-figure-continuous.png" width="100%"/>

### Compare three methods (i.e., NI, MC and WS) for calculating posterior probabilities for continuous endpoint

``` r
# Load packages
library(BayesianQDM)
library(dplyr)

# Calculate bayesian posterior probabilities
set.seed(1)
comparison.CalcMethod = as_tibble(
  expand.grid(
    q = c(0, 1, 2),
    mu.t1 = seq(0, 8, by = 2),
    mu.t2 = seq(0, 8, by = 2),
    sd.t1 = c(1, 5, 10),
    sd.t2 = c(1, 5, 10),
    nu.t1 = c(5, 10, 20),
    nu.t2 = c(5, 10, 20)
  )
) %>% 
  filter(mu.t1 >= mu.t2) %>% 
  group_by(q, nu.t1, nu.t2) %>% 
  mutate(
    NI = pNIdifft(
      q = unique(q), mu.t1 = mu.t1, mu.t2 = mu.t2, 
      sd.t1 = sd.t1, sd.t2 = sd.t2, 
      nu.t1 = unique(nu.t1), nu.t2 = unique(nu.t2)
    ),
    MC = pMCdifft(
      nMC = 1e+5, q = unique(q), mu.t1 = mu.t1, mu.t2 = mu.t2, 
      sd.t1 = sd.t1, sd.t2 = sd.t2, 
      nu.t1 = unique(nu.t1), nu.t2 = unique(nu.t2)
    ),
    WS = pWSdifft(
      q = unique(q), mu.t1 = mu.t1, mu.t2 = mu.t2, 
      sd.t1 = sd.t1, sd.t2 = sd.t2, 
      nu.t1 = unique(nu.t1), nu.t2 = unique(nu.t2)
    )
  ) %>% 
  ungroup()

# Summarizing comparison of three methods
MC.vs.NI = comparison.CalcMethod %>%
  mutate(
    abs.diff = abs(MC - NI)
  ) %>% 
  filter(
    abs.diff == max(abs.diff)
  ) %>% 
  pull(abs.diff)
MC.vs.WS = comparison.CalcMethod %>%
  mutate(
    abs.diff = abs(MC - WS)
  ) %>% 
  filter(
    abs.diff == max(abs.diff)
  ) %>% 
  pull(abs.diff)
NI.vs.WS = comparison.CalcMethod %>%
  mutate(
    abs.diff = abs(NI - WS)
  ) %>% 
  filter(
    abs.diff == max(abs.diff)
  ) %>% 
  pull(abs.diff)
print(tibble(MC.vs.NI, MC.vs.WS, NI.vs.WS))

# # A tibble: 1 × 3
#   MC.vs.NI MC.vs.WS NI.vs.WS
#      <dbl>    <dbl>    <dbl>
# 1  0.00502   0.0261   0.0248
```

## References

Kang et al.(20XX)). Title
