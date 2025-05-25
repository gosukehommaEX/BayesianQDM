# BayesianQDM

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
library(ggplot2)

# Calculate bayesian decision probabilities
results = BayesDecisionProbBinary(
  prob = 'predictive', design = 'controlled', theta.TV = NULL, theta.MAV = NULL, theta.NULL = 0, gamma1 = 0.9, gamma2 = 0.3,
  pi1 = seq(0, 1, l = 101), pi2 = rep(0.2, 101), n1 = 12, n2 = 12, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5, z = NULL,
  m1 = 30, m2 = 30, ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
)

# Display a figure showing bayesian decision making
figure = results %>% 
  as_tibble() %>% 
  mutate(theta = pi1 - pi2) %>% 
  pivot_longer(
    cols = c(Go, NoGo, Gray), names_to = 'Decision', values_to = 'Prob'
  ) %>% 
  mutate(
    Decision = factor(Decision, levels = c('Go', 'Gray', 'NoGo'))
  ) %>% 
  ggplot(aes(x = theta, y = Prob, group = Decision)) +
  geom_line(aes(color = Decision, linetype = Decision), linewidth = 2) +
  theme_bw() +
  scale_color_manual(
    values = c('Go' = '#658D1B', 'Gray' = '#939597', 'NoGo' = '#D91E49'),
    labels =  c('Go', 'Gray', 'NoGo')
  ) +
  scale_x_continuous(
    #expand=c(0, 0),
    limits = c(0 - 0.2, 1 - 0.2),
    breaks = seq(0 - 0.2, 1 - 0.2, l = 6)
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
    text = element_text(size = 30),
    legend.background = element_rect(fill = 'white', color = 'black'),
    legend.key.width = unit(4, 'cm'),
    legend.text = element_text(size = 30),
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
library(scales)
library(patchwork)

# Calculate bayesian decision probabilities
results = tibble(
  Approach = c('Convolution', 'WS.approx'),
  approx = c(FALSE, TRUE)
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
            nsim = 10000, prob = prob, design = design, prior = prior, approx = approx,
            theta.TV = 2, theta.MAV = 0, theta.NULL = 0.5, gamma1 = 0.8, gamma2 = 0.3,
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
    Approach = factor(Approach, levels = c('Convolution', 'WS.approx')),
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
  geom_line(aes(colour = Decision, linetype = Approach), linewidth = 1) +
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

## References

Kang et al.(20XX)). Title
