## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  fig.alt = "Visualization of statistical methods and decision probabilities for continuous endpoints"
)

## ----setup--------------------------------------------------------------------
library(BayesianQDM)
library(dplyr)
library(tidyr)
library(ggplot2)

## ----basic-example------------------------------------------------------------
# Calculate posterior probability that treatment effect > 1.5
posterior_prob <- BayesPostPredContinuous(
  prob = 'posterior',
  design = 'controlled', 
  prior = 'vague',
  CalcMethod = 'WS',  # Fast method for vignette
  theta0 = 1.5,
  n1 = 12, n2 = 12,
  bar.y1 = 5.2, bar.y2 = 3.1,
  s1 = 1.8, s2 = 1.6
)

cat("Posterior probability P(μ₁ - μ₂ > 1.5):", round(posterior_prob, 4))

## ----method-comparison--------------------------------------------------------
# Set common parameters
params <- list(
  prob = 'posterior',
  design = 'controlled',
  prior = 'vague', 
  theta0 = 1,
  n1 = 15, n2 = 15,
  bar.y1 = 4.2, bar.y2 = 2.8,
  s1 = 1.5, s2 = 1.4
)

# Compare NI and WS methods (skip MC and MCMC for vignette speed)
ni_result <- do.call(BayesPostPredContinuous, c(params, list(CalcMethod = 'NI')))
ws_result <- do.call(BayesPostPredContinuous, c(params, list(CalcMethod = 'WS')))

# Results comparison
results_df <- data.frame(
  Method = c("Numerical Integration", "Welch-Satterthwaite"),
  Probability = c(ni_result, ws_result),
  Speed = c("Slower", "Faster"),
  Accuracy = c("Exact", "Approximation")
)

knitr::kable(results_df, digits = 4, caption = "Method Comparison")

## ----decision-framework-------------------------------------------------------
# Calculate Go/NoGo/Gray probabilities
decision_result <- BayesDecisionProbContinuous(
  nsim = 30,  # Reduced for vignette speed
  prob = 'posterior',
  design = 'controlled',
  prior = 'vague',
  CalcMethod = 'WS',  # Fast method
  theta.TV = 2.0,   # Target value for Go
  theta.MAV = 0.5,  # Minimum acceptable value for NoGo
  theta.NULL = NULL,
  gamma1 = 0.8,     # Go threshold
  gamma2 = 0.3,     # NoGo threshold
  n1 = 12, n2 = 12,
  mu1 = 4.5, mu2 = 2.0,
  sigma1 = 1.5, sigma2 = 1.3,
  seed = 123
)

print(decision_result)

## ----decision-visualization---------------------------------------------------
# Create operating characteristics
scenarios <- expand.grid(
  true_diff = c(0, 0.5, 1.0, 1.5, 2.0, 2.5),
  sample_size = c(10, 15, 20)
) %>%
  rowwise() %>%
  mutate(
    decision = list(BayesDecisionProbContinuous(
      nsim = 20,  # Reduced for vignette speed
      prob = 'posterior',
      design = 'controlled',
      prior = 'vague',
      CalcMethod = 'WS',
      theta.TV = 1.5,
      theta.MAV = 0.5,
      theta.NULL = NULL,
      gamma1 = 0.8,
      gamma2 = 0.3,
      n1 = sample_size, n2 = sample_size,
      mu1 = 2.0 + true_diff, mu2 = 2.0,
      sigma1 = 1.5, sigma2 = 1.5,
      seed = 123
    ))
  ) %>%
  unnest_wider(decision)

# Plot operating characteristics  
scenarios_long <- scenarios %>%
  select(true_diff, sample_size, Go, NoGo, Gray) %>%
  pivot_longer(cols = c(Go, NoGo, Gray), names_to = "Decision", values_to = "Probability")

ggplot(scenarios_long, aes(x = true_diff, y = Probability, color = Decision)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~paste("N =", sample_size, "per arm")) +
  scale_color_manual(values = c('Go' = '#2E8B57', 'Gray' = '#808080', 'NoGo' = '#DC143C')) +
  labs(
    title = "Operating Characteristics for Continuous Endpoints",
    x = "True Treatment Difference",
    y = "Decision Probability",
    caption = "Note: Reduced simulation size for vignette speed"
  ) +
  theme_bw()

## ----informative-priors-------------------------------------------------------
# Example with informative priors
informative_result <- BayesPostPredContinuous(
  prob = 'posterior',
  design = 'controlled',
  prior = 'N-Inv-Chisq',
  CalcMethod = 'NI',
  theta0 = 1.0,
  n1 = 10, n2 = 10,
  bar.y1 = 3.2, bar.y2 = 2.1,
  s1 = 1.4, s2 = 1.2,
  # Prior parameters
  kappa01 = 2, kappa02 = 2,
  nu01 = 4, nu02 = 4, 
  mu01 = 3.0, mu02 = 2.0,
  sigma01 = 1.5, sigma02 = 1.5
)

cat("Posterior probability with informative priors:", round(informative_result, 4))

## ----external-data, eval=FALSE------------------------------------------------
# # Note: This example is not evaluated in the vignette due to computation time
# # In practice, use higher nMCMCsample values (e.g., 10000)
# 
# external_result <- pMCMCdiff(
#   nMCMCsample = 1000,  # Reduced for demonstration
#   q = 1.5,
#   mu.n1 = 4.2, mu.n2 = 2.8,
#   sd.n1 = 1.5, sd.n2 = 1.4,
#   n1 = 12, n2 = 12,
#   ne1 = 24, ne2 = 24,  # External data sample sizes
#   alpha01 = 0.5, alpha02 = 0.5  # Power prior parameters
# )
# 
# cat("Probability with external data:", round(external_result, 4))

