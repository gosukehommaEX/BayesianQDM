## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  fig.alt = "Visualization of Go/NoGo/Gray decision probabilities for binary endpoints"
)

## ----setup--------------------------------------------------------------------
library(BayesianQDM)
library(dplyr)
library(tidyr)
library(ggplot2)

## ----posterior_basic----------------------------------------------------------
# Calculate posterior probability
posterior_prob <- BayesPostPredBinary(
  prob = 'posterior', 
  design = 'controlled', 
  theta0 = 0.15,
  n1 = 12, n2 = 15, y1 = 7, y2 = 9, 
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  m1 = NULL, m2 = NULL, # Not needed for posterior probability
  ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL # No external data
)

cat("Posterior probability that treatment effect > 0.15:", round(posterior_prob, 4))

## ----predictive_basic---------------------------------------------------------
# Calculate posterior predictive probability
predictive_prob <- BayesPostPredBinary(
  prob = 'predictive', 
  design = 'controlled', 
  theta0 = 0.1,
  n1 = 12, n2 = 15, y1 = 7, y2 = 7, 
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  m1 = 30, m2 = 30,    # Required for predictive probability
  ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL # No external data
)

cat("Posterior predictive probability that future trial effect > 0.1:", round(predictive_prob, 4))

## ----external_design----------------------------------------------------------
# Posterior probability with external control data
external_prob <- BayesPostPredBinary(
  prob = 'posterior', 
  design = 'external', 
  theta0 = 0.15,
  n1 = 12, n2 = 15, y1 = 7, y2 = 9, 
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  m1 = NULL, m2 = NULL, # Not needed for posterior probability
  ne1 = 12, ne2 = 12, ye1 = 6, ye2 = 6, 
  ae1 = 0.5, ae2 = 0.5
)

cat("Posterior probability with external control:", round(external_prob, 4))

## ----decision_posterior-------------------------------------------------------
# Calculate Go/NoGo/Gray probabilities using posterior probability
decision_posterior <- BayesDecisionProbBinary(
  prob = 'posterior', 
  design = 'controlled', 
  theta.TV = 0.3,    # Go threshold
  theta.MAV = 0.1,   # NoGo threshold
  theta.NULL = NULL, # Not needed for posterior probability
  gamma1 = 0.8,      # Minimum probability for Go
  gamma2 = 0.2,      # Maximum probability for NoGo
  pi1 = c(0.2, 0.4, 0.6, 0.8), 
  pi2 = rep(0.2, 4), 
  n1 = 12, n2 = 12,
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  z = NULL,          # Not needed for controlled design
  m1 = NULL, m2 = NULL, # Not needed for posterior probability
  ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL # No external data
)

print(decision_posterior)

## ----decision_predictive------------------------------------------------------
# Calculate Go/NoGo/Gray probabilities using posterior predictive probability
decision_predictive <- BayesDecisionProbBinary(
  prob = 'predictive', 
  design = 'controlled', 
  theta.TV = NULL, theta.MAV = NULL, # Not needed for predictive probability
  theta.NULL = 0.15,  # Threshold for future trial
  gamma1 = 0.9,       # Minimum probability for Go
  gamma2 = 0.3,       # Maximum probability for NoGo
  pi1 = c(0.2, 0.4, 0.6, 0.8), 
  pi2 = rep(0.2, 4), 
  n1 = 12, n2 = 12,
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  z = NULL,           # Not needed for controlled design
  m1 = 30, m2 = 30,   # Required for predictive probability
  ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL # No external data
)

print(decision_predictive)

## ----visualization, fig.width=10, fig.height=6, fig.alt="Comparison of posterior vs posterior predictive decision probabilities showing Go, Gray, and NoGo zones across different treatment effect values"----
# Create comprehensive decision probability data
# Split into separate scenarios to avoid vectorization issues
posterior_results <- BayesDecisionProbBinary(
  prob = 'posterior', design = 'controlled', 
  theta.TV = 0.3, theta.MAV = 0.1, theta.NULL = NULL, 
  gamma1 = 0.8, gamma2 = 0.2,
  pi1 = seq(0.1, 0.9, by = 0.05), pi2 = rep(0.2, length(seq(0.1, 0.9, by = 0.05))), 
  n1 = 12, n2 = 15, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  z = NULL, m1 = NULL, m2 = NULL,
  ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
) %>%
  mutate(prob = 'posterior', theta = pi1 - pi2)

predictive_results <- BayesDecisionProbBinary(
  prob = 'predictive', design = 'controlled', 
  theta.TV = NULL, theta.MAV = NULL, theta.NULL = 0.15, 
  gamma1 = 0.9, gamma2 = 0.3,
  pi1 = seq(0.1, 0.9, by = 0.05), pi2 = rep(0.2, length(seq(0.1, 0.9, by = 0.05))), 
  n1 = 12, n2 = 15, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  z = NULL, m1 = 30, m2 = 30,
  ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
) %>%
  mutate(prob = 'predictive', theta = pi1 - pi2)

# Combine results
results <- bind_rows(posterior_results, predictive_results) %>%
  mutate(Gray = pmax(Gray, 0))  # Ensure non-negative Gray probabilities

# Create visualization
figure <- results %>%
  pivot_longer(
    cols = c(Go, NoGo, Gray), names_to = 'Decision', values_to = 'Prob'
  ) %>%
  mutate(
    prob = factor(prob, levels = c('posterior', 'predictive')),
    Decision = factor(Decision, levels = c('Go', 'Gray', 'NoGo'))
  ) %>%
  ggplot(aes(x = theta, y = Prob)) +
  facet_wrap(~ prob, labeller = labeller(prob = c(
    'posterior' = 'Posterior Probability',
    'predictive' = 'Posterior Predictive Probability'
  ))) +
  geom_line(aes(colour = Decision, linetype = Decision), linewidth = 1.2) +
  theme_bw() +
  scale_color_manual(
    values = c('Go' = '#2E8B57', 'Gray' = '#808080', 'NoGo' = '#DC143C'),
    labels = c('Go', 'Gray', 'NoGo')
  ) +
  scale_linetype_manual(
    values = c('Go' = 'solid', 'Gray' = 'dashed', 'NoGo' = 'dotted'),
    labels = c('Go', 'Gray', 'NoGo')
  ) +
  scale_x_continuous(
    name = expression(paste("Treatment Effect (", pi[1] - pi[2], ")")),
    breaks = seq(-0.1, 0.7, by = 0.2)
  ) +
  scale_y_continuous(
    name = "Probability",
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1)
  ) +
  labs(
    title = "Go/Gray/NoGo Decision Probabilities for Binary Endpoints",
    subtitle = "Comparison of posterior vs. posterior predictive approaches"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 11),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 11)
  )

print(figure)

## ----uncontrolled-------------------------------------------------------------
# Uncontrolled design example
uncontrolled_result <- BayesDecisionProbBinary(
  prob = 'posterior', 
  design = 'uncontrolled', 
  theta.TV = 0.3, 
  theta.MAV = 0.1,
  theta.NULL = NULL, # Not needed for posterior probability
  gamma1 = 0.8, 
  gamma2 = 0.2,
  pi1 = c(0.3, 0.5, 0.7), 
  pi2 = rep(0.2, 3),  # Historical control rate
  n1 = 20, 
  n2 = 20,    # Provide n2 even for uncontrolled design
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  z = 3,      # Fixed number of historical control responders
  m1 = NULL, m2 = NULL, # Not needed for posterior probability
  ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL # No external data
)

print(uncontrolled_result)

## ----external_comparison------------------------------------------------------
# Compare different power prior weights
power_weights <- c(0.25, 0.5, 0.75, 1.0)
external_results <- purrr::map_dfr(power_weights, ~ {
  result <- BayesDecisionProbBinary(
    prob = 'posterior', 
    design = 'external', 
    theta.TV = 0.25, 
    theta.MAV = 0.05,
    theta.NULL = NULL, # Not needed for posterior probability
    gamma1 = 0.8, 
    gamma2 = 0.2,
    pi1 = c(0.3, 0.5), 
    pi2 = rep(0.2, 2), 
    n1 = 15, n2 = 15,
    a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
    z = NULL,           # Not needed for external design
    m1 = NULL, m2 = NULL, # Not needed for posterior probability
    ne1 = 20, ne2 = 20, ye1 = 8, ye2 = 4, 
    ae1 = .x, ae2 = .x
  )
  result$power_weight <- .x
  result
})

print(external_results)

## ----prior_sensitivity--------------------------------------------------------
# Compare different prior specifications
prior_specs <- list(
  "Vague (0.5, 0.5)" = list(a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5),
  "Uniform (1, 1)" = list(a1 = 1, a2 = 1, b1 = 1, b2 = 1),
  "Informative (2, 8)" = list(a1 = 2, a2 = 2, b1 = 8, b2 = 8)
)

sensitivity_results <- purrr::map_dfr(names(prior_specs), function(prior_name) {
  specs <- prior_specs[[prior_name]]
  
  result <- BayesPostPredBinary(
    prob = 'posterior', design = 'controlled', theta0 = 0.15,
    n1 = 20, n2 = 20, y1 = 12, y2 = 8,
    a1 = specs$a1, a2 = specs$a2, b1 = specs$b1, b2 = specs$b2,
    m1 = NULL, m2 = NULL,
    ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  )
  
  tibble(prior_type = prior_name, posterior_prob = result)
})

print(sensitivity_results)

## ----sample_size_comparison---------------------------------------------------
# Compare different sample sizes
sample_sizes <- c(10, 20, 30, 50)
oc_results <- purrr::map_dfr(sample_sizes, function(n) {
  result <- BayesDecisionProbBinary(
    prob = 'posterior', design = 'controlled', 
    theta.TV = 0.3, theta.MAV = 0.1, theta.NULL = NULL,
    gamma1 = 0.8, gamma2 = 0.2,
    pi1 = 0.5, pi2 = 0.2,  # True effect of 0.3
    n1 = n, n2 = n,
    a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
    z = NULL, m1 = NULL, m2 = NULL,
    ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  )
  result$sample_size <- n
  result
})

print(oc_results)

