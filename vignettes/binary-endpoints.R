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
# Calculate posterior probability for a single trial outcome
posterior_prob <- pPPsinglebinary(
  prob = 'posterior',
  design = 'controlled',
  theta0 = 0.15,  # Minimum clinically meaningful difference
  n1 = 20, n2 = 20,  # Sample sizes
  y1 = 12, y2 = 6,   # Observed responders
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,  # Beta prior parameters
  m1 = NULL, m2 = NULL,
  ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
)

cat("Posterior probability that treatment effect > 0.15:", round(posterior_prob, 4))

## ----predictive_basic---------------------------------------------------------
# Calculate posterior predictive probability for future trial
predictive_prob <- pPPsinglebinary(
  prob = 'predictive',
  design = 'controlled',
  theta0 = 0.1,
  n1 = 20, n2 = 20,     # Current trial sample sizes
  y1 = 12, y2 = 6,      # Current trial outcomes
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  m1 = 100, m2 = 100,   # Future trial sample sizes
  ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
)

cat("Probability of future trial success:", round(predictive_prob, 4))

## ----operating_characteristics------------------------------------------------
# Define true response rates to evaluate
true_rates <- expand.grid(
  pi1 = c(0.2, 0.3, 0.4, 0.5, 0.6),
  pi2 = 0.2
)

# Calculate Go/NoGo/Gray probabilities
oc_results <- pGNGsinglebinary(
  prob = 'posterior',
  design = 'controlled',
  theta.TV = 0.25,   # Target value for Go
  theta.MAV = 0.05,  # Minimum acceptable value for NoGo
  theta.NULL = NULL,
  gamma1 = 0.8,      # Go threshold
  gamma2 = 0.2,      # NoGo threshold
  pi1 = true_rates$pi1,
  pi2 = true_rates$pi2,
  n1 = 30, n2 = 30,
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  z = NULL, m1 = NULL, m2 = NULL,
  ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
)

print(oc_results)

## ----decision_plot, fig.width=8, fig.height=5---------------------------------
# Reshape data for plotting
plot_data <- oc_results %>%
  select(pi1, pi2, Go, NoGo, Gray) %>%
  pivot_longer(cols = c(Go, NoGo, Gray), 
               names_to = "Decision", 
               values_to = "Probability")

# Create visualization
ggplot(plot_data, aes(x = pi1, y = Probability, color = Decision, group = Decision)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c('Go' = '#2E8B57', 'Gray' = '#808080', 'NoGo' = '#DC143C'),
    labels = c('Go', 'Gray', 'NoGo')
  ) +
  labs(
    title = "Operating Characteristics for Binary Endpoint Trial",
    subtitle = "Control response rate = 20%",
    x = "True Treatment Response Rate",
    y = "Decision Probability",
    color = "Decision"
  ) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_x_continuous(labels = scales::percent) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

## ----external_design----------------------------------------------------------
# Calculate probability with external control data
external_result <- pPPsinglebinary(
  prob = 'posterior',
  design = 'external',
  theta0 = 0.15,
  n1 = 20, n2 = 20,     # Current trial
  y1 = 12, y2 = 6,      # Current outcomes
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  m1 = NULL, m2 = NULL,
  ne1 = 30, ne2 = 30,   # Historical trial sizes
  ye1 = 15, ye2 = 8,    # Historical outcomes
  ae1 = 0.5, ae2 = 0.5  # Power prior parameters (50% borrowing)
)

cat("Posterior probability with external data:", round(external_result, 4))

## ----borrowing_sensitivity----------------------------------------------------
# Test different levels of borrowing
alpha_values <- seq(0, 1, by = 0.2)

borrowing_results <- data.frame(
  alpha = alpha_values,
  probability = sapply(alpha_values, function(a) {
    pPPsinglebinary(
      prob = 'posterior', design = 'external', theta0 = 0.15,
      n1 = 20, n2 = 20, y1 = 12, y2 = 6,
      a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
      m1 = NULL, m2 = NULL,
      ne1 = 30, ne2 = 30, ye1 = 15, ye2 = 8,
      ae1 = a, ae2 = a
    )
  })
)

print(borrowing_results)

## ----borrowing_plot, fig.width=7, fig.height=4--------------------------------
ggplot(borrowing_results, aes(x = alpha, y = probability)) +
  geom_line(linewidth = 1.2, color = '#2E8B57') +
  geom_point(size = 3, color = '#2E8B57') +
  labs(
    title = "Sensitivity to Historical Data Borrowing",
    x = "Power Prior Parameter (α)",
    y = "Posterior Probability",
    subtitle = "α = 0: no borrowing, α = 1: full borrowing"
  ) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

