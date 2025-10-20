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
library(purrr)

## ----posterior_basic----------------------------------------------------------
# Calculate posterior probability from trial data
posterior_prob <- pPPsinglecontinuous(
  prob = 'posterior',
  design = 'controlled',
  prior = 'vague',
  CalcMethod = 'NI',
  theta0 = 1.5,  # Minimum clinically meaningful difference
  n1 = 25, n2 = 25,
  bar.y1 = 5.2, bar.y2 = 3.5,  # Observed means
  s1 = 2.1, s2 = 1.9           # Observed SDs
)

cat("Posterior probability that treatment effect > 1.5:", round(posterior_prob, 4))

## ----predictive_basic---------------------------------------------------------
# Calculate probability of future trial success
predictive_prob <- pPPsinglecontinuous(
  prob = 'predictive',
  design = 'controlled',
  prior = 'vague',
  CalcMethod = 'WS',
  theta0 = 1.0,
  n1 = 25, n2 = 25,      # Current trial
  m1 = 100, m2 = 100,    # Future trial
  bar.y1 = 5.2, bar.y2 = 3.5,
  s1 = 2.1, s2 = 1.9
)

cat("Probability of future trial success:", round(predictive_prob, 4))

## ----method_comparison--------------------------------------------------------
# Compare NI and WS methods
cat("=== Method Comparison ===\n\n")

# NI method (exact)
time_ni <- system.time({
  prob_ni <- pPPsinglecontinuous(
    prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
    theta0 = 1.5, n1 = 20, n2 = 20, bar.y1 = 5, bar.y2 = 3, s1 = 2, s2 = 2
  )
})

cat("NI method:\n")
cat("  Probability:", round(prob_ni, 6), "\n")
cat("  Time:", round(time_ni[3], 4), "seconds\n\n")

# WS method (approximation)
time_ws <- system.time({
  prob_ws <- pPPsinglecontinuous(
    prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'WS',
    theta0 = 1.5, n1 = 20, n2 = 20, bar.y1 = 5, bar.y2 = 3, s1 = 2, s2 = 2
  )
})

cat("WS method:\n")
cat("  Probability:", round(prob_ws, 6), "\n")
cat("  Time:", round(time_ws[3], 4), "seconds\n\n")

cat("Difference:", round(abs(prob_ni - prob_ws), 6), "\n")
cat("Speed-up:", round(time_ni[3] / time_ws[3], 2), "x\n")

## ----operating_characteristics------------------------------------------------
# Evaluate decision probabilities across different scenarios
true_effects <- data.frame(
  mu1 = c(3, 4, 5, 6, 7),
  mu2 = 2
)

oc_results <- do.call(rbind, lapply(1:nrow(true_effects), function(i) {
  pGNGsinglecontinuous(
    nsim = 50,  # Increased for better estimates
    prob = 'posterior',
    design = 'controlled',
    prior = 'vague',
    CalcMethod = 'WS',
    theta.TV = 2.0,   # Target value for Go
    theta.MAV = 0.5,  # Minimum acceptable value for NoGo
    theta.NULL = NULL,
    nMC = NULL,
    gamma1 = 0.8,     # Go threshold
    gamma2 = 0.2,     # NoGo threshold
    n1 = 30, n2 = 30,
    m1 = NULL, m2 = NULL,
    kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
    mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
    mu1 = true_effects$mu1[i],
    mu2 = true_effects$mu2[i],
    sigma1 = 2, sigma2 = 2,
    r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
    bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
    seed = 123 + i
  )
}))

print(oc_results)

## ----decision_plot, fig.width=8, fig.height=5---------------------------------
# Calculate true treatment effect
oc_results$true_effect <- oc_results$mu1 - oc_results$mu2

# Reshape for plotting
plot_data <- oc_results %>%
  select(true_effect, Go, NoGo, Gray) %>%
  pivot_longer(cols = c(Go, NoGo, Gray), 
               names_to = "Decision", 
               values_to = "Probability")

# Create visualization
ggplot(plot_data, aes(x = true_effect, y = Probability, color = Decision)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c('Go' = '#2E8B57', 'Gray' = '#808080', 'NoGo' = '#DC143C')
  ) +
  geom_vline(xintercept = 2.0, linetype = "dashed", alpha = 0.5, color = "blue") +
  geom_vline(xintercept = 0.5, linetype = "dashed", alpha = 0.5, color = "red") +
  annotate("text", x = 2.0, y = 0.95, label = "θ_TV", color = "blue") +
  annotate("text", x = 0.5, y = 0.95, label = "θ_MAV", color = "red") +
  labs(
    title = "Operating Characteristics for Continuous Endpoint Trial",
    subtitle = "Sample size = 30 per arm, σ = 2",
    x = "True Treatment Effect (μ₁ - μ₂)",
    y = "Decision Probability",
    color = "Decision"
  ) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

## ----prior_comparison---------------------------------------------------------
# Vague prior
prob_vague <- pPPsinglecontinuous(
  prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
  theta0 = 1.5, n1 = 15, n2 = 15, bar.y1 = 5, bar.y2 = 3, s1 = 2, s2 = 2
)

# Informative prior (Normal-Inverse-Chi-squared)
prob_informative <- pPPsinglecontinuous(
  prob = 'posterior', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'NI',
  theta0 = 1.5, n1 = 15, n2 = 15,
  kappa01 = 5, kappa02 = 5, nu01 = 10, nu02 = 10,
  mu01 = 5, mu02 = 3, sigma01 = 2, sigma02 = 2,
  bar.y1 = 5, bar.y2 = 3, s1 = 2, s2 = 2
)

cat("Vague prior:       ", round(prob_vague, 4), "\n")
cat("Informative prior: ", round(prob_informative, 4), "\n")

## ----external_design----------------------------------------------------------
# Incorporating external control data
external_result <- pPPsinglecontinuous(
  prob = 'posterior',
  design = 'external',
  prior = 'vague',
  CalcMethod = 'WS',
  theta0 = 1.5,
  n1 = 20, n2 = 20,           # Current trial
  bar.y1 = 5.5, bar.y2 = 3.2,
  s1 = 2.1, s2 = 1.9,
  ne1 = NULL, ne2 = 40,       # External control only
  alpha01 = NULL, alpha02 = 0.5,  # 50% borrowing
  bar.ye1 = NULL, bar.ye2 = 3.0,  # Historical control mean
  se1 = NULL, se2 = 2.0           # Historical control SD
)

cat("Posterior probability with external control:", round(external_result, 4))

## ----power_prior_sensitivity--------------------------------------------------
# Test different levels of borrowing
alpha_values <- seq(0, 1, by = 0.2)

borrowing_results <- data.frame(
  alpha = alpha_values,
  probability = sapply(alpha_values, function(a) {
    pPPsinglecontinuous(
      prob = 'posterior', design = 'external', prior = 'vague', CalcMethod = 'WS',
      theta0 = 1.5, n1 = 20, n2 = 20, bar.y1 = 5.5, bar.y2 = 3.2, s1 = 2.1, s2 = 1.9,
      ne1 = NULL, ne2 = 40, alpha01 = NULL, alpha02 = a,
      bar.ye1 = NULL, bar.ye2 = 3.0, se1 = NULL, se2 = 2.0
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

## ----power_analysis-----------------------------------------------------------
# Evaluate optimal sample size for desired operating characteristics
sample_sizes <- seq(10, 50, by = 10)

oc_results_n <- do.call(rbind, lapply(sample_sizes, function(n) {
  result <- pGNGsinglecontinuous(
    nsim = 30,  # Reduced for vignette speed
    prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'WS',
    theta.TV = 1.5, theta.MAV = 0.5, theta.NULL = NULL,
    nMC = NULL, gamma1 = 0.8, gamma2 = 0.3,
    n1 = n, n2 = n, m1 = NULL, m2 = NULL,
    kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
    mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
    mu1 = 5, mu2 = 3,  # Assume true effect of 2.0
    sigma1 = 2, sigma2 = 2,
    r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
    bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
    seed = 123
  )
  
  data.frame(
    sample_size = n,
    go_prob = result$Go,
    nogo_prob = result$NoGo,
    gray_prob = result$Gray
  )
}))

print(oc_results_n)

## ----power_plot, fig.width=8, fig.height=5------------------------------------
# Reshape for plotting
power_plot_data <- oc_results_n %>%
  pivot_longer(cols = c(go_prob, gray_prob, nogo_prob),
               names_to = "Decision",
               values_to = "Probability") %>%
  mutate(Decision = recode(Decision,
                           "go_prob" = "Go",
                           "gray_prob" = "Gray",
                           "nogo_prob" = "NoGo"))

ggplot(power_plot_data, aes(x = sample_size, y = Probability, color = Decision)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c('Go' = '#2E8B57', 'Gray' = '#808080', 'NoGo' = '#DC143C')
  ) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "blue", alpha = 0.5) +
  annotate("text", x = 45, y = 0.85, label = "Target: 80%", color = "blue") +
  labs(
    title = "Power Analysis: Sample Size vs Decision Probabilities",
    x = "Sample Size (per arm)",
    y = "Decision Probability",
    color = "Decision",
    subtitle = "Assuming true treatment effect = 2.0, σ = 2"
  ) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

