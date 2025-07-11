## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5
)

## ----setup--------------------------------------------------------------------
library(BayesianQDM)

## ----binary_example-----------------------------------------------------------
# Calculate decision probabilities for binary endpoint
result_binary <- BayesDecisionProbBinary(
  prob = 'posterior', 
  design = 'controlled', 
  theta.TV = 0.3, theta.MAV = 0.1,
  gamma1 = 0.8, gamma2 = 0.2,
  pi1 = c(0.3, 0.5, 0.7), 
  pi2 = rep(0.2, 3), 
  n1 = 15, n2 = 15,
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  z = NULL, m1 = NULL, m2 = NULL,
  ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
)

print(result_binary)

## ----binary_posterior---------------------------------------------------------
# Calculate posterior probability
posterior_prob <- BayesPostPredBinary(
  prob = 'posterior', 
  design = 'controlled', 
  theta0 = 0.15,
  n1 = 12, n2 = 15, y1 = 7, y2 = 5, 
  a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
  m1 = NULL, m2 = NULL,
  ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
)

cat("Posterior probability:", round(posterior_prob, 4))

## ----continuous_example-------------------------------------------------------
# Calculate decision probabilities for continuous endpoint  
result_continuous <- BayesDecisionProbContinuous(
  nsim = 50,  # Small nsim for vignette speed
  prob = 'posterior', 
  design = 'controlled', 
  prior = 'vague', 
  CalcMethod = 'WS',  # Fast method for vignette
  theta.TV = 1.5, theta.MAV = 0.5, theta.NULL = NULL,  # All required thresholds
  nMC = NULL, nMCMCsample = NULL,  # Not needed for WS method
  gamma1 = 0.8, gamma2 = 0.3,
  n1 = 12, n2 = 12, m1 = NULL, m2 = NULL,
  mu1 = 4.5, mu2 = 2.0, sigma1 = 1.5, sigma2 = 1.3,
  seed = 123
)

print(result_continuous)

