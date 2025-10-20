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
result_binary <- pGNGsinglebinary(
  prob = 'posterior', 
  design = 'controlled', 
  theta.TV = 0.3, theta.MAV = 0.1, theta.NULL = NULL,
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
posterior_prob <- pPPsinglebinary(
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
result_continuous <- pGNGsinglecontinuous(
  nsim = 50,  # Small nsim for vignette speed
  prob = 'posterior', 
  design = 'controlled', 
  prior = 'vague', 
  CalcMethod = 'WS',  # Fast method for vignette
  theta.TV = 1.5, theta.MAV = 0.5, theta.NULL = NULL,
  nMC = NULL, gamma1 = 0.8, gamma2 = 0.3,
  n1 = 12, n2 = 12, m1 = NULL, m2 = NULL,
  kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL,
  mu01 = NULL, mu02 = NULL, sigma01 = NULL, sigma02 = NULL,
  mu1 = 4.5, mu2 = 2.0, sigma1 = 1.5, sigma2 = 1.3,
  r = NULL, ne1 = NULL, ne2 = NULL, alpha01 = NULL, alpha02 = NULL,
  bar.ye1 = NULL, bar.ye2 = NULL, se1 = NULL, se2 = NULL,
  seed = 123
)

print(result_continuous)

## ----method_comparison--------------------------------------------------------
# Compare NI and WS methods
prob_ni <- pPPsinglecontinuous(
  prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
  theta0 = 1, n1 = 12, n2 = 12, bar.y1 = 3, bar.y2 = 1, s1 = 1.5, s2 = 1.2
)

prob_ws <- pPPsinglecontinuous(
  prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'WS',
  theta0 = 1, n1 = 12, n2 = 12, bar.y1 = 3, bar.y2 = 1, s1 = 1.5, s2 = 1.2
)

cat("NI method:", round(prob_ni, 4), "\n")
cat("WS method:", round(prob_ws, 4), "\n")
cat("Difference:", round(abs(prob_ni - prob_ws), 4), "\n")

