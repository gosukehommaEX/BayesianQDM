# tibble(mu.t = c(1, 2, 4)) %>%
#   group_by_all() %>%
#   mutate(
#     Go.and.NoGo.prob.continuous(
#       nsim = 1000, msim = 1000, n.t = 12, n.c = 12, m.t = NULL, m.c = NULL, r = NULL,
#       mu.0t = 5, mu.0c = 5, k.0t = 5, k.0c = 5, nu.0t = 5, nu.0c = 5, sigma.0t = sqrt(5), sigma.0c = sqrt(5),
#       mu.t = mu.t, mu.c = 0, sigma.t = 1, sigma.c = 1,
#       theta0 = c(2, 0), gamma1 = 0.8, gamma2 = 0.3, design = 'controlled', prior = 'N-Inv-Chisq', prob.type = 'posterior', seed = 1
#     )
#   )
BayesDecisionProbContinuous = function(nsim, msim, n.t, n.c, m.t, m.c, r,
                                       mu.0t, mu.0c, k.0t, k.0c, nu.0t, nu.0c, sigma.0t, sigma.0c,
                                       mu.t, mu.c, sigma.t, sigma.c,
                                       theta0, gamma1, gamma2, design, prior, prob.type, seed) {
  # Set seed number
  set.seed(seed)
  # Random numbers of treatment group in PoC study
  y.t = matrix(rnorm(nsim * n.t, mu.t, sigma.t), nrow = nsim)
  # Sample mean of treatment group
  bar.y.t = rowSums(y.t) / n.t
  # Standard deviation of treatment group
  s.t = sqrt(rowSums((y.t - bar.y.t) ^ 2) / (n.t - 1))
  ## For controlled design
  if(design == 'controlled') {
    # Random numbers of control group in PoC study
    y.c = matrix(rnorm(nsim * n.c, mu.c, sigma.c), nrow = nsim)
    # Sample mean of control group
    bar.y.c = rowSums(y.c) / n.c
    # Standard deviation of control group
    s.c = sqrt(rowSums((y.c - bar.y.c) ^ 2) / (n.c - 1))
  } else if(design == 'uncontrolled') {
    # Sample mean of control group
    bar.y.c = NULL
    # Standard deviation of control group
    s.c = NULL
  }
  # Probability of success
  prob.success = sapply(seq(nsim), function(i) {
    h.theta0(
      msim, n.t, n.c, m.t, m.c, r,
      mu.0t, mu.0c, k.0t, k.0c, nu.0t, nu.0c, sigma.0t, sigma.0c,
      bar.y.t = bar.y.t[i], bar.y.c = bar.y.c[i], s.t = s.t[i], s.c = s.c[i],
      theta0, design, prior, prob.type, seed
    )
  })
  # Posterior probability (Note: theta0 should be two values, e.g., theta0 = c(2, 0))
  if(prob.type == 'posterior') {
    # Go, NoGo and Grey probabilities
    Go = sum(prob.success[1, ] >= gamma1) / msim
    NoGo = sum(prob.success[2, ] <= gamma2) / msim
    Grey = 1 - Go - NoGo
    results = data.frame(Go, NoGo, Grey)
    ## Predictive probability (Note: theta0 should be one value, e.g., theta0 = 0)
  } else if(prob.type == 'predictive') {
    # Go, NoGo and Grey probabilities
    Go = sum(prob.success >= gamma1) / msim
    NoGo = sum(prob.success <= gamma2) / msim
    Grey = 1 - Go - NoGo
    results = data.frame(Go, NoGo, Grey)
  }
  # Return output
  return(results)
}
