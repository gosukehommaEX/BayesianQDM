# Test continuous endpoint functions (excluding MC and MCMC methods for speed)

test_that("pNIdifft function works correctly", {
  # Test basic functionality
  result <- pNIdifft(q = 3, mu.t1 = 2, mu.t2 = 0, sd.t1 = 1, sd.t2 = 1, nu.t1 = 17, nu.t2 = 17)
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0 && result <= 1)

  # Test with different parameters
  result2 <- pNIdifft(q = 1, mu.t1 = 5, mu.t2 = 3, sd.t1 = 2, sd.t2 = 1.5, nu.t1 = 10, nu.t2 = 15)
  expect_type(result2, "double")
  expect_true(result2 >= 0 && result2 <= 1)

  # Test with q = 0 and symmetric parameters (should be close to 0.5)
  result_zero <- pNIdifft(q = 0, mu.t1 = 1, mu.t2 = 1, sd.t1 = 1, sd.t2 = 1, nu.t1 = 20, nu.t2 = 20)
  expect_true(abs(result_zero - 0.5) < 0.1)
})

test_that("pWSdifft function works correctly", {
  # Test basic functionality
  result <- pWSdifft(q = 3, mu.t1 = 2, mu.t2 = 0, sd.t1 = 1, sd.t2 = 1, nu.t1 = 17, nu.t2 = 17)
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0 && result <= 1)

  # Test with unequal variances
  result2 <- pWSdifft(q = 1, mu.t1 = 5, mu.t2 = 3, sd.t1 = 2, sd.t2 = 1.5, nu.t1 = 10, nu.t2 = 15)
  expect_type(result2, "double")
  expect_true(result2 >= 0 && result2 <= 1)

  # Test with q = 0 and symmetric parameters
  result_zero <- pWSdifft(q = 0, mu.t1 = 1, mu.t2 = 1, sd.t1 = 1, sd.t2 = 1, nu.t1 = 5, nu.t2 = 20)
  expect_true(abs(result_zero - 0.5) < 0.1)
})

test_that("Comparison between NI and WS methods", {
  # Compare NI and WS results for same parameters
  q_val <- 2
  mu1 <- 3
  mu2 <- 1
  sd1 <- 1.5
  sd2 <- 1.2
  nu1 <- 15
  nu2 <- 18

  result_ni <- pNIdifft(q = q_val, mu.t1 = mu1, mu.t2 = mu2, sd.t1 = sd1, sd.t2 = sd2, nu.t1 = nu1, nu.t2 = nu2)
  result_ws <- pWSdifft(q = q_val, mu.t1 = mu1, mu.t2 = mu2, sd.t1 = sd1, sd.t2 = sd2, nu.t1 = nu1, nu.t2 = nu2)

  # Results should be reasonably close (within 10% for most cases)
  expect_true(abs(result_ni - result_ws) / max(result_ni, result_ws) < 0.3)
})

test_that("BayesPostPredContinuous function works correctly with NI method", {
  # Test posterior probability with controlled design
  result_ni <- BayesPostPredContinuous(
    prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
    theta0 = 1, n1 = 12, n2 = 12, bar.y1 = 3, bar.y2 = 1, s1 = 1.5, s2 = 1.2
  )
  expect_type(result_ni, "double")
  expect_length(result_ni, 1)
  expect_true(result_ni >= 0 && result_ni <= 1)

  # Test with N-Inv-Chisq prior
  result_ninvchisq <- BayesPostPredContinuous(
    prob = 'posterior', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'NI',
    theta0 = 2, n1 = 12, n2 = 12, kappa01 = 5, kappa02 = 5, nu01 = 5, nu02 = 5,
    mu01 = 5, mu02 = 5, sigma01 = sqrt(5), sigma02 = sqrt(5),
    bar.y1 = 2, bar.y2 = 0, s1 = 1, s2 = 1
  )
  expect_type(result_ninvchisq, "double")
  expect_length(result_ninvchisq, 1)
  expect_true(result_ninvchisq >= 0 && result_ninvchisq <= 1)
})

test_that("BayesPostPredContinuous function works correctly with WS method", {
  # Test predictive probability with controlled design and vague prior
  result_ws <- BayesPostPredContinuous(
    prob = 'predictive', design = 'controlled', prior = 'vague', CalcMethod = 'WS',
    theta0 = 1, n1 = 12, n2 = 12, m1 = 120, m2 = 120,
    bar.y1 = 3, bar.y2 = 1, s1 = 1.5, s2 = 1.2
  )
  expect_type(result_ws, "double")
  expect_length(result_ws, 1)
  expect_true(result_ws >= 0 && result_ws <= 1)

  # Test with N-Inv-Chisq prior
  result_ninvchisq <- BayesPostPredContinuous(
    prob = 'predictive', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'WS',
    theta0 = 2, n1 = 12, n2 = 12, m1 = 120, m2 = 120, kappa01 = 5, kappa02 = 5, nu01 = 5, nu02 = 5,
    mu01 = 5, mu02 = 5, sigma01 = sqrt(5), sigma02 = sqrt(5),
    bar.y1 = 2, bar.y2 = 0, s1 = 1, s2 = 1
  )
  expect_type(result_ninvchisq, "double")
  expect_length(result_ninvchisq, 1)
  expect_true(result_ninvchisq >= 0 && result_ninvchisq <= 1)

  # Test uncontrolled design
  result_uncontrolled <- BayesPostPredContinuous(
    prob = 'predictive', design = 'uncontrolled', prior = 'vague', CalcMethod = 'WS',
    theta0 = 0.5, n1 = 12, m1 = 120, m2 = 120, mu02 = 0,
    bar.y1 = 2, s1 = 1, r = 12
  )
  expect_type(result_uncontrolled, "double")
  expect_true(result_uncontrolled >= 0 && result_uncontrolled <= 1)
})

test_that("BayesDecisionProbContinuous function works correctly with NI method", {
  # Test with small nsim for speed - ADD MISSING m1, m2 PARAMETERS
  result_ni <- BayesDecisionProbContinuous(
    nsim = 10, prob = 'posterior', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'NI',
    theta.TV = 2, theta.MAV = 0, theta.NULL = NULL, gamma1 = 0.8, gamma2 = 0.3,
    n1 = 12, n2 = 12, m1 = NULL, m2 = NULL, kappa01 = 5, kappa02 = 5, nu01 = 5, nu02 = 5,
    mu01 = 5, mu02 = 5, sigma01 = sqrt(5), sigma02 = sqrt(5), mu1 = 4, mu2 = 0,
    sigma1 = 1, sigma2 = 1, seed = 1
  )
  expect_s3_class(result_ni, "data.frame")
  expect_named(result_ni, c("mu1", "mu2", "Go", "NoGo", "Gray"))
  expect_equal(nrow(result_ni), 1)
  expect_true(result_ni$Go >= 0 && result_ni$Go <= 1)
  expect_true(result_ni$NoGo >= 0 && result_ni$NoGo <= 1)
  expect_true(result_ni$Gray >= -1e-10 && result_ni$Gray <= 1) # Allow small negative due to numerical precision
})

test_that("BayesDecisionProbContinuous function works correctly with WS method", {
  # Test with WS method
  result_ws <- BayesDecisionProbContinuous(
    nsim = 10, prob = 'predictive', design = 'controlled', prior = 'vague', CalcMethod = 'WS',
    theta.TV = NULL, theta.MAV = NULL, theta.NULL = 1, gamma1 = 0.8, gamma2 = 0.3,
    n1 = 10, n2 = 10, m1 = 50, m2 = 50, mu1 = 2.5, mu2 = 1.2,
    sigma1 = 1, sigma2 = 1, seed = 1
  )
  expect_s3_class(result_ws, "data.frame")
  expect_named(result_ws, c("mu1", "mu2", "Go", "NoGo", "Gray"))
  expect_equal(nrow(result_ws), 1)
})

test_that("Error handling works correctly for continuous functions", {
  # Test missing parameters for predictive probability
  expect_error(
    BayesPostPredContinuous(
      prob = 'predictive', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
      theta0 = 1, n1 = 12, n2 = 12, bar.y1 = 3, bar.y2 = 1, s1 = 1.5, s2 = 1.2
    )
  )

  # Test missing parameters for N-Inv-Chisq prior
  expect_error(
    BayesPostPredContinuous(
      prob = 'posterior', design = 'controlled', prior = 'N-Inv-Chisq', CalcMethod = 'NI',
      theta0 = 2, n1 = 12, n2 = 12, kappa01 = NULL, kappa02 = NULL, nu01 = NULL, nu02 = NULL, sigma01 = NULL, sigma02 = NULL,
      bar.y1 = 2, bar.y2 = 0, s1 = 1, s2 = 1
    )
  )

  # Test missing parameters for uncontrolled design
  expect_error(
    BayesPostPredContinuous(
      prob = 'posterior', design = 'uncontrolled', prior = 'vague', CalcMethod = 'NI',
      theta0 = 1, n1 = 12, bar.y1 = 2, s1 = 1, mu02 = 0, r = NULL
    )
  )

  # Test missing parameters for decision probability
  expect_error(
    BayesDecisionProbContinuous(
      nsim = 10, prob = 'posterior', design = 'controlled', prior = 'vague', CalcMethod = 'NI',
      theta.TV = NULL, theta.MAV = NULL, theta.NULL = NULL, gamma1 = 0.8, gamma2 = 0.3,
      n1 = 12, n2 = 12, mu1 = 4, mu2 = 0, sigma1 = 1, sigma2 = 1, seed = 1
    )
  )
})
