# Test binary endpoint functions

test_that("AppellsF1 function works correctly", {
  # Test basic functionality
  result <- AppellsF1(0.5, 0.5, 0, 1, 0.96, 1.2)
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(is.finite(result))

  # Test with known values
  result2 <- AppellsF1(1, 1, 1, 2, 0.3, 0.4)
  expect_type(result2, "double")
  expect_true(result2 > 0)
})

test_that("pBetadiff function works correctly", {
  # Test basic functionality with symmetric parameters
  result <- pBetadiff(0.2, 0.5, 0.5, 0.5, 0.5)
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0 && result <= 1)

  # Test with q = 0 (should give 0.5 for symmetric case)
  result_zero <- pBetadiff(0, 1, 1, 1, 1)
  expect_equal(result_zero, 0.5, tolerance = 1e-6)

  # Test with different parameters
  result2 <- pBetadiff(-0.1, 2, 3, 1, 4)
  expect_type(result2, "double")
  expect_true(result2 >= 0 && result2 <= 1)
})

test_that("pBetaBinomdiff function works correctly", {
  # Test basic functionality
  result <- pBetaBinomdiff(0.2, 12, 12, 0.5, 0.5, 0.5, 0.5)
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0 && result <= 1)

  # Test with different sample sizes
  result2 <- pBetaBinomdiff(0.1, 20, 15, 1, 1, 1, 1)
  expect_type(result2, "double")
  expect_true(result2 >= 0 && result2 <= 1)

  # Test with q = 0 for symmetric case
  result_zero <- pBetaBinomdiff(0, 10, 10, 2, 2, 3, 3)
  expect_type(result_zero, "double")
  expect_true(result_zero >= 0 && result_zero <= 1)
})

test_that("BayesPostPredBinary function works correctly", {
  # Test posterior probability with controlled design
  result_post <- BayesPostPredBinary(
    prob = 'posterior', design = 'controlled', theta0 = 0.15,
    n1 = 12, n2 = 15, y1 = 7, y2 = 9, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
    m1 = NULL, m2 = NULL, ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  )
  expect_type(result_post, "double")
  expect_length(result_post, 1)
  expect_true(result_post >= 0 && result_post <= 1)

  # Test predictive probability with controlled design
  result_pred <- BayesPostPredBinary(
    prob = 'predictive', design = 'controlled', theta0 = 0.5,
    n1 = 12, n2 = 15, y1 = 7, y2 = 7, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
    m1 = 12, m2 = 12, ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  )
  expect_type(result_pred, "double")
  expect_length(result_pred, 1)
  expect_true(result_pred >= 0 && result_pred <= 1)

  # Test external design
  result_ext <- BayesPostPredBinary(
    prob = 'posterior', design = 'external', theta0 = 0.15,
    n1 = 12, n2 = 15, y1 = 7, y2 = 9, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
    m1 = NULL, m2 = NULL, ne1 = 12, ne2 = 12, ye1 = 6, ye2 = 6, ae1 = 0.5, ae2 = 0.5
  )
  expect_type(result_ext, "double")
  expect_length(result_ext, 1)
  expect_true(result_ext >= 0 && result_ext <= 1)
})

test_that("BayesDecisionProbBinary function works correctly", {
  # Test posterior probability decision making
  result_post <- BayesDecisionProbBinary(
    prob = 'posterior', design = 'controlled', theta.TV = 0.4, theta.MAV = 0.2, theta.NULL = NULL,
    gamma1 = 0.5, gamma2 = 0.2, pi1 = c(0.2, 0.4, 0.6, 0.8), pi2 = rep(0.2, 4), n1 = 12, n2 = 12,
    a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5, z = NULL, m1 = NULL, m2 = NULL, ne1 = NULL, ne2 = NULL,
    ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  )
  expect_s3_class(result_post, "data.frame")
  expect_named(result_post, c("pi1", "pi2", "Go", "Gray", "NoGo"))
  expect_equal(nrow(result_post), 4)
  expect_true(all(result_post$Go >= 0 & result_post$Go <= 1))
  expect_true(all(result_post$NoGo >= 0 & result_post$NoGo <= 1))
  expect_true(all(result_post$Gray >= -1e-10 & result_post$Gray <= 1)) # Allow small negative due to numerical precision

  # Test predictive probability decision making
  result_pred <- BayesDecisionProbBinary(
    prob = 'predictive', design = 'controlled', theta.TV = NULL, theta.MAV = NULL, theta.NULL = 0,
    gamma1 = 0.9, gamma2 = 0.3, pi1 = c(0.2, 0.4), pi2 = rep(0.2, 2), n1 = 12, n2 = 12,
    a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5, z = NULL, m1 = 30, m2 = 30, ne1 = NULL, ne2 = NULL,
    ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  )
  expect_s3_class(result_pred, "data.frame")
  expect_named(result_pred, c("pi1", "pi2", "Go", "Gray", "NoGo"))
  expect_equal(nrow(result_pred), 2)
})

test_that("Error handling works correctly for binary functions", {
  # Test missing parameters for predictive probability
  expect_error(
    BayesPostPredBinary(
      prob = 'predictive', design = 'controlled', theta0 = 0.5,
      n1 = 12, n2 = 15, y1 = 7, y2 = 7, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
      m1 = NULL, m2 = NULL, ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
    ),
    "m1 and m2 should be non-null"
  )

  # Test missing parameters for external design
  expect_error(
    BayesPostPredBinary(
      prob = 'posterior', design = 'external', theta0 = 0.15,
      n1 = 12, n2 = 15, y1 = 7, y2 = 9, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
      m1 = NULL, m2 = NULL, ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
    ),
    "ne1, ne2, ye1, ye2, ae1 and ae2 should be non-null"
  )

  # Test missing parameters for decision probability
  expect_error(
    BayesDecisionProbBinary(
      prob = 'posterior', design = 'controlled', theta.TV = NULL, theta.MAV = 0.2, theta.NULL = NULL,
      gamma1 = 0.5, gamma2 = 0.2, pi1 = c(0.2, 0.4), pi2 = rep(0.2, 2), n1 = 12, n2 = 12,
      a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5, z = NULL, m1 = NULL, m2 = NULL, ne1 = NULL, ne2 = NULL,
      ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
    ),
    "theta.TV and theta.MAV should be non-null"
  )
})
