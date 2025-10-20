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

test_that("d2betadiff function works correctly", {
  # Test basic functionality
  result <- d2betadiff(0.2, 0.5, 0.5, 0.5, 0.5)
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0)

  # Test with different parameters
  result2 <- d2betadiff(-0.1, 2, 3, 1, 4)
  expect_type(result2, "double")
  expect_true(result2 >= 0)
})

test_that("p2betadiff function works correctly", {
  # Test basic functionality with symmetric parameters
  result <- p2betadiff(0.2, 0.5, 0.5, 0.5, 0.5)
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0 && result <= 1)

  # Test with q = 0 (should give 0.5 for symmetric case)
  result_zero <- p2betadiff(0, 1, 1, 1, 1)
  expect_equal(result_zero, 0.5, tolerance = 1e-6)

  # Test with different parameters
  result2 <- p2betadiff(-0.1, 2, 3, 1, 4)
  expect_type(result2, "double")
  expect_true(result2 >= 0 && result2 <= 1)
})

test_that("p2betabinomdiff function works correctly", {
  # Test basic functionality
  result <- p2betabinomdiff(0.2, 12, 12, 0.5, 0.5, 0.5, 0.5)
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0 && result <= 1)

  # Test with different sample sizes
  result2 <- p2betabinomdiff(0.1, 20, 15, 1, 1, 1, 1)
  expect_type(result2, "double")
  expect_true(result2 >= 0 && result2 <= 1)

  # Test with q = 0 for symmetric case
  result_zero <- p2betabinomdiff(0, 10, 10, 2, 2, 3, 3)
  expect_type(result_zero, "double")
  expect_true(result_zero >= 0 && result_zero <= 1)
})

test_that("pPPsinglebinary function works correctly", {
  # Test posterior probability with controlled design
  result_post <- pPPsinglebinary(
    prob = 'posterior', design = 'controlled', theta0 = 0.15,
    n1 = 12, n2 = 15, y1 = 7, y2 = 9, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
    m1 = NULL, m2 = NULL, ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  )
  expect_type(result_post, "double")
  expect_length(result_post, 1)
  expect_true(result_post >= 0 && result_post <= 1)

  # Test predictive probability with controlled design
  result_pred <- pPPsinglebinary(
    prob = 'predictive', design = 'controlled', theta0 = 0.5,
    n1 = 12, n2 = 15, y1 = 7, y2 = 7, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
    m1 = 12, m2 = 12, ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  )
  expect_type(result_pred, "double")
  expect_length(result_pred, 1)
  expect_true(result_pred >= 0 && result_pred <= 1)

  # Test external design
  result_ext <- pPPsinglebinary(
    prob = 'posterior', design = 'external', theta0 = 0.15,
    n1 = 12, n2 = 15, y1 = 7, y2 = 9, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
    m1 = NULL, m2 = NULL, ne1 = 12, ne2 = 12, ye1 = 6, ye2 = 6, ae1 = 0.5, ae2 = 0.5
  )
  expect_type(result_ext, "double")
  expect_length(result_ext, 1)
  expect_true(result_ext >= 0 && result_ext <= 1)
})

test_that("pGNGsinglebinary function works correctly", {
  # Test posterior probability decision making with small n for speed
  result_post <- pGNGsinglebinary(
    prob = 'posterior', design = 'controlled', theta.TV = 0.4, theta.MAV = 0.2, theta.NULL = NULL,
    gamma1 = 0.5, gamma2 = 0.2, pi1 = c(0.3, 0.6), pi2 = rep(0.2, 2), n1 = 10, n2 = 10,
    a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5, z = NULL, m1 = NULL, m2 = NULL, ne1 = NULL, ne2 = NULL,
    ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  )
  expect_s3_class(result_post, "data.frame")
  expect_named(result_post, c("pi1", "pi2", "Go", "Gray", "NoGo", "Miss"))
  expect_equal(nrow(result_post), 2)
  expect_true(all(result_post$Go >= 0 & result_post$Go <= 1))
  expect_true(all(result_post$NoGo >= 0 & result_post$NoGo <= 1))

  # Test predictive probability decision making
  result_pred <- pGNGsinglebinary(
    prob = 'predictive', design = 'controlled', theta.TV = NULL, theta.MAV = NULL, theta.NULL = 0,
    gamma1 = 0.9, gamma2 = 0.3, pi1 = c(0.3, 0.5), pi2 = rep(0.2, 2), n1 = 10, n2 = 10,
    a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5, z = NULL, m1 = 20, m2 = 20, ne1 = NULL, ne2 = NULL,
    ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
  )
  expect_s3_class(result_pred, "data.frame")
  expect_named(result_pred, c("pi1", "pi2", "Go", "Gray", "NoGo", "Miss"))
  expect_equal(nrow(result_pred), 2)
})

test_that("Error handling works correctly for binary functions", {
  # Test missing parameters for predictive probability
  expect_error(
    pPPsinglebinary(
      prob = 'predictive', design = 'controlled', theta0 = 0.5,
      n1 = 12, n2 = 15, y1 = 7, y2 = 7, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
      m1 = NULL, m2 = NULL, ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
    )
  )

  # Test missing parameters for external design
  expect_error(
    pPPsinglebinary(
      prob = 'posterior', design = 'external', theta0 = 0.15,
      n1 = 12, n2 = 15, y1 = 7, y2 = 9, a1 = 0.5, a2 = 0.5, b1 = 0.5, b2 = 0.5,
      m1 = NULL, m2 = NULL, ne1 = NULL, ne2 = NULL, ye1 = NULL, ye2 = NULL, ae1 = NULL, ae2 = NULL
    )
  )
})
