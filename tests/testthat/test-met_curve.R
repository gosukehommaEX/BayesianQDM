test_that("met_curve() gives the probability of meeting a threshold", {
  p <- c(0.2, 0.5, 0.8)
  w <- c(0.5, 0.3, 0.2)
  expect_equal(met_curve(p, w, c(0.1, 0.2, 0.3, 0.5, 0.9)),
               c(1, 1, 0.5, 0.5, 0))
  # Ties within the tolerance meet the threshold
  expect_equal(met_curve(0.5 - 1e-12, 1, 0.5), 1)
})
