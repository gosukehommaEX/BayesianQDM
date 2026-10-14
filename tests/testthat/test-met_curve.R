test_that("met_curve() gives the probability of meeting a threshold", {
  p <- c(0.2, 0.5, 0.8)
  w <- c(0.5, 0.3, 0.2)
  expect_equal(met_curve(p, w, c(0.1, 0.2, 0.3, 0.5, 0.9)),
               c(1, 1, 0.5, 0.5, 0))
  # Ties within the tolerance meet the threshold
  expect_equal(met_curve(0.5 - 1e-12, 1, 0.5), 1)
})

test_that("met_curve() handles pairs of outcomes and marginal rules", {
  # Probabilities over 2 x 2 pairs of treatment and control outcomes
  p <- matrix(c(0.2, 0.6, 0.4, 0.8), 2L)
  w <- list(t = c(0.3, 0.7), c = c(0.5, 0.5))
  expect_equal(met_curve(p, w, c(0.1, 0.5, 0.7)),
               c(1, 0.7 * 0.5 + 0.7 * 0.5, 0.7 * 0.5))
  # Marginal rules: the same threshold for both endpoints
  pl <- list(c(0.2, 0.6, 0.9), c(0.7, 0.5, 0.4))
  wv <- c(0.2, 0.3, 0.5)
  expect_equal(met_curve(pl, wv, c(0.3, 0.5), "and"), c(0.3 + 0.5, 0.3))
  expect_equal(met_curve(pl, wv, c(0.3, 0.5), "or"), c(1, 1))
})
